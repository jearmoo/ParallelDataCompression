/*
 * ByteReader.cpp
 *
 * author: tto@andrew.cmu.edu and cdaruwal@andrew.cmu.edu
 *
 */
#include "ByteReader.h"

#include "util.h"

#include <fstream>
#include <iostream>
#include <string>

#include <assert.h>
#include <omp.h>
#include <stdlib.h>

using namespace std;

bool ByteReader::load_from_file(const string &file_name) {
    ifstream infile_stream(file_name, ios::in | ios::binary | ios::ate);

    streampos infile_size;

    if (infile_stream.is_open()) {
        infile_size = infile_stream.tellg();

        data = new unsigned char[infile_size];
        size = infile_size;

        infile_stream.seekg(0, ios::beg);
        infile_stream.read((char *)data, infile_size);

        infile_stream.close();
        return true;
    }

    return false;
}

bool ByteReader::load_from_file_parallel(const string &file_name, int num_threads) {
    omp_set_num_threads(min(num_threads,omp_get_max_threads()));

    ifstream infile_stream(file_name, ios::in | ios::binary | ios::ate);
    streampos infile_size;

    if (infile_stream.is_open()) {
        infile_size = infile_stream.tellg();
        infile_stream.close();

        data = new unsigned char[infile_size];
        size = infile_size;

#pragma omp parallel
        {
            int num_threads = omp_get_num_threads();
            uint64_t chunk_size = get_chunksize(infile_size, num_threads);
            int thread_num = omp_get_thread_num();

            ifstream thread_infile(file_name, ios::in | ios::binary);

            uint64_t thread_start_pos = chunk_size * thread_num;

            if (thread_infile.is_open()) {
                thread_infile.seekg(thread_start_pos, ios::beg);
                thread_infile.read((char *)data + thread_start_pos,
                                   min((thread_num + 1) * chunk_size,
                                       static_cast<uint64_t>(infile_size)) -
                                       thread_start_pos);
                thread_infile.close();
            }
        }
        return true;
    }

    return false;
}

bool ByteReader::read(void *buf, uint64_t num_bytes, uint64_t manual_pos) {
    if (manual_pos != -1) {
        if (!(0 <= manual_pos && manual_pos + num_bytes <= size))
            return false;
        for (uint64_t i = 0; i < num_bytes; ++i) {
            static_cast<unsigned char *>(buf)[i] = data[manual_pos++];
        }
        return true;
    }
    if (pos + num_bytes > size) {
        return false;
    }
    for (uint64_t i = 0; i < num_bytes; ++i) {
        static_cast<unsigned char *>(buf)[i] = data[pos++];
    }
    return true;
}

void ByteReader::seek(uint64_t new_pos) {
    D(assert(0 <= new_pos && new_pos < size);)
    pos = new_pos;
}

void ByteReader::rewind() { pos = 0; }

bool ByteReader::at_end_of_data() { return pos == size; }

uint64_t ByteReader::get_size() { return size; }

unsigned char *ByteReader::get_data() { return data; }

uint64_t ByteReader::get_pos() { return pos; }
