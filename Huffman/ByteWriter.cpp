/*
 * ByteWriter.cpp
 *
 * author: tto@andrew.cmu.edu and cdaruwal@andrew.cmu.edu
 *
 */
#include "ByteWriter.h"

#include "util.h"

#include <fstream>
#include <string>
#include <vector>

#include <assert.h>
#include <omp.h>
#include <stdint.h>
#include <stdlib.h>

void ByteWriter::write(void *buf, uint64_t num_bytes, uint64_t manual_pos) {
    if (manual_pos != -1L) {
        D(assert(0 <= manual_pos && manual_pos + num_bytes <= size);)
        for (uint64_t i = 0; i < num_bytes; ++i) {
            data[manual_pos + i] = static_cast<unsigned char *>(buf)[i];
        }
        return;
    }
    D(assert(pos + num_bytes <= size);)
    for (uint64_t i = 0; i < num_bytes; ++i) {
        data[pos++] = static_cast<unsigned char *>(buf)[i];
    }
}

void ByteWriter::write_byte(unsigned char c, uint64_t manual_pos) {
    if (manual_pos != -1L) {
        D(assert(0 <= manual_pos && manual_pos < size);)
        data[manual_pos] = c;
        return;
    }
    D(assert(pos + 1 <= size);)
    data[pos++] = c;
}

void ByteWriter::write_to_file(const string &file_name) {
    ofstream outfile_stream(file_name, ios::out | ios::binary);
    if (outfile_stream.is_open()) {
        outfile_stream.write((char *)data, size);
    }
}

void ByteWriter::write_to_file_parallel(const string &file_name) {
#pragma omp parallel
    {
        int num_threads = omp_get_num_threads();
        uint64_t chunk_size = get_chunksize(size, num_threads);
        int thread_num = omp_get_thread_num();

        ofstream thread_outfile(file_name, ios::out | ios::binary);

        uint64_t thread_start_pos = chunk_size * thread_num;

        if (thread_outfile.is_open()) {
            thread_outfile.seekp(thread_start_pos, ios::beg);
            thread_outfile.write((char *)data + thread_start_pos,
                                 min((thread_num + 1) * chunk_size, size) -
                                     thread_start_pos);
            thread_outfile.close();
        }
    }
}

void ByteWriter::seek(uint64_t new_pos) {
    D(assert(0 <= new_pos && new_pos <= size);)
    pos = new_pos;
}

void ByteWriter::rewind() { pos = 0; }

uint64_t ByteWriter::get_size() { return size; }

void ByteWriter::resize(uint64_t size) {
    delete[] data;
    this->size = size;
    data = new unsigned char[size];
}

uint64_t ByteWriter::get_pos() { return pos; }

void ByteWriter::assign(unsigned char *data, uint64_t size) {
    this->data = data;
    this->size = size;
}

void ByteWriter::set_size(uint64_t size) { this->size = size; }

unsigned char *ByteWriter::get_data() { return data; }
