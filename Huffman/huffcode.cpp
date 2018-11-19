/*
 * huffcode.cpp
 *
 * author: tto@andrew.cmu.edu and cdaruwal@andrew.cmu.edu
 *
 */
#include "ByteReader.h"
#include "ByteWriter.h"
#include "huffman_par.h"
#include "huffman_seq.h"
#include "util.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <omp.h>
#include <stdio.h>
#include <unistd.h>

using namespace std;

static void usage(FILE *out) {
    // Sample usage to run benchmarking. ./huffmancode -b -i input_file
    fputs(
        "Usage:\n"
        "  huffcode -c -i <input file> -o <output file> [-t <num threads>]\n"
        "  huffcode -d -i <input file> -o <output file> [-t <num threads>]\n\n"
        "  -h - print usage information\n"
        "  -i - input file\n"
        "  -o - output file\n"
        "  -c - compress\n"
        "  -d - decompress\n"
        "  -t - number of threads (defaults to 1)\n",
        out);
}

int main(int argc, char **argv) {

    char *infile, *outfile;
    int c;

    bool compress_flag = false;
    bool decompress_flag = false;
    bool sequential_freq = false;

    int num_threads = 1;

    while ((c = getopt(argc, argv, "cdhi:o:t:f")) != -1)
        switch (c) {
        case 'f':
            sequential_freq = true;
            break;
        case 't':
            num_threads = atoi(optarg);
            break;
        case 'c':
            compress_flag = true;
            break;
        case 'd':
            decompress_flag = true;
            break;
        case 'h':
            usage(stdout);
            return 0;
        case 'i':
            infile = optarg;
            break;
        case 'o':
            outfile = optarg;
            break;
        default:
            return 1;
        }

    if (infile == nullptr || outfile == nullptr) {
        cerr << "AN INPUT AND OUTPUT FILE IS REQUIRED" << endl;
        return 1;
    }

    cout << "--------------------------------------------" << endl;
    cout << "Using " << num_threads << " threads" << endl;

    // This sets up things for offloading to the MIC
    #ifdef RUN_MIC
    #pragma offload_transfer target(mic)
    #endif

    ByteReader input;
    if (!input.load_from_file_parallel(infile, num_threads)) {
        cerr << "UNABLE TO READ THE FILE " << infile << endl;
        return 1;
    }

    ByteWriter output;

    // Run on MIC if compiled to do so
    if (compress_flag) {
        #ifdef RUN_MIC
        // Running both sequential and parallel on the xeon phi to test speedup
        // vs num_threads
        huffman_encode_par_offload(input, output, num_threads, sequential_freq);
        #else
        huffman_encode_par(input, output, num_threads, sequential_freq, true);
        #endif

        cout << "Compression ratio is "
             << (double)output.get_size() / input.get_size() << endl;
    } else {
        #ifdef RUN_MIC
        // Running both sequential and parallel on the xeon phi to test speedup
        // vs num_threads
        huffman_decode_par_offload(input, output, num_threads);
        #else
        huffman_decode_par(input, output, num_threads, true);
        #endif
    }

    output.write_to_file(outfile);

    return 0;
}