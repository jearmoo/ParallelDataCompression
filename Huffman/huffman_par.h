/*
 * huffman_par.h
 *
 * author: tto@andrew.cmu.edu and cdaruwal@andrew.cmu.edu
 *
 */
#pragma once

#include "ByteReader.h"
#include "ByteWriter.h"

#include <vector>

using namespace std;

/* if alloc is true, allocates space for output */
void huffman_encode_par(ByteReader &input, ByteWriter &output, int num_threads,
                        bool sequential_freq, bool alloc);

/* offload to the xeon phi */
void huffman_encode_par_offload(ByteReader &input, ByteWriter &output,
                                int num_threads, bool sequential_freq);

/* if alloc is true, allocates space for output */
void huffman_decode_par(ByteReader &input, ByteWriter &output, int num_threads,
                        bool alloc);

/* offload to the xeon phi */
void huffman_decode_par_offload(ByteReader &input, ByteWriter &output,
                                int num_threads);