/*
 * huffman_seq.h
 *
 * author: tto@andrew.cmu.edu and cdaruwal@andrew.cmu.edu
 *
 */
#pragma once

#include "ByteReader.h"
#include "ByteWriter.h"

#include <vector>

using namespace std;

/* true if succeeds, false if fails */
void huffman_encode_seq(ByteReader &input, ByteWriter &output);

/* true if succeeds, false if fails */
void huffman_decode_seq(ByteReader &input, ByteWriter &output);