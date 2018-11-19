/*
 * ByteWriter.h
 *
 * author: tto@andrew.cmu.edu and cdaruwal@andrew.cmu.edu
 *
 */
#pragma once

#include <string>
#include <vector>

#include <stdint.h>

using namespace std;

class ByteWriter {
  private:
    unsigned char *data;
    uint64_t size;
    uint64_t pos;

  public:
    ByteWriter() : data(nullptr), size(0), pos(0){};
    ByteWriter(unsigned char *data, uint64_t size)
        : data(data), size(size), pos(0){};

    void set_size(uint64_t size);

    void assign(unsigned char *data, uint64_t size);

    /* write num_bytes bytes at pos, update pos */
    /* if manual_pos is specified, update there instead */
    void write(void *buf, uint64_t num_bytes, uint64_t manual_pos = -1L);

    /* append write a byte at pos, update pos */
    /* if manual_pos is specified, update there instead */
    void write_byte(unsigned char c, uint64_t manual_pos = -1);

    /* write everything in the byte store to the file file_name */
    void write_to_file(const string &file_name);
    void write_to_file_parallel(const string &file_name);

    /* change current position */
    void seek(uint64_t new_pos);

    /* reset current position to the beginning */
    void rewind();

    /* returns the number of bytes stored */
    size_t get_size();

    /* changes the size of the buffer */
    void resize(uint64_t size);

    /* get the position of the writer */
    uint64_t get_pos();

    /* get the underlying data */
    unsigned char *get_data();
};