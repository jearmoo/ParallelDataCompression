/*
 * ByteReader.h
 *
 * author: tto@andrew.cmu.edu and cdaruwal@andrew.cmu.edu
 *
 */
#pragma once

#include <string>
#include <vector>

#include <stdint.h>

using namespace std;

class ByteReader {
  private:
    unsigned char *data;
    uint64_t size;
    uint64_t pos;

  public:
    ByteReader() : data(nullptr), size(0), pos(0){};
    ByteReader(unsigned char *data, uint64_t size)
        : data(data), size(size), pos(0){};

    /* initialize the bytereader with the content of a file */
    /* returns true on success, false on failure*/
    bool load_from_file(const string &file_name);
    bool load_from_file_parallel(const string &file_name, int num_threads);

    /* load buf with num_bytes from the data at the current position pos */
    /* returns true on success, false on failure */
    bool read(void *buf, uint64_t num_bytes, uint64_t manual_pos = -1L);

    /* change current position */
    void seek(uint64_t new_pos);

    /* reset current position to the beginning */
    void rewind();

    /* returns true if at the end of the data, false otherwise */
    bool at_end_of_data();

    /* returns the number of bytes stored */
    uint64_t get_size();

    /* get a pointer to the underlying datastructure holding the data */
    unsigned char *get_data();

    /* returns the position of the reader */
    uint64_t get_pos();
};