/*
 * huffman_par.cpp
 *
 * author: tto@andrew.cmu.edu and cdaruwal@andrew.cmu.edu
 *
 * adapted from:
 *
 *  huffcode - Encode/Decode files using Huffman encoding.
 *  http://huffman.sourceforge.net
 *  Copyright (C) 2003  Douglas Ryan Richardson
 */

#include "huffman_par.h"

#include "mic.h"
#include "util.h"

#include <chrono>
#include <iostream>

#include <arpa/inet.h>
#include <assert.h>
#include <omp.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

typedef struct huffman_node_tag {
    unsigned char isLeaf;
    unsigned long count;
    struct huffman_node_tag *parent;

    union {
        struct {
            struct huffman_node_tag *zero, *one;
        };
        unsigned char symbol;
    };
} huffman_node;

typedef struct huffman_code_tag {
    /* The length of this code in bits. */
    unsigned long numbits;

    /* The bits that make up this code. The first
       bit is at position 0 in bits[0]. The second
       bit is at position 1 in bits[0]. The eighth
       bit is at position 7 in bits[0]. The ninth
       bit is at position 0 in bits[1]. */
    unsigned char *bits;
} huffman_code;

static unsigned long numbytes_from_numbits(unsigned long numbits) {
    return numbits / 8 + (numbits % 8 ? 1 : 0);
}

/*
 * get_bit returns the ith bit in the bits array
 * in the 0th position of the return value.
 */
static unsigned char get_bit(unsigned char *bits, unsigned long i) {
    return (bits[i / 8] >> i % 8) & 1;
}

static void reverse_bits(unsigned char *bits, unsigned long numbits) {
    unsigned long numbytes = numbytes_from_numbits(numbits);
    unsigned char *tmp = (unsigned char *)alloca(numbytes);
    unsigned long curbit;
    long curbyte = 0;

    memset(tmp, 0, numbytes);

    for (curbit = 0; curbit < numbits; ++curbit) {
        unsigned int bitpos = curbit % 8;

        if (curbit > 0 && curbit % 8 == 0)
            ++curbyte;

        tmp[curbyte] |= (get_bit(bits, numbits - curbit - 1) << bitpos);
    }

    memcpy(bits, tmp, numbytes);
}

#define MAX_SYMBOLS 256
typedef huffman_node *SymbolFrequencies[MAX_SYMBOLS];
typedef huffman_code *SymbolEncoder[MAX_SYMBOLS];

static huffman_node *new_leaf_node(unsigned char symbol) {
    huffman_node *p = (huffman_node *)malloc(sizeof(huffman_node));
    p->isLeaf = 1;
    p->symbol = symbol;
    p->count = 0;
    p->parent = 0;
    return p;
}

static huffman_node *new_nonleaf_node(unsigned long count, huffman_node *zero,
                                      huffman_node *one) {
    huffman_node *p = (huffman_node *)malloc(sizeof(huffman_node));
    p->isLeaf = 0;
    p->count = count;
    p->zero = zero;
    p->one = one;
    p->parent = 0;

    return p;
}

static void free_huffman_tree(huffman_node *subtree) {
    if (subtree == NULL)
        return;

    if (!subtree->isLeaf) {
        free_huffman_tree(subtree->zero);
        free_huffman_tree(subtree->one);
    }

    free(subtree);
}

static void init_frequencies(SymbolFrequencies *pSF) {
    memset(*pSF, 0, sizeof(SymbolFrequencies));
}

static void get_symbol_frequencies(SymbolFrequencies *pSF, ByteReader &input,
                                   vector<uint64_t> &chunk_frequency_counts,
                                   int num_threads, bool sequential) {
    /* Set all frequencies to 0. */
    init_frequencies(pSF);

    uint64_t chunk_size = get_chunksize(input.get_size(), num_threads);

    if (sequential) {

        /* Count the frequency of each symbol in the input file. */
        unsigned char uc;
        uint64_t counter = 0;
        while (input.read(&uc, 1)) {
            if (!(*pSF)[uc])
                (*pSF)[uc] = new_leaf_node(uc);
            ++(*pSF)[uc]->count;

            ++chunk_frequency_counts[MAX_SYMBOLS * (counter / chunk_size) + uc];

            ++counter;
        }
        input.rewind();
        return;
    }

    uint64_t char_chunk_size = get_chunksize(MAX_SYMBOLS, num_threads);

    vector<uint64_t> char_counts(MAX_SYMBOLS, 0);

    #pragma omp parallel
    {
        int thread_num = omp_get_thread_num();
        uint64_t start_index = chunk_size * thread_num;
        uint64_t end_index = min(start_index + chunk_size, input.get_size());

        uint64_t freq_counts_offset = MAX_SYMBOLS * thread_num;

        unsigned char uc;
        for (int i = start_index; i < end_index; ++i) {
            input.read(&uc, 1, i);

            ++chunk_frequency_counts[freq_counts_offset + uc];
        }

        #pragma omp barrier
        uint64_t start_char = thread_num * char_chunk_size;
        uint64_t end_char =
            min(start_char + char_chunk_size, (uint64_t)MAX_SYMBOLS);

        for (int i = start_char; i < end_char; ++i) {
            for (int j = 0; j < num_threads; ++j) {
                char_counts[i] += chunk_frequency_counts[MAX_SYMBOLS * j + i];
            }
        }
    }

    for (int i = 0; i < MAX_SYMBOLS; ++i) {
        if (char_counts[i] > 0) {
            (*pSF)[i] = new_leaf_node(i);
            (*pSF)[i]->count = char_counts[i];
        }
    }
}

/*
 * When used by qsort, SFComp sorts the array so that
 * the symbol with the lowest frequency is first. Any
 * NULL entries will be sorted to the end of the list.
 */
static int SFComp(const void *p1, const void *p2) {
    const huffman_node *hn1 = *(const huffman_node **)p1;
    const huffman_node *hn2 = *(const huffman_node **)p2;

    /* Sort all NULLs to the end. */
    if (hn1 == NULL && hn2 == NULL)
        return 0;
    if (hn1 == NULL)
        return 1;
    if (hn2 == NULL)
        return -1;

    if (hn1->count > hn2->count)
        return 1;
    else if (hn1->count < hn2->count)
        return -1;

    return 0;
}

/*
 * new_code builds a huffman_code from a leaf in
 * a Huffman tree.
 */
static huffman_code *new_code(const huffman_node *leaf) {
    /* Build the huffman code by walking up to
     * the root node and then reversing the bits,
     * since the Huffman code is calculated by
     * walking down the tree. */
    unsigned long numbits = 0;
    unsigned char *bits = NULL;
    huffman_code *p;

    while (leaf && leaf->parent) {
        huffman_node *parent = leaf->parent;
        unsigned char cur_bit = (unsigned char)(numbits % 8);
        unsigned long cur_byte = numbits / 8;

        /* If we need another byte to hold the code,
           then allocate it. */
        if (cur_bit == 0) {
            uint64_t newSize = cur_byte + 1;
            bits = (unsigned char *)realloc(bits, newSize);
            bits[newSize - 1] = 0; /* Initialize the new byte. */
        }

        /* If a one must be added then or it in. If a zero
         * must be added then do nothing, since the byte
         * was initialized to zero. */
        if (leaf == parent->one)
            bits[cur_byte] |= 1 << cur_bit;

        ++numbits;
        leaf = parent;
    }

    if (bits)
        reverse_bits(bits, numbits);

    p = (huffman_code *)malloc(sizeof(huffman_code));
    p->numbits = numbits;
    p->bits = bits;
    return p;
}

/*
 * build_symbol_encoder builds a SymbolEncoder by walking
 * down to the leaves of the Huffman tree and then,
 * for each leaf, determines its code.
 */
static void build_symbol_encoder(huffman_node *subtree, SymbolEncoder *pSF) {
    if (subtree == NULL)
        return;

    if (subtree->isLeaf)
        (*pSF)[subtree->symbol] = new_code(subtree);
    else {
        build_symbol_encoder(subtree->zero, pSF);
        build_symbol_encoder(subtree->one, pSF);
    }
}

/*
 * calculate_huffman_codes turns pSF into an array
 * with a single entry that is the root of the
 * huffman tree. The return value is a SymbolEncoder,
 * which is an array of huffman codes index by symbol value.
 */
static SymbolEncoder *calculate_huffman_codes(SymbolFrequencies *pSF) {
    unsigned int i = 0;
    unsigned int n = 0;
    huffman_node *m1 = NULL, *m2 = NULL;
    SymbolEncoder *pSE = NULL;

    /* Sort the symbol frequency array by ascending frequency. */
    qsort((*pSF), MAX_SYMBOLS, sizeof((*pSF)[0]), SFComp);

    /* Get the number of symbols. */
    for (n = 0; n < MAX_SYMBOLS && (*pSF)[n]; ++n)
        ;

    /*
     * Construct a Huffman tree. This code is based
     * on the algorithm given in Managing Gigabytes
     * by Ian Witten et al, 2nd edition, page 34.
     * Note that this implementation uses a simple
     * count instead of probability.
     */
    for (i = 0; i < n - 1; ++i) {
        /* Set m1 and m2 to the two subsets of least probability. */
        m1 = (*pSF)[0];
        m2 = (*pSF)[1];

        /* Replace m1 and m2 with a set {m1, m2} whose probability
         * is the sum of that of m1 and m2. */
        (*pSF)[0] = m1->parent = m2->parent =
            new_nonleaf_node(m1->count + m2->count, m1, m2);
        (*pSF)[1] = NULL;

        /* Put newSet into the correct count position in pSF. */
        qsort((*pSF), n - i, sizeof((*pSF)[0]), SFComp);
    }

    /* Build the SymbolEncoder array from the tree. */
    pSE = (SymbolEncoder *)malloc(sizeof(SymbolEncoder));
    memset(pSE, 0, sizeof(SymbolEncoder));
    build_symbol_encoder((*pSF)[0], pSE);
    return pSE;
}

static void free_code(huffman_code *p) {
    free(p->bits);
    free(p);
}

static void free_encoder(SymbolEncoder *pSE) {
    unsigned long i;
    for (i = 0; i < MAX_SYMBOLS; ++i) {
        huffman_code *p = (*pSE)[i];
        if (p)
            free_code(p);
    }

    free(pSE);
}

/*
 * Write the huffman code table. The format is:
 * 4 byte code count in network byte order.
 * 4 byte number of bytes encoded per section
 *   (if you decode the data, you should get this number of bytes)
 * 4 byte number of processors it was encoded with
 * code1
 * ...
 * codeN, where N is the count read at the begginning of the file.
 * Each codeI has the following format:
 * 1 byte symbol, 1 byte code bit length, code bytes.
 * Each entry has numbytes_from_numbits code bytes.
 * The last byte of each code may have extra bits, if the number of
 * bits in the code is not a multiple of 8.
 */
static int write_code_table(ByteWriter &output, SymbolEncoder *se,
                            uint64_t symbol_count, int num_threads) {
    uint32_t i, count = 0;

    /* Determine the number of entries in se. */
    for (i = 0; i < MAX_SYMBOLS; ++i) {
        if ((*se)[i])
            ++count;
    }

    /* Write the number of entries in network byte order. */
    i = htonl(count);
    output.write(&i, sizeof(uint32_t));

    /* Write the number of bytes that will be encoded per chunk. */
    symbol_count = htonl(symbol_count);
    output.write(&symbol_count, sizeof(uint32_t));

    /* Write the entries. */
    for (i = 0; i < MAX_SYMBOLS; ++i) {
        huffman_code *p = (*se)[i];
        if (p) {
            unsigned int numbytes;
            /* Write the 1 byte symbol. */
            output.write_byte(i);
            /* Write the 1 byte code bit length. */
            output.write_byte(p->numbits);
            /* Write the code bytes. */
            numbytes = numbytes_from_numbits(p->numbits);
            output.write(p->bits, numbytes);
        }
    }

    return 0;
}

static int do_file_encode(ByteReader &input, ByteWriter &output,
                          SymbolEncoder *se, int num_threads,
                          vector<uint64_t> &chunk_offsets) {
    unsigned char *data = input.get_data();
    uint64_t chunk_size = get_chunksize(input.get_size(), num_threads);

    uint64_t output_offset = output.get_pos();
    #pragma omp parallel
    {
        int num_threads = omp_get_num_threads();
        int thread_num = omp_get_thread_num();
        uint64_t thread_start = chunk_size * thread_num;
        uint64_t thread_end = min(thread_start + chunk_size, input.get_size());
        unsigned char curbyte = 0;
        unsigned char curbit = 0;

        output.write(&chunk_offsets[thread_num], sizeof(uint64_t),
                     output_offset + sizeof(uint64_t) * thread_num);

        uint64_t chunk_offset =
            output_offset + 8 * num_threads + chunk_offsets[thread_num];

        uint64_t index = thread_start;

        uint64_t num_bytes_written = 0;

        while (index < thread_end) {
            huffman_code *code = (*se)[data[index]];
            unsigned long i;

            for (i = 0; i < code->numbits; ++i) {
                /* Add the current bit to curbyte. */
                curbyte |= get_bit(code->bits, i) << curbit;

                /* If this byte is filled up then write it
                 * out and reset the curbit and curbyte. */
                if (++curbit == 8) {
                    output.write_byte(curbyte,
                                      chunk_offset + num_bytes_written);
                    curbyte = 0;
                    curbit = 0;
                    ++num_bytes_written;
                }
            }

            index++;
        }

        /*
         * If there is data in curbyte that has not been
         * output yet, which means that the last encoded
         * character did not fall on a byte boundary,
         * then output it.
         */
        if (curbit > 0) {
            output.write_byte(curbyte, chunk_offset + num_bytes_written);
            ++num_bytes_written;
        }

        D(if (thread_num != num_threads - 1) {
            assert(num_bytes_written ==
                   chunk_offsets[thread_num + 1] - chunk_offsets[thread_num]);
        } else {
            assert(output.get_size() == chunk_offset + num_bytes_written);
        })
    }

    /*uint64_t pos = 0;
    for (int i = 0; i < num_threads; ++i) {
        // append the starting location of each section
        output[num_threads].append(&pos, sizeof(uint64_t));
        pos += output[i].get_size();
    }*/

    return 0;
}

/*
 * read_code_table builds a Huffman tree from the code
 * in the in file. This function returns NULL on error.
 * The returned value should be freed with free_huffman_tree.
 */
static huffman_node *read_code_table(ByteReader &input, uint32_t *pDataBytes,
                                     vector<uint64_t> &chunk_offsets,
                                     int num_threads) {
    huffman_node *root = new_nonleaf_node(0, NULL, NULL);
    uint32_t count;

    /* Read the number of entries.
       (it is stored in network byte order). */
    if (!input.read(&count, sizeof(uint32_t))) {
        free_huffman_tree(root);
        return NULL;
    }

    count = ntohl(count);

    /* Read the number of data bytes this encoding represents. */
    if (!input.read(pDataBytes, sizeof(uint32_t))) {
        free_huffman_tree(root);
        return NULL;
    }

    *pDataBytes = ntohl(*pDataBytes);

    /* Read the entries. */
    while (count-- > 0) {
        unsigned int curbit;
        unsigned char symbol;
        unsigned char numbits;
        unsigned char numbytes;
        unsigned char *bytes;
        huffman_node *p = root;

        if (!input.read(&symbol, 1)) {
            free_huffman_tree(root);
            return NULL;
        }

        if (!input.read(&numbits, 1)) {
            free_huffman_tree(root);
            return NULL;
        }

        numbytes = (unsigned char)numbytes_from_numbits(numbits);
        bytes = (unsigned char *)malloc(numbytes);
        if (!(input.read(bytes, numbytes))) {
            free(bytes);
            free_huffman_tree(root);
            return NULL;
        }

        /*
         * Add the entry to the Huffman tree. The value
         * of the current bit is used switch between
         * zero and one child nodes in the tree. New nodes
         * are added as needed in the tree.
         */
        for (curbit = 0; curbit < numbits; ++curbit) {
            if (get_bit(bytes, curbit)) {
                if (p->one == NULL) {
                    p->one = curbit == (unsigned char)(numbits - 1)
                                 ? new_leaf_node(symbol)
                                 : new_nonleaf_node(0, NULL, NULL);
                    p->one->parent = p;
                }
                p = p->one;
            } else {
                if (p->zero == NULL) {
                    p->zero = curbit == (unsigned char)(numbits - 1)
                                  ? new_leaf_node(symbol)
                                  : new_nonleaf_node(0, NULL, NULL);
                    p->zero->parent = p;
                }
                p = p->zero;
            }
        }

        free(bytes);
    }

    for (int i = 0; i < num_threads; ++i) {
        uint64_t offset;
        input.read(&offset, sizeof(uint64_t));
        chunk_offsets[i] = offset;
    }

    return root;
}

uint64_t get_output_size(ByteReader &input, SymbolEncoder *se, int num_threads,
                         vector<uint64_t> &chunk_frequency_counts,
                         vector<uint64_t> &chunk_offsets) {
    // 4 bytes for # entries in se and 4 bytes for # bytes that will be encoded
    // numthreads * 8 bytes to encode the offsets for decompression
    uint64_t res = 8 + 8 * num_threads;

    vector<uint64_t> thread_byte_counts(num_threads);

    // find the encoded length of what each thread encodes
    #pragma omp parallel
    {
        int thread_num = omp_get_thread_num();
        uint64_t bits_counter = 0;

        for (int i = 0; i < MAX_SYMBOLS; ++i) {
            huffman_code *p = (*se)[i];
            if (p) {
                /* number of bytes used to encode in code table */
                bits_counter +=
                    chunk_frequency_counts[thread_num * MAX_SYMBOLS + i] *
                    p->numbits;

                // the root thread will account for the size of the code table
                if (thread_num == 0) {
                    // 1 byte symbol and 1 byte code bit length
                    res += 2;

                    /* number of bytes used to encode in code table */
                    res += numbytes_from_numbits(p->numbits);
                }
            }
        }

        thread_byte_counts[thread_num] = numbytes_from_numbits(bits_counter);
    }

    chunk_offsets[0] = 0;
    for (int i = 0; i < num_threads; ++i) {
        res += thread_byte_counts[i];
        if (i != num_threads - 1) {
            chunk_offsets[i + 1] = chunk_offsets[i] + thread_byte_counts[i];
        }
    }

    return res;
}

void huffman_encode_par(ByteReader &input, ByteWriter &output, int num_threads,
                        bool sequential_freq, bool alloc) {
    if (alloc) {
        output.resize(input.get_size());
    }

    cout << "max threads is " << omp_get_max_threads() << endl;
    omp_set_num_threads(min(num_threads, omp_get_max_threads()));
    SymbolFrequencies sf;
    SymbolEncoder *se;
    huffman_node *root = NULL;

    uint64_t symbol_count = input.get_size();

    auto start_time = Clock::now();
    double init_time;
    double prev_time = 0;

    vector<uint64_t> chunk_frequency_counts(num_threads * MAX_SYMBOLS, 0);

    /* Get the frequency of each symbol in the input file. */
    // PARALLEL
    get_symbol_frequencies(&sf, input, chunk_frequency_counts, num_threads,
                           sequential_freq);
    init_time = duration_cast<dsec>(Clock::now() - start_time).count();

    cout << "Frequency time:\t" << init_time - prev_time << endl;
    cout << "Frequency time (cum):\t" << init_time << endl;

    prev_time = init_time;

    /* Build an optimal table from the symbolCount. */
    // SEQUENTIAL
    se = calculate_huffman_codes(&sf);
    root = sf[0];

    init_time = duration_cast<dsec>(Clock::now() - start_time).count();

    cout << "Build Tree time:\t" << init_time - prev_time << endl;
    cout << "Build Tree time (cum):\t" << init_time << endl;

    prev_time = init_time;

    /* Scan the file again and, using the table
    previously built, encode it into the output file. */
    // SEQUENTIAL
    write_code_table(output, se, symbol_count, num_threads);
    init_time = duration_cast<dsec>(Clock::now() - start_time).count();

    cout << "Write Code Table Time:\t" << init_time - prev_time << endl;
    cout << "Write Code Table time (cum):\t" << init_time << endl;

    prev_time = init_time;

    vector<uint64_t> chunk_offsets(num_threads);
    /* Calculate the output size of the encoded file */
    // PARALLEL
    uint64_t output_size = get_output_size(
        input, se, num_threads, chunk_frequency_counts, chunk_offsets);

    output.set_size(output_size);

    // PARALLEL
    do_file_encode(input, output, se, num_threads, chunk_offsets);
    init_time = duration_cast<dsec>(Clock::now() - start_time).count();

    cout << "Encode time:\t" << init_time - prev_time << endl;
    cout << "Encode time (cum total):\t" << init_time << endl;

    prev_time = init_time;

    /* Free the Huffman tree. */
    free_huffman_tree(root);
    free_encoder(se);

    cout << "Output size is " << output_size << endl;
}

void huffman_encode_par_offload(ByteReader &input, ByteWriter &output,
                                int num_threads, bool sequential_freq) {
    unsigned char *in_data = input.get_data();
    uint64_t input_size = input.get_size();

    unsigned char *out_data = new unsigned char[input_size];

    uint64_t output_size;

    #ifdef RUN_MIC
    #pragma offload target(mic) \
        in(in_data: length(input_size) INOUT) \
        out(out_data: length(input_size) INOUT) \
        out(output_size)
    #endif
    {
        ByteReader mic_reader(in_data, input_size);
        ByteWriter mic_writer(out_data, input_size);
        huffman_encode_par(mic_reader, mic_writer, num_threads, sequential_freq,
                           false);
        output_size = mic_writer.get_size();
    }

    output.assign(out_data, output_size);
}

void huffman_decode_par(ByteReader &input, ByteWriter &output, int num_threads,
                        bool alloc) {

    huffman_node *root;
    int c;
    unsigned int data_count;

    cout << "max threads is " << omp_get_max_threads() << endl;
    omp_set_num_threads(min(num_threads, omp_get_max_threads()));

    vector<uint64_t> chunk_offsets(num_threads);

    auto start_time = Clock::now();

    /* Read the Huffman code table. */
    root = read_code_table(input, &data_count, chunk_offsets, num_threads);

    if (alloc) {
        cout << "Output size is " << data_count << endl;
        output.resize(data_count);
    }

    uint64_t chunk_size = get_chunksize(data_count, num_threads);

    uint64_t input_offset = input.get_pos();

    /* Decode the file. */
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        uint64_t start_index = chunk_size * tid;
        uint64_t end_index =
            min(start_index + chunk_size, (uint64_t)data_count);

        uint64_t chunk_offset = chunk_offsets[tid];

        uint64_t thread_data_count = end_index - start_index;
        uint64_t read_bytes = 0;
        uint64_t written_bytes = 0;

        huffman_node *p = root;

        unsigned char byte;
        while (written_bytes < thread_data_count) {
            assert(
                input.read(&byte, 1, input_offset + chunk_offset + read_bytes));
            ++read_bytes;

            unsigned char mask = 1;
            while (written_bytes < thread_data_count && mask) {
                p = byte & mask ? p->one : p->zero;
                mask <<= 1;

                if (p->isLeaf) {
                    output.write_byte(p->symbol, start_index + written_bytes);
                    p = root;
                    ++written_bytes;
                }
            }
        }
    }

    cout << "Decompression Time: "
         << duration_cast<dsec>(Clock::now() - start_time).count() << endl;

    free_huffman_tree(root);

} /* Set all frequencies to 0. */

void huffman_decode_par_offload(ByteReader &input, ByteWriter &output,
                                int num_threads) {

    unsigned char *in_data = input.get_data();
    uint64_t input_size = input.get_size();

    uint32_t num_bytes;
    input.read(&num_bytes, sizeof(uint32_t), 4);

    num_bytes = ntohl(num_bytes);
    cout << "Output size is " << num_bytes << endl;

    output.resize(num_bytes);
    unsigned char *out_data = output.get_data();

    #ifdef RUN_MIC
    #pragma offload target(mic) \
        in(in_data: length(input_size) INOUT) \
        out(out_data: length(num_bytes) INOUT)
    #endif
    {
        ByteReader mic_reader(in_data, input_size);
        ByteWriter mic_writer(out_data, num_bytes);
        huffman_decode_par(mic_reader, mic_writer, num_threads, false);
    }
}
