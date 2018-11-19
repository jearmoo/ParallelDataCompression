/*
 * huffman_seq.cpp
 *
 * author: tto@andrew.cmu.edu and cdaruwal@andrew.cmu.edu
 *
 * adapted from:
 *
 *  huffcode - Encode/Decode files using Huffman encoding.
 *  http://huffman.sourceforge.net
 *  Copyright (C) 2003  Douglas Ryan Richardson
 */
#include "huffman_seq.h"

#include "util.h"

#include <chrono>
#include <iostream>

#include <arpa/inet.h>
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
    /* number of times code appears*/
    unsigned long count;

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

static void get_symbol_frequencies(SymbolFrequencies *pSF, ByteReader &input) {
    /* Set all frequencies to 0. */
    init_frequencies(pSF);

    /* Count the frequency of each symbol in the input file. */
    unsigned char uc;
    while (input.read(&uc, 1)) {
        if (!(*pSF)[uc])
            (*pSF)[uc] = new_leaf_node(uc);
        ++(*pSF)[uc]->count;
    }
    input.rewind();
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

    int count = leaf->count;

    while (leaf && leaf->parent) {
        huffman_node *parent = leaf->parent;
        unsigned char cur_bit = (unsigned char)(numbits % 8);
        unsigned long cur_byte = numbits / 8;

        /* If we need another byte to hold the code,
           then allocate it. */
        if (cur_bit == 0) {
            size_t newSize = cur_byte + 1;
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
    p->count = count;
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
 * 4 byte number of bytes encoded
 *   (if you decode the data, you should get this number of bytes)
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
                            uint32_t symbol_count) {
    uint32_t i, count = 0;

    /* Determine the number of entries in se. */
    for (i = 0; i < MAX_SYMBOLS; ++i) {
        if ((*se)[i])
            ++count;
    }

    /* Write the number of entries in network byte order. */
    i = htonl(count);
    output.write(&i, sizeof(uint32_t));

    /* Write the number of bytes that will be encoded. */
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
                          SymbolEncoder *se) {
    unsigned char curbyte = 0;
    unsigned char curbit = 0;
    char c;

    while (input.read(&c, 1)) {
        unsigned char uc = (unsigned char)c;
        huffman_code *code = (*se)[uc];
        unsigned long i;

        for (i = 0; i < code->numbits; ++i) {
            /* Add the current bit to curbyte. */
            curbyte |= get_bit(code->bits, i) << curbit;

            /* If this byte is filled up then write it
             * out and reset the curbit and curbyte. */
            if (++curbit == 8) {
                output.write_byte(curbyte);
                curbyte = 0;
                curbit = 0;
            }
        }
    }
    input.rewind();

    /*
     * If there is data in curbyte that has not been
     * output yet, which means that the last encoded
     * character did not fall on a byte boundary,
     * then output it.
     */
    if (curbit > 0)
        output.write_byte(curbyte);

    return 0;
}

/*
 * read_code_table builds a Huffman tree from the code
 * in the in file. This function returns NULL on error.
 * The returned value should be freed with free_huffman_tree.
 */
static huffman_node *read_code_table(ByteReader &input, uint32_t *pDataBytes) {
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

    return root;
}

uint64_t get_output_size(ByteReader &input, SymbolEncoder *se) {
    // 4 bytes for # entries in se and 4 bytes for # bytes that will be encoded
    uint64_t res = 8;

    uint64_t bits_counter = 0;
    for (int i = 0; i < MAX_SYMBOLS; ++i) {
        huffman_code *p = (*se)[i];
        if (p) {
            // 1 byte symbol and 1 byte code bit length
            res += 2;

            /* number of bytes used to encode in code table */
            res += numbytes_from_numbits(p->numbits);

            /* account for number of bits in the actual content */
            bits_counter += p->count * p->numbits;
        }
    }

    res += numbytes_from_numbits(bits_counter);

    return res;
}

void huffman_encode_seq(ByteReader &input, ByteWriter &output) {
    SymbolFrequencies sf;
    SymbolEncoder *se;
    huffman_node *root = NULL;
    uint64_t symbol_count = input.get_size();

    output.resize(input.get_size());

    auto start_time = Clock::now();
    double init_time;
    double prev_time = 0;
    /* Get the frequency of each symbol in the input file. */
    get_symbol_frequencies(&sf, input);

    init_time = duration_cast<dsec>(Clock::now() - start_time).count();

    cout << "Frequency time:\t" << init_time - prev_time << endl;
    cout << "Frequency time (cum):\t" << init_time << endl;

    prev_time = init_time;

    /* Build an optimal table from the symbolCount. */
    se = calculate_huffman_codes(&sf);
    root = sf[0];

    init_time = duration_cast<dsec>(Clock::now() - start_time).count();

    cout << "Build Tree time:\t" << init_time - prev_time << endl;
    cout << "Build Tree time (cum):\t" << init_time << endl;

    prev_time = init_time;

    /* Scan the file again and, using the table
       previously built, encode it into the output file. */
    write_code_table(output, se, symbol_count);

    init_time = duration_cast<dsec>(Clock::now() - start_time).count();

    cout << "Write Code Table Time:\t" << init_time - prev_time << endl;
    cout << "Write Code Table time (cum):\t" << init_time << endl;

    prev_time = init_time;

    /* Calculate the output size of the encoded file */
    uint64_t output_size = get_output_size(input, se);

    output.set_size(output_size);

    do_file_encode(input, output, se);

    init_time = duration_cast<dsec>(Clock::now() - start_time).count();

    cout << "Encode time:\t" << init_time - prev_time << endl;
    cout << "Encode time (cum total):\t" << init_time << endl;

    prev_time = init_time;

    /* Free the Huffman tree. */
    free_huffman_tree(root);
    free_encoder(se);
}

void huffman_decode_seq(ByteReader &input, ByteWriter &output) {
    huffman_node *root, *p;
    int c;
    unsigned int data_count;

    /* Read the Huffman code table. */
    root = read_code_table(input, &data_count);

    cout << "Output size is " << data_count << endl;
    output.resize(data_count);

    /* Decode the file. */
    p = root;
    while (data_count > 0 && (input.read(&c, 1))) {
        unsigned char byte = (unsigned char)c;
        unsigned char mask = 1;
        while (data_count > 0 && mask) {
            p = byte & mask ? p->one : p->zero;
            mask <<= 1;

            if (p->isLeaf) {
                output.write_byte(p->symbol);
                p = root;
                --data_count;
            }
        }
    }

    free_huffman_tree(root);
}