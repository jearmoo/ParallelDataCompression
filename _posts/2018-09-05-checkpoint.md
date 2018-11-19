---
title: "checkpoint"
bg: tuscan  #defined in _config.yml, can use html color like '#0fbfcf'
color: white   #text color
style: left
---

## Summary of What We Completed

We have written and optimized the sequential version of the Huffman encoding and
decoding algorithms, and tested it. For the parallel CPU version of this, we
were debating between SIMD intrinsics and ISPC, and OpenMP.

However, huffman coding compression and decompression doesn’t seem to have a
workload that can appropriately use SIMD. This is because there is no elegant
way of dealing with bits instead of bytes in SIMD. Moreover, different bytes
compress to a different number of bits (there is no fixed mapping of input
vector size to output vector size), which makes byte alignment in SIMD very
difficult (for example, the compressed form for a random 4 byte input could
range from 2 to 4 bytes). This is much worse for decompression, where resolving
bit-level conflicts (where a specific encoding spreads over 2 bytes) is almost
impossible and might actually result in the algorithm being slower than the
sequential version. Therefore, we decided to focus on OpenMP.

For compression, we first sort the array in parallel, to minimize number of
concurrent updates to the shared frequency dictionary, reducing contention and
false sharing. We also use fine-grained locking for the frequency dictionary,
individually locking each key-value pair. Once the symbol codes have been
determined, each symbol is replaced by its code, and all symbols are so
processed in parallel.

Decompression is inherently sequential, and hence much harder to parallelize.
In this case, we take advantage of the self-synchronizing property of Huffman
coding, which allows us to start at an arbitrary point in the encoded bits, and
assume that at some point, the offset in bits will correct itself, resulting in
the correct output thereafter.

We read about the LZ77 algorithm and explored the different variants of the
algorithm. We also explored different ways to parallelize LZ77. One naive
approach is running the LZ77 algorithm along different segments of the data.
This approach could output the same result as the sequential implementation if
we use a fixed size sliding window and reread over some of the data. Another
approach is the one outlined in Practical Parallel Lempel-Ziv Factorization
which uses an unbounded sliding window and employs the use of prefix sums and
segment trees to calculate the Lempel-Ziv factorization in parallel.

## Update on Deliverables

Our sequential implementations are close to finished, and we have some idea of how to parallelize the algorithms. Our goal for the checkpoint was to have both of these parts finished, but we have not completely met the goal. We may pivot and work on parallelizing the compression and decompression of the huffman coding algorithm and drop the LZ77 part of the project altogether.

**Our new goals:**
1. Parallelize the Huffman Coding compression.
2. Parallelize the Huffman Coding decompression or LZ77 compression

**Hope to achieve:**
1. Both parts of part 2 in our new goals.




## What We Will Bring to the Poster Session

At the poster session we plan to present graphs that show the speedup vs number
of processors and compression ratio vs number of processors for the Huffman
Encoding and LZ77 algorithms. We will also have a demo of our compression
algorithms with some sample data files to use them on.

## Issues

We are stuck at deciding what our encoding for LZ77 will be. There are so many
different implementations of LZ77 and encoding schemes for the algorithm.
We have not yet settled on one encoding scheme. *Practical Parallel Lempel-Ziv
Factorization* doesn't specify what kind of encoding scheme to use. The fact
that they use unbounded windows complicates things.