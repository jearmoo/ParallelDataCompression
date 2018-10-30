---
title: "background"
bg: palm-leaf  #defined in _config.yml, can use html color like '#0fbfcf'
color: white   #text color
fa-icon: pied-piper-alt
---

So far, the main compression algorithm that we have looked at is Huffman Coding

The parallelizable components are as follows
1. Reading from input file to make a frequency dictionary
2. Building a Huffman tree.
3. Traversing the Huffman Tree and build the prefix code table.
4. Encoding the input file using the prefix code table (maximum potential for parallelism since it is the slowest sequentially)