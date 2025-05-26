# rsexh

Reed-Solomon & Extended Hamming codec.

An arbitrary linear block code other than the Extended Hamming code can also be set.

## Description

Two codes are used connected in cascade.

The first code is a binary linear block code (n, k, d) with the code distance d, where k - the number of input symbols. The input symbol is an arbitrary vector of integers (your choice). Binary code is selected due to simplicity, but one symbol is a vector of integers. The code works in erasures recovery mode for the simplicity and efficiency. The time complexity for the full decode procedure is O(nr) + O(e^2 * r) = O(nr), where e - the number of erasured symbols.

The second code is a Reed-Solomon (RS) code in the field GF(2^4) with an arbitrary code distance (your choice). The code length is n = 15. Such choice was made to use LUT for 1- and 2-error corrections for better performance. An RS codeword is the one symbol for the linear block code (the first code). If RS-code cannot correct the error (but detect it), the output symbol will marked as erased. So, two codes complement each other if you tune the codec properly. Two cases are shown in main.cpp file.

For 1-error corrections full LUT is used, but for 2-error corrections the followed trick is applied to reduce memory consumption.

We believe a 2-error has the form (x, 0, 0, ..., x, ...), i.e. the first broken symbol has always zero index. LUT is filled for all such errors. When 2-errored codeword has been received, calculated cyndrome is looked up in the LUT. If not found, the cyndrome is recalculated for left-shifted (by one symbol) error vector, etc. No more then (n-1) checks is required for the described procedure, where n - the code length.

The cyndrome calculation has O(nr) time complexity. Total time complexity for 2-error correction is O(n) + O(nr) = O(nr). For LUT implementation a hash table is used, therefore for 1-error correction the total time complexity is O(nr) + O(1) = O(nr). So, error correction procedure has O(nr) time complexity. The space complexity is O(r * n^3) due to LUT usage.

The decoding procedure (after error correction) has O(n^2) complexity, but can be translated to O(n log n) due to Fast Fourier Transformation property. The code length n = 15 is small enough, so no need to make the decoding algorithm more complex: it can be implemented in hardware (for example: GPU, FPGA, ASIC) in parallel.
