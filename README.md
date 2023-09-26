# Cache-Controller-Design
The goal of this project is to determine the best architecture for a cache controller that interfaces to a 32- bit microprocessor with a 32-bit data bus. While the microprocessor 
is general purpose in design and performs a number of functions, it is desired to speed up certain signal processing functions, such as a fast Fourier transform (FFT) routine.
The project allows 4GiB of SDRAM to be interfaced using a 32-bit data bus (arranged as 230 x 32-bits) and the cache should be limited to 256 KiB in size. The size of the FFT will be 32768 points (512KiB) in normal operation.
The metric being used to evaluate the best performance is the total access time required to run the Radix2FFT(), including the time required to flush the cache after the FFT is computed.
Project describes the set size (N-way), number of cache lines (L), block size (B), write strategy (write-back allocate, write-through allocate, and write-through non-allocate), and cache line block replacement strategy. (round robin or LRU)
