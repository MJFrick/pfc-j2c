Accelerated pipeline for quick prototyping PFC code by writing in Julia and generating an equivalent high performance, mpi parallelized C code. 

It is worth noting that there are a few bugs and edge cases that have not been worth the time to develop fixes for as they can be easily corrected after the fact.

Additionally, there is a conflict between the native handling of fourier transforms in Julia and C. Julia automatically normalizes the fourier transforms while C does not. As I have used this primarily to code directly into HPC I generally encode a normalization into the julia code which means that it would not correctly simulate if run in Julia.

Known bugs:
- When replacing the caret symbol for powers with the power function, the base is not recursively passed through replacements. No current plans to correct this.
- Not currently defining all the various functions in the header. Again, I may fix this later.

Known edge cases:
- The calculation of the square of the wave vector is not properly transposed such that kx is calculated from the second index and ky by the first index. This must be done manually. This does not seem worth the dev time to correct.
- Similarly, the code has no way of knowing when a loop index requires global indexing. Since this is rare, generally only needed for the calculation of the wavevector and potentially the seed I've chosen to simply do this manually.
- I have not written a convenience function for reading real valued fourier space arrays. Given that only real space arrays are required to restart a simulation and data file size can scale quickly with PFC simulations I doubt I will correct this.

