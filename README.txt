**************************************************************************************
PA4  Implementation
**************************************************************************************
THIS IS A GROUP SUBMISSION:
**************************************************************************************
Harshit Khaitan : hkhaitan@stanford.edu
Avinash Parthasarathy : avinash1@stanford.edu
**************************************************************************************
The fourier transform functions on X-Axis and Y-Axis on real image is
performed by using OPENMP with number of threads set to 8. We are seeing
Max performance with 8 threads since we have 8 cores on Amazon EC2
machines.

We are seeing a speedup of 50X using Makefile. The parallel threads are spawned on Inner
for-loop rather than outside for-loop to acheive maximum performance and
also to utilize shared memory cache. The outside for-loop is optimized serially to run 1024/8
times by copying the subsequent 8 rows/columns depending upon X-axis or
Y-Axis computation. This copying of original image to buffer expoits
shared memory so that required array elements are present for individual 
worker threads to work parallely. The Reading of array element happen 
parallely from buffer and writing after computation happens parallely to 
original image. The thread scheduling is set as dynamic to maximize
parallelization so that all threads are busy until completion of
computation. The number of chunks to schedule clause is set to size_x/OMP_NUM_THREADS.
