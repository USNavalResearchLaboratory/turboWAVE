Parallelism
===========

In order to run turboWAVE efficiently, a basic understanding of the parallel programming model is necessary.  TurboWAVE uses a combination of distributed memory and shared memory methods.  The distributed memory method is called "domain decomposition", and corresponds to physically partitioning the simulation region into chunks that can be worked on by different processors.  The memory that holds each chunk is distributed in the sense that it may be connected only by a network cable.  The shared memory method is called "fork-join".  In this approach independent software threads are running on different cores, but in a setting where each core has equal access to the memory (e.g., the threads could be on different cores of the same multi-core processor).

The distributed memory model is implemented using software called MPI.  The shared memory model is implemented using software called OpenMP.

As the user, you have to choose how to partition the problem.  In particular, you must choose the number of MPI processes, and the number of threads.  If you opt to have more than one MPI process, you may also want to control exactly how the physical simulation domain is partitioned into chunks.

On HPC systems in particular, the term "node" refers to a set of physical processor cores that have access to the same shared memory.  Often you must specify resources requested in terms of nodes.
