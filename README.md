This project implements 3 flavors of parallel google pagerank as described follows. 

The first two flavors (mpi_parallel and pthread_pagerank) use the same concept of considering the input web graph whose pagerank needs to 
be calculated as a Compressed Row Storage matrix(useful for sparse matrices), anda normalized pagerank vector. The graphs  are of size 
200 K and 4 M, and graphs are partitioned according to the number of pthreads/cores.Used METIS (set of serial programs for partitioning 
graphs) The pagerank is then computed through iterative multiplication with the matrix (power iterations). The process continues until the
pageranks do not change significantly (the two-norm of the pagerank vector changes by less than 10^{-5}. 

In the MPI approach (message passing), one processor reads the file and communicates equal number of nodes to all other processes,The processes
initialize their pagerank values and computer pagerank vector by communicating the necessary elements of the vector to processors 
that need them using MPI's send recv constructors, and not broadcasting the entire intermediate pagerank vector to all processors.

In the PThreads approach(shared space),  each thread takes care of it's own part of pagerank vector and matrix, by modifying the matrix 
and the vector stored in memory by accessing the data concurrently taking barriers to avoid data inconsistency. 

The third flavor(pthread_randomwalk),uses pthreads, where each thread has a set of walkers, on each iteration,it picks its walker, 
chooses a random surrounding node and jumps to it, and then updates a global visited counter for the node it lands up on, the pagerank 
is the number of visits by a walker in a fixed num of iterations


    
