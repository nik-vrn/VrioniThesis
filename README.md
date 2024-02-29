# VrioniThesis -Update of Structures for Finding Optimal Routes in Dynamic Road Networks

A Dynamic Contraction Hierarchy Weight Decrease algorithm is designed to efficiently handle edge weight reduction in dynamic road networks. This algorithm includes two main components: preprocessing and dynamic handling of queries and updates in microseconds.

Requirements
Before running the code, ensure you have the following dependencies installed:

PGL: Library for graph algorithms and data structures (https://www.ceid.upatras.gr/webpages/faculty/zaro/software/pgl/index.html).

Routingkit: (https://github.com/RoutingKit/RoutingKit).

Boost library: (https://www.boost.org/)

DIMACS9: Using graphs in DIMACS9 format (http://www.diag.uniroma1.it/challenge9/download.shtml).

The Makefile contains the correct path to the pgl/include folder
The Makefile contains the correct path to the boost header files folder (you have to install boost first, speciffically libboost-dev and libboost-program-options)

Run the preprocessing script to generate the contraction hierarchies (CH) data file (ch.dat) for the desired graph and use it for dynamic weight decreade updates and queries.
