# LocalSignallingEntropy
Ongoing Project


Entropy matrix should be produced first in order to produce bootstrapping, perm, and pertubation testings.

Permutation test has two similiar for faster results. You can run either of them. One is _FAST — single-threaded but optimized. Pre-computes the neighbor list once and reuses it. Runs on one CPU core. The other one is _PARALLEL — uses doParallel and foreach to distribute the 100 permutations across multiple cores (default 4). Each core runs ~25 permutations simultaneously.
Requirement (minimum) for this script to work: 2+ Cores CPU, 32 GB RAM.
