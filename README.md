# Succinct Dictionary

This is a succinct dictionary data structure that supports rank and select in O(1) time  while still using space that is very close to the information theoretic limit. The concept of the succinct data structure was introduced by Jacobson to encode bit vectors. 
This provides a two level indexing system for the dictionary. This data structure is modeled after the one introduced by Sebastiano Vigna for 64 bit systems: Ref: Broadword Implementation of Rank/Select Queries: Vigna, Sebastiano. "Broadword implementation of rank/select queries." International Workshop on Experimental and Efficient Algorithms. Springer, Berlin, Heidelberg, 2008. Not theoretically optimal by much more practical than other algorithms Runs much faster.
