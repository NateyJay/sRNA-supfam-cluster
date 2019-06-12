# sRNA-supfam-cluster
Clusters sRNA sequences by edit distance, nominally into superfamilies

#### Motivation

Micro RNAs (miRNAs) and other small regualtory RNAs (sRNAs) can be divided into families, usually defined by near perfect sequence similarity (edit distance maximum of 1 or 2). However, more distant relationships can exist between sRNAs (Xia et al. 2013 - Plant Cell). This script allows to group sRNA sequences into superfamilies, based on larger user-defined distance cutoffs. The script takes a fasta-file as input and produces files grouping the sRNAs into superfamilies, based on 1 or more cutoff distances.

The function used for calculating edit distance is a modified version of Hamming distance, which allows for base mutations and shift errors when measuring distance, but heavily penalizes insertions and deletions. We expect that shifts or base mutations may commonly occur with related sRNAs, but indels likely do not as they would be expected to invalidate the sRNAs targeting relationship.

#### Usage

    01-cluster_by_mod-hamming.py \
      -file [input.fa/input.fasta] \
      -cutoff [integer_1] [integer_2] [...] # allows the user to supply multiple cutoffs
      
Example:

    python 01-cluster_by_mod-hamming.py \
      -file test_data.fasta \
      -cutoff 1 2 5 9


#### Output

Fields indicate:
1. sRNA name
2. sRNA sequence
3. assigned superfamily

