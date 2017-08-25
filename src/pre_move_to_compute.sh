#! /bin/bash

# move tries from warm storage to compute
# warm storage profile for tries: /warm/leslie/pritykin/projects/CRISPR/<organism>/<organism>_all/<organism>_all_trie.dat
mv /warm/leslie/zamparol/crisprML/tries/mm10_all_trie.dat /data/leslie/zamparol/crisprML/tries/mm10_all_trie.dat

# move kmers from warm storage to compute
# warm storage profile for kmer files: /warm/leslie/zamparol/crisprML/kmers/mm10/Cas9/mm10_all_kmers.txt
mv /warm/leslie/zamparol/crisprML/kmers/mm10/Cas9/mm10_all_kmers.txt /data/leslie/zamparol/crisprML/kmers/mm10/Cas9/mm10_all_kmers.txt

