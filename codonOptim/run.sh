#!/bin/bash
#mCherry
#mRuby2
#mVenus
#mCerulean
#lacZ

python optimise.py -s PCA -v 5 --pca-groups "bg,0.49,-0.03,1.5;rb,-2.1,3.4,1.5;ph,-2.8,-1.073,1.25" stats/B.Sub_168 data/*.fasta --amino -o out/
#python optimise.py -s demo -v 5 stats/B.Sub_168 data/*.fasta --amino -o out/

