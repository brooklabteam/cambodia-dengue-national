# This folder contains the script to generate Figure 5 of the main text.

The script "make-transmission-trees.R" loads the phylogenies from Fig4 and creates transmission chains based on phylogenetic distance and most recent common ancestor. It results in the file "AllDENVtransTreeDat.csv" which is stored in the "data" folder of the main repo.

The script "fig5.R" loads this csv file to produce Figure 5 of the main text.
