# Building a Maximum Likelihood Tree

To prep the ML tree, first make list to pull from GenBank, which includes our sequences. Here is the call we use to download all sequences for DENV-1 and DENV-2:


```
rm(list=ls())

homewd="/Users/carabrook/Developer/cambodia-dengue-national/"
setwd(homewd)

dat.new <- read.csv(file = paste0(homewd, "figure-development/FigS4/ML-Sequences.csv"))
denv1 <- paste(dat.new$Accession[dat.new$Serotype=="DENV-1"], collapse=", ")
denv2 <- paste(dat.new$Accession[dat.new$Serotype=="DENV-2"], collapse=", ")

denv1 <- paste0("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=text&id=", denv1)
denv2 <- paste0("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=text&id=", denv2)


```
Then, make sure you attach one "outgroup" to each tree. Root them both in the DENV4 reference sequence (Accession # NC_002640)


```
ref <- " NC_002640"
denv1 <- paste0(denv1, ", ", ref)
denv2 <- paste0(denv2, ", ", ref)


```

After downloading the files for DENV1 and DENV2 into our browsers, we store each compiled fasta in the 'ML-tree' folder, then load and rename:

```
rm(list=ls())
library(seqinr)

homewd="/Users/carabrook/Developer/cambodia-dengue-national/"
setwd(homewd)

denv1.seq <- read.fasta(file = paste0(homewd, "figure-development/FigS4/all-DENV1-ML.fasta"), as.string=T, forceDNAtolower=F)

denv2.seq <- read.fasta(file = paste0(homewd, "figure-development/FigS4/all-DENV2-ML.fasta"), as.string=T, forceDNAtolower=F)

dat.new <- read.csv(file = paste0(homewd, "figure-development/FigS4/ML-Sequences.csv"))
dat.new$rename <- paste(paste(dat.new$Accession, dat.new$Locality, sep = "_"), dat.new$Year, sep="_")

#and rename
denv1.names <- dat.new$rename[dat.new$Serotype=="DENV-1"]
denv2.names <- dat.new$rename[dat.new$Serotype=="DENV-2"]

denv1.names <- c(denv1.names, "NC_002640_DENV4")
denv2.names <- c(denv2.names, "NC_002640_DENV4")

write.fasta(denv1.seq, names = denv1.names, file.out = paste0(homewd, "figure-development/FigS4/all-DENV1-ML.fasta"), as.string=T)

write.fasta(denv2.seq, names = denv2.names, file.out = paste0(homewd, "figure-development/FigS4/all-DENV2-ML.fasta"), as.string=T)

```


Then send to [MAFFT](https://mafft.cbrc.jp/alignment/server/) for alignment to produce the files "align-DENV1-ML.fasta" and "align-DENV2-ML.fasta". Then send these to ModelTest-NG to determine the best nucleotide substitution model for the data. Here's the script to use on the computing cluster:

```

#!/bin/bash
#SBATCH --job-name=denv2ML
#SBATCH --account=fc_sirmodel
#SBATCH --partition=savio3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00


module load vim/7.4 
module load java/1.8.0_121
module load emacs/25.1 
module load cmake/3.15.1
module load python/3.6 
module load gcc/7.4.0
module load openmpi/4.0.1-gcc

export CC=`which gcc`
export CXX=`which c++`


./bin/modeltest-ng -i dengueML/align-DENV2-ML.fasta -t ml -p 8

```

I saved the results of ModelTest-NG in subfolder "dengueML" in the FigS4 directory.

Once ModelTest-NG finishes (it can take several hours), it is time to build a maximum likelihood tree using RAxML. See documentation on their website for how to get this running on your home computer and/or computing cluster. In my case, ModelTest-NG told me that the best support was recovered for a GTR+I+G4 model for both DENV-1 and DENV-2 phylogenies, so that is what I ran in raxml. I followed the RAxML tutorial to build a tree from my MSA, first checking that RAxML could read the alignment:

```
/global/home/users/cbrook/raxml-ng/raxml-ng-mpi --check --msa align-DENV2-ML.fasta --model GTR+I+G4 --prefix T1

```
...then parsing the alignment to find the appropriate number of threads (12) with which to run RAxML:

```
/global/home/users/cbrook/raxml-ng/raxml-ng-mpi --parse --msa align-DENV2-ML.fasta --model GTR+I+G4 --prefix T2
```
Finally, I kicked off RAxML (including bootstraps) with the following script:

```
#!/bin/bash
#SBATCH --job-name=denv1
#SBATCH --account=fc_sirmodel
#SBATCH --partition=savio3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00

module load vim/7.4
module load emacs/25.1
module load python/3.6
module load java/1.8.0_121
module load cmake/3.15.1
module load gcc/7.4.0
module load openmpi/4.0.1-gcc

/global/home/users/cbrook/raxml-ng/raxml-ng-mpi --all --msa align-DENV1-ML.fasta --model GTR+I+G4 --prefix T3  --seed 12 --threads 8 --bs-metric fbp,tbe

```

Once RAxML finished (a few hours later), I imported the resulting tree into R and built a phylogenetic tree for Fig S4. Note that for quick viewing, you can easily view the tree in FigTree. See the script, FigS4 for instructions on how to build the tree.



