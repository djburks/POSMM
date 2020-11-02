<p align="center"><img src="assets/posmm.svg" width="500" title="POSMM Logo"></p>

# POSMM
Python-Optimized Single-Order Markov Model 
Metagenomic Read Classifier for GNU/Linux

POSMM (pronounced 'possum') is an alignment-free metagenomic read classifier.  Intended for short reads (e.g. 30-500nt), POSMM uses single-order Markov models and logistic regression models to assign taxa with an associated confidence at Phylum, Class, Order, Family, Genus, and Species ranks.  Much like its namesake, POSMM excels at being a highly-adaptive, downstream application that can take advantage of the 'leftovers' from ultrafast upstream classifiers such as Kraken or CLARK, providing higher sensitivity than alignment-based classification programs while not wasting time on reads already classified.


## Installation

Download the POSMM-1.0.tar.gz archive from the dist folder, and simply install with pip.
Dependencies are taken care of by pip, but you will need a relatively modern GCC version for compiling SMM.

```
pip install POSMM-1.0.tar.gz
```
## Features

- **Database Free**:  POSMM does not have a database building stage.  Markov-chains are built directly from genome fasta files, so each run can be tailored to a specific analysis.
- **Confidence Scores**: Confidence values ranging from 0 to 1 are produced for each taxonomic rank.  
- **Genome Model Downloader**:  Customized and quick-setup options are available for downloading genomes directly from NCBI RefSeq.  Sharing model groups is also extremely straightforward using simple POSMM taxlists.
- **Lineage Assignment**:  Lineage assignments are provided directly from NCBI Taxonomy and based on the Taxid of each organism.
- **Multiple Group Management**: POSMM makes it easy to switch genome sets, and create new ones.
- **Multi-thread Support**:  POSMM can take advantage of multi-core / multi-thread systems for both model comparison and confidence scoring.

## Options

To see all available options and descriptions for their use, just ask for help.
```
POSMM -h
```

## Setting up Genome Sets
### Building a Quick-Setup Genome Set from RefSeq

Before analyzing reads, you need to set up at least one set of genomes for POSMM.  By using the setup runmode, you can specify whether you would like to download all available prokaryotic reference or representative genomes.  By default, genome sets are written to ~/POSMM/Genomes, but this can be changed with the --gdir parameter.
```
POSMM --runmode setup --gtype reference --gdir ~/POSMM/Reference_Genomes
```
### Building a Custom Genome Set from RefSeq
By providing POSMM with a list of GCF numbers, you can control which genomes are pulled from RefSeq.  This can be useful if an upstream classifier narrows down a set of reads to a particular taxonomic group (phylum, class, etc), or if you only care about isolating reads to a specific set of organisms.
```
POSMM --runmode setup --gtype custom --gdir ~/POSMM/Custom_Genomes  --taxlist GCFList.txt
```
The taxlist needs to be a text-file with one GCF annotation per line:
```
GCF_000006785.2
GCF_000006885.1
GCF_000007265.1
```
### Skipping Downloads from RefSeq

Already have your genomes downloaded?  POSMM assigns reads to genomes based on their taxid, so just make sure the first part of the name is the correct taxid followed by a period (e.g. 1311.GCF_000007265.1.fna).  The rest of the name can be whatever you like.  Just rename and drop them in a directory of your choosing, and use the --gdir parameter to make sure POSMM finds them.

## Analyzing MetaFasta Reads

Once you have a set of genomes ready, running POSMM is easy.
```
POSMM -f metagenomic_fasta_reads.fna -out posmm_analyze_reads.txt
```
By default, POSMM will use the genome set in ~/POSMM/Genomes.  This can always be changed with the --gdir parameter, just as when downloading genomes.
POSMM also uses the highest order (12) for generating SMMs for highest accuracy.  While the most accurate, it is also the slowest, and can be dropped to 10 for even faster runs.

## Multi-Threading Support

The SMM analysis and confidence scoring stages are spread across all available logical cores.  Please note that more threads == more RAM use, so either split your metagenomic fasta files before analysis, or use a lower number of threads with the --thread/-t parameter.  Lower orders also use less RAM, such as --order 10, but note that this leads to less accurate classifications.

POSMM is heavily CPU-bound, and will use 100% of the CPU if you let it.  As such, we have found that setting -t to the number of physical cores sometimes results in better performance on some systems.

## Output

Output is a tab-delimited file of taxa and confidence score pairs, separated by ":::". 
Given the size of metagenomic datasets, output is kept as minimal as possible.  The line number corresponds to the read order of the original metagenomic fasta file.
You can also export a list of taxonomic assignments and the raw Markov model score with --runmode raw.

## Notes on Confidence Scores

Confidence scores help avoid false-positive classifications inherent to Markov model methods of classification.  As a general rule, we recommend using a cutoff >0.50 for your own analysis.  Using a cutoff of 0.25 is similar to using no cutoff, as this will remove obvious false-positives such as assignment of human contaminant reads to bacterial genomes.  Think of this like using Kraken without a --confidence cutoff. The reliability of using low (i.e. < 0.5) cutoffs depends on the representative quality of your genome set.  

## Future Updates
POSMM is far from fully evolved.  Some of the upcoming plans include:
- Custom (Non-RefSeq) Genome Set Support
- More Output Options 
- PyPI Integration / Conda Packages
