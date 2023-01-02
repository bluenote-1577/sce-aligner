# sce-aligner 
A basic seed-chain-extend aligner with linear-gap cost chaining and quadratic time extension for experiments. This is the tool used for our paper "Seed-chain-extend alignment is accurate and runs in close to1
O(m log n) time for similar sequences: a rigorous average-case analysis" for aligning nanopore reads to references.

Also included in this repository are the scripts for generating our figures. 

## Installing the aligner

```
git clone https://github.com/bluenote-1577/sce-aligner
cd sce-aligner
cargo build --release
./target/release/sce_aligner -h
```

To use the aligner on real sequences, do

```
./target/release/sce_aligner reference.fa query.fastq
```

The output will look like
```
extension_time  chaining_time read_length aligned_fraction  read_name
```
for each read. The extension time and chaining time will be -1 if the aligned_fraction is less than 90%. It will be -2 if there is a large gap in the chaining (> 10kb). it will be -3 if there are too many anchors. 

## Plotting 

### get_passing_reads.py

This script takes in a list of sam files and outputs reads that have gap-compressed identity close to 95%. 

### good_plot.py

This script plots the extension and chaining times. Preprocessed results from sce_aligner are already present in this folder, so it can be run without
regenerating intermediate files. 

To regenerate the intermediate files, do `sce_aligner human_ref.fa human_reads.fastq > human_results.txt` etc. The reads and references used in our study can be found
in the supplementary table of our paper. 

### recov_plot.py

This script plots the aligned fraction as a function of read length. 
