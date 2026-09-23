---
layout: post
title: "How to Build and Run a Genomic Sequencing Demultiplexer"
date: 2026-09-22
category: technical-writing
---

## Tool Overview
In this tutorial, you'll build a Python-based demultiplexing pipeline for paired-end Illumina sequencing data with dual indexes. The pipeline will classify reads based on their index combinations and write them to separate output files. It will also handle instances of index hopping and invalid indexes. 

### What You'll Learn
* Parse FASTQ files
* Match dual indexes to samples
* Identify index hopping and invalid indexes
* Write classified reads to output files
* Package the pipeline as a command-line tool with `argparse`
* Test the pipeline with `pytest`

## Background
Genetic sequencing is expensive both in time and resources. If we were to sequence the fragments from one sample at a time, it would be very inefficient. Thus, it’s common to multiplex–-that is, pool fragments from tens, hundreds, even thousands of different sample groups and sequence them in a single run. This also limits slight variations in reagents, temperature, or the sequencing machine between runs. Normally, these would result in batch effects, differentiating the samples based on minor experimental/protocol changes rather than real, biological differences.

We still need a way to separate the reads back into the sample groups they originated from after the sequencing process is complete to analyze our data by comparing across and within samples and pursue further experimental explorations. This is accomplished with **indexing**, the process of tagging every nucleic fragment with a short DNA sequence prior to sequencing that corresponds to its sample. Thus, the reads generated from the fragments will contain their indexes as well. We can use these indexes to group the sequenced reads by sample. This is the process we call **“demultiplexing”**. There are a few different ways that indexing is implemented.

## Illumina Indexing
Every library preparation kit handles indexing differently. For our purposes, we will focus on the Illumina protocol. Illumina indexes are located between the DNA insert and the sequencing adapters. They are sequenced in separate runs from the insert. There is a further distinction between single-indexed and dual-indexed libraries.

Single-indexed libraries contain one index sequence appended to the end of the fragment. In Illumina sequencing, this would be the side that contains the P7 adapter. In this case, there would be one additional sequencing run that encodes Index 1 (I1) after the run for R1.

Dual-indexed libraries contain two index sequences appended on either side of the DNA fragment. In this case, there will be two additional sequencing runs that will each capture I1 and Index 2 (I2) separately after R1. There are two approaches to dual-indexing which differ primarily on how many unique index sequences are used. **Combinatorial dual indexing** is designed such that individual index sequences may repeat across sample groups, but the combination of I1 and I2 sequences are unique to each sample group. **Unique dual indexing (UDI)** uses separate, non-redundant pairs of I1 and I2 for each sample group. The index of one pair never reoccurs as part of another pair. A major advantage of UDI is that it makes it much easier to detect and filter cases of **index hopping**, which is a rare phenomenon where a fragment gets attached with the wrong index on either or both ends. Usually, it occurs as a result of excess free adapters during library prep. If gone undetected, it can cause sample contamination and skew downstream analyses.

## Understand the Input Data
### FASTQ Files
FASTQ is one of many bioinformatics file formats. It holds base calls from the sequencer along with a quality score indicating its confidence in the base call at each position. Each fragment that was sequenced is stored as a record, which have the following structure:

```
@header
Base calls for this fragment
+
Quality scores (encoded as ASCII symbols; read more about Phred scores).
```

Note that the header is the first line of the record, and always begins with the "@" character. The second line contains the sequence as read by the sequencer. The third line usually just contains the "+" character, though it can also contain other metadata. The fourth line is the per base quality score.

In the case of paired-end, dual indexed data, the different components of one fragment are stored across four FASTQ files. They are stored in the order that they are sequenced. The forward read (read 1; R1) is stored in the first file. Next, index 1 (I1). Then index 2 (I2), and finally, the reverse read (read 2; R2). 

### Index File
We need to tell the tool the expected index sequences. This information will be stored in a text file, which lists the sample name and its corresponding I1 sequence, tab-separated:

```
A1  ATTGCACC
B2  TGGCTACA
...
```

## Assumptions
Here are the assumptions we’re making for this tool:
### Input
* Input is from paired-end Illumina sequencing data.
* UDI where the reverse complement of I2 = I1 was used to label sample groups.
* Data is contained across 4 FASTQ files: R1, I1, I2, and R2.
* List of expected indexes is provided. 

### Demultiplexing logic
If the reverse complement of I2 = I1, then that reads belongs to the group encoded by that index. 
Index hopping is when I2 is not the reverse complement of I1, but both indexes exist in the set of expected indexes.
Indexes that contain “N” nucleotide calls or are not in the set of expected indexes are considered unknown/invalid.

## Set Up the Project
Directory structure:
```
demultiplex/
├── src/
│   └── dmux/
├── tests/
├── data/
├── pyproject.toml
└── README.md
```

## Build Demux Tools
As a rule of thumb, we want to modularize the pipeline into functions that accomplish subprocesses. This makes it much easier to test and debug the pipeline.

### Implement a Reverse Complement Generator
We will be reverse complementing DNA sequences many times in order to classify reads, so let's make it a function.

```python
def reverse_complement(sequence):
  # hashmap of complementary bases
  comp_bases = {"A":"T", ... "N":"N"}

  rev_comp = ""

  # check if the base is a key in the hashmap; only continue if it does, else raise error
  for base in sequence:
    if base not in comp_bases:
      raise ValueError

    rev_comp += comp_bases[base]

  return rev_comp

```

### Implement a FASTQ Parser
Rather than loading an entire FASTQ file into memory, the pipeline will load records one at a time using a generator.

```python
def fastq_parser(fastq_file):
  while True:
    header = fastq_file.readline().strip("\n")
    # break the loop at the end of the file, at which point header will be an empty string, which equates to 'False'.
    if not header:
      break
    seq = fastq_file.readline().strip("\n")
    plus = fastq_file.readline().strip("\n")
    qscore = fastq_file.readline().strip("\n")

    # use of yield makes the product of the function a generator, so the position of the pointer is not reset at the end of the run; the next time it runs, it will pick up from where it left off last time,      generating the next record.
    yield header, seq, plus, qscore

```

The majority of memory consumption for this tool will be through reading the FASTQ files. These files can get very large, so we need to build our tool such that regardless of the input data size, a controlled amount of memory is used. Processing records incrementally keeps the memory usage low and allows more efficient handling of large sequencing files. 

### Implement Index Matching
In order to classify reads, we need to compare I1 and I2. As outlined earlier in demultiplexing rules, we use if-else statements to identify which condition the read meets.

```python
index1 = i1_record[1]
index2 = i2_record[1]

# get the reverse complement of index 2
rc_index2 = reverse_complement(index2)

# now the classification process using index1 and rc_index2
if (
    "N" in index_seq_1
    or "N" in index_seq_2
    or index_seq_1 not in indexes
    or rc_index_seq_2 not in indexes
):
  ...

elif index1 == rc_index2:
  ...

elif index 1 != rc_index2:
  ...

```

### Write Classified Reads to Output Files
Before the classification process, the pipeline will create all output FASTQ files. Each expected index will have output FASTQ files named after it. Additionally, there will be hopped.fastq files to hold index hopped reads and unknown.fastq for reads with unknown reads.

```python
pass
```

Notice that each class has two output files: R1 and R2. R1 and R2 must be synchronized; these will hold the forward and reverse read, respectively, of that fragment. In other words, the record from input R1.fastq should be written to its output R1 FASTQ and its counterpart from input R2.fastq will be written to the corresponding output R2 FASTQ.

Further to this, we will ensure the index information is not lost after demultiplexing by appending I1-I2 pairs to the headers of each record. This is useful particularly in hopped and unknown instances, so that the user can use this information to make decisions on a case-by-case basis of these records.

### Implement a New Header Generator
We need a function that handles appending the indexes to the header.

```python
def create_new_header(header, index1, rc_index2):
  pass
```

## Add a Command Line Interface (CLI)
We will use the Python package `argparse` to build a CLI for this tool. 

```bash
demultiplex \
  -r1 R1.fastq.gz \
  -i1 R2.fastq.gz \
  -i2 R3.fastq.gz \
  -r2 R4.fastq.gz \
  -i indexes.txt \
  -o path/to/output/directory/
```

## Add Tests
### Unit Tests
Test individual functions to ensure subprocesses are being conducted as expected:

```python
def test_reverse_complement():
  pass

def test_fastq_parser():
  pass

def test_demultiplex():
  pass
```

### Input Validation Test 
We should also create tests that ensure input FASTQ files are not malformed or corrupted. By default, it will run on the small, synthetic test dataset provided in the project repo, but it can also be applied to user-specified data.

```python
def test_inputs():
  pass
```

### Integration Test
Test the synchronicity of the entire pipeline end-to-end, from parsing user input to generating output:

```python
def test_cli():
  pass
```

### Validate the Results
Using our test dataset, we should also have an expected output, which we will ensure is the output produced by the tool when the test data is provided as input.


