Scripts for using ARGs sampled by *ARGweaver* to analyze patterns of species differentiation
=======

This directory includes all scripts used to process *ARGweaver* output files (.SMC) and to analyze ARG-based measures along the genome. All scripts are in ``R`` of ``BASH``. The analysis pipeline is split into two components, described below.



------------


smc_to_stats
------------

This directory contains scripts that take *ARGweaver* ouptut files and create windowed summary stats.  Before running *ARGweaver*, long genomic scaffolds are typically split into overlapping ARG blocks. So the input to this pipeline is an *ARGweaver* output file (.SMC) per ARG block.The final product is a statistic table (.stat.bed) per genomic scaffold, with one line of stats per window. The pipeline summary:

**Input:**   SMC files (*.smc.gz)							                    - one per ARG block

**Step 1:**  Filter and trim the ARG blocks

**Step 2:**  Generate BED files  (*.bed.gz + *.bed.gz.tbi)	- one per ARG block, after filtering for length

**Step 3:**  Generate TREE files (*.tre.gz)							   - one per ARG block, after trimming

**Step 4:**  Generate STAT files (*.stat.gz)							  - one per ARG block (same # lines as TREE file)

**Step 5:**  Generate STAT BED file (*.stat.bed)					- one per genomic scaffold, line per window (20 kb)

For more details, see *smc_to_stats/README.md*


------------

summarize_stats
------------

This directory contains scripts that take windowed stat summary files (.STAT.BED) and intermediate TREE files and analyze patterns of species differentiation along the genome. This implements **tests 1**, **2**, and **3**, used in **Table 1** of our manuscript. Additional scripts are given for plotting summaries and trees.

For more details, see *summarize_stats/README.md*


------------

sampleFiles
------------

This directory contains examples of input files that we used in the analysis. These are referenced in the appropriate part of each pipeline.


------------

<span style="color:darkblue">Prerequisites:</span>
------------

The R scripts in this repository use the following packages:

* ``ape``
* ``phytools``
* ``castor``
* ``plyr``

Additional programs used in the pipeline:

* ``tabix`` (for creating index for BED files)

* ``gzip``and ``bgzip`` (for compressing BED, TREE, and STAT files)
* ``smc2bed`` (for converting *ARGweaver* output from SMC to bed format; from *ARGweaver* package)


------------