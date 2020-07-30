Scripts for using ARGs sampled by *ARGweaver* to analyze patterns of species differentiation - smc_to_stats
=======

This directory contains scripts that take *ARGweaver* ouptut files and create windowed summary stats.  Before running *ARGweaver*, long genomic scaffolds are typically split into overlapping ARG blocks. So the input to this pipeline is an *ARGweaver* output file (.SMC) per ARG block. The final product is a statistic table (.stat.bed) per genomic scaffold, with one line of stats per window. The pipeline summary:

**Input:**   SMC files (*.smc.gz)							                    - one per ARG block

**Step 1:**  Filter and trim the ARG blocks

**Step 2:**  Generate BED files  (*.bed.gz + *.bed.gz.tbi)	- one per ARG block, after filtering for length

**Step 3:**  Generate TREE files (*.tre.gz)							   - one per ARG block, after trimming

**Step 4:**  Generate STAT files (*.stat.gz)							  - one per ARG block (same # lines as TREE file)

**Step 5:**  Generate STAT BED file (*.stat.bed)					- one per genomic scaffold, line per window (20 kb)


------------

Setup
------------

- Run *ARGweaver* on your data sets. When splitting long scaffolds / chromosomes into smaller ARG blocks (~2Mb), make sure to mainain an overlap of ~100 kb, so that 50 kb can be trimmed from the edges.
- Put all *ARGweaver* output files  (.SMC.GZ ) in a single smc source dir
- Put all *ARGweaver* output log files  (.LOG) in a single log source dir
- Create directory *infoTables/*
- Prepare sequencing species-individual key file (*infoTables/individual-species-key-modified.txt*). This file contains a table with a row per haploid sample in the data set (120 in our analysis) with the following columns: *Species*  (species name)   *Individual_hap* (haploid id)    *Color*  (RGB color code for species). See example in *sampleFiles/individual-species-key-modified-sporophila.txt* in root directory. 
- Follow the specific setup instructions of each step in the pipeline below


------------

Step 1 - Filter and trim ARG blocks
------------

**Input:**

- *ARGweaver* output files  (.SMC.GZ ) 

**Output:**

- *infoTables/ARGblock-coordinates.txt*   -   (columns: <SMC-file-path> <scaffold> <MCMCiter> <startIndex> <endIndex>)  
- *infoTables/ARGblock-coordinates-trimmed.txt*   -   (same format but coordinates after trimming and contains filtered set of ARG-blocks)  
- *infoTables/ARGblock-coordinates-info.txt*   -   (free text file with information on filtering) 

**Steps:**

``==> bash create_arg_block_region_file.sh``

- **Setup:**  set smc source dir in script
- Creates ARG block coordinate file based on SMC files.

`==> Rscript trim_arg_blocks.R `

- **Setup:**  set filtering and trimming parameters in script `trim_arg_blocks.R ` : (1) MCMC sampling iteration for filtering (default=1000); (2) minimum ARG block length for filtering (default=100,000); (3) window length (default =20,000); (4) trim size in each end of an ARG block (default=50,0000). The trim size should be at most half the length of the overlap between consecutive ARG blocks in the same scaffold to ensure complete coverage.
- Filters out short ARG blocks and SMC files generated in early MCMC iterations
  - Trims edges of each ARG block (according to overlap)
  - Computes a new region for each ARG block so that it contains an integer number of windows of given length (20kb in our default settings)


------------

Step 2 - Generate ARG bed files
------------

**Input:**

- *ARGweaver* output files  (.SMC.GZ ) 
- *infoTables/ARGblock-coordinates.txt*   (generated in **Step 1** above)

**Output:**

- *argBedFiles/<ARGblock>_out.<MCMCiter>.bed.gz*   -   (ARG BED files, one per ARG block)

  [ format for line in file: <scaffold> <startPos> <endPos> <MCMCiter=0> <Newick_tree> ]  

- *argBedFiles/<ARGblock>_out.<MCMCiter>.bed.gz.tbi*   -   (index files for ARG BED files)

- *infoTables/arg_blocks_with_no_bed.txt*   -   (list of ARG blocks that are filtered out and will not have a bed file) 

**Steps:**

``==> Rscript make_command_lines_for_beds.R``

- **Setup:**  set smc source dir & log source dir in script; set MCMC sampling iteration (default=1000) & minimum ARG block length for filtering (default=100,000); set the number of parallel jobs that can run on your machine (default=20)
- Creates bash script ``create_bed_files.sh`` in current directory with ``Rscript`` calls to ``smc_to_bed.R`` for all ARG blocks after filtering. Adds ``wait`` command every *numJobs* jobs.

``==> bash ./create_bed_files.sh`` 

* **Setup:**  set appropriate path for  ``smc2bed`` program in the system call in script  ``smc_to_bed.R``.
* Runs ``smc_to_bed.R`` on all SMC files. Each run applies the  ``smc2bed`` program + compression with ``bgzip`` + indexing with ``tabix``.


------------

Step 3 - Generate Tree files
------------

**Input:**

- ARG BED files  (.BED.GZ ; generated in **Step 2 above) 
- *infoTables/ARGblock-coordinates-trimmed.txt*  - (generated in **Step 1** above)

**Output:**

- *argTreeFiles/<scaffold>.<startInd>-<endInd>.<MCMCiter>.tre.gz*   -   (ARG TREE files, one per ARG block, ~15x smaller than BED)

  [ format for line in file: <position> <tree> ] 

**Steps:**

``==> Rscript make_command_lines_for_trees.R``

- **Setup:**  set interval between sampled trees (default=500); set the number of parallel jobs that can run on your machine (default=20)
- Creates bash script ``create_tree_files.sh`` in current directory with ``Rscript`` calls to ``bed_to_tre.R`` for all ARG blocks with BED files (unless empty after trimming). Adds ``wait`` command every *numJobs* jobs.

``==> bash ./create_tree_files.sh``

* Runs ``bed_to_tre.R`` on all BED files. Each run applies ``tabix`` to extract trimmed regions in BED file, and then extracts trees in specified intervals (default=500). Resulting file is compressed with ``gzip``.


------------

Step 4 - Generate Stat files
------------

**Input:**

- ARG tree files  (.TRE.GZ ; generated in **Step 3** above)
- *infoTables/individual-species-key-modified.txt* - to map haploid samples to species

**Output:**

- *argStats/<scaffold>.<startInd>-<endInd>.<MCMCiter>.stat.gz*   -   (ARG STAT files, one per ARG block, 10-20 times smaller than TRE)

  [ format for line in file (one line per tree): <TMRCAH-all> <RT12> <spe1-RTH>  . . . <spe1-enrich>  . . .] 

  [  these are all statistics used in test 1 (enrichment scores) & test 2 (RTH'); Cross coalescence times are not windowed ]

**Steps:**

``==> Rscript make_command_lines_for_stats.R``

- **Setup:**  set the number of parallel jobs that can run on your machine (default=20)
- Creates bash script ``create_stat_files.sh`` in current directory with ``Rscript`` calls to ``tre_to_stat.R`` for all ARG blocks with TRE files. Adds ``wait`` command every *numJobs* jobs.

``==> bash ./create_stat_files.sh``

* **Setup:**  set total number of haploid samples (*total_n*) and number of haploid samples per species (*species_n*) in  ``tre_to_stat.R`` to fit your data set. Note that RT12 will be set according half of *species_n*. If different species have different *species_n*, then you may need to further adjust script.
* Runs ``tre_to_stat.R`` on all TRE files. Each run applies ``ape::read.tree()`` to read TREE file, and then computes stats for each tree. Resulting file is compressed with ``gzip``.
* Uses auxiliary functions in script ``treeStatFunctions.R``. 
* The stats: [ some stats are coded with specific settings for the input set ]
  * TMRCAH_all - the TMRCAH of the entire sample set (TMRCA of 1/2 of all samples). 
  * RT12 - the smallest TMRCA of *species_n*/2 samples (the 12 originates from the fact that in our data set *species_n*=24).
  * <species>_RTH - the RTH' of each species (one per species), which is the TMRCAH of that species divided by TMRCAH_all
  * <species>_enrich - the enrichment score of each species (one per species), determined using a hypergeometric test.


------------

Step 5 - Generate Stat BED files
------------

**Input:**

- STAT files  (.STAT.GZ ; generated in **Step 4 above)

**Output:**

- *argStats-windowed/<scaffold>.stat.bed*   -   (STAT BED file; one per scaffold/chromosome)

  [ format for line identical to STAT file (see above) with one line per window and start and end indices of window specified ] 

* *infoTables/scaffold-window-coverage.txt* - (table containing the ranges covered for every scaffold)

**Steps:**

``==> Rscript make_all_window_stats.R``

- **Setup:**  set the length of non-overlapping windows (default=20,000; should be consistent with setting in **Step 1**)
- Runs script ``window_stats.R`` on all scaffolds with specified window length.
- Writes repost on window coverage in *infoTables/scaffold-window-coverage.txt*. In some scaffolds have missing ARG blocks, this will be apparent in report.

``==> bash scripts/check_missing_blocks.sh``  (optional)

* **Setup:** set the path to a file containing a table with lengths of all scaffolds (to check for missing blocks); set various length parameters based on trimming length and window length (see guidelines in comments).
* Prints a list of missing blocks in scaffolds. If doesn't print anything, then all scaffolds are fully covered.


------------