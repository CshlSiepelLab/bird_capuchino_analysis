Scripts for using ARGs sampled by ARGweaver to analyze patterns of species differentiation - summarize_stats
=======

This directory contains scripts that take windowed stat summary files (.STAT.BED) and intermediate TREE files and analyze patterns of species differentiation along the genome. This implements **tests 1**, **2**, and **3**, used in **Table 1** of our manuscript. Additional scripts are given for plotting summaries and trees. The pipeline summary:

**Input:**   Summary files computed in smc_to_stats + table with Fst peaks.

**Step 1:**  Compute control windows and blocks.

**Step 2:**  Summarize enrichment scores, RTH', and RT12 in Fst peaks and control windows / blocks (**test1** and **test 2**)

**Step 3:**  Summarize cross-coalescence (CC) times in Fst peaks, flanking regions, and control blocks ( part of **test3**)


------------

Setup
------------

- Run the pipeline described in *smc_to_stats/README.md*  to generate STAT BED files and ARG Tree files with summaries
- Prepare table with Fst peaks (*infoTables/FST-table-base.txt*). This file contains a table with a row per Fst peak (25 in our analysis) with the following columns: *peakID* (id for peak to be used in tables)	*Scaffold*  *startPos* *endPos*  (coordinates of peak). Additional columns are optional and not used by script (such as *Chromosome* indicating mapping to species with complete genome build). See example in *sampleFiles/sporophila-FST-table-base.txt* in root directory. 
- Follow the specific setup instructions of each step in the pipeline below


------------

Step 1 - Determine control blocks
------------

**Input:**

- *infoTables/FST-table-base.txt* - (Fst peak source table)
- *infoTables/scaffold-window-coverage.txt* - (table containing the ranges covered for every scaffold; grnerated in **Step 5** of smc_to_stats pipleine)

**Output:**

- *infoTables/control-scaffolds.txt*   -   (line per scaffold that does not contain an Fst peak; columns: <chrom> <startPos> <endPos>)  
- *infoTables/control-blocks.txt*   -   (same format but line per control block) 

**Steps:**

``==> bash create_control_region_list.sh``

- **Setup:**  set the length of control blocks (default=500,000)
- Uses  Fst peak table and scaffold coverage file to determine a list of scaffolds that do not contain Fst peaks, and then tiles these scaffolds with blocks of a given length.


------------

Step 2 - Summarize enrichment scores, RTH', and RT12 in Fst peaks and control windows / blocks
------------

**Input:**

- *infoTables/FST-table-base.txt* - (Fst peak source table)
- *argStats-windowed/<scaffold>.stat.bed*   -   (STAT BED file; one per scaffold/chromosome; grnerated in **Step 5** of smc_to_stats pipleine)
- *infoTables/control-scaffolds.txt*   -   (generated in **Step 1** above)  
- *infoTables/control-blocks.txt*   -   (generated in **Step 1** above) 

**Output:**

- *infoTables/fst_peak_arg_stats.tsv*   -   (tab-delimited table with maximum species enrichment and minimum species-RTH' and RT12 for every Fst peak; Columns: <scaffold> <startPos> <endPos> <stat1> <stat2> . . .)  ***These are values used in Tables S2 & S3***
- *infoTables/control_block_arg_stats.tsv*   -   (tab-delimited table with maximum species enrichment and minimum species-RTH' and RT12 for every control block; same format as above)
- *infoTables/control_regions_summary.txt*   -   (free text summary of stats measured in control scaffolds and blocks; see below).

**Steps:**

``==> Rscript summarize_stats_fst_n_control.R``

- **Setup:**  set the desired empirical p-value for species enrichment (**test 1**) and RTH' (**test 2**) (default=0.0001), and a separate (more relaxed) p-value for RT12 (default=0.001). You can also set the fractionof values that will be printed from the tails of the statistic distributions in the report file. This value should be greater than the p-value (default=0.0005). Additinally, if you wish only to analyze the Fst peaks or the control blocks or windows, you can toggle the three flags in the top of the script file from **TRUE** to **FALSE**.
- Run script with ``compute_FST_stats=TRUE`` to generate min/max stat table for Fst peaks
- Run script with ``compute_control_block_stats=TRUE`` to generate smin/max stat table for control blocks.
- Run script with ``compute_control_scaffold_stats=TRUE`` to compute distribution of stat values in all windows in control scaffolds and determine significance thresholds based on empirical p-values (see **Setup**). The report contains three parts:
  - The number of windows in control scaffolds and a list of tail values for all species enrichments, RTH' and RT12. 
  - A table of significance thresholds for each stat (thresholds used in **test1** and **test 2**)
  - (if ``compute_control_block_stats=TRUE``) A report of the number of control blocks that exceed the significance thresdholds determined above.


------------

Step 3 - Summarize cross-coalescence (CC) times in Fst peaks, flanking regions, and control blocks
------------

**Input:**

- *infoTables/FST-table-base.txt* - (Fst peak source table)
- *infoTables/individual-species-key-modified.txt* - to map haploid samples to species
- *argTreeFiles/<scaffold>.<startInd>-<endInd>.<MCMCiter>.tre.gz*   -   (ARG TREE files; grnerated in **Step 3** of smc_to_stats pipleine)
- *infoTables/control-blocks.txt*   -   (generated in **Step 1** above) 

**Output:**

- *infoTables/cross_coals_fst_peaks.tsv*   -   (tab-delimited table with quantile differences in CC times for every species pair and every Fst peak; Columns: <peakId> <pair1> <pair2> . . .)  ***These are values used in Table S5***
- *infoTables/cross_coals_control_blocks.tsv*   -   (tab-delimited table with quantile differences in CC times for every species pair and every control block; same format as above)
- *infoTables/cross_coals_control_blocks_summary.txt*   -   (free text summary of quantile differences measured in control blocks; see below).

**Steps:**

``infoTables/FST-table-with-center.txt``

* Prepare  ***infoTables/FST-table-with-center.txt***  - a version of the Fst peak table with a "center" specified for every peak. This should be a <= 200 kb region that is contained in the peak and contains the most species-enriched windows (the entire peak if it's shorter than 200 kb). This is done manually by copying  *infoTables/FST-table-base.txt* to file named *infoTables/FST-table-with-center.txt* and adding two columns (<peakCenterStart> and <peakCenterEnd>). See example in *sampleFiles/sporophila-FST-table-with-center.txt* in root directory. 

``==> bash cross_coal_fst_or_control.sh \``
        ``infoTables/FST-table-with-center.txt  infoTables/cross_coals_fst_peaks.tsv   fst  2``

- **Setup:** set the number of parallel jobs that can run on your machine (default=20); you may also set the start and end rows that you which to analyze in the Fst table file (default start = $4=2; default end = end-of-file).
- For each row processed in the Fst table file, the script runs script ``cross_coal_per_pregion.R`` with appropriate arguments (see below for description of this script and arguments). This produces a temporary one-row table file for this row.
- Temporary table files for individual rows are appended to *infoTables/cross_coals_fst_peaks.tsv* every time a parallel batch ends.

``==> bash cross_coal_fst_or_control.sh \``
        ``infoTables/control-blocks.txt  infoTables/cross_coals_control_blocks.tsv   control  1``

- **Setup:** set the number of parallel jobs that can run on your machine (default=20); you may also set the start and end rows that you which to analyze in the control block file (default start = $4=1; default end = end-of-file). Script makes hard assumption that control blocks are 500 kb long.
- For each row processed in the control block file, the script runs script ``cross_coal_per_pregion.R`` with appropriate arguments (see below for description of this script and arguments). This produces a temporary one-row table file for this row.
- Temporary table files for individual rows are appended to *infoTables/cross_coals_control_blocks.tsv* every time a parallel batch ends.

``==> Rscript cross_coal_control_summary.R``

- **Setup:** set p-value used for determining the threshold for significant elevation / reduction in quantile differences (default=0.01). By default, the script compute threshold for elevation, and if you wish to analyze lower tail of quantile differences, toggle *upperOrLowerTail* from **TRUE** to **FALSE**.
- Analyzes the distribution of quantile differences in *infoTables/cross_coals_control_blocks.tsv*, and reports the upper/lower tail of values and the resulting threshold. Reports a global threshold (across all pairwise comparisons) as well as pairwise-specific thresholds. We use the global threshold in determning significant elevation / reduction in **Table S5** and in **test 3**.

**Auxiliary scripts:**



````R
crossCoalFunctions.R
````

- Contains various utility functions for computing the youngest CCs in a given tree, the distribution of normalized young CC times in a given list of trees, and the quantile difference between distributions.



````
cross_coal_per_pregion.R
````

- Input arguments:
  - *treeDir*                 - directory where TREE files are found (*./argTreeFiles/*)
  - *regionID*               - an ID string for region (fst peak ID or <scaffold>.control_<i> for control block) 
  - *scaffold*                - scaffold ID
  - *leftFlankEnd*         - rightmost position of left flanking region (see below)
  - *rightFlankStart*    -  leftmost position of right flanking region (see below)
  - *startPos*                - leftmost positoin of region of interest
  - *endPos*                  - rightmost positoin of region of interest
  - *outFile*                   - name of (temporary) file where the one-line table is written

* Total length of flanking regions is set to 200 kb in code.
* For Fst peaks, *leftFlankEnd* and *rightFlankStart* are set according to the peak boundaries, and *startPos* and *endPos* are set based on the peak's center.
* For control blocks, *leftFlankEnd*=*startPos* is set 150 kb from start of block (100 kb left of center), and *rightFlankStart*=endPos* is set 350 kb from start of block (100 kb right of center). 
* The script computes the flanking regions by trying to take 100 kb on each side starting from *leftFlankEnd* and *rightFlankStart*. In the control blocks, we are guaranteed to have such flanking regions, and the first and last 50 kb in each block are not analyzed. If an Fst peak is near the edge of the scaffold and there are fewer than 100 kb on one side, we elongate the flanking region on the other side s.t. the total length is 200 kb. If the total number of bases flanking an Fst peak on the scaffold is less than 200 kb, the script prints a message and outputs all NA values to the output table.
* Once the region of interest and flanking regions are determined, the script reads trees from the relevant TREE files in *treeDir*, and for every species-pair, it applies the function ``getCrossCoalDistribution()`` from the auxiliary script ``crossCoalFunctions.R`` to get the distribution of normalized recent CC times. Then, the quantile difference between these two distributions is computed and printed to the output table.


------------
