Scripts for using ARGs sampled by ARGweaver to analyze patterns of species differentiation - summarize_stats
=======

This directory contains scripts that take windowed stat summary files (.STAT.BED) and intermediate TREE files and analyze patterns of species differentiation along the genome. This implements **tests 1**, **2**, and **3**, used in **Table 1** of our manuscript. Additional scripts are given for plotting summaries and trees. The pipeline summary:

**Input:**   Summary files computed in smc_to_stats + table with Fst peaks.

**Step 1:**  Compute control windows and blocks.

**Step 2:**  Summarize enrichment scores, RTH', and RT12 in Fst peaks and control windows / blocks (**test1** and **test 2**)

**Step 3:**  Summarize cross-coalescence (CC) times in Fst peaks, flanking regions, and control blocks (part of **test3**)

**Step 4:**  Plotting and visualization (not inherent part of analysis)


------------

Setup
------------

- Run the pipeline described in *smc_to_stats/README.md*  to generate STAT BED files and ARG Tree files with summaries in directories *argTrees/* and *argStats-windowed/* under the current directory. The pipeline also generates the file *scaffold-window-coverage.txt*  in the *infoTables/* subdirectory, which is used here.
- Prepare table with Fst peaks (*infoTables/FST-table-base.txt*). This file contains a table with a row per Fst peak (25 in our analysis) with the following columns: *peakID* (id for peak to be used in tables)	*Scaffold*  *startPos* *endPos*  (coordinates of peak). Additional columns are optional and not used by script (such as *Chromosome* indicating mapping to species with complete genome build). See example in *sampleFiles/sporophila-FST-table-base.txt* in root directory. 
- Set bash variable `$scriptDir` the path of **this** directory (the one containing these scripts).
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

``==> bash $scriptDir/create_control_region_list.sh``

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

``==> Rscript $scriptDir/summarize_stats_fst_n_control.R``

- **Setup:**  set the desired empirical p-value for species enrichment (**test 1**) and RTH' (**test 2**) (default=0.0001), and a separate (more relaxed) p-value for RT12 (default=0.001). You can also set the fraction of values that will be printed from the tails of the statistic distributions in the report file. This value should be greater than the p-value (default=0.0005). Additinally, if you wish only to analyze the Fst peaks or the control blocks or windows, you can toggle the three flags in the top of the script file from **TRUE** to **FALSE**.
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

``==> bash $scriptDir/cross_coal_fst_or_control.sh \``
        ``infoTables/FST-table-with-center.txt  infoTables/cross_coals_fst_peaks.tsv   fst  2``

- **Setup:** set the number of parallel jobs that can run on your machine (default=20); you may also set the start and end rows that you which to analyze in the Fst table file (default start = $4=2; default end = end-of-file). If the flag `doFSTstatFile` is set to `true`, then this script will also create a windowed stats file with average cross coalescence times across windows. This is not used in analysis, but could be used in exploratory plots (see **Step 4** below).
- For each row processed in the Fst table file, the script runs script ``cross_coal_per_region.R`` with appropriate arguments (see below for description of this script and arguments). This produces a temporary one-row table file for this row.
- Temporary table files for individual rows are appended to *infoTables/cross_coals_fst_peaks.tsv* every time a parallel batch ends.
- If the flag `doFSTstatFile` is set to `true`, then the script ``window_cross_coal.R`` is executed to produce a stats file with average values per 20 kb window. This file can be used to plot values in peak and in flanking 200 kb. It is not used in analysis.

``==> bash $scriptDir/cross_coal_fst_or_control.sh \``
        ``infoTables/control-blocks.txt  infoTables/cross_coals_control_blocks.tsv   control  1``

- **Setup:** set the number of parallel jobs that can run on your machine (default=20); you may also set the start and end rows that you which to analyze in the control block file (default start = $4=1; default end = end-of-file). Script makes hard assumption that control blocks are 500 kb long.
- For each row processed in the control block file, the script runs script ``cross_coal_per_region.R`` with appropriate arguments (see below for description of this script and arguments). This produces a temporary one-row table file for this row.
- Temporary table files for individual rows are appended to *infoTables/cross_coals_control_blocks.tsv* every time a parallel batch ends.

``==> Rscript $scriptDir/cross_coal_control_summary.R``

- **Setup:** set p-value used for determining the threshold for significant elevation / reduction in quantile differences (default=0.01). By default, the script compute threshold for elevation, and if you wish to analyze lower tail of quantile differences, toggle *upperOrLowerTail* from **TRUE** to **FALSE**.
- Analyzes the distribution of quantile differences in *infoTables/cross_coals_control_blocks.tsv*, and reports the upper/lower tail of values and the resulting threshold. Reports a global threshold (across all pairwise comparisons) as well as pairwise-specific thresholds. We use the global threshold in determning significant elevation / reduction in **Table S5** and in **test 3**. All reports are directoed to the standard output.

**Auxiliary scripts:**

````
cross_coal_per_region.R
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
* Once the region of interest and flanking regions are determined, the script reads trees from the relevant TREE files in *treeDir*, and for every species-pair, it applies the function ``getCrossCoalDistribution()`` from the auxiliary script ``crossCoalFunctions.R`` in ``smc_to_stats`` drectory to get the distribution of normalized recent CC times. Then, the quantile difference between these two distributions is computed and printed to the output table.

````
crossCoalFunctions.R 
````

* Contains various utility functions for computing the youngest cross coalescence in a given tree, the distribution of normalized young CC times in a given list of trees, and the quantile difference between distributions. This script is located in the ``smc_to_stats`` directory.



------------

Step 4 - Plotting and visualization
------------

#### A. Plot average statics across a given region

**Input:**

- *argStats-windowed/<scaffold>.stat.bed*   -   (STAT BED file; one per scaffold/chromosome; grnerated in **Step 5** of smc_to_stats pipleine)
- *argStats-windowed/<scaffold>.<range>.crosscoal.bed*   -   (STAT BED file with average normalized cross coalescence times in genomic windows; needed for plotting cross coalescence times; see below)

**Output:**

- PDF file with plots describing tree-based statistics in region of interest.

``==> Rscript $scriptDir/plot_segment_stats_per_region.R \``
        ``./argStats-windowed <plotFilePath.pdf> <scaffold> <startPos> <endPos> <regionStart> <regionEnd>``

For example, to plot statistics in the region surrounding the Fst peak on Contig252 (200 kb flanking in each direction) into a PDF file named *peak_stats_Contig252.pdf* in directory *plots/*, we apply:

``==> Rscript $scriptDir/plot_segment_stats_per_region.R \``
        ``./argStats-windowed plots/peak_stats_Contig252.pdf Contig252 220000 710000 420000 510000``

- **Setup:** 
  - This script requires installing additional libraries: ``ggplot2``,  ``cowplot``, and  ``reshape2``.
  - If you wish for cross coalescence times to be plotted, then make sure that the average values have been computed and saved in the appropriate file in *./argStats-windowed/*. This is done for the Fst peaks when you run the script ``cross_coal_fst_or_control.sh`` with the flag `doFSTstatFile` turned on (see **Step 3** above). For other regions, you need to execute the script ``window_cross_coal.R``  with appropriate parameters.
  - Set the type of plots you wish to be plotted: species enrichment / RTH / normalized cross coalescence times. You may turn each type off or on, and they will be plotted in the above order.
  - Set significance thresholds for each type, if you wish for a horizontal dashed line to indicate it (set to -1 if you do not wish a threshold to be plotted).
  - Set the dimensions of each plot (panel). 
- The script generates a single PDF file with all relevant plots. The plots are determined by the flags set in the script and are presented from top to bottom in this order: species enrichment --> RTH --> normalized corss coalescence times. Colors for species enrichment and RTH are determined by the colors specified in the species key table (*infoTables/individual-species-key-modified.txt*). The RT12 statistic is plotted in gray and the pairwise cross coalescences are plotted in colors from the rainbow palette.

#### B. Plot trees with extreme values for species statistics

**Input:**

- *argStats/\*.tre.gz*   -   (TREE files  (.TRE.GZ)  generated in **Step 3** of ``smc_to_stats`` pipeline)
- *argStats/\*.stat.gz*   -   (STAT files  (.STAT.GZ)  generated in **Step 4** of ``smc_to_stats`` pipeline)

**Output:**

- PDF file with demonstrative local trees extracted from the infrerred ARGs (see details below).

``==> Rscript $scriptDir/plot_local_tree.R <scaffold> <startPos> <endPos> <outFile>``

For example, to plot trees that show various species separations in the Fst peak on Contig252 into a PDF file named *demo_trees_Contig252.pdf* in directory *plots/*, we apply:

``==> Rscript $scriptDir/plot_local_tree.R  Contig252 420000 510000 plots/demo_trees_Contig252.pdf``

- **Setup:** this script requires installing [the ``argweaver`` R library](https://github.com/CshlSiepelLab/argweaver/tree/master/R), for the main procedure that plots the tree. If you wish for a subset of individuals to be plotted, you can set the ``inds_to_keep`` variable to a list of appropriate haploid sample names.
- The script generates a single PDF file with all relevant plots.  Two trees are plotted per species - one with the highest enrichment score, and one with the lowest RTH. Each tree is plotted in a separate page and the species are ordered in decreasing order of the highest enrichment score achieved for that species in the specified genomic region. The colors of tips in each tree are colored based on the color scheme defined for the species in the *infoTables/individual-species-key-modified.txt* file. A short caption is written in the bottom of every page indicating the position of the tree, the statistic it demonstrates and the value of that statistic.


------------
