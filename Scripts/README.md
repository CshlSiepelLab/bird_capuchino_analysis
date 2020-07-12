This directory contains analysis and plotting scripts for reproducing the analysis in the paper "Genomic islands of differentiation in a rapid avian radiation have been driven by recent selective sweeps" by Hussein A. Hejase, Ayelet Salman-Minkov, Leonardo Campagna, Melissa J. Hubisz, Irby J. Lovette, Ilan Gronau, and Adam Siepel (2020).

<h4>Our pipeline assumes you've downloaded the following scripts (used in the S/HIC software [1]):</h4>
1. splitMsOutputIntoWindows.py - https://github.com/kr-colab/shIC <br />
2. combineWinsIntoFeatureVec.py - https://github.com/kr-colab/shIC <br />
3. niceStats - https://github.com/kr-colab/shIC <br />

<h4>Example run of pipeline:</h4>
Let's assume you've run hard sweep simulations (hardp2.slim; 10,000 replicates) in species #1  and stored the output in directory hardp2. You can run the pipeline using the following command: sh pipeline.sh hard 10000 hardp2 1 5 48. The arguments are: type of simulation = hard; number of replicates = 10,000; directory = hardp2; id = 1 (or any random integer); number of windows = 5; sample size = 48.

<h4>References:</h4>
[1] Schrider, Daniel R., and Andrew D. Kern. "S/HIC: robust identification of soft and hard sweeps using machine learning." PLoS genetics 12.3 (2016): e1005928.
