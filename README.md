# Relapse to cocaine-seeking is regulated by medial habenula Nr4a2 in mice
Childs & Morabito et al. 2024 (Cell Reports, in press)

This Repository contains the code used for data processing and analysis in our manuscript
titled "Relapse to cocaine-seeking is regulated by medial habenula Nr4a2 in mice". In this study
we used single-nucleus RNA-seq (snRNA-seq) to profile the habenula in four different groups of mice to
study the molecular changes following a cocaine reinstatement behavioral experiment and manipulation
of the transcription factor Nr4a2. Here we list the major sections of the paper and provide links to the data
analysis steps in each part.

## Code to program the cocaine reinstatement behavioral experiment

We supply the script that was used during the cocaine reinstatement behavioral experiment.

* [Cocaine reinstatement experiment script (.mpc file)](other/behavioral_test.mpc)

## Processing sequencing data and quantifying gene expression

snRNA-seq was performed using the 10X Genomics kit, and we used CellRanger to quantify gene expression 
from the raw sequencing reads. After running CellRanger, ambient RNA was removed using cellbender. 

* [Data preprocessing scripts](preprocessing/)

## Clustering analysis (Figure 2)

In our snRNA-seq dataset we had four different groups of mice: NURR2C and GFP mice that were
behaviorally experienced and behaviorally naive. The behaviorally naive and behaviorally experienced 
groups were analyzed separately before performing an integrated analysis. 

* [Behaviorally naive clustering analysis (Jupyter Notebook)](clustering/scanpy_processing_naive.ipynb)
* [Behaviorally experienced clustering analysis (Jupyter Notebook)](clustering/scanpy_processing_experienced.ipynb)
* [Integrated clustering analysis (Jupyter Notebook)](clustering/integrated_clustering.ipynb)
* [Plotting clustering results and QC metrics (R markdown)](clustering/plotting_for_paper.Rmd)
* [Cluster marker gene script (bash script)](differential_expression/cluster_markers.sub)
* [Cell type marker gene script (bash script)](differential_expression/celltype_markers.sub)
* [Misc. plotting utility script (R script)](other/plotting_utils.R)

To compare our snRNA-seq with a [previously published dataset of the mouse habenula](https://pubmed.ncbi.nlm.nih.gov/32272058/), 
we performed an additional integrated analysis.

* [Integration with Hashikawa et al (R Markdown)](clustering/integrate_hashikawa.Rmd)

Due to our experimental manipulation of Nr4a2 in the habenula, we were interested in investigating 
snRNA-seq read pileup at the Nr4a2 locus. For this analysis, we processed the sequencing data so that it could 
be viewed on the genome browser. 

* [Script to make genome browser trackhubs (R markdown)](trackhubs/make_trackhubs.Rmd)
* [Alternative splicing analysis with Swan (Jupyter Notebook)](other/swan.ipynb)
* [Additional helper scripts are in this directory](trackhubs/)

## Transcription factor (TF) network analysis (Figure 3)

As part of this study, we developed a custom strategy for TF network analysis. These functions have been formally added to the 
[hdWGCNA R package](https://smorabit.github.io/hdWGCNA/), but we provide the code as is here from before these functions
were incuded in hdWGCNA.

* [TF network analysis (R markdown)](TF_nets/TF_analysis.Rmd)
* [Individual functions to construct the TF nets (R script)](TF_nets/TF_functions.R)
* [Script to run TF net construction (R script)](TF_nets/parallel_TFnets.R)
* [Script to launch jobs to build separate TF nets on the HPC (bash script)](TF_nets/run_ConstructTFNets.sub)

## Differential expression analysis (Figures 4 and 5)

We performed differential expression analysis in each cell type and cell cluster to compare gene expression signatures
between behaviorally naive (Figure 4) and behaviorally experienced (Figure 5) NURR2C and GFP mice. We also compared 
these DEGs to each other (Figure 5). Note that the results plotting R markdown file contains plotting code 
for Figures 4 and 5, which summarize the DEG results for the naive and experienced mice.

* [DEG Results plotting code (R markdown)](differential_expression/DEG_plotting.Rmd)
* [Script to run DEG tests in parallel](other/parallel_DEGs.R)

Behaviorally naive (Figure 4)
* [Cluster DEG script (bash script)](differential_expression/cluster_Naive_Nurr2c_vs_GFP.sub)
* [Cell type DEG script (bash script)](differential_expression/celltype_Naive_Nurr2c_vs_GFP.sub)

Behaviorally experienced (Figure 5)
* [Cluster DEG script (bash script)](differential_expression/cluster_Nurr2c_vs_GFP.sub)
* [Cell type DEG script (bash script)](differential_expression/celltype_Nurr2c_vs_GFP.sub)

We also performed relative likelihood analysis with [MELD](https://github.com/KrishnaswamyLab/MELD)
in the behaviorally experienced mice to quantify the transcriptome-wide perturbation effects of Nr4a2 
manipulation.
* [MELD Relative likelihood analysis script (Jupyter Notebook)](other/run_MELD_revision.ipynb)

## Co-expression network analysis in medial habenula neurons (Figure 6)

We used [hdWGCNA](https://github.com/smorabit/hdWGCNA/) to perform co-expression network analysis to 
specifically study systems-level transcriptome changes in the medial habenula neurons, 
and to link TF regulatory networks to co-expression modules. 

* [Co-expression network analysis with hdWGCNA](hdWGCNA/hdWGCNA.Rmd)
