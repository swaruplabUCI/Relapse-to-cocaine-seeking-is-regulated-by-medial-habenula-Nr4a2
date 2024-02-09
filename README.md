# Relapse to cocaine-seeking is regulated by medial habenula Nr4a2 in mice
Childs & Morabito et al. 2024 (Cell Reports, in press)

This Repository contains the code used for data processing and analysis in our manuscript
titled "Relapse to cocaine-seeking is regulated by medial habenula Nr4a2 in mice". In this study
we used single-nucleus RNA-seq (snRNA-seq) to profile the habenula in four different groups of mice to
study the molecular changes following a cocaine reinstatement behavioral experiment and manipulation
of the transcription factor Nr4a2. Here we list the major sections of the paper and provide links to the data
analysis steps in each part.


## Processing sequencing data and quantifying gene expression

snRNA-seq was performed using the 10X Genomics kit, and we used CellRanger to quantify gene expression 
from the raw sequencing reads. After running CellRanger, ambient RNA was removed using cellbender. 

* [preprocessing scripts](preprocessing/)

## Clustering analysis (Figure 2)

In our snRNA-seq dataset we had four different groups of mice: NURR2C and GFP mice that were
behaviorally experienced and behaviorally naive. The behaviorally naive and behaviorally experienced 
groups were analyzed separately before performing an integrated analysis. 

* [Behaviorally naive clustering analysis](clustering/scanpy_processing_april2023.ipynb)

## Transcription factor (TF) network analysis (Figure 3)

## Differential expression analysis in behaviorally naive mice (Figure 4)

## Differential expression analysis in behaviorally experienced mice (Figure 5)

## Co-expression network analysis in medial habenula neurons (Figure 6)

