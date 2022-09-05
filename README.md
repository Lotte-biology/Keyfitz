# Keyfitz
Code for figures from paper "A note on discretising Keyfitz entropy"
The matlab file create_figure1.m can be run to create Figure 1, and to calculate the entropies in Table 1. Change Deltat to 0.01 or 1 to get the two columns of Table 1. 

The R files (mainEntropyAnalysis_ver1.R; stageConversion_ver1.R; KeyfitzVisuals_ver1.R) can be used to create Figure 2:
- stageConversion_ver1.R has code to convert stage based matrix models to age-based matrix models
- mainEntropyAnalysis_ver1.R calculates the original and the proposed Keyfitz metrics for matrix models from COMPADRE and COMADRE
- KeyfitzVisuals_ver1.R creates plots
