# ATAUC Supplementary Materials

This repository provides supplementary data and R scripts associated with the paper:

**Quantifying Ruminal Health: A Statistical Review and Application of Area and Time Under the Curve in Animal Science**  
Luis O. Tedeschi

## Description

The ATAUC (Area and Time Above and Under the Curve) software is a proprietary tool developed for analyzing time-series data related to ruminal pH. This repository **does not contain the ATAUC program or source code**, but includes:

- Example datasets (`data.xlsx`)
- The `getauc` R function used to calculate AUC, AAC, TUC, and TAC
- Selected ATAUC output files to illustrate how the software works with real data
- 
These materials are sufficient to **reproduce the statistical analyses and visualizations** reported in the manuscript.

## Files

- `Data.xlsx`: Excel file used in the manuscript, with example data (see Figure 1).
- `getauc.R`: R function to calculate area and time metrics with spline interpolation.
- `atauc_output.RAR`: Example output file generated by ATAUC (read-only illustration).

## How to Use

1. Install R (version 4.4.1 or newer).
2. Install required R packages:
   ```R
   install.packages(c("MESS", "rootSolve"))
