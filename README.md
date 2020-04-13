# Module Analysis Scripts
## Author: Pasan Fernando

This repository contains the scripts used for the project titled: "Studying the changes in gene modules associated with fin to limb transition using quality-enhanced gene networks"

The scripts are arranged into following folders

**Gene_extractor**: A script to extract genes annotated directly to a given anatomical entity. Used to prepare the genes with original annotations.

**Anatomy_part_extractor**: A script to extract genes annotated to the parts of a given anatomical entity. Used to prepare the genes with original annotations.

**Precision_threshold_variation**: Using the trial and error method to select a precision threshold for module extraction.

**Threshold_mdetection**: Selects candidate genes for a given anatomical entity and extracts the module based on a given precision threshold score.

**Module_net_prediction**: Predict scores for each unknown gene based on Hishigaki network-based candidate gene prediction method. Moreover, contains a script for evaluating the network-based candidate gene prediction performance using leave-one-out cross-validation to generate ROC and precision-recall curves

**Module_comparison**: Compares two given modules (e.g., pelvic fin vs hindlimb). Also, contains a script to distinguish whether the compared genes are directly annotated or annotated to parts or precursors.

**Weighted_degree_comparison**: Performs comparisons of weighted degrees for predicted vs genes with original annotations and module-specific vs conserved genes

**Uberon_enrichment_analysis**: Perform enrichment analysis using Uberon terms and Fisher's exact test for a given gene list

