# Welcome to the ClusteredCollapsing wiki!

This wiki is a step-by-step guide for running analyses. For granular detail regarding ATAV functions, please refer to the ATAV wiki (https://redmine.igm.cumc.columbia.edu/projects/atav/wiki)

Setup
* [Repository Setup](https://github.com/igm-team/ClusteredCollapsing/wiki/Repository-Setup) - Instructions on cloning/downloading GitHub repo. Best located on /projects when possible. Each cohort will require its own copy of the repo.

* [Yaml Setup](https://github.com/igm-team/ClusteredCollapsing/wiki/Yaml-Setup) - The code uses a configuration file (.yaml). This is the only part of the code that is edited by the user. Once done, you're ready to begin!

Analysis
* [cohortSelection](https://github.com/igm-team/ClusteredCollapsing/wiki/Creating-Cohort) - In this step, you create your ped file. This is a file of cases and controls. 
* [SampleFileCheck](https://github.com/igm-team/ClusteredCollapsing/wiki/SampleFileCheck) - Once the cohort is created, this quick steps makes sure that each row is formatted correctly for ATAV
* [initialCoverage](https://github.com/igm-team/ClusteredCollapsing/wiki/Coverage-of-CaseControl-Cohort) - Next, coverage is performed on the case/control cohort. The output for this step is used in the step when removing of related individuals (given a choice, the program keep the better covered individual)
* [CohortPCA](https://github.com/igm-team/ClusteredCollapsing/wiki/Creation-of-Principle-Components-and-Removal-of-Related-Samples) - In this step, we create ancestry based principle components and remove related individuals. it uses the coverage data from initialCoverage to keep the best covered samples.
* [clustering](https://github.com/igm-team/ClusteredCollapsing/wiki/Create-Clusters) - This step uses the PCAs generated in the prior step to create ancestry based case/control clusters.
* [clusteringMixedAncestry](https://github.com/igm-team/ClusteredCollapsing/wiki/Create-Mixed-Ancestry-Cluster) - For those interested in analyzing a single cluster with all ancestries, this step facilitates that. It can be used in Regenie analyses.

The following steps are then performed in each cluster
* [pcaCluster](https://github.com/igm-team/ClusteredCollapsing/wiki/pcaCluster) - Re-runs ancestry PCA on each cluster. In doing so, excludes ancestry outliers. This step can be skipped if desired.
* [clusterCoverage](https://github.com/igm-team/ClusteredCollapsing/wiki/Coverage-by-Cluster) - Coverage harmonization is performed within each cluster. This step can be skipped if desired. 
* [createMasterGeno](https://github.com/igm-team/ClusteredCollapsing/wiki/Master-Variant-File) - This step creates large variant files for each cluster. Subsequent models will be drawn from these files so only one trip to the database is necessary.
* [models](https://github.com/igm-team/ClusteredCollapsing/wiki/Create-Models) - This step filters the variants in the master genotype files created in the preceding step in order to create specific models

The following steps aggregate the analysis across included clusters
* [cmh](https://github.com/igm-team/ClusteredCollapsing/wiki/Gene-Based-CMH)
* [CMH Permutations](https://github.com/igm-team/ClusteredCollapsing/wiki/CMH-Permutations)
* [Create QQ Plots](https://github.com/igm-team/ClusteredCollapsing/wiki/Create-QQ-Plots)

The following steps can be performed for the mixed ancestry cluster
* pcaMixedAncestry
* clusterCoverageMixedAncestry
* createMasterGenoMixedAncestry
* modelsMixedAncestry


The following are optional steps usually necessary for methods or supplemental data for collapsing papers/presentations
* [Create Cluster Map](https://github.com/igm-team/ClusteredCollapsing/wiki/Create-Cluster-Map)
* [Track Sample Removal](https://github.com/igm-team/ClusteredCollapsing/wiki/Track-Sample-Removal)

Other analyses
* [Run Regenie](https://github.com/igm-team/ClusteredCollapsing/wiki/Run-Regenie) - Run regenie 3.1.3 on your data
* [Run NHC](https://github.com/igm-team/ClusteredCollapsing/wiki/Run-NHC)
* [Forest Plot](https://github.com/igm-team/ClusteredCollapsing/wiki/Forest-Plot) - Generate forest plot if provided a csv with ORs, CIs, P values, etc.
* [findThreshold](https://github.com/igm-team/ClusteredCollapsing/wiki/findThreshold) - Detects optimal threshold for gene set analysis. Currently only working for LOEUF.

Steps without documentation
* clusterCoverageMale
* clusterCoverageFemale


[Frequently Asked Questions](https://github.com/igm-team/ClusteredCollapsing/wiki/Frequently-Asked-Questions)

[Troubleshooting](https://github.com/igm-team/ClusteredCollapsing/wiki/Troubleshooting)