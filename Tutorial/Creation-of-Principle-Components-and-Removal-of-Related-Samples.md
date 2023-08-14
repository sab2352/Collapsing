## Navigation
### Go back to <a href = "/Tutorial/Home.md">Wiki</a>
<hr>

# Kinship and FlashPCA

At this stage, we run an ATAV command that accomplishes two objects.

* We need to remove all related samples from our cohort. We use KING to ensure that only unrelated (up to third-degree) individuals are retained in the sample list.

* We need to create principal components based on the SNPs associated with geographic ancestry. This will be used for clustering in the next step.
Note: by default --flashpca will not perform PLINK nearest-neighbor pruning anymore use --flashpca-plink-pruning to trigger running plink pruning step if you want and add --flashpca-num-nearest-neighbor 200 

**Must be done on qs1**

Command: `./Scripts/collapsing_lclust.sh CohortPCA`

### Input:
* Automatically detects the file created in the SampleFileCheck step found in /Results/SampleFileCheck/*_existing.sample.txt
* Detects kinship-relatedness-threshold from the yaml file. The default value used by ATAV is 0.0884, or third degree relatedness. This option removes samples if they have a relatedness score above this threshold (ie. If the threshold is 0.0884, then samples with a kinship relatedness score greater than 0.0884 will be removed).

### Output:
* /Results/KinshipFlashPCA$dir
  * txt files of kinship results and pruned samples
 

### Kinship Ranges:
* [0.177, 0.354] : 1st-degree
* [0.0884, 0.177] : 2nd-degree
* [0.0442, 0.0884] : 3rd-degree
