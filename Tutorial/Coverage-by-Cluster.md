## Navigation
### Go back to <a href = "/Tutorial/Home.md">Wiki</a>
<hr>

# Coverage by Cluster
* We next run coverage on all of the sample files created in the clustering step.
* Coverage accomplishes the following: To alleviate issues coming from differential coverage between the case and control cohorts, we next prune out noisy nucleotide sites among the consensus coding sequence (CCDS release 20 [I think]) public transcript or their 2 base pair intronic extensions:
* We have full clusters, male only, female only 
* Note: generally, you need only run the first coverage loop. Only run the sex-specific loops if you intend to run sex-specific analyses. (TODO: discuss with Gundula circumstances for running sex-specific analyses) 

**Must be run on qs1**

Command: `./Scripts/collapsing_lclust.sh clusterCoverage`

### Input:
* Uses cluster files found at /Results/KinshipFlashPCA

### Output:
* Each cluster gets it own folder found at /Results/Coverage
  * txt files of ATAV logs
  * csv files of coverage summary and sample summary
  * txt files of existing samples and exon site list