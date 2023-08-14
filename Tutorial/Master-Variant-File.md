## Navigation
### Go back to <a href = "/Tutorial/Home.md">Wiki</a>
<hr>

# Master Genotype File

* The next step is to create a master variant file for each cluster. This file will subsequently be filtered to create models so this file should be very inclusive of any model you may want to create in the future.

**Must be run on qs1**

Command: `./Scripts/collapsing_lclust.sh createMasterGeno`

### Input:
* Uses clusters found at /Results/Coverage

### Output:
* Each cluster gets it own folder, found at /Results/Collapsing
  * txt files of ATAV logs
  * csv files of existing samples and genotypes