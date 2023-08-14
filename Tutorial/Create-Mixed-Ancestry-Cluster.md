## Navigation
### Go back to <a href = "/Tutorial/Home.md">Wiki</a>
<hr>

# Create Mixed Ancestry Cluster

## This will create a large combined cluster with multiple ancestries but keeps ratios intact

**Must run on dev1-dev4**

Command: `$ ./Scripts/collapsing_lclust.sh clusteringMixedAncestry`

```
Usage: 
  lclust_Flash_combo_cluster.R --resolution_var=<resolution_var> [--min_sample=<min_sample>]  [--case_group=<case_group>] [--max_ratio=<max_ratio>]
  
  Options:
  -h --help
  --resolution_var=<resolution_var> resolution for clustering. Typically 0.1-0.4
  --case_group=<case_group> label for analysis. example might "picu" or "IPF" [default: temp_case]
  --min_sample=<min_sample> minimum size of case/control allowed to include cluster into combined cluster [default: 0]
  --max_ratio=<max_ratio> desired case/control ratio. Will downsample controls to achieve ratio. If blank, program will find the max value in among clusters and then choose that value. If users sets this, will ignore clusters that are greater than user value and downsample controsl for other clusters [default: NA] 
```
### YAML Variables Used:
```
USER_VARIABLE:
  case_group: project_name
COLLAPSING_VARIABLE:
  resolution: "0_2"
  min_sample: 20
  max_ratio: NA
  mixed_ancestry_cluster: NA
```
* case_group: project name, added to output file(s) name
* resolution: used to identify correct files
* min_sample: minimum size of case/control allowed to include cluster into combined cluster 
* max_ratio: desired case/control ratio
  * enter max_ratio in decimal form
* mixed_ancestry_cluster: needs to be set to a unique cluter number. Usually nubmer of clusters + 1

### Output:
* /Results/KinshipFlashPCA
  * txt file: a large combined cluster with multiple ancestries
