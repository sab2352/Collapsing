## Navigation
### Go back to <a href = "/Tutorial/Home.md">Wiki</a>
<hr>

# Track Sample Removal 

* Creates a table with cases and controls at each step

**Must be run on dev1-dev4**

Command: `$ ./Scripts/collapsing_lclust.sh sampleTracking`

```
Usage: 
  table_sample_counts.R --resolution_var=<resolution_var> --min_sample=<min_sample> [--case_group=<case_group>]
  
  Options:
  -h --help
  --min_sample=<min_sample> minimum number of cases/controls of included clusters
  --resolution_var=<resolution_var> resolution for clustering.
  --case_group=<case_group> label for analysis. example might "picu" or "IPF" [default: temp_case]
```
If you want to view help message in terminal, run command: `$ ./Scripts/collapsing_lclust.sh sampleTracking --help`
### Yaml Variables Used:
```
USER_VARIABLE:
  case_group: project_name
COLLAPSING_VARIABLE:
  resolution: "0_2"
  min_sample: 20
```
* case_group: project name, added to output file(s) name
* resolution: used to identify correct files
* min_sample: minimum number of cases/controls of included clusters

### Input:
  * Output form previous steps (Creating Cohort - Creating Clusters)

### Output:
  * file with cases and controls at each step 