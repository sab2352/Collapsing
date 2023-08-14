## Navigation
### Go back to <a href = "/Tutorial/Home.md">Wiki</a>
<hr>

# Creating Clusters

* Using the principal components created in the prior step, the following code uses Louvain clustering to create geographic ancestry based clustering. 
* The --resolution_var can be changed to change the size of the clusters. 
* Optimal clusters balance adequate cluster size while ensuring clusters are as uniform as possible.

**Runs on QS1**

Command: `./Scripts/collapsing_lclust.sh clustering`

```
 Usage: 
  lclust_Flash.R [--case_group=<case_group>] [--resolution_var=<resolution_var>]
  
  Options:
  -h --help
  --case_group=<case_group> label for analysis. example might "picu" or "IPF" [default: temp_case]
  --resolution_var=<resolution_var> resolution for clustering. Typically 0.1-0.4 [default: 0.3]
```
### Yaml Variables Used:
```
USER_VARIABLE:
  case_group: project_name
COLLAPSING_VARIABLE:
  resolution: "0_2"
```
* case_group: project name, added to output file(s) name
* resolution: resolution variable used for clustering
  * Can be changed to change the size of the clusters. Most likely will need to change resolution in order to create optimal clusters. Must use underscore when changing resolution (0.05 = 0_05)

### Output:
* /Results/Plots
  * txt files of the cluster sizes and ratios
  * png files of cluster maps