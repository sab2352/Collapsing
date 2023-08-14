## Navigation
### Go back to <a href = "/Tutorial/Home.md">Wiki</a>
<hr>

# Runs Network-based Heterogeneity Clustering (NHC)

**Must Run on dev1-dev4**

Command: `./Scripts/collapsing_lclust.sh runNHC`

```
Usage: 
  run_NHC.R --cluster=<cluster> --resolution_var=<resolution_var> [--path_to_nhc_git=<path_to_nhc_git>] [--gene_list_case_path=<gene_list_case_path>] [--gene_list_ctrl_path=<gene_list_ctrl_path>]  [--min_sample=<min_sample>]  [--case_group=<case_group>] [--max_ratio=<max_ratio>] [--num_pcas=<num_pcas>] [--model=<model>] [--mixed_ancestry_logical=<mixed_ancestry_logical>]
  
  Options:
  -h --help
  --cluster=<cluster> define cluster to run NHC on
  --path_to_nhc_git=<path_to_nhc_git> define path to NHC git hub repo. See https://github.com/casanova-lab/NHC for NHC documentation [default: /nfs/labs/goldstein/jm4279/github_repos/NHC/]
  --model=<model> will run on a model if specified. If specified, no need for gene_list_case_path or gene_list_ctrl_path [default: NA]
  --gene_list_case_path=<--gene_list_case_path> path to list of cases and genes harboring QV. Best practice is for this list to live in a folder /Results/NHC [default: NA]
  --resolution_var=<resolution_var> resolution for clustering. Typically 0.1-0.4
  --gene_list_ctrl_path=<--gene_list_ctrl_path> path to list of controls and genes harboring QV. [default: NA]
  --case_group=<case_group> [default: temp_group]
  --min_sample=<min_sample> minimum size of case/control allowed to include cluster into combined cluster. Need to overide default if using combined cluster. Otherwise, will  [default: NA]
  --max_ratio=<max_ratio> this is added to allow loading of combined ancestry cluster [default: NA]
  --mixed_ancestry_logical=<mixed_ancestry_logical> [default: FALSE]
```
If you want to view help message in terminal, run command: `$ ./Scripts/collapsing_lclust.sh runNHC --help`
### Yaml Variables Used:
```
COLLAPSING_VARIABLE:
  resolution: "0_2"
NHC_VARIABLE:
  model: URPTV_pext9_igm_af
  cluster: 1
```
* resolution: used to identify correct files
* model: model to run on
* cluster: cluster to run NHC on

### Input:
  *  Previously generated cluster

### Output:
  * Results/NHC/
    * .txt files of results