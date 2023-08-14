# Fast Gene Set Enrichment Analysis (FGSEA)

This is based off of Fast Gene Set Enrichment Analysis: https://github.com/ctlab/fgsea

```
Usage: 
  fgsea_run.R --resolution_var=<resolution_var> --model=<model> --min_sample=<min_sample> [--case_group=<case_group>] [--debug=<debug>] 
  
  Options:
  -h --help
  --resolution_var=<resolution_var> resolution for clustering.
  --model=<model> will run on a model if specified. If specified, no need for gene_list_case_path or gene_list_ctrl_path
  --min_sample=<min_sample> minimum size of case/control allowed to include cluster into combined cluster.
  --case_group=<case_group> [default: temp_group]
  --debug=<debug> [default: FALSE]
  ```