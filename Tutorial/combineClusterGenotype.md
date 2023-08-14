# Combine Cluster Genotype 

## This combines the genotype files from a model across clusters.
* Has the option to save a file for total, case-only, and ctrl-only genotypes. 

* This looks at the models from the YAML file under the variable: perm_models

* The default sets combine_ctrl_output to FALSE. This means that only the cases file will be produced. If the combine_ctrl_output is set to TRUE, then there will be three files that are outputted: a case-only file, a ctrl-only file and a total fie.

```
Usage: 
  review_genotypes_from_model.R --resolution=<resolution> [--min_sample=<min_sample>] [--case_group=<case_group>] [--combine_ctrl_output=<combine_ctrl_output>] [--debug=<debug>]


  Options:
  -h --help
    --resolution=<resolution> resolution used to create clusters
    --min_sample=<min_sample> [default: NA]
    --models_to_exclude=<models_to_exclude> [default: NA]
    --case_group=<case_group> [default: temp]
    --combine_ctrl_output=<combine_ctrl_output> [default: False]
    --debug=<debug> [default: FALSE]
```