## Navigation
### Go back to <a href = "/Tutorial/Home.md">Wiki</a>
<hr>


# CMH

* In this step, an odds ratio and p-value is calculated for each gene to assess association

**Run on QS1**

Command: `$ ./Scripts/collapsing_lclust.sh cmh`

```
Usage: 
  cmh_lclust.R --resolution_var=<resolution_var> --min_sample=<min_sample> [--models_to_exclude=<models_to_exclude>] [--write_excel_only=<write_excel_only>]
  
  Options:
  -h --help
  --min_sample=<min_sample> minimum number of cases/controls of included clusters
  --resolution_var=<resolution_var> resolution for clustering.
  --models_to_exclude=<models_to_exclude> exclude these models. Typically the master genotype model. [default: dominantNoneMAF]
  --write_excel_only=<write_excel_only> If all models are done but excel not correctly written, change to TRUE [default: FALSE] 
```
If you want to view help message in terminal, run command: `$ ./Scripts/collapsing_lclust.sh cmh --help`
### Yaml Variables Used:
```
COLLAPSING_VARIABLE:
  resolution: "0_2"
  min_sample: 20
  models_to_exclude: "dominantNoneMAF"
```
* resolution: used to identify correct files
* min_sample: minimum number of cases/controls of included clusters
* models_to_exclude: models that will be excluded. Typically just the master genotype model
* write_excel_only: state if you need the excel output only. If all models are done but excel not correctly written, change to TRUE

### Input:
  *  Previously generated models

### Output:
  * /Results/CMH
    * .xlsx file that contains an odds ratio and p-value for each gene
    * .RDS file