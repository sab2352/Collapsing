## Navigation
### Go back to <a href = "/Tutorial/Home.md">Wiki</a>
<hr>

# CMH Permutations 
* Runs permutations

**Run on QS1**

Command: `./Scripts/collapsing_lclust.sh cmhPerm`

```
 Usage: 
  cmh_permute_lclust_short_script.R --project project --res res --minsample minsample --model model --permstart permstart --permend permend --cores cores
  
  Options:
    -p, --project     path to Results directory of project
    -d, --dir         name of group
    -r, --res         default="0_3", resolution
    -n, --minsample   default=5, minimum number of cases/controls per cluster
    -m, --model       default="dominantSynonymous", model to permute
    -s, --permstart   default=1, first permutation
    -e, --permend     default=1, last permutation
    -c, --cores       default=5, number of cores
```
If you want to view help message in terminal, run command: `$ ./Scripts/collapsing_lclust.sh cmh --help`

### Yaml Variables Used:
```
COLLAPSING_VARIABLE:
  resolution: "0_2"
  min_sample: 20
PERMUTATION_VARIABLE:
  permstart: 1
  permend: 100
  cores: 15
  perm_models: "dominantSynonymous"
```

* resolution: used to identify correct files
* min_sample: minimum number of cases/controls of included clusters
* permstart: where the permutation starts, first permutation
* permend: where the permutation ends, last permutation
* cores: number of cores to be used
* perm_models: model(s) to be used for permutations 

### Input:
* Uses previously generated models

### Output:
* /Results/Permutations
  * RDS files of the permutations 

