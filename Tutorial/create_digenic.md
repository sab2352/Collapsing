## Navigation
### Go back to <a href = "/Tutorial/Home.md">Wiki</a>
<hr>

# Create Digenic 
### Creates a digenic model out of a single model 
* Will take in an individual model and then create a model where the columns are Sample.Name and Digene which is gene1_gene2 if the Sample.Name has a QV in gene1 and a QV in gene2. 
* Will allow for qq plots in specific gene1_gene2 pairs but probably will be more helpful in gene set analyses.

Command:
`./Scripts/collapsing_lclust.sh create_gene_set_model`

### Usage: 
```
  digenic_create_model.R --resolution_var=<resolution_var> --model=<model>
  -h --help
  --resolution_var=<resolution_var> resolution for clustering. Typically 0.1-0.4
  --model=<model> 
```
  
