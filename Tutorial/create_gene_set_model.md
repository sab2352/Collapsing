## Navigation
### Go back to <a href = "/Tutorial/Home.md">Wiki</a>
<hr>

# Create Gene Set Model 

Command:
`./Scripts/collapsing_lclust.sh create_gene_set_model`

### Usage: 
```
fp_create_gene_set.R --model=<model> --resolution=<resolution> --cluster=<cluster> --gene_set_path=<gene_set_path> --gene_set_name=<gene_set_name> --ccds_genes_path=<ccds_genes_pathj> [--debug=<debug>]
  -h --help
    --model=<model> model_to
    --resolution=<resolution> resolution used to create clusters
    --cluster=<cluster> 
    --gene_set_path=<gene_set_path> path to gene set file. should be a .txt file. no header. First line of gene_set text file should start with # and be a description of the gene set. This line will be skipped when reading in file.
    --gene_set_name=<gene_set_name> name of gene set. no spaces.
    --ccds_genes_path=<ccds_genes_path> path to ccds genes in the atav database. will reject gene list and error if some genes are not in database
    --debug=<debug> [default: FALSE]
```