## Navigation
### Go back to <a href = "/Tutorial/Home.md">Wiki</a>
<hr>

# Create Cluster Map

Run the following code snippet to create cluster map.

```
echo "Creating umap with clusters and tables plots ..."

$Rscript_var $umap_clusters_fig --help

$Rscript_var $umap_clusters_fig --resolution_var $resolution --min_sample 20 --case_group $case_group

exit

```