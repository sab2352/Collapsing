## Navigation
### Go back to <a href = "/Tutorial/Home.md">Wiki</a>
<hr>

# Create QQ Plots

* Produces QQ plots: This page is under construction

**Runs on QS1**

Command: `./Scripts/collapsing_lclust.sh qq`

```
 Usage: 
  cmh_qq.R --resolution_var=<resolution_var> --min_sample=<min_sample> --nperm_end=<nperm_end> [--top_plot=<top_plot>] [--sig_line=<sig_line>] [--text_label_color=<text_label_color>] [--y_lim_low=<y_lim_low>] [--title_size=<title_size>] [--case_group=<case_group>]
  
  Options:
  -h --help
  --min_sample=<min_sample> minimum number of cases/controls of included clusters
  --resolution_var=<resolution_var> resolution for clustering.
  --nperm_end=<nperm_end> number of permutations to incorporate 
  --top_plot=<top_plot> max y value displayed in qq plot [default: 8.9]
  --sig_line=<sig_line> p value to add horizontal line for significance [default: NA]
  --text_label_color=<text_label_color> color of text for label [default: black]
  --y_lim_low=<y_lim_low> i think this may never get used. will remove in future versions. [default: 6]
  --title_size=<title_size> Font size for title [default: 9]
  --case_group=<case_group> label for analysis. example might "picu" or "IPF" [default: temp_case]
 
```
If you want to view help message in terminal, run command: `$ ./Scripts/collapsing_lclust.sh qq --help`
### Yaml Variables Used:
```
USER_VARIABLE:
  case_group: project_name
COLLAPSING_VARIABLE:
  resolution: "0_2"
  min_sample: 20
  permend: 100
ADDITIONAL_PARAMETERS:
  cmh_qq: ""
```
* case_group: project name, added to output file(s) name
* resolution: used to identify correct files
* min_sample: minimum number of cases/controls of included clusters
* permend: number of permutations to incorporate

### Input:
* .xlsx cmh excel file
* .txt cluster sizes

### Output:
* /Results/CMH
  * .png qq plots

When running permutations, typically 100 permutations is sufficient data to generate preliminary QQ plots.  For publication-quality plots with more precise confidence intervals, 1000 permutations can be run.
