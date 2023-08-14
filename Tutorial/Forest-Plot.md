## Navigation
### Go back to <a href = "/Tutorial/Home.md">Wiki</a>
<hr>

# Forest Plot 

* Creates a forest plot
* For when we run forest plots after geneset analyses


**Must be run on dev1-dev4**

### Input 
The input file for this function is a .csv file
The file has the following columns: 
1. path_to_git (absolute path to repo)
2. model 
3. resolution
4. min_sample
5. analysis 
6. row_name

Command: `$ ./Scripts/collapsing_lclust.sh forestPlot_display`
```
# fp_forest_plot.R: creates a forest plot, Created by Josh.

'Usage: 
  fp_forest_plot.R --output_directory=<output_directory> --path_to_fp_data=<path_to_fp_data> [--title_str=<title_str>] [--left_border=<left_border>] [--left_label_width=<left_label_width>] [--num_groups_var=<num_groups_var>] [--color_var=<color_var>] [--pub_textsize=<pub_textsize>] [--pub_breaks_log=<pub_breaks_log>] [--pub_leg_pos=<pub_leg_pos>] [--pub_breaks=<pub_breaks>] [--label_height_factor=<label_height_factor>] [--group_label_left_border=<group_label_left_border>] [--pub_width=<pub_width>] [--pub_height=<pub_height>] [--log_flag=<log_flag>]  [--group_label_alpha=<group_label_alpha>] [--x_val_factor=<x_val_factor>] [--pub_text_start_buf=<pub_text_start_buf>] [--axis_text_size=<axis_text_size>] [--row_text_size=<row_text_size>] [--group_text_size=<group_text_size>] [--p_val_textsize=<p_val_textsize>]  [--p_line_space=<p_line_space>]     [--y_axis_alpha=<y_axis_alpha>]     [--palette_index=<palette_index>]     [--group_color_palette_index=<group_color_palette_index>] [--dashed_alpha=<dashed_alpha>]     [--y_axis_width=<y_axis_width>]     [--x_axis_height=<x_axis_height>]      [--show_x_axis_title=<show_x_axis_title>]     [--color_margin=<color_margin>]     [--y_top_buffer=<y_top_buffer>]     [--y_bottom_buffer=<y_bottom_buffer>]    [--p_alpha_p_val=<p_alpha_p_val>]     [--p_val_hjust=<p_val_hjust>]    [--pointrange_size=<pointrange_size>]     [--arrow_size=<arrow_size>]     [--x_axis_line_end_buffer=<x_axis_line_end_buffer>]     [--sig_line_alpha=<sig_line_alpha>] [--show_per_w_p=<show_per_w_p>]     [--margin_vector=<margin_vector>]  [--y_max_pub=<y_max_pub>]  [--group_label_vert_axis_alpha=<group_label_vert_axis_alpha>] [--pub_width_table_mod=<pub_width_table_mod>] [--show_ci_w_p=<show_ci_w_p>]

  Options:
  -h --help
    --output_directory=<output_directory> path to output directory
    --path_to_fp_data=<path_to_fp_data> path to data for forest plot. Requires CSV with headers 
    --title_str=<title_str> string for labelling output files [default: temp]
    --left_border=<left_border> the left border of the plot. More negative gives more room for low ORs. can adjust left border of publication output if need image wider. More negative makes image wider. [default: -2]
    --left_label_width=<left_label_width> distance from left border to left part of page. Larger number gives more room for label. [default: 2]
    --num_groups_var=<num_groups_var> Variable used to determine how many groups for color shading. [default: Group_1]
    --color_var=<color_var> Variable used to determine variable for coloring points. [default: Group_1]
    --pub_textsize=<pub_textsize> Text size in publication figure. [default: 5]
    --pub_breaks_log=<pub_breaks_log> Breaks for axis in publication figure. [default: -2,0,2]
    --pub_leg_pos=<pub_leg_pos> Controls position of legend. [default: none]
    --pub_breaks=<pub_breaks> TODO [default: 1,1.5,2,4,6,8]
    --label_height_factor=<label_height_factor> [default: 0.95]
    --group_label_left_border=<group_label_left_border> [default: 4]
    --pub_width=<pub_width> [default: 8]
    --pub_height=<pub_height> [default: 8]
    --log_flag=<log_flag> [default: TRUE]
    --group_label_alpha=<group_label_alpha> set to 0 if no group label desired [default: 1]
    --x_val_factor=<x_val_factor> compresses or expands along the original x (then Y) axis [default: 1]
    --pub_text_start_buf=<pub_text_start_buf> [default: 0.1]
    --axis_text_size=<axis_text_size> [default: 10]
    --row_text_size=<row_text_size> [default: 4]
    --group_text_size=<group_text_size> [default: 10]
    --p_val_textsize=<p_val_textsize> [default: 3]
    --p_line_space=<p_line_space> [default: 0.3]
    --y_axis_alpha=<y_axis_alpha> [default: 1]
    --palette_index=<palette_index>  index of color to be used for each line. must be a list of integers >=1 with , separation and no spaces [default: 1,2,3,4,5,6,7,8]
    --group_color_palette_index=<group_color_palette_index>] must be a list of integers >=1 with , separation and no spaces [default: 1,2,3,4,5,6,7,8]
    --dashed_alpha=<dashed_alpha> translucency of dashed line which divides groups along vertical axis. Set to 0 to hide. [default: 0.5]
    --y_axis_width=<y_axis_width> [default: 0.8]
    --x_axis_height=<x_axis_height> [default: 1.5]
    --show_x_axis_title=<show_x_axis_title> [default: TRUE]
    --color_margin=<color_margin> [default: 0.5]
    --y_top_buffer=<y_top_buffer> [default: 0.5]
    --y_bottom_buffer=<y_bottom_buffer> [default: 0]
    --p_alpha_p_val=<p_alpha_p_val> value which will filter out point labels with OR and P. Set to -1 to see none, set to 2 to see all [default: 0.05]
    --p_val_hjust=<p_val_hjust> Determines position of p value label relative to line [default: 0.15]
    --pointrange_size=<pointrange_size> [default: 0.5]
    --arrow_size=<arrow_size> [default: 0.5]
    --x_axis_line_end_buffer=<x_axis_line_end_buffer> [default: 0]
    --sig_line_alpha=<sig_line_alpha> [default: 0.2]
    --show_per_w_p=<show_per_w_p> [default: TRUE]
    --margin_vector=<margin_vector> [default: 0.5,0,0,0]
    --y_max_pub=<y_max_pub> [default: 3.1]
    --group_label_vert_axis_alpha=<group_label_vert_axis_alpha> set to 0 to hide labels on vertical axis [default: 1]
    --pub_width_table_mod=<pub_width_table_mod> adjust table width relative to graph wide. e.g. 1 means will be one inch wider than forest plot [default: 0]
    --show_ci_w_p show=<show> confidence intervals over points in forest plot [default: TRUE]
    
```
If you want to view help message in terminal, run command: `$ ./Scripts/collapsing_lclust.sh forestPlot_display --help`
### Yaml Variables Used:
```
ADDITIONAL_PARAMETERS:
  forestPlot: ""
```
* forestPlot: enter any additional parameters you wish to specify between the quotes
  * should be look like forestPlot: "--y_max_pub 3.1 --y_top_buffer 0.5"

### Output:
  * forestplot