# Power Collapsing 

## performs power or odds ratio calculations 

Usage: 
  power_collapsing.R [--calc=<calc>] [--alpha=<alpha>] [--odds_ratio=<odds_ratio>] [--case_group=<case_group>]
  
  Options:
  -h --help
  --calc=<calc> what we are calculating, can be "power" or "es" (odds ratio) [default: power]
  --alpha=<alpha> type 1 error rate, or "Alpha" [default: 0.08]
  --odds_ratio=<odds_ratio> the detectable effect size (or odds ratio in the case of a binary outcome variable)[default: 3]
  --case_group=<case_group> label for analysis. example might "picu" or "IPF" [default: temp_case]