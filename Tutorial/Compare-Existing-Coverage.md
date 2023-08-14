## Navigation
### Go back to <a href = "/Tutorial/Home.md">Wiki</a>
<hr>

# Compare Existing Coverage
* This compares ped.txt, coverage_existing.sample.txt and the case list to identify what samples have been removed
* Must be run after initalCoverage fails and before running existingCoverage

**Must be run on dev1-dev4**

Command: `./Scripts/collapsing_lclust.sh compareCoverage`
```
Usage:
  identify_removed_samples.py [-h] --caselist CASELIST --ped_txt PED_TXT --coverage_existing COVERAGE_EXISTING --out OUT

Program to review samples removed to create ped.txt and coverage_exisitng

Arguments:
  -h, --help            show this help message and exit
  --caselist CASELIST   file path for case list
  --ped_txt PED_TXT     file path for ped.txt
  --coverage_existing COVERAGE_EXISTING
                        file path for coverage_existing
  --out OUT             output directory
```

### Input:
* ped.txt
* coverage_existing.sample.txt
* case list txt file

### Output:
* /Sample_check
  * log file
  * txt file case samples not in ped.txt
  * txt file samples removed from ped.txt to create coverage_existing.sample.txt
  * txt file case samples among ones removed from ped.txt

Next: [Existing Coverage of CaseControl Cohort](https://github.com/igm-team/ClusteredCollapsing/wiki/Existing-Coverage-of-CaseControl-Cohort)
