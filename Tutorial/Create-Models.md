## Navigation
### Go back to <a href = "/Tutorial/Home.md">Wiki</a>
<hr>

# Create Model

* We next filter the master variant file to create specific models
* Enter the model num(s) you would like to run in the Yaml file

**Must be run on qs1**

Command: `./Scripts/collapsing_lclust.sh models`
### YAML Variables Used:
```
MODEL_VARIABLE:
  model_num: "1 4 5 7-10 13"
  model_1:...
```
* model_num: model numbers you wish to run
  * Enter model number(s)/range of models in YAML file you wish to run
  * For example, running model_1-model_4 would be model_num: "1-4"

### Input:
* Uses the genotypes of each cluster found at /Results/Coverage

### Output:
* In /Results/Collapsing, within each cluster's directory there will be a folder for each model

### Special Note: 
* When you name your models, make sure that each name is less than 30 characters in length. This is because model names will be used for the excel sheet tabs and it has a 30 character limit. If this rule is not met, the cmh step will not run properly! 