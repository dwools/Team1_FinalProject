# TIGRESS 


**Trustful Inference of Gene REgulation using Stability Selection (TIGRESS)** is a Gene Regulatory Network inferencing algorithm that uses least angle regression. 


## Running the notebook: 

Simply make sure that an expression matri file, its metadata file, and a TF list file, are located in the resources section of the project, and make sure to change the file names if using the non-default data. There should also be an sc.data and sc.metadata file for VIPER to compare with.

In the section where TIGRESS is run, if you dont want to run TIGRESS on the full dataset, ensure that use.truncated is set to TRUE. 




## Output files: 

One file for each cell_type is output while running TIGRESS, there is a 'full' output file that gets plugged into VIPER, and then there is a VIPER plot file. 


## Documentation: 

Any other TIGRESS documentation can be found at the [TIGRESS GitHub](https://github.com/jpvert/tigress)