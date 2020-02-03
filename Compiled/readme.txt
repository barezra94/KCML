KCML Executable

1. Prerequisites for Deployment 

Verify that version 9.7 (R2019b) of the MATLAB Runtime is installed.   
If not, you can run the MATLAB Runtime installer.
To find its location, enter
  
    >>mcrinstaller
      
at the MATLAB prompt.
NOTE: You will need administrator rights to run the MATLAB Runtime installer. 

Alternatively, download and install the Windows version of the MATLAB Runtime for R2019b 
from the following link on the MathWorks website:

    https://www.mathworks.com/products/compiler/mcr/index.html
   
For more information about the MATLAB Runtime and the MATLAB Runtime installer, see 
"Distribute Applications" in the MATLAB Compiler documentation  
in the MathWorks Documentation Center.

2. Files to Deploy and Package

Files to Package for Standalone 
================================
-KCML.exe
-MCRInstaller.exe 
    Note: if end users are unable to download the MATLAB Runtime using the
    instructions in the previous section, include it when building your 
    component by clicking the "Runtime included in package" link in the
    Deployment Tool.
-This readme file 


3. Other requirements to run KCML please provide 
- desired values in config.txt including experiment name, sigma values, and number of cross validation during training.

- annotations files in annotations folder. The first column should have gene names and the second column should have functional annotation name as in the included example file. The name of this file should follow the following format: $expName$Annotation. For example if expName=Demo then annotation file name should be DemoAnnotation.csv.

- data files in annotations folder. The first column should have gene names as in the included example file.The name of this file should follow the following format: $expName$Data. For example if expName=Demo then data file name should be DemoData.csv.

KCML Output:
- The prediction of KCML on various functional terms and its confidence in these predictions 
- The performance of the classifiers
- The selected features by different classifiers

If you use this program, please cite
Salem H., Rittscher J., Pelkmans L., 2019, KCML: a machine-learning framework for inference of multi-scale gene functions from genetic perturbation screens, bioRxiv 761106


4. Definitions

For information on deployment terminology, go to
https://www.mathworks.com/help and select MATLAB Compiler >
Getting Started > About Application Deployment >
Deployment Product Terms in the MathWorks Documentation
Center.




