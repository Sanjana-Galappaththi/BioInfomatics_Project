# Master script: runs the full analysis in order.
# Set the working directory to the project root first.
source("C:/Users/User/Desktop/4.1sem/BioInfomatics/ICA03_project/ICA03_Bioinformatics_Dexamethasone/code/01_preprocessing.R")
source("C:/Users/User/Desktop/4.1sem/BioInfomatics/ICA03_project/ICA03_Bioinformatics_Dexamethasone/code/02_exploratory_analysis.R")
source("C:/Users/User/Desktop/4.1sem/BioInfomatics/ICA03_project/ICA03_Bioinformatics_Dexamethasone/code/03_differential_expression.R")
source("C:/Users/User/Desktop/4.1sem/BioInfomatics/ICA03_project/ICA03_Bioinformatics_Dexamethasone/code/04_visualisations.R")
writeLines(capture.output(sessionInfo()), "results/session_info.txt")