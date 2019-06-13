INSTRUCTIONS ON HOW TO RUN GWIZ.

Go to folder “Data”
Here there are two example csv files with GWAS data.
Using the format of those files as a template, load your own data into a csv file in the folder “data”. Take special care to ensure the column names are written exactly as they are in the example csv files. 

Now load “GWAS_call_batch.R”
Go to line 3 and change the working directory to the appropriate location.
Select all and run the code.
A table called aucout should appear with your results.

After running the code two folders will also be created in your working directory. “Results data” contains the coefficients of the logistic regression model created, ROC curve plots and the AUROC.
“Resampled data” contains the simulated population created from your GWAS summary level data.

Congrats you are now an expert in using GWIZ!