# This folder contains the R script to generate Figure2 of the main text and Figure S2 of the supplement.

First, the script "fig2-model-fitting.R" fits an FOI model to the 18-year age-structured national dengue data for Cambodia ("DENV-Nat-Aged.csv") which is located in the "data" subfolder of the main repo. It results in an output file of FOI estimates and corresponding confidence intervals which is also stored in the "data" folder. This file is called "foi-fit-national.csv".

Then, the script, "fig2.R", loads that csv with its foi estimates, runs the model, and compares with data to generate figure 2 and S2. You can find the final figures "fig2.png" and "figS2.png" in the folder "final-figures".
