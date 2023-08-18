# This folder contains the script to generate Figure 3 of the main text and Fig S3 of the supplement.

First, running the fitting scripts to generate the FOI estimates. The script, "fig3-model-fitting-cohort.R", fits an FOI model to 2 years (2019 and 2020) of age-structured incidence data from our cohort study in Kampong Speu province. The parameter estimates generated from this scripts (as well as the dataset) are located in the "data" folder. The FOI estimates for the 2019 and 2020 cohort data are called "foi-kampong-speu-2019-cohort.csv".

The script "fig3.R" pulls these FOI estimates, runs the age-prevalence model, and compares with the corresponding data to produce Figure 3 of the main text.
