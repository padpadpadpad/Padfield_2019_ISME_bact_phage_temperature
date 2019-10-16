# Analyses and data of:

_Padfield et al. (2019) Temperature-dependent changes to host-parasite interactions alter the thermal performance of a bacterial host. ISME_ 

DOI of paper:

DOI of preprint:

[doi: 10.1101/554717](https://doi.org/10.1101/554717)

DOI of analyses and dataset:


### Outline

This repository contains the final datasets, analyses and figures of the above-mentioned paper. It can recreate the all of the analyses and figures in the main text (Figure 1 to Figure 4) and all but one of the Supplementary Information (Figures S1-S15).

### Feedback

- Please report any problems or bugs in the code in the [Issues](https://github.com/padpadpadpad/Padfield_2019_ISME_bact_phage_temperature) tab of the GitHub repository. Alternatively, please email _d.padfield@exeter.ac.uk_.

### Licensing

This code is licensed under GPL-3.

### Running the scripts and analyses

- The project can be `cloned` or for those not familiar with GitHub, a zip file of this project can be downloaded using the "Clone or download" button at the top right of this page.
- Open the R project file in the downloaded folder. [R projects](https://support.rstudio.com/hc/en-us/articles/200526207-Using-Projects) automatically assigns the root directory to the directory in which the project resides. Consequently all of the analyses should be runnable without altering paths. These are very easy to open using RStudio. All of the scripts for the analyses can be found in `scripts/`.
- `bacteria_with_without_phage_growth_models.R` does the analysis of the logistic growth curve models of the bacteria grown in the presence and absence of phage. Recreates Figures S1:S8, Figure S10, and Figure S15.
- `bacteria_phage_tpcs.R` runs the TPC models for phage replication and bacterial growth in the presence and absence of phage. Recreates Figure 1, Figure 2 and Figure S9.
- `cost_of_resistance_logistic_growth_model.R` analyses the logistic growth models of the susceptible and resistant bacterial clones. Recreates Figure S11, Figure S12, Figure S13 and Figure S14.
- `cost_of_resistance_rolling_regression.R` analyses calculates growth rates using a rolling regression and runs the TPC models for susceptible and resistant clones. Recreates Figure 4.
- `resistance_assays.R` runs the code for the phage streak assays. Recreates Figure 3.
- `growth_curve_functions.R` contains the functions for modelling growth rates.
- All of the data needed to run the analyses are stored in `data`.
- All figures produced by the analysis are saved in `plots/` and are labelled as they are in the main text.
- 

__All analyses are done in R version 3.5.3, on macOS High Sierra 10.13.6. I am unsure whether some older version of R will support all of the packages and whether the analyses will run exactly the same.__
