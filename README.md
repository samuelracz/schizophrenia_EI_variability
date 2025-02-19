# Supplementary code for "Reduced temporal variability of cortical excitation/inhibition ratio in schizophrenia" by Racz et al. (2025).

This folder contains the complete analysis pipeline for Racz et al. (2025) [Reduced temporal variability of cortical excitation/inhibition ratio in schizophrenia](https://www.nature.com/articles/s41537-025-00568-3) (DOI: https://doi.org/10.1038/s41537-025-00568-3). The scripts reproduce the analysis outcomes and figures published in the article. Specifically,

- 'script_01_analyze_IRASA_bimodal.m' loads *pre-processed* EEG data and conducts the bimodal IRASA analysis, saving outputs (per participant) into results_IRASA_ec/ 
   as Matlab workspaces.
- 'script_02_statistics_IRASA_bimodal.m' conducts all statistical analyses presented in the article.
- 'script_03_plot_results.m' re-creates all the figures in the article and its supplementary material.

Additionally, 'script_00_preproc_ec_mara.m' provides the means to reproduce the results from raw EEG data. Raw EEG recordings are available at Zenodo.org in the repository [Resting-state EEG, clinical, and demographics data from schizophrenia patients and age-matched healthy controls](https://zenodo.org/records/14808296) (DOI: 10.5281/zenodo.14808296). The list of included participants and the specific position of the 30s of EEG data selected for analysis is contained in the matlab workspace miscellaneous/fnames_times_ec.mat.

Note that the code loads the pre-processed EEG data from Matlab workspaces (in data_ec_mara_avg_256Hz_mat/); however, pre-processed EEG is also provided in .edf format (in data_ec_mara_avg_256Hz/).

**When using this resource, please cite** Racz, F.S., Farkas, K., Becske, M. *et al*. Reduced temporal variability of cortical excitation/inhibition ratio in schizophrenia. *Schizophr* **11**, 20 (2025). https://doi.org/10.1038/s41537-025-00568-3

Frigyes Samuel Racz

The University of Texas at Austin

email: fsr324@austin.utexas.edu

2025.
