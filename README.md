# UKBMetabolomics
This is work of the Shah Group assessing the value of 1H-NMR serum metabolites for incident HF risk stratification within UK Biobank. A preprint, also containing detailed methodological information, can be found on medrXiv: https://doi.org/10.1101/2023.08.21.23294202

# Code
The code was run in the following order:
1. adf
   - Reformats ICD tables
   - Searches for earliest diagnoses
2. Earliest_date_ICD_DR.R
3. Earliest_date_OPCS4.R
4. Earliest_date_VI.R
5. Endpoints.R
6. Metabolomics_Models_Analysis.R (main analysis)
   - Reformatting, exclusion
   - Baseline characteristics
   - Individual metabolite associations
   - Partition and model building
   - Model assessment
        - Absolute & relative
        - Discrimination & calibration
        - Survival stratification
        - Insights into feature utilisation
