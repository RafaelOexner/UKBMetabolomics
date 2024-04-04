# UKBMetabolomics
This is work of the Shah Group assessing the value of 1H-NMR serum metabolites for incident HF risk stratification within UK Biobank. A preprint, also containing detailed methodological information, can be found on medrXiv: https://doi.org/10.1101/2023.08.21.23294202

# Code
The code was run in the following order:
1. Earliest_date_ICD_DR.R
   - Calculates earliest dates for ICD diagnosis & disease status from Death records
   - Creates the file 'ICD_earliest_date_sum.tsv' and 'Death_record_status.tsv'
2. Earliest_date_OPCS4.R
   - Calculates earliest dates for OPCS4 procedures and disease status 
   - Creates file 'OPCS4_disease_status.tsv'
3. Earliest_date_VI.R
   - Determines disease status with verbal interview data (NCIC and OPC)
   - Creates file 'Verbal_interview_disease_status.tsv'
4. Endpoints.R
   - Calls in all earliest dates and disease status files and merges them 
   - Calculates ultimate earliest date, time to event and disease status 
   - Creates file 'Time_to_follow_up.tsv'
5. Reformatting_dataset.R
   - Calls in additional stratification information
   - Reformats categorical data
   - Creates disease status columns for DM, HTN and DM at baseline.
   - Merges with endpoint data created above
7. Metabolomics_Models_Analysis.R (main analysis)
   - Performs reformatting and exclusion
   - Calculates baseline characteristics
   - Assesses individual metabolite associations
   - Partition and model training
   - Model assessment:
        - Absolute & relative
        - Discrimination & calibration
        - Survival stratification
        - Insights into feature utilisation
