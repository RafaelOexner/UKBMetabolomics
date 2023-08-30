# UKBMetabolomics
This is work of the Shah Group assessing the value of 1H-NMR serum metabolites for incident HF risk stratification within UK Biobank. A preprint, also containing detailed methodological information, can be found on medrXiv: https://doi.org/10.1101/2023.08.21.23294202

# Code
The code was run in the following order:
1. Earliest_date_ICD_DR.R
   - Calculates earliest dates for ICD diagnosis & disease status from Death records
   - Creates the file 'ICD_earliest_date_sum.tsv' and 'Death_record_status.tsv'
2. Earliest_date_OPCS4.R
   - Calculates earliest dates for OPCS4 procedures and disease status 'Yes' 'No' and 'Yes at baseline'
   - Creates file 'OPCS4_disease_status.tsv'
3. Earliest_date_VI.R
   - Determines disease status with verbal interview data (NCIC and OPC), as 'TRUE' and 'FALSE'
   - Creates file 'Verbal_interview_disease_status.tsv'
4. Endpoints.R
   - Calls in all earliest dates and disease status files, merges them into data frame TF2
   - Calculates ultimate earliest date, time to event and disease status (Yes, No, Yes at baseline)
   - Creates file 'Time_to_follow_up.tsv'
5. Reformatting_dataset.R
   - Calls in different .tsv files from 'TableDownloads' folder, reformats categorical data mostly into 'TRUE' or 'FALSE'          or NA.
   - Also creates disease status columns for DM, HTN and DM at baseline.
   - Merges all data to create 'Before_exclusion.tsv'.
   - Contains lots of additional information, e.g. diet
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
