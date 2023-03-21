# BreastCancerSNP
Supplementary Data and Codes for "A nomogram based on genotypic and clinicopathologic factors to predict the non-sentinel lymph node metastasis in Chinese women breast cancer patients".

Abstract
Background: Sentinel lymph node biopsy (SLNB) is the standard therapy for patients with breast cancer with clinically negative axilla. However, axillary lymph node dissection (ALND) remains the standard of care for patients with sentinel lymph node (SLN) positivity. Clinical data reveals that approximately 40–75% of patients do not develop non-sentinel lymph node (NSLN) metastasis after ALND.
Methods: We retrospectively collected data from 1,879 women with breast cancer enrolled in multiple centers. Genotypic features include 96 single nucleotide polymorphisms (SNPs) associated with breast cancer susceptibility, therapy, and prognosis.
Results: Among 229 patients who underwent SLNB followed by ALND without neo-adjuvant therapy, 79 (34%) had positive axillary NSLN metastasis. LDA-RFE identified characteristics including lymphovascular invasion, number of positive and negative SLNs, and two SNP clusters as significant NSLN metastasis predictors. Furthermore, SVM-RFE identified 29 significant SNPs for predicting NSLN metastasis. In the internal validation, the median AUCs of the clinical and all SNPs combined, the clinical and 29 significant SNPs combined, and the clinical models were 0.837, 0.795, and 0.708, respectively. Conversely, in the external validation, the AUCs of the three models were 0.817, 0.815, and 0.745, respectively.
Conclusion: We present a new nomogram that combines genotypic and clinicopathological factors to achieve higher sensitivity and specificity than the traditional clinicopathological factors in predicting NSLN metastasis in Chinese women with breast cancer. However, more validation is required in prospective studies among different patient populations.


The preprocessed data and model building R codes are listed here.

EDC_code.R: The R codes processing the input files

input.zip: Compressed input files
