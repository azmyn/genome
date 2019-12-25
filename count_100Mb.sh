CANCER_TYPE=$1

Rscript count_mut_100Mb.R /Users/azumi/Genome/mc3/mc3.v0.2.8.PUBLIC.maf.gz /Users/azumi/Genome/${CANCER_TYPE}_tcga/data_bcr_clinical_data_patient.txt ${CANCER_TYPE}

# Rscript PCA_plot.R lusc.mut_count.txt /Users/azumi/Genome/lusc_tcga/data_bcr_clinical_data_patient.txt
# Rscript count_mut.R /Users/azumi/Genome/mc3/mc3.v0.2.8.PUBLIC.maf.gz /Users/azumi/Genome/lusc_tcga/data_bcr_clinical_data_patient.txt > lusc.mut_count.txt