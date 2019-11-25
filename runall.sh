CANCER_TYPE = $1

wget http://download.cbioportal.org/${CANCER_TYPE}_tcga.tar.gz -P /Users/azumi/Genome/
mkdir -m 777 /Users/azumi/Genome/${CANCER_TYPE}_tcga
tar -zxvf /Users/azumi/Genome/${CANCER_TYPE}_tcga.tar.gz -C /Users/azumi/Genome/${CANCER_TYPE}_tcga
Rscript count_mut.R /Users/azumi/Genome/mc3/mc3.v0.2.8.PUBLIC.maf.gz /Users/azumi/Genome/${CANCER_TYPE}_tcga/data_bcr_clinical_data_patient.txt > ${CANCER_TYPE}.mut_count.txt
Rscript PCA_plot.R ${CANCER_TYPE}.mut_count.txt /Users/azumi/Genome/${CANCER_TYPE}_tcga/data_bcr_clinical_data_patient.txt

# Rscript PCA_plot.R lusc.mut_count.txt /Users/azumi/Genome/lusc_tcga/data_bcr_clinical_data_patient.txt
# Rscript count_mut.R /Users/azumi/Genome/mc3/mc3.v0.2.8.PUBLIC.maf.gz /Users/azumi/Genome/lusc_tcga/data_bcr_clinical_data_patient.txt > lusc.mut_count.txt