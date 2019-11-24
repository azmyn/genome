CANCER_TYPE=$1

#python count_mut.py mc3.v0.2.8.PUBLIC.maf.gz ${CANCER_TYPE} > ${CANCER_TYPE}.mut_count.txt
wget http://download.cbioportal.org/${CANCER_TYPE}_tcga.tar.gz
# python process_clinical_info.py ${CANCER_TYPE}_tcga.tar.gz > ${CANCER_TYPE}.clinical_info.txt
Rscript count_mut.R /Users/azumi/Genome/mc3/mc3.v0.2.8.PUBLIC.maf.gz > ${CANCER_TYPE}.mut_count.txt

# Rscript PCA_plot.R ${CANCER_TYPE}.mut_count.txt ${CANCER_TYPE}.clinical_info.txt 