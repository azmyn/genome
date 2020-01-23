# genome
codes for master thesis

# 中身

* **PCAWG_screening.py** 前処理用のpythonのコード
* **common_preprocessing.R** 実験1,2に共通する前処理部分
* **LDA_tSNE_preprocessing.R** 実験2のみの前処理部分
* **genome_analysis_3matrix_tSNE.R** 実験1の本体
* **genome_analysis_LDA_tSNE.R** 実験2の本体

# 動かす前にやること

1. https://dcc.icgc.org/releases/PCAWG/consensus_snv_indel から final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz" をダウンロードする.
2. https://dcc.icgc.org/releases/PCAWG/driver_mutations から "TableS3_panorama_driver_mutations_ICGC_samples.public.tsv.gz" をダウンロードし解凍する.
3. https://dcc.icgc.org/releases/PCAWG/clinical_and_histology から "pcawg_donor_clinical_August2016_v9.xlsx" をダウンロードする.

# 動かす順番

1. "final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz" をダウンロードした状態で最初にPCAWG_screening.pyを動かし, "f_c_p_screened.csv"を得る.
2. "f_c_p_screened.csv"を得た状態でcommon_preprocessing.Rを動かす. 
3. (実験1の場合)  genome_analysis_3matrix_tSNE.Rを動かす. 2.で作った "PCAWG_matrix_labels.csv","PCAWG_matrix_position.csv", "PCAWG_matrix_type.csv"とダウンロードした"TableS3_panorama_driver_mutations_ICGC_samples.public.tsv"が必要.
4. (実験2の場合)  LDA_tSNE_preprocessing.Rを動かす. 2.で作った "PCAWG_matrix_labels.csv","mutation_small.csv"とダウンロードした "TableS3_panorama_driver_mutations_ICGC_samples.public.tsv"が必要.
5. (実験2の場合)  genome_analysis_LDA_tSNE.Rを動かす. 4.で作った "position_lda_vector.csv", "type_lda_vector.csv", "gene_lda_vector.csv"と"pcawg_donor_clinical_August2016_v9.xlsx"が必要.
