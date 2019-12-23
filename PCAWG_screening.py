import csv
import pandas as pd

# ����������
path_w = '/Users/azumi/Dropbox/KU/shimolab_2019/genome/test.csv'


# �i���z�ߤ�SNP�γ��
df = pd.read_csv('/Users/azumi/Genome/final_consensus_passonly.snv_mnv_indel.icgc.public.maf', sep='\t', index_col=0,usecols=[1,2,3,5,6,7,8,9,12,15,41,42])
df = df.query('Variant_Type == "SNP"')

# �v���ζ��x
alelle_comp = {"A":"T", "T":"A", "G":"C", "C":"G"}


# �䮐�η����yһ
list(df.reset_index().query('Reference_Allele == "C"').index)
# with open(path_w, mode='w') as f:
df.to_csv(path_w)



# with codecs.open(gzip_path, 'r','utf-8', 'ignore') as f:
#     s = f.readlines()
#     print(s[1])
