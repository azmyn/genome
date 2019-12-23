import csv
import pandas as pd

# 書き出し先
path_w = '/Users/azumi/Dropbox/KU/shimolab_2019/genome/test.csv'


# 読み込みとSNPの抽出
df = pd.read_csv('/Users/azumi/Genome/final_consensus_passonly.snv_mnv_indel.icgc.public.maf', sep='\t', index_col=0,usecols=[1,2,3,5,6,7,8,9,12,15,41,42])
df = df.query('Variant_Type == "SNP"')

# 関数の定義
alelle_comp = {"A":"T", "T":"A", "G":"C", "C":"G"}


# 変異の方向を統一
list(df.reset_index().query('Reference_Allele == "C"').index)
# with open(path_w, mode='w') as f:
df.to_csv(path_w)



# with codecs.open(gzip_path, 'r','utf-8', 'ignore') as f:
#     s = f.readlines()
#     print(s[1])
