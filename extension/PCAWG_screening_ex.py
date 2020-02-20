import csv
import pandas as pd

# write
path_w = '/Users/azumi/Dropbox/KU/shimolab_2019/genome/f_c_p_screened_withDNP.csv'


# read and exract SNP
df = pd.read_csv('/Users/azumi/Genome/final_consensus_passonly.snv_mnv_indel.icgc.public.maf', sep='\t', index_col=0,usecols=[1,2,3,5,6,7,8,9,12,15,41,42])
# df = df.query(Variant_Type in ('SNP','DNP'))
df = df[(df["Variant_Type"] == "SNP") | (df["Variant_Type"] == "DNP")| (df["Variant_Type"] == "INS")| (df["Variant_Type"] == "DEL")]


# u = df['Variant_Type'].unique()
# print(u)
d = df["Variant_Type"].value_counts().to_dict()
print(d)

df.to_csv(path_w)

