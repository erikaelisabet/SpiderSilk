import pandas as pd
import os
import json
import numpy as np

orthoG = "C:/Users/46705/Documents/SpiderSilk/data/raw_data/Orthogroups.GeneCount.tsv"
file_connectoion = "C:/Users/46705/Documents/SpiderSilk/data/raw_data/S1-S4/data_s1.csv"


df_counts = pd.read_csv(orthoG, sep='\t')


df_connections = pd.read_csv(file_connectoion, dtype=str, sep=";")

sub_sub = df_counts.iloc[:10, :]
# # sub_sub = df_orthogroups

# exploded_df['ID'] = exploded_df['Species'].str.replace('_longest_ORFs_aa$', '', regex=True)

df_DRR_species = df_connections[["ID", 'DRR', 'species', 'Genus', 'Family']]
df_DRR_species["idv_id"] = df_connections['ID'].str.extract(r'(\d+)-W', expand=False).astype(str)

# merged_df = pd.merge(exploded_df, df_DRR_species,  on='ID', how='left')

# # ANNOTATION

# # path_info  = "/crex/proj/uppstore2019013/miguel_analysis/bridgeSpider/Lsc.v1.1.genes.AA.signalP.lascID.info.csv"
path_info  = "C:/Users/46705/Documents/SpiderSilk/data/annotations/Lsc.v1.1.genes.AA.signalP.lascID.info.csv"

df_info = pd.read_csv(path_info, sep=',')

df_info["Group"].unique()
sp = len(df_info[df_info["Group"] == "Spidroin"])
print(f"nr of Spidroin genes: {sp} ")
spiCE = len(df_info[df_info["Group"] == "SpiCE"])
print(f"nr of SpiCE genes: {spiCE} ")
part_df = sub_sub[["Orthogroup", "Lsc.v1.1.genes.AA.signalP.lascID"]]

part_df=part_df.dropna()
melted_ann = pd.melt(part_df, id_vars='Orthogroup', var_name='Lsc.v1.1.genes.AA.signalP.lascID', value_name='Genes')
melted_ann['Genes'] = melted_ann['Genes'].str.split(', ')

exploded_df = melted_ann.explode('Genes')

full_ortho_ann = pd.merge(exploded_df, df_info, how='left', left_on='Genes', right_on='lascID')

# grouped_df = full_ortho_ann.groupby('Orthogroup').agg({'product': list, 'geneName': list,'Group': lambda x: list(x.unique())}).reset_index()
def join_with_nan(x):
    non_null_values = x.dropna() 
    if len(non_null_values) > 0:  
        return ','.join(non_null_values.astype(str).unique())  
    else:
        return np.nan  
    
grouped_df = full_ortho_ann.groupby('Orthogroup').agg({'product': join_with_nan,
                                                       'geneName': join_with_nan,
                                                       'Group': join_with_nan, 'lascID': join_with_nan}).reset_index()
df_merged_anno = pd.merge(merged_df, grouped_df, on="Orthogroup", how="left")
