import pandas as pd
import os
import json
import numpy as np

def print_file_size(file_path):
    """Prints the size of a file."""
    if os.path.exists(file_path):
        file_size = os.path.getsize(file_path)
        print(f"File size of {file_path}: {file_size} bytes")
    else:
        print(f"File {file_path} does not exist.")



def get_dataframe_memory_usage(df):
    """Returns the memory usage of a Pandas DataFrame."""
    memory = df.memory_usage(deep=True).sum() / (1024 ** 2) 
    return memory


def remove_OG(df, abundance): 
    frac = (100 - abundance) / 100
    threshold = len(df.columns) * frac
    empty_row_counts = df.isna().sum(axis=1)  
    rows_to_drop_empty_cell = empty_row_counts[empty_row_counts > threshold].index
    df_filtered_empty_cell = df.drop(index=rows_to_drop_empty_cell)
    return df_filtered_empty_cell

# orthoG = "C:/Users/46705/Documents/SpiderSilk/annotation/Orthogroups_ls.tsv"
# reads = "C:/Users/46705/Documents/SpiderSilk/big_data/full_body_rc"
# file_connectoion = "C:/Users/46705/Documents/SpiderSilk/data/raw_data/S1-S4/data_s1.csv"

orthoG = "/proj/naiss2023-6-13/miguel_analysis/orthol/orthofinder_mech_240125/Results_Jan25/WorkingDirectory/OrthoFinder/Results_Feb20/Orthogroups/Orthogroups.tsv"
reads = "/crex/proj/uppstore2019013/nobackup/private/data/bulkRNAseq/kSpiders/kallisto_quantification/individual/"
file_connectoion = "/proj/naiss2023-6-13/erika_analysis/data_s1.csv"

abundance_filename = "abundance.tsv"
info_filename = "run_info.json" 


df_orthogroups = pd.read_csv(orthoG, sep='\t')
import sys

df_connections = pd.read_csv(file_connectoion, dtype=str, sep=";")

# sub_sub = df_orthogroups.iloc[:1000, :]
sub_sub = df_orthogroups


melted_df = pd.melt(sub_sub, id_vars='Orthogroup', var_name='Species', value_name='Genes')
melted_df['Genes'] = melted_df['Genes'].str.split(', ')


exploded_df = melted_df.explode('Genes')
exploded_df['ID'] = exploded_df['Species'].str.replace('_longest_ORFs_aa$', '', regex=True)


df_DRR_species = df_connections[["ID", 'DRR', 'species', 'Genus', 'Family']]
df_DRR_species["idv_id"] = df_connections['ID'].str.extract(r'(\d+)-W', expand=False).astype(str)

merged_df = pd.merge(exploded_df, df_DRR_species,  on='ID', how='left')


# ANNOTATION

path_info  = "/crex/proj/uppstore2019013/miguel_analysis/bridgeSpider/Lsc.v1.1.genes.AA.signalP.lascID.info.csv"
# path_info  = "C:/Users/46705/Documents/SpiderSilk/annotation/Lsc.v1.1.genes.AA.signalP.lascID.info.csv"
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
    
# Group by 'Orthogroup' 
grouped_df = full_ortho_ann.groupby('Orthogroup').agg({'product': join_with_nan,
                                                       'geneName': join_with_nan,
                                                       'Group': join_with_nan}).reset_index()


df_merged_anno = pd.merge(merged_df, grouped_df, on="Orthogroup", how="left")

spidroin_df = df_merged_anno[df_merged_anno['Group']== "Spidroin"]


non_spidroin_df = df_merged_anno[(df_merged_anno['Group'] != "Spidroin") | (df_merged_anno['Group'].isna())]


df_reads = pd.DataFrame(columns=["target_id", "Cpm", "Tpm"])


# DDR_list = ['DRR296821', 'DRR296993', 'DRR297745']
DDR_list = merged_df['DRR'].unique()
unmapped_list = []
num=0

for sub_dir in DDR_list:
    file_path = f'{reads}/{sub_dir}/{abundance_filename}'
    json_path = f'{reads}/{sub_dir}/{info_filename}'
    
    if os.path.exists(file_path) or os.path.exists(json_path):
        with open(json_path, 'r') as f:
            data = json.load(f)

        n_processed = data['n_processed']
        chunk_size = 10000 
        
        for chunk in pd.read_csv(file_path, sep='\t', chunksize=chunk_size):
            chunk['Tpm'] = chunk['tpm']
            chunk['Cpm'] = (chunk['est_counts'] / n_processed) * 1000000
            df_reads = pd.concat([df_reads, chunk[['target_id', 'Cpm', "Tpm"]]], ignore_index=True)
    else:
        unmapped_list.append(sub_dir)
        

    num+=1
    print(num)

    
with open("unmapped_list.txt", "w") as file:

    file.write("\n".join(map(str, unmapped_list)))



full_df_spidroin = pd.merge(spidroin_df, df_reads, how='left', left_on='Genes', right_on='target_id')
print("merge complete")

full_df_non_spidroin = pd.merge(non_spidroin_df, df_reads, how='left', left_on='Genes', right_on='target_id')

# Adding CPM values
filled_df_non_spidroin = full_df_non_spidroin.fillna({'Group': 'NaN', 'geneName': 'NaN', 'product': 'NaN'})
filled_df_spidroin = full_df_spidroin.fillna({'Group': 'NaN', 'geneName': 'NaN', 'product': 'NaN'})

pivot_df_non_spidroin = filled_df_non_spidroin.pivot_table(
    index=['ID', 'idv_id', 'DRR', 'species', 'Genus', 'Family'], 
    columns=['Orthogroup', 'Group', "geneName", "product"], 
    values='Tpm', 
    aggfunc='sum'
).fillna(0)

pivot_df_non_spidroin.replace('NaN', np.nan, inplace=True)


pivot_df_spidroin = full_df_spidroin.pivot_table(index=['ID', 'idv_id', 'DRR', 'species', 'Genus', 'Family'], 
                               columns=['Orthogroup', 'Group', "geneName", "product"], values='Cpm', aggfunc='sum').fillna(0)

pivot_df_spidroin.replace('NaN', np.nan, inplace=True)

memory_usage = get_dataframe_memory_usage(pivot_df_spidroin)
print(f"Memory usage of pivot_df_spidroin: {memory_usage:.2f} MB")

pivot_df_non_spidroin.reset_index(inplace=True)
pivot_df_spidroin.reset_index(inplace=True)

pivot_df_spidroin.to_csv("spidroins_cpm_annotated.csv", index=True)
pivot_df_non_spidroin.to_csv("non_spidroins_tpm_annotated.csv", index=True) 
pivot_df_all = pd.merge(pivot_df_non_spidroin, pivot_df_spidroin)
pivot_df_all.to_csv("all_annotated.csv", index=True) 