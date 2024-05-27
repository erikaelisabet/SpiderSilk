import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

all_path = "C:/Users/46705/Documents/SpiderSilk/data/pre_filtering/all_annotated.csv" 
properties_path = 'C:/Users/46705/Documents/SpiderSilk/data/raw_data/mechanical_properties.csv'


df = pd.read_csv(all_path, sep=',')


def make_bool(df, mech_prop, families): 
       spidroin_conts_path = "C:/Users/46705/Documents/SpiderSilk/data/pre_filtering/df_grouped.csv"

       spidroin_conts = pd.read_csv(spidroin_conts_path, sep= ",")

       spidroins = spidroin_conts[['ID']].copy()  
       spidroins = pd.concat([spidroins, spidroin_conts.iloc[:, 12:-19]], axis=1)  

       properties_df = pd.read_csv(properties_path, sep=",")

       df_mer = pd.merge(df, spidroins, on='ID', how='left')

       df = pd.merge(df_mer, properties_df, on='idv_id', how='left')

       metadata = df.iloc[:3]
       data = df.iloc[3:]

       #######################
       filtered_df = data[data["Family"].isin(families)].dropna(subset=mech_prop)

       data_df_filled = filtered_df.iloc[:, 7:-29].fillna(0)

       df_boolean = data_df_filled.map(lambda x: 1 if x != 0 else 0)

       filtered_df.iloc[:, 7: -29] = df_boolean
       ##############3

       common_columns = metadata.columns.intersection(filtered_df.columns)
       concatenated_df_boolean = pd.concat([metadata[common_columns], filtered_df], axis=0)


       return concatenated_df_boolean

fam =  ["Araneidae", "Salticidae", "Tetragnathidae","Theridiidae" ]
# fam =  ["Araneidae"]
bool_df = make_bool(df, ['toughness',"young's_modulus", 'tensile_strength', 'strain_at_break'], fam)


bool_df.to_csv("boolean_df_multip_sp.csv", index= False)
 

 









































# def print_file_size(file_path):
#     """Prints the size of a file."""
#     if os.path.exists(file_path):
#         file_size = os.path.getsize(file_path)
#         print(f"File size of {file_path}: {file_size} bytes")
#     else:
#         print(f"File {file_path} does not exist.")

# # Example usage:
# # Assuming 'file_path' is the path to your file
# # and you want to print memory usage before and after loading the file.

# def get_dataframe_memory_usage(df):
#     """Returns the memory usage of a Pandas DataFrame."""
#     memory = df.memory_usage(deep=True).sum() / (1024 ** 2)  # Convert bytes to megabytes
#     return memory

# # Filter columns where the count of empty lists exceeds the threshold
# # Calculate the threshold for 20% of the rows
# def remove_OG(df, abundance): 
#     frac = (100 - abundance) / 100
#     threshold = len(df.columns) * frac
#     empty_row_counts = df.isna().sum(axis=1)  # Count NaN values in each row
#     rows_to_drop_empty_cell = empty_row_counts[empty_row_counts > threshold].index

#     # Drop the filtered rows
#     df_filtered_empty_cell = df.drop(index=rows_to_drop_empty_cell)

#     # Print the filtered DataFrame
#     return df_filtered_empty_cell

# # Example usage:
# # Assuming 'df' is your DataFrame
# # For example:
# # df = pd.read_csv('your_file.csv')

# orthoG = "C:/Users/46705/Documents/skräp_spider_silk/annotation/Orthogroups_ls.tsv"
# # reads = "C:/Users/46705/Documents/SpiderSilk/big_data/full_body_rc"
# file_connectoion = "C:/Users/46705/Documents/skräp_spider_silk/data/raw_data/S1-S4/data_s1.csv"

# # orthoG = "/proj/naiss2023-6-13/miguel_analysis/orthol/orthofinder_mech_240125/Results_Jan25/WorkingDirectory/OrthoFinder/Results_Feb20/Orthogroups/Orthogroups.tsv"
# # reads = "/crex/proj/uppstore2019013/nobackup/private/data/bulkRNAseq/kSpiders/kallisto_quantification/individual/"
# # file_connectoion = "/proj/naiss2023-6-13/erika_analysis/data_s1.csv"

# abundance_filename = "abundance.tsv"
# info_filename = "run_info.json" 

# # Read orthofile
# df_orthogroups = pd.read_csv(orthoG, sep='\t')
# import sys

# # Read the CSV file with dtype options for the identifier column
# df_connections = pd.read_csv(file_connectoion, dtype=str, sep=";")

# # Choose to include full or not
# sub_sub = df_orthogroups.iloc[:10, :]
# # sub_sub = df_orthogroups

# # Separate genes --> make long df
# melted_df = pd.melt(sub_sub, id_vars='Orthogroup', var_name='Species', value_name='Genes')
# melted_df['Genes'] = melted_df['Genes'].str.split(', ')

# # Explode the lists in the "Genes" column into separate rows
# exploded_df = melted_df.explode('Genes')
# exploded_df['ID'] = exploded_df['Species'].str.replace('_longest_ORFs_aa$', '', regex=True)

# # Add metadata and DRR for read data and orthologs to be connected.
# df_DRR_species = df_connections[["ID", 'DRR', 'species', 'Genus', 'Family']]
# df_DRR_species["idv_id"] = df_connections['ID'].str.extract(r'(\d+)-W', expand=False).astype(str)

# merged_df = pd.merge(exploded_df, df_DRR_species,  on='ID', how='left')

# # ANNOTATION

# # path_info  = "/crex/proj/uppstore2019013/miguel_analysis/bridgeSpider/Lsc.v1.1.genes.AA.signalP.lascID.info.csv"
# path_info  = "C:/Users/46705/Documents/skräp_spider_silk/annotation/Lsc.v1.1.genes.AA.signalP.lascID.info.csv"
# # Read orthofile
# df_info = pd.read_csv(path_info, sep=',')

# df_info["Group"].unique()
# sp = len(df_info[df_info["Group"] == "Spidroin"])
# print(f"nr of Spidroin genes: {sp} ")
# spiCE = len(df_info[df_info["Group"] == "SpiCE"])
# print(f"nr of SpiCE genes: {spiCE} ")
# part_df = sub_sub[["Orthogroup", "Lsc.v1.1.genes.AA.signalP.lascID"]]

# part_df=part_df.dropna()
# melted_ann = pd.melt(part_df, id_vars='Orthogroup', var_name='Lsc.v1.1.genes.AA.signalP.lascID', value_name='Genes')
# melted_ann['Genes'] = melted_ann['Genes'].str.split(', ')

# exploded_df = melted_ann.explode('Genes')

# full_ortho_ann = pd.merge(exploded_df, df_info, how='left', left_on='Genes', right_on='lascID')

# # grouped_df = full_ortho_ann.groupby('Orthogroup').agg({'product': list, 'geneName': list,'Group': lambda x: list(x.unique())}).reset_index()
# def join_with_nan(x):
#     non_null_values = x.dropna()  # Drop NaN values
#     if len(non_null_values) > 0:  # If there are non-null values
#         return ','.join(non_null_values.astype(str).unique())  # Join non-null unique values
#     else:
#         return np.nan  # Return NaN if all values are NaN
    
# # Group by 'Orthogroup' and aggregate columns
# grouped_df = full_ortho_ann.groupby('Orthogroup').agg({'product': join_with_nan,
#                                                        'geneName': join_with_nan,
#                                                        'Group': join_with_nan}).reset_index()


# df_merged_anno = pd.merge(merged_df, grouped_df, on="Orthogroup", how="left")

# # Make pivot table and add CPM values
# # Fill NaN values with a placeholder value before pivoting
# filled_df = df_merged_anno.fillna({'Group': 'NaN', 'geneName': 'NaN', 'product': 'NaN'})
# filled_df

# # Pivot the DataFrame and fill NaN values with 0
# pivot_df = filled_df.pivot_table(
#     index=['ID', 'idv_id', 'DRR', 'species', 'Genus', 'Family'], 
#     columns=['Orthogroup', 'Group', "geneName", "product"], 
#     values='Tpm', 
#     aggfunc='size'
# ).fillna(0)

# # Replace the placeholder values back to NaN
# pivot_df.replace('NaN', np.nan, inplace=True)

# # pivot_df_non_spidroin.reset_index(inplace=True)
# # pivot_df_spidroin.reset_index(inplace=True)

# # pivot_df.to_csv("nr_genes_full_df.csv", index=True)



