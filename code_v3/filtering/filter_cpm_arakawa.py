import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

sp_path = "C:/Users/46705/Documents/SpiderSilk/data/pre_filtering/spidroins_cpm_annotated_lascID.csv" # OG and from the newest orthogroup with l.s included. however that species (bridgespider) can not be included due to no counts assigned. . 
non_sp_path = "C:/Users/46705/Documents/SpiderSilk/data/pre_filtering/non_spidroins_tpm_annotated_lascID.csv" 
all_path = "C:/Users/46705/Documents/SpiderSilk/data/pre_filtering/all_annotated_lascID.csv" 
properties_path = 'C:/Users/46705/Documents/SpiderSilk/data/raw_data/mechanical_properties.csv'

path_gland = "C:/Users/46705/Documents/SpiderSilk/data/pre_filtering/gland_og.csv"

path_gland_sac = "C:/Users/46705/Documents/SpiderSilk/data/pre_filtering/sac_gland_og.csv"

path_gland_tail = "C:/Users/46705/Documents/SpiderSilk/data/pre_filtering/tail_gland_og.csv"

path_tail_sac = "C:/Users/46705/Documents/SpiderSilk/data/pre_filtering/tail_sac_5_og_log.csv"

# df = pd.read_csv(all_path, sep=',')

# df_gland_og = pd.read_csv(path_gland, sep=',')
df_tail_sac_filtered = pd.read_csv(path_tail_sac, sep=',')

# df_sac = pd.read_csv(path_gland_sac, sep=',')
# df_tail = pd.read_csv(path_gland_tail, sep=',')

print("df read")
# df_small = pd.read_csv(sp_path, sep=',')

def initial_filtering(df, mech_prop, families, filtering_threshold, zero_prop_threshold, log=False): 
       spidroin_conts_path = "C:/Users/46705/Documents/SpiderSilk/data/pre_filtering/df_grouped.csv"

       spidroin_conts = pd.read_csv(spidroin_conts_path, sep= ",")

       spidroins = spidroin_conts[['ID']].copy()  # Copy the "ID" column
       spidroins = pd.concat([spidroins, spidroin_conts.iloc[:, 12:-19]], axis=1)  # Concatenate with desired columns from spidroin_conts
       
       properties_df = pd.read_csv(properties_path, sep=",")
       print(spidroins)
       print(properties_df)


       df_mer = pd.merge(df, spidroins, on='ID', how='left')
       df = pd.merge(df_mer, properties_df, on='idv_id', how='left')
       cols_to_drop_ = df.iloc[:, 7:-43].columns
       df = df.drop(columns=cols_to_drop_)

       # df = pd.merge(df, properties_df, on='idv_id', how='left')
       print(df["Family"].unique())

       metadata = df.iloc[:4]
       data = df.iloc[4:]

       filtered_df = data[data["Family"].isin(families)].dropna(subset=mech_prop)

       subset_early = filtered_df.iloc[:, 7:-29].apply(pd.to_numeric, errors='coerce').fillna(0)

       zero_prop_early = (subset_early == 0).mean()

       print(subset_early)


       sns.kdeplot((abs(zero_prop_early-1)))
       plt.title('')
       plt.axvline(x=zero_prop_threshold, color='r', linestyle='--')
       plt.xlabel('Expression coverage across species ALL')
       plt.ylabel('Frequency')
       plt.grid(True)
       plt.show()


       col_sums = filtered_df.iloc[:, 7:-29].astype(float).sum()
       # cols_to_drop = col_sums[col_sums == 0].index
       cols_to_drop = col_sums[col_sums < filtering_threshold].index
       print(len(cols_to_drop))
       filtered_df_zero_GG_removed = filtered_df.drop(columns=cols_to_drop)

       if log: 
              filtered_df_zero_GG_removed.iloc[:,7:-29] = filtered_df_zero_GG_removed.iloc[:,7:-29].astype(float).apply(lambda x: np.log2(x + 1)) 
       
       col_sums_2 = filtered_df_zero_GG_removed.iloc[:, 7:-29].astype(float).sum()

       cols_to_drop_2 = col_sums_2[col_sums_2 < filtering_threshold].index
       print(len(cols_to_drop_2))
       filtered_df_zero_GG_removed_2 = filtered_df_zero_GG_removed.drop(columns=cols_to_drop_2)


       common_columns = metadata.columns.intersection(filtered_df_zero_GG_removed_2.columns)
       concatenated_df = pd.concat([metadata[common_columns], filtered_df_zero_GG_removed_2], axis=0)

       plt.figure(figsize=(20, 12))
       # Plot histogram
       plt.hist(col_sums_2, bins=200, color='skyblue', edgecolor='black')
       plt.title('reads per Gene Group distribution')
       plt.axvline(x=filtering_threshold, color='r', linestyle='--')
       plt.xlabel('Sum of Columns')
       plt.ylabel('Frequency')
       plt.show()

       concatenated_df
       subset = concatenated_df.iloc[4 :, 7:-29]


       mask = subset.isna().any(axis=1)
       concatenated_df = concatenated_df.drop(concatenated_df.iloc[4:,:][mask].index)

#########################
       
       subset = concatenated_df.iloc[4:, 7:-29]

       zero_prop = (subset == 0).mean()
       columns_to_drop = zero_prop[zero_prop > zero_prop_threshold].index
       print(len(columns_to_drop))

       # Drop columns from main_df
       concatenated_df.drop(columns=columns_to_drop, inplace=True)

#########################
       sns.kdeplot((abs(zero_prop-1)))
       plt.title('')
       plt.axvline(x=zero_prop_threshold, color='r', linestyle='--')
       plt.xlabel('Expression coverage across species')
       plt.ylabel('Frequency')
       plt.grid(True)
       plt.show()


       return concatenated_df, abs(zero_prop_early-1)

fam =  ["Araneidae", "Nephilidae", "Tetragnathidae", "Theridiidae" ]
# fam =  ["Araneidae"]
# filtered_sac = initial_filtering(df_sac, ['toughness',"young's_modulus", 'tensile_strength', 'strain_at_break'], fam, 20, 0.9, log=True)

 
# filtered_tail = initial_filtering(df_tail, ['toughness',"young's_modulus", 'tensile_strength', 'strain_at_break'], fam, 20, 0.9, log=True)



filtered_arakawa, coverage_across_species = initial_filtering(df_tail_sac_filtered, ['toughness',"young's_modulus", 'tensile_strength', 'strain_at_break'], fam, 20, 0.5, log=True)


# filtered_arakawa.to_csv("C:/Users/46705/Documents/SpiderSilk/data/post_filtering/filtered_full_multisp_log_20_lascID.csv", index=False)
filtered_arakawa.to_csv("C:/Users/46705/Documents/SpiderSilk/data/post_filtering/arakawa_anno_log_multisp.csv", index=False)




