import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


## Load data

df = pd.read_csv("data/presence.csv", sep="\t")
# Dataset already filtered by year and continent through GBIF
df_filtered = df.copy()
# df.head()


## Filter by temporal scope (year)

(min_year, max_year) = (2015, 2025)
df_filtered = df_filtered[
    (df_filtered["year"] >= min_year) & 
    (df_filtered["year"] <= max_year)
]
# df_filtered[["year"]].describe()
# 246,362 remaining records


## Filter by geographic scope (latitude and longitude)

(min_lat, max_lat, min_lon, max_lon) = (5, 85, -180, -50)
df_filtered = df_filtered[
    (df_filtered["decimalLatitude"] > min_lat) & 
    (df_filtered["decimalLatitude"] < max_lat) & 
    (df_filtered["decimalLongitude"] > min_lon) & 
    (df_filtered["decimalLongitude"] < max_lon)
]
# df_filtered[["decimalLatitude", "decimalLongitude"]].describe()
# 246,206 remaining records


## Filter by taxonomic scope (genus)

selected_genera = ["Argia", "Aeshna", "Calopteryx", "Hetaerina", 
                   "Libellula", "Sympetrum", "Plathemis"]
df_filtered = df_filtered[
    df_filtered['genus'].isin(selected_genera)
]
# df_filtered.groupby("genus").count()[["gbifID"]]
# df_filtered[["genus"]].describe()
# 79,859 remaining records


## Aggregate records into species; filter by number of records

species_counts_full = df_filtered.groupby(['genus', 'species']).size()
# species_counts_full.count()
# 122 remaining species

species_counts_filtered = species_counts_full[species_counts_full >= 30]
# species_counts_filtered.count()
# 84 remaining species


## Select records from at most 5 species per genera

selected_species = species_counts_filtered.groupby('genus').nlargest(5)\
    .reset_index(level=0, drop=True).index.get_level_values('species').tolist()

df_final = df_filtered[df_filtered['species'].isin(selected_species)]
species_counts_final = df_final.groupby("species").count()["gbifID"].sort_index()

# species_counts_final.sum()
# 58,900 selected records

# species_counts_final.count()
# 30 selected species

# species_counts_final
# See Table 1 in paper


## Save cleaned data and species names for SDM

df_final[["species", "decimalLatitude", "decimalLongitude"]]\
    .to_csv("output/cleaning/presence.csv", index=False)

with open("output/cleaning/selected_species.json", 'w') as file:
    json.dump(selected_species, file)

