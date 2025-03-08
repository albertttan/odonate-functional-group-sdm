import json
import pandas as pd

with open("../data/odonata/selected_species.json") as file:
    selected_species = json.load(file)
data = []

for species in selected_species:
    with open(f"output/maxent_files/{species}/2020_final.json") as file:
        data.append(json.load(file))
    for year in ["2050", "2075", "2100"]: 
        for ssp in ["ssp1", "ssp2", "ssp5"]: 
            with open(f"output/maxent_files/{species}/{year}_{ssp}.json") as file:
                data.append(json.load(file))

df = pd.DataFrame(data).to_csv("output/results.csv", index=False)
