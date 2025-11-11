from Bio import Entrez
import pandas as pd
from tqdm import tqdm

Entrez.email = "288847@student.pwr.edu.pl"

accessions = [
    'NC_005089.1', 'NC_012374.1', 'NC_000884.1', 'NC_037509.1', 'NC_034314.1', 
    'NC_027684.1', 'NC_027683.1', 'NC_048490.1', 'NC_018367.1', 'NC_025902.1', 
    'NC_025316.1', 'NC_023780.1', 'NC_053822.1', 'NC_031802.1', 'NC_001913.1', 
    'NC_011120.1', 'NC_011137.1', 'NC_012920.1', 'NC_002082.1', 'NC_014042.1', 
    'NC_001643.1', 'Y18001.1', 'NC_020009.2', 'NC_002083.1', 'NC_014047.1', 
    'MT711860.1', 'NC_028306.1', 'NC_001700.1', 'OR095102.1', 'OR095103.1', 
    'NC_027083.1', 'NC_022842.1', 'NC_010642.1', 'OR777682.1', 'MW257216.1', 
    'NC_016470.1', 'NC_026529.1', 'NC_008434.1', 'AY729880.1', 'NC_067757.1', 
    'NC_009686.1', 'NC_013700.1', 'NC_028427.1', 'NC_036369.1', 'NC_008093.1', 
    'NC_013445.1', 'NC_019591.1', 'NC_019590.1', 'NC_019589.1', 'NC_019588.1', 
    'NC_064558.1', 'NC_012059.1', 'NC_012058.1', 'NC_019577.1', 'NC_019578.2', 
    'NC_083222.1', 'NC_083094.1', 'NC_012062.1', 'NC_012051.1', 'NC_012053.1', 
    'NC_012061.1', 'NC_012057.1', 'NC_083148.1', 'OR490498.1', 'NC_082320.1', 
    'PQ064114.1', 'NC_057092.1', 'NC_056111.1', 'NC_041638.1', 'NC_041160.1', 
    'NC_027237.1', 'NC_034227.1', 'NC_029849.1', 'NC_029422.1', 'NC_029346.1', 
    'NC_029191.1', 'MZ457524.1', 'NC_007596.2', 'NC_015529.1', 'NC_005129.2', 
    'NC_000934.1', 'OL628830.1', 'JN673263.1', 'MK360903.1', 'NC_001794.1', 
    'NC_061372.1', 'NC_069653.1', 'NC_057520.1', 'NC_057519.1', 'NC_008133.1', 
    'NC_018788.1', 'NC_039717.1', 'NC_011944.1', 'NC_008447.1', 'NC_008145.1', 
    'NC_007631.1', 'FJ515782.1', 'KY996500.1'
]

GROUP_RULES = {
    "Kotowate": ["Felidae"],
    "Psowate": ["Canidae"],
    "Słoniowate": ["Proboscidea", "Elephantidae"],
    "Mroczkowate": ["Chiroptera", "Vespertilionidae"],
    "Naczelne": ["Primates"],
    "Gryzonie": ["Rodentia"],
    "Walenie": ["Cetacea", "Artiodactyla"],
    "Torbacze": ["Marsupialia", "Diprotodontia"]
}

LABEL_TO_INT = {
    "Kotowate": 0,
    "Psowate": 1,
    "Słoniowate": 2,
    "Mroczkowate": 3,
    "Naczelne": 4,
    "Gryzonie": 5,
    "Walenie": 6,
    "Torbacze": 7,
    "inne": 8
}

def classify_group(taxonomy_list):
    for group, keys in GROUP_RULES.items():
        for key in keys:
            if key in taxonomy_list:
                return group
    return "inne"

def fetch_taxonomy(acc):
    try:
        handle = Entrez.efetch(db="nuccore", id=acc, retmode="xml")
        record = Entrez.read(handle)[0]
        handle.close()

        sci_name = record["GBSeq_organism"]
        taxonomy_list = record["GBSeq_taxonomy"].split("; ")

        group = classify_group(taxonomy_list)

        return sci_name, "; ".join(taxonomy_list), group

    except:
        return None, None, "inne"

results = []
for acc in tqdm(accessions, desc="NCBI metadata download"):
    sci, tax, grp = fetch_taxonomy(acc)
    results.append([acc, sci, tax, grp])

df = pd.DataFrame(results, columns=["Accession", "Scientific_name", "Taxonomy", "Group"])
df["Label_ID"] = df["Group"].map(LABEL_TO_INT)

df.to_csv("taxonomy_labels_FIXED.csv", index=False)
print("✅ Zapisano taxonomy_labels_FIXED.csv")
print(df[df["Group"] != "inne"].head())
