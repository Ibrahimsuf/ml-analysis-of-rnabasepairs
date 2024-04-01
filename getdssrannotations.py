import json
import pandas as pd
import os
import requests



def main():
    filenames  = pd.read_csv("3_5/UU_beaddistances.csv.gz")["pdb_id"].unique().tolist()
    output_dir = "dssrpairs3_5"
    missing_files = []
    for pdb_id in filenames:
        try:
            get_dssr_annotations(pdb_id.replace(".pdb", "").replace(".cif", ""), output_dir)
        except:
            missing_files.append(pdb_id)

    print(missing_files)


def get_dssr_annotations(pdb_id, output_dir, url = None):
    url = f"http://skmatic.x3dna.org/pdb/{pdb_id.lower()}/{pdb_id.lower()}.json" if url is None else url
    print(url)
    annotations_json = requests.get(url).text
    try:
        pd.DataFrame(json.loads(annotations_json)["pairs"]).to_csv(f"{output_dir}/{pdb_id}.csv")
    except:
        print(f"No DSSR Pairs for {pdb_id}")
        pd.DataFrame().to_csv(f"{output_dir}/{pdb_id}.csv")



if __name__ == "__main__":
   main()