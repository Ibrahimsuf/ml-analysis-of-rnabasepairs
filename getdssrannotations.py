import json
import pandas as pd
import os
import requests



def main():
    filenames  = os.listdir("test_pdbs")
    missing_files = []
    for pdb_id in filenames:
        try:
            get_dssr_annotations(pdb_id.replace(".pdb", "").replace(".cif", ""))
        except:
            missing_files.append(pdb_id)

    print(missing_files)


def get_dssr_annotations(pdb_id):
    print(f"http://skmatic.x3dna.org/pdb/{pdb_id}/{pdb_id}.json")
    annotations_json = requests.get(f"http://skmatic.x3dna.org/pdb/{pdb_id.lower()}/{pdb_id.lower()}.json").text
    try:
        pd.DataFrame(json.loads(annotations_json)["pairs"]).to_csv(f"dssrstacks_pairs_test_pdbs/{pdb_id}.csv")
    except:
        print(f"No DSSR Pairs for {pdb_id}")
        pd.DataFrame().to_csv(f"dssr_annotations/{pdb_id}.csv")

main()