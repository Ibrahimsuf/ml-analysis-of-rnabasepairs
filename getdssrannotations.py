import json
import pandas as pd
import os
import requests



def main():
    filenames  = os.listdir("test_pdbs2")
    missing_files = []
    for pdb_id in filenames:
        try:
            get_dssr_annotations(pdb_id.replace(".pdb", "").replace(".cif", ""), "dssrpairs_test_pdbs2")
        except:
            missing_files.append(pdb_id)

    print(missing_files)


def get_dssr_annotations(pdb_id, output_dir):
    print(f"http://skmatic.x3dna.org/pdb/{pdb_id}/{pdb_id}.json")
    annotations_json = requests.get(f"http://skmatic.x3dna.org/pdb/{pdb_id.lower()}/{pdb_id.lower()}.json").text
    try:
        pd.DataFrame(json.loads(annotations_json)["pairs"]).to_csv(f"{output_dir}/{pdb_id}.csv")
    except:
        print(f"No DSSR Pairs for {pdb_id}")
        pd.DataFrame().to_csv(f"{output_dir}/{pdb_id}.csv")

main()