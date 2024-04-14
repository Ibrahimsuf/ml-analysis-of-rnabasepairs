import requests
import csv
import os
import argparse
import pandas as pd
def getdssrstacks(pdb_id, output_path, url = None):
    url = f"http://skmatic.x3dna.org/pdb/{pdb_id.lower()}/{pdb_id.lower()}.json" if not url else url
    data = requests.get(url).json()
    with open(f"{output_path}/{pdb_id}.csv", "w") as f:
        writer = csv.writer(f)
        writer.writerow(["BasePairType", "nt1", "nt2"])
        for stack in data["stacks"]:
            for i in range(stack["num_nts"] - 1):
                nt1, nt2 = stack["nts_long"].replace("/", "").split(",")[i:i+2]
                base_pair_type = stack["nts_short"][i:i+2]
                writer.writerow([base_pair_type, nt1, nt2])
        if "helices" not in data:
            return
        for helix in data["helices"]:
            for pair in helix["pairs"]:
                nt1 = pair["nt1"].replace("/", "")
                nt2 = pair["nt2"].replace("/", "")
                base_pair_type = pair["bp"].replace("-", "")
                writer.writerow([base_pair_type, nt1, nt2])


def main(pdbs_list_file, output_path):
    pdb_list = pd.read_csv(pdbs_list_file)["PDB_ID"].unique().tolist()
    # with open("/Users/ibrahims/Downloads/pdb_ids_3_5.csv") as f:
    #     # first line of the file is the header
    #     next(f)
    #     for line in f:
    #         pdb_list.append(line.strip())
    
    os.makedirs(output_path, exist_ok=True)
    
    for pdb_id in pdb_list:
        print(pdb_id)
        getdssrstacks(pdb_id, output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("pdbs_list_file")
    parser.add_argument("output_path")
    args = parser.parse_args()
    main(args.pdbs_list_file, args.output_path)