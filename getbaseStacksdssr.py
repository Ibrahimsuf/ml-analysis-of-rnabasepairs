import requests
import csv
import os
def getdssrstacks(pdb_id, output_path):
    data = requests.get(f"http://skmatic.x3dna.org/pdb/{pdb_id.lower()}/{pdb_id.lower()}.json").json()
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


def main():
    pdb_list = []
    with open("/Users/ibrahims/Downloads/pdb_ids_3_5.csv") as f:
        # first line of the file is the header
        next(f)
        for line in f:
            pdb_list.append(line.strip())
    
    output_path = f"dssrstacks_3_5"
    os.makedirs(output_path, exist_ok=True)
    
    for pdb_id in pdb_list:
        print(pdb_id)
        getdssrstacks(pdb_id, "dssrstacks_3_5")


if __name__ == "__main__":
    main()