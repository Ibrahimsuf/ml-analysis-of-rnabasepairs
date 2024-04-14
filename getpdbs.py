import requests
import gzip
import shutil
import os
import pandas as pd
import argparse


def main(pdb_list_file, output_dir):
    pdb_ids = pd.read_csv(pdb_list_file)["PDB_ID"].values
    missing_files = []
    for pdb_id in pdb_ids:
        try:
            download_pbd(pdb_id, output_dir)
        except(Exception) as e:
            print(e)
            missing_files.append(pdb_id, output_dir)

    print(missing_files)

def download_pbd(pdb_id, output_dir):
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
    print(url)

    # Create the output folder if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_path = os.path.join(output_dir, pdb_id.upper() + ".cif")

    # Download the file
    response = requests.get(url, stream=True)
    with open(output_path + '.gz', 'wb') as f:
        shutil.copyfileobj(response.raw, f)

    # Unzip the file
    with gzip.open(output_path + '.gz', 'rb') as f_in:
        with open(output_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    # Remove the gzipped file if needed
    os.remove(output_path + '.gz')

    print(f"File downloaded and unzipped as {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb_list_file")
    parser.add_argument("output_dir")
    args = parser.parse_args()
    main(args.pdb_list_file, args.output_dir) 