import pandas as pd

resolution = "3_5"
split_percentage = 0.2
basepairtypes = ["AA", "AC", "AG", "AU", "CC", "CG", "CU", "GG", "GU", "UU"]

for basepairtype in basepairtypes:
  df = pd.read_csv(f"{resolution}/{basepairtype}_stacks.csv")
  df["BasePair"].fillna(False, inplace=True)
  

  num_test_examples_of_each_class = int(df["BasePair"].sum() * split_percentage)
  pairs = df[df["BasePair"]]
  nonpairs = df[~df["BasePair"]]

  test_pairs = pairs.sample(num_test_examples_of_each_class)
  train_pairs=  pairs.drop(test_pairs.index)

  test_nonpairs = nonpairs.sample(num_test_examples_of_each_class)
  train_nonpairs = nonpairs.drop(test_nonpairs.index)

  train = pd.concat([train_pairs, train_nonpairs]).sample(frac=1).reset_index(drop=True)  
  test = pd.concat([test_pairs, test_nonpairs]).sample(frac=1).reset_index(drop=True)

  train.to_csv(f"{resolution}/train/{basepairtype}_pairs.csv", index = False)
  test.to_csv(f"{resolution}/test/{basepairtype}_pairs.csv", index = False)

