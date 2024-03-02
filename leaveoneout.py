import pandas as pd
from XGBModel import XGBModel
from sklearn.model_selection import LeaveOneOut
import warnings
from sklearn.metrics import confusion_matrix, classification_report
import json
import numpy as np
from numpyencoder import NumpyEncoder

warnings.filterwarnings("ignore")

metrics = {"AA": {}, "AC": {}, "AG": {}, "AU": {}, "CC": {}, "CG": {}, "CU": {}, "GG": {}, "GU": {}, "UU": {}}


def run_analysis_for_basepairtype(basepairtype):
  data = pd.read_csv(f"0_3/{basepairtype}_stacks.csv")
  data["Class"] = data.apply(getclass, axis=1)

  class_num_sample_mapping = num_samples_for_each_class(data)
  samples_to_test_on = data.groupby("Class").apply(lambda x: x.sample(class_num_sample_mapping[x["Class"].iloc[0]])).sample(frac=1).index.droplevel(0).values
  data = pd.concat([data.iloc[samples_to_test_on], data.drop(samples_to_test_on)])
  cv = LeaveOneOut()

  count = 0
  y_true = []
  y_pred = []
  for train_ix, test_ix in cv.split(data):
    count += 1
    model = XGBModel(data.iloc[train_ix], target="Class", split=False)
    model.fit(100)
    y_pred.append(model.predict(data.iloc[test_ix]))
    y_true.append(data.iloc[test_ix]["Class"])

    if count % 10 == 0:
      print(f"{count} it/s")
    
    if count == 1_000:
      break
  
  metrics[basepairtype]["Classification Report"] = classification_report(y_true, y_pred, output_dict=True, labels=[0, 1, 2])
  confusion_matrix_for_basepairtype = confusion_matrix(y_true, y_pred, labels=[0, 1, 2])
  pd.DataFrame(confusion_matrix_for_basepairtype, columns=["Neither", "BasePair", "BaseStack"], index=["Neither", "BasePair", "BaseStack"]).to_csv(f"metrics/confusion_matrix/{basepairtype}.csv")





def getclass(row):
  if row["BasePair"] == True:
    return 1
  elif row["BaseStack"] == True:
    return 2
  else:
    return 0

def num_samples_for_each_class(df):
  num_total_samples_left = 1_000
  class_num_sample_mapping = {0: 0, 1: 0, 2: 0}
  
  class_num_sample_mapping[2] = min(df[df["Class"] == 2]["Class"].shape[0], num_total_samples_left // 3)
  num_total_samples_left -= class_num_sample_mapping[2]

  class_num_sample_mapping[1] = min(df[df["Class"] == 1]["Class"].shape[0], num_total_samples_left // 2)
  num_total_samples_left -= class_num_sample_mapping[1]

  class_num_sample_mapping[0] = num_total_samples_left
  return class_num_sample_mapping


def main():
  for basepairtype in ["CU", "GG", "GU", "UU"]:
    run_analysis_for_basepairtype(basepairtype)

  with open("metrics/leave_one_out_accuracies.json", "w") as f:
    json.dump(metrics, f, cls=NumpyEncoder)

if __name__ == "__main__":
  main()