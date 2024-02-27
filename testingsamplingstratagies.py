import pandas as pd
from imblearn.under_sampling import RandomUnderSampler
from sklearn.decomposition import PCA
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from imblearn.under_sampling import NearMiss
from imblearn.under_sampling import CondensedNearestNeighbour
from imblearn.under_sampling import TomekLinks
from imblearn.combine import SMOTETomek
from imblearn.over_sampling import ADASYN
import matplotlib.pyplot as plt
from XGBModel import XGBModel
import json
from tqdm import tqdm
from numpyencoder import NumpyEncoder

def main(outputfile="metrics/sampling_strategies_accuracies.json"):
  accuracies = {}
  for basepairtype in tqdm(["AA", "AC", "AG", "AU", "CC", "CG", "CU", "GG", "GU", "UU"]):
    accuracies[basepairtype] = get_sampling_stratagies_accuracies(basepairtype)
  with open(outputfile, "w") as f:
    json.dump(accuracies, f, cls=NumpyEncoder)

  


def get_sampling_stratagies_accuracies(basepairtype):
  data = pd.read_csv(f"0_3/{basepairtype}_stacks.csv")
  data_3_5 = pd.read_csv(f"3_5/{basepairtype}_stacks.csv")
  train, test = create_test_set_with_even_pairs_and_non_pairs(data)
  test_3_5 = create_even_test_set(data_3_5)
  
  del data
  del data_3_5

  features = [feature for feature in train.columns if feature not in ["pdb_id", "nt1", "nt2", "BasePair", 'BaseStack']]
  pca = make_pipeline(StandardScaler(), PCA(n_components=2))
  pca.fit(train[features])

  sampling_stratagies = {
    "RandomUnderSampler1": RandomUnderSampler(random_state=42, replacement=False),
    "RandomUnderSampler2": RandomUnderSampler(random_state=43, replacement=False),
    "NearMiss": NearMiss(),
    "CondensedNearestNeighbour": CondensedNearestNeighbour(),
    "TomekLinks": TomekLinks(),
    "SMOTETomek": SMOTETomek(),
    "ADASYN": ADASYN()
  }
  sampling_strategy_metrics = {}
  for sampling_method_name, method in tqdm(sampling_stratagies.items()):
    X_resampled, y_resampled = method.fit_resample(train[features], train["BasePair"])
    X_resampled.reset_index(drop=True, inplace=True)
    y_resampled.reset_index(drop=True, inplace=True)
    pca_plot = make_pca(train, pd.concat([X_resampled, y_resampled], axis=1), features, sampling_method_name, pca)
    pca_plot.savefig(f"PCAplots/{sampling_method_name}_{basepairtype}.png")


    metrics_0_3 = train_and_test_resampling(X_resampled, y_resampled, test, features, "0_3")
    metrics_3_5 = train_and_test_resampling(X_resampled, y_resampled, test_3_5, features, "3_5")
    metrics = {**metrics_0_3, **metrics_3_5}
    sampling_strategy_metrics[sampling_method_name] = metrics

  return sampling_strategy_metrics



def train_and_test_resampling(X_resampled, y_resampled, test, features, resolution):
    model = XGBModel(pd.concat([X_resampled, y_resampled], axis=1), target = "BasePair", split = None)
    model.set_test_data(test[features], test["BasePair"])
    model.fit(100)
    classification_report = model.get_classifcation_report(output_dict=True)
    tn, fp, fn, tp = model.get_confusion_matrix()

    metrics = {}
    metrics[f"Train True Examples"] = y_resampled.sum()
    metrics[f"Train False Examples"] = len(y_resampled) - y_resampled.sum()
    metrics[f"True Precision_{resolution}"] = (classification_report["True"]["precision"])
    metrics[f"True Recall_{resolution}"] = (classification_report["True"]["recall"])
    metrics[f"True F1_{resolution}"] = (classification_report["True"]["f1-score"])
    metrics[f"True Support_{resolution}"] = (classification_report["True"]["support"])
    metrics[f"TP_{resolution}"] = (tp)
    metrics[f"FP_{resolution}"] = (fp)
    metrics[f"FN_{resolution}"] = (fn)
    metrics[f"TN_{resolution}"] = (tn)
    metrics[f"False Precision_{resolution}"] = classification_report["False"]["precision"]
    metrics[f"False Recall_{resolution}"] = classification_report["False"]["recall"]
    metrics[f"False F1_{resolution}"] = classification_report["False"]["f1-score"]
    metrics[f"False Support_{resolution}"] = classification_report["False"]["support"]
    return metrics

def create_test_set_with_even_pairs_and_non_pairs(data):
    # Drop column if it exists
    if "Unnamed: 0" in data.columns:
        data.drop(["Unnamed: 0"], axis=1, inplace=True)

    data["BasePair"].fillna(False, inplace=True)

    pairs = data[data["BasePair"]]
    nonpairs = data[~data["BasePair"]]

    num_test_examples_of_each_class = int(pairs["BasePair"].sum() * 0.2)

    # Sample directly from the original data
    test_pairs = pairs.sample(num_test_examples_of_each_class, random_state=42)
    test_nonpairs = nonpairs.sample(num_test_examples_of_each_class, random_state=42)

    # Use iloc for faster indexing
    train_pairs = pairs.drop(test_pairs.index)
    train_nonpairs = nonpairs.drop(test_nonpairs.index)

    # Concatenate dataframes once
    train = pd.concat([train_pairs, train_nonpairs]).sample(frac=1).reset_index(drop=True)
    test = pd.concat([test_pairs, test_nonpairs]).sample(frac=1).reset_index(drop=True)

    return train, test

def create_even_test_set(data):
  if "Unnamed: 0" in data.columns:
    data.drop(["Unnamed: 0"], axis=1, inplace=True)
  data["BasePair"].fillna(False, inplace=True)
  sample_size = int(data["BasePair"].sum() * 0.2)
  return data.groupby('BasePair').apply(lambda x: x.sample(sample_size))



def make_pca(real_data, resampled_data, features, sampling_method, pca):
  fig, axes = plt.subplots(2,2,figsize=(10, 10))
  ax = axes[1, 0]
  ax.set_xlabel('Principal Component 1', fontsize=15)
  ax.set_ylabel('Principal Component 2', fontsize=15)
  for data, color_pairs, plot in zip([real_data, resampled_data], [('r', 'b'), ('g', 'm')], [axes[0,0], axes[0,1]]):
    principalComponents = pca.transform(data[features])
    principalDf = pd.DataFrame(data=principalComponents)

    finalDf = pd.concat([principalDf, data['BasePair']], axis=1)
    targets = ["Non Base Pair", "Base Pair"]
    for target, color in zip(targets, color_pairs):
        if target == "Base Pair":
            indicesToKeep = finalDf['BasePair'] == True
        else:
            indicesToKeep = finalDf['BasePair'] == False

        ax.scatter(finalDf.loc[indicesToKeep, 0], finalDf.loc[indicesToKeep, 1], c=color, s=50)
        plot.scatter(finalDf.loc[indicesToKeep, 0], finalDf.loc[indicesToKeep, 1], c=color, s=50)
        plot.legend(targets, prop={'size': 15})


  ax.legend(["Real Non Base Pair", "Real Base Pair", "Resampled Non Base Pair", "Resampled Base Pair"])
  fig.suptitle(f"PCA for orginal and resampled data with {sampling_method}", fontsize=20)
  axes[0,0].title.set_text("Real Data")
  axes[0,1].title.set_text("Resampled Data")
  axes[1,0].title.set_text("Combined Real and Resampled Data")
  axes[1,1].remove()
  return fig



if __name__ == "__main__":
  main()