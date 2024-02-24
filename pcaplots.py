import json
import pandas as pd
from sklearn.decomposition import PCA
import numpy as np
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt


def make_pca(df, features, basepair_type, log_plot = False, resolution = None, features_to_use = "None", plot_stacks = True):
    X = df.loc[:, features].values
    X = StandardScaler().fit_transform(X)

    #perform pca
    pca = PCA()
    principalComponents = pca.fit_transform(X)
    principalDf = pd.DataFrame(data = principalComponents)

    if plot_stacks:
        finalDf = pd.concat([principalDf, df[['BasePair', "BaseStack"]]], axis = 1)
    else:
        finalDf = pd.concat([principalDf, df['BasePair']], axis = 1)

    if log_plot:
        fig, axes = plt.subplots(3, 1, figsize=(10, 30))
    else:
        fig, axes = plt.subplots(2, 1, figsize=(10, 20))

    # Plot 1: PCA 1,2 v Base Pair type
    ax1 = axes[0]
    ax1.set_xlabel('Principal Component 1', fontsize=15)
    ax1.set_ylabel('Principal Component 2', fontsize=15)
    ax1.set_title(f'2 component PCA for {basepair_type}', fontsize=20)

    if plot_stacks:
        targets = ["Non Base Pair/Non Base Stack", "Base Pair", "Base Stack"]
        colors = ["r", "b", "g", "y"]
        for target, color in zip(targets, colors):
            if target == "Base Stack":
                indicesToKeep = (finalDf['BaseStack'] == True)
            elif target == "Non Base Pair/Non Base Stack":
                indicesToKeep = (finalDf['BasePair'] == False) & (finalDf['BaseStack'] == False)
            else:
                indicesToKeep = (finalDf['BasePair'] == True)
            ax1.scatter(finalDf.loc[indicesToKeep, 0], finalDf.loc[indicesToKeep, 1], c=color, s=50)
        ax1.legend(targets)
        ax1.grid()
    else:
        targets = ["Non Base Pair", "Base Pair"]
        colors = ["r", "b"]
        for target, color in zip(targets, colors):
            if target == "Base Pair":
                indicesToKeep = (finalDf['BasePair'] == True)
            else:
                indicesToKeep = (finalDf['BasePair'] == False)
            ax1.scatter(finalDf.loc[indicesToKeep, 0], finalDf.loc[indicesToKeep, 1], c=color, s=50)
        ax1.legend(targets)
        ax1.grid()

    # Plot 2: Full Scree Plot
    if log_plot:
        ax2 = axes[1]
        PC_values = np.arange(pca.n_components_) + 1
        ax2.plot(PC_values, pca.explained_variance_ratio_, 'ro-', linewidth=2)
        ax2.set_title(f'Full Scree Plot {basepair_type}')
        ax2.set_xlabel('Principal Component')
        ax2.set_ylabel('Proportion of Variance Explained')
        ax2.set_xscale("log")

    # Plot 3: Top 10 Scree Plot
    ax3 = axes[2] if log_plot else axes[1]
    PC_values = np.arange(10) if pca.explained_variance_ratio_.shape[0] > 10 else np.arange(pca.explained_variance_ratio_.shape[0])
    ax3.plot(PC_values, pca.explained_variance_ratio_[0:10], 'ro-', linewidth=2)
    ax3.set_title(f'Top 10 Scree Plot for {basepair_type}')
    ax3.set_xlabel('Principal Component')
    ax3.set_ylabel('Proportion of Variance Explained')
    annotations = f"The first 3 pca components explain {(pca.explained_variance_ratio_[0:2].sum() * 100):0.2f}% of the variance"
    annotations += f"\nThe first 6 pca components explain {(pca.explained_variance_ratio_[0:5].sum() * 100):0.2f}% of the variance" if features_to_use == "all" else ""
    plt.text(0, 0, annotations)
    plt.savefig(f"PCAPlots/{basepair_type}_{resolution}_pca_{features_to_use}_features.png")



def plot_pcas(resolution, features_to_use, log_plot = False, ):
    basepair_types = ["AA", "AC", "AG", "AU", "CC", "CG", "CU", "GG", "GU", "UU"]
    for basepair_type in basepair_types:
        df = pd.read_csv(f"{resolution}/{basepair_type}_stacks.csv")
        df["BasePair"].fillna(False, inplace=True)
        df["BaseStack"].fillna(False, inplace=True)

        if features_to_use == "all":
            features = [feature for feature in list(df.columns) if feature not in ['pdb_id', "nt1", "nt2", "BasePair", "BaseStack"]]
        elif features_to_use == "top 5":
            features = pd.read_csv(f"feature_importances/feature_importances_{basepair_type}.csv").sort_values("importance", ascending=True)[0:5]["Unnamed: 0"].values

        make_pca(df, features, basepair_type, resolution=resolution, features_to_use=features_to_use, log_plot=log_plot)        
        




def main():
    # plot_pcas("3_5", "top 5")
    plot_pcas("0_3", "top 5")

if __name__ == "__main__":
    main()