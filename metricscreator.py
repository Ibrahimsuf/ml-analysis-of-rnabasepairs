from collections import defaultdict
from XGBModel import XGBModel
import xgboost as xgb
import pandas as pd
import numpy as np
from sklearn.metrics import confusion_matrix, classification_report
class MetricDataFrameCreator:
    def __init__(self, params, features, target = "BasePair"):
        self.params = params
        if not features in ["all", "top 5"]:
            raise ValueError("features must be 'all' or 'top 5'")
        if not target in ["BasePair", "BaseStack"]:
            raise ValueError("target must be 'BasePair' or 'BaseStack'")
        self.target = target
        self.features = features
        self.basepairtypes = ["AA", "AC", "AG", "AU", "CC", "CG", "CU", "GG", "GU", "UU"]
        self.metrics = defaultdict(dict)
        self.models = {}
    def train_model(self, basepairtype):
        data_df = pd.read_csv(f"0_3/{basepairtype}_stacks.csv")

        data_df[self.target].fillna(False, inplace=True)

        if self.features == "all":
            self.models[basepairtype] = XGBModel(data_df, target = self.target, parms = self.params)
        elif self.features == "top 5":
            features = pd.read_csv(f"feature_importances/feature_importances_{basepairtype}.csv").sort_values("importance", ascending=True)[0:5]["Unnamed: 0"].values
            self.models[basepairtype] = XGBModel(data_df[list(features) + [self.target]], target = self.target, parms = self.params)
        
        self.models[basepairtype].fit(100)

    def get_metrics(self, basepairtype, resolution, test_pdbs_dir = None):
        model = self.models[basepairtype]
        if resolution == "0_3":
            tn, fp, fn, tp = model.get_confusion_matrix()
            claissification_report_dict = model.get_classifcation_report(output_dict=True)
        elif resolution == "3_5":
            data_test = pd.read_csv(f"3_5/{basepairtype}_stacks.csv")
            features = model.bst.feature_names

            dtest = xgb.DMatrix(data = data_test[features])
            y_preds = model.bst.predict(dtest).round()
            y_true = data_test[self.target].fillna(False)
            tn, fp, fn, tp = confusion_matrix(y_true, y_preds).ravel()
            claissification_report_dict = classification_report(y_true, y_preds, output_dict=True)
        elif resolution == "test_pdbs":
            data_test = pd.read_csv(f"{test_pdbs_dir}/{basepairtype}.csv")
            features = model.bst.feature_names
            dtest = xgb.DMatrix(data = data_test[features])
            y_preds = model.bst.predict(dtest).round()
            y_true = data_test[self.target].fillna(False)
            tn, fp, fn, tp = confusion_matrix(y_true, y_preds, labels=[False, True]).ravel()
            claissification_report_dict = classification_report(y_true, y_preds, output_dict=True, labels=[False, True], zero_division=np.nan)

        record = {}
        record[f"True Precision_{resolution}"] = (claissification_report_dict["True"]["precision"])
        record[f"True Recall_{resolution}"] = (claissification_report_dict["True"]["recall"])
        record[f"True F1_{resolution}"] = (claissification_report_dict["True"]["f1-score"])
        record[f"True Support_{resolution}"] = (claissification_report_dict["True"]["support"])
        record[f"TP_{resolution}"] = (tp)
        record[f"FP_{resolution}"] = (fp)
        record[f"FN_{resolution}"] = (fn)
        record[f"TN_{resolution}"] = (tn)
        record[f"False Precision_{resolution}"] = claissification_report_dict["False"]["precision"]
        record[f"False Recall_{resolution}"] = claissification_report_dict["False"]["recall"]
        record[f"False F1_{resolution}"] = claissification_report_dict["False"]["f1-score"]
        record[f"False Support_{resolution}"] = claissification_report_dict["False"]["support"]
        
        if resolution == "0_3":
            data_df = pd.read_csv(f"0_3/{basepairtype}_stacks.csv")
            data_df[self.target].fillna(False, inplace=True)
            record["% True_0_3"] = (data_df[self.target].fillna(False).mean())
        elif resolution == "3_5":
            record["% True_3_5"] = (data_test[self.target].fillna(False).mean())

        self.metrics[basepairtype].update(record)

    def get_metrics_for_all(self):
        for basepairtype in self.basepairtypes:
            self.train_model(basepairtype)
            self.get_metrics(basepairtype, "0_3")
            self.get_metrics(basepairtype, "3_5")
        return self.metrics
    def get_metrics_test_pdbs(self, test_pdbs_dir):
        for basepairtype in self.basepairtypes:
            self.train_model(basepairtype)
            self.get_metrics(basepairtype, "test_pdbs", test_pdbs_dir)
        return self.metrics


def main():
    param_weighted = {'max_depth': 2, 'eta': 1, 'objective': 'binary:logistic', "scale_pos_weight": 100}
    param_weighted['nthread'] = 4
    param_weighted['eval_metric'] = ['logloss', "error", 'pre']
    # # metric_df_creator = MetricDataFrameCreator(param_weighted, "all")
    # metric_df_creator = MetricDataFrameCreator(param_weighted, "top 5")


    # metrics = metric_df_creator.get_metrics_for_all()
    # # pd.DataFrame(metrics).to_csv("WeightedParams_All_feautures_metrics.csv")
    # pd.DataFrame(metrics).to_csv("WeightedParams_Top5_feautures_metrics.csv")

    # metric_df_creator = MetricDataFrameCreator(param_weighted, "top 5", target = "BaseStack")
    # metrics = metric_df_creator.get_metrics_for_all()
    # pd.DataFrame(metrics).to_csv("WeightedParams_Top5_feautures_metrics_stacks.csv")

    # metric_df_creator_all = MetricDataFrameCreator(param_weighted, "all", target = "BaseStack")
    # metrics_all = metric_df_creator_all.get_metrics_for_all()
    # pd.DataFrame(metrics_all).to_csv("WeightedParams_All_feautures_metrics_stacks.csv")

    metric_df_creator_all = MetricDataFrameCreator(param_weighted, "all")
    metrics_all = metric_df_creator_all.get_metrics_test_pdbs("test_pdbs2")
    pd.DataFrame(metrics_all).to_csv("WeightedParams_All_feautures_metrics_test_pdbs2.csv")

    metric_df_creator_all = MetricDataFrameCreator(param_weighted, "top 5")
    metrics_all = metric_df_creator_all.get_metrics_test_pdbs("test_pdbs2")
    pd.DataFrame(metrics_all).to_csv("WeightedParams_Top5_feautures_metrics_test_pdbs2.csv")

if __name__ == "__main__":
    main()