from collections import defaultdict
from XGBModel import XGBModel
import xgboost as xgb
import pandas as pd
from sklearn.metrics import confusion_matrix, classification_report
class MetricDataFrameCreator:
    def __init__(self, params, features):
        self.params = params
        if not features in ["all", "top 5"]:
            raise ValueError("features must be 'all' or 'top 5'")
        self.features = features
        self.basepairtypes = ["AA", "AC", "AG", "AU", "CC", "CG", "CU", "GG", "GU", "UU"]
        self.metrics = defaultdict(dict)
        self.models = {}
    def train_model(self, basepairtype):
        data_df = pd.read_csv(f"0_3/{basepairtype}.csv")

        data_df["BasePair"].fillna(False, inplace=True)

        if self.features == "all":
            self.models[basepairtype] = XGBModel(self.params, data_df)
        elif self.features == "top 5":
            features = pd.read_csv(f"feature_importances_{basepairtype}.csv").sort_values("importance", ascending=True)[0:5]["Unnamed: 0"].values
            self.models[basepairtype] = XGBModel(self.params, data_df[list(features) + ["BasePair"]])
        
        self.models[basepairtype].fit(100)

    def get_metrics(self, basepairtype, resolution):
        model = self.models[basepairtype]
        if resolution == "0_3":
            tp, fp, fn, tn = model.get_confusion_matrix()
            claissification_report_dict = model.get_classifcation_report(output_dict=True)
        elif resolution == "3_5":
            data_test = pd.read_csv(f"3_5/{basepairtype}.csv")
            features = model.bst.feature_names

            dtest = xgb.DMatrix(data = data_test[features])
            y_preds = model.bst.predict(dtest).round()
            y_true = data_test["BasePair"].fillna(False)
            tn, fp, fn, tp = confusion_matrix(y_true, y_preds).ravel()
            claissification_report_dict = classification_report(y_true, y_preds, output_dict=True)
                

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
            data_df = pd.read_csv(f"0_3/{basepairtype}.csv")
            data_df["BasePair"].fillna(False, inplace=True)
            record["% True_0_3"] = (data_df["BasePair"].fillna(False).mean())
        elif resolution == "3_5":
            record["% True_3_5"] = (data_test["BasePair"].fillna(False).mean())

        self.metrics[basepairtype].update(record)

    def get_metrics_for_all(self):
        for basepairtype in self.basepairtypes:
            self.train_model(basepairtype)
            self.get_metrics(basepairtype, "0_3")
            self.get_metrics(basepairtype, "3_5")
        return self.metrics


def main():
    param_weighted = {'max_depth': 2, 'eta': 1, 'objective': 'binary:logistic', "scale_pos_weight": 100}
    param_weighted['nthread'] = 4
    param_weighted['eval_metric'] = ['logloss', "error", 'pre']
    metric_df_creator = MetricDataFrameCreator(param_weighted, "all")


    metrics = metric_df_creator.get_metrics_for_all()
    pd.DataFrame(metrics).to_csv("WeightedParams_All_feautures_metrics.csv")


if __name__ == "__main__":
    main()