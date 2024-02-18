from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
import xgboost as xgb
from numpyencoder import NumpyEncoder
import json
import pandas as pd
class XGBModel:
    def __init__(self, parms, data, target = "BasePair"):
        self.parms = parms
        self.bst = None
        if not target in ["BasePair", "BaseStack"]:
            raise ValueError("target must be 'BasePair' or 'BaseStack'")
        self.target = target
        if isinstance(data, dict):
            self.data = pd.DataFrame(data)
            self.data[self.target].fillna(False, inplace=True)
        else:
            self.data = data


        features = [feature for feature in self.data.columns if feature not in ["Unnamed: 0", "pdb_id", "nt1", "nt2", "BasePair", "BaseStack"]]
        self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(self.data[features], self.data[target], test_size=0.2, random_state=42)
        self.dtrain = xgb.DMatrix(self.X_train, label=self.y_train)
        self.dval = xgb.DMatrix(self.X_test, label=self.y_test)
        self.dtrain = self.dtrain
        self.evallist = [(self.dtrain, 'train'), (self.dval, 'eval')]

    def fit(self, num_round):
        bst = xgb.train(self.parms, self.dtrain, num_round, evals = self.evallist, verbose_eval=False)
        self.bst = bst

    
    def fit_save_and_evaluate(self, name, num_round, write_to_file=False):
        self.fit(num_round)
        self.save_model(name)
        if write_to_file:
            classification_report = self.get_classifcation_report(output_dict=True)
            tn, fp, fn, tp = self.get_confusion_matrix()
            classification_report["tn"] = tn
            classification_report["fp"] = fp
            classification_report["fn"] = fn
            classification_report["tp"] = tp
            with open(f"{name}_classifcation_report", "w") as f:
                json.dump(classification_report, f, cls=NumpyEncoder)
        else: 
            self.get_classifcation_report()
            self.get_confusion_matrix()
    
    def save_model(self, name):
        self.bst.save_model(name)
    
    def get_classifcation_report(self, output_dict=False):
        y_preds = self.bst.predict(self.dval).round()
        if output_dict:
            return classification_report(self.y_test, y_preds, output_dict=True)
        else:
            print(classification_report(self.y_test, y_preds, output_dict=False))
    def get_confusion_matrix(self):
        y_preds = self.bst.predict(self.dval).round()
        tn, fp, fn, tp = confusion_matrix(self.y_test, y_preds).ravel()
        print(f" True Negative: {tn}, False Positive: {fp}, False Negative: {fn}, True Positive: {tp}")
        return tn, fp, fn, tp
    
    def eval(self, test_data):
        dtest = xgb.DMatrix(test_data.drop(["BasePair", "BaseStack"], axis=1), label=test_data[self.target])
        return self.bst.eval(dtest)
    def predict(self, dtest):
        return self.bst.predict(dtest)