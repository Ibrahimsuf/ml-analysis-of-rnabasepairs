from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
import xgboost as xgb
from numpyencoder import NumpyEncoder
import json
import pandas as pd
class XGBModel:
    def __init__(self, data, target = "BasePair", split = 0.2, parms = None):
        if parms is None:
          self.parms = {'max_depth': 2, 'eta': 1, 'objective': 'binary:logistic', "scale_pos_weight": 100}
          self.parms['nthread'] = 4
          self.parms['eval_metric'] = ['logloss', "error", 'pre']
        else:
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


        self.features = [feature for feature in self.data.columns if feature not in ["Unnamed: 0", "pdb_id", "nt1", "nt2", "BasePair", "BaseStack"]]
        
        
        if split:
            self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(self.data[self.features], self.data[target], test_size=split, random_state=42)
            self.dval = xgb.DMatrix(self.X_test, label=self.y_test)
            self.dtrain = xgb.DMatrix(self.X_train, label=self.y_train)
            self.evallist = [(self.dtrain, 'train'), (self.dval, 'eval')]
        else:
            self.X_train = self.data[self.features]
            self.y_train = self.data[target]
            self.dval = None
            self.dtrain = xgb.DMatrix(self.X_train, label=self.y_train)
            self.evallist = [(self.dtrain, 'train')]  

    def set_test_data(self, X_test, y_test):
        self.X_test = X_test
        self.y_test = y_test
        self.dval = xgb.DMatrix(self.X_test, label=self.y_test)

    def fit(self, num_round):
        bst = xgb.train(self.parms, self.dtrain, num_round, evals = self.evallist, verbose_eval=False)
        self.bst = bst

    
    def fit_save_and_evaluate(self, num_round, write_to_file=False, name=None):
        self.fit(num_round)
        if name:
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
        dtest = xgb.DMatrix(dtest[self.features], label=dtest[self.target])
        return self.bst.predict(dtest)