from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
import xgboost as xgb
from numpyencoder import NumpyEncoder
import json
import pandas as pd
import numpy as np
class XGBModel:
    def __init__(self, data, target = "BasePair", split = 0.2, parms = None):
        if parms is None:
          self.parms = {'max_depth': 2, 'eta': 1, "scale_pos_weight": 100}
          self.parms['nthread'] = 4
          self.parms['objective'] = 'multi:softmax' if target == "Class" else 'binary:logistic'
          self.parms['num_class'] = 3 if target == "Class" else 2
        else:
          self.parms = parms

        self.bst = None
        if not target in ["BasePair", "BaseStack", "Class"]:
            raise ValueError("target must be 'BasePair' or 'BaseStack' or 'Class'")
        self.target = target
        if isinstance(data, dict):
            self.data = pd.DataFrame(data)
            self.data[self.target].fillna(False, inplace=True)
        else:
            self.data = data


        self.features = [feature for feature in self.data.columns if feature not in ["Unnamed: 0", "pdb_id", "nt1", "nt2", "BasePair", "BaseStack", "Class"]]
        
        
        if split:
            self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(self.data[self.features], self.data[target], test_size=split, random_state=42)
            self.dval = xgb.DMatrix(self.X_test, label=self.y_test)
            if target == "Class":
                weights = XGBModel.get_weight_from_class_distribution(self.y_train)
                self.dtrain = xgb.DMatrix(self.X_train, label=self.y_train, weight=weights)
            else: 
              self.dtrain = xgb.DMatrix(self.X_train, label=self.y_train)
            self.evallist = [(self.dtrain, 'train'), (self.dval, 'eval')]
        else:
            self.X_train = self.data[self.features]
            self.y_train = self.data[target]
            self.dval = None
            if target == "Class":
                weights = XGBModel.get_weight_from_class_distribution(self.y_train)
                self.dtrain = xgb.DMatrix(self.X_train, label=self.y_train, weight=weights)
            else: 
              self.dtrain = xgb.DMatrix(self.X_train, label=self.y_train, weight=weights)
            self.evallist = [(self.dtrain, 'train')]  

    def set_test_data(self, X_test, y_test):
        self.X_test = X_test
        self.y_test = y_test
        self.dval = xgb.DMatrix(self.X_test, label=self.y_test)

    def fit(self, num_round):
        bst = xgb.train(self.parms, self.dtrain, num_round,  evals = self.evallist, verbose_eval=False)
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
        if self.target == "Class":
          print(confusion_matrix(self.y_test, y_preds, labels = [0, 1, 2]))
        else:
          tn, fp, fn, tp = confusion_matrix(self.y_test, y_preds).ravel()
          print(f" True Negative: {tn}, False Positive: {fp}, False Negative: {fn}, True Positive: {tp}")
          return tn, fp, fn, tp
    
    def eval(self, test_data):
        dtest = xgb.DMatrix(test_data.drop(["BasePair", "BaseStack"], axis=1), label=test_data[self.target])
        return self.bst.eval(dtest)
    def predict(self, dtest):
        dtest = xgb.DMatrix(dtest[self.features], label=dtest[self.target])
        return self.bst.predict(dtest)
    
    @staticmethod
    def get_weight_from_class_distribution(y_train):
        proportions = np.bincount(y_train) / y_train.shape[0]
        smallest_propotion = np.min(proportions)

        class_weights =  smallest_propotion / proportions
        return [class_weights[y] for y in y_train]