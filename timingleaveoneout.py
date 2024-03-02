import pandas as pd
from sklearn.model_selection import LeaveOneOut
import xgboost as xgb
import time
import psutil
import os
def getclass(row):
  if row["BasePair"] == True:
    return 1
  elif row["BaseStack"] == True:
    return 2
  else:
    return 0
  
CG = pd.read_csv("CG_0_3_data.tar.gz")

CG["Class"] = CG.apply(getclass, axis=1, threads=4)
params = {'max_depth': 2, 'eta': 1, "scale_pos_weight": 100}
params['nthread'] = 4
params["device"] = "cuda"
params['objective'] = 'multi:softmax' 
params['num_class'] = 3

cv = LeaveOneOut()
y_true = []
y_pred = []
features = [feature for feature in CG.columns if feature not in ["Unnamed: 0", "pdb_id", "nt1", "nt2", "BasePair", "BaseStack", "Class"]]
for train_index, test_index in [next(cv.split(CG))]:
  start_time = time.time()
  X_train, X_test = CG.iloc[train_index], CG.iloc[test_index]
  y_train, y_test = X_train["Class"], X_test["Class"]

  dtrain = xgb.DMatrix(X_train[features], label=y_train)
  bst = xgb.train(params, dtrain, 100)

  y_pred.append(bst.predict(xgb.DMatrix(X_test[features])))
  y_true.append(y_test)
  print("--- %s seconds ---" % (time.time() - start_time))
  # Get the current process ID
  pid = os.getpid()
  # Get the number of threads used by the process
  threads = psutil.Process(pid).num_threads()
  break

print("Number of threads:", threads)


