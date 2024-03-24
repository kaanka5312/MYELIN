
# -------------------- GENERATED WITH PHOTON WIZARD (beta) ------------------------------
# PHOTON Project Folder: E:/EIB/photonai

import pandas as pd
import numpy as np
from photonai.base import Hyperpipe, PipelineElement, OutputSettings, Switch, Preprocessing
from photonai.optimization import Categorical, IntegerRange, FloatRange
from sklearn.model_selection import StratifiedKFold
             
# Specify how results are going to be saved
# Define hyperpipe
hyperpipe = Hyperpipe('selfnonselfclassification',
                      project_folder = 'E:/EIB/photonai',
                      optimizer="random_grid_search",
                      optimizer_params={'n_configurations': 30},
                      metrics=['accuracy', 'balanced_accuracy', 'f1_score', 'auc'],
                      best_config_metric="auc",
                      outer_cv = StratifiedKFold(n_splits=10,shuffle=True),
                      inner_cv = StratifiedKFold(n_splits=10, shuffle=True))
        
# Add transformer elements
hyperpipe += PipelineElement("StandardScaler", hyperparameters={}, 
                             test_disabled=False, with_mean=True, with_std=True)

hyperpipe += PipelineElement("MinMaxScaler", hyperparameters={}, 
                             test_disabled=False)

# Add estimator
estimator_switch = Switch('EstimatorSwitch')
#estimator_switch += PipelineElement("SVC", hyperparameters={'C': FloatRange(0.5, 2), 'kernel': ['linear', 'rbf', 'poly']}, gamma='scale', max_iter=1000000)
estimator_switch += PipelineElement("RandomForestClassifier", hyperparameters={'n_estimators': IntegerRange(5, 20), 'min_samples_split': IntegerRange(2,5), 'min_samples_leaf': IntegerRange(1,3)}, criterion='gini', max_depth=None)
estimator_switch += PipelineElement("DecisionTreeClassifier", hyperparameters={'min_samples_leaf': IntegerRange(1,3), 'min_samples_split': IntegerRange(2,5)}, criterion='gini', max_depth=None)
hyperpipe += estimator_switch                

# Load data
df = pd.read_excel('E:/EIB/DATA/DAT.xlsx')
X = np.asarray(df.iloc[:, [0,1,2,4] ])
y = np.asarray(df.iloc[:, 5])

# Fit hyperpipe
hyperpipe.fit(X, y)

handler=hyperpipe.results_handler

performance_table = handler.get_performance_table()
performance_table.to_csv('E:/EIB/DATA/performance_table.csv', index=False)

with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    print(performance_table)
print(" ")
