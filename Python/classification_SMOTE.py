
# -------------------- GENERATED WITH PHOTON WIZARD (beta) ------------------------------
# PHOTON Project Folder: E:/EIB/

import pandas as pd
import numpy as np
from photonai.base import Hyperpipe, PipelineElement, OutputSettings, Switch, Preprocessing
from photonai.optimization import Categorical, IntegerRange, FloatRange
from sklearn.model_selection import StratifiedKFold
             
# Specify how results are going to be saved
# Define hyperpipe
hyperpipe = Hyperpipe('classificationsmote',
                      project_folder = 'E:/EIB/',
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
hyperpipe += PipelineElement("ImbalancedDataTransformer", hyperparameters={'method_name': ['RandomUnderSampler', 'ClusterCentroids']}, 
                             test_disabled=False)

# Add estimator
estimator_switch = Switch('EstimatorSwitch')
estimator_switch += PipelineElement("SVC", hyperparameters={'C': FloatRange(0.5, 2), 'kernel': ['linear', 'rbf']}, gamma='scale', max_iter=1000000)
estimator_switch += PipelineElement("RandomForestClassifier", hyperparameters={'n_estimators': IntegerRange(5, 20), 'min_samples_split': IntegerRange(2,5), 'min_samples_leaf': IntegerRange(1,3)}, criterion='gini', max_depth=None)
estimator_switch += PipelineElement("DecisionTreeClassifier", hyperparameters={'min_samples_leaf': IntegerRange(1,3), 'min_samples_split': IntegerRange(2,5)}, criterion='gini', max_depth=None)

hyperpipe += estimator_switch                

# Load data
df = pd.read_csv('E:/EIB/DATA/MED_2.csv')
X = np.asarray(df.iloc[:, 0:2])
y = np.asarray(df.iloc[:, 3])

# Fit hyperpipe
hyperpipe.fit(X, y)

handler=hyperpipe.results_handler

performance_table = handler.get_performance_table()
performance_table.to_csv('E:/EIB/DATA/performance_table_RegWise.csv', index=False)

