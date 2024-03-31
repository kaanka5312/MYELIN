# This script produces the figure 1C. It basically performs # # 10x10 nested CV and produces ROC curves.

import numpy as np
from sklearn.datasets import load_breast_cancer
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import cross_val_score, KFold, train_test_split
from sklearn.metrics import roc_curve, auc
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, MinMaxScaler
import matplotlib.pyplot as plt
import pandas as pd

# Load data
# data = load_breast_cancer()
# X, y = data.data, data.target

df = pd.read_excel('E:/EIB/DATA/DAT.xlsx')
X = np.asarray(df.iloc[:, [0,1,2,4] ])
y = np.asarray(df.iloc[:, 5])
y = y - 1

# Initialize variables for plotting
tprs = []
aucs = []
mean_fpr = np.linspace(0, 1, 100)

fig, ax = plt.subplots()

# Define a pipeline with preprocessing steps and a classifier
pipeline = Pipeline([
    ('scaler', StandardScaler()),
    ('minmax', MinMaxScaler()),
    ('classifier', DecisionTreeClassifier())
])

'''pipeline = Pipeline([
    ('classifier', DecisionTreeClassifier())
]) '''


# Outer CV
outer_cv = KFold(n_splits=10, shuffle=True, random_state=42)

for i, (train_index, test_index) in enumerate(outer_cv.split(X, y)):
    X_train, X_test = X[train_index], X[test_index]
    y_train, y_test = y[train_index], y[test_index]

    # Fit the pipeline on the training data
    pipeline.set_params(classifier__random_state=i)  # Ensure reproducibility
    pipeline.fit(X_train, y_train)
    y_score = pipeline.predict_proba(X_test)[:, 1]

    # Compute ROC curve and AUC for this fold
    fpr, tpr, thresholds = roc_curve(y_test, y_score)
    tprs.append(np.interp(mean_fpr, fpr, tpr))
    tprs[-1][0] = 0.0
    roc_auc = auc(fpr, tpr)
    aucs.append(roc_auc)
    ax.plot(fpr, tpr, lw=1, alpha=0.3, label=f'ROC fold {i+1} (AUC = {roc_auc:.2f})')

# Plot the mean ROC curve
mean_tpr = np.mean(tprs, axis=0)
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
std_auc = np.std(aucs)
ax.plot(mean_fpr, mean_tpr, color='blue',
        label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
        lw=2, alpha=0.8)

# Plot chance level
ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='red', label='Chance', alpha=0.8)

# Formatting the plot
ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05],
       title="Receiver Operating Characteristic - Subject-wise")
ax.legend(loc="lower right")
plt.show()
