# Regionwise classification of Self and nonself. The code produces ROC curves in figures
import numpy as np
from sklearn.datasets import make_classification
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from imblearn.pipeline import make_pipeline as make_pipeline_imb
from imblearn.under_sampling import RandomUnderSampler
from imblearn.under_sampling import ClusterCentroids
import matplotlib.pyplot as plt
import pandas as pd

# Generate synthetic data (replace this with your actual data loading code)
# X, y = make_classification(n_samples=1000, n_features=20, n_classes=2, weights=[0.5, 0.5], random_state=42)

df = pd.read_csv('./DATA/REG_BP_DT.csv')
X = np.asarray(df.iloc[:, 0:2])
y = np.asarray(df.iloc[:, 3])
y = y - 1

# Define the outer cross-validation procedure
cv_outer = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

# Collect the ROC curve data
tprs = []
aucs = []
mean_fpr = np.linspace(0, 1, 100)
fig, ax = plt.subplots()

for train_ix, test_ix in cv_outer.split(X, y):
    # Split data
    X_train, X_test = X[train_ix], X[test_ix]
    y_train, y_test = y[train_ix], y[test_ix]
    
    # Define the pipeline
    pipeline = make_pipeline_imb(
        StandardScaler(),
        MinMaxScaler(),
        ClusterCentroids(random_state=42),
        RandomUnderSampler(random_state=42),
        DecisionTreeClassifier(random_state=42)
    )
    
    # Fit the model
    pipeline.fit(X_train, y_train)
    y_proba = pipeline.predict_proba(X_test)[:, 1]
    
    # Compute ROC curve and AUC for this fold
    fpr, tpr, _ = roc_curve(y_test, y_proba)
    roc_auc = auc(fpr, tpr)
    aucs.append(roc_auc)
    
    # Interpolate the ROC curve
    tprs.append(np.interp(mean_fpr, fpr, tpr))
    tprs[-1][0] = 0.0
    
    # Plot the ROC curve for the current fold
    ax.plot(fpr, tpr, lw=1, alpha=0.3, label=f'Fold ROC (AUC = {roc_auc:.2f})')

# Plot the mean ROC curve
mean_tpr = np.mean(tprs, axis=0)
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
std_auc = np.std(aucs)
ax.plot(mean_fpr, mean_tpr, color='blue', label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc), lw=2, alpha=0.8)

# Plot the chance line
ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='red', label='Chance', alpha=0.8)

# Finalize the plot
ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05], title="Receiver Operating Characteristic - Region-Wise")
ax.legend(loc="lower right")
plt.show()
