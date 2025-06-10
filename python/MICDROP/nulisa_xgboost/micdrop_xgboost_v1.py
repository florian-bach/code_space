from xgboost import XGBClassifier, plot_tree
from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import confusion_matrix, accuracy_score, roc_curve, auc
from sklearn.model_selection import GridSearchCV
import pandas as pd
import numpy as np
import time
import matplotlib.pyplot as plt
# read data
# 52 week z score slightly worse than conc for malaria & treatment arm
# # fold change from 8 to 52 AUC ~0.7 for both anymalar and treatmentarm
wide_data = pd.read_csv('~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/wide_data.csv')
X_train, X_test, y_train, y_test = train_test_split(wide_data.iloc[:,3:len(wide_data.columns)], wide_data['anymalar'], test_size=.2)
le = LabelEncoder()
y_train = le.fit_transform(y_train)

param_grid = {
    'n_estimators': [50, 100, 200, 400],
    'learning_rate': [0.01, 0.05, 0.1, 0.3],
    'max_depth': [2, 3, 5, 9]
}

model = XGBClassifier(eval_metric='logloss')

grid = GridSearchCV(
    model, 
    param_grid,
    cv=5,
    scoring='roc_auc',
    verbose=1,
    n_jobs=-1
)
start_time = time.time()
grid.fit(X_train, y_train)
print('Fit time : ', time.time() - start_time)

print(grid.best_params_)



# for treatmentarm: n_estimators=100, max_depth=3, learning_rate=0.3
# for anymalar:  'n_estimators': 500, 'max_depth': 9, learning_rate=0.05
bst = XGBClassifier(n_estimators=100, max_depth=3, learning_rate=0.3, objective='binary:logistic')
# fit model
bst.fit(X_train, y_train)
# make predictions
preds = bst.predict(X_test)

preds = le.inverse_transform(preds)
cm = confusion_matrix(y_test, preds)
print(cm)
accuracy_score(y_test, preds)




plot_tree(bst, num_trees=0)
plt.show()
y_pred_proba = bst.predict_proba(X_test)[:, 1]
fpr, tpr, thresholds = roc_curve(le.fit_transform(y_test), y_pred_proba)

# Calculate the area under the ROC curve (AUC)
roc_auc = auc(fpr, tpr)

# Plot the ROC curve
plt.figure(figsize=(8, 6))
plt.plot(fpr, tpr, color='blue', label=f'ROC curve (AUC = {roc_auc:.2f})')
plt.plot([0, 1], [0, 1], color='red', linestyle='--', label='Random guess')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic (ROC) Curve')
plt.legend(loc="lower right")
plt.show()
