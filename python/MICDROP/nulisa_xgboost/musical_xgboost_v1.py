from xgboost import XGBClassifier, plot_tree
from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import confusion_matrix, accuracy_score, roc_curve, auc
from sklearn.model_selection import GridSearchCV
from random import sample
import pandas as pd
import numpy as np
import time
import random
import matplotlib.pyplot as plt
# read data
wide_data = pd.read_csv('~/postdoc/stanford/plasma_analytes/MUSICAL/combo/wide_data.csv')

# 15 singletons 31 paired data; 36A & 41S;
# pick 8 individuals with both A&S data as test set

#get random list of 8 individuals with both timepoints for test set
#Count number of occurrences of each id

#id_counts = wide_data['id'].value_counts()
#Filter IDs that appear more than once

#duplicate_ids = id_counts[id_counts > 1].index
test_ids= random.sample(list(duplicate_ids), k=8)

#test_ids=[268, 662, 355, 225, 132, 463, 385, 218]


#tilde negates..!
X_train=wide_data[~wide_data['id'].isin(test_ids)].iloc[:, 2:len(wide_data.columns)]
X_test=wide_data[wide_data['id'].isin(test_ids)].iloc[:, 2:len(wide_data.columns)]
y_train=wide_data[~wide_data['id'].isin(test_ids)].iloc[:, 1]
y_test = wide_data[wide_data['id'].isin(test_ids)].iloc[:, 1]

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
#start_time = time.time()
grid.fit(X_train, y_train)
#print('Fit time : ', time.time() - start_time)

print(grid.best_params_)



bst = XGBClassifier(n_estimators=50, max_depth=3, learning_rate=0.05, objective='binary:logistic')
# fit model
bst.fit(X_train, y_train)
# make predictions
preds = bst.predict(X_test)

preds = le.inverse_transform(preds)
cm = confusion_matrix(y_test, preds)
print(cm)
accuracy_score(y_test, preds)


y_pred_proba = bst.predict_proba(X_test)[:, 1]
fpr, tpr, thresholds = roc_curve(le.fit_transform(y_test), y_pred_proba)

# Calculate the area under the ROC curve (AUC)
roc_auc = auc(fpr, tpr)

plot_tree(bst, num_trees=0)
plt.show()
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




