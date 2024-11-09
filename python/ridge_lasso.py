# load and summarize housing dataset

import pandas as pd
from matplotlib import pyplot
from numpy import mean
from numpy import std
from numpy import absolute
from numpy import arange
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedKFold
from sklearn.linear_model import RidgeCV
from sklearn.linear_model import Ridge
from sklearn.metrics import roc_curve



# load the mtcars dataset; header=0 means infer
dataframe = read_csv("~/Downloads/nulisa_for_ridge.csv", header=0)
# [:, :-1] is all rows, and all columns except the last one
# [:, -1]  is all rows and ONLY the last column
X = dataframe.iloc[:, 3:254].values
y = dataframe['class'].values

# define model
#model = Ridge(alpha=1.0)

# define model evaluation method; 10 fold crossvalidation, report average mean absolute error
cv = RepeatedKFold(n_splits=10, n_repeats=3, random_state=1)
model = RidgeCV(alphas=arange(0, 1, 0.01), cv=cv, scoring='neg_mean_absolute_error')
# fit model
model.fit(X, y)
# summarize chosen configuration
model.coef_

# making a binary classifier

from sklearn.linear_model import RidgeClassifier
clf = RidgeClassifier(alpha=1.0) # Adjust alpha for regularization strength
clf.fit(X, y)
y_pred = clf.predict(X)
confidence = clf.decision_function(X)
clf.roc_curve

