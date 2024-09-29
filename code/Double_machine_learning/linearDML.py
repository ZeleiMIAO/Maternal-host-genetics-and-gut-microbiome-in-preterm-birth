"""
linearDML
train tsbc, test webirth
"""


# Some imports to get us started
# Utilities
import os
import urllib.request
import numpy as np
import pandas as pd
import shap
from sklearn.model_selection import train_test_split

# Generic ML imports
from sklearn.preprocessing import PolynomialFeatures
from sklearn.ensemble import GradientBoostingRegressor

# EconML imports
from econml.dml import LinearDML, CausalForestDML, SparseLinearDML
from econml.cate_interpreter import SingleTreeCateInterpreter, SingleTreePolicyInterpreter
from econml.sklearn_extensions.linear_model import WeightedLassoCVWrapper, WeightedLasso, WeightedLassoCV

from itertools import product
from sklearn.linear_model import Lasso, LassoCV, LogisticRegression, LogisticRegressionCV
import matplotlib.pyplot as plt

from lightgbm import LGBMClassifier, LGBMRegressor
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from scipy.stats import ttest_ind
from sklearn.model_selection import GridSearchCV

import matplotlib.pyplot as plt

%matplotlib inline
##############################################################################################################
# PRS and MRS

##############################################################################################################

# Load the data
df = pd.read_csv('thsbc_data.csv')

# Defining Y, T, X, W, etc. based on bootstrap_df
Y_train = df['preterm']  # Outcome of interest
T_train = df['prs']  # Intervention
X_train = df[['mrs']]  # Features
W_train = df[['parity', 'sexbaby', 'age']]  # Covariates
ids_train = df['id']  # Identifiers


df1 = pd.read_csv('webirth_data.csv')

# Defining Y, T, X, W, etc. based on bootstrap_df
Y_test = df1['preterm']  # Outcome of interest
T_test = df1['prs']  # Intervention
X_test = df1[['mrs']]  # Features
W_test = df1[['parity', 'sexbaby', 'age']]  # Covariates
ids_test = df1['spid']  # Identifiers

# Modeling and estimating effects
est = LinearDML(random_state=20,model_y=LogisticRegressionCV())
est.fit(Y_train, T_train, X=X_train, W=W_train, inference='bootstrap')

#### summary the results   
g_linear = est.summary()
ate = est.ate(X_train)
ate = est.ate_inference(X_train)

df_summary = pd.DataFrame(g_linear.tables[0])
df_summary.to_csv('res_linearDML_p.csv', index=False)

