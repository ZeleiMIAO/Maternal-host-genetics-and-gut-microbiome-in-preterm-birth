"""
forestDML
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

##############################################################################################################
# CLR transformed data
##############################################################################################################

df = pd.read_csv('data.csv')

# Ensure there are no missing values before proceeding
assert not df.isna().any().any(), "There are NaN values in the DataFrame"

# Extract microbial feature columns
cols = list(df.columns)
mic_cols = cols[18:221]  

# Extract original X and IDs for inference and saving results
original_X = df[mic_cols]

original_ids = df['id']

# Parameters
n_bootstrap_samples = 1000

# Storing individual bootstrap results
shap_vals_bootstrap = []

# Main Loop for Bootstrap Sampling
for i in range(n_bootstrap_samples):
    # Perform bootstrap sampling
    bootstrap_indices = np.random.choice(df.index, size=len(df), replace=True)
    bootstrap_df = df.loc[bootstrap_indices].reset_index(drop=True)
    
    # Splitting into train/test for the bootstrap sample
    train_indices, test_indices = train_test_split(np.arange(len(bootstrap_df)), test_size=0.2, random_state=i)
    
    # Defining Y, T, X, W, etc. based on bootstrap_df
    Y = bootstrap_df['preterm']  # Outcome of interest
    T = bootstrap_df['prs']  # Intervention
    X = bootstrap_df[mic_cols]  # Features
    W = bootstrap_df[['parity', 'sexbaby', 'age']]  # Covariates
    ids = bootstrap_df['id']  # Identifiers
    
    X_train, X_test = X.iloc[train_indices], X.iloc[test_indices]
    Y_train, Y_test = Y.iloc[train_indices], Y.iloc[test_indices]
    W_train, W_test = W.iloc[train_indices], W.iloc[test_indices]
    T_train, T_test = T.iloc[train_indices], T.iloc[test_indices]
    id_train, id_test = ids.iloc[train_indices], ids.iloc[test_indices]

    # Modeling and estimating effects
    est = CausalForestDML(model_y=LGBMRegressor(n_estimators=300),
                          model_t=LGBMRegressor(n_estimators=300),
                          n_estimators=600,
                          min_samples_leaf=5,
                          max_depth=100,
                          max_samples=0.2,
                          discrete_treatment=False,
                          random_state=i)
    
    
    est.fit(Y_train, T_train, X=X_train, W=W_train)
   
    # shap values for original data (gut microbiota)
    shap_values = est.shap_values(original_X, feature_names=original_X.columns)['preterm']['prs']
    shap_vals_bootstrap.append(shap_values.values)

    print(f"Iteration {i+1}/{n_bootstrap_samples} Completed")

# Averaging shap values over all bootstrap samples
shap_vals_avg = np.mean(shap_vals_bootstrap, axis=0)

# Saving averaged results
shap_avg_df = pd.DataFrame(shap_vals_avg, columns=X.columns)
shap_avg_df.to_csv('shap.csv', index=False)

