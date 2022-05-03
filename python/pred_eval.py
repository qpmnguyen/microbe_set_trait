import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.calibration import CalibratedClassifierCV
from sklearn.model_selection import KFold, cross_val_score, cross_validate
from sklearn.preprocessing import FunctionTransformer
from sklearn.metrics import accuracy_score, brier_score_loss, roc_auc_score
from sklearn.inspection import permutation_importance
from skbio.stats.composition import clr, multiplicative_replacement
from sklearn.pipeline import Pipeline
from sklearn import preprocessing
import pickle 

"""
This function performs some apriori preprocessing based on data set required 
"""
def prior_preprocess(feat, meta, outcome_label, pos_class):
    X = feat.to_numpy()
    Y = meta.loc[:, outcome_label].to_numpy()
    Y = Y[~np.all(X == 0, axis = 1)]
    y = [1 if i == pos_class else 0 for i in Y]
    X = X[~np.all(X == 0, axis=1)]
    return(X, y)

"""
This function performs the clr transformation with imputation 
"""
def clr_transform(arr, pseudocount=True):
    if pseudocount == False:
        # impute with multiplicative replacement
        arr = multiplicative_replacement(arr)
    else:
        arr = arr + 10e-5
    # clr transformation 
    arr = clr(arr)
    return(arr)

"""
This function creates a scikit-learn pipeline 
Here we use a base default random forest classifier and calibrated it using 5-fold cross-validation. 
The remainder calibrated classifier is an ensemble of the individual classifiers. If the sample size is low, 
we're going to use sigmoid calibration approach. 
For preprocessing, we created a custom preprocessing function using the CLR transformation using scikit-bio
"""
def create_pipeline(clr_trans=True):
    rf = RandomForestClassifier(n_estimators=500, max_features='sqrt')
    estimators = [("calib_rf", CalibratedClassifierCV())]
    if clr_trans == True:
        transf = FunctionTransformer(clr_transform, validate=True, kw_args={'pseudocount': True})
        estimators.insert(0, ("clr_transformer", transf))
    pipe = Pipeline(estimators)
    pipe.set_params(calib_rf__base_estimator=rf, calib_rf__cv=5, calib_rf__ensemble=True, 
                    calib_rf__method='sigmoid')
    return(pipe)

if __name__ == "__main__":
    np.random.seed(210595)
    print("Load and preprocess data")
    feat_path = snakemake.input[0]
    meta_path = snakemake.input[1]
    
    # parameters to inputs
    outcome_label = snakemake.params[0]
    pos_class = snakemake.params[1]
    
    # clarification for outputs
    condition = snakemake.params[2]
    sequencing = snakemake.params[3]
    input_type = snakemake.params[4]
    
    feat = pd.read_csv(feat_path, index_col=0)
    meta = pd.read_csv(meta_path, index_col=0)
    
    X, y = prior_preprocess(feat, meta, outcome_label, pos_class)
    print("Create pipeline")
    pipe = create_pipeline(clr_trans=snakemake.params[5])
    
    print("Fit pipeline via 10-fold cross validation")
    scoring = ["neg_brier_score", "roc_auc"]
    scores = cross_validate(pipe, X, y, scoring = scoring, cv=10)
    
    out = pd.DataFrame({
        "condition": np.repeat(np.array([condition]), 10),
        "sequencing": np.repeat(np.array([sequencing]), 10),
        "input_type": np.repeat(np.array([input_type]), 10),
        "fold": range(1,11),
        "brier": -1 * scores['test_neg_brier_score'],
        "roc_auc": scores["test_roc_auc"]
    })
    
    out.to_csv(snakemake.output[0])
    
    
    
    
