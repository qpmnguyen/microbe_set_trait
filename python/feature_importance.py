import numpy as np 
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.calibration import CalibratedClassifierCV
from sklearn.model_selection import KFold, cross_val_score, cross_validate, train_test_split
from sklearn.preprocessing import FunctionTransformer
from sklearn.metrics import accuracy_score, brier_score_loss, roc_auc_score
from sklearn.inspection import permutation_importance
from skbio.stats.composition import clr, multiplicative_replacement
from sklearn.pipeline import Pipeline
from sklearn import preprocessing
from pred_eval import prior_preprocess, clr_transform, create_pipeline 
import timeit

if __name__ == "__main__":
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
    
    print("Using threads:" + str(snakemake.threads))
    
    print("Test set feature importance using permutational test")
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=160497)
    fitted_mod = pipe.fit(X_train, y_train)
    
    y_pred = fitted_mod.predict(X_test)
    perf_value = roc_auc_score(y_test, y_pred)
    
    start = timeit.default_timer()
    r = permutation_importance(fitted_mod, X_test, y_test, n_repeats = 20, 
                                n_jobs = snakemake.threads, 
                               scoring = "roc_auc", random_state=160497)
    end = timeit.default_timer()
    
    print("Time: ", end - start)
    # top 10 features 
    sorted_idx = r.importances_mean.argsort()[range(0,11)]
    print(sorted_idx)
    importances = r.importances[sorted_idx].T
    labels = feat.columns[sorted_idx]
    out = pd.DataFrame(importances, columns =labels)
    out["performance"] = np.repeat(perf_value, out.shape[0])
    
    out.to_csv(snakemake.output[0])

