from pred_eval import prior_preprocess, clr_transform, create_pipeline 

def get_feat():
    


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

