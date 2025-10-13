This is the accompanying code of the paper "Bayesian Causal Effect Estimation using Staged Trees". The most relevant files are:

 - 'ppmx_Hamming.R' which implements the Metropolis-Hastings algorithm to learn staged trees using a Hamming distance prior;
 - 'ppx_Hamming_Z.R' which extends the previous algorithm to include information from continuous covariates;
 - 'bayesian_cluster_summary.R' which implements Bayesian clustering techniques to summarize a sample of posterior staged trees into a point estimate;
 - 'effect_sevt.R' which implements estimators for the Average Treatment Effect (ATE) based on as staged tree.

The other files either implement the simulation studies of the paper or the data applications.
