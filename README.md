# R_Parafac2_CV_Example_Code

A set of R code scripts that I wrote for an analytic research project currently in development. 

This code implements a novel dimensionality reduction method (<ins>**Par**</ins>allel <ins>**Fac**</ins>tor Analysis - 2, or Parafac2) applied to a four-dimensional multi-nested array of human EEG data (~2 million entries). Parafac2 can explain variation across all four dimensions (modes) of the data *simultaneously* in a single model. Unlike bimodal methods such as PCA and ICA, this method requires much fewer free parameters to model the same data, allows us to use all data without unfolding the data and thus violating the multiway structure of the data, and provides a unique solution (e.g., no rotational indeterminacy like in PCA/ICA) that aids in interpreting the factors.

The optimal model is chosen using cross-validation of the change in mean square error, and is then validated via a split-half analysis. The resulting model parsimoniously explains over 60% of the variance across all four dimensions *simultaneously* with only five factors.

The code is organized into several scripts with different purposes:

1) Structure and preprocess the raw data into a four-way data array,
2) Visualize the raw data, 
3) Conduct Parafac2 modeling of he four-way array,
4) Use the cross-validated mean square error to decide the optimal number of factors to retain,
5) Conduct split-half analyses and calculate Tucker congruence coefficients to validate the chosen model, 
6) Visualize the final Parafac2-derived model, and
7) Conduct statistical analyses using linear mixed models to understand the real-world significance of the Parafac2 model factors.
