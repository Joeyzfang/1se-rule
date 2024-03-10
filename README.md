# 1se-rule
## This R program is used to estimate the performance of 1 se rule in cross validation when tuning logistic model.
The estimation would contain 3 main parts: i) standard error of cross validation error;
ii) accuracy of variable selection; iii) out-of-sample prediction error. 
The fourth part is meta-model building.

In standard error of CV error section, simulation involves both lasso and ridge. We only conducted estimation analysis on lasso's performance.
Due to the massive amount of computing time and too many raw data files, we provided data that could be visualized directly.

