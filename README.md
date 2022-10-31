# run shell scripts



## lasso

```shell
feature_file_num=/d/python_learning/project/ASD_metal/manuscript/Scripts/all_diff_features_num.tsv
outdir=/d/python_learning/project/ASD_metal/manuscript/Scripts/out
mkdir -p $outdir
Rscript lasso.R --feature_profile=$feature_file_num --out=$outdir
```



## ROC

```shell
importance_profile=$outdir"/lasso/lasso_coefs_top.txt"
importance_type=LASSO
Rscript ROC.R --feature_profile=$feature_file_num --importance_profile=$importance_profile --importance_type=$importance_type --out=$outdir
```



## 

