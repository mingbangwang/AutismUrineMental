rm(list=ls())

library(optparse)

option_list = list(
        make_option("--feature_profile", type="character", default=NULL, help="Input feature_profile file, the first two columns are Sample_names, Group_names. default: %default [required]" ),
        make_option("--importance_profile", type="character", default=NULL, help="Input importance_profile file, the first column is feature name., second column is feature coef/IncNodePurity default: %default [required]" ),
        make_option("--importance_type", type="character", default=NULL, help="Input importance_type,lasso or randomforest: %default [required]" ),
        make_option("--out",type="character", default=NULL, help="output dir, default: %default

  Example:
  Rscript ROC.R --feature_profile=feature_profile --importance_profile=importance_profile --importance_type --out=out_dir")

)

opt = parse_args(OptionParser(usage="%prog[option]", option_list=option_list, add_help_option=TRUE, prog="ROC.R"), print_help_and_exit=TRUE, positional_arguments=FALSE)


feature_profile = opt$feature_profile
importance_profile = opt$importance_profile
importance_type = opt$importance_type
workdir = opt$out


options(digits = 4)
q_975 <- qnorm(0.975)


## 
library(pROC)
library(rms)
library(randomForest)
library(glue)

setwd(workdir) 
dir.create("ROC") 
outdir <- 'ROC/' 

# 
profile <- read.csv(feature_profile,sep ="\t",header=TRUE) 

sample_col <- names(profile)[1] # "Sample_names" # 
group_col <- names(profile)[2] # "Group_names"

# 
profile_samplecol <- names(profile) %in% c(sample_col)
profile <-profile[!profile_samplecol]

dd=datadist(profile)
options(datadist="dd") 


## formula
imp <- read.csv(importance_profile,sep ="\t",header=TRUE) # 

featureIDs <-imp[ ,names(imp)[1]]
featureCount <- length(featureIDs)
# featureCount

tops <- seq(3,featureCount,3)
topLabel1 <- glue("Top{tops[1]}_importance_features")
topLabel2 <- glue("Top{tops[2]}_importance_features")
topLabel3 <- glue("all_{featureCount}_importance_features")

formula_1 <- as.formula(paste(group_col,"~",paste(featureIDs[1:tops[1]],collapse=" + "))) # top 5
formula_2 <- as.formula(paste(group_col,"~",paste(featureIDs[1:tops[2]],collapse=" + "))) # top 10
formula_3 <- as.formula(paste(group_col,"~",paste(featureIDs[1:featureCount],collapse=" + "))) # all


if (identical(toupper(importance_type),"LASSO")){
  print("LASSO")
  # LASSO 
  fit1 <- glm(formula=formula_1,data = profile,family = binomial())
  fit2 <- glm(formula=formula_2,data = profile,family = binomial())
  fit3 <- glm(formula=formula_3,data = profile,family = binomial())
}  else {
  print("randomForest")
  fit1 <- randomForest(formula_1,data=profile) # Mortality, Group_names
  fit2 <- randomForest(formula_2,data=profile) # Mortality, Group_names
  fit3 <- randomForest(formula_3,data=profile) # Mortality, Group_names
}


profile$predvalue1 <- predict(fit1)
profile$predvalue2 <- predict(fit2)
profile$predvalue3 <- predict(fit3)

# AUC
ROC1 <- roc(profile[,group_col],profile$predvalue1)
ROC2 <- roc(profile[,group_col],profile$predvalue2)
ROC3 <- roc(profile[,group_col],profile$predvalue3)


# plot
pdf(file=glue("{outdir}/ROC_{importance_type}.pdf"),width = 5,height = 4) # width = 8,height = 8

plot(1-ROC1$specificities,ROC1$sensitivities,type = "l",col="red",lty=1,xlab = "1-specificities",ylab = "sensitivities",lwd=2)
lines(1-ROC2$specificities,ROC2$sensitivities,col="green",lty=1,lwd=2)
lines(1-ROC3$specificities,ROC3$sensitivities,col="blue",lty=1,lwd=2)
abline(0,1)
legend(0.01,0.55,
       c(paste("top",tops[1], " features auc= ",round(auc(ROC1),3),": (",round(ci(auc(ROC1))[1],3),"-",round(ci(auc(ROC1))[3],3),")",sep = ""),
         paste("top",tops[2], " features auc= ",round(auc(ROC2),3),": (",round(ci(auc(ROC2))[1],3),"-",round(ci(auc(ROC2))[3],3),")",sep = ""),
         paste("all",featureCount, " features auc= ",round(auc(ROC3),3),": (",round(ci(auc(ROC3))[1],3),"-",round(ci(auc(ROC3))[3],3),")",sep = "")),
       lty = c(1,1,1),
       lwd = c(2,2,2),
       col = c("red","green","blue"),
       bty = "n" 
)
title(glue("{importance_type}_importance_features"))

dev.off()

