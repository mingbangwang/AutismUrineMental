rm(list=ls())

library(optparse)

option_list = list(
        make_option("--feature_profile", type="character", default=NULL, help="Input feature_profile file, the first two columns are Sample_names, Group_names. default: %default [required]" ),
        make_option("--trait_profile", type="character", default=NULL, help="Input trait_profile file, the first two columns are Sample_names, Group_names. default: %default [required]" ),
        make_option("--importance_profile", type="character", default=NULL, help="Input importance_profile file, the first column is feature name., second column is feature coef/IncNodePurity default: %default [required]" ),
        make_option("--out",type="character", default=NULL, help="output dir, default: %default

  Example:
  Rscript calibration.R --feature_profile=feature_profile --trait_profile=trait_profile --importance_profile=importance_profile --out=out_dir")

)

opt = parse_args(OptionParser(usage="%prog[option]", option_list=option_list, add_help_option=TRUE, prog="calibration.R"), print_help_and_exit=TRUE, positional_arguments=FALSE)

options(digits = 3)
q_975 <- qnorm(0.975)

feature_profile = opt$feature_profile
trait_profile <- opt$trait_profile
importance_profile = opt$importance_profile
workdir = opt$out
merge_trait <- FALSE # merge trait


## 载入R包
library(riskRegression)
library(rms)
library(glue)

setwd(workdir) 
dir.create("calibration") 
outdir <- 'calibration/' 

profile <- read.csv(feature_profile,sep ="\t",header=TRUE) 

sample_col <- names(profile)[1] # "Sample_names" # 
group_col <- names(profile)[2] # "Group_names"

trait <- read.csv(trait_profile,sep ="\t",header=TRUE) # 
traitCols <- names(trait)[3:length(trait)] #

df <- merge(profile,trait,by=c(sample_col,group_col))
df_samplecol <- names(df) %in% c(sample_col)
df <-df[!df_samplecol]


dd=datadist(df)
options(datadist="dd") 

## formula
imp <- read.csv(importance_profile,sep ="\t",header=TRUE) 
featureIDs <-imp[ ,names(imp)[1]]
featureCount <- length(featureIDs)
featureCount

tops <- seq(3,featureCount,3) # 
tops
topLabel1 <- glue("Top{tops[1]}_importance_features")
topLabel2 <- glue("Top{tops[2]}_importance_features")
topLabel3 <- glue("All{featureCount}_importance_features")

if (merge_trait==TRUE) {
        print("merge_trait")
        # merge trait
        formula_1 <- as.formula(paste(group_col,"~",paste(c(featureIDs[1:tops[1]],traitCols),collapse=" + "))) # top 5
        formula_2 <- as.formula(paste(group_col,"~",paste(c(featureIDs[1:tops[2]],traitCols),collapse=" + "))) # top 10
        formula_3 <- as.formula(paste(group_col,"~",paste(c(featureIDs[1:featureCount],traitCols),collapse=" + ")))
} else {
        print("not merge_trait")
        # no merge trait
        formula_1 <- as.formula(paste(group_col,"~",paste(featureIDs[1:tops[1]],collapse=" + "))) # top 5
        formula_2 <- as.formula(paste(group_col,"~",paste(featureIDs[1:tops[2]],collapse=" + "))) # top 10
        formula_3 <- as.formula(paste(group_col,"~",paste(featureIDs[1:featureCount],collapse=" + "))) # all
}

formula_3


# 
calibration_plot <- function(formula,profile,formulaName,outdir) {
        fit <- lrm(formula,data = profile,x=TRUE,y=TRUE)
        cal <- calibrate(fit,method="boot",B=500)
        
        pdf(file=glue("{outdir}/calibrated_method_1_{formulaName}.pdf"),width = 6.5,height = 4.5)
        plot(cal,
             xlim = c(0,1),
             ylab = "predicted Probability",
             legend = FALSE,
             subtitles = FALSE
        )
        abline(0,1,col="blue",lty=1,lwd=2)
        lines(cal[,c("predy","calibrated.orig")],type="l",lwd=2,col="red",pch=16)
        lines(cal[,c("predy","calibrated.corrected")],type="l",lwd=2,col="green",pch=16)
        legend(0.15,0.35,
               c("Apparent","Ideal","Bias-corrected"),
               lty=c(2,1,1),
               lwd=c(2,1,1),
               col = c("blue","red","green"),
               bty="n" # "O"为加边框
        )
        dev.off()
        
}

# method 1
calibration_plot(formula_1,df,topLabel1,outdir)
calibration_plot(formula_2,df,topLabel2,outdir)
calibration_plot(formula_3,df,topLabel3,outdir)

#
fit1 <- lrm(formula_1,data = df,x=TRUE,y=TRUE)
fit2 <- lrm(formula_2,data = df,x=TRUE,y=TRUE)
fit3 <- lrm(formula_3,data = df,x=TRUE,y=TRUE)

cal1 <- calibrate(fit1,method="boot",B=500)
cal2 <- calibrate(fit2,method="boot",B=500)
cal3 <- calibrate(fit3,method="boot",B=500)

pdf(file=paste(outdir,glue("calibrated_method_1_combine.pdf"),sep = ""),width = 6.5,height = 4.5)
plot(1,type="n",
     xlim = c(0,1),
     ylim = c(0,1),
     xlab = "Predicted Probability",
     ylab = "Observed Probability",
     legend = FALSE,
     subtitles = FALSE
)

abline(0,1,col="black",lty=1,lwd=2)
lines(cal1[,c("predy","calibrated.corrected")],type="l",lwd=2,col="green",pch=16)
lines(cal2[,c("predy","calibrated.corrected")],type="l",lwd=2,col="blue",pch=16)
lines(cal3[,c("predy","calibrated.corrected")],type="l",lwd=2,col="red",pch=16)
legend(0.15,0.35,
       c(topLabel1,topLabel2,topLabel3),
       lty=c(2,2,2),
       lwd=c(2,2,2),
       col = c("green","blue","red"),
       bty="n" 
)
dev.off()


# method 2
# library(riskRegression)
fit1 <- glm(formula_1,data = df,family = binomial())
fit2 <- glm(formula_2,data = df,family = binomial())
fit3 <- glm(formula_3,data = df,family = binomial())

xb <- Score(list("top_1_3_importance_features"=fit1,
                 "top_2_3_importance_features"=fit2,
                 "All_importance_features"=fit3
                 ),
            formula=Group_names~1,
            null.model=FALSE,
            conf.int=TRUE,
            plots=c("calibration","ROC"),
            metrics = c("auc","brier"),
            B=1000,
            M=50,
            data=df
            )
pdf(file=paste(outdir,"calibrated_method_2_combine.pdf",sep = ""),width = 6.5,height = 4.5)
plotCalibration(xb,col=c("blue","red","green"))
dev.off()

rm(list=ls())



