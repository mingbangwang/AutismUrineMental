rm(list=ls())

library(optparse)

option_list = list(
        make_option("--feature_profile", type="character", default=NULL, help="Input feature_profile file, the first two columns are Sample_names, Group_names. default: %default [required]" ),
        make_option("--trait_profile", type="character", default=NULL, help="Input trait_profile file, the first two columns are Sample_names, Group_names. default: %default [required]" ),
        make_option("--importance_profile", type="character", default=NULL, help="Input importance_profile file, the first column is feature name., second column is feature coef/IncNodePurity default: %default [required]" ),
        make_option("--traitCol", type="character", default=NULL, help="traitCol: %default [required]" ),
        make_option("--out",type="character", default=NULL, help="output dir, default: %default

  Example:
  Rscript clinical_decision.R --feature_profile=feature_profile --trait_profile=trait_profile --importance_profile=importance_profile --traitCol=traitCol --out=out_dir")

)

opt = parse_args(OptionParser(usage="%prog[option]", option_list=option_list, add_help_option=TRUE, prog="clinical_decision.R"), print_help_and_exit=TRUE, positional_arguments=FALSE)

options(digits = 4)
q_975 <- qnorm(0.975)

feature_profile = opt$feature_profile
trait_profile = opt$trait_profile
importance_profile = opt$importance_profile
traitCol = opt$traitCol
workdir = opt$out
merge_trait <- FALSE


library(rmda)
library(rms)
library(glue)

setwd(workdir) 
dir.create("clinical_decision") 
outdir <- 'clinical_decision/' 

profile <- read.csv(feature_profile,sep ="\t",header=TRUE) 

sample_col <- names(profile)[1] # "Sample_names" # 
group_col <- names(profile)[2] # "Group_names"

trait <- read.csv(trait_profile,sep ="\t",header=TRUE) 
traitCols <- names(trait)[3:length(trait)] #
traitCols

df <- merge(profile,trait,by=c(sample_col,group_col))
df_samplecol <- names(df) %in% c(sample_col)
df <-df[!df_samplecol]

dd=datadist(df)
options(datadist="dd") 

# 
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

model_1 <- decision_curve(formula_1,
                          data = df,
                          family = binomial(link = "logit"),
                          thresholds = seq(0,1,by = 0.01),
                          confidence.intervals = 0.95,
                          study.design = 'case-control',
                          population.prevalence = 0.3
                         )

model_2 <- decision_curve(formula_2,
                          data = df,
                          family = binomial(link = "logit"),
                          thresholds = seq(0,1,by = 0.01),
                          confidence.intervals = 0.95,
                          study.design = 'case-control',
                          population.prevalence = 0.3
                          )

model_3 <- decision_curve(formula_3,
                          data = df,
                          family = binomial(link = "logit"),
                          thresholds = seq(0,1,by = 0.01),
                          confidence.intervals = 0.95,
                          study.design = 'case-control',
                          population.prevalence = 0.3
                          )

plot_decision_curve(model_1,curve.names = c(topLabel1),xlim = c(0,0.8),
                    cost.benefit.axis = FALSE,col = c('red'),
                    confidence.intervals = FALSE,standardize = FALSE
                    )

plot_decision_curve(model_2,curve.names = c(topLabel2),xlim = c(0,0.8),
                    cost.benefit.axis = FALSE,col = c('red'),
                    confidence.intervals = FALSE,standardize = FALSE
                    )

plot_decision_curve(model_3,curve.names = c(topLabel3),xlim = c(0,0.8),
                    cost.benefit.axis = FALSE,col = c('red'),
                    confidence.intervals = FALSE,standardize = FALSE
                    )

# merge dca plot
model_all <- list(model_1,model_2,model_3)
# save fig
pdf(file=paste(outdir,"clinical_decision_curve_combine.pdf",sep = ""),width = 7,height = 5)
plot_decision_curve(model_all,curve.names = c(topLabel1,
                                              topLabel2,
                                              topLabel3),
                    xlim = c(0,0.8),
                    cost.benefit.axis = FALSE,col = c('red','green','blue','pink'),
                    confidence.intervals = FALSE,standardize = FALSE
                    ) 
dev.off()

# clinical impact curve
# model_1
pdf(file=paste(outdir,topLabel1,"_clinical_impact_model.pdf",sep = ""),width = 6,height = 4)
plot_clinical_impact(model_1,population.size = 1000,cost.benefit.axis = T,
                     n.cost.benefits = 8,col = c('red','blue'),
                     confidence.intervals = F
                     ) + title(topLabel1) # 
dev.off()

# model_2
pdf(file=paste(outdir,topLabel2,"_clinical_impact_model.pdf",sep = ""),width = 6,height = 4)
plot_clinical_impact(model_2,population.size = 1000,cost.benefit.axis = T,
                     n.cost.benefits = 8,col = c('red','blue'),
                     confidence.intervals = F
                     ) + title(topLabel2) #
dev.off()

# model_3
pdf(file=paste(outdir,topLabel3,"_clinical_impact_model.pdf",sep = ""),width = 6,height = 4)
plot_clinical_impact(model_3,population.size = 1000,cost.benefit.axis = T,
                     n.cost.benefits = 8,col = c('red','blue'),
                     confidence.intervals = F
                     ) + title(topLabel3) # 
dev.off()

