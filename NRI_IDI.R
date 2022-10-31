rm(list=ls())

# https://www.jianshu.com/p/3995e3312a90
# https://www.plob.org/article/12496.html

library(optparse)

option_list = list(
  make_option("--feature_profile", type="character", default=NULL, help="Input feature_profile file, the first two columns are Sample_names, Group_names. default: %default [required]" ),
  make_option("--trait_profile", type="character", default=NULL, help="Input trait_profile file, the first two columns are Sample_names, Group_names. default: %default [required]" ),
  make_option("--importance_profile", type="character", default=NULL, help="Input importance_profile file, the first column is feature name., second column is feature coef/IncNodePurity default: %default [required]" ),
  make_option("--traitCol", type="character", default=NULL, help="traitCol: %default [required]" ),
  # make_option(c("-m", "--merge_trait"), type = "logical", default = FALSE,
  #             help = "merge trait [default]"),
  make_option("--out",type="character", default=NULL, help="output dir, default: %default

  Example:
  Rscript NRI_IDI.R --feature_profile=feature_profile --importance_profile=importance_profile --importance_profile=importance_profile --traitCol=traitCol --out=out_dir")

)

opt = parse_args(OptionParser(usage="%prog[option]", option_list=option_list, add_help_option=TRUE, prog="TestHeatMap.R"), print_help_and_exit=TRUE, positional_arguments=FALSE)


feature_profile = opt$feature_profile
trait_profile = opt$trait_profile
importance_profile = opt$importance_profile
traitCol = opt$traitCol
workdir = opt$out
# merge_trait =opt$merge_trait
merge_trait <- FALSE

# # test 1
# feature_profile <- "D:/R_learning/Rscript/example/data/MRI/feature_num.tsv" # feature_num.tsv,all_diff_features_num.tsv
# trait_profile <- "D:/R_learning/Rscript/example/data/MRI/trait_clean_num.tsv"
# importance_profile <- "D:/R_learning/project/WES_301/lasso/lasso_coefs_top.txt" # lasso/lasso_coefs_top.txt,randomForest/importance.txt
# workdir <- 'D:/python_learning/biosoft/clinical_analysis/example/asd_MRI/clinical_analysis/' # "D:/python_learning/biosoft/clinical_analysis/example/asd_MRI/clincial_analysis/",'D:/R_learning/project/WES_301/
# traitCol <- "IQ" # "AI_Output"   "Age"         "Gender"      "IQ"          "DSM"
# merge_trait <- FALSE # 是否合并trait,FALSE,TRUE


# # test 2
# feature_profile <- "D:/python_learning/biosoft/clinical_analysis/example/data/asd_metagenome/all_diff_features_num.tsv" # feature_num.tsv,all_diff_features_num.tsv
# trait_profile <- "D:/python_learning/biosoft/clinical_analysis/example/data/asd_metagenome/trait_clean_num.tsv"
# importance_profile <- "D:/python_learning/biosoft/clinical_analysis/example/asd_metagenome/clinical_analysis/lasso/lasso_coefs_top.txt" # lasso/lasso_coefs_top.txt,randomForest/importance.txt
# workdir <- 'D:/python_learning/biosoft/clinical_analysis/example/asd_metagenome/clinical_analysis/' # "D:/python_learning/biosoft/clinical_analysis/example/asd_MRI/clincial_analysis/",'D:/R_learning/project/WES_301/
# traitCol <- "IgA" #  "Age"         "Gender"
# merge_trait <- FALSE # 是否合并trait,FALSE,TRUE


# 设置有效数字位数
options(digits = 4)
q_975 <- qnorm(0.975)


## 载入R包
library(nricens) # NRI
library(PredictABEL)
library(glue)
library(rms)

# 设置工作文件夹
setwd(workdir) #进入文件夹
dir.create("nri_idi") #创建文件夹
outdir <- 'nri_idi/' #输出文件夹

# 读取数据
profile <- read.csv(feature_profile,sep ="\t",header=TRUE) # profile_df_trait_diversity_num.csv
# head(profile)

# 获取样本名称和组名列
sample_col <- names(profile)[1] # "Sample_names" # 
group_col <- names(profile)[2] # "Group_names"

trait <- read.csv(trait_profile,sep ="\t",header=TRUE) # profile_df_trait_diversity_num.csv
traitCols <- names(trait)[3:length(trait)] #
# traitCols <- c("AI_Output","Age" ,"Gender" ,"IQ")
traitCols

df <- merge(profile,trait,by=c(sample_col,group_col))
# 使用过滤不要的"Sample_names"列
df_samplecol <- names(df) %in% c(sample_col)
df <-df[!df_samplecol]

## 关键步骤：按照nomogram要求“打包”数据，可以输入??datadist查看详细说明
dd=datadist(df)
options(datadist="dd") 

# 读取重要的feature
imp <- read.csv(importance_profile,sep ="\t",header=TRUE) # profile_df_trait_diversity_num.csv
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
    # 引入trait
    formula_1 <- as.formula(paste(group_col,"~",paste(c(featureIDs[1:tops[1]],traitCols),collapse=" + "))) # top 5
    formula_2 <- as.formula(paste(group_col,"~",paste(c(featureIDs[1:tops[2]],traitCols),collapse=" + "))) # top 10
    formula_3 <- as.formula(paste(group_col,"~",paste(c(featureIDs[1:featureCount],traitCols),collapse=" + ")))
} else {
    print("not merge_trait")
    # 不引入trait
    formula_1 <- as.formula(paste(group_col,"~",paste(featureIDs[1:tops[1]],collapse=" + "))) # top 5
    formula_2 <- as.formula(paste(group_col,"~",paste(featureIDs[1:tops[2]],collapse=" + "))) # top 10
    formula_3 <- as.formula(paste(group_col,"~",paste(featureIDs[1:featureCount],collapse=" + "))) # all
}

formula_3

## 设置生存公式
# Surv两个变量
timeCol <-  traitCol  # 时间，或者重要的表型
eventCol <- group_col  # 发生的事件，0，1变量

# 设置生存公式的Y
Y <- glue("Surv({timeCol},{eventCol})")
# Y <- glue("Surv(Precentral_R,Group_names)")
# Y 

if (merge_trait==TRUE) {
    print("merge_trait")
    # 引入trait
    # 去掉trait里面作为timeCol的变量
    traitList <- as.list(traitCols) # traitCols = c("AI_Output","IQ")
    for (i in 1:length(traitList)) {
        if (traitList[i] == timeCol){
            print (traitList[i])
            traitList[i] = NULL
        }
    }
    traitList
    X1 <- paste(c(featureIDs[1:tops[1]],traitList),collapse=" + ") # traitList,traitCols
    X2 <- paste(c(featureIDs[1:tops[2]],traitList),collapse=" + ") # traitList,traitCols
    X3 <- paste(c(featureIDs[1:featureCount],traitList),collapse=" + ") # traitList,traitCols
    
} else {
    print("not merge_trait")
    # 设置其他协变量,不引入trait,成功
    X1 <- paste(featureIDs[1:tops[1]],collapse=" + ") # length(featureIDs)
    X2 <- paste(featureIDs[1:tops[2]],collapse=" + ") # length(featureIDs)
    X3 <- paste(featureIDs[1:featureCount],collapse=" + ") # length(featureIDs)
}

X3

# 得到生存公式公式
survFormula1 <- as.formula(glue("{Y} ~ {X1}"))
survFormula2 <- as.formula(glue("{Y} ~ {X2}"))
survFormula3 <- as.formula(glue("{Y} ~ {X3}"))
survFormula3

# surv_formula
model1 <- coxph(survFormula1, data = df)
model2 <- coxph(survFormula2, data = df)
model3 <- coxph(survFormula3, data = df)


## 针对上述三个模型进行模型自身及模型间的评价与验证
## 净重新分类指数(Net ReclassificationIndex, NRI)
## 手动计算P值（两两比较）

# 记录输出
con_nri <- file(glue("{workdir}/{outdir}/compare_model_log_NRI.txt")) # 创建一个.log文件
sink(con_nri, append=TRUE) # 记录output
sink(con_nri, append=TRUE, type="message") # 记录message
# 所有的output和message都会记录到test.log中，而控制台中不在有信息显示

# eventCol <- group_col # "Group_names"
event <- df[ , c(eventCol)]
z_std <- as.matrix(subset(profile,select = featureIDs[1:tops[1]])) # 
z_new1 <- as.matrix(subset(profile,select = featureIDs[1:tops[2]]))
z_new2 <- as.matrix(subset(profile,select = featureIDs[1:featureCount]))

mstd <- glm(event~.,binomial(logit),data.frame(event,z_std),x=TRUE)
mnew1 <- glm(event~.,binomial(logit),data.frame(event,z_new1),x=TRUE)
mnew2 <- glm(event~.,binomial(logit),data.frame(event,z_new2),x=TRUE)

# 截断值的选择至关重要
set.seed(1)
print(glue("compare top {tops[1]} vs top {tops[2]}"))
result <- nribin(mdl.std = mstd,# mstd
                 mdl.new = mnew1,# mnew1
                 cut = c(0.2,0.4),niter = 1000,updown = "category")

# Point & Interval estimates:
#   Estimate  Std.Error      Lower      Upper
# NRI           0.34415584 0.14426800 0.12362177 0.68558764
# NRI+          0.07142857 0.04735891 0.00000000 0.18518519
# NRI-          0.27272727 0.12277085 0.09523810 0.56097561
# Pr(Up|Case)   0.12500000 0.04430605 0.03571429 0.21153846
# Pr(Down|Case) 0.05357143 0.02295963 0.00000000 0.08474576
# Pr(Down|Ctrl) 0.34090909 0.10943827 0.18604651 0.60495716
# Pr(Up|Ctrl)   0.06818182 0.04113631 0.00000000 0.15047170

# z <- abs(0.34415584  /0.14426800  ) # Estimate/Std.Error
# pvalue <- (1-pnorm(z))*2 # 0.01705447,top3 vs top 10, 说明模型改善大
# print(pvalue) # 

print(
    "z <- abs(Estimate/Std.Error) \n
    pvalue <- (1-pnorm(z))*2 \n
    "    
)

print("--------------------------")

print(glue("compare top {tops[1]} vs all"))
nribin(mdl.std = mstd, # mstd
       mdl.new = mnew2, # mnew2
       cut = c(0.2,0.4),niter = 1000,updown = "category")
# pvalue <- (1-pnorm(z))*2 # 2.293883e-09,top3 vs all,说明模型改善大
# print(pvalue) # 2.293883e-09，说明模型改善大

# z <- abs(0.7159091/0.14359608) # Estimate/Std.Error
# pvalue <- (1-pnorm(z))*2 # 6.177765e-07-09,top5 vs all,说明模型改善大
# print(pvalue) #
print("--------------------------")

print(glue("compare top {tops[2]} vs all"))
nribin(mdl.std = mnew1, # mnew1
       mdl.new = mnew2, # mnew2
       cut = c(0.2,0.4),niter = 1000,updown = "category")

# z <- abs(0.4935065/0.11965690) # Estimate/Std.Error
# pvalue <- (1-pnorm(z))*2 # 3.717891e-05,top10 vs all,说明模型改善大
# print(pvalue) #
print("--------------------------")


# 连续净重新分类指数（cfNRI）
print("cfNRI")
set.seed(1)
print(glue("compare top {tops[1]} vs top {tops[2]}"))
nribin(mdl.std = mstd,
       mdl.new = mnew1,
       cut = 0,niter = 1000,updown = "diff")
print("--------------------------")

print(glue("compare top {tops[1]} vs all"))
nribin(mdl.std = mstd,
       mdl.new = mnew2,
       cut = 0,niter = 1000,updown = "diff")
print("--------------------------")

print(glue("compare top {tops[2]} vs all"))
nribin(mdl.std = mnew1,
       mdl.new = mnew2,
       cut = 0,niter = 1000,updown = "diff")
print("--------------------------")

# NRI记录完毕后，重置output和message的记录，运行完一下两行，后续的输入命令重新显示到控制台中
sink()
sink(type="message")


## IDI（Integrated Discrimination Improvement，综合判别改善指数）
# 记录输出
con <- file(glue("{workdir}/{outdir}/compare_model_log_IDI.txt")) # 创建一个.log文件
sink(con, append=TRUE) # 记录output
sink(con, append=TRUE, type="message") # 记录message
# 所有的output和message都会记录到test.log中，而控制台中不在有信息显示
print("Integrated Discrimination Improvement")

# 同时有NRI
# IDI的有P值
pstd <- mstd$fitted.values
pnew1 <- mnew1$fitted.values
pnew2 <- mnew1$fitted.values

demo <- as.matrix(df) # 
print(glue("compare top {tops[1]} vs top {tops[2]}"))
reclassification(data = demo,
                 cOutcome = eventCol,
                 predrisk1 = pstd,
                 predrisk2 = pnew1,
                 cutoff = c(0,0.2,0.4,1))
# NRI(Categorical) [95% CI]: 0.2581 [ 0.062 - 0.4542 ] ; p-value: 0.00989 # 改善大
# NRI(Continuous) [95% CI]: 0.9545 [ 0.6071 - 1.302 ] ; p-value: 0 
# IDI [95% CI]: 0.1703 [ 0.0951 - 0.2455 ] ; p-value: 1e-05 
print("--------------------------")

print(glue("compare top {tops[1]} vs all"))
reclassification(data = demo,
                 cOutcome = eventCol,
                 predrisk1 = pstd,
                 predrisk2 = pnew2,
                 cutoff = c(0,0.2,0.4,1))
# NRI(Categorical) [95% CI]: 0.2581 [ 0.062 - 0.4542 ] ; p-value: 0.00989 # 改善大
# NRI(Continuous) [95% CI]: 0.9545 [ 0.6071 - 1.302 ] ; p-value: 0 
# IDI [95% CI]: 0.1703 [ 0.0951 - 0.2455 ] ; p-value: 1e-05 
print("--------------------------")

print(glue("compare top {tops[2]} vs all"))
reclassification(data = demo,
                 cOutcome = eventCol,
                 predrisk1 = pnew1,
                 predrisk2 = pnew2,
                 cutoff = c(0,0.2,0.4,1))
# NRI(Categorical) [95% CI]: 0 [ 0 - 0 ] ; p-value: NaN # 改善大
# NRI(Continuous) [95% CI]: 0 [ 0 - 0 ] ; p-value: NaN 
# IDI [95% CI]: 0 [ 0 - 0 ] ; p-value: NaN 
print("--------------------------")

# 记录完毕后，重置output和message的记录，运行完一下两行，后续的输入命令重新显示到控制台中
sink()
sink(type="message")

