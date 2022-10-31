rm(list=ls())


library(optparse)

option_list = list(
        make_option("--feature_profile", type="character", default=NULL, help="Input feature_profile file, the first two columns are Sample_names, Group_names. default: %default [required]" ),
        make_option("--trait_profile", type="character", default=NULL, help="Input trait_profile file, the first two columns are Sample_names, Group_names. default: %default [required]" ),
        make_option("--importance_profile", type="character", default=NULL, help="Input importance_profile file, the first column is feature name., second column is feature coef/IncNodePurity default: %default [required]" ),
        make_option("--traitCol", type="character", default=NULL, help="traitCol: %default [required]" ),
        make_option("--out",type="character", default=NULL, help="output dir, default: %default

  Example:
Rscript nomogram_reg.R --feature_profile=feature_profile --trait_profile=trait_profile --importance_profile=importance_profile --traitCol=traitCol --out=out_dir")

)

opt = parse_args(OptionParser(usage="%prog[option]", option_list=option_list, add_help_option=TRUE, prog="nomogram_reg.R"), print_help_and_exit=TRUE, positional_arguments=FALSE)

feature_profile = opt$feature_profile
trait_profile = opt$trait_profile
importance_profile = opt$importance_profile
traitCol = opt$traitCol
workdir = opt$out
merge_trait <- TRUE


options(digits = 2) # 4
q_975 <- qnorm(0.975)


library(rms)
library(survival)
library(survminer)
library(regplot)
library(nomogramEx) # for total points
library(glue)

setwd(workdir) #
dir.create("nomogram") #
outdir <- 'nomogram/' #

profile <- read.csv(feature_profile,sep ="\t",header=TRUE) 

sample_col <- names(profile)[1] # "Sample_names" # 
group_col <- names(profile)[2] # "Group_names"

trait <- read.csv(trait_profile,sep ="\t",header=TRUE) 
traitCols <- names(trait)[3:length(trait)] #


df <- merge(profile,trait,by=c(sample_col,group_col))
# 使用过滤不要的"Sample_names"列
df_samplecol <- names(df) %in% c(sample_col)
df <-df[!df_samplecol]
        
dd=datadist(df)
options(datadist="dd") 

## for formula
imp <- read.csv(importance_profile,sep ="\t",header=TRUE) # profile_df_trait_diversity_num.csv
featureIDs <-imp[ ,names(imp)[1]]
featureCount <- length(featureIDs)

tops <- seq(3,featureCount,3) # 
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


## define function of logistic regression model for nomogram
nomogram_lrm_plot <- function(formula,data,formulaName,outdir){
        
        ## 建logistic regression model
        model_lrm <- lrm(formula, data)
        nom_lrm <- nomogram(model_lrm, # model name
                            fun=plogis, 
                            lp=TRUE, # 
                            fun.at=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),# 
                            funlabel="Risk") # 
        
        # 计算total points
        table <- nomogramEx(nom_lrm)
        results <- print(table,
                          showAllLevels=FALSE) # FALSE
        write.csv(results,glue("{outdir}/total_points_nom_lrm_{formulaName}.csv"))
        
        # save fig
        pdf(file=glue("{outdir}/nomogram_{formulaName}.pdf"),width = 8,height = 8)
        plot(nom_lrm,
             lplabel = 'Linear Predictor',
             col.grid = c("blue","yellow"), # 
             fun.side = c(1,3,1,3,1,3,1,3,1), # 
             xfrac=.3,
             
        )
        
        dev.off()
        
}


nomogram_lrm_plot(formula=formula_1,data=df,formulaName=topLabel1,outdir=outdir)
nomogram_lrm_plot(formula=formula_2,data=df,formulaName=topLabel2,outdir=outdir)
nomogram_lrm_plot(formula=formula_2,data=df,formulaName=topLabel3,outdir=outdir)


## define function of COXPH for nomogram
nomogram_coxph_plot <- function(formula,profile,formulaName,time,outdir){
        model_coxph <- psm(formula=formula,data =profile, dist='lognormal') # 
        
        med <- Quantile(model_coxph) # 
        surv <- Survival(model_coxph) # 
        nom_cox <- nomogram(model_coxph,
                            fun=function(x) med(lp=x),
                            funlabel=paste("Median",time, "value",sep = " ")) # Median Survival Time
        # for total points
        table<- nomogramEx(nom_cox)
        results <- print(table,
                         showAllLevels=FALSE) # FALSE
        write.csv(results,glue("{outdir}/total_points_nom_coxPH_{formulaName}.csv"))
        
        # save fig
        pdf(file=glue("{outdir}/nomogram_coxph_{formulaName}.pdf"),width = 8,height = 8)
        plot(nom_cox,
             lplabel = 'Linear Predictor',
             col.grid = c("blue","yellow"), 
             xfrac=.3,
             
        )
        
        dev.off()
          
}


## 构建COX比例风险模型
# COX回归中位生存时间的Nomogram
# Surv两个变量
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

# X3


# 得到生存公式公式
survFormula1 <- as.formula(glue("{Y} ~ {X1}"))
survFormula2 <- as.formula(glue("{Y} ~ {X2}"))
survFormula3 <- as.formula(glue("{Y} ~ {X3}"))
# survFormula3


nomogram_coxph_plot(formula=survFormula1,
                    profile=df, # 
                    formulaName=glue("{topLabel1}_{timeCol}"),
                    time=timeCol,
                    outdir=outdir)

nomogram_coxph_plot(formula=survFormula2,
                    profile=df, #
                    formulaName=glue("{topLabel2}_{timeCol}"),
                    time=timeCol,
                    outdir=outdir)

nomogram_coxph_plot(formula=survFormula3,
                    profile=df, #
                    formulaName=glue("{topLabel3}_{timeCol}"),
                    time=timeCol,
                    outdir=outdir)



## 方法二、regplot: Plots a regression nomogram showing covariate distribution，需要手工完成
# library(regplot)

# 先计算timeCol的取值区间
timeCol.min <- min(df[,timeCol])
timeCol.1st_Qu <- quantile(df[,timeCol],1/4)
timeCol.med <- median(df[,timeCol])
timeCol.mean <- mean(df[,timeCol])
timeCol.3rd_Qu <- quantile(df[,timeCol],3/4)
timeCol.max <- max(df[,timeCol])

# A lrm
lrm_model <- lrm(formula_3,data =  df) #

norm_lrm <- regplot(lrm_model,
        observation = df[1,], # 指定第几个观测值，选择一个患者的列名试试（如2,11）
        title = glue("nomogram of {timeCol} vs {eventCol} using lrm"), # logistic regression model
        points = TRUE, # Point的最大刻度为100
        odds = FALSE, # 设置是否设置OR
        rank = "sd", # 按照回归系数的SD变量排序
        clickable=FALSE, # FALSE,是否可以点击，进行交互
        center = TRUE, # TRUE,将各个变量设置不从0开始
        showP = TRUE, # 显示变量是否存在统计学意义
        # interval = "confidence", # 置信区间
        )



# B coxPH
coxph_model <-  coxph(formula =survFormula3,data =df) #

regplot(coxph_model,
        observation=df[37,], # 指定第几个观测值，选择一个患者的列名试试（如71）
        clickable=FALSE, # TRUE
        points=TRUE,
        rank="sd",
        showP = TRUE, # 显示变量是否存在统计学意义
        center = TRUE, # TRUE,将各个变量设置不从0开始
        # interval = "confidence", # 置信区间
        title = glue("nomogram of {timeCol} vs {eventCol}  CoxPH model"), # Proportional Hazards model
        failtime = c(timeCol.3rd_Qu,timeCol.med,timeCol.1st_Qu), # 重要变量的阈值,failtime should be in range of observed: 0.3218 0.6644
        )

rm(list=ls())
