rm(list=ls())

library(optparse)

option_list = list(
  make_option("--feature_profile", type="character", default=NULL, help="Input feature_profile file, the first two columns are Sample_names, Group_names. default: %default [required]" ),
  make_option("--out",type="character", default=NULL, help="output dir, default: %default

  Example:
  Rscript lasso.R --feature_profile=feature_profile --out=out_dir")

)

opt = parse_args(OptionParser(usage="%prog[option]", option_list=option_list, add_help_option=TRUE, prog="TestHeatMap.R"), print_help_and_exit=TRUE, positional_arguments=FALSE)

options(digits = 22)
q_975 <- qnorm(0.975)

feature_profile = opt$feature_profile
workdir = opt$out

library(corrplot)
library(glmnet)
library(glue)
library(survival)

setwd(workdir) 
dir.create("lasso") 
outdir <- 'lasso/' 

profile <- read.csv(feature_profile,sep='\t', header=TRUE) 
profile <- na.omit(profile)

sample_col <- names(profile)[1]
group_col <- names(profile)[2] 


profile_samplecol <- names(profile) %in% sample_col
profile <-profile[!profile_samplecol]


profile_groupcol <- names(profile) %in% c(group_col)
profile_nogroup <-profile[!profile_groupcol]
X <- as.matrix(profile_nogroup)
Y <- profile[ , c(group_col)]

lambdas <- seq(0,0.5,length.out = 200)
set.seed(123)
cv.lasso <- cv.glmnet(X,Y,alpha = 1,lambda = lambdas,nfolds = 3,family = "binomial")
lasso_lse <- cv.lasso$lambda.min # min
print(glue("best lambda: {lasso_lse}")) 

# plot cv lasso 
pdf(file=paste(outdir,"lasso_cv.pdf",sep = ""),width = 5,height = 4)
plot(cv.lasso)
abline(v=log(c(cv.lasso$lambda.min,cv.lasso$lambda.1se)),col=c("blue","red"),lty=c(2,2),lwd=c(2,2))
dev.off()

# plot lambda
pdf(file=paste(outdir,"lasso_cv_lambda.pdf",sep = ""),width = 5,height = 4)
plot(cv.lasso$glmnet.fit,xvar = "lambda",label = F)  
abline(v=log(cv.lasso$lambda.min),lwd=2,col="blue",lty=2)
dev.off()


# save lasso formula
lasso.coef <- coef(cv.lasso$glmnet.fit, s = lasso_lse,exact = F) # 


lasso_coef <- as.matrix(lasso.coef)
lasso_coef <- as.data.frame(lasso_coef)

lasso_coef_file <- paste(outdir,"lasso_coefs.txt",sep="")
lasso_coef_top_file <- paste(outdir,"lasso_coefs_top.txt",sep="")

write.table(lasso_coef,file = lasso_coef_file,quote = F,sep = '\t', row.names = T, col.names = F)

coefs <- read.csv(lasso_coef_file,sep='\t', header=T,col.names = c("id","value"))

coefs_top <- coefs[coefs$value >0 |coefs$value <0 ,]
coefs_top <- coefs_top[order(abs(coefs_top[,"value"]),decreasing = T), ] # 

write.table(coefs_top,file = lasso_coef_top_file,quote = F,sep = '\t', row.names = T, col.names = T)

barplot_lasso <- function(data_file,idCol, valueCol,name,outdir){
        df <- read.csv(data_file,sep = "\t",header=TRUE) # profile_df_trait_diversity_num.csv
        values <- df[,valueCol]
        cols <-c()
        for (i in values) {
                if (i < 0) { 
                        cols <-c(cols,"blue")
                }
                else {
                        cols <-c(cols,"red")
                }
        }

        pdf(file=glue("{outdir}/barplot_{name}.pdf"),width = 7,height = 5)
        p <- barplot(values, # valueCol
                     col=cols,
                     beside=TRUE,
                     horiz=TRUE # TRUE
        )
        ids <- as.vector(df[,idCol])
        text(p, ids, label = tmp, pos = 4) #                                                               
        dev.off()
        
}
# barplot
barplot_lasso(data_file = glue("{outdir}/lasso_coefs_top.txt"),
              idCol="id", 
              valueCol="value",
              name="lasso_coefs_top",
              outdir)

# formula_lasso
formula_lasso <- as.formula(paste(group_col,"~",paste(coefs_top$id,collapse=" + ")))


# evaluate formula
model_lasso <- glm(formula_lasso,data =  profile,family = binomial()) #formula
formula_evaluate1 <- summary(model_lasso)$coefficients
formula_evaluate2 <- exp(cbind("OR"=coef(model_lasso),confint(model_lasso)))
write.table(cbind(formula_evaluate1,formula_evaluate2),file = paste(outdir,"formula_evaluate.tsv",sep = ""),quote = F,sep = '\t', row.names = T, col.names = T)


rm(list=ls())

