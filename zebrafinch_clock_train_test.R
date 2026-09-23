################ loading packages #############

library(caret) #package for machine learning training
library(glmnet) #package for elastic net/ridge/lasso regression
library(doMC) #package for running loops in parallel on different cores
library(glmnetUtils) #additional package for elastic net regression
library(tidyverse)
library(readxl)
library(ggplot2)
library(GenomicRanges)
library(purrr)
library(dplyr)


############## set up the environment ############

args <- commandArgs(trailingOnly = TRUE) # otherwise Eddie server will bump into errors
sample_file <- args[1] 
r2 <- as.numeric(args[2]) #0  default should have been 0 
prefix_text <- args[3] #''
training_pct <- 1
adult_cutoff <- 0
out_dir <- args[4] 
meth_data_dir <- './meth_corr/' 
correlation_dir <- './meth_corr/' 
file <- ''
cutoff <- as.numeric(args[5])
surplus <- as.numeric(args[6])
sample_file_train <- args[7]
sample_file_test <- args[8]

############## setting parameters ##############

lambda_grid <- 10 ^ seq(-2,2,length = 81) 
alpha_grid <-  0.5 
grid2 <- expand.grid(alpha = alpha_grid, lambda = lambda_grid)

training_percentage <- training_pct  
prefix <- prefix_text  
r2_cutoff <- r2

prefix2 <- gsub('r201','no_filter',prefix)
meth_data2 <- readRDS(paste0(meth_data_dir,'/',prefix2,"_meth_data.RDS"))


################## load metadata ################
sample_data <- read.csv(sample_file)
sample_data <- sample_data[order(sample_data$Sample), ]
sample_data$Age <- as.numeric(sample_data$Age)
sample_data$Bird <- gsub("_.*", "",sample_data$Sample.ID)


################## processing sample data ###############

prefix3 <- gsub('no_filter','r201',prefix)
if(r2_cutoff > 0){
 correlation_results <- readRDS(paste0(correlation_dir,'/',prefix3,'_correlation_results.RDS'))
  keep_clusters <- correlation_results[correlation_results$R_squared >= r2_cutoff,]$cluster
  keep_clusters <- intersect(colnames(meth_data2),keep_clusters) 
}else{
  keep_clusters <- colnames(meth_data2)
}
keep_clusters <- keep_clusters[keep_clusters != 'Age']
keep_clusters <- keep_clusters[keep_clusters != 'Sex_dummy'] 

meth_data3 <- meth_data2[,colnames(meth_data2) %in% keep_clusters]

meth_data3 <- meth_data3[sample_data$Sample.ID,]
identical(sample_data$Sample.ID, rownames(meth_data3))
meth_data <- meth_data3


######################## logarithmic functions ################

young_transform <- function(x,adult.age= cutoff, surplus = surplus) {
  x=(x+surplus)/(adult.age + surplus) 
  y=ifelse(x<=1, log(x),x-1)
  return(y)
}

relative_age <- function(x,adult.age= cutoff, surplus = surplus) {
  x=(x+surplus)/(adult.age + surplus) 
  y=ifelse(x<=1, log(x),x-1)
  y = y*(adult.age +surplus) + adult.age 
  return(y)
}

relative_age_f <- function(f, adult.age= cutoff, surplus = surplus){
  y = f*(adult.age +surplus) + adult.age 
  return(y)
}

anti_young_transform <- function(f, adult.age= cutoff, surplus = surplus){
  y=  ifelse(f<0, (adult.age + surplus)*exp(f)-surplus, (adult.age + surplus)*f+adult.age)
  return(y)
}


############## load sample files ################

sample_data_train  <- read.csv(sample_file_train)
sample_data_train <- sample_data_train[order(sample_data_train$Sample.ID), ]
sample_data_train$Age <- as.numeric(sample_data_train$Age)
  
surplus_grid <- surplus


############## before training, do age transformation ####################

results <- data.frame(file=character(), RMSE_cv = numeric(), MAE_cv = numeric(),RMSE_cv_relative = numeric(),cutoff = numeric(), surplus = numeric(),
                      MAE_train_F=numeric(), RMSE_train_F = numeric(), r_sq_train_F=numeric(),
                      MAE_train_age=numeric(), RMSE_train_age = numeric(), r_sq_train_age=numeric(),
                      MAE_train_relative_age=numeric(), RMSE_train_relative_age = numeric(), r_sq_train_relative_age=numeric(),
                      MAE_test=numeric(),  RMSE_test =numeric(), r_sq_test=numeric(), 
                      alpha = numeric(), lambda = numeric(), total_par = numeric(), non_zero_par =numeric(),stringsAsFactors=FALSE)


##################### model training and testing ###################

rm(meth_data2, meth_data3)

if(file.exists(sample_file_train)){}else{ 
  
  meth_data_train <- meth_data[trainIndex,]
  meth_data_test <- meth_data[-trainIndex,]
}

meth_data_train <- meth_data[sample_data_train$Sample.ID,]
identical(rownames(meth_data_train), sample_data_train$Sample.ID )
meth_data_test <- meth_data[sample_data_test$Sample.ID,]

meth_data_train_full_age <- meth_data_train
sample_data_train_full_age <- sample_data_train

sample_data_train <- sample_data_train[sample_data_train$Age >= adult_cutoff,] # age filtering
adult_samples <- sample_data_train$Sample.ID
meth_data_train <- meth_data_train[sample_data_train$Sample.ID,]
identical(rownames(meth_data_train), sample_data_train$Sample.ID )

sample_data_train$Bird <- gsub("_.*", "",sample_data_train$Sample.ID)
out_index <- split(1:nrow(sample_data_train), sample_data_train$Bird)
cv_index <- list()
for (i in seq_along(out_index)) {
  # Get the complement of the current split
  complement <- setdiff(1:nrow(sample_data_train), out_index[[i]])
  
  # Store the complement in cv_index
  cv_index[[i]] <- complement
  names(cv_index)[[i]] <- names(out_index)[[i]]
}


### Run elastic net model and find the best parameters for the lasso regression ###

sample_data_train$Fstat <- young_transform (sample_data_train$Age,adult.age= cutoff, surplus = surplus) 
sample_data_train$relative_age <- relative_age (sample_data_train$Age,adult.age= cutoff, surplus = surplus) 

control <- trainControl(method = "repeatedcv",number = length(cv_index),
                        index = cv_index,indexOut = out_index,
                        repeats = 1, 
                        search = "grid",
                        verboseIter = TRUE,
                        allowParallel = TRUE,
                        returnData = TRUE,
                        trim = TRUE)

garbage <- capture.output(elastic_model_train2 <- train(meth_data_train, sample_data_train$Fstat,
                                                        data = NULL,
                                                        method = "glmnet",
                                                        tuneLength = 15,#nrow(grid), #how many combinations to try - already has grid
                                                        tuneGrid = grid2,
                                                        trControl = control))
pdf(paste0(out_dir,'/',prefix,'_cutoff',cutoff,'_surplus',surplus,'_model2_diagnostic.pdf'))
print(plot(elastic_model_train2,  scales = list(x = list(log = 10)),  main = list(paste0('model2_diagnostic ', prefix,'_cutoff',cutoff,'_surplus',surplus), cex=0.8) )) # this package uses lattice plot package to do the
dev.off()

write.csv(elastic_model_train2[["results"]],file = paste0(out_dir,'/',prefix,'_cutoff',cutoff,'_surplus',surplus,'_model2.csv'))
saveRDS(elastic_model_train2, file = paste0(out_dir,'/',prefix,'_cutoff',cutoff,'_surplus',surplus,'_model2.RDS'))
best_alpha2 <- elastic_model_train2$bestTune$alpha
best_lambda2 <- elastic_model_train2$bestTune$lambda

cv_res <- elastic_model_train2[["results"]]
cv_res_sorted <- cv_res[order(cv_res$RMSE, -cv_res$lambda), ]
RMSE_cv <-  cv_res_sorted$RMSE[1]
MAE_cv <-    cv_res_sorted$MAE[1]
RMSE_cv_relative <- RMSE_cv * (cutoff + surplus)
write.csv(cv_res, file = paste0(out_dir,'/',prefix,'_cutoff',cutoff,'_surplus',surplus,"_cv_results.csv"))


# extract weights of cpg sites

weights1 <- data.frame(as.matrix( coef(elastic_model_train2$finalModel, s =best_lambda2 ))) 
saveRDS(weights1, file = paste0(out_dir,'/',prefix,'_cutoff',cutoff,'_surplus',surplus,"_weights_interpolated.RDS"))
write.csv(weights1, file = paste0(out_dir,'/',prefix,'_cutoff',cutoff,'_surplus',surplus,"_weights_interpolated.csv"))

model3 = glmnet(x = meth_data_train, y= sample_data_train$Fstat, lambda=best_lambda2, alpha=best_alpha2)#running GLMNET using lamda from the cv procedure
weights = data.frame(coef.name = dimnames(coef(model3))[[1]], coef.value = matrix(coef(model3)))#These functions weight the variable x by a specific vector of weights.

saveRDS(weights, file = paste0(out_dir,'/',prefix,'_cutoff',cutoff,'_surplus',surplus,"_weights.RDS"))
write.csv(weights, file = paste0(out_dir,'/',prefix,'_cutoff',cutoff,'_surplus',surplus,"_weights.csv"))
saveRDS(model3, file = paste0(out_dir,'/',prefix,'_cutoff',cutoff,'_surplus',surplus,"_model3.RDS"))

ggplot(weights[-1,], aes(x = coef.value))+geom_histogram(fill = 'pink3')+ scale_y_continuous(trans = 'log10')+theme_classic()+ggtitle(paste0("Model parameters ",prefix, ' ', file))
ggsave(paste0(out_dir,'/',prefix,'_cutoff',cutoff,'_surplus',surplus,"_model_parameters.pdf"), height=7, width=10)

total_par <- ncol(meth_data)
non_zero_par <- sum(weights$coef.value != 0) -1



 ################################# evaluation of model ################################################

predicted_train_F <- model3 %>% predict(as.matrix(meth_data_train))
corr<-cor.test(sample_data_train$Fstat, predicted_train_F, alternative = c("two.sided", "less", "greater"),
               method = c("pearson", "kendall", "spearman"),
               exact = NULL, conf.level = 0.95, continuity = FALSE)

r_sq_train_F <- (corr[["estimate"]][["cor"]])^2
RMSE_train_F <- RMSE(sample_data_train$Fstat, predicted_train_F)
MAE_train_F <- MAE(sample_data_train$Fstat, predicted_train_F)

predicted_train_age <- anti_young_transform(f = predicted_train_F, adult.age= cutoff, surplus = surplus)
corr<-cor.test(sample_data_train$Age, predicted_train_age, alternative = c("two.sided", "less", "greater"),
               method = c("pearson", "kendall", "spearman"),
               exact = NULL, conf.level = 0.95, continuity = FALSE)

r_sq_train_age <- (corr[["estimate"]][["cor"]])^2
RMSE_train_age <- RMSE(sample_data_train$Age, predicted_train_age)
MAE_train_age <- MAE(sample_data_train$Age, predicted_train_age)

predicted_train_relative_age <- relative_age_f(f = predicted_train_F, adult.age= cutoff, surplus = surplus)
corr<-cor.test(sample_data_train$relative_age, predicted_train_relative_age, alternative = c("two.sided", "less", "greater"),
               method = c("pearson", "kendall", "spearman"),
               exact = NULL, conf.level = 0.95, continuity = FALSE)

r_sq_train_relative_age <- (corr[["estimate"]][["cor"]])^2
RMSE_train_relative_age <- RMSE(sample_data_train$relative_age, predicted_train_relative_age)
MAE_train_relative_age <- MAE(sample_data_train$relative_age, predicted_train_relative_age)

predicted_ages_df <- data.frame(predicted_train_age = predicted_train_age, predicted_train_F = predicted_train_F, predicted_train_relative_age = predicted_train_relative_age)
predicted_ages_df$actual_age <- sample_data_train$Age
predicted_ages_df$relative_age <- sample_data_train$relative_age
predicted_ages_df$Fstat <- sample_data_train$Fstat
predicted_ages_df$Sample <- sample_data_train$Sample.ID
predicted_ages_df$run <- prefix
write.csv(predicted_ages_df, file=paste0(out_dir,'/',prefix,'_cutoff',cutoff,'_surplus',surplus,"_predicted_ages_train.csv"), row.names=FALSE)


r_sq_test <- NA
RMSE_test  <- NA
MAE_test <- NA

#add to the results dataframe
results <- rbind(results, data.frame(file=file,RMSE_cv = RMSE_cv , MAE_cv = MAE_cv, RMSE_cv_relative = RMSE_cv_relative,cutoff =cutoff, surplus = surplus,
                                     MAE_train_F=MAE_train_F,RMSE_train_F=RMSE_train_F,
                                     r_sq_train_F=r_sq_train_F,
                                     MAE_train_age=MAE_train_age,RMSE_train_age=RMSE_train_age,
                                     r_sq_train_age=r_sq_train_age,
                                     MAE_train_relative_age=MAE_train_relative_age,
                                     RMSE_train_relative_age=RMSE_train_relative_age,
                                     r_sq_train_relative_age=r_sq_train_relative_age,
                                     MAE_test=MAE_test, RMSE_test=RMSE_test, r_sq_test=r_sq_test,
                                     alpha = best_alpha2, lambda = best_lambda2,total_par = total_par, non_zero_par = non_zero_par, stringsAsFactors=FALSE))

ggplot(NULL, aes(x=sample_data_train$Fstat, y=predicted_train_F, color = sample_data_train$Sex)) +
  geom_point() +  scale_colour_manual(values = c('M' = 'cyan3', 'F' = 'deeppink2'))+
  xlab("Transformed age") +
  ylab("Predicted Transformed") +
  geom_abline(intercept=0, slope=1, alpha = 0.5, color = 'green3') +
  geom_smooth(method = "lm", se = FALSE, color = "orchid1", alpha = 0.5)+
  # ylim(NA,3500) +
  annotate("text", x=Inf, y=Inf, size=2, vjust =1,hjust =1, label=paste0(" MAE = ",round(MAE_train_F,3)))+
  annotate("text", x=Inf, y=Inf, size=2, vjust =3, hjust =1,label=paste0(" R-squared = ",round(r_sq_train_F,3)))+
  annotate("text", x=Inf, y=Inf, size=2, vjust =5,hjust =1, label=paste0(" RMSE = ",round(RMSE_train_F,3)))+
  ggtitle(paste0("Training data ",prefix, '_cutoff',cutoff,'_surplus',surplus))+theme_classic()+coord_fixed() +guides(colour = guide_legend(title = "Sex"))
ggsave(paste0(out_dir,'/',prefix,'_cutoff',cutoff,'_surplus',surplus,"_elastic_net_train_plot_F.pdf"), height=7, width=10)


ggplot(NULL, aes(x=sample_data_train$Age, y=predicted_train_age, color = sample_data_train$Sex)) +
  geom_point() +  scale_colour_manual(values = c('M' = 'cyan3', 'F' = 'deeppink2'))+
  xlab("Chronological Age") +
  ylab("Predicted Age") +
  geom_abline(intercept=0, slope=1, alpha = 0.5, color = 'green3') +
  geom_smooth(method = "lm", se = FALSE, color = "orchid1", alpha = 0.5)+
  # ylim(NA,3500) +
  annotate("text", x=Inf, y=Inf, size=2, vjust =1,hjust =1, label=paste0(" MAE = ",round(MAE_train_age,3)))+
  annotate("text", x=Inf, y=Inf, size=2, vjust =3, hjust =1,label=paste0(" R-squared = ",round(r_sq_train_age,3)))+
  annotate("text", x=Inf, y=Inf, size=2, vjust =5,hjust =1, label=paste0(" RMSE = ",round(RMSE_train_age,3)))+
  ggtitle(paste0("Training data ",prefix, '_cutoff',cutoff,'_surplus',surplus))+theme_classic()+coord_fixed() +guides(colour = guide_legend(title = "Sex"))
ggsave(paste0(out_dir,'/',prefix,'_cutoff',cutoff,'_surplus',surplus,"_elastic_net_train_plot_Age.pdf"), height=7, width=10)


ggplot(NULL, aes(x=sample_data_train$relative_age, y=predicted_train_relative_age, color = sample_data_train$Sex)) +
  geom_point() +  scale_colour_manual(values = c('M' = 'cyan3', 'F' = 'deeppink2'))+
  xlab("Relative Age") +
  ylab("Predicted Relative Age") +
  geom_abline(intercept=0, slope=1, alpha = 0.5, color = 'green3') +
  geom_smooth(method = "lm", se = FALSE, color = "orchid1", alpha = 0.5)+
  # ylim(NA,3500) +
  annotate("text", x=Inf, y=Inf, size=2, vjust =1,hjust =1, label=paste0(" MAE = ",round(MAE_train_relative_age,3)))+
  annotate("text", x=Inf, y=Inf, size=2, vjust =3, hjust =1,label=paste0(" R-squared = ",round(r_sq_train_relative_age,3)))+
  annotate("text", x=Inf, y=Inf, size=2, vjust =5,hjust =1, label=paste0(" RMSE = ",round(RMSE_train_relative_age,3)))+
  ggtitle(paste0("Training data ",prefix,'_cutoff',cutoff,'_surplus',surplus))+theme_classic()+coord_fixed() +guides(colour = guide_legend(title = "Sex"))
ggsave(paste0(out_dir,'/',prefix,'_cutoff',cutoff,'_surplus',surplus,"_elastic_net_train_plot_relative_age.pdf"), height=7, width=10)


sample_data_test$Fstat <- young_transform (sample_data_test$Age,adult.age= cutoff, surplus = surplus) 
sample_data_test$relative_age <- relative_age (sample_data_test$Age,adult.age= cutoff, surplus = surplus) 


predicted_test_F <- model3 %>% predict(as.matrix(meth_data_test))
corr<-cor.test(sample_data_test$Fstat, predicted_test_F, alternative = c("two.sided", "less", "greater"),
               method = c("pearson", "kendall", "spearman"),
               exact = NULL, conf.level = 0.95, continuity = FALSE)

r_sq_test_F <- (corr[["estimate"]][["cor"]])^2
RMSE_test_F <- RMSE(sample_data_test$Fstat, predicted_test_F)
MAE_test_F <- MAE(sample_data_test$Fstat, predicted_test_F)

predicted_test_age <- anti_young_transform(f = predicted_test_F, adult.age= cutoff, surplus = surplus)
corr<-cor.test(sample_data_test$Age, predicted_test_age, alternative = c("two.sided", "less", "greater"),
               method = c("pearson", "kendall", "spearman"),
               exact = NULL, conf.level = 0.95, continuity = FALSE)

r_sq_test_age <- (corr[["estimate"]][["cor"]])^2
RMSE_test_age <- RMSE(sample_data_test$Age, predicted_test_age)
MAE_test_age <- MAE(sample_data_test$Age, predicted_test_age)

predicted_test_relative_age <- relative_age_f(f = predicted_test_F, adult.age= cutoff, surplus = surplus)
corr<-cor.test(sample_data_test$relative_age, predicted_test_relative_age, alternative = c("two.sided", "less", "greater"),
               method = c("pearson", "kendall", "spearman"),
               exact = NULL, conf.level = 0.95, continuity = FALSE)

r_sq_test_relative_age <- (corr[["estimate"]][["cor"]])^2
RMSE_test_relative_age <- RMSE(sample_data_test$relative_age, predicted_test_relative_age)
MAE_test_relative_age <- MAE(sample_data_test$relative_age, predicted_test_relative_age)

predicted_ages_df <- data.frame(predicted_test_age = predicted_test_age, predicted_test_F = predicted_test_F, predicted_test_relative_age = predicted_test_relative_age)
predicted_ages_df$actual_age <- sample_data_test$Age
predicted_ages_df$relative_age <- sample_data_test$relative_age
predicted_ages_df$Fstat <- sample_data_test$Fstat
predicted_ages_df$Sample <- sample_data_test$Sample.ID
predicted_ages_df$run <- prefix
write.csv(predicted_ages_df, file=paste0(out_dir,'/',prefix,'_cutoff',cutoff,'_surplus',surplus,"_predicted_ages_test.csv"), row.names=FALSE)


#save the results
saveRDS(results, file=paste0(out_dir,'/',prefix,'_cutoff',cutoff,'_surplus',surplus,"_results_table.RDS"))
write.csv(results, file=paste0(out_dir,'/',prefix,'_cutoff',cutoff,'_surplus',surplus,"_results_table.csv"), row.names=FALSE)

