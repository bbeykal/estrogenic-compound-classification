########################################
# Burcu Beykal
# May 31, 2020
########################################

library(dendextend)
library(ggplot2)
library(ggdendro)
library(caret)
library(boot)
library(randomForest)
library(statip)
library(sm)
library(pROC)
library(stats)
library(reshape)

# Souce Function Files
source("./Project4_Functions_v3.r")

# Read Selected Experimental Data Set
data_sets <- NULL
for (j in 1:6) {
	data_sets[[j]] <- read.csv(paste("EPA_Reps_CSV/20190531A_Rep_",(j+3),"_Well_Data_updated_cmpnd.csv", sep=""), header = T) 
}

for (j in 7:12) {
        data_sets[[j]] <- read.csv(paste("EPA_Reps_CSV/20190531B_Rep_",(j+3),"_Well_Data_updated_cmpnd.csv", sep=""), header = T) 
}

for (j in 13:18) {
        data_sets[[j]] <- read.csv(paste("EPA_Reps_CSV/20190531C_Rep_",(j+3),"_Well_Data_updated_cmpnd.csv", sep=""), header = T)
}


############################## DATA PRE-PROCESSING ######################################################
# Call pre-processing function
processed_data     <- NULL
for (i in 1:length(data_sets)) {
        processed_data[[i]] <- ERData_process(data_sets[[i]])
}

# Replicate 8 (index =5) is selected for analysis 
new_frame <- data.frame(processed_data[[5]]$clean_set, Preferred.Name = processed_data[[5]]$expand_names$Preferred.Name)

# Average biological replicates raw biological replicates for outlier detection
aggregated_res <- aggregate(new_frame[,-dim(new_frame)[2]], by=list(new_frame$Preferred.Name), FUN = mean)
rownames(aggregated_res) <- aggregated_res$Group.1

aggregate_numeric <- aggregated_res[,-1]
# Cluster results and visualize it with a dendrogram
hc <- hclust(dist(aggregate_numeric, method = "euclidean"), method = "complete")
dd.row <- as.dendrogram(hc)
ddata_x <- dendro_data(dd.row,type ="rectangle")

# Vertical Dendrogram
ggplot(segment(ddata_x)) + geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
       geom_text(data=label(ddata_x),aes(label=label, x=x, y=0),hjust = 0, nudge_y = max(segment(ddata_x)$y) * 0.01, size = 5) +
       labs(x="",y="") +
       coord_flip() + 
	scale_y_reverse(expand=c(1,1))  +
        theme(legend.position="none",
	axis.text.y = element_blank(), 
	axis.line.y=element_blank(), 
	axis.ticks.y=element_blank(), 
	axis.ticks.x=element_blank(), 
	panel.background=element_rect(fill="white"),
	panel.grid=element_blank(),
	axis.text.x = element_blank(), 
	axis.title.x = element_blank())
	ggsave("Dendrogram_Outliers_P4.pdf", width = 10, height = 10, dpi = 300)

# Horizontal Dendrogram
#x11(width=12, height=10)
#ggplot(segment(ddata_x)) + geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
#        geom_text(data=label(ddata_x),aes(label=label, x=x, y=0), hjust = 1.03, angle=90, size = 4.5) +
#        labs(x="",y="") +
#        scale_y_continuous(expand=c(1,1)) +
#	theme(legend.position="none",
#        axis.text.y = element_blank(), 
#        axis.line.y=element_blank(), 
#        axis.ticks.y=element_blank(), 
#        axis.ticks.x=element_blank(), 
#        panel.background=element_rect(fill="white"),
#        panel.grid=element_blank(),
#        axis.text.x = element_blank(), 
#        axis.title.x = element_blank())

#    ggsave("Dendrogram_Outliers_P4.png", width = 12, height = 10, dpi = 300)



# Normalize data set
scale_to_mad <- NULL
clear_active <- NULL
for (i in 1:length(data_sets)) {
	clear_active[[i]] <- normalize_mad(processed_data[[i]]$clean_set, processed_data[[i]]$e2_set)
}

########################################### END OF DATA PRE-PROCESSING #############################################

########################################## FEATURE SELECTION ##################################################
# Cluster features using Pearson Correlation only for replicate 8 (index = 5)
similarity_pearson <- data.frame(matrix(NA, nrow = dim(clear_active[[5]])[2], ncol = dim(clear_active[[5]])[2]))
for(i in 1:dim(clear_active[[5]])[2]){
        for(j in 1:dim(clear_active[[5]])[2]){
            similarity_pearson[i,j] <- cor(as.numeric(as.vector(clear_active[[5]][,i])),as.numeric(as.vector(clear_active[[5]][,j])),method="pearson")
        }
}
rownames(similarity_pearson) <- colnames(similarity_pearson) <- colnames(clear_active[[5]])

hc1 <- hclust(as.dist(1-abs(similarity_pearson)), method = "complete")
dend <- as.dendrogram(hc1)
#png("./Dendo_Feat_Sel2.png",units='in',width = 10, height = 10, res = 600)
pdf("./Dendo_Feat_Sel2.pdf",width = 12, height = 15)
par(cex = 2, lwd=1.7, mar = c(4,0.05,0.05, 14) + 0.1,
    xpd = NA)

plot(dend, horiz = TRUE, xlab="Height")
abline(v = max(hc1$height)*0.05, col = 'red')
dev.off()

# Cut tree to 5% similarity
feat_sel <- cutree(hc1, h = max(hc1$height)*0.05)

# Select the topmost bilogically relevant features from independent clusters and reduce the data matrix size
col_select 	<- c("Array Area","Array PI Variance","Array Mean PI","Array Total PI","Array to Nucleoplasm Intensity Ratio")
reduced_data 	<- NULL
temp_class 	<- NULL
antagonists 	<- NULL

for (i in 1:length(data_sets)) {
	reduced_data[[i]] <- cbind(clear_active[[i]][,col_select],processed_data[[i]]$expand_names) # Append compound names and ER Activity
	temp_class[[i]]   <- data.frame(reduced_data[[i]], Class.Info = 0) # Append Class info 0 Agonist, +1 Antagonist
	antagonists[[i]] <- temp_class[[i]][which(temp_class[[i]]$ER.Activity %in% "Antagonist"),]
 	temp_class[[i]]$Class.Info[as.numeric(rownames(antagonists[[i]]))] <- 1
	temp_class[[i]]$Class.Info <- as.factor(temp_class[[i]]$Class.Info) # Make class info factor
}
class_data <- temp_class[[5]]

########################################### END OF FEATURE SELECTION #############################################


########################################### CLASSIFICATION MODEL TRAINING ##########################################
# Setting the seed for reproducibility
#RNGkind(sample.kind = "Rounding")
set.seed(1)

# Split Training Compounds
subset_data <- rbind(class_data[which(class_data$ER.Activity %in% "Antagonist"),],class_data[which(class_data$Preferred.Name %in% "Dicofol"),],class_data[which(class_data$Preferred.Name %in% "Diethylstilbestrol"),], class_data[which(class_data$Preferred.Name %in% "Estrone"),],class_data[which(class_data$Preferred.Name %in% "Fenarimol"),],class_data[which(class_data$Preferred.Name %in% "o,p'-DDT"),])

#Reduce Training Data to Numeric input-output
train_data <- data.frame(subset_data[,1:5], Class.Info = subset_data$Class.Info)

# Assign remaining compunds to test/validation set
test_data <- class_data[-as.numeric(rownames(subset_data)),]

test_data_numeric <- test_data[,1:5]
test_data_class	<- data.frame(test_data[,1:5], Class.Info = test_data$Class.Info)

########################## WILCOXON RANK SUM TEST ######################################
# Create combinations of 4 chemicals out of 5. This is done to get equal sizes of agonist and antagonist vectors for the wilcoxon test. There are 5 different combinations.
antagonist_group <- train_data[train_data$Class.Info == 1,] # all antagonists are selected
agonist_Maingroup <- train_data[train_data$Class.Info == 0,]

wilcox_p <- data.frame(pvalue = matrix(NA, nrow = 5, ncol = 1))
for (j in 1:(dim(train_data)[2]-1)) {
	wil_res <- wilcox.test(antagonist_group[,j], agonist_Maingroup[,j], paired=FALSE)
	wilcox_p[j,] <- wil_res$p.value
	rownames(wilcox_p)[j] <- colnames(agonist_Maingroup)[j] 
}

wilcox_p_round <- round(wilcox_p, 3)

write.csv(wilcox_p_round, file = "Wilcoxon_Test_pvalues.csv", row.names = TRUE)

######################### END WILCOXON RANK SUM TEST #######################################


########################## QUICK PCA ANALYSIS BEFORE MODELING #########################
# Show the PCA analysis and the distribution of training and test sets over 2 principal components
# For this part we will be using the processed unscaled data and we will center and scale during within the PCA function.

pca_data <- processed_data[[5]]$clean_set[,col_select]

# Call PCA function using prcomp
experiment.pca <- prcomp(pca_data, center = TRUE, scale = TRUE)

# Calculating percentage proportion of each PC 
percentage <- round(experiment.pca$sdev^2/sum(experiment.pca$sdev^2)*100, 2)
# Save results into a data frame with class information
df_out <- as.data.frame(experiment.pca$x)
df_out$group2 <- rep("Agonist", times = dim(df_out)[1])
df_out$group2[as.numeric(rownames(class_data[which(class_data$ER.Activity %in% "Antagonist"),]))] <- "Antagonist"

# For second type of visualization, mark the training and testing sets
df_out$group <- rep("Testing", times = dim(df_out)[1])

# Use training indices to remark the training rows
df_out$group[as.numeric(rownames(train_data))] <- "Training"

percentage <- paste(colnames(df_out),"(",paste(as.character(percentage),"% explained var.",")",sep=""))

# Visualizing PCA Analysis with biplots
ggplot(df_out, aes(x=PC1,y=PC2,color=group)) + geom_point() + stat_ellipse(level = 0.95) + 
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
	theme_bw() + xlab(percentage[1]) + ylab(percentage[2]) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position= c(0.1,0.1), legend.title = element_blank(), legend.background = element_blank()) + 
	scale_color_manual(values=c("#E69F00", "dodgerblue3")) + xlim(-7,5) + ylim(-3.5,4.5)

        ggsave("Train-test-PCA.pdf", width = 5.5, height = 4, dpi = 300)

ggplot(df_out, aes(x=PC1,y=PC2,color=group2)) + geom_point() + stat_ellipse(level = 0.95) + 
	geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
	theme_bw() + xlab(percentage[1]) + ylab(percentage[2]) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.12, 0.1), legend.title = element_blank(), legend.background = element_blank()) + 
        scale_color_manual(values=c("seagreen4", "mediumpurple3")) + xlim(-7,5) + ylim(-3.5, 4.5)

        ggsave("Agonist-Antagonist-PCA.pdf", width = 5.5, height = 4, dpi = 300)

######################## END PCA ANALYSIS ###################################


################## START LOGISTIC REGRESSION ##
# Fit Logistic Regression for Single Features (Total of 5 models)
cname_train <- colnames(train_data)[1:(length(colnames(train_data)) -1)]

# Initialization
model_param 	<- NULL
train_accuracy 	<- NULL
AIC_logit 	<- NULL
cv_error	<- NULL
cv_accuracy 	<- NULL
test_accuracy 	<- NULL
conf_interval 	<- NULL
confusion_test	<- NULL
pdf("./ROC_curve_training.pdf",height = 5, width = 6)
for (i in 1:length(cname_train)) {
	fmla <- paste("Class.Info ~", cname_train[i])
        fmla <- as.formula(fmla)
	fit_logit <- glm(fmla, family=binomial(link="logit"), data=train_data)

	# Calculate Training Accuracy
	predict_logit 	  <- predict(fit_logit,type = "response")
	confusion_train   <- confusionMatrix(factor(round(predict_logit)),train_data$Class.Info, positive = "0", dnn = c("Predicted", "True"))
        train_accuracy[i] <- round(as.numeric(confusion_train$overall[1]),2)
	AIC_logit[i] 	  <- round(AIC(fit_logit),2)
	model_param[[i]]  <- coef(fit_logit) 

	# Calculate 5-fold Cross Validation Error and Accuracy
	cv_error[i]	<- cv.glm(train_data, fit_logit, K =5)$delta[1]
        cv_accuracy[i]	<- round(1-cv_error[i],2)	

	# Calculate Testing Accuracy
        predict_test 	<- round(predict(fit_logit,test_data_numeric,type = "response"),3)
        confusion_test[[i]]  <- confusionMatrix(factor(round(predict_test)),test_data_class$Class.Info, positive = "0", dnn = c("Predicted", "True"))
        test_accuracy[i]   <- round(as.numeric(confusion_test[[i]]$overall[1]),2)
	
	# Calculate the Confidence Intervals via Bootstrapping
	out_boot <- boot(data=train_data, statistic=bs, R=400, formula=fmla)
	conf_interval[[i]] <- round(boot.ci(out_boot, type ="bca", index=1+1)$bca,2)

	# Plotting the ROC curves for the training data set
	color_list <- c("#1F78B4", "#33A02C", "#E31A1C", "#FF7F00", "#6A3D9A")
	if (i == 1) {
#		par(pty="s", mar=c(2,2,2,2))
		roc(train_data$Class.Info, predict_logit, plot = TRUE, legacy.axes = TRUE, xlab = "False Positive Rate (1 - Specificity)", ylab = "True Positive Rate (Sensitivity)", col = color_list[i], lwd = 4, print.auc = TRUE, asp = NA)
	} else {
		plot.roc(train_data$Class.Info, predict_logit, col = color_list[i], lwd = 4, print.auc = TRUE,  add = TRUE, print.auc.y = 0.5 - 0.05*(i-1), asp = NA)
	}
	legend("bottomright", legend=col_select, col = color_list, lwd = 2, cex = 0.75, box.lty=0)
}
dev.off()

# Bind results into 1 data frame and sort with respect to AIC Values
results_logit <- as.data.frame(cbind(cname_train, AIC_logit, cv_accuracy = cv_accuracy, test_accuracy = test_accuracy))
results_logit <- results_logit[order(test_accuracy, decreasing = TRUE),]


####### RESULTS VALIDATION ####
# Validate the best performing logistic regression model using other experimental replicates
validation_set <- temp_class
# Remove train/test replicate (index = 5)
validation_set[[5]] <- NULL

# Remove trained compounds from the set for unbiased estimate of prediction performance
subset_validation 	<- NULL
split_calculate 	<- NULL
valid_data_numeric <- NULL
valid_data_class <- NULL

for (i in 1:length(validation_set)) {
	subset_validation[[i]] <- rbind(validation_set[[i]][which(validation_set[[i]]$ER.Activity %in% "Antagonist"),],
					validation_set[[i]][which(validation_set[[i]]$Preferred.Name %in% "Dicofol"),],
					validation_set[[i]][which(validation_set[[i]]$Preferred.Name %in% "Diethylstilbestrol"),],
					validation_set[[i]][which(validation_set[[i]]$Preferred.Name %in% "Estrone"),],
					validation_set[[i]][which(validation_set[[i]]$Preferred.Name %in% "Fenarimol"),],
					validation_set[[i]][which(validation_set[[i]]$Preferred.Name %in% "o,p'-DDT"),])

	# Assign remaining compunds to test/validation set
	split_calculate[[i]] <- validation_set[[i]][-as.numeric(rownames(subset_validation[[i]])),]
	
	# Separate numerics and categorical outputs in validation data
	valid_data_numeric[[i]] <- split_calculate[[i]][,1:5]
	valid_data_class[[i]] <- data.frame(split_calculate[[i]][,1:5], Class.Info = split_calculate[[i]]$Class.Info)
}

# Validate also with full data
valid_full_numeric <- NULL
valid_full_class  <- NULL

for (i in 1:length(validation_set)) {
        # Separate numerics and categorical outputs in validation data
        valid_full_numeric[[i]] <- validation_set[[i]][,1:5]
        valid_full_class[[i]] <- data.frame(validation_set[[i]][,1:5], Class.Info = validation_set[[i]]$Class.Info)
}

# Call fit and predict functions
# Array PI Variance

fmla <- paste("Class.Info ~", cname_train[2])
fmla <- as.formula(fmla)
fit_logit_PI_Var <- glm(fmla, family=binomial(link="logit"), data=train_data)

temp_valid <- NULL
blind_PI_Var <- NULL
temp_valid2 <- NULL
full_PI_Var <- NULL
for (i in 1:length(validation_set)) {
	temp_valid <- mperform(fit_logit_PI_Var, valid_data_numeric[[i]], valid_data_class[[i]], type = "logit")
	blind_PI_Var <- rbind(blind_PI_Var,temp_valid)

	temp_valid2 <- mperform(fit_logit_PI_Var, valid_full_numeric[[i]], valid_full_class[[i]], type = "logit")
	full_PI_Var <- rbind(full_PI_Var,temp_valid2)
}

# Save results to csv file
write.csv(blind_PI_Var, file = "BlindValidation_Results_Array_PI_Variance_Logit.csv", row.names = FALSE)
write.csv(full_PI_Var, file = "FullValidation_Results_Array_PI_Variance_Logit.csv", row.names = FALSE)


# Array to Nucleoplasm Intensity Ratio
fmla <- paste("Class.Info ~", cname_train[5])
fmla <- as.formula(fmla)
fit_logit_NIR <- glm(fmla, family=binomial(link="logit"), data=train_data)

temp_valid <- NULL
blind_NIR <- NULL
temp_valid2 <- NULL
full_NIR <- NULL

for (i in 1:length(validation_set)) {
        temp_valid <- mperform(fit_logit_NIR, valid_data_numeric[[i]], valid_data_class[[i]], type = "logit")
        blind_NIR <- rbind(blind_NIR,temp_valid)

        temp_valid2 <- mperform(fit_logit_NIR, valid_full_numeric[[i]], valid_full_class[[i]], type = "logit")
        full_NIR <- rbind(full_NIR,temp_valid2)
}

# Save results to csv file
write.csv(blind_NIR, file = "BlindValidation_Results_NIR_Logit.csv", row.names = FALSE)
write.csv(full_NIR, file = "FullValidation_Results_NIR_Logit.csv", row.names = FALSE)

######### END OF VALIDATION ####


###################### END LOGISTIC REGRESSION ##

###################### START RANDOM FOREST ######
set.seed(123)
#set.seed(44)
fit_rf <- randomForest(formula = Class.Info ~ ., data=train_data, importance = TRUE, ntree = 130)

# Calculate Training Accuracy
predict_rf  	<- predict(fit_rf)
confusion_RF_train <- confusionMatrix(predict_rf,train_data$Class.Info, positive = "0", dnn = c("Predicted", "True"))
train_RF_accuracy  <- round(as.numeric(confusion_RF_train$overall[1]),2)

# Calculate Testing Accuracy
predict_rf_test    <- predict(fit_rf,test_data_numeric)
confusion_RF_test  <- confusionMatrix(predict_rf_test,test_data_class$Class.Info, positive = "0", dnn = c("Predicted", "True"))
test_RF_accuracy   <- round(as.numeric(confusion_RF_test$overall[1]),2)

# Calculate Mean Decrease in Gini Index
feature_imp <- as.data.frame(fit_rf$importance)
feature_rank <- as.matrix(feature_imp[order(feature_imp$MeanDecreaseGini, decreasing = T), 4])
rownames(feature_rank) <- rownames(feature_imp)[order(feature_imp$MeanDecreaseGini, decreasing = T)]
colnames(feature_rank) <- c("Mean Decrease in Gini Index")
feature_rank


write.csv(feature_rank, file = "Feauture_Ranking_RF.csv", row.names = TRUE)


# Calculate Validation Accuracies

temp_valid <- NULL
blind_RF <- NULL
temp_valid2 <- NULL
full_RF <- NULL

for (i in 1:length(validation_set)) {
        temp_valid <- mperform(fit_rf, valid_data_numeric[[i]], valid_data_class[[i]], type = "rf")
        blind_RF <- rbind(blind_RF,temp_valid)

        temp_valid2 <- mperform(fit_rf, valid_full_numeric[[i]], valid_full_class[[i]], type ="rf")
        full_RF <- rbind(full_RF,temp_valid2)
}


# Save results to csv file
write.csv(blind_RF, file = "BlindValidation_Results_RF.csv", row.names = FALSE)
write.csv(full_RF, file = "FullValidation_Results_RF.csv", row.names = FALSE)

###################### END RANDOM FOREST ##

##################### BEGIN FEATURE DENSITY VISUALIZATION ########################################
agonists <- NULL
for (i in 1:length(temp_class)) {
	agonists[[i]] <- temp_class[[i]][which(temp_class[[i]]$Class.Info == 0),]
}

# Calculate Hellinger Distance between Agonists and Antagonists for Array PI Variance

hdist_PI_Var <-NULL
pdf("./Histogram_PI_Var.pdf", width = 8, height = 16)
par(mfrow=c(6,3), cex.axis = 2.2, cex.lab = 2.2, lwd =2, mar = c(4.2, 4.8, 2,1.5)+ 0.1)

for (i in 1:length(agonists)) {
	hdist_PI_Var[[i]] <- round(hellinger(agonists[[i]]$Array.PI.Variance, antagonists[[i]]$Array.PI.Variance,-Inf, Inf, method = 1),2)
        sm.density.compare(temp_class[[i]]$Array.PI.Variance, temp_class[[i]]$Class.Info, xlab=paste("Exp",i,", HD=", hdist_PI_Var[[i]], sep=""), lty = c(1, 1),
             lwd = c(2, 2), col = c("dodgerblue3","red3"))
}

dev.off()

# Calculate Hellinger Distance between Agonists and Antagonists for Array to Nucleoplasm Intensity Ratio

hdist_NIR <-NULL
pdf("./Histogram_NIR.pdf", width = 8, height = 16)
par(mfrow=c(6,3), cex.axis = 2.2, cex.lab = 2.2, lwd =2, mar = c(4.2, 4.8, 2,1.5)+ 0.1)

for (i in 1:length(agonists)) {
        hdist_NIR[[i]] <- round(hellinger(agonists[[i]]$Array.to.Nucleoplasm.Intensity.Ratio, antagonists[[i]]$Array.to.Nucleoplasm.Intensity.Ratio,-Inf, Inf, method = 1),2)
	sm.density.compare(temp_class[[i]]$Array.to.Nucleoplasm.Intensity.Ratio, temp_class[[i]]$Class.Info, xlab=paste("Exp",i,", HD=", hdist_NIR[[i]],sep=""), lty = c(1, 1),
             lwd = c(2, 2), col = c("dodgerblue3","red3"))
}
dev.off()

##### BEGIN BOXPLOTS FOR PERFORMANCE METRICS ############
# Append model names to data frames
full_NIR$Model <- rep("Logistic Regression 'Array to Nucleoplasm Intensity Ratio'",times = dim(full_NIR)[1])  
full_PI_Var$Model <- rep("Logistic Regression 'Array PI Variance'", times=dim(full_PI_Var)[1])
full_RF$Model <- rep("Random Forest", times = dim(full_RF)[1])

# Merge all results
merged_data <- rbind(full_PI_Var,full_NIR,full_RF)

#Drop confidence intervals
merged_data$CI.Lower <- NULL
merged_data$CI.Upper <- NULL 
colnames(merged_data)[4] <- "Balanced Accuracy"
# Melt data for plotting with respect to model type
melted_data <- melt(merged_data, id.vars = "Model")

ggplot(melted_data, aes(x=variable, y =value, fill=Model)) + geom_boxplot() + theme_bw() + 
theme(panel.grid.minor = element_blank(), 
	axis.title.x=element_blank(),
	axis.text=element_text(size=10),
	legend.justification=c("left","bottom"),
	legend.position = c(.01, .05),
	legend.text=element_text(size=8),
	legend.title = element_blank()) + 
scale_fill_brewer(palette ="Blues") + 
labs(y="Model Performance")

ggsave("Boxplot_all_active_validation.pdf", width = 6.5, height = 3.16, dpi = 600)


# Removing high-density experimental replicates from the analysis
full_NIR_red <- full_NIR[-c(2,4,5,7,9,11,12,13,15,17),]
full_PI_Var_red <- full_PI_Var[-c(2,4,5,7,9,11,12,13,15,17),]
full_RF_red <- full_RF[-c(2,4,5,7,9,11,12,13,15,17),]

# Merge all results
merged_data2 <- rbind(full_PI_Var_red,full_NIR_red,full_RF_red)

#Drop confidence intervals
merged_data2$CI.Lower <- NULL
merged_data2$CI.Upper <- NULL 
colnames(merged_data2)[4] <- "Balanced Accuracy"
# Melt data for plotting with respect to model type
melted_data2 <- melt(merged_data2, id.vars = "Model")

ggplot(melted_data2, aes(x=variable, y =value, fill=Model)) + geom_boxplot() + theme_bw() + 
theme(panel.grid.minor = element_blank(), 
	axis.title.x=element_blank(),
	axis.text=element_text(size=10),
	legend.justification=c("left","bottom"),
	legend.position = c(.01, .05),
	legend.text=element_text(size=8),
	legend.title = element_blank()) +
scale_fill_brewer(palette ="Reds") + ylim(0,1) +
labs(y="Model Performance")

ggsave("Boxplot_all_active_after_data_quality_check.pdf", width = 6.5, height = 3.16, dpi = 600)
