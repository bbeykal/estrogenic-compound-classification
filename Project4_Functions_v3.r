########################################
# Burcu Beykal
# May 31, 2020
########################################

# Read Name File
data_names <- read.csv("Reference_ER_Activity.csv", header = T)
feature_names <- read.csv("Feature_Names.csv", header = T)

####################################### DATA PREPROCESSING ###########################################
ERData_process <- function(data_set) {
	# Order Data Sets based on compound ID
	data_set <- data_set[order(data_set$compound, decreasing=F),]

	# Drop columns with channel 00
	data_set <- data_set[,-grep("ch00", colnames(data_set))]

	# Separate out Negative Control - Media
	media <- data_set[which(data_set$compound %in% "DMSO"),]
	data_set <- data_set[-which(data_set$compound %in% "DMSO"),]

	# Separate out Positive Control - E2 (Estradiol)
	e2 <- data_set[which(data_set$compound %in% "E2"),]
	data_set <- data_set[-which(data_set$compound %in% "E2"),]

	# Separate out 4HT
	ht4 <-  data_set[which(data_set$compound %in% "4HT"),]
	data_set <- data_set[-which(data_set$compound %in% "4HT"),]

	# Remove Inactive Compounds
	a02 <- data_set[which(data_set$compound %in% "A02"),]
	data_set <- data_set[-which(data_set$compound %in% "A02"),]

	a03 <- data_set[which(data_set$compound %in% "A03"),]
        data_set <- data_set[-which(data_set$compound %in% "A03"),]

	b06 <- data_set[which(data_set$compound %in% "B06"),]
        data_set <- data_set[-which(data_set$compound %in% "B06"),]

	c02 <- data_set[which(data_set$compound %in% "C02"),]
        data_set <- data_set[-which(data_set$compound %in% "C02"),]

	c05 <- data_set[which(data_set$compound %in% "C05"),]
        data_set <- data_set[-which(data_set$compound %in% "C05"),]

	c06 <- data_set[which(data_set$compound %in% "C06"),]
        data_set <- data_set[-which(data_set$compound %in% "C06"),]

	d03 <- data_set[which(data_set$compound %in% "D03"),]
        data_set <- data_set[-which(data_set$compound %in% "D03"),]

	d06 <- data_set[which(data_set$compound %in% "D06"),]
        data_set <- data_set[-which(data_set$compound %in% "D06"),]

	e03 <- data_set[which(data_set$compound %in% "E03"),]
        data_set <- data_set[-which(data_set$compound %in% "E03"),]

	f01 <- data_set[which(data_set$compound %in% "F01"),]
        data_set <- data_set[-which(data_set$compound %in% "F01"),]

	f02 <- data_set[which(data_set$compound %in% "F02"),]
        data_set <- data_set[-which(data_set$compound %in% "F02"),]

	h05 <- data_set[which(data_set$compound %in% "H05"),]
        data_set <- data_set[-which(data_set$compound %in% "H05"),]

	g02 <- data_set[which(data_set$compound %in% "G02"),]
        data_set <- data_set[-which(data_set$compound %in% "G02"),]



	# Merge masked compund IDs with correct names
	expand_names <- data.frame(matrix(NA, ncol = 2, nrow = dim(data_set)[1]))
	for (j in 1:length(data_names$Short_ID)) {
		for (i in 1:dim(data_set)[1]) {
			if (data_names$Short_ID[j] == data_set$compound[i]) {
				expand_names[i,1] <- toString(data_names$Preferred_Name[j])
				expand_names[i,2] <- toString(data_names$ER.Activity[j])
			}
		}
	}
	colnames(expand_names) <- c('Preferred.Name', 'ER.Activity')

	# Create a new data frame with experimental features only
	clean_set <- data_set[,-(1:6)]
	media_set <- media[,-(1:6)]
	e2_set 	  <- e2[,-(1:6)]

	for (k in 1:length(colnames(clean_set))) {
		if (colnames(clean_set)[k] == feature_names$Original.Feature.Names[k]) {
			colnames(clean_set)[k] <- as.character(feature_names$Even.Reduced.Names[k])
		}
	} 

	return(list(clean_set = clean_set, expand_names = expand_names, e2_set=e2_set, media_set=media_set))
}

# Average biological replicates
#avg_normal_control <- function(processed_data) {
#	shift <- as.numeric(rownames(unique(processed_data$expand_names)))
#	avg_normal_set <- data.frame(matrix(NA, nrow = dim(unique(processed_data$expand_names))[1], ncol = dim(processed_data$normal_control)[2]))
#	save_name <- NULL
#	for (i in 1:dim(avg_normal_set)[1]) {
 #       	for (j in 1:dim(avg_normal_set)[2]) {
  #              	avg_normal_set[i,j] <- mean(processed_data$normal_control[shift[i]:(shift[i]+7),j])
   #             	save_name[i] <- unique(processed_data$expand_names[shift[i]:(shift[i]+7),1])
   #     	}
#	}
#	rownames(avg_normal_set) <- save_name

#	return(avg_normal_set)
#}


# Data Normalization with respect to mean absolute deviation
normalize_mad <- function(clean_set, e2_set) {
        mad_e2 <- NULL
        bm_e2 <- data.frame(matrix(NA, nrow = dim(e2_set)[1], ncol = dim(clean_set)[2]))
        for (i in 1:dim(e2_set)[1]) {
                for (j in 1:dim(e2_set)[2]) {
                        bm_e2[i,j] <- abs(e2_set[i,j] - mean(e2_set[,j]))
                        mad_e2[j] <- mean(bm_e2[,j])
                }
        }

        # Built in function can also be used: careful with the definition of mad in documentation
        # normal_mad <- data.frame(matrix(NA, nrow = dim(clean_set)[1], ncol = dim(clean_set)[2]))
        # for (i in 1:dim(clean_set)[1]) {
        #        for (j in 1:dim(clean_set)[2]) {
        #               normal_mad[i,j] <- (clean_set[i,j] - median(e2_set[,j]))/mad(e2_set[,j], center = mean(e2_set[,j]), constant = 1)
        #        }
        # }

	df_mad <- data.frame(matrix(NA, nrow = dim(clean_set)[1], ncol = dim(clean_set)[2]))
        for (i in 1:dim(clean_set)[1]) {
                for (j in 1:dim(clean_set)[2]) {
                        df_mad[i,j] <- (clean_set[i,j] - median(e2_set[,j]))/mad_e2[j]
                }
        }
	colnames(df_mad) <- colnames(clean_set)
	return(df_mad)
}
 

# Data normalization with respect to controls
#normalize_control <- function(clean_set, e2_set, media_set){
#        df_control <- data.frame(matrix(NA, nrow = dim(clean_set)[1], ncol = dim(clean_set)[2]))
#        for (i in 1:dim(clean_set)[1]) {
#                for (j in 1:dim(clean_set)[2]) {
#                        df_control[i,j] <- (clean_set[i,j] - median(media_set[,j]))/(median(e2_set[,j]) - median(media_set[,j]))
#                }
#        }
#	colnames(df_control) <- colnames(clean_set)
#	return(df_control)
#}

# Regression fit function for CI calculation
bs <- function(formula, data, indices) {
	d <- data[indices,] # allows boot to select sample
	fit <- glm(formula, family=binomial(link="logit"), data=d)
	return(coef(fit))
}

# Model Performance Calculation
mperform <- function(fit_model, valid_data_numeric, valid_data_class, type ="logit") {
	if (type =="logit") {
 		predict_test   <- round(predict(fit_model,valid_data_numeric,type = "response"),3)
        	conf_valid     <- confusionMatrix(factor(round(predict_test)),valid_data_class$Class.Info, positive = "0", dnn = c("Predicted", "True"))
		Accuracy     <- round(as.numeric(conf_valid$overall[1]),2)
        	Sensitivity  <- round(as.numeric(conf_valid$byClass[1]),2)
        	Specificity  <- round(as.numeric(conf_valid$byClass[2]),2)
        	BA           <- round(as.numeric(conf_valid$byClass[11]),2)
        	CI_Lower     <- round(as.numeric(conf_valid$overall[3]),2)
        	CI_Upper     <- round(as.numeric(conf_valid$overall[4]),2)
	
	} else if (type == "rf") {
		predict_test <- predict(fit_model,valid_data_numeric)
                conf_valid   <- confusionMatrix(predict_test, valid_data_class$Class.Info, positive = "0", dnn = c("Predicted", "True"))
                Accuracy     <- round(as.numeric(conf_valid$overall[1]),2)
                Sensitivity  <- round(as.numeric(conf_valid$byClass[1]),2)
                Specificity  <- round(as.numeric(conf_valid$byClass[2]),2)
                BA           <- round(as.numeric(conf_valid$byClass[11]),2)
                CI_Lower     <- round(as.numeric(conf_valid$overall[3]),2)
                CI_Upper     <- round(as.numeric(conf_valid$overall[4]),2)
	} 
		summary.res  <- data.frame(Accuracy = Accuracy, CI.Lower = CI_Lower, CI.Upper = CI_Upper, Sensitivity = Sensitivity, Specificity = Specificity, Balanced.Accuracy = BA)
       return(summary.res)
}


