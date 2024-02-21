# Verena Rubel 
# RPTU Kaiserslautern Landau
# 21.02.2024
# find indicators for dPCR assay using top 200 ASVs stored in tab_relab2
# first do tree approximation using error plot
# RF model using 1000 trees, mtry default, LOO x 10 CV
# indicators: ASVs which in no model show a MDA<0

# load required packages 
library(tidyr)
library(vegan)
library(tibble)
library(ggplot2)
library(dplyr)
library(caret)
library(randomForest)

# load env with info containing grouping (good/bad)
mapping <- read.csv("../data/Dunstaffnage_TS_metadata.csv", sep="\t") %>% filter(Seq.number!="X") %>%
  rename(Sample=Sample.code) %>% rename(EQ_class=Sample.name)
mapping_eq <- mapping %>% select(Sample, EQ_class)

# read target relative abundance matrix
tab_relab<- read.csv("tab_relab2.csv", row.names = 1)

# RF
# add label for rf prediction
rf_input <- tab_relab %>% as.data.frame() %>% rownames_to_column("Sample") %>%
  left_join(mapping_eq) %>% column_to_rownames("Sample")
rf_input$EQ_class = factor(rf_input$EQ_class) 

# errorplot to estimate required number of trees 
set.seed(42)
model_acc <- randomForest(EQ_class ~ ., data=rf_input, ntree=5000, importance = TRUE, proximity = TRUE)
saveRDS(model_acc, "model_acc.rds")
model_acc<-readRDS("model_acc.rds")

oob.error.data <- data.frame(
  Trees=rep(1:nrow(model_acc$err.rate), times=4),
  Type=rep(c("OOB", "CE", "AZE", "REF"), each=nrow(model_acc$err.rate)),
  Error=c(model_acc$err.rate[,"OOB"],
          model_acc$err.rate[,"CE"],
          model_acc$err.rate[,"AZE"],
          model_acc$err.rate[,"REF"]))
errorrate_acc <- ggplot(data=oob.error.data, aes(x=Trees, y=Error)) +
  geom_line(aes(color=Type))+
  scale_color_manual(values=c("orange", "brown1", "black", "forestgreen"))
errorrate_acc
ggsave("css_errorplot.pdf",errorrate_acc,height=5,width=6)


# Random forest leave one out approach
# classificatoin in CE/AZE/REF
# RANDOM FOREST CLASSIFICATION

### Random Forest LOO
train <- rf_input
uniqueIDs <- rownames(train)          
nruns <- length(uniqueIDs)    # number of cross validation runs: one for each unique ID           
crossclass <- match(rownames(train) , uniqueIDs)  
nobs <- nrow(na.omit(train))
crossPredict <- rep(NA, nobs)
choose_wdh <- 10 #number of repeats
result_rf_matrix <- matrix(nrow = choose_wdh, ncol = nruns)
colnames(result_rf_matrix) <- uniqueIDs
pb <- txtProgressBar(min = 0, max = choose_wdh, style = 3)

for (j in 1:choose_wdh) {
  crossPredict <- rep(NA, nobs)
  for (i in 1:nruns) {
    indtrain <- which(crossclass != i)
    indvalidate <- setdiff(1:nobs, indtrain)
    cat("Run", i, ": training only on observations with ID not", uniqueIDs[i], "\n")
    #IMPORTANT: set seed only for LOO without repetition!!!
    #set.seed(666)
    rf_df_CV <- randomForest(EQ_class ~ ., data = train[indtrain,],
                             ntree = 1000, na.action = "na.omit", importance=TRUE)
    crossPredict[indvalidate] <- predict(rf_df_CV, train[indvalidate,])
    result_rf_matrix[j,i] <- crossPredict[indvalidate]
    saveRDS(rf_df_CV, 
            paste0("/work/vdully/BA_time_series/RF_top200/models/model_without_", uniqueIDs[i], "_modelrepeat_", j, "X", ".rds"))
  }
  Sys.sleep(0.1)
  setTxtProgressBar(pb, j)
}

#end
close(pb)
result_rf_matrix

result_rf_matrix2 <- apply(result_rf_matrix, 2, function(x) ifelse(x == 1, "AZE", ifelse(x == 2, "CE", ifelse(x == 3, "REF", x))))
print(t(result_rf_matrix2))
write.csv(t(result_rf_matrix2), "/work/vdully/BA_time_series/RF_top200/predictions1.csv")
#result_rf_matrix2<- read.csv("/work/vdully/BA_time_series/RF_top200/predictions1.csv", row.names = 1) %>% t()
df <- result_rf_matrix2 %>% t() %>% as.data.frame()
df2 <- df %>% mutate("majority"=apply(df,1,function(x) names(which.max(table(x)))))

vgl_table <- df2 %>% select(majority) %>%
  rownames_to_column("Sample") %>%
  left_join(mapping_eq) %>%
  droplevels()
colnames(vgl_table) <- c("Sample", "pred", "truth")

conf_mat_final <- confusionMatrix(as.factor(vgl_table$truth), as.factor(vgl_table$pred))
print(conf_mat_final)
saveRDS(conf_mat_final, "/work/vdully/BA_time_series/RF_top200/confusion_matrix1.rds")

# read all models
# get all variable importances
# and combine


# MB NOR 
path <- "/work/vdully/BA_time_series/RF_top200/models/"
temp = list.files(path, pattern="*rds")
temp2 <- paste0(path, temp)
sams <- sapply(strsplit(basename(temp), ".rds"), `[`, 1)
sams2 <- sapply(strsplit(basename(sams), "without_"), `[`, 2)
for (i in 1:length(temp2)) assign(sams2[i], readRDS(temp2[i]))

data <- get(sams2[1])
model_0 <- importance(data) %>% as.data.frame() %>% rownames_to_column("ASV")

data2 <- get(sams2[2])
importance(get(sams2[2]))
model_1 <- importance(data2) %>% as.data.frame() %>% rownames_to_column("ASV")

model_comb <- model_0 %>% full_join(model_1, by="ASV")

sams3 <- paste0(sams2, "_VARIMP")
for (i in 1:length(sams2)) assign(sams3[i], importance(get(sams2[i])))

# construct empty matrix for each feature and run
result_matrix <- matrix(ncol= length(sams3), nrow=200)
colnames(result_matrix) <- sams3

anfang <- get(sams3[1])
result_matrix[,1] <- anfang[,3]
rownames(result_matrix) <- rownames(anfang)

for (i in 1:length(sams3)) {
  vals_tmp <- get(sams3[i])
  result_matrix[,i] <- vals_tmp[,3]
}

write.csv(result_matrix, "/work/vdully/BA_time_series/RF_top200/varimp_pro_model1.csv")

find_indicators <-  as.data.frame(result_matrix)
find_indicators$min <- do.call(pmin, find_indicators)
find_indicators_clean <- find_indicators %>% filter(min>0)
dim(find_indicators_clean)
write.csv(find_indicators_clean, "/work/vdully/BA_time_series/RF_top200/MIN0_indicators1.csv")

in_new_plot3 <- as.data.frame(result_matrix) %>% rownames_to_column("ASV") %>% gather("a", "b", -ASV) %>%
  filter(ASV %in% rownames(find_indicators_clean))

means_df <- in_new_plot3 %>% group_by(ASV) %>% summarise(mean=mean(b)) %>% arrange(-mean)  %>% top_n(25)
plot_varimps2 <- ggplot(means_df, aes(x=reorder(ASV,mean), y=mean)) + 
  geom_bar(stat="identity") +
  theme_light() +
  theme(axis.text.x = element_text(angle=0), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())+
  labs(y="Mean Variable Importance", x="")+
  coord_flip()
plot_varimps2
ggsave("/work/vdully/BA_time_series/RF_top200/varimp1.pdf", plot_varimps2, width=5, height =7)

