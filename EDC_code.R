###########################################
##SNP聚类##
###########################################

rm(list=ls())
setwd("E:/")
path="E:/"

data_input <- read.csv(paste(path, "EDC_1879_96_num.csv", sep=""), header=T)
data_snp <- data_input[,2:1880]
rownames(data_snp) <- data_input[,1]

# 统计缺失
library(VIM)
miss <- aggr(data_snp, prop=T, numbers=T)
write.table(miss[["missings"]], "missing_sample_snp.txt", sep="\t", row.names=F, quote=F)

# 去掉缺失10%以上的样本
data_snp <- data_snp[,c(-50,-128,-177,-186,-203,-242,-255,-262,-278,-281,-292,-364,-381,-382,-408,-449,-468,-470,-471,-472,-473,-482,-858,-859,-860,-873,-874,-875,-876,-877,-890,-891,-897,-907,-908,-918,-941,-942,-943,-951,-976,-977,-990,-1042,-1043,-1058,-1068,-1069,-1106,-1128,-1131,-1132,-1133,-1151,-1166,-1176,-1195,-1204,-1206,-1233,-1236,-1237,-1238,-1299,-1300,-1302,-1303,-1304,-1334,-1335,-1353,-1354,-1355,-1356,-1360,-1362,-1363,-1364,-1365,-1366,-1368,-1369,-1370,-1371,-1372,-1373,-1374,-1375,-1393,-1398,-1399,-1400,-1416,-1419,-1420,-1421,-1427,-1487,-1488,-1579,-1659,-1696,-1781,-1784,-1787,-1791)]

# Identifying Zero-Variance Predictors
library(caret)
zv <- preProcess(t(data_snp), method=c("zv"))
zv[["method"]][["remove"]]

# 插值众数
MS <- function(x){return(as.numeric(names(table(x))[table(x)==max(table(x))]))}
data_snp <- as.data.frame(t(data_snp))
n=ncol(data_snp)
for(i in 1:n){
  zs <- MS(data_snp[,i])
  data_snp[is.na(data_snp[,i]),i] <- zs
  print(i)
}
data_snp <- as.data.frame(t(data_snp))
data_snp[,] <- lapply(data_snp[,], as.factor)

# 随机森林插补缺失值（可以处理混合型数据）
library(missForest)
data_snp <- as.data.frame(t(data_snp))
data_snp[,] <- lapply(data_snp[,], as.factor)
set.seed(101)
# process <- missForest(data_snp, verbose=T)
process <- missForest(data_snp, maxiter=2, ntree=1000, mtry=44, verbose=T)
data_snp <- process$ximp
data_snp <- as.data.frame(t(data_snp))

# mtry selection
mtry_used <- c()
OOB_NRMSE <- c()
OOB_PFC <- c()
i <- 1
for (loop_seed in 1:5){
  for(mtry in 1:96){
    print(paste0("seed = ", loop_seed, ", mtry = ", mtry))
    set.seed(loop_seed)
    process <- missForest(data_snp, maxiter=10, ntree=100, mtry=mtry, verbose=T)
    OOB_NRMSE[i] <- process$OOBerror[2]
    OOB_PFC[i] <- process$OOBerror[1]
    mtry_used[i] <- mtry
    i <- i+1
  }
}
library(dplyr)
library(tidyr)
mtry_df <- data.frame(mtry_used, OOB_NRMSE, OOB_PFC) %>% 
  group_by(mtry_used) %>% 
  summarize(`OOB NRMSE`=mean(OOB_NRMSE), `OOB PFC`=mean(OOB_PFC)) %>% 
  gather("metric", "error", -mtry_used)
ggplot(mtry_df, aes(x=mtry_used, y=error, col=factor(metric))) + 
  geom_point() + 
  geom_line() + 
  scale_x_continuous(breaks=seq(0, 96, 2)) + 
  scale_color_manual(values=c("deepskyblue", "mediumseagreen")) + 
  facet_wrap(~ metric, scales="free_y") + 
  theme(legend.position="none", 
        axis.title.y=element_blank()) + 
  labs(x="mtry")

# 计算矩阵相关性
data_snp[,] <- lapply(data_snp[,], as.numeric)
data_snp[,] <- data_snp[,]-1
diag(cor(data_snp, data_snp2))

# km聚类
set.seed(101)
km <- kmeans(data_snp[,1:1773], 2)
data_snp <- cbind(data_snp, km_cluster=km$cluster)

# km_cluster颜色
library(ggsci)
library(scales)
colors <- pal_igv("default")(2)
names(colors) <- unique(data_snp$km_cluster)

# t-SNE 2D
library(Rtsne)
library(pca3d)
perp <- (nrow(data_snp)-1)/3
set.seed(101)
tsne <- Rtsne(data_snp[,1:1773], dims=2, perplexity=perp, check_duplicates=FALSE)
plot(tsne$Y, t='n', main="SNP")
text(tsne$Y, labels=data_snp$km_cluster, col=colors[data_snp$km_cluster])
tsneY <- tsne[["Y"]]
pca2d(tsneY, components=1:2, group=data_snp$km_cluster, show.ellipses=TRUE)

# 基于tsne的km聚类(相当于dbscan聚类)
set.seed(101)
ds <- kmeans(tsne$Y, 2)
data_snp <- cbind(data_snp, ds_cluster=ds$cluster)

# data_snp <- read.table(paste(path, "EDC_1773_96_num_km_ds.txt", sep=""), header=TRUE, row.names=1)

# ds_cluster颜色
library(ggsci)
library(scales)
colors <- pal_igv("default")(2)
names(colors) <- unique(data_snp$ds_cluster)

# t-SNE 2D
library(Rtsne)
library(pca3d)
perp <- (nrow(data_snp)-1)/3
set.seed(101)
tsne <- Rtsne(data_snp[,1:1773], dims=2, perplexity=perp, check_duplicates=FALSE)
plot(tsne$Y, t='n', main="SNP")
text(tsne$Y, labels=data_snp$ds_cluster, col=colors[data_snp$ds_cluster])
tsneY <- tsne[["Y"]]
pca2d(tsneY, components=1:2, group=data_snp$ds_cluster, show.ellipses=TRUE)
# or
library(Rtsne)
perp <- (nrow(data_snp)-1)/3
set.seed(101)
tsne <- Rtsne(data_snp[,1:1773], dims=2, perplexity=27, check_duplicates=FALSE)
pch <- c(19,17)
names(pch) <- unique(data_snp$ds_cluster)
plot(tsne$Y, type="p", main="SNP", xlab="", ylab="", cex.axis=1.2, pch=pch[data_snp$ds_cluster], lwd=2, col=colors[data_snp$ds_cluster], axes=F)
legend(-8, 4, c("cluster1","cluster2"), cex=1.2, x.intersp=0.5, y.intersp=0.5, pch=c(19,17), pt.lwd=2, col=pal_igv("default")(2), box.col="white")
axis(side=1, lwd=2, font=2)
axis(side=2, lwd=2, font=2)

data_snp <- cbind(SNP=rownames(data_snp), data_snp)
write.table(data_snp, "EDC_1773_96_num_km_ds.txt", sep="\t", row.names=F, quote=F)

###########################################
##SNP筛选##
###########################################

rm(list=ls())

setwd("E:/")
path="E:/"

data_input <- read.table(paste(path, "EDC_1773_96_num_km_ds.txt", sep=""), header=TRUE, row.names=1)
data_snp <- as.data.frame(t(data_input[,1:1773]))
data_y <- read.csv(paste(path, "EDC_1879_133_num.csv", sep=""), header=TRUE)
data_y <- data_y[c(-50,-128,-177,-186,-203,-242,-255,-262,-278,-281,-292,-364,-381,-382,-408,-449,-468,-470,-471,-472,-473,-482,-858,-859,-860,-873,-874,-875,-876,-877,-890,-891,-897,-907,-908,-918,-941,-942,-943,-951,-976,-977,-990,-1042,-1043,-1058,-1068,-1069,-1106,-1128,-1131,-1132,-1133,-1151,-1166,-1176,-1195,-1204,-1206,-1233,-1236,-1237,-1238,-1299,-1300,-1302,-1303,-1304,-1334,-1335,-1353,-1354,-1355,-1356,-1360,-1362,-1363,-1364,-1365,-1366,-1368,-1369,-1370,-1371,-1372,-1373,-1374,-1375,-1393,-1398,-1399,-1400,-1416,-1419,-1420,-1421,-1427,-1487,-1488,-1579,-1659,-1696,-1781,-1784,-1787,-1791),c(98,99)]
data_snp <- cbind(data_snp, data_y)

data_snp <- data_snp[is.na(data_snp$是否发生前哨淋巴结转移)==FALSE & is.na(data_snp$是否发生非前哨淋巴结转移)==FALSE,]
data_snp <- data_snp[data_snp$是否发生前哨淋巴结转移==1,]
data_snp <- data_snp[,-97]

# Train and test sets
set.seed(3)
index_test_0 <- sample(which(data_snp$是否发生非前哨淋巴结转移==0), 30)
set.seed(3)
index_test_1 <- sample(which(data_snp$是否发生非前哨淋巴结转移==1), 16)
index_test <- c(index_test_0,index_test_1)
data_train <- data_snp[-index_test,]
data_test <- data_snp[index_test,]
# write.table(rownames(data_train), "data_train.txt", sep="\t", row.names=F, col.names=F, quote=F)
# write.table(rownames(data_test), "data_test.txt", sep="\t", row.names=F, col.names=F, quote=F)

# 查看比例
table(data_snp$是否发生非前哨淋巴结转移)
prop.table(table(data_snp$是否发生非前哨淋巴结转移))
table(data_train$是否发生非前哨淋巴结转移)
prop.table(table(data_train$是否发生非前哨淋巴结转移))
table(data_test$是否发生非前哨淋巴结转移)
prop.table(table(data_test$是否发生非前哨淋巴结转移))

# Identifying Zero-Variance Predictors
library(caret)
zv <- preProcess(data_train, method=c("zv"))
zv[["method"]][["remove"]]
idx <- !(colnames(data_train) %in% zv[["method"]][["remove"]])
data_train <- data_train[,idx]

# wilcox.test
library(data.table)
data_all <- data_train
data_all_wilcox <- data.table(mut=colnames(data_all[,1:(ncol(data_all)-1)]), p.value=numeric(ncol(data_all)-1))
for(i in 1:(ncol(data_all)-1)){
  p.value <- wilcox.test(data_all[,i]~是否发生非前哨淋巴结转移, data=data_all, alternative='two.sided', conf.level=0.95)[["p.value"]]
  data_all_wilcox[mut==colnames(data_all)[i], "p.value"] <- p.value
}
rm(p.value, i)
write.table(data_all_wilcox, "snp_wilcox_3.txt", sep="\t", row.names=F, quote=F)

# feature selection(svm&wilcox)
library(ROSE)
library(e1071)
data_all <- data_train
result_pvals <- matrix(ncol=3)
for(k in 99:101){
  set.seed(k)
  folds <- createFolds(data_all[,"是否发生非前哨淋巴结转移"], k=5)
  for(i in (ncol(data_all)-1):1){
    data_cv <- data_all[,c(order(data_all_wilcox$p.value)[1:i],ncol(data_train))]
    for(j in 1:5){
      fold_test <- data_cv[folds[[j]],]
      fold_train <- data_cv[-folds[[j]],]
      svm_train <- svm(是否发生非前哨淋巴结转移~., data=fold_train, kernel="sigmoid") # kernel=c("linear","radial"(default),"sigmoid")
      svm_test <- predict(svm_train, newdata=fold_test)
      auc <- roc.curve(fold_test$是否发生非前哨淋巴结转移, svm_test)$auc
      line <- c(i, j, auc)
      result_pvals <- rbind(result_pvals, line)
    }
  }
}
result_pvals <- as.data.frame(result_pvals[-1,])
colnames(result_pvals) <- c("Features","Fold","AUC")
write.table(result_pvals, "pvals_auc_non-sentinel_sigmoid.txt", sep="\t", row.names=F, quote=F)
ggplot(result_pvals, aes(x=as.factor(Features), y=AUC)) + 
  geom_boxplot(colour="#0073C2FF") + 
  scale_y_continuous(limits=c(0.45,0.95), breaks=seq(0.45,0.95,0.05)) + 
  # stat_summary(fun=mean, shape=19, col='red', geom='point') + 
  xlab("Number of features") + 
  ylab("AUC (cross-validation)") + 
  theme(
    axis.text.x=element_text(size=rel(0.8),colour="black"),
    axis.text.y=element_text(size=rel(1.2),colour="black"),
    axis.title.x=element_text(size=rel(1.2),colour="black"),
    axis.title.y=element_text(size=rel(1.2),colour="black"),
    axis.line.x=element_line(size=0.5,colour="black"),
    axis.line.y=element_line(size=0.5,colour="black"),
    panel.border=element_blank(),
    panel.grid=element_blank(),
    panel.background=element_rect(fill="transparent",color=NA)
  )

###########################################
##PRS计算##
###########################################

rm(list=ls())

setwd("E:/")
path="E:/"

library(data.table)
data_input <- read.csv(paste(path, "EDC_1879_133_num.csv", sep=""), header=TRUE)
snp_cluster <- fread(paste(path, "EDC_1773_96_num_km_ds.txt", sep=""), sep="\t")
data_snp <- cbind(as.data.frame(t(snp_cluster[,2:1774])),data_input[c(-50,-128,-177,-186,-203,-242,-255,-262,-278,-281,-292,-364,-381,-382,-408,-449,-468,-470,-471,-472,-473,-482,-858,-859,-860,-873,-874,-875,-876,-877,-890,-891,-897,-907,-908,-918,-941,-942,-943,-951,-976,-977,-990,-1042,-1043,-1058,-1068,-1069,-1106,-1128,-1131,-1132,-1133,-1151,-1166,-1176,-1195,-1204,-1206,-1233,-1236,-1237,-1238,-1299,-1300,-1302,-1303,-1304,-1334,-1335,-1353,-1354,-1355,-1356,-1360,-1362,-1363,-1364,-1365,-1366,-1368,-1369,-1370,-1371,-1372,-1373,-1374,-1375,-1393,-1398,-1399,-1400,-1416,-1419,-1420,-1421,-1427,-1487,-1488,-1579,-1659,-1696,-1781,-1784,-1787,-1791),c(98,99)])
colnames(data_snp) <- c(t(snp_cluster[,1]),"是否发生前哨淋巴结转移","是否发生非前哨淋巴结转移")
rownames(data_snp) <- data_input[c(-50,-128,-177,-186,-203,-242,-255,-262,-278,-281,-292,-364,-381,-382,-408,-449,-468,-470,-471,-472,-473,-482,-858,-859,-860,-873,-874,-875,-876,-877,-890,-891,-897,-907,-908,-918,-941,-942,-943,-951,-976,-977,-990,-1042,-1043,-1058,-1068,-1069,-1106,-1128,-1131,-1132,-1133,-1151,-1166,-1176,-1195,-1204,-1206,-1233,-1236,-1237,-1238,-1299,-1300,-1302,-1303,-1304,-1334,-1335,-1353,-1354,-1355,-1356,-1360,-1362,-1363,-1364,-1365,-1366,-1368,-1369,-1370,-1371,-1372,-1373,-1374,-1375,-1393,-1398,-1399,-1400,-1416,-1419,-1420,-1421,-1427,-1487,-1488,-1579,-1659,-1696,-1781,-1784,-1787,-1791),1]

data_snp <- data_snp[is.na(data_snp$是否发生前哨淋巴结转移)==FALSE & is.na(data_snp$是否发生非前哨淋巴结转移)==FALSE,]
data_snp <- data_snp[data_snp$是否发生前哨淋巴结转移==1,]
data_snp <- data_snp[,c(-97)]

# Train and test sets
set.seed(3)
index_test_0 <- sample(which(data_snp$是否发生非前哨淋巴结转移==0), 30)
set.seed(3)
index_test_1 <- sample(which(data_snp$是否发生非前哨淋巴结转移==1), 16)
index_test <- c(index_test_0,index_test_1)
data_train <- data_snp[-index_test,]
data_test <- data_snp[index_test,]

# coef计算
data_snp_all <- data_input[,c(2:97,98,99)]
rownames(data_snp_all) <- data_input[,1]
data_snp_all <- data_snp_all[is.na(data_snp_all$是否发生非前哨淋巴结转移)==FALSE,]
# data_snp_all <- data_snp_all[is.na(data_snp_all$是否发生前哨淋巴结转移)==FALSE,]
# data_snp_all <- data_snp_all[is.na(data_snp_all$是否发生非前哨淋巴结转移)==FALSE | data_snp_all$是否发生前哨淋巴结转移==0,]
# data_snp_all$是否发生非前哨淋巴结转移[is.na(data_snp_all$是否发生非前哨淋巴结转移)] <- 0
# idx <- !(rownames(data_snp_all) %in% rownames(data_test))
# data_snp_all <- data_snp_all[idx,]
data_snp_all <- data_snp_all[,-97]

# 单因素二分类逻辑回归(binomial)
library(data.table)
data_logit_all <- data.table(mut=colnames(data_snp_all)[-ncol(data_snp_all)], coe=numeric(ncol(data_snp_all)-1))
for(i in 1:(ncol(data_snp_all)-1)){
  data_logit <- data_snp_all[is.na(data_snp_all[,i])==FALSE,]
  model_logit <- glm(是否发生非前哨淋巴结转移 ~ data_logit[,i], family=binomial(link="logit"), data=data_logit, control=list(maxit=100))
  coe <- coef(model_logit)[2]
  data_logit_all[mut==colnames(data_snp_all)[i], c("coe")] <- c(coe)
}
rm(data_logit, model_logit, coe, i)
# write.table(data_logit_all, "snp_coef_non-sentienel_392.txt", sep="\t", row.names=F, quote=F)

# Identifying Zero-Variance Predictors
library(caret)
data_train <- data_train[,c(-97)]
data_snp <- data_snp[,c(-97)]
zv <- preProcess(data_train, method=c("zv"))
zv[["method"]][["remove"]]
idx <- !(colnames(data_train) %in% zv[["method"]][["remove"]])
data_snp <- data_snp[,idx]

# all SNP
snp_cluster <- snp_cluster[idx,]
data_logit <- data_logit_all[idx,]
data_logit$ds_cluster <- snp_cluster$ds_cluster

# selected SNP
snp_cluster <- snp_cluster[idx,]
data_logit <- data_logit_all[idx,]
data_logit$ds_cluster <- snp_cluster$ds_cluster
snp_order <- fread(paste(path, "snp_wilcox_3.txt", sep=""), sep="\t")
idx_sig1 <- colnames(data_snp) %in% sapply(snp_order[order(snp_order$p.value)[1:29],1], as.character)
data_snp <- data_snp[,idx_sig1]
idx_sig2 <- sapply(data_logit[,1], as.character) %in% sapply(snp_order[order(snp_order$p.value)[1:29],1], as.character)
data_logit <- data_logit[idx_sig2,]

# cluster PRS
for(i in 1:nrow(data_snp)){
  cluster1 <- 0
  cluster2 <- 0
  for(j in 1:ncol(data_snp)){
    if(data_logit$ds_cluster[j]==1 & data_snp[i,j]==1){cluster1 <- cluster1+data_logit$coe[j]}
    if(data_logit$ds_cluster[j]==2 & data_snp[i,j]==1){cluster2 <- cluster2+data_logit$coe[j]}
  }
  data_snp$cluster1[i] <- 1/(1+exp(-cluster1))
  data_snp$cluster2[i] <- 1/(1+exp(-cluster2))
}
data_snp <- cbind(序号=rownames(data_snp), data_snp)
write.table(data_snp, "snp_cluster_non-sentinel_229_3_29.txt", sep="\t", row.names=F, quote=F)
data_snp <- data_snp[,-1]

###########################################
##特征选择##
###########################################

rm(list=ls())

setwd("E:/")
path="E:/"

data_input <- read.csv(paste(path, "EDC_229_35_num_non-sentinel_3_84.csv", sep=""), header=TRUE)
data <- data_input[,c(-1,-21,-23,-34)]

# Train and test sets
set.seed(3)
index_test_0 <- sample(which(data$是否发生非前哨淋巴结转移==0), 30)
set.seed(3)
index_test_1 <- sample(which(data$是否发生非前哨淋巴结转移==1), 16)
index_test <- c(index_test_0,index_test_1)
data_train <- data[-index_test,]
data_test <- data[index_test,]

# Identifying Zero- and Near Zero-Variance Predictors
library(caret)
nzv <- nearZeroVar(data_train, saveMetrics=T)
nzv[nzv$nzv,]
nzv_index <- nearZeroVar(data_train)
data_train <- data_train[,-nzv_index]

# Identifying Correlated Predictors
library(caret)
data.cor <- cor(data_train[,1:(ncol(data_train)-1)], method="spearman")
summary(data.cor[upper.tri(data.cor)])
high.cor <- findCorrelation(data.cor, cutoff=0.90)
data_train <- data_train[,-high.cor]

# Normalization
library(caret)
standard <- preProcess(data_train[,1:(ncol(data_train)-1)], method=c("center","scale"))
data_train[,1:(ncol(data_train)-1)] <- predict(standard, data_train[,1:(ncol(data_train)-1)])

# LDA-RFE
library(MASS)
library(ROSE)
# 获得15个list
data_all <- data_train
result <- matrix(ncol=3)
least <- matrix(ncol=4)
least <- as.data.frame(least[-1,])
colnames(least) <- c("Features","Seed","Fold","Feature")
for(k in 99:101){
  set.seed(k)
  folds <- createFolds(data_all[,"是否发生非前哨淋巴结转移"], k=5)
  for(j in 1:5){
    fold_test <- data_all[folds[[j]],]
    fold_train <- data_all[-folds[[j]],]
    for(i in (ncol(data_all)-1):1){
      lda_train <- lda(是否发生非前哨淋巴结转移~., data=fold_train[,c(1:i,ncol(fold_train))])
      lda_test <- predict(lda_train, newdata=fold_test)
      auc <- roc.curve(fold_test$是否发生非前哨淋巴结转移, lda_test$posterior[,2])$auc
      line <- c(i, j, auc)
      result <- rbind(result, line)
      line <- data.frame(Features=i, Seed=k, Fold=j, Feature=rownames(which(abs(lda_train$scaling)==min(abs(lda_train$scaling)), arr.ind=TRUE)))
      least <- rbind(least, line)
      n <- which(abs(lda_train$scaling)==min(abs(lda_train$scaling)))
      fold_train <- fold_train[,-n]
    }
  }
}
result <- as.data.frame(result[-1,])
colnames(result) <- c("Features","Fold","AUC")
write.table(result, "lda_auc_non-sentinel.txt", sep="\t", row.names=F, quote=F)
write.table(least, "lda_importance_non-sentinel.txt", sep="\t", row.names=F, quote=F)
# 箱型图
ggplot(result, aes(x=as.factor(Features), y=AUC)) + 
  geom_boxplot(colour="#0073C2FF") + 
  scale_y_continuous(limits=c(0.50,1.00), breaks=seq(0.50,1.00,0.05)) + 
  stat_summary(fun=mean, shape=19, col='red', geom='point') +   
  xlab("Number of features") + 
  ylab("AUC (cross-validation)") + 
  theme(
    axis.text.x=element_text(size=rel(0.8),colour="black"),
    axis.text.y=element_text(size=rel(1.2),colour="black"),
    axis.title.x=element_text(size=rel(1.2),colour="black"),
    axis.title.y=element_text(size=rel(1.2),colour="black"),
    axis.line.x=element_line(size=0.5,colour="black"),
    axis.line.y=element_line(size=0.5,colour="black"),
    panel.border=element_blank(),
    panel.grid=element_blank(),
    panel.background=element_rect(fill="transparent",color=NA)
  )
model_lda_sig <- lda(是否发生非前哨淋巴结转移~cluster2+cluster1+前哨淋巴结阳性个数分组+前哨淋巴结阴性个数分组, data=data)
plot(model_lda_sig, col="#0073C2FF", type="both")

###########################################
##Logit Regression##
###########################################

rm(list=ls())

setwd("E:/")
path="E:/"

# 外部验证之前运行
result <- matrix(ncol=9)

data_input <- read.csv(paste(path, "EDC_229_35_num_non-sentinel_3_29.csv", sep=""), header=TRUE)
# cluster2+cluster1+前哨淋巴结阳性个数分组+前哨淋巴结阴性个数分组
data <- data_input[,c(3,2,31,32,35)]
data$是否发生非前哨淋巴结转移 <- factor(data$是否发生非前哨淋巴结转移, levels=c(1,0), labels=c("yes","no"))
data$是否发生非前哨淋巴结转移 <- factor(data$是否发生非前哨淋巴结转移, levels=c("no","yes"), labels=c("no","yes"))

# Train and test sets
set.seed(3)
index_test_0 <- sample(which(data$是否发生非前哨淋巴结转移=="no"), 30)
set.seed(3)
index_test_1 <- sample(which(data$是否发生非前哨淋巴结转移=="yes"), 16)
index_test <- c(index_test_0,index_test_1)
data_train <- data[-index_test,]
data_test <- data[index_test,]
data_train[,c(3:5)] <- lapply(data_train[,c(3:5)], as.factor)
data_test[,c(3:5)] <- lapply(data_test[,c(3:5)], as.factor)

# clinic only
data_train <- data_train[,c(-1,-2)]
data_test <- data_test[,c(-1,-2)]

# external validation
library(pROC)
lr_train <- glm(是否发生非前哨淋巴结转移~., family=binomial(link=logit), data=data_train)
lr_test <- predict(lr_train, type='response', newdata=data_test)
# 使用pROC包自动标出最优临界点
modelroc <- roc(data_test$是否发生非前哨淋巴结转移, lr_test, levels=c("no","yes"))
plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1,0.2), grid.col=c("green","red"), max.auc.polygon=TRUE, auc.polygon.col="skyblue", print.thres=TRUE)
thres <- as.numeric(coords(modelroc, "best", best.method="closest.topleft", ret=c("threshold", "specificity", "sensitivity"))[1,1])
# thres <- 0.5
direction <- modelroc$direction
auc <- modelroc$auc
# 评估模型的预测效果
if(direction=="<"){
  predict.results <- ifelse(lr_test>thres, 1, 0)
}
if(direction==">"){
  predict.results <- ifelse(lr_test<thres, 1, 0)
}
predict.results <- factor(predict.results, levels=c("0","1"))
table <- table(predict.results, data_test$是否发生非前哨淋巴结转移)
print(table)
tp <- table[2,2]
fn <- table[1,2]
fp <- table[2,1]
tn <- table[1,1]
acc <- (tp+tn)/(tp+fn+fp+tn)
sen <- tp/(fn+tp)
f0r <- fn/(fn+tn)
spe <- tn/(fp+tn)
fdr <- fp/(fp+tp)
rec <- tp/(tp+fn)
pre <- tp/(tp+fp)
f1 <- 2*(rec*pre)/(rec+pre)
print(paste("auc=", auc, sep=""))
print(paste("acc=", acc, sep=""))
print(paste("sen=", sen, sep=""))
print(paste("for=", f0r, sep=""))
print(paste("spe=", spe, sep=""))
print(paste("fdr=", fdr, sep=""))
print(paste("rec=", rec, sep=""))
print(paste("pre=", pre, sep=""))
print(paste("f1=", f1, sep=""))
model <- c(auc, acc, sen, f0r, spe, fdr, rec, pre, f1)
result <- rbind(result, model)

# 外部验证之后运行
result <- result[-1,]
colnames(result) <- c("AUC","ACC","SEN","FOR","SPE","FDR","REC","PRE","F1")
write.csv(result, "ex_logit_non-sentinel.csv", quote=FALSE)

# cross validation
library(pROC)
library(caret)
data_all <- data_train
result <- matrix(ncol=9)
for(j in 99:101){
  set.seed(j)
  folds <- createFolds(data_all[,"是否发生非前哨淋巴结转移"], k=5)
  for(i in 1:5){
    fold_test <- data_all[folds[[i]],] # 取folds[[i]]作为测试集
    fold_train <- data_all[-folds[[i]],] # 剩下的数据作为训练集
    lr_train <- glm(是否发生非前哨淋巴结转移~., family=binomial(link=logit), data=fold_train)
    lr_test <- predict(lr_train, type='response', newdata=fold_test)
    # 使用pROC包自动标出最优临界点
    modelroc <- roc(fold_test$是否发生非前哨淋巴结转移, lr_test, levels=c("no","yes"))
    plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1,0.2), grid.col=c("green","red"), max.auc.polygon=TRUE, auc.polygon.col="skyblue", print.thres=TRUE)
    thres <- as.numeric(coords(modelroc, "best", best.method="closest.topleft", ret=c("threshold", "specificity", "sensitivity"))[1,1])
    # thres <- 0.5
    direction <- modelroc$direction
    auc <- modelroc$auc
    # 评估模型的预测效果
    if(direction=="<"){
      predict.results <- ifelse(lr_test>thres, 1, 0)
    }
    if(direction==">"){
      predict.results <- ifelse(lr_test<thres, 1, 0)
    }
    predict.results <- factor(predict.results, levels=c("0","1"))
    table <- table(predict.results, fold_test$是否发生非前哨淋巴结转移)
    print(table)
    tp <- table[2,2]
    fn <- table[1,2]
    fp <- table[2,1]
    tn <- table[1,1]
    acc <- (tp+tn)/(tp+fn+fp+tn)
    sen <- tp/(fn+tp)
    f0r <- fn/(fn+tn)
    spe <- tn/(fp+tn)
    fdr <- fp/(fp+tp)
    rec <- tp/(tp+fn)
    pre <- tp/(tp+fp)
    f1 <- 2*(rec*pre)/(rec+pre)
    print(paste("auc=", auc, sep=""))
    print(paste("acc=", acc, sep=""))
    print(paste("sen=", sen, sep=""))
    print(paste("for=", f0r, sep=""))
    print(paste("spe=", spe, sep=""))
    print(paste("fdr=", fdr, sep=""))
    print(paste("rec=", rec, sep=""))
    print(paste("pre=", pre, sep=""))
    print(paste("f1=", f1, sep=""))
    fold <- c()
    fold <- c(fold, auc, acc, sen, f0r, spe, fdr, rec, pre, f1)
    result <- rbind(result, fold)
  }
}
result <- result[-1,]
rownames(result) <- c(1:15)
colnames(result) <- c("AUC","ACC","SEN","FOR","SPE","FDR","REC","PRE","F1")
result <- as.data.frame(result)
result <- as.data.frame(lapply(result, function(x){paste(x)}))
write.csv(result, "in_logit_non-sentinel.csv", quote=FALSE)

###########################################
##MSKCC##
###########################################

rm(list=ls())

setwd("E:/")
path="E:/"

data_input <- read.csv(paste(path, "EDC_229_35_num_non-sentinel_3_84.csv", sep=""), header=TRUE)

# 提取数据
# 手术病理大小分组+ER状态+术后病理类型+前哨淋巴结阳性个数分组+前哨淋巴结阴性个数分组
data <- data_input[,c(22,24,29,31,32,35)]
data$是否发生非前哨淋巴结转移 <- factor(data$是否发生非前哨淋巴结转移, levels=c(1,0), labels=c("yes","no"))
data$是否发生非前哨淋巴结转移 <- factor(data$是否发生非前哨淋巴结转移, levels=c("no","yes"), labels=c("no","yes"))

# Train and test sets
set.seed(3)
index_test_0 <- sample(which(data$是否发生非前哨淋巴结转移=="no"), 30)
set.seed(3)
index_test_1 <- sample(which(data$是否发生非前哨淋巴结转移=="yes"), 16)
index_test <- c(index_test_0,index_test_1)
data_train <- data[-index_test,]
data_test <- data[index_test,]
data_train[,c(1:6)] <- lapply(data_train[,c(1:6)], as.factor)
data_test[,c(1:6)] <- lapply(data_test[,c(1:6)], as.factor)

# external validation
library(pROC)
result <- matrix(ncol=9)
lr_train <- glm(是否发生非前哨淋巴结转移~手术病理大小分组+ER状态+术后病理类型+前哨淋巴结阳性个数分组+前哨淋巴结阴性个数分组, family=binomial(link=logit), data=data_train)
lr_test <- predict(lr_train, type='response', newdata=data_test)
# 使用pROC包自动标出最优临界点
modelroc <- roc(data_test$是否发生非前哨淋巴结转移, lr_test, levels=c("no","yes"))
plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1,0.2), grid.col=c("green","red"), max.auc.polygon=TRUE, auc.polygon.col="skyblue", print.thres=TRUE)
thres <- as.numeric(coords(modelroc, "best", best.method="closest.topleft", ret=c("threshold", "specificity", "sensitivity"))[1,1])
# thres <- 0.5
direction <- modelroc$direction
auc <- modelroc$auc
# 评估模型的预测效果
if(direction=="<"){
  predict.results <- ifelse(lr_test>thres, 1, 0)
}
if(direction==">"){
  predict.results <- ifelse(lr_test<thres, 1, 0)
}
predict.results <- factor(predict.results, levels=c("0","1"))
table <- table(predict.results, data_test$是否发生非前哨淋巴结转移)
print(table)
tp <- table[2,2]
fn <- table[1,2]
fp <- table[2,1]
tn <- table[1,1]
acc <- (tp+tn)/(tp+fn+fp+tn)
sen <- tp/(fn+tp)
f0r <- fn/(fn+tn)
spe <- tn/(fp+tn)
fdr <- fp/(fp+tp)
rec <- tp/(tp+fn)
pre <- tp/(tp+fp)
f1 <- 2*(rec*pre)/(rec+pre)
print(paste("auc=", auc, sep=""))
print(paste("acc=", acc, sep=""))
print(paste("sen=", sen, sep=""))
print(paste("for=", f0r, sep=""))
print(paste("spe=", spe, sep=""))
print(paste("fdr=", fdr, sep=""))
print(paste("rec=", rec, sep=""))
print(paste("pre=", pre, sep=""))
print(paste("f1=", f1, sep=""))
model <- c(auc, acc, sen, f0r, spe, fdr, rec, pre, f1)
result <- rbind(result, model)
colnames(result) <- c("AUC","ACC","SEN","FOR","SPE","FDR","REC","PRE","F1")
write.csv(result, "ex_mskcc.csv", quote=FALSE)

# cross validation
library(pROC)
library(caret)
data_all <- data_train
result <- matrix(ncol=9)
for(j in 99:101){
  set.seed(j)
  folds <- createFolds(data_all[,"是否发生非前哨淋巴结转移"], k=5)
  for(i in 1:5){
    fold_test <- data_all[folds[[i]],] # 取folds[[i]]作为测试集
    fold_train <- data_all[-folds[[i]],] # 剩下的数据作为训练集
    lr_train <- glm(是否发生非前哨淋巴结转移~手术病理大小分组+ER状态+术后病理类型+前哨淋巴结阳性个数分组+前哨淋巴结阴性个数分组, family=binomial(link=logit), data=fold_train)
    lr_test <- predict(lr_train, type='response', newdata=fold_test)
    # 使用pROC包自动标出最优临界点
    modelroc <- roc(fold_test$是否发生非前哨淋巴结转移, lr_test, levels=c("no","yes"))
    plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1,0.2), grid.col=c("green","red"), max.auc.polygon=TRUE, auc.polygon.col="skyblue", print.thres=TRUE)
    thres <- as.numeric(coords(modelroc, "best", best.method="closest.topleft", ret=c("threshold", "specificity", "sensitivity"))[1,1])
    # thres <- 0.5
    direction <- modelroc$direction
    auc <- modelroc$auc
    # 评估模型的预测效果
    if(direction=="<"){
      predict.results <- ifelse(lr_test>thres, 1, 0)
    }
    if(direction==">"){
      predict.results <- ifelse(lr_test<thres, 1, 0)
    }
    predict.results <- factor(predict.results, levels=c("0","1"))
    table <- table(predict.results, fold_test$是否发生非前哨淋巴结转移)
    print(table)
    tp <- table[2,2]
    fn <- table[1,2]
    fp <- table[2,1]
    tn <- table[1,1]
    acc <- (tp+tn)/(tp+fn+fp+tn)
    sen <- tp/(fn+tp)
    f0r <- fn/(fn+tn)
    spe <- tn/(fp+tn)
    fdr <- fp/(fp+tp)
    rec <- tp/(tp+fn)
    pre <- tp/(tp+fp)
    f1 <- 2*(rec*pre)/(rec+pre)
    print(paste("auc=", auc, sep=""))
    print(paste("acc=", acc, sep=""))
    print(paste("sen=", sen, sep=""))
    print(paste("for=", f0r, sep=""))
    print(paste("spe=", spe, sep=""))
    print(paste("fdr=", fdr, sep=""))
    print(paste("rec=", rec, sep=""))
    print(paste("pre=", pre, sep=""))
    print(paste("f1=", f1, sep=""))
    fold <- c()
    fold <- c(fold, auc, acc, sen, f0r, spe, fdr, rec, pre, f1)
    result <- rbind(result, fold)
  }
}
result <- result[-1,]
rownames(result) <- c(1:15)
colnames(result) <- c("AUC","ACC","SEN","FOR","SPE","FDR","REC","PRE","F1")
result <- as.data.frame(result)
result <- as.data.frame(lapply(result, function(x){paste(x)}))
write.csv(result, "in_mskcc.csv", quote = FALSE)

###########################################
##MOU##
###########################################

rm(list=ls())

setwd("E:/")
path="E:/"

data_input <- read.csv(paste(path, "EDC_229_35_num_non-sentinel_3_84.csv", sep=""), header=TRUE)

# 提取数据
# 手术病理大小分组+术后病理类型+前哨淋巴结阳性个数分组+前哨淋巴结阴性个数分组
data <- data_input[,c(22,29,31,32,35)]
data$是否发生非前哨淋巴结转移 <- factor(data$是否发生非前哨淋巴结转移, levels=c(1,0), labels=c("yes","no"))
data$是否发生非前哨淋巴结转移 <- factor(data$是否发生非前哨淋巴结转移, levels=c("no","yes"), labels=c("no","yes"))

# Train and test sets
set.seed(3)
index_test_0 <- sample(which(data$是否发生非前哨淋巴结转移=="no"), 30)
set.seed(3)
index_test_1 <- sample(which(data$是否发生非前哨淋巴结转移=="yes"), 16)
index_test <- c(index_test_0,index_test_1)
data_train <- data[-index_test,]
data_test <- data[index_test,]
data_train[,c(1:5)] <- lapply(data_train[,c(1:5)], as.factor)
data_test[,c(1:5)] <- lapply(data_test[,c(1:5)], as.factor)

# external validation
library(pROC)
result <- matrix(ncol=9)
lr_train <- glm(是否发生非前哨淋巴结转移~手术病理大小分组+术后病理类型+前哨淋巴结阳性个数分组+前哨淋巴结阴性个数分组, family=binomial(link=logit), data=data_train)
lr_test <- predict(lr_train, type='response', newdata=data_test)
# 使用pROC包自动标出最优临界点
modelroc <- roc(data_test$是否发生非前哨淋巴结转移, lr_test, levels=c("no","yes"))
plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1,0.2), grid.col=c("green","red"), max.auc.polygon=TRUE, auc.polygon.col="skyblue", print.thres=TRUE)
thres <- as.numeric(coords(modelroc, "best", best.method="closest.topleft", ret=c("threshold", "specificity", "sensitivity"))[1,1])
# thres <- 0.5
direction <- modelroc$direction
auc <- modelroc$auc
# 评估模型的预测效果
if(direction=="<"){
  predict.results <- ifelse(lr_test>thres, 1, 0)
}
if(direction==">"){
  predict.results <- ifelse(lr_test<thres, 1, 0)
}
predict.results <- factor(predict.results, levels=c("0","1"))
table <- table(predict.results, data_test$是否发生非前哨淋巴结转移)
print(table)
tp <- table[2,2]
fn <- table[1,2]
fp <- table[2,1]
tn <- table[1,1]
acc <- (tp+tn)/(tp+fn+fp+tn)
sen <- tp/(fn+tp)
f0r <- fn/(fn+tn)
spe <- tn/(fp+tn)
fdr <- fp/(fp+tp)
rec <- tp/(tp+fn)
pre <- tp/(tp+fp)
f1 <- 2*(rec*pre)/(rec+pre)
print(paste("auc=", auc, sep=""))
print(paste("acc=", acc, sep=""))
print(paste("sen=", sen, sep=""))
print(paste("for=", f0r, sep=""))
print(paste("spe=", spe, sep=""))
print(paste("fdr=", fdr, sep=""))
print(paste("rec=", rec, sep=""))
print(paste("pre=", pre, sep=""))
print(paste("f1=", f1, sep=""))
model <- c(auc, acc, sen, f0r, spe, fdr, rec, pre, f1)
result <- rbind(result, model)
colnames(result) <- c("AUC","ACC","SEN","FOR","SPE","FDR","REC","PRE","F1")
write.csv(result, "ex_mou.csv", quote=FALSE)

# cross validation
library(pROC)
library(caret)
data_all <- data_train
result <- matrix(ncol=9)
for(j in 99:101){
  set.seed(j)
  folds <- createFolds(data_all[,"是否发生非前哨淋巴结转移"], k=5)
  for(i in 1:5){
    fold_test <- data_all[folds[[i]],] # 取folds[[i]]作为测试集
    fold_train <- data_all[-folds[[i]],] # 剩下的数据作为训练集
    lr_train <- glm(是否发生非前哨淋巴结转移~手术病理大小分组+术后病理类型+前哨淋巴结阳性个数分组+前哨淋巴结阴性个数分组, family=binomial(link=logit), data=fold_train)
    lr_test <- predict(lr_train, type='response', newdata=fold_test)
    # 使用pROC包自动标出最优临界点
    modelroc <- roc(fold_test$是否发生非前哨淋巴结转移, lr_test, levels=c("no","yes"))
    plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1,0.2), grid.col=c("green","red"), max.auc.polygon=TRUE, auc.polygon.col="skyblue", print.thres=TRUE)
    thres <- as.numeric(coords(modelroc, "best", best.method="closest.topleft", ret=c("threshold", "specificity", "sensitivity"))[1,1])
    # thres <- 0.5
    direction <- modelroc$direction
    auc <- modelroc$auc
    # 评估模型的预测效果
    if(direction=="<"){
      predict.results <- ifelse(lr_test>thres, 1, 0)
    }
    if(direction==">"){
      predict.results <- ifelse(lr_test<thres, 1, 0)
    }
    predict.results <- factor(predict.results, levels=c("0","1"))
    table <- table(predict.results, fold_test$是否发生非前哨淋巴结转移)
    print(table)
    tp <- table[2,2]
    fn <- table[1,2]
    fp <- table[2,1]
    tn <- table[1,1]
    acc <- (tp+tn)/(tp+fn+fp+tn)
    sen <- tp/(fn+tp)
    f0r <- fn/(fn+tn)
    spe <- tn/(fp+tn)
    fdr <- fp/(fp+tp)
    rec <- tp/(tp+fn)
    pre <- tp/(tp+fp)
    f1 <- 2*(rec*pre)/(rec+pre)
    print(paste("auc=", auc, sep=""))
    print(paste("acc=", acc, sep=""))
    print(paste("sen=", sen, sep=""))
    print(paste("for=", f0r, sep=""))
    print(paste("spe=", spe, sep=""))
    print(paste("fdr=", fdr, sep=""))
    print(paste("rec=", rec, sep=""))
    print(paste("pre=", pre, sep=""))
    print(paste("f1=", f1, sep=""))
    fold <- c()
    fold <- c(fold, auc, acc, sen, f0r, spe, fdr, rec, pre, f1)
    result <- rbind(result, fold)
  }
}
result <- result[-1,]
rownames(result) <- c(1:15)
colnames(result) <- c("AUC","ACC","SEN","FOR","SPE","FDR","REC","PRE","F1")
result <- as.data.frame(result)
result <- as.data.frame(lapply(result, function(x){paste(x)}))
write.csv(result, "in_mou.csv", quote = FALSE)

###########################################
##Ljubljana##
###########################################

rm(list=ls())

setwd("E:/")
path="E:/"

data_input <- read.csv(paste(path, "EDC_229_35_num_non-sentinel_3_84.csv", sep=""), header=TRUE)

# 提取数据
# 手术病理大小分组+前哨淋巴结阳性个数分组+前哨淋巴结阴性个数分组
data <- data_input[,c(22,31,32,35)]
data$是否发生非前哨淋巴结转移 <- factor(data$是否发生非前哨淋巴结转移, levels=c(1,0), labels=c("yes","no"))
data$是否发生非前哨淋巴结转移 <- factor(data$是否发生非前哨淋巴结转移, levels=c("no","yes"), labels=c("no","yes"))

# Train and test sets
set.seed(3)
index_test_0 <- sample(which(data$是否发生非前哨淋巴结转移=="no"), 30)
set.seed(3)
index_test_1 <- sample(which(data$是否发生非前哨淋巴结转移=="yes"), 16)
index_test <- c(index_test_0,index_test_1)
data_train <- data[-index_test,]
data_test <- data[index_test,]
data_train[,c(1:4)] <- lapply(data_train[,c(1:4)], as.factor)
data_test[,c(1:4)] <- lapply(data_test[,c(1:4)], as.factor)

# external validation
library(pROC)
result <- matrix(ncol=9)
lr_train <- glm(是否发生非前哨淋巴结转移~手术病理大小分组+前哨淋巴结阳性个数分组+前哨淋巴结阴性个数分组, family=binomial(link=logit), data=data_train)
lr_test <- predict(lr_train, type='response', newdata=data_test)
# 使用pROC包自动标出最优临界点
modelroc <- roc(data_test$是否发生非前哨淋巴结转移, lr_test, levels=c("no","yes"))
plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1,0.2), grid.col=c("green","red"), max.auc.polygon=TRUE, auc.polygon.col="skyblue", print.thres=TRUE)
thres <- as.numeric(coords(modelroc, "best", best.method="closest.topleft", ret=c("threshold", "specificity", "sensitivity"))[1,1])
# thres <- 0.5
direction <- modelroc$direction
auc <- modelroc$auc
# 评估模型的预测效果
if(direction=="<"){
  predict.results <- ifelse(lr_test>thres, 1, 0)
}
if(direction==">"){
  predict.results <- ifelse(lr_test<thres, 1, 0)
}
predict.results <- factor(predict.results, levels=c("0","1"))
table <- table(predict.results, data_test$是否发生非前哨淋巴结转移)
print(table)
tp <- table[2,2]
fn <- table[1,2]
fp <- table[2,1]
tn <- table[1,1]
acc <- (tp+tn)/(tp+fn+fp+tn)
sen <- tp/(fn+tp)
f0r <- fn/(fn+tn)
spe <- tn/(fp+tn)
fdr <- fp/(fp+tp)
rec <- tp/(tp+fn)
pre <- tp/(tp+fp)
f1 <- 2*(rec*pre)/(rec+pre)
print(paste("auc=", auc, sep=""))
print(paste("acc=", acc, sep=""))
print(paste("sen=", sen, sep=""))
print(paste("for=", f0r, sep=""))
print(paste("spe=", spe, sep=""))
print(paste("fdr=", fdr, sep=""))
print(paste("rec=", rec, sep=""))
print(paste("pre=", pre, sep=""))
print(paste("f1=", f1, sep=""))
model <- c(auc, acc, sen, f0r, spe, fdr, rec, pre, f1)
result <- rbind(result, model)
colnames(result) <- c("AUC","ACC","SEN","FOR","SPE","FDR","REC","PRE","F1")
write.csv(result, "ex_ljubljana.csv", quote=FALSE)

# cross validation
library(pROC)
library(caret)
data_all <- data_train
result <- matrix(ncol=9)
for(j in 99:101){
  set.seed(j)
  folds <- createFolds(data_all[,"是否发生非前哨淋巴结转移"], k=5)
  for(i in 1:5){
    fold_test <- data_all[folds[[i]],] # 取folds[[i]]作为测试集
    fold_train <- data_all[-folds[[i]],] # 剩下的数据作为训练集
    lr_train <- glm(是否发生非前哨淋巴结转移~手术病理大小分组+前哨淋巴结阳性个数分组+前哨淋巴结阴性个数分组, family=binomial(link=logit), data=fold_train)
    lr_test <- predict(lr_train, type='response', newdata=fold_test)
    # 使用pROC包自动标出最优临界点
    modelroc <- roc(fold_test$是否发生非前哨淋巴结转移, lr_test, levels=c("no","yes"))
    plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1,0.2), grid.col=c("green","red"), max.auc.polygon=TRUE, auc.polygon.col="skyblue", print.thres=TRUE)
    thres <- as.numeric(coords(modelroc, "best", best.method="closest.topleft", ret=c("threshold", "specificity", "sensitivity"))[1,1])
    # thres <- 0.5
    direction <- modelroc$direction
    auc <- modelroc$auc
    # 评估模型的预测效果
    if(direction=="<"){
      predict.results <- ifelse(lr_test>thres, 1, 0)
    }
    if(direction==">"){
      predict.results <- ifelse(lr_test<thres, 1, 0)
    }
    predict.results <- factor(predict.results, levels=c("0","1"))
    table <- table(predict.results, fold_test$是否发生非前哨淋巴结转移)
    print(table)
    tp <- table[2,2]
    fn <- table[1,2]
    fp <- table[2,1]
    tn <- table[1,1]
    acc <- (tp+tn)/(tp+fn+fp+tn)
    sen <- tp/(fn+tp)
    f0r <- fn/(fn+tn)
    spe <- tn/(fp+tn)
    fdr <- fp/(fp+tp)
    rec <- tp/(tp+fn)
    pre <- tp/(tp+fp)
    f1 <- 2*(rec*pre)/(rec+pre)
    print(paste("auc=", auc, sep=""))
    print(paste("acc=", acc, sep=""))
    print(paste("sen=", sen, sep=""))
    print(paste("for=", f0r, sep=""))
    print(paste("spe=", spe, sep=""))
    print(paste("fdr=", fdr, sep=""))
    print(paste("rec=", rec, sep=""))
    print(paste("pre=", pre, sep=""))
    print(paste("f1=", f1, sep=""))
    fold <- c()
    fold <- c(fold, auc, acc, sen, f0r, spe, fdr, rec, pre, f1)
    result <- rbind(result, fold)
  }
}
result <- result[-1,]
rownames(result) <- c(1:15)
colnames(result) <- c("AUC","ACC","SEN","FOR","SPE","FDR","REC","PRE","F1")
result <- as.data.frame(result)
result <- as.data.frame(lapply(result, function(x){paste(x)}))
write.csv(result, "in_ljubljana.csv", quote = FALSE)

###########################################
##Cambridge##
###########################################

rm(list=ls())

setwd("E:/")
path="E:/"

data_input <- read.csv(paste(path, "EDC_229_35_num_non-sentinel_3_84.csv", sep=""), header=TRUE)

# 提取数据
# 前哨淋巴结阳性比例分组
data <- data_input[,c(33,35)]
data$是否发生非前哨淋巴结转移 <- factor(data$是否发生非前哨淋巴结转移, levels=c(1,0), labels=c("yes","no"))
data$是否发生非前哨淋巴结转移 <- factor(data$是否发生非前哨淋巴结转移, levels=c("no","yes"), labels=c("no","yes"))

# Train and test sets
set.seed(3)
index_test_0 <- sample(which(data$是否发生非前哨淋巴结转移=="no"), 30)
set.seed(3)
index_test_1 <- sample(which(data$是否发生非前哨淋巴结转移=="yes"), 16)
index_test <- c(index_test_0,index_test_1)
data_train <- data[-index_test,]
data_test <- data[index_test,]
data_train[,c(1:2)] <- lapply(data_train[,c(1:2)], as.factor)
data_test[,c(1:2)] <- lapply(data_test[,c(1:2)], as.factor)

# external validation
library(pROC)
result <- matrix(ncol=9)
lr_train <- glm(是否发生非前哨淋巴结转移~前哨淋巴结阳性比例分组, family=binomial(link=logit), data=data_train)
lr_test <- predict(lr_train, type='response', newdata=data_test)
# 使用pROC包自动标出最优临界点
modelroc <- roc(data_test$是否发生非前哨淋巴结转移, lr_test, levels=c("no","yes"))
plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1,0.2), grid.col=c("green","red"), max.auc.polygon=TRUE, auc.polygon.col="skyblue", print.thres=TRUE)
thres <- as.numeric(coords(modelroc, "best", best.method="closest.topleft", ret=c("threshold", "specificity", "sensitivity"))[1,1])
# thres <- 0.5
direction <- modelroc$direction
auc <- modelroc$auc
# 评估模型的预测效果
if(direction=="<"){
  predict.results <- ifelse(lr_test>thres, 1, 0)
}
if(direction==">"){
  predict.results <- ifelse(lr_test<thres, 1, 0)
}
predict.results <- factor(predict.results, levels=c("0","1"))
table <- table(predict.results, data_test$是否发生非前哨淋巴结转移)
print(table)
tp <- table[2,2]
fn <- table[1,2]
fp <- table[2,1]
tn <- table[1,1]
acc <- (tp+tn)/(tp+fn+fp+tn)
sen <- tp/(fn+tp)
f0r <- fn/(fn+tn)
spe <- tn/(fp+tn)
fdr <- fp/(fp+tp)
rec <- tp/(tp+fn)
pre <- tp/(tp+fp)
f1 <- 2*(rec*pre)/(rec+pre)
print(paste("auc=", auc, sep=""))
print(paste("acc=", acc, sep=""))
print(paste("sen=", sen, sep=""))
print(paste("for=", f0r, sep=""))
print(paste("spe=", spe, sep=""))
print(paste("fdr=", fdr, sep=""))
print(paste("rec=", rec, sep=""))
print(paste("pre=", pre, sep=""))
print(paste("f1=", f1, sep=""))
model <- c(auc, acc, sen, f0r, spe, fdr, rec, pre, f1)
result <- rbind(result, model)
colnames(result) <- c("AUC","ACC","SEN","FOR","SPE","FDR","REC","PRE","F1")
write.csv(result, "ex_cambridge.csv", quote=FALSE)

# cross validation
library(pROC)
library(caret)
data_all <- data_train
result <- matrix(ncol=9)
for(j in 99:101){
  set.seed(j)
  folds <- createFolds(data_all[,"是否发生非前哨淋巴结转移"], k=5)
  for(i in 1:5){
    fold_test <- data_all[folds[[i]],] # 取folds[[i]]作为测试集
    fold_train <- data_all[-folds[[i]],] # 剩下的数据作为训练集
    lr_train <- glm(是否发生非前哨淋巴结转移~前哨淋巴结阳性比例分组, family=binomial(link=logit), data=fold_train)
    lr_test <- predict(lr_train, type='response', newdata=fold_test)
    # 使用pROC包自动标出最优临界点
    modelroc <- roc(fold_test$是否发生非前哨淋巴结转移, lr_test, levels=c("no","yes"))
    plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1,0.2), grid.col=c("green","red"), max.auc.polygon=TRUE, auc.polygon.col="skyblue", print.thres=TRUE)
    thres <- as.numeric(coords(modelroc, "best", best.method="closest.topleft", ret=c("threshold", "specificity", "sensitivity"))[1,1])
    # thres <- 0.5
    direction <- modelroc$direction
    auc <- modelroc$auc
    # 评估模型的预测效果
    if(direction=="<"){
      predict.results <- ifelse(lr_test>thres, 1, 0)
    }
    if(direction==">"){
      predict.results <- ifelse(lr_test<thres, 1, 0)
    }
    predict.results <- factor(predict.results, levels=c("0","1"))
    table <- table(predict.results, fold_test$是否发生非前哨淋巴结转移)
    print(table)
    tp <- table[2,2]
    fn <- table[1,2]
    fp <- table[2,1]
    tn <- table[1,1]
    acc <- (tp+tn)/(tp+fn+fp+tn)
    sen <- tp/(fn+tp)
    f0r <- fn/(fn+tn)
    spe <- tn/(fp+tn)
    fdr <- fp/(fp+tp)
    rec <- tp/(tp+fn)
    pre <- tp/(tp+fp)
    f1 <- 2*(rec*pre)/(rec+pre)
    print(paste("auc=", auc, sep=""))
    print(paste("acc=", acc, sep=""))
    print(paste("sen=", sen, sep=""))
    print(paste("for=", f0r, sep=""))
    print(paste("spe=", spe, sep=""))
    print(paste("fdr=", fdr, sep=""))
    print(paste("rec=", rec, sep=""))
    print(paste("pre=", pre, sep=""))
    print(paste("f1=", f1, sep=""))
    fold <- c()
    fold <- c(fold, auc, acc, sen, f0r, spe, fdr, rec, pre, f1)
    result <- rbind(result, fold)
  }
}
result <- result[-1,]
rownames(result) <- c(1:15)
colnames(result) <- c("AUC","ACC","SEN","FOR","SPE","FDR","REC","PRE","F1")
result <- as.data.frame(result)
result <- as.data.frame(lapply(result, function(x){paste(x)}))
write.csv(result, "in_cambridge.csv", quote = FALSE)

###########################################
##MDA##
###########################################

rm(list=ls())

setwd("E:/")
path="E:/"

data_input <- read.csv(paste(path, "EDC_229_35_num_non-sentinel_3_84.csv", sep=""), header=TRUE)

# 提取数据
# 手术病理大小分组+前哨淋巴结活检总数分组
data <- data_input[,c(22,30,35)]
data$是否发生非前哨淋巴结转移 <- factor(data$是否发生非前哨淋巴结转移, levels=c(1,0), labels=c("yes","no"))
data$是否发生非前哨淋巴结转移 <- factor(data$是否发生非前哨淋巴结转移, levels=c("no","yes"), labels=c("no","yes"))

# Train and test sets
set.seed(3)
index_test_0 <- sample(which(data$是否发生非前哨淋巴结转移=="no"), 30)
set.seed(3)
index_test_1 <- sample(which(data$是否发生非前哨淋巴结转移=="yes"), 16)
index_test <- c(index_test_0,index_test_1)
data_train <- data[-index_test,]
data_test <- data[index_test,]
data_train[,c(1:3)] <- lapply(data_train[,c(1:3)], as.factor)
data_test[,c(1:3)] <- lapply(data_test[,c(1:3)], as.factor)

# external validation
library(pROC)
result <- matrix(ncol=9)
lr_train <- glm(是否发生非前哨淋巴结转移~手术病理大小分组+前哨淋巴结活检总数分组, family=binomial(link=logit), data=data_train)
lr_test <- predict(lr_train, type='response', newdata=data_test)
# 使用pROC包自动标出最优临界点
modelroc <- roc(data_test$是否发生非前哨淋巴结转移, lr_test, levels=c("no","yes"))
plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1,0.2), grid.col=c("green","red"), max.auc.polygon=TRUE, auc.polygon.col="skyblue", print.thres=TRUE)
thres <- as.numeric(coords(modelroc, "best", best.method="closest.topleft", ret=c("threshold", "specificity", "sensitivity"))[1,1])
# thres <- 0.5
direction <- modelroc$direction
auc <- modelroc$auc
# 评估模型的预测效果
if(direction=="<"){
  predict.results <- ifelse(lr_test>thres, 1, 0)
}
if(direction==">"){
  predict.results <- ifelse(lr_test<thres, 1, 0)
}
predict.results <- factor(predict.results, levels=c("0","1"))
table <- table(predict.results, data_test$是否发生非前哨淋巴结转移)
print(table)
tp <- table[2,2]
fn <- table[1,2]
fp <- table[2,1]
tn <- table[1,1]
acc <- (tp+tn)/(tp+fn+fp+tn)
sen <- tp/(fn+tp)
f0r <- fn/(fn+tn)
spe <- tn/(fp+tn)
fdr <- fp/(fp+tp)
rec <- tp/(tp+fn)
pre <- tp/(tp+fp)
f1 <- 2*(rec*pre)/(rec+pre)
print(paste("auc=", auc, sep=""))
print(paste("acc=", acc, sep=""))
print(paste("sen=", sen, sep=""))
print(paste("for=", f0r, sep=""))
print(paste("spe=", spe, sep=""))
print(paste("fdr=", fdr, sep=""))
print(paste("rec=", rec, sep=""))
print(paste("pre=", pre, sep=""))
print(paste("f1=", f1, sep=""))
model <- c(auc, acc, sen, f0r, spe, fdr, rec, pre, f1)
result <- rbind(result, model)
colnames(result) <- c("AUC","ACC","SEN","FOR","SPE","FDR","REC","PRE","F1")
write.csv(result, "ex_mda.csv", quote=FALSE)

# cross validation
library(pROC)
library(caret)
data_all <- data_train
result <- matrix(ncol=9)
for(j in 99:101){
  set.seed(j)
  folds <- createFolds(data_all[,"是否发生非前哨淋巴结转移"], k=5)
  for(i in 1:5){
    fold_test <- data_all[folds[[i]],] # 取folds[[i]]作为测试集
    fold_train <- data_all[-folds[[i]],] # 剩下的数据作为训练集
    lr_train <- glm(是否发生非前哨淋巴结转移~手术病理大小分组+前哨淋巴结活检总数分组, family=binomial(link=logit), data=fold_train)
    lr_test <- predict(lr_train, type='response', newdata=fold_test)
    # 使用pROC包自动标出最优临界点
    modelroc <- roc(fold_test$是否发生非前哨淋巴结转移, lr_test, levels=c("no","yes"))
    plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1,0.2), grid.col=c("green","red"), max.auc.polygon=TRUE, auc.polygon.col="skyblue", print.thres=TRUE)
    thres <- as.numeric(coords(modelroc, "best", best.method="closest.topleft", ret=c("threshold", "specificity", "sensitivity"))[1,1])
    # thres <- 0.5
    direction <- modelroc$direction
    auc <- modelroc$auc
    # 评估模型的预测效果
    if(direction=="<"){
      predict.results <- ifelse(lr_test>thres, 1, 0)
    }
    if(direction==">"){
      predict.results <- ifelse(lr_test<thres, 1, 0)
    }
    predict.results <- factor(predict.results, levels=c("0","1"))
    table <- table(predict.results, fold_test$是否发生非前哨淋巴结转移)
    print(table)
    tp <- table[2,2]
    fn <- table[1,2]
    fp <- table[2,1]
    tn <- table[1,1]
    acc <- (tp+tn)/(tp+fn+fp+tn)
    sen <- tp/(fn+tp)
    f0r <- fn/(fn+tn)
    spe <- tn/(fp+tn)
    fdr <- fp/(fp+tp)
    rec <- tp/(tp+fn)
    pre <- tp/(tp+fp)
    f1 <- 2*(rec*pre)/(rec+pre)
    print(paste("auc=", auc, sep=""))
    print(paste("acc=", acc, sep=""))
    print(paste("sen=", sen, sep=""))
    print(paste("for=", f0r, sep=""))
    print(paste("spe=", spe, sep=""))
    print(paste("fdr=", fdr, sep=""))
    print(paste("rec=", rec, sep=""))
    print(paste("pre=", pre, sep=""))
    print(paste("f1=", f1, sep=""))
    fold <- c()
    fold <- c(fold, auc, acc, sen, f0r, spe, fdr, rec, pre, f1)
    result <- rbind(result, fold)
  }
}
result <- result[-1,]
rownames(result) <- c(1:15)
colnames(result) <- c("AUC","ACC","SEN","FOR","SPE","FDR","REC","PRE","F1")
result <- as.data.frame(result)
result <- as.data.frame(lapply(result, function(x){paste(x)}))
write.csv(result, "in_mda.csv", quote = FALSE)

###########################################
##SNUH##
###########################################

rm(list=ls())

setwd("E:/")
path="E:/"

data_input <- read.csv(paste(path, "EDC_229_35_num_non-sentinel_3_84.csv", sep=""), header=TRUE)

# 提取数据
# T分期+前哨淋巴结阳性个数分组+前哨淋巴结阴性个数分组
data <- data_input[,c(28,31,32,35)]
data$是否发生非前哨淋巴结转移 <- factor(data$是否发生非前哨淋巴结转移, levels=c(1,0), labels=c("yes","no"))
data$是否发生非前哨淋巴结转移 <- factor(data$是否发生非前哨淋巴结转移, levels=c("no","yes"), labels=c("no","yes"))

# Train and test sets
set.seed(3)
index_test_0 <- sample(which(data$是否发生非前哨淋巴结转移=="no"), 30)
set.seed(3)
index_test_1 <- sample(which(data$是否发生非前哨淋巴结转移=="yes"), 16)
index_test <- c(index_test_0,index_test_1)
data_train <- data[-index_test,]
data_test <- data[index_test,]
data_train[,c(1:4)] <- lapply(data_train[,c(1:4)], as.factor)
data_test[,c(1:4)] <- lapply(data_test[,c(1:4)], as.factor)

# external validation
library(pROC)
result <- matrix(ncol=9)
lr_train <- glm(是否发生非前哨淋巴结转移~T分期+前哨淋巴结阳性个数分组+前哨淋巴结阴性个数分组, family=binomial(link=logit), data=data_train)
lr_test <- predict(lr_train, type='response', newdata=data_test)
# 使用pROC包自动标出最优临界点
modelroc <- roc(data_test$是否发生非前哨淋巴结转移, lr_test, levels=c("no","yes"))
plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1,0.2), grid.col=c("green","red"), max.auc.polygon=TRUE, auc.polygon.col="skyblue", print.thres=TRUE)
thres <- as.numeric(coords(modelroc, "best", best.method="closest.topleft", ret=c("threshold", "specificity", "sensitivity"))[1,1])
# thres <- 0.5
direction <- modelroc$direction
auc <- modelroc$auc
# 评估模型的预测效果
if(direction=="<"){
  predict.results <- ifelse(lr_test>thres, 1, 0)
}
if(direction==">"){
  predict.results <- ifelse(lr_test<thres, 1, 0)
}
predict.results <- factor(predict.results, levels=c("0","1"))
table <- table(predict.results, data_test$是否发生非前哨淋巴结转移)
print(table)
tp <- table[2,2]
fn <- table[1,2]
fp <- table[2,1]
tn <- table[1,1]
acc <- (tp+tn)/(tp+fn+fp+tn)
sen <- tp/(fn+tp)
f0r <- fn/(fn+tn)
spe <- tn/(fp+tn)
fdr <- fp/(fp+tp)
rec <- tp/(tp+fn)
pre <- tp/(tp+fp)
f1 <- 2*(rec*pre)/(rec+pre)
print(paste("auc=", auc, sep=""))
print(paste("acc=", acc, sep=""))
print(paste("sen=", sen, sep=""))
print(paste("for=", f0r, sep=""))
print(paste("spe=", spe, sep=""))
print(paste("fdr=", fdr, sep=""))
print(paste("rec=", rec, sep=""))
print(paste("pre=", pre, sep=""))
print(paste("f1=", f1, sep=""))
model <- c(auc, acc, sen, f0r, spe, fdr, rec, pre, f1)
result <- rbind(result, model)
colnames(result) <- c("AUC","ACC","SEN","FOR","SPE","FDR","REC","PRE","F1")
write.csv(result, "ex_snuh.csv", quote=FALSE)

# cross validation
library(pROC)
library(caret)
data_all <- data_train
result <- matrix(ncol=9)
for(j in 99:101){
  set.seed(j)
  folds <- createFolds(data_all[,"是否发生非前哨淋巴结转移"], k=5)
  for(i in 1:5){
    fold_test <- data_all[folds[[i]],] # 取folds[[i]]作为测试集
    fold_train <- data_all[-folds[[i]],] # 剩下的数据作为训练集
    lr_train <- glm(是否发生非前哨淋巴结转移~T分期+前哨淋巴结阳性个数分组+前哨淋巴结阴性个数分组, family=binomial(link=logit), data=fold_train)
    lr_test <- predict(lr_train, type='response', newdata=fold_test)
    # 使用pROC包自动标出最优临界点
    modelroc <- roc(fold_test$是否发生非前哨淋巴结转移, lr_test, levels=c("no","yes"))
    plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1,0.2), grid.col=c("green","red"), max.auc.polygon=TRUE, auc.polygon.col="skyblue", print.thres=TRUE)
    thres <- as.numeric(coords(modelroc, "best", best.method="closest.topleft", ret=c("threshold", "specificity", "sensitivity"))[1,1])
    # thres <- 0.5
    direction <- modelroc$direction
    auc <- modelroc$auc
    # 评估模型的预测效果
    if(direction=="<"){
      predict.results <- ifelse(lr_test>thres, 1, 0)
    }
    if(direction==">"){
      predict.results <- ifelse(lr_test<thres, 1, 0)
    }
    predict.results <- factor(predict.results, levels=c("0","1"))
    table <- table(predict.results, fold_test$是否发生非前哨淋巴结转移)
    print(table)
    tp <- table[2,2]
    fn <- table[1,2]
    fp <- table[2,1]
    tn <- table[1,1]
    acc <- (tp+tn)/(tp+fn+fp+tn)
    sen <- tp/(fn+tp)
    f0r <- fn/(fn+tn)
    spe <- tn/(fp+tn)
    fdr <- fp/(fp+tp)
    rec <- tp/(tp+fn)
    pre <- tp/(tp+fp)
    f1 <- 2*(rec*pre)/(rec+pre)
    print(paste("auc=", auc, sep=""))
    print(paste("acc=", acc, sep=""))
    print(paste("sen=", sen, sep=""))
    print(paste("for=", f0r, sep=""))
    print(paste("spe=", spe, sep=""))
    print(paste("fdr=", fdr, sep=""))
    print(paste("rec=", rec, sep=""))
    print(paste("pre=", pre, sep=""))
    print(paste("f1=", f1, sep=""))
    fold <- c()
    fold <- c(fold, auc, acc, sen, f0r, spe, fdr, rec, pre, f1)
    result <- rbind(result, fold)
  }
}
result <- result[-1,]
rownames(result) <- c(1:15)
colnames(result) <- c("AUC","ACC","SEN","FOR","SPE","FDR","REC","PRE","F1")
result <- as.data.frame(result)
result <- as.data.frame(lapply(result, function(x){paste(x)}))
write.csv(result, "in_snuh.csv", quote = FALSE)

###########################################
##Tenon##
###########################################

rm(list=ls())

setwd("E:/")
path="E:/"

data_input <- read.csv(paste(path, "EDC_229_35_num_non-sentinel_3_84.csv", sep=""), header=TRUE)

# 提取数据
# 手术病理大小分组+前哨淋巴结阳性比例分组
data <- data_input[,c(22,33,35)]
data$是否发生非前哨淋巴结转移 <- factor(data$是否发生非前哨淋巴结转移, levels=c(1,0), labels=c("yes","no"))
data$是否发生非前哨淋巴结转移 <- factor(data$是否发生非前哨淋巴结转移, levels=c("no","yes"), labels=c("no","yes"))

# Train and test sets
set.seed(3)
index_test_0 <- sample(which(data$是否发生非前哨淋巴结转移=="no"), 30)
set.seed(3)
index_test_1 <- sample(which(data$是否发生非前哨淋巴结转移=="yes"), 16)
index_test <- c(index_test_0,index_test_1)
data_train <- data[-index_test,]
data_test <- data[index_test,]
data_train[,c(1:3)] <- lapply(data_train[,c(1:3)], as.factor)
data_test[,c(1:3)] <- lapply(data_test[,c(1:3)], as.factor)

# external validation
library(pROC)
result <- matrix(ncol=9)
lr_train <- glm(是否发生非前哨淋巴结转移~手术病理大小分组+前哨淋巴结阳性比例分组, family=binomial(link=logit), data=data_train)
lr_test <- predict(lr_train, type='response', newdata=data_test)
# 使用pROC包自动标出最优临界点
modelroc <- roc(data_test$是否发生非前哨淋巴结转移, lr_test, levels=c("no","yes"))
plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1,0.2), grid.col=c("green","red"), max.auc.polygon=TRUE, auc.polygon.col="skyblue", print.thres=TRUE)
thres <- as.numeric(coords(modelroc, "best", best.method="closest.topleft", ret=c("threshold", "specificity", "sensitivity"))[1,1])
# thres <- 0.5
direction <- modelroc$direction
auc <- modelroc$auc
# 评估模型的预测效果
if(direction=="<"){
  predict.results <- ifelse(lr_test>thres, 1, 0)
}
if(direction==">"){
  predict.results <- ifelse(lr_test<thres, 1, 0)
}
predict.results <- factor(predict.results, levels=c("0","1"))
table <- table(predict.results, data_test$是否发生非前哨淋巴结转移)
print(table)
tp <- table[2,2]
fn <- table[1,2]
fp <- table[2,1]
tn <- table[1,1]
acc <- (tp+tn)/(tp+fn+fp+tn)
sen <- tp/(fn+tp)
f0r <- fn/(fn+tn)
spe <- tn/(fp+tn)
fdr <- fp/(fp+tp)
rec <- tp/(tp+fn)
pre <- tp/(tp+fp)
f1 <- 2*(rec*pre)/(rec+pre)
print(paste("auc=", auc, sep=""))
print(paste("acc=", acc, sep=""))
print(paste("sen=", sen, sep=""))
print(paste("for=", f0r, sep=""))
print(paste("spe=", spe, sep=""))
print(paste("fdr=", fdr, sep=""))
print(paste("rec=", rec, sep=""))
print(paste("pre=", pre, sep=""))
print(paste("f1=", f1, sep=""))
model <- c(auc, acc, sen, f0r, spe, fdr, rec, pre, f1)
result <- rbind(result, model)
colnames(result) <- c("AUC","ACC","SEN","FOR","SPE","FDR","REC","PRE","F1")
write.csv(result, "ex_tenon.csv", quote=FALSE)

# cross validation
library(pROC)
library(caret)
data_all <- data_train
result <- matrix(ncol=9)
for(j in 99:101){
  set.seed(j)
  folds <- createFolds(data_all[,"是否发生非前哨淋巴结转移"], k=5)
  for(i in 1:5){
    fold_test <- data_all[folds[[i]],] # 取folds[[i]]作为测试集
    fold_train <- data_all[-folds[[i]],] # 剩下的数据作为训练集
    lr_train <- glm(是否发生非前哨淋巴结转移~手术病理大小分组+前哨淋巴结阳性比例分组, family=binomial(link=logit), data=fold_train)
    lr_test <- predict(lr_train, type='response', newdata=fold_test)
    # 使用pROC包自动标出最优临界点
    modelroc <- roc(fold_test$是否发生非前哨淋巴结转移, lr_test, levels=c("no","yes"))
    plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1,0.2), grid.col=c("green","red"), max.auc.polygon=TRUE, auc.polygon.col="skyblue", print.thres=TRUE)
    thres <- as.numeric(coords(modelroc, "best", best.method="closest.topleft", ret=c("threshold", "specificity", "sensitivity"))[1,1])
    # thres <- 0.5
    direction <- modelroc$direction
    auc <- modelroc$auc
    # 评估模型的预测效果
    if(direction=="<"){
      predict.results <- ifelse(lr_test>thres, 1, 0)
    }
    if(direction==">"){
      predict.results <- ifelse(lr_test<thres, 1, 0)
    }
    predict.results <- factor(predict.results, levels=c("0","1"))
    table <- table(predict.results, fold_test$是否发生非前哨淋巴结转移)
    print(table)
    tp <- table[2,2]
    fn <- table[1,2]
    fp <- table[2,1]
    tn <- table[1,1]
    acc <- (tp+tn)/(tp+fn+fp+tn)
    sen <- tp/(fn+tp)
    f0r <- fn/(fn+tn)
    spe <- tn/(fp+tn)
    fdr <- fp/(fp+tp)
    rec <- tp/(tp+fn)
    pre <- tp/(tp+fp)
    f1 <- 2*(rec*pre)/(rec+pre)
    print(paste("auc=", auc, sep=""))
    print(paste("acc=", acc, sep=""))
    print(paste("sen=", sen, sep=""))
    print(paste("for=", f0r, sep=""))
    print(paste("spe=", spe, sep=""))
    print(paste("fdr=", fdr, sep=""))
    print(paste("rec=", rec, sep=""))
    print(paste("pre=", pre, sep=""))
    print(paste("f1=", f1, sep=""))
    fold <- c()
    fold <- c(fold, auc, acc, sen, f0r, spe, fdr, rec, pre, f1)
    result <- rbind(result, fold)
  }
}
result <- result[-1,]
rownames(result) <- c(1:15)
colnames(result) <- c("AUC","ACC","SEN","FOR","SPE","FDR","REC","PRE","F1")
result <- as.data.frame(result)
result <- as.data.frame(lapply(result, function(x){paste(x)}))
write.csv(result, "in_tenon.csv", quote = FALSE)

###########################################
##Louisville##
###########################################

rm(list=ls())

setwd("E:/")
path="E:/"

data_input <- read.csv(paste(path, "EDC_229_35_num_non-sentinel_3_84.csv", sep=""), header=TRUE)

# 提取数据
# T分期+前哨淋巴结阳性个数分组+前哨淋巴结阳性比例分组
data <- data_input[,c(28,31,33,35)]
data$是否发生非前哨淋巴结转移 <- factor(data$是否发生非前哨淋巴结转移, levels=c(1,0), labels=c("yes","no"))
data$是否发生非前哨淋巴结转移 <- factor(data$是否发生非前哨淋巴结转移, levels=c("no","yes"), labels=c("no","yes"))

# Train and test sets
set.seed(3)
index_test_0 <- sample(which(data$是否发生非前哨淋巴结转移=="no"), 30)
set.seed(3)
index_test_1 <- sample(which(data$是否发生非前哨淋巴结转移=="yes"), 16)
index_test <- c(index_test_0,index_test_1)
data_train <- data[-index_test,]
data_test <- data[index_test,]
data_train[,c(1:4)] <- lapply(data_train[,c(1:4)], as.factor)
data_test[,c(1:4)] <- lapply(data_test[,c(1:4)], as.factor)

# external validation
library(pROC)
result <- matrix(ncol=9)
lr_train <- glm(是否发生非前哨淋巴结转移~T分期+前哨淋巴结阳性个数分组+前哨淋巴结阳性比例分组, family=binomial(link=logit), data=data_train)
lr_test <- predict(lr_train, type='response', newdata=data_test)
# 使用pROC包自动标出最优临界点
modelroc <- roc(data_test$是否发生非前哨淋巴结转移, lr_test, levels=c("no","yes"))
plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1,0.2), grid.col=c("green","red"), max.auc.polygon=TRUE, auc.polygon.col="skyblue", print.thres=TRUE)
thres <- as.numeric(coords(modelroc, "best", best.method="closest.topleft", ret=c("threshold", "specificity", "sensitivity"))[1,1])
# thres <- 0.5
direction <- modelroc$direction
auc <- modelroc$auc
# 评估模型的预测效果
if(direction=="<"){
  predict.results <- ifelse(lr_test>thres, 1, 0)
}
if(direction==">"){
  predict.results <- ifelse(lr_test<thres, 1, 0)
}
predict.results <- factor(predict.results, levels=c("0","1"))
table <- table(predict.results, data_test$是否发生非前哨淋巴结转移)
print(table)
tp <- table[2,2]
fn <- table[1,2]
fp <- table[2,1]
tn <- table[1,1]
acc <- (tp+tn)/(tp+fn+fp+tn)
sen <- tp/(fn+tp)
f0r <- fn/(fn+tn)
spe <- tn/(fp+tn)
fdr <- fp/(fp+tp)
rec <- tp/(tp+fn)
pre <- tp/(tp+fp)
f1 <- 2*(rec*pre)/(rec+pre)
print(paste("auc=", auc, sep=""))
print(paste("acc=", acc, sep=""))
print(paste("sen=", sen, sep=""))
print(paste("for=", f0r, sep=""))
print(paste("spe=", spe, sep=""))
print(paste("fdr=", fdr, sep=""))
print(paste("rec=", rec, sep=""))
print(paste("pre=", pre, sep=""))
print(paste("f1=", f1, sep=""))
model <- c(auc, acc, sen, f0r, spe, fdr, rec, pre, f1)
result <- rbind(result, model)
colnames(result) <- c("AUC","ACC","SEN","FOR","SPE","FDR","REC","PRE","F1")
write.csv(result, "ex_louisville.csv", quote=FALSE)

# cross validation
library(pROC)
library(caret)
data_all <- data_train
result <- matrix(ncol=9)
for(j in 99:101){
  set.seed(j)
  folds <- createFolds(data_all[,"是否发生非前哨淋巴结转移"], k=5)
  for(i in 1:5){
    fold_test <- data_all[folds[[i]],] # 取folds[[i]]作为测试集
    fold_train <- data_all[-folds[[i]],] # 剩下的数据作为训练集
    lr_train <- glm(是否发生非前哨淋巴结转移~T分期+前哨淋巴结阳性个数分组+前哨淋巴结阳性比例分组, family=binomial(link=logit), data=fold_train)
    lr_test <- predict(lr_train, type='response', newdata=fold_test)
    # 使用pROC包自动标出最优临界点
    modelroc <- roc(fold_test$是否发生非前哨淋巴结转移, lr_test, levels=c("no","yes"))
    plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1,0.2), grid.col=c("green","red"), max.auc.polygon=TRUE, auc.polygon.col="skyblue", print.thres=TRUE)
    thres <- as.numeric(coords(modelroc, "best", best.method="closest.topleft", ret=c("threshold", "specificity", "sensitivity"))[1,1])
    # thres <- 0.5
    direction <- modelroc$direction
    auc <- modelroc$auc
    # 评估模型的预测效果
    if(direction=="<"){
      predict.results <- ifelse(lr_test>thres, 1, 0)
    }
    if(direction==">"){
      predict.results <- ifelse(lr_test<thres, 1, 0)
    }
    predict.results <- factor(predict.results, levels=c("0","1"))
    table <- table(predict.results, fold_test$是否发生非前哨淋巴结转移)
    print(table)
    tp <- table[2,2]
    fn <- table[1,2]
    fp <- table[2,1]
    tn <- table[1,1]
    acc <- (tp+tn)/(tp+fn+fp+tn)
    sen <- tp/(fn+tp)
    f0r <- fn/(fn+tn)
    spe <- tn/(fp+tn)
    fdr <- fp/(fp+tp)
    rec <- tp/(tp+fn)
    pre <- tp/(tp+fp)
    f1 <- 2*(rec*pre)/(rec+pre)
    print(paste("auc=", auc, sep=""))
    print(paste("acc=", acc, sep=""))
    print(paste("sen=", sen, sep=""))
    print(paste("for=", f0r, sep=""))
    print(paste("spe=", spe, sep=""))
    print(paste("fdr=", fdr, sep=""))
    print(paste("rec=", rec, sep=""))
    print(paste("pre=", pre, sep=""))
    print(paste("f1=", f1, sep=""))
    fold <- c()
    fold <- c(fold, auc, acc, sen, f0r, spe, fdr, rec, pre, f1)
    result <- rbind(result, fold)
  }
}
result <- result[-1,]
rownames(result) <- c(1:15)
colnames(result) <- c("AUC","ACC","SEN","FOR","SPE","FDR","REC","PRE","F1")
result <- as.data.frame(result)
result <- as.data.frame(lapply(result, function(x){paste(x)}))
write.csv(result, "in_louisville.csv", quote = FALSE)

###########################################
##Stanford##
###########################################

rm(list=ls())

setwd("E:/")
path="E:/"

data_input <- read.csv(paste(path, "EDC_229_35_num_non-sentinel_3_84.csv", sep=""), header=TRUE)

# 提取数据
# 手术病理大小分组
data <- data_input[,c(22,35)]
data$是否发生非前哨淋巴结转移 <- factor(data$是否发生非前哨淋巴结转移, levels=c(1,0), labels=c("yes","no"))
data$是否发生非前哨淋巴结转移 <- factor(data$是否发生非前哨淋巴结转移, levels=c("no","yes"), labels=c("no","yes"))

# Train and test sets
set.seed(3)
index_test_0 <- sample(which(data$是否发生非前哨淋巴结转移=="no"), 30)
set.seed(3)
index_test_1 <- sample(which(data$是否发生非前哨淋巴结转移=="yes"), 16)
index_test <- c(index_test_0,index_test_1)
data_train <- data[-index_test,]
data_test <- data[index_test,]
data_train[,c(1:2)] <- lapply(data_train[,c(1:2)], as.factor)
data_test[,c(1:2)] <- lapply(data_test[,c(1:2)], as.factor)

# external validation
library(pROC)
result <- matrix(ncol=9)
lr_train <- glm(是否发生非前哨淋巴结转移~手术病理大小分组, family=binomial(link=logit), data=data_train)
lr_test <- predict(lr_train, type='response', newdata=data_test)
# 使用pROC包自动标出最优临界点
modelroc <- roc(data_test$是否发生非前哨淋巴结转移, lr_test, levels=c("no","yes"))
plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1,0.2), grid.col=c("green","red"), max.auc.polygon=TRUE, auc.polygon.col="skyblue", print.thres=TRUE)
thres <- as.numeric(coords(modelroc, "best", best.method="closest.topleft", ret=c("threshold", "specificity", "sensitivity"))[1,1])
# thres <- 0.5
direction <- modelroc$direction
auc <- modelroc$auc
# 评估模型的预测效果
if(direction=="<"){
  predict.results <- ifelse(lr_test>thres, 1, 0)
}
if(direction==">"){
  predict.results <- ifelse(lr_test<thres, 1, 0)
}
predict.results <- factor(predict.results, levels=c("0","1"))
table <- table(predict.results, data_test$是否发生非前哨淋巴结转移)
print(table)
tp <- table[2,2]
fn <- table[1,2]
fp <- table[2,1]
tn <- table[1,1]
acc <- (tp+tn)/(tp+fn+fp+tn)
sen <- tp/(fn+tp)
f0r <- fn/(fn+tn)
spe <- tn/(fp+tn)
fdr <- fp/(fp+tp)
rec <- tp/(tp+fn)
pre <- tp/(tp+fp)
f1 <- 2*(rec*pre)/(rec+pre)
print(paste("auc=", auc, sep=""))
print(paste("acc=", acc, sep=""))
print(paste("sen=", sen, sep=""))
print(paste("for=", f0r, sep=""))
print(paste("spe=", spe, sep=""))
print(paste("fdr=", fdr, sep=""))
print(paste("rec=", rec, sep=""))
print(paste("pre=", pre, sep=""))
print(paste("f1=", f1, sep=""))
model <- c(auc, acc, sen, f0r, spe, fdr, rec, pre, f1)
result <- rbind(result, model)
colnames(result) <- c("AUC","ACC","SEN","FOR","SPE","FDR","REC","PRE","F1")
write.csv(result, "ex_stanford.csv", quote=FALSE)

# cross validation
library(pROC)
library(caret)
data_all <- data_train
result <- matrix(ncol=9)
for(j in 99:101){
  set.seed(j)
  folds <- createFolds(data_all[,"是否发生非前哨淋巴结转移"], k=5)
  for(i in 1:5){
    fold_test <- data_all[folds[[i]],] # 取folds[[i]]作为测试集
    fold_train <- data_all[-folds[[i]],] # 剩下的数据作为训练集
    lr_train <- glm(是否发生非前哨淋巴结转移~手术病理大小分组, family=binomial(link=logit), data=fold_train)
    lr_test <- predict(lr_train, type='response', newdata=fold_test)
    # 使用pROC包自动标出最优临界点
    modelroc <- roc(fold_test$是否发生非前哨淋巴结转移, lr_test, levels=c("no","yes"))
    plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1,0.2), grid.col=c("green","red"), max.auc.polygon=TRUE, auc.polygon.col="skyblue", print.thres=TRUE)
    thres <- as.numeric(coords(modelroc, "best", best.method="closest.topleft", ret=c("threshold", "specificity", "sensitivity"))[1,1])
    # thres <- 0.5
    direction <- modelroc$direction
    auc <- modelroc$auc
    # 评估模型的预测效果
    if(direction=="<"){
      predict.results <- ifelse(lr_test>thres, 1, 0)
    }
    if(direction==">"){
      predict.results <- ifelse(lr_test<thres, 1, 0)
    }
    predict.results <- factor(predict.results, levels=c("0","1"))
    table <- table(predict.results, fold_test$是否发生非前哨淋巴结转移)
    print(table)
    tp <- table[2,2]
    fn <- table[1,2]
    fp <- table[2,1]
    tn <- table[1,1]
    acc <- (tp+tn)/(tp+fn+fp+tn)
    sen <- tp/(fn+tp)
    f0r <- fn/(fn+tn)
    spe <- tn/(fp+tn)
    fdr <- fp/(fp+tp)
    rec <- tp/(tp+fn)
    pre <- tp/(tp+fp)
    f1 <- 2*(rec*pre)/(rec+pre)
    print(paste("auc=", auc, sep=""))
    print(paste("acc=", acc, sep=""))
    print(paste("sen=", sen, sep=""))
    print(paste("for=", f0r, sep=""))
    print(paste("spe=", spe, sep=""))
    print(paste("fdr=", fdr, sep=""))
    print(paste("rec=", rec, sep=""))
    print(paste("pre=", pre, sep=""))
    print(paste("f1=", f1, sep=""))
    fold <- c()
    fold <- c(fold, auc, acc, sen, f0r, spe, fdr, rec, pre, f1)
    result <- rbind(result, fold)
  }
}
result <- result[-1,]
rownames(result) <- c(1:15)
colnames(result) <- c("AUC","ACC","SEN","FOR","SPE","FDR","REC","PRE","F1")
result <- as.data.frame(result)
result <- as.data.frame(lapply(result, function(x){paste(x)}))
write.csv(result, "in_stanford.csv", quote = FALSE)

###########################################
##Mayo##
###########################################

rm(list=ls())

setwd("E:/")
path="E:/"

data_input <- read.csv(paste(path, "EDC_229_35_num_non-sentinel_3_84.csv", sep=""), header=TRUE)

# 提取数据
# 年龄分组+手术病理大小分组+ER状态+前哨淋巴结阳性个数分组+前哨淋巴结阴性个数分组
data <- data_input[,c(4,22,24,31,32,35)]
data$是否发生非前哨淋巴结转移 <- factor(data$是否发生非前哨淋巴结转移, levels=c(1,0), labels=c("yes","no"))
data$是否发生非前哨淋巴结转移 <- factor(data$是否发生非前哨淋巴结转移, levels=c("no","yes"), labels=c("no","yes"))

# Train and test sets
set.seed(3)
index_test_0 <- sample(which(data$是否发生非前哨淋巴结转移=="no"), 30)
set.seed(3)
index_test_1 <- sample(which(data$是否发生非前哨淋巴结转移=="yes"), 16)
index_test <- c(index_test_0,index_test_1)
data_train <- data[-index_test,]
data_test <- data[index_test,]
data_train[,c(1:6)] <- lapply(data_train[,c(1:6)], as.factor)
data_test[,c(1:6)] <- lapply(data_test[,c(1:6)], as.factor)

# external validation
library(pROC)
result <- matrix(ncol=9)
lr_train <- glm(是否发生非前哨淋巴结转移~年龄分组+手术病理大小分组+ER状态+前哨淋巴结阳性个数分组+前哨淋巴结阴性个数分组, family=binomial(link=logit), data=data_train)
lr_test <- predict(lr_train, type='response', newdata=data_test)
# 使用pROC包自动标出最优临界点
modelroc <- roc(data_test$是否发生非前哨淋巴结转移, lr_test, levels=c("no","yes"))
plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1,0.2), grid.col=c("green","red"), max.auc.polygon=TRUE, auc.polygon.col="skyblue", print.thres=TRUE)
thres <- as.numeric(coords(modelroc, "best", best.method="closest.topleft", ret=c("threshold", "specificity", "sensitivity"))[1,1])
# thres <- 0.5
direction <- modelroc$direction
auc <- modelroc$auc
# 评估模型的预测效果
if(direction=="<"){
  predict.results <- ifelse(lr_test>thres, 1, 0)
}
if(direction==">"){
  predict.results <- ifelse(lr_test<thres, 1, 0)
}
predict.results <- factor(predict.results, levels=c("0","1"))
table <- table(predict.results, data_test$是否发生非前哨淋巴结转移)
print(table)
tp <- table[2,2]
fn <- table[1,2]
fp <- table[2,1]
tn <- table[1,1]
acc <- (tp+tn)/(tp+fn+fp+tn)
sen <- tp/(fn+tp)
f0r <- fn/(fn+tn)
spe <- tn/(fp+tn)
fdr <- fp/(fp+tp)
rec <- tp/(tp+fn)
pre <- tp/(tp+fp)
f1 <- 2*(rec*pre)/(rec+pre)
print(paste("auc=", auc, sep=""))
print(paste("acc=", acc, sep=""))
print(paste("sen=", sen, sep=""))
print(paste("for=", f0r, sep=""))
print(paste("spe=", spe, sep=""))
print(paste("fdr=", fdr, sep=""))
print(paste("rec=", rec, sep=""))
print(paste("pre=", pre, sep=""))
print(paste("f1=", f1, sep=""))
model <- c(auc, acc, sen, f0r, spe, fdr, rec, pre, f1)
result <- rbind(result, model)
colnames(result) <- c("AUC","ACC","SEN","FOR","SPE","FDR","REC","PRE","F1")
write.csv(result, "ex_mayo.csv", quote=FALSE)

# cross validation
library(pROC)
library(caret)
data_all <- data_train
result <- matrix(ncol=9)
for(j in 99:101){
  set.seed(j)
  folds <- createFolds(data_all[,"是否发生非前哨淋巴结转移"], k=5)
  for(i in 1:5){
    fold_test <- data_all[folds[[i]],] # 取folds[[i]]作为测试集
    fold_train <- data_all[-folds[[i]],] # 剩下的数据作为训练集
    lr_train <- glm(是否发生非前哨淋巴结转移~年龄分组+手术病理大小分组+ER状态+前哨淋巴结阳性个数分组+前哨淋巴结阴性个数分组, family=binomial(link=logit), data=fold_train)
    lr_test <- predict(lr_train, type='response', newdata=fold_test)
    # 使用pROC包自动标出最优临界点
    modelroc <- roc(fold_test$是否发生非前哨淋巴结转移, lr_test, levels=c("no","yes"))
    plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1,0.2), grid.col=c("green","red"), max.auc.polygon=TRUE, auc.polygon.col="skyblue", print.thres=TRUE)
    thres <- as.numeric(coords(modelroc, "best", best.method="closest.topleft", ret=c("threshold", "specificity", "sensitivity"))[1,1])
    # thres <- 0.5
    direction <- modelroc$direction
    auc <- modelroc$auc
    # 评估模型的预测效果
    if(direction=="<"){
      predict.results <- ifelse(lr_test>thres, 1, 0)
    }
    if(direction==">"){
      predict.results <- ifelse(lr_test<thres, 1, 0)
    }
    predict.results <- factor(predict.results, levels=c("0","1"))
    table <- table(predict.results, fold_test$是否发生非前哨淋巴结转移)
    print(table)
    tp <- table[2,2]
    fn <- table[1,2]
    fp <- table[2,1]
    tn <- table[1,1]
    acc <- (tp+tn)/(tp+fn+fp+tn)
    sen <- tp/(fn+tp)
    f0r <- fn/(fn+tn)
    spe <- tn/(fp+tn)
    fdr <- fp/(fp+tp)
    rec <- tp/(tp+fn)
    pre <- tp/(tp+fp)
    f1 <- 2*(rec*pre)/(rec+pre)
    print(paste("auc=", auc, sep=""))
    print(paste("acc=", acc, sep=""))
    print(paste("sen=", sen, sep=""))
    print(paste("for=", f0r, sep=""))
    print(paste("spe=", spe, sep=""))
    print(paste("fdr=", fdr, sep=""))
    print(paste("rec=", rec, sep=""))
    print(paste("pre=", pre, sep=""))
    print(paste("f1=", f1, sep=""))
    fold <- c()
    fold <- c(fold, auc, acc, sen, f0r, spe, fdr, rec, pre, f1)
    result <- rbind(result, fold)
  }
}
result <- result[-1,]
rownames(result) <- c(1:15)
colnames(result) <- c("AUC","ACC","SEN","FOR","SPE","FDR","REC","PRE","F1")
result <- as.data.frame(result)
result <- as.data.frame(lapply(result, function(x){paste(x)}))
write.csv(result, "in_mayo.csv", quote = FALSE)

###########################################
##Between-Models##
###########################################

# 如若报错可能需要重启R
rm(list=ls())

setwd("E:/")
path="E:/"

data_input <- read.csv(paste(path, "seed3/in_prediction_non-sentinel.csv", sep=""), header=TRUE)
data <- data_input[is.na(data_input$Value)==FALSE,]
data$Model <- factor(data$Model, levels=c("clinical+genotype(29 SNPs)","MSKCC","MOU","Ljubljana","SNUH","Mayo","Tenon","Louisville","Cambridge","Stanford","MDA"))
data$Index <- factor(data$Index, levels=c("AUC","Accuracy","Sensitivity","FOR","Specificity","FDR","F1 score"))

library(ggsci)
library(scales)
library(ggplot2)
p <- ggplot(data) + 
  geom_boxplot(aes(x=Index, y=Value, fill=Model), width=0.6, position=position_dodge(0.8)) + 
  scale_fill_manual(values=pal_igv("default")(12)[2:12], 
                    breaks=c("clinical+genotype(29 SNPs)","MSKCC","MOU","Ljubljana","SNUH","Mayo","Tenon","Louisville","Cambridge","Stanford","MDA"), 
                    labels=c("clinical+genotype(29 SNPs)","MSKCC","MOU","Ljubljana","SNUH","Mayo","Tenon","Louisville","Cambridge","Stanford","MDA")) + 
  xlab("") + 
  ylab("") + 
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.05)) + 
  theme_bw() + 
  theme(
    legend.position="top",
    legend.background=element_blank(),
    legend.key = element_blank(),
    legend.margin=margin(0,0,0,0,"mm"),
    axis.text.x=element_text(size=rel(2),colour="black"),
    axis.text.y=element_text(size=rel(2),colour="black"),
    axis.line.x=element_line(size=0.5,colour="black"),
    axis.line.y=element_line(size=0.5,colour="black"),
    legend.text=element_text(size=rel(1.8)),
    legend.title=element_blank(),
    panel.border=element_blank(),
    panel.grid=element_blank(),
    panel.background=element_rect(fill="transparent",color=NA),
    plot.background=element_rect(fill="transparent",color=NA)
  ) + 
  guides(color="none")
png('in_prediction_non-sentinel.png', width=1200, height=600, bg="transparent")
print(p)
dev.off()

###########################################
##Between-Models##SNP##
###########################################

# 如若报错可能需要重启R
rm(list=ls())

setwd("E:/")
path="E:/"

data_input <- read.csv(paste(path, "seed3/in_prediction_non-sentinel_snp.csv", sep=""), header=TRUE)
data <- data_input[is.na(data_input$Value)==FALSE,]
data$Model <- factor(data$Model, levels=c("clinical+genotype","clinical+genotype(29 SNPs)","clinical+genotype(20 SNPs)","clinical+genotype(15 SNPs)","clinical"))
data$Index <- factor(data$Index, levels=c("AUC","Accuracy","Sensitivity","FOR","Specificity","FDR","F1 score"))

library(ggsci)
library(scales)
library(ggplot2)
p <- ggplot(data) + 
  geom_boxplot(aes(x=Index, y=Value, fill=Model), width=0.6, position=position_dodge(0.8)) + 
  scale_fill_manual(values=pal_igv("default")(12), 
                    breaks=c("clinical+genotype","clinical+genotype(29 SNPs)","clinical+genotype(20 SNPs)","clinical+genotype(15 SNPs)","clinical"), 
                    labels=c("clinical+genotype","clinical+genotype(29 SNPs)","clinical+genotype(20 SNPs)","clinical+genotype(15 SNPs)","clinical")) + 
  xlab("") + 
  ylab("") + 
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.05)) + 
  theme_bw() + 
  theme(
    legend.position="top",
    legend.background=element_blank(),
    legend.key = element_blank(),
    legend.margin=margin(0,0,0,0,"mm"),
    axis.text.x=element_text(size=rel(2),colour="black"),
    axis.text.y=element_text(size=rel(2),colour="black"),
    axis.line.x=element_line(size=0.5,colour="black"),
    axis.line.y=element_line(size=0.5,colour="black"),
    legend.text=element_text(size=rel(1.8)),
    legend.title=element_blank(),
    panel.border=element_blank(),
    panel.grid=element_blank(),
    panel.background=element_rect(fill="transparent",color=NA),
    plot.background=element_rect(fill="transparent",color=NA)
  ) + 
  guides(color="none")
png('in_prediction_non-sentinel_snp.png', width=1200, height=600, bg="transparent")
print(p)
dev.off()

###########################################
##nomogram##
###########################################

rm(list=ls())

setwd("E:/")
path="E:/"

data_input <- read.csv(paste(path, "EDC_229_35_num_non-sentinel_3_29.csv", sep=""), header=TRUE)
# cluster2+cluster1+前哨淋巴结阳性个数分组+前哨淋巴结阴性个数分组
data <- data_input[,c(3,2,31,32,35)]
colnames(data) <- c("cluster2","cluster1","number_of_positive_SLNs","number_of_negative_SLNs","NSLN_metastasis")

# Train and test sets
set.seed(3)
index_test_0 <- sample(which(data$NSLN_metastasis==0), 30)
set.seed(3)
index_test_1 <- sample(which(data$NSLN_metastasis==1), 16)
index_test <- c(index_test_0,index_test_1)
data_train <- data[-index_test,]
data_test <- data[index_test,]
data_train[,c(3:5)] <- lapply(data_train[,c(3:5)], as.factor)
data_test[,c(3:5)] <- lapply(data_test[,c(3:5)], as.factor)

# 设置因子变量（label按照分组顺序）
data_train$number_of_positive_SLNs <- factor(data_train$number_of_positive_SLNs, labels=c('1','2','≥3'))
relevel(data_train$number_of_positive_SLNs, ref='1')
data_train$number_of_negative_SLNs <- factor(data_train$number_of_negative_SLNs, labels=c('0','[1,2]','≥3'))
relevel(data_train$number_of_negative_SLNs, ref='0')

# 绘制nomogram
library(rms)
# 打包数据
dd <- datadist(data_train)
options(datadist="dd")
# 逻辑回归
model_logit <- lrm(NSLN_metastasis~., data=data_train, x=T, y=T)
model_logit
# 绘制LASSO回归的风险预测值的nomogram
nom <- nomogram(model_logit, fun=function(x)1/(1+exp(-x)), # or fun=plogis
                fun.at=c(0.001,0.01,0.05,seq(0.1,0.9,by=0.1),0.95,0.99,0.999),
                lp=F, funlabel="Risk of NSLN metastasis")
plot(nom)

###########################################
##AF检验##
###########################################

rm(list=ls())
library(ggpubr)

setwd("E:/")
path="E:/"

data_input <- read.csv(paste(path, "SNP_AF.csv", sep=""), header=TRUE)

# Wilcoxon rank-sum test
ggbarplot(data_input, x="cluster", y="BC", width=0.5, size=1, 
               add=c("mean_se","jitter"), shape="cluster", 
               color="cluster", palette=c("#5050FFFF","#CE3D32FF")) + 
  xlab("") + 
  ylab("Allele Frequency") + 
  # ggtitle("Allele Frequency") + 
  theme(
    legend.position="none",
    plot.title = element_text(size=rel(1.4),hjust=0.5),
    axis.text.x=element_text(size=rel(1.4),colour="black"),
    axis.text.y=element_text(size=rel(1.2),colour="black"),
    axis.line.x=element_line(size=1,colour="black"),
    axis.line.y=element_line(size=1,colour="black"),
    axis.title.y=element_text(size=rel(1.4),color="black")
  ) + 
  scale_y_continuous(limits=c(0,1.05), breaks=seq(0,1,0.25)) + 
  stat_compare_means(aes(label=paste0("p ", ..p.format.., " (", ..p.signif.., ")")), 
                     label.y=1, size=5)

# Wilcoxon signed-rank test
# data <- data_input[is.na(data_input$gnomAD_EAS)==FALSE & is.na(data_input$X1000g2015aug_eas)==FALSE,]
data <- data_input
data[is.na(data)] <- 0
data_1 <- data.frame(data[,c(1,2)], database="gnomAD_EAS", AF=data[,c(4)])
data_2 <- data.frame(data[,c(1,2)], database="1000g_EAS", AF=data[,c(5)])
data_3 <- data.frame(data[,c(1,2)], database="BC_cohort", AF=data[,c(8)])
data_db <- rbind(data_3, data_2, data_1)

# cluster
# ggpaired(data_db[c(1:150),], x="database", y="AF", 
ggpaired(data_db[c(1:192),], x="database", y="AF", 
              color="cluster", palette=c("#5050FFFF","#CE3D32FF"), 
              line.color="gray", line.size=0.4, 
              facet.by="cluster", short.panel.labs=TRUE) + 
  xlab("") + 
  ylab("Allele Frequency") + 
  # ggtitle("Allele Frequency") + 
  theme(
    legend.position="none",
    plot.title = element_text(size=rel(1),hjust=0.5),
    axis.text.x=element_text(size=rel(0.8),colour="black",face="bold"),
    axis.text.y=element_text(size=rel(0.8),colour="black",face="bold"),
    axis.title.y=element_text(size=rel(0.8),color="black",face="bold"),
    strip.text.x=element_text(size=rel(1.2),color="black",face="bold"),
  ) + 
  scale_y_continuous(limits=c(0,1.05), breaks=seq(0,1,0.25)) + 
  stat_compare_means(aes(label=paste0("p = ", ..p.format.., " (", ..p.signif.., ")")), 
                     label.y=1, size=3.5, fontface=2, paired=TRUE)

# ggpaired(data_db[c(1:75,151:225),], x="database", y="AF", 
ggpaired(data_db[c(1:96,193:288),], x="database", y="AF", 
              color="cluster", palette=c("#5050FFFF","#CE3D32FF"), 
              line.color="gray", line.size = 0.4, 
              facet.by="cluster", short.panel.labs=TRUE) + 
  xlab("") + 
  ylab("Allele Frequency") + 
  # ggtitle("Allele Frequency") + 
  theme(
    legend.position="none",
    plot.title = element_text(size=rel(1),hjust=0.5),
    axis.text.x=element_text(size=rel(0.8),colour="black",face="bold"),
    axis.text.y=element_text(size=rel(0.8),colour="black",face="bold"),
    axis.title.y=element_text(size=rel(0.8),color="black",face="bold"),
    strip.text.x=element_text(size=rel(1.2),color="black",face="bold"),
  ) + 
  scale_y_continuous(limits=c(0,1.05), breaks=seq(0,1,0.25)) + 
  stat_compare_means(aes(label=paste0("p = ", ..p.format.., " (", ..p.signif.., ")")), 
                     label.y=1, size=3.5, fontface=2, paired=TRUE)

# ggpaired(data_db[c(76:225),], x="database", y="AF", 
ggpaired(data_db[c(97:288),], x="database", y="AF", 
              color="cluster", palette=c("#5050FFFF","#CE3D32FF"), 
              line.color="gray", line.size = 0.4, 
              facet.by="cluster", short.panel.labs=TRUE) + 
  xlab("") + 
  ylab("Allele Frequency") + 
  # ggtitle("Allele Frequency") + 
  theme(
    legend.position="none",
    plot.title = element_text(size=rel(1),hjust=0.5),
    axis.text.x=element_text(size=rel(0.8),colour="black",face="bold"),
    axis.text.y=element_text(size=rel(0.8),colour="black",face="bold"),
    axis.title.y=element_text(size=rel(0.8),color="black",face="bold"),
    strip.text.x=element_text(size=rel(1.2),color="black",face="bold"),
  ) + 
  scale_y_continuous(limits=c(0,1.05), breaks=seq(0,1,0.25)) + 
  stat_compare_means(aes(label=paste0("p = ", ..p.format.., " (", ..p.signif.., ")")), 
                     label.y=1, size=3.5, fontface=2, paired=TRUE)

# all
# ggpaired(data_db[c(1:150),], x="database", y="AF", 
ggpaired(data_db[c(1:192),], x="database", y="AF", 
         color="cluster", palette=c("#5050FFFF","#CE3D32FF"), 
         line.color="gray", line.size=0.4) + 
  stat_compare_means(paired=TRUE)

# ggpaired(data_db[c(1:75,151:225),], x="database", y="AF", 
ggpaired(data_db[c(1:96,193:288),], x="database", y="AF", 
         color="cluster", palette=c("#5050FFFF","#CE3D32FF"), 
         line.color="gray", line.size=0.4) + 
  stat_compare_means(paired=TRUE)

# ggpaired(data_db[c(76:225),], x="database", y="AF", 
ggpaired(data_db[c(97:288),], x="database", y="AF", 
         color="cluster", palette=c("#5050FFFF","#CE3D32FF"), 
         line.color="gray", line.size=0.4) + 
  stat_compare_means(paired=TRUE)

#######################################
# 基因功能分析
#######################################

rm(list=ls())
setwd("E:/")
path="E:/"

library(data.table)
data_input <- fread(paste(path, "EDC_1773_96_num_km_ds_gene.txt", sep=""), sep="\t")
cluster1 <- unique(data_input[data_input$ds_cluster==1,]$gene)
cluster2 <- unique(data_input[data_input$ds_cluster==2,]$gene)

# 转ID
library(clusterProfiler)
library(org.Hs.eg.db)
eg1 <- bitr(cluster1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
eg2 <- bitr(cluster2, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
write.table(eg1, "cluster1.txt", sep="\t", row.names=F, col.names=F, quote=F)
write.table(eg2, "cluster2.txt", sep="\t", row.names=F, col.names=F, quote=F)
eg <- eg1
eg <- eg2

# KEGG分析
ekegg <- enrichKEGG(gene          = eg$ENTREZID,
                    organism      = 'hsa',
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05)
# head(ekegg)
dim(ekegg)
# barplot(ekegg, showCategory=20)
dotplot(ekegg, showCategory=20, font.size=15)
# emapplot(ekegg)
# upsetplot(ekegg)
browseKEGG(ekegg, ekegg@result$ID[10])

###########################################
##特征统计##
###########################################

rm(list=ls())

setwd("E:/")
path="E:/"

data_input <- read.csv(paste(path, "EDC_229_35_num_non-sentinel_3_84.csv", sep=""), header=TRUE)
data <- data_input[,c(4:6,9:10,12,15,18:20,22,24:27,29:33,35)]

# Train and test sets
set.seed(3)
index_test_0 <- sample(which(data$是否发生非前哨淋巴结转移==0), 30)
set.seed(3)
index_test_1 <- sample(which(data$是否发生非前哨淋巴结转移==1), 16)
index_test <- c(index_test_0,index_test_1)
data_train <- data[-index_test,]
data_test <- data[index_test,]

library(data.table)
data_all <- data_train
# which(data_all$体格检查淋巴结个数==2, arr.ind=TRUE)
# data_all <- data_all[-3,]
data_all_stat <- data.table("Category"="","0"="", "1"="")
data_all_stat <- data_all_stat[-1,]
for(i in 1:(ncol(data_all)-1)){
  temp <- as.data.frame.matrix(table(data_all[,c(i,21)]))
  temp <- cbind(Category=rownames(temp), temp)
  data_all_stat <- rbind(data_all_stat, temp)
}
write.table(data_all_stat, "statistics_feature_non-sentinel.txt", sep="\t", row.names=F, quote=F)
