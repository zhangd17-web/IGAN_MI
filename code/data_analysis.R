library(caret)
library(e1071)
library(randomForest)
####################
####data loading####
####################
library(readxl)
new_data <- read_xlsx(path = "GWAS.0314.xlsx",sheet = 1,na = "NA",
                      col_types = c(rep("guess",69),"numeric",rep("guess",12)))
imputation_TG = read_xlsx("IgA-2022.3.13.xlsx",sheet = 1,col_names = TRUE,na = "NA")

imputation_TCHO = read_xlsx("IgA-2022.3.13.xlsx",sheet = 2,col_names = TRUE,na="NA")

imputation_LDL = read_xlsx("IgA-2022.3.13.xlsx",sheet = 3,col_names = TRUE,na="NA")

imputation_IgA = read_xlsx("IgA-2022.3.13.xlsx",sheet = 4,col_names = TRUE,na="NA")

new_data$`TG mmol/L`[match(imputation_TG$病理号8位,new_data$病理号8位)] = imputation_TG$TG
new_data$`TCHO mmol/L`[match(imputation_TCHO$病理号8位,new_data$病理号8位)] = imputation_TCHO$TCHO
new_data$`LDL mmol/L`[match(imputation_LDL$病理号8位,new_data$病理号8位)] = imputation_LDL$LDL
new_data$`IgA g/L`[match(imputation_IgA$病理号8位,new_data$病理号8位)] = imputation_IgA$IgA
new_data_withT <- filter(new_data,!is.na(T))
new_data$composite<- 0
for (i in 1:nrow(new_data)){
  if ((new_data$END7[i]==1) || (new_data$END8[i]==1) || (new_data$RRT[i]==1) ){
    new_data$composite[i]<- 1
  }}
new_data$composite_time <- 0 
for (i in 1:nrow(new_data)){
  new_data$composite_time[i]<-min(new_data$END7T[i],new_data$END8T[i],new_data$RRTT[i],na.rm = TRUE)
}

data_follow <-read.csv("origin_data_follow.csv")
id_forT <- data_follow$数据库DNA号
id_forT <- clean_data$FID[!(clean_data$FID%in%combined_data$FID)]
test_data <- filter(new_data,(数据库DNA号%in%id_forT)&T!="NA")%>%
  mutate(uacid=0)
test_data <- test_data%>%mutate(uacid=0)
test_data$uacid[which(test_data$GENDER==1&test_data$`尿酸(umol/L)`>420)]=1
test_data$uacid[which(test_data$GENDER==2&test_data$`尿酸(umol/L)`>360)]=1

test_T<-test_data%>%select(数据库DNA号,T,age_Y,sbp,dbp,`GFR-baseline cMDRD`,`IgA g/L`,`baseline UTP(g/d)`,`尿酸(umol/L)`,uacid)
test_T <- test_T[,-c(1,9)] 
test_T <- mutate(test_T,biT = ifelse(T!=0,1,0))

test_T <- test_T[,-1]
colnames(test_T) <- c("age","SBP","DBP","eGFR","IgA...g.L.","upro","uacid","biT")

for (i in c(1,2,3,4,6)) {
  test_T[is.na(test_T[,i]),i] = mean(as.matrix(test_T[,i]),na.rm=TRUE) 
}

test_T <- test_T%>%tidyr::drop_na()
training_dat.T <- read.csv("training_dat.T.csv")
training_dat.T = WPRS_data.1[which(!WPRS_data$FID%in%follow_data$FID),]
training_dat.T$biT=as.factor(training_dat.T$biT)
##########################
###training model for T###
##########################

set.seed(35726)
cross.index = createFolds(1:nrow(training_dat.T),k=10)
WRF_score = rep(0,nrow(training_dat.T))
WSVM_score = rep(0,nrow(training_dat.T))
GLM_score = rep(0,nrow(training_dat.T))
P_weights = rep(1,nrow(training_dat.T))
P_weights[which(training_dat.T$biT==1)]=2.5

for (i in 1:10) {
  k =  cross.index[[i]]
  modsvm.cv = svm(as.factor(biT)~.,data= training_dat.T[setdiff(1:nrow(training_dat.T),k),],class.weights=c("0"=1,"1"=2.5),probability=TRUE)
  predsvm = attr(predict(modsvm.cv,training_dat.T[k,],probability = TRUE),"probabilities")[,2]
  modrf.cv = randomForest(as.factor(biT)~., classwt = c(1,2.5) ,data = training_dat.T[setdiff(1:nrow(training_dat.T),k),])
  predrf = predict(modrf.cv,training_dat.T[k,],type="prob")[,2]
  modglm.cv = glm(biT~.,family = "binomial",weights =P_weights[setdiff(1:nrow(training_dat.T),k)],data = training_dat.T[setdiff(1:nrow(training_dat.T),k),])
  predglm = predict(modglm.cv,training_dat.T[k,],type="response")
  WRF_score[k] = predrf
  WSVM_score[k] = predsvm
  GLM_score[k] = predglm
  
}
####################################
####second layer model training#####
####################################
train_data2 = data.frame(biT=training_dat.T$biT,WRF=WRF_score,WSVM =WSVM_score,GLM=GLM_score)
sec_model = glm(biT~.,data=train_data2,family = "binomial")
svm_train = svm(as.factor(biT)~.,data=training_dat.T,class.weights=c("0"=1,"1"=2.5),probability=TRUE)
rf_train = randomForest(as.factor(biT)~.,classwt = c(1,2.5),data = training_dat.T)
glm_train = glm(biT~.,family="binomial",weights = P_weights,data = training_dat.T)
########################

svm_pred = attr(predict(svm_train,newdata = test_T,probability = TRUE),"probabilities")[,2]
rf_pred = predict(rf_train,test_T,type="prob")[,2]
glm_pred = predict(glm_train,test_T,type="response")
test_data.sec = data.frame(biT=test_T$biT,WRF = rf_pred,WSVM =svm_pred,GLM= glm_pred)

test_fit = predict(sec_model,test_data.sec,type="response")
library(pROC)
roc1 = roc(test_T$biT~test_fit)
print(roc1)
ci.auc(roc1)
coords(roc1,"best",transpose=FALSE)
p1<-ggroc(roc1,legacy.axes = TRUE,color="blue")+geom_text(aes(x=0.3,y=0.7,label="AUC=0.823"),size=5,check_overlap = TRUE)+theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+geom_abline(intercept = 0,slope=1)
ggsave(filename = "T_prediction.png",plot = p1,device = "png",dpi = 320)

####################################################################
##############final predicting model training###############
####################################################################

train_temp <- read.csv("train_temp.csv")
train_temp <- mutate(train_temp,biT=ifelse(train_temp$T>0,1,0))
test_temp <- read.csv("test_temp.csv")
test_temp<- mutate(test_temp,biT=ifelse(test_temp$T>0,1,0))
mod_temp <- randomForest::randomForest(as.factor(composite)~.,data=train_temp,mtry=4,ntree=200)
temp_rT[id_temp] <- predict(mod_temp,test_temp,type="prob")[,2]
print(pROC::roc(train_dat_w$composite[id_temp]~predict(mod_temp,test_temp,type="prob")[,2]))

coords(pROC::roc(train_dat_w$composite[id_temp]~predict(mod_temp,test_temp,type="prob")[,2]),"best",transpose=FALSE)
plot.roc(pROC::roc(train_dat_w$composite[id_temp]~predict(mod_temp,test_temp,type="prob")[,2]),print.auc=TRUE, auc.polygon = TRUE,grid=c(0.1,0.2),grid.col=c("green","red"),auc.polygon.col="lightblue",print.thres=TRUE)
temp = predict(mod_temp,test_temp,type="prob")[,2]
sum(test_temp$composite==(temp>0.087))
ci.auc(roc(train_dat_w$composite[id_temp]~predict(mod_temp,test_temp,type="prob")[,2]))
glm_cv <- glmnet::cv.glmnet(x = as.matrix(train_temp[,-12]),y=as.matrix(train_temp$composite),family = "binomial")
mod_logi <- glmnet::glmnet(x = train_temp[,-12],y=train_temp$composite,family = "binomial",lambda = 0.2,alpha = 0)
mod_ANN <- ANN2::neuralnetwork(X = as.matrix(train_temp[,-12]),y=as.matrix(train_temp$composite),hidden.layers = c(100,30,30),optim.type = 'adam',learn.rates = 1e-04,val.prop = 0,n.epochs = 35,regression = FALSE,activ.functions = 'relu')
print(roc(test_temp$composite~predict(mod_ANN,newdata = test_temp[,-12])$probabilities[,2]))
temp = predict(mod_ANN,newdata = test_temp[,-12])$probabilities[,2]
sum(test_temp$composite==(temp>0.064))
coords(roc(test_temp$composite~predict(mod_ANN,newdata = test_temp[,-12])$probabilities[,2]),"best",transpose=FALSE)
print(roc(test_temp$composite~predict(mod_logi,as.matrix(test_temp[,-12]),type="response"),ci=TRUE))
coords(roc(test_temp$composite~predict(mod_logi,as.matrix(test_temp[,-12]),type="response")),"best",transpose=FALSE)
temp = predict(mod_logi,as.matrix(test_temp[,-12]),type="response")
sum(test_temp$composite==(temp>0.068))

library(xgboost)
test_temp <- train_dat_w[id_temp,]
test_temp<- mutate(test_temp,biT=ifelse(test_temp$T>0,1,0))
X <-as.matrix(train_temp[,-c(12,13,16:19)])
dtrain <- xgb.DMatrix(X,label=as.numeric(train_temp$composite))
param <- list(max_depth=5, eta = 0.05, objective = "binary:logistic")
watchlist <- list(train = dtrain)
mod_xgb <- xgb.train(param, dtrain, nrounds=110, watchlist,verbose = FALSE)
print(roc(test_temp$composite~predict(mod_xgb,as.matrix(test_temp[,-c(12,13,16:19)]))))
a.1 <- predict(mod_xgb,as.matrix(test_temp[,-c(12,13,16:19)]))
coords(roc(test_temp$composite~predict(mod_xgb,as.matrix(test_temp[,-c(12,13,16:19)]))),"best",transpose = FALSE)
roc1 = roc(test_temp$composite~predict(mod_xgb,as.matrix(test_temp[,-c(12,13,16:19)])))
plot.roc(roc1,print.auc=TRUE, auc.polygon = FALSE,print.thres=TRUE)
ci.auc(roc1)


X <-as.matrix(train_temp[,-c(12,13,16:20)])
dtrain <- xgb.DMatrix(X,label=as.numeric(train_temp$composite))
param <- list(max_depth=5, eta = 0.05, objective = "binary:logistic")
watchlist <- list(train = dtrain)
mod_xgb <- xgb.train(param, dtrain, nrounds=110, watchlist,verbose = FALSE)
print(roc(test_temp$composite~predict(mod_xgb,as.matrix(test_temp[,-c(12,13,16:20)]))))
a.2 <- predict(mod_xgb,as.matrix(test_temp[,-c(12,13,16:20)]))
coords(roc(test_temp$composite~predict(mod_xgb,as.matrix(test_temp[,-c(12,13,16:20)]))),"best",transpose = FALSE)
roc2 = roc(test_temp$composite~predict(mod_xgb,as.matrix(test_temp[,-c(12,13,16:20)])))

test_temp <- train_dat_w_pT[id_temp,]
X <-as.matrix(train_temp[,-c(12,13,16:19)])
dtrain <- xgb.DMatrix(X,label=as.numeric(train_temp$composite))
param <- list(max_depth=5, eta = 0.05, objective = "binary:logistic")
watchlist <- list(train = dtrain)
mod_xgb <- xgb.train(param, dtrain, nrounds=110, watchlist,verbose = FALSE)
print(roc(test_temp$composite~predict(mod_xgb,as.matrix(test_temp[,-c(12,13,16:19)]))))
a.3<-predict(mod_xgb,as.matrix(test_temp[,-c(12,13,16:19)]))
coords(roc(test_temp$composite~predict(mod_xgb,as.matrix(test_temp[,-c(12,13,16:19)]))),"best",transpose = FALSE)
roc3 = roc(test_temp$composite~predict(mod_xgb,as.matrix(test_temp[,-c(12,13,16:19)])))

roc4 = roc(test_temp$composite~a.2+a.1+a.3)
roc4 %>% 
  map(~tibble(AUC = .x$auc)) %>% 
  bind_rows(.id = "name") -> data.auc
data.auc %>% 
  mutate(label_long=paste0(name," , AUC = ",paste(round(AUC,2))),
         label_AUC=paste0("AUC = ",paste(round(AUC,2)))) -> data.labels
data.labels$label_long<- c("base model, AUC=0.86","base model+Tbio, AUC=0.92","base model+Tpre, AUC=0.90")
p1<-ggroc(roc4,legacy.axes = TRUE)+theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+guides(fill="none")+scale_color_discrete(name="",labels=data.labels$label_long)+
  geom_segment(aes(x=1,xend=0,y=1,yend=0),color="darkgrey",linetype="dashed")
ggsave(filename = "Multiple_roc.png",plot = p1,device = "png",dpi = 320)

cp_1<-ggplot(data=data.frame(obs=test_temp$composite,pre=a.2),aes(pre,obs)) +
  geom_segment(aes(x=1,xend=0,y=1,yend=0),color="darkgrey",linetype="dashed")+
  geom_smooth(method = stats::loess, se = FALSE) +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  xlab("Estimated Prob.") +
  ylab("Empirical Prob.") +
  theme(plot.title = element_text(hjust = 0.5)) +theme_bw()+
  theme(panel.grid.major=element_blank(),panel.border = element_rect(size=1.2),panel.grid.minor=element_blank(),axis.title.x = element_text(size=16),axis.title.y=element_text(size=16))+guides(fill="none")+scale_color_discrete(name="",labels=data.labels$label_long)
ggsave(filename = "calibration_plot_realT.png",plot = cp_1,device = "png",dpi = 320)
cp_2<-ggplot(data=data.frame(obs=test_temp$composite,pre=a.1),aes(pre,obs)) +
  geom_segment(aes(x=1,xend=0,y=1,yend=0),color="darkgrey",linetype="dashed") +
  geom_smooth(method = stats::loess, se = FALSE) +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  xlab("Estimated Prob.") +
  ylab("Empirical Prob.") +
  theme(plot.title = element_text(hjust = 0.5)) +theme_bw()+
  theme(panel.grid.major=element_blank(),panel.border = element_rect(size=1.2),panel.grid.minor=element_blank(),axis.title.x = element_text(size=16),axis.title.y=element_text(size=16))+guides(fill="none")+scale_color_discrete(name="",labels=data.labels$label_long)
ggsave(filename = "calibration_plot_noT.png",plot = cp_2,device = "png",dpi = 320)
cp_4<- ggroc(roc4,legacy.axes = TRUE)+theme_bw()+
  theme(panel.grid.major=element_blank(),panel.border = element_rect(size=1.2),panel.grid.minor=element_blank(),axis.title.x = element_text(size=16),axis.title.y=element_text(size=16))+guides(fill="none")+scale_color_discrete(name="",labels=data.labels$label_long)+
  geom_segment(aes(x=1,xend=0,y=1,yend=0),color="darkgrey",linetype="dashed")
ggsave(filename = "roc_plot.png",plot = cp_4,device = "png",dpi = 320)
cp_3<-ggplot(data=data.frame(obs=test_temp$composite,pre=a.3),aes(pre,obs)) +
  geom_segment(aes(x=1,xend=0,y=1,yend=0),color="darkgrey",linetype="dashed") +
  geom_smooth(method = stats::loess, se = FALSE) +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  xlab("Estimated Prob.") +
  ylab("Empirical Prob.") +
  theme(plot.title = element_text(hjust = 0.5)) +theme_bw()+
  theme(panel.grid.major=element_blank(),panel.border = element_rect(size=1.2),panel.grid.minor=element_blank(),axis.title.x = element_text(size=16),axis.title.y=element_text(size=16))+guides(fill="none")+scale_color_discrete(name="",labels=data.labels$label_long)
ggsave(filename = "calibration_plot_predictT.png",plot = cp_3,device = "png",dpi = 320)
p1 <-ggarrange(cp_1,cp_2,cp_3,labels = c("A","B","C"),ncol = 2,nrow = 2,)
ggsave(filename = "calibration_plot.png",plot = p1,device = "png",dpi = 320)

hoslem.test(x=test_temp$composite,y=predict(mod_xgb,as.matrix(test_temp[,-c(12,13,16:19)])))


p2<-ggroc(roc1,legacy.axes = TRUE,color="blue")+geom_text(aes(x=0.3,y=0.7,label="AUC=0.925"),size=5,check_overlap = TRUE)+theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+geom_abline(intercept = 0,slope=1)
temp=predict(mod_xgb,as.matrix(test_temp[,-12]))
sum(test_temp$composite==(temp>0.068))

X <-as.matrix(train_temp[,-c(12,13,16:19)])
dtrain <- xgb.DMatrix(X,label=as.numeric(train_temp$composite)-1)
param <- list(max_depth=5, eta = 0.05, silent = 1, objective = "binary:logistic")
watchlist <- list(train = dtrain)
mod_xgb <- xgb.train(param, dtrain, nrounds=110, watchlist,verbose = FALSE)
print(roc(test_temp$composite~predict(mod_xgb,as.matrix(test_temp[,-c(12,13,16:19)]))))
coords(roc(test_temp$composite~predict(mod_xgb,as.matrix(test_temp[,-c(12,13,16:19)]))),"best",transpose = FALSE)
roc2= roc(test_temp$composite~predict(mod_xgb,as.matrix(test_temp[,-c(12,13,16:19)])))
ggplot()+geom_line(aes(x=roc1$specificities,y=roc1$sensitivities))
p1<-ggroc(roc1,legacy.axes = TRUE,color="blue")+geom_text(aes(x=0.3,y=0.7,label="AUC=0.860"),size=5,check_overlap = TRUE)+theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+geom_abline(intercept = 0,slope=1)
ggsave(filename = "/UsersnoT_predicting.png",plot = p1,device = "png",dpi = 320)


roc.test(roc2,roc3)

test_temp1<- train_dat_w[id_temp,]
test_temp1 <- mutate(test_temp1,biT=ifelse(test_temp1$T>0,1,0))
X <-as.matrix(train_temp[,-c(12,13,16:19)])
dtrain <- xgb.DMatrix(X,label=as.numeric(train_temp$composite)-1)
param <- list(max_depth=5, eta = 0.05, silent = 1, objective = "binary:logistic")
watchlist <- list(train = dtrain)
mod_xgb <- xgb.train(param, dtrain, nrounds=110, watchlist,verbose = FALSE)
print(roc(test_temp$composite~predict(mod_xgb,as.matrix(test_temp1[,-c(12,13,16:19)]))))
coords(roc(test_temp$composite~predict(mod_xgb,as.matrix(test_temp1[,-c(12,13,16:19)]))),"best",transpose = FALSE)
roc3= roc(test_temp$composite~predict(mod_xgb,as.matrix(test_temp1[,-c(12,13,16:19)])))


temp=predict(mod_xgb,as.matrix(test_temp[,-c(12,13)]))

temp = predict(mod_ANN,newdata = test_temp[,-12])$probabilities[,2]

mod_svm <- svm(composite~.,data=train_temp,type="C-classification",probability=TRUE)
temp = predict(mod_svm,newdata = test_temp[,-12],probability=TRUE)
temp = attr(temp,"probabilities")[,2]
print(roc(test_temp$composite~temp))
coords(roc(test_temp$composite~temp),"best",transpose = FALSE)
###importance
train_dat_w <- mutate(train_dat_w,T=test_T$T)
train_temp <- train_dat_w[-id_temp,-13]
X <-as.matrix(train_temp[,-12])
dtrain <- xgb.DMatrix(X,label=as.numeric(train_temp$composite)-1)
param <- list(max_depth=5, eta = 0.05, silent = 1, objective = "binary:logistic")
watchlist <- list(train = dtrain)
mod_xgb <- xgb.train(param, dtrain, nrounds=110, watchlist,verbose = FALSE)
write.csv(xgb.importance(model=mod_xgb),"importance_clinicl_MESTC.csv")

X <-as.matrix(train_temp[,-(1:14)])
dtrain <- xgb.DMatrix(X,label=as.numeric(train_temp$composite)-1)
param <- list(max_depth=5, eta = 0.05, silent = 1, objective = "binary:logistic")
watchlist <- list(train = dtrain)
mod_xgb <- xgb.train(param, dtrain, nrounds=110, watchlist,verbose = FALSE)

write.csv(xgb.importance(model=mod_xgb),"importance_MESTC.csv") ###importance calculation

####################################################
########survival machine learning model training####
####################################################
library(survival)
colnames(follow_up_train)[c(17,18)] = c("jisu_mianyi","ACEI_ARB")
train_sur <- read.csv("train_sur.csv")
test_sur <- read.csv("test_sur.csv")
X <- as.matrix(train_sur[,c(1,2,3,4,5,6,7,8,9,10,11,14,17,18)])
surv.time <- train_sur$composite_time
surv.time[train_sur$composite==0] = -train_sur$composite_time[train_sur$composite==0]
dtrain <- xgb.DMatrix(X,label=surv.time)

param <- list(max_depth=5, eta = 0.05, silent = 1, objective = "survival:cox", eval_metric = "cox-nloglik")
watchlist <- list(train = dtrain)
bst.cli <- xgb.train(param, dtrain, nrounds=77, watchlist,verbose = FALSE)
Hmisc::rcorr.cens(-log(predict(bst.cli,as.matrix(test_sur[,c(1,2,3,4,5,6,7,8,9,10,11,14,17,18)]))),Surv(test_sur$composite_time,test_sur$composite))
Hmisc::rcorr.cens(-log(predict(bst.cli,as.matrix(train_sur[,c(1,2,3,4,5,6,7,8,9,10,11,14,17,18)]))),Surv(train_sur$composite_time,train_sur$composite))

follow_up_train_pT<-follow_up_train
follow_up_train_pT$biT <- train_predicted_T
test_sur_pT <- follow_up_train_pT[id_temp,]
compareC(test_sur$composite_time,test_sur$composite,predict(bst.cli,as.matrix(test_sur[,c(1,2,3,4,5,6,7,8,9,10,11,14,17,18)])),predict(bst.cli,as.matrix(test_sur_pT[,c(1,2,3,4,5,6,7,8,9,10,11,14,17,18)])))


train_sur <- follow_up_train[-id_temp,]
test_sur <- follow_up_train[id_temp,]
test_sur$biT <- train_predicted_T
X <- as.matrix(train_sur[,c(1,2,3,4,5,6,7,8,9,10,11,14,17,18)])
surv.time <- train_sur$composite_time
surv.time[train_sur$composite==0] = -train_sur$composite_time[train_sur$composite==0]
dtrain <- xgb.DMatrix(X,label=surv.time)

param <- list(max_depth=5, eta = 0.05, silent = 1, objective = "survival:cox", eval_metric = "cox-nloglik")
watchlist <- list(train = dtrain)
bst.cli <- xgb.train(param, dtrain, nrounds=80, watchlist,verbose = FALSE)
Hmisc::rcorr.cens(-log(predict(bst.cli,as.matrix(test_sur[,c(1,2,3,4,5,6,7,8,9,10,11,14,17,18)]))),Surv(test_sur$composite_time,test_sur$composite))
Hmisc::rcorr.cens(-log(predict(bst.cli,as.matrix(train_sur[,c(1,2,3,4,5,6,7,8,9,10,11,14,17,18)]))),Surv(train_sur$composite_time,train_sur$composite))
###no T
train_sur <- follow_up_train[-id_temp,]
test_sur <- follow_up_train[id_temp,]
X <- as.matrix(train_sur[,c(1,2,3,4,5,6,7,8,9,10,11,17,18)])
surv.time <- train_sur$composite_time
surv.time[train_sur$composite==0] = -train_sur$composite_time[train_sur$composite==0]
dtrain <- xgb.DMatrix(X,label=surv.time)

param <- list(max_depth=5, eta = 0.05, silent = 1, objective = "survival:cox", eval_metric = "cox-nloglik")
watchlist <- list(train = dtrain)
bst.cli <- xgb.train(param, dtrain, nrounds=80, watchlist,verbose = FALSE)
Hmisc::rcorr.cens(-log(predict(bst.cli,as.matrix(test_sur[,c(1,2,3,4,5,6,7,8,9,10,11,17,18)]))),Surv(test_sur$composite_time,test_sur$composite))
Hmisc::rcorr.cens(-log(predict(bst.cli,as.matrix(train_sur[,c(1,2,3,4,5,6,7,8,9,10,11,17,18)]))),Surv(train_sur$composite_time,train_sur$composite))















