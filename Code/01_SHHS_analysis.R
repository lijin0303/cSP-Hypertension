set.seed(123)
setwd("~/Documents/HTN_SDB/")
source("code/00_Utilities.R")
#===== SHHS:cSDB PSG =====
##### Penalized regression #####
SHHS_process_PSG <- SHHSprocess(exclusionKEY = ">=[ ]{0,1}[245]{1}%",assothres = 0.2,misthres = 0.3,Drange ="complete",lung=T)
winsorize_examination(D = SHHS_process_PSG$Dat,vec =SHHS_process_PSG$candi,DR ="complete")
WR <-  c("AvDNOP3","dMxBROH","HREMt1P",
         "MnHROA","MnHROP","MxDNBP",
         "OANOA","OARDNBA","OARDNBP",
         "REMt1P","SlpLatP")
WL <- NULL
SHHS_PSG_elastic  <- penalLogit(D_list=SHHS_process_PSG,Alpha=0.5,sdt=T,winthorize=T,winRight=WR,winLeft=WL,DR = "complete")
##### Table 1 #####
shhs_ANAset <- readRDS("Result/Fixed_rlt_0530/PSG/WinsorD_PSG.rds") #use the result generated from previous code in 5/30
shhs_ANAset[,SmokingStatus:=ifelse(SmokingStatus==0,"Never",ifelse(SmokingStatus==1,"Current","Former"))]
tableshhs_complete <- CreateTableOne(vars = colnames(shhs_ANAset), data = shhs_ANAset, 
                            factorVars = c("Sex","Race","SmokingStatus"),strata = "IncHT")

tableOneCSV(tableOne = tableshhs_complete,savefile = "~/Documents/HTN_SDB/Result/Fixed_rlt_0530/tableOne_SHHS_PSG.csv")
##### Clustering #####
colorVals <- c("#FA8072", "#808000", "#3CB371", "#4682B4","#BA55D3","#FFD500")
shhs_cluster_PSG<- SHHS_PSG_elastic$penalizedPool[,10:51]%>%
  t()%>%
  set_colnames(SHHS_process_PSG$Dat%>%pull(SHHSid))
# Diagnosis plot
pdf("Result/Fixed_rlt_0530/diagnosis_PSG_elbow.pdf",width = 5,height = 5)
fviz_nbclust(shhs_cluster_PSG, FUN = hcut, method = "wss")
dev.off()
cls_PSG <- 6
# Hierarchical clustering
clustered_SHHS_PSG<- daisyHIER(k=cls_PSG,labelColors =colorVals[1:cls_PSG],data =shhs_cluster_PSG,jpgname = "hier_PSG")

ggplot(prcomp(shhs_cluster_PSG, scale = FALSE)$x%>%
         as.data.frame()%>%
         mutate(traits = rownames(shhs_cluster_PSG),
                hier = as.factor(clustered_SHHS_PSG)),
       aes(x=PC1,y=PC2,label=traits,color=hier))+
  scale_color_manual(values = 
                       colorVals[1:cls_PSG])+
  geom_point()+
  geom_label(cex=1.5)+
  theme_bw()+
  theme(axis.text=element_text(size=5),
        axis.title=element_text(size=5,face="bold"))+
  theme(legend.position="none",
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank())+
  ggsave("./Result/Fixed_rlt_0530/pca_PSG.jpg",
         width = 4, height = 3, units = "in",
         dpi = 2000)
##### cSDB PSG description #####
SHHS_PSG_elastic$coef[,.(`Selected Traits`=Trait,Weight=Effect,Description = display_name)]%>%
  mutate(Weight = formatC(Weight,digits = 2,format = "E"))%>%
  openxlsx::write.xlsx("Result/Fixed_rlt_0530/cSDB_PSG.xlsx")

shhs_ANAset$nSDB <- as.matrix(shhs_ANAset[,SHHS_PSG_elastic$coef$Trait,with=F])%*%SHHS_PSG_elastic$coef$Effect
ggplot(shhs_ANAset, aes(x=nSDB)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") +
  theme_bw()+
  theme(axis.text=element_text(size=5),
        axis.title=element_text(size=5,face="bold"))+
  theme(legend.position="none",
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank())+
  labs(x="cSP")+
  ggsave(filename="Result/Fixed_rlt_0530/PSG/histSHHS_PSG.jpg",
         width = 3, height = 3, units = "in",
         dpi = 2000)

ggplot(shhs_ANAset, aes(x=nSDB,color=IncHT,fill=IncHT)) + 
  #geom_histogram(aes(y=..density..), position="identity", alpha=0.5)+
  geom_density(alpha=0.6) +
  scale_color_brewer(palette="Dark2",name = "Incident HTN",label=c("No","Yes"))+
  scale_fill_brewer(palette="Dark2",name = "Incident HTN",label=c("No","Yes"))+
  theme(
    plot.title = element_text(color='black', hjust = 0.5),
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    panel.grid = element_blank(),
    axis.text = element_text(color='black')
    )+labs(x="cSP")+
  ggsave(filename="~/Desktop/histSHHS_PSG_HTN.jpg",
         width = 6, height = 6, units = "in",
         dpi = 2000)


shhs_ANAset%>%mutate(q_csp = case_when(nSDB<=quantile(nSDB,0.25)~"Q1",
                                       nSDB>quantile(nSDB,0.25) &nSDB<=quantile(nSDB,0.5)~"Q2",
                                       nSDB>quantile(nSDB,0.5) &nSDB<=quantile(nSDB,0.75)~"Q3",
                                       TRUE~"Q4"))%>%
  group_by(q_csp,IncHT)%>%
  summarise(n=n())%>%
  tidyr::spread(IncHT,n)%>%
  mutate(total = `0`+`1`)%>%
  mutate(noHTN = round(`0`*100/total,2),
         withHTN = round(`1`*100/total,2))%>%
  select(q_csp,noHTN,withHTN)%>%
  tidyr::gather(HTNstatus,Percent,-q_csp)%>%
  mutate(label = paste0(round(Percent),"%"))%>%
  ggplot(., aes(fill=HTNstatus, y=Percent, x=q_csp)) + 
  geom_bar(position="stack", stat="identity")+
  viridis::scale_fill_viridis(discrete = T,alpha = 0.95,
                              name="Incident HTN",
                              label=c("No","Yes"))+
  geom_text(
    data = .%>%filter(HTNstatus!="noHTN"),
    aes(x=q_csp,y = rep(12,4), label = label), size =3,
    family = "Times New Roman", fontface = "bold")+
  theme(
    plot.title = element_text(color='black', hjust = 0.5),
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    panel.grid = element_blank(),
    axis.text = element_text(color='black')
  )+labs(x="Quantiles of cSP",y="Percentage")+
ggsave(width = 6,height = 6,filename = "~/Desktop/StackedPercent_all.png",
       units = "in",dpi=1500)

#===== SHHS:cSDB HST =====
SHHS_process_HST<- SHHSprocess(exclusionKEY = ">=[ ]{0,1}[245]{1}%",assothres = 0.2,misthres = 0.3,Drange ="nonEGG")
winsorize_examination(D = SHHS_process_HST$Dat,vec =SHHS_process_HST$candi,DR ="nonEGG")
SHHS_HST_elastic <- penalLogit(D_list=SHHS_process_HST,Alpha=0.5,sdt=T,winthorize=T,winRight=c("mnhop"),winLeft=NULL,DR = "nonEGG")
##### Table 1 #####
shhs_ANAset_hst <- readRDS("Result/Fixed_rlt_0530/HST/WinsorD_HST.rds") #use the result generated from previous code in 5/30
shhs_ANAset_hst[,SmokingStatus:=ifelse(SmokingStatus==0,"Never",ifelse(SmokingStatus==1,"Current","Former"))]
tableshhs_HST <- CreateTableOne(vars = colnames(shhs_ANAset_hst), data = shhs_ANAset_hst, 
                                     factorVars = c("Sex","Race","SmokingStatus"),strata = "IncHT")

tableOneCSV(tableOne = tableshhs_HST,savefile = "~/Documents/HTN_SDB/Result/Fixed_rlt_0530/HST/tableOne_SHHS_HST.csv")
##### Clustering #####
shhs_cluster_HST <-  SHHS_HST_elastic$penalizedPool[,10:ncol(SHHS_HST_elastic$penalizedPool)]%>%
  t()%>%
  set_colnames(SHHS_process_HST$Dat%>%pull(SHHSid))
pdf("./Result/Fixed_rlt_0530/HST/diagnosis_HST_wss.pdf",width = 5,height = 5)
fviz_nbclust(shhs_cluster_HST, FUN = hcut, method = "wss")
dev.off()
cls_HST <- 5
# Hierarchical clustering
clustered_SHHS_HST<- daisyHIER(k=cls_HST,labelColors =colorVals[1:cls_HST],data =shhs_cluster_HST,jpgname = "hier_HST")

ggplot(prcomp(shhs_cluster_HST, scale = FALSE)$x%>%
         as.data.frame()%>%
         mutate(traits = rownames(shhs_cluster_HST),
                hier = as.factor(clustered_SHHS_HST)),
       aes(x=PC1,y=PC2,label=traits,color=hier))+
  scale_color_manual(values = 
                       colorVals[1:cls_HST])+
  geom_point()+
  geom_label(cex=1.5)+
  theme_bw()+
  theme(axis.text=element_text(size=5),
        axis.title=element_text(size=5,face="bold"))+
  theme(legend.position="none",
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank())+
  ggsave("./Result/Fixed_rlt_0530/HST/pca_HST.jpg",
         width = 4, height = 3, units = "in",
         dpi = 2000)
##### cSDB HST description #####
SHHS_HST_elastic$coef[,.(`Selected Traits`=Trait,Weight=Effect,Description = display_name)]%>%
  mutate(Weight = formatC(Weight,digits = 2,format = "E"))%>%
  openxlsx::write.xlsx("Result/Fixed_rlt_0530/HST/cSDB_HST.xlsx")

shhs_ANAset_hst$nSDB <- as.matrix(shhs_ANAset_hst[,SHHS_HST_elastic$coef$Trait,with=F])%*%SHHS_HST_elastic$coef$Effect
ggplot(shhs_ANAset_hst, aes(x=nSDB)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") +
  theme_bw()+
  theme(axis.text=element_text(size=5),
        axis.title=element_text(size=5,face="bold"))+
  theme(legend.position="none",
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank())+
  labs(x="cSP")+
  ggsave(filename="Result/Fixed_rlt_0530/HST/histSHHS_HST.jpg",
         width = 3, height = 3, units = "in",
         dpi = 2000)

#===== SHHS: covariates recall for follow-up analysis =====
shhs1SDB <- fread("Data/shhs-data-dictionary-0.15.0-variables.csv",header = T)
interested_comp <- shhs1SDB%>%filter(grepl("AHI.*>= 3%.*all apneas",display_name))%>%pull(id)
shhs1 <- fread("Data/shhs1-dataset-0.15.0.csv",header = T)
shhs1%<>%
  mutate(WHratio = waist/Hip)%>%
  select(SHHSid = nsrrid,
         Age = age_s1,
         Sex = gender,
         Race = race,
         Height = height,
         BaseHT = HTNDerv_s1,
         BMI = bmi_s1,
         WHratio = WHratio,
         NeckGrith = NECK20,
         PkYr = CgPkYr,
         SmokingStatus = smokstat_s1,
         AHI_3oxydesat = ahi_a0h3,
         AHI_4oxydesat = ahi_a0h4,
         OAHI4P = oahi,
         AHI_3oxydesat_arousal = ahi_a0h3a,
         sleepiness_start = ESS_s1)
shhs2 <- fread("Data/shhs2-dataset-0.15.0.csv",
               header = T)[,.(nsrrid,HTNDerv_s2,ess_s2)]
colnames(shhs2)<-c("SHHSid","IncHT","sleepiness_end")
shhs <- merge(shhs1,shhs2,by="SHHSid")
#===== Analysis with fixed data results =====
PSG_comb <- readRDS("~/Downloads/本月钥匙/cSDB_manuscript/fixed_rlt/PSGcomb_SHHS.rds")
HST_comb <- readRDS("~/Downloads/本月钥匙/cSDB_manuscript/fixed_rlt/HST_comb_SHHS.rds")
shhs <- readRDS("~/Downloads/本月钥匙/cSDB_manuscript/fixed_rlt/Compare_cov_shhs.rds")
shhsHP <- shhs%>%filter(BaseHT==0)
shhsHP%<>%
  mutate_at(c("Sex","Race","IncHT","SmokingStatus"),as.factor)%>%
  mutate(AHI_cat5 = case_when(AHI_3oxydesat>5~1,TRUE~0),
         AHI_cat15 = case_when(AHI_3oxydesat>15~1,TRUE~0))
glmAHI1 <- fit_model_shhs(paste("IncHT~","AHI_3oxydesat",
                     "+ Age + as.factor(Sex) + as.factor(Race) + BMI + WHratio + NeckGrith + as.factor(SmokingStatus) + PkYr",sep=""),
                     shhsHP)
ModelMetrics::auc(glmAHI1)
coef(summary(glmAHI1))[2,4]

fit_model_shhs(paste("IncHT~","AHI_cat5",
                     "+ Age + as.factor(Sex) + as.factor(Race) + BMI + WHratio + NeckGrith + as.factor(SmokingStatus) + PkYr",sep=""),
               shhsHP)%>%summary()%>%coef()%>%.[2,]

fit_model_shhs(paste("IncHT~","AHI_cat15",
                     "+ Age + as.factor(Sex) + as.factor(Race) + BMI + WHratio + NeckGrith + as.factor(SmokingStatus) + PkYr",sep=""),
               shhsHP)%>%summary()%>%coef()%>%.[2,]

PSG_glm <- fit_model_shhs(paste("IncHT~","nSDB",
                     "+ Age + as.factor(Sex) + as.factor(Race) + BMI + WHratio + NeckGrith + as.factor(SmokingStatus) + PkYr",sep=""),
               PSG_comb$Data)
ModelMetrics::auc(PSG_glm)
coef(summary(PSG_glm))[2,]

HST_glm <- fit_model_shhs(paste("IncHT~","nSDB",
                     "+ Age + as.factor(Sex) + as.factor(Race) + BMI + WHratio + NeckGrith + as.factor(SmokingStatus) + PkYr",sep=""),
                     HST_comb$Data)
ModelMetrics::auc(HST_glm)
coef(summary(HST_glm))[2,]

p_cSP <- ggplot(PSG_comb$Data, aes(x=nSDB)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") +
  theme_bw()+
  theme(axis.text=element_text(size=5),
        axis.title=element_text(size=5,face="bold"))+
  theme(legend.position="none",
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank())+
  labs(x="cSP")

p_AHI <- ggplot(shhsHP, aes(x=AHI_3oxydesat)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="dodgerblue3") +
  theme_bw()+
  theme(axis.text=element_text(size=5),
        axis.title=element_text(size=5,face="bold"))+
  theme(legend.position="none",
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank())+
  labs(x="AHI")

require(ggpubr)
ggarrange(p_cSP,p_AHI,ncol=2)%>%
  ggsave(.,filename="~/Desktop/histogram_ahi_cSP.jpg",
         width = 6, height = 3, units = "in",
         dpi = 2000)

#===== Incident Sleepiness investigation =====
ESS_dat <- shhs%>%
  mutate(base_catESS = case_when(sleepiness_start>10~1,
                                 TRUE~0),
         Inc_catESS =  case_when(sleepiness_end>10~1,
                                 TRUE~0))%>%
  filter(base_catESS==0)%>%
  inner_join(PSG_comb$Data%>%select(SHHSid,nSDB),by="SHHSid")


ESS_dat%>%select(nSDB,Inc_catESS)%>%
  ggplot(., aes(x=nSDB,color=factor(Inc_catESS),fill=factor(Inc_catESS))) + 
  geom_density(alpha=0.6) +
  scale_color_brewer(palette="Dark2",name = "Incident Sleepiness",label=c("No","Yes"))+
  scale_fill_brewer(palette="Dark2",name = "Incident Sleepiness",label=c("No","Yes"))+
  theme(
    plot.title = element_text(color='black', hjust = 0.5),
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    panel.grid = element_blank(),
    axis.text = element_text(color='black')
  )+labs(x="cSP")+
  ggsave(filename="~/Desktop/histSHHS_PSG_ESS.jpg",
         width = 6, height = 6, units = "in",
         dpi = 2000)

PSG_ESS <- fit_model_shhs(paste("Inc_catESS~","nSDB",
                     "+ Age + as.factor(Sex) + as.factor(Race) + BMI + WHratio + NeckGirth + as.factor(SmokingStatus) + PkYr",sep=""),
               ESS_dat)

AHI_ESS <- fit_model_shhs(paste("Inc_catESS~","AHI_3oxydesat",
                     "+ Age + as.factor(Sex) + as.factor(Race) + BMI + WHratio + NeckGirth + as.factor(SmokingStatus) + PkYr",sep=""),
               ESS_dat)

auc(PSG_ESS)
auc(AHI_ESS)
coef(summary(PSG_ESS))[2,4]
coef(summary(AHI_ESS))[2,4]

#===== Other lung traits investigation =====
require(nricens)
PSG_comb <- readRDS("~/Downloads/本月钥匙/cSDB_manuscript/fixed_rlt/PSGcomb_SHHS.rds")
HST_comb <- readRDS("~/Downloads/本月钥匙/cSDB_manuscript/fixed_rlt/HST_comb_SHHS.rds")
shhs <- readRDS("~/Downloads/本月钥匙/cSDB_manuscript/fixed_rlt/Compare_cov_shhs.rds")
shhsHP <- shhs%>%filter(BaseHT==0)
shhsHP%<>%
  mutate_at(c("Sex","Race","IncHT","SmokingStatus"),as.factor)
SHHS_process_PSG <- SHHSprocess(exclusionKEY = ">=[ ]{0,1}[245]{1}%",assothres = 0.2,misthres = 0.3,Drange ="complete",lung=T)
model_compare <- function(D,predictor){
  D <- na.omit(D[,c("SHHSid","IncHT",predictor,"Age","Sex","Race","BMI","WHratio","NeckGrith","SmokingStatus","PkYr"),with=F])%>%
    as.data.frame()
  sd_pred <- sd(D%>%pull(predictor),na.rm = T)
  expr <- paste("IncHT~",predictor,
                "+ Age + as.factor(Sex) + as.factor(Race) + BMI + WHratio + NeckGrith + as.factor(SmokingStatus) + PkYr",sep="")
  glmodel <-  glm(as.formula(expr),binomial(logit), data = D, x=TRUE)
  cut_glm <- Youden_max_cut(D = data.frame(prob = predict(object = glmodel, type = 'response'),
                                           label = D$IncHT))
  CI <- exp(confint(glmodel)[2,]*sd_pred)%>%round(.,2)
  return(list(stats = data.frame(N = nrow(D),
                                 OR =exp(coef(summary(glmodel))[2,1])%>%round(.,2),
                                 OR_sd =exp(coef(summary(glmodel))[2,1]*sd_pred)%>%round(.,2),
                                 l_OR = CI[1],
                                 u_OR = CI[2],
                    pval = coef(summary(glmodel))[2,4],
                    AUC = ModelMetrics::auc(glmodel),
                    cut = cut_glm),
                    model = glmodel,
                    newD = D,sd = sd_pred)
         )
  
}
Youden_max_cut <- function(D,l=seq(0,1,0.001)){
  colnames(D) <- c("prob","label")
  youden_D <- sapply(l,function(s){
    true_pos <- D%>%filter(prob>s&label==1)%>%nrow()
    true_neg <- D%>%filter(prob<s&label==0)%>%nrow()
    false_pos <- D%>%filter(prob>s&label==0)%>%nrow()
    false_neg <- D%>%filter(prob<s&label==1)%>%nrow()
    youden <- true_pos/(true_pos+false_neg) +  true_neg/(true_neg+false_pos) -1
    return(youden)
  })
  cutoff_choice <- l[youden_D==max(youden_D)]
  return(cutoff_choice)}
PSG_model <- model_compare(D = PSG_comb$Data,predictor="nSDB")
fev1pp_model<- model_compare(D = SHHS_process_PSG$wholeD%>%filter(SHHSid %in% PSG_model$newD$SHHSid),predictor="FEV1pp")
NRI_fev1pp <- nribin(mdl.std = glm(IncHT ~ AHI_3oxydesat + Age + as.factor(Sex) + as.factor(Race) + BMI + WHratio + PkYr + as.factor(SmokingStatus) + NeckGrith,
                                binomial(logit), data = shhsHP%>%filter(SHHSid %in% PSG_model$newD$SHHSid), x=TRUE),
                  mdl.new =fev1pp_model$model, 
                  cut =fev1pp_model$stats$cut, updown = 'category',niter=0)

fvcpp_model<- model_compare(D = SHHS_process_PSG$wholeD%>%filter(SHHSid %in% PSG_model$newD$SHHSid),predictor="FVCpp")
NRI_fvcpp <- nribin(mdl.std = glm(IncHT ~ AHI_3oxydesat + Age + as.factor(Sex) + as.factor(Race) + BMI + WHratio + PkYr + as.factor(SmokingStatus) + NeckGrith,
                                   binomial(logit), data = shhsHP%>%filter(SHHSid %in% PSG_model$newD$SHHSid), x=TRUE),
                     mdl.new =fvcpp_model$model, 
                     cut =fvcpp_model$stats$cut, updown = 'category',niter=0)
PSG_model <- model_compare(D = PSG_comb$Data,predictor="nSDB")
NRI_PSG <- nribin(mdl.std = glm(IncHT ~ AHI_3oxydesat + Age + as.factor(Sex) + as.factor(Race) + BMI + WHratio + PkYr + as.factor(SmokingStatus) + NeckGrith,
                                  binomial(logit), data = shhsHP%>%filter(SHHSid %in% PSG_model$newD$SHHSid), x=TRUE),
                    mdl.new =PSG_model$model, 
                    cut =PSG_model$stats$cut, updown = 'category',niter=0)
HST_model <- model_compare(D = HST_comb$Data%>%filter(SHHSid %in% PSG_model$newD$SHHSid),predictor="nSDB")
NRI_HST <- nribin(mdl.std = glm(IncHT ~ AHI_3oxydesat + Age + as.factor(Sex) + as.factor(Race) + BMI + WHratio + PkYr + as.factor(SmokingStatus) + NeckGrith,
                                binomial(logit), data = shhsHP%>%filter(SHHSid %in% PSG_model$newD$SHHSid), x=TRUE),
                  mdl.new =HST_model$model, 
                  cut =HST_model$stats$cut, updown = 'category',niter=0)

WR <-  c("AvDNOP3","dMxBROH","HREMt1P",
         "MnHROA","MnHROP","MxDNBP",
         "OANOA","OARDNBA","OARDNBP",
         "REMt1P","SlpLatP")
WL <- NULL
SHHS_process_PSG_nolung <- SHHSprocess(exclusionKEY = ">=[ ]{0,1}[245]{1}%",assothres = 0.2,misthres = 0.3,Drange ="complete",lung=F)
SHHS_PSG_elastic_nolung  <- penalLogit(D_list=SHHS_process_PSG_nolung,Alpha=0.5,sdt=T,winthorize=T,winRight=WR,winLeft=WL,DR = "complete")
PSG_nolungD <- SHHS_process_PSG_nolung$Dat
PSG_nolungD[,SmokingStatus:=ifelse(SmokingStatus==0,"Never",ifelse(SmokingStatus==1,"Current","Former"))]
PSG_nolungD$nSDB <- as.matrix(PSG_nolungD[,SHHS_PSG_elastic_nolung$coef$Trait,with=F])%*%SHHS_PSG_elastic_nolung$coef$Effect
PSG_nolung_model <- model_compare(D = PSG_nolungD%>%filter(SHHSid %in% PSG_model$newD$SHHSid),predictor="nSDB")
NRI_PSG_nolung <- nribin(mdl.std = glm(IncHT ~ AHI_3oxydesat + Age + as.factor(Sex) + as.factor(Race) + BMI + WHratio + PkYr + as.factor(SmokingStatus) + NeckGrith,
                                binomial(logit), data = shhsHP%>%filter(SHHSid %in% PSG_model$newD$SHHSid), x=TRUE),
                  mdl.new =PSG_nolung_model$model, 
                  cut =PSG_nolung_model$stats$cut, updown = 'category',niter=0)

SHHS_process_HST_nolung <- SHHSprocess(exclusionKEY = ">=[ ]{0,1}[245]{1}%",assothres = 0.2,misthres = 0.3,Drange ="nonEGG",lung=F)
SHHS_HST_elastic_nolung  <- penalLogit(D_list=SHHS_process_HST_nolung,Alpha=0.5,sdt=T,winthorize=T,winRight=c("mnhop"),winLeft=NULL,DR = "nonEGG")
HST_nolungD <- SHHS_process_HST_nolung$Dat
HST_nolungD[,SmokingStatus:=ifelse(SmokingStatus==0,"Never",ifelse(SmokingStatus==1,"Current","Former"))]
HST_nolungD$nSDB <- as.matrix(HST_nolungD[,SHHS_HST_elastic_nolung$coef$Trait,with=F])%*%SHHS_HST_elastic_nolung$coef$Effect
HST_nolung_model <- model_compare(D = HST_nolungD%>%filter(SHHSid %in% PSG_model$newD$SHHSid),predictor="nSDB")
NRI_HST_nolung <- nribin(mdl.std = glm(IncHT ~ AHI_3oxydesat + Age + as.factor(Sex) + as.factor(Race) + BMI + WHratio + PkYr + as.factor(SmokingStatus) + NeckGrith,
                                       binomial(logit), data = shhsHP%>%filter(SHHSid %in% PSG_model$newD$SHHSid), x=TRUE),
                         mdl.new =HST_nolung_model$model, 
                         cut =HST_nolung_model$stats$cut, updown = 'category',niter=0)


ahi3_model <- model_compare(D = shhsHP%>%filter(SHHSid %in% PSG_model$newD$SHHSid),predictor="AHI_3oxydesat")
ahi3_model_arousal <- model_compare(D = shhsHP%>%filter(SHHSid %in% PSG_model$newD$SHHSid),predictor="AHI_3oxydesat_arousal")
ahi4_model <- model_compare(D = shhsHP%>%filter(SHHSid %in% PSG_model$newD$SHHSid),predictor="AHI_4oxydesat")
oahi4_model <- model_compare(D = shhsHP%>%filter(SHHSid %in% PSG_model$newD$SHHSid),predictor="OAHI4P")


sd_collect <- sapply(list(fev1pp_model,fvcpp_model,PSG_model,HST_model,PSG_nolung_model,HST_nolung_model,ahi3_model),function(s) s$sd)
sd_SHHS <- data.frame(sd = sd_collect,
           method = c("FEV1","FVC","PSG","HST","PSG_nolung","HST_nolung","AHI"))
saveRDS(sd_SHHS,"~/Desktop/model_pheno_sd.rds")

#===== Correlation plot =====
require(RColorBrewer)
save_pheatmap_pdf <- function(x, filename, width=4.5, height=7.5) {
  pdf(filename, width = width, height = height)
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png <- function(x, filename, width=4.5, height=7.5) {
  png(filename, width = width, height = height,res=1200,units="in")
  grid::grid.draw(x$gtable)
  dev.off()
}
pheatmap::pheatmap(cor(shhs_ANAset[,SHHS_PSG_elastic$coef$Trait,with=F]),
                   color = colorRampPalette(
                     brewer.pal(n = 9, name ="GnBu"))(100),
                   cluster_cols = T,
                   cluster_rows = T,
                   cellwidth = 20,
                   cellheight =20,
                   fontsize_row = 5.3,
                   fontsize_col = 5.3,
                   show_colnames= TRUE,angle_col = "45")%>%
  save_pheatmap_png(., "~/Desktop/cor_PSG_pool.png",height = 8,width = 8)
pheatmap::pheatmap(cor(shhs_ANAset_hst[,SHHS_HST_elastic$coef$Trait,with=F]),
                   color = colorRampPalette(
                     brewer.pal(n = 9, name ="GnBu"))(100),
                   cluster_cols = T,
                   cluster_rows = T,
                   cellwidth = 20,
                   cellheight = 20,
                   fontsize_row = 4.5,
                   fontsize_col = 4.5,
                   show_colnames= TRUE,angle_col = "45")%>%
  save_pheatmap_png(., "~/Desktop/cor_HST_pool.png",height = 4.7,width = 4.7)


