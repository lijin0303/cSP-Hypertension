set.seed(123)
source("Code/00_Utilities.R")
#####  Compose cSP traits in MESA  ##### 
MESA_PSG<- MESAprocess(penalcomb=SHHS_PSG_elastic,impuNO=10)
MESA_HST<- MESAprocess(penalcomb=SHHS_HST_elastic,impuNO=10)
MESA_PSG_nolung <- MESAprocess(penalcomb=SHHS_PSG_elastic_nolung,impuNO=10)
MESA_HST_nolung <- MESAprocess(penalcomb=SHHS_HST_elastic_nolung,impuNO=10)
##### Figure S1 #####
mesaD_PSG <- MESA_PSG$OriDat
colnames(mesaD_PSG)[grep("WHr",colnames(mesaD_PSG))] <- "WHR"
mesaD_HST <- MESA_HST$OriDat
colnames(mesaD_HST)[grep("WHr",colnames(mesaD_HST))] <- "WHR"
mesaD_cSDB <- mesaD_HST%>%
  inner_join(mesaD_PSG,by = intersect(colnames(mesaD_HST),colnames(mesaD_PSG)))
jpeg('Result/Fixed_rlt_0530/missing_MESA.jpg',
     width = 6, height = 5.5, units = 'in', res = 400)
gg_miss_upset(mesaD_cSDB,nsets = 8,nintersects=15)
dev.off()
##### Table 3 #####
covsub <- mesaD_PSG[,.(Age,Sex,Race,BMI,WHR,HT,PkYr,SmokingStatus)]
covsub[,SmokingStatus:=ifelse(SmokingStatus==0,"Never",ifelse(SmokingStatus==1,"Former","Current"))]
covsub[,Race:=ifelse(Race==1,"White, Caucasian",
                     ifelse(Race==3,"Black, African-American",ifelse(Race==4,"Hispanic","Asian")))]
covsub[,Sex:=ifelse(Sex==1,"Male","Female")]
tablecov <- CreateTableOne(vars = colnames(covsub), data = covsub, 
                           factorVars = c("Sex","Race","HT","SmokingStatus"),
                           strata = "HT")
tableOneCSV(tableOne = tablecov,savefile = "Result/Fixed_rlt_0530/PSG/tableOne_MESA.csv")
##### Table S3 #####
statcSDB_PSG <- data.table(ImputationNO=c(1:10),
                           Mean= round(unlist(with(MESA_PSG$PSGIMPU,mean(CombSDB))$analyses),4),
                           Sd = round(unlist(with(MESA_PSG$PSGIMPU,sd(CombSDB))$analyses),4),
                           Quartile25 = round(unlist(with(MESA_PSG$PSGIMPU,quantile(CombSDB)[2])$analyses),4),
                           Quartile75 = round(unlist(with(MESA_PSG$PSGIMPU,quantile(CombSDB)[4])$analyses),4))

statcSDB_HST <- data.table(ImputationNO=c(1:10),
                           Mean= round(unlist(with(MESA_HST$PSGIMPU,mean(CombSDB))$analyses),4),
                           Sd = round(unlist(with(MESA_HST$PSGIMPU,sd(CombSDB))$analyses),4),
                           Quartile25 = round(unlist(with(MESA_HST$PSGIMPU,quantile(CombSDB)[2])$analyses),4),
                           Quartile75 = round(unlist(with(MESA_HST$PSGIMPU,quantile(CombSDB)[4])$analyses),4))

statcSDB_PSG%>%mutate(source="PSG")%>%rbind(statcSDB_HST%>%mutate(source="HST"))%>%
  openxlsx::write.xlsx("Result/Fixed_rlt_0530/cSDB_MESA_10impu.xlsx")
##### Figure 3 #####
mesa6 <- read.csv("~/Documents/HTN_SDB/Data/mesae6_finallabel_20200811_htn_variables.csv")
PSG_mesa6 <- complete(MESA_PSG$PSGIMPU, action = "long", include = TRUE)%>%
  inner_join(mesaD_PSG%>%
               mutate(.id= 1:n())%>%
               inner_join(mesa6,by="idno")%>%
               mutate(mesa6_HT = case_when(htnmed6c==1|htn6c==1~1,TRUE~0))%>%
               select(.id,mesa6_HT,Age,Sex,Race,BMI,WHR,PkYr,SmokingStatus,HT),
             by=".id")
imp_avg_comb <- PSG_mesa6%>%
  filter(.imp!=0)%>%
  select(.id,.imp,HT,mesa6_HT,CombSDB)%>%
  group_by(.id,HT,mesa6_HT)%>%
  summarise(imp_avg = mean(CombSDB))%>%
  mutate(HT_status = factor(case_when(HT==1~"MESA 5 cases",
                                      HT==0&mesa6_HT==1~"MESA 6 cases",
                                      TRUE~"Non-cases"),levels = c("Non-cases",
                                                                   "MESA 6 cases",
                                                                   "MESA 5 cases")))

perc_Q10 <- imp_avg_comb%>%
  as.data.frame()%>%
  dplyr::mutate(cutcSP = gtools::quantcut(imp_avg_comb$imp_avg, q=10))%>%
  group_by(cutcSP,HT_status)%>%
  summarise(n=n())%>%
  tidyr::spread(HT_status,n)%>%
  arrange(cutcSP)%>%
  as.data.frame()%>%
  mutate(Quantile = factor(sapply(1:10,function(s) paste0("Q",s)),
                           levels = sapply(1:10,function(s) paste0("Q",s))))%>%
  mutate(total = `MESA 5 cases` + `MESA 6 cases` + `Non-cases`)%>%
  mutate(`MESA 5 cases` = round(`MESA 5 cases`*100/total,2),
         `MESA 6 cases` = round(`MESA 6 cases`*100/total,2),
         `Non-cases` = round(`Non-cases`*100/total,2))%>%
  select(Quantile,`MESA 5 cases`,`MESA 6 cases`,`Non-cases`)%>%
  tidyr::gather(HTNstatus,Percent,-Quantile)%>%
  dplyr::mutate(label = paste0(sprintf("%.0f",Percent), "%"))%>%
  mutate(HTNstatus = factor(HTNstatus,
                            levels = c("Non-cases",
                                       "MESA 6 cases",
                                       "MESA 5 cases")))
q_mesa5 <- perc_Q10%>%filter(HTNstatus=="MESA 5 cases")%>%pull(Percent)
ggplot(perc_Q10, aes(fill=HTNstatus, y=Percent, x=Quantile)) + 
  geom_bar(position="stack", stat="identity")+
  geom_text(
    data = .%>%filter(HTNstatus!="Non-cases"),
    aes(x=Quantile,y = c(rep(25,10),q_mesa5+5), label = label), size =3,
    family = "Times New Roman", fontface = "bold") +
  viridis::scale_fill_viridis(discrete = T,alpha = 0.95,
                              name="HTN status")+
  theme(
    plot.title = element_text(color='black', hjust = 0.5),
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    panel.grid = element_blank(),
    axis.text = element_text(color='black'),
    axis.title= element_blank())+
  ggsave(filename = "~/Desktop/stackedBar_MESA5_cSP_PSG_quantile10.png",
         width=10,height=6)
##### Table 5 #####
mesa_cov_HT <- fread("Data/mesa_sdbtraits_20190827.txt",header=T,data.table = F)
mesa_cov_demgrh <- fread("Data/mesa_ms5527_20190130.csv",header=T,data.table = F)
mesa_cov_demgrh <- mesa_cov_demgrh[,grepALL("idno|bmi|race|age|gender|ppfvc|ppfev1|pkyrs5c|cig5c",colnames(mesa_cov_demgrh),names = T)]%>%
  magrittr::set_colnames(c("idno","BMI","Race","Age","Sex","FVCpp","FEV1pp","PkYr","SmokingStatus"))
mesa_cov_HT <- mesa_cov_HT[,grepALL("idno|waist|hip|med|sbp|dbp",colnames(mesa_cov_HT),names = T)]%>%
  magrittr::set_colnames(c("idno","WAIST","HIP","HTMed","SBP","DBP"))%>%
  mutate(WHratio = WAIST/HIP,
         HT = case_when(SBP>130|DBP>80|HTMed==1~1,TRUE~0))
mesa_cov_lung<- mesa_cov_HT%>%inner_join(mesa_cov_demgrh,by="idno") 
mesa_extended <- readr::read_tsv("Data/MESA_combinedtraits.txt",col_names =T)
mesadata <- read.csv("Data/ruitong-mesa-full-20200423.csv",header = T,stringsAsFactors = F)%>%
  filter(!grepl("Changed the",idno))%>%
  mutate_at(vars(idno),as.integer)%>%
  left_join(mesa_cov_lung,by="idno")%>%
  left_join(mesa_extended,by="idno")%>%
  mutate(FVCpp = FVCpp/100,FEV1pp = FEV1pp/100)%>%
  mutate(WHratio = case_when(is.infinite(WHratio)~1,TRUE~WHratio))
model_compare_mesa <- function(D,predictor){
  D1 <- na.omit(D[,c("idno","HT",predictor,"Age","Sex","Race","BMI","WHratio","SmokingStatus","PkYr")])
  expr <- paste("HT~",predictor,
                "+ Age + as.factor(Sex) + as.factor(Race) + BMI + WHratio + as.factor(SmokingStatus) + PkYr",sep="")
  glmodel <-  glm(as.formula(expr),binomial(logit), data = D1, x=TRUE)
  cut_glm <- Youden_max_cut(D = data.frame(prob = predict(object = glmodel, type = 'response'),
                                           label = D1$HT))
  ahi_expr <- "HT ~ AHI3 + Age + as.factor(Sex) + as.factor(Race) + BMI + WHratio + PkYr + as.factor(SmokingStatus)"
  D2 <- na.omit(D[,c("idno","HT","AHI3","Age","Sex","Race","BMI","WHratio","SmokingStatus","PkYr",predictor)])
  ahi_model <- glm(as.formula(ahi_expr),binomial(logit), data = D2, x=TRUE)
  glmodel_zoom <-  glm(as.formula(expr),binomial(logit), data = D2, x=TRUE)
  cat_nri <- nribin(mdl.std = ahi_model,
                    mdl.new =glmodel_zoom, 
                    cut =cut_glm,
                    updown = 'category',niter=0)
  CI <- exp(confint(glmodel)[2,])%>%round(.,2)
  return(list(stats = data.frame(N = nrow(D1),
                                 OR =exp(coef(summary(glmodel))[2,1])%>%round(.,2),
                                 coef = coef(summary(glmodel))[2,1],
                                 coef_low = confint(glmodel)[2,1],
                                 coef_high = confint(glmodel)[2,2],
                                 l_OR = CI[1],
                                 u_OR = CI[2],
                                 pval = coef(summary(glmodel))[2,4],
                                 AUC = ModelMetrics::auc(glmodel)%>%round(.,2),
                                 cat_nri = cat_nri$nri$Estimate[1]%>%round(.,2)),
              model = glmodel,
              newD = D1,
              cat_nri = cat_nri$nri)
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
AHI_model <- model_compare_mesa(D = mesadata,predictor="AHI3")

message("Average AUC for imputed FVCpp is ",round(sapply(1:10,function(s){
  D_impu <- complete(MESA_PSG$PSGIMPU, action = "long", include = TRUE)%>%
    filter(.imp==s)%>%cbind(mesaD_PSG%>%select(idno,HT,Age,Sex,Race,BMI,WHR,PkYr,SmokingStatus))%>%
    na.omit()
  glm(as.formula(
    "HT ~ FVCpp + Age + as.factor(Sex) + as.factor(Race) + BMI + WHR + PkYr + as.factor(SmokingStatus)"), 
    family=binomial(link="logit"),data=D_impu)%>%auc(.)
})%>%mean(),2))
message("Average AUC for imputed FEV1pp is ",round(sapply(1:10,function(s){
  D_impu <- complete(MESA_PSG$PSGIMPU, action = "long", include = TRUE)%>%
    filter(.imp==s)%>%cbind(mesaD_PSG%>%select(idno,HT,Age,Sex,Race,BMI,WHR,PkYr,SmokingStatus))%>%
    na.omit()
  glm(as.formula(
    "HT ~ FEV1pp + Age + as.factor(Sex) + as.factor(Race) + BMI + WHR + PkYr + as.factor(SmokingStatus)"), 
    family=binomial(link="logit"),data=D_impu)%>%auc(.)
})%>%mean(),2))
message("Average AUC for imputed cSP PSG is ",round(sapply(1:10,function(s){
  D_impu <- complete(MESA_PSG$PSGIMPU, action = "long", include = TRUE)%>%
    filter(.imp==s)%>%cbind(mesaD_PSG%>%select(idno,HT,Age,Sex,Race,BMI,WHR,PkYr,SmokingStatus))%>%
    na.omit()
  glm(as.formula(
    "HT ~ CombSDB + Age + as.factor(Sex) + as.factor(Race) + BMI + WHR + PkYr + as.factor(SmokingStatus)"), 
    family=binomial(link="logit"),data=D_impu)%>%auc(.)
})%>%mean(),2))
message("Average AUC for imputed cSP HST is ",round(sapply(1:10,function(s){
  D_impu <- complete(MESA_HST$PSGIMPU, action = "long", include = TRUE)%>%
    filter(.imp==s)%>%cbind(mesaD_PSG%>%select(idno,HT,Age,Sex,Race,BMI,WHR,PkYr,SmokingStatus))%>%
    na.omit()
  glm(as.formula(
    "HT ~ CombSDB + Age + as.factor(Sex) + as.factor(Race) + BMI + WHR + PkYr + as.factor(SmokingStatus)"), 
    family=binomial(link="logit"),data=D_impu)%>%auc(.)
})%>%mean(),2))

psg_model <- with(data = mice::cbind(MESA_PSG$PSGIMPU,mesaD_PSG[,.(HT,Age,Sex,Race,BMI,WHR,PkYr,SmokingStatus)]), 
                  exp = glm(as.formula(
                    "HT ~ CombSDB + Age + as.factor(Sex) + as.factor(Race) + BMI + WHR + PkYr + as.factor(SmokingStatus)"), 
                    family=binomial(link="logit")))%>%
  pool()

hst_model <- with(data = mice::cbind(MESA_HST$PSGIMPU,mesaD_HST[,.(HT,Age,Sex,Race,BMI,WHR,PkYr,SmokingStatus)]), 
                  exp = glm(as.formula(
                    "HT ~ CombSDB + Age + as.factor(Sex) + as.factor(Race) + BMI + WHR + PkYr + as.factor(SmokingStatus)"), 
                    family=binomial(link="logit")))%>%
  pool()

summary(psg_model, conf.int = TRUE, exponentiate = TRUE)[2,c(2,7,8)]%>%round(.,2)
summary(hst_model, conf.int = TRUE, exponentiate = TRUE)[2,c(2,7,8)]%>%round(.,2)

nri_mesa_psg <- sapply(1:10,function(s){
         D_impu <- complete(MESA_PSG$PSGIMPU, action = "long", include = TRUE)%>%
           filter(.imp==s)%>%
           cbind(mesaD_PSG%>%select(idno,HT,Age,Sex,Race,BMI,WHR,PkYr,SmokingStatus))%>%
           rename(WHratio = WHR)%>%
           inner_join(mesadata%>%select(idno,AHI3),by="idno")%>%
           na.omit()
         r <- model_compare_mesa(D = D_impu,predictor = "CombSDB")
         return(r$stats$cat_nri)
       })%>%mean()
nri_mesa_hst <- sapply(1:10,function(s){
  D_impu <- complete(MESA_HST$PSGIMPU, action = "long", include = TRUE)%>%
    filter(.imp==s)%>%
    cbind(mesaD_HST%>%select(idno,HT,Age,Sex,Race,BMI,WHR,PkYr,SmokingStatus))%>%
    rename(WHratio = WHR)%>%
    inner_join(mesadata%>%select(idno,AHI3),by="idno")%>%
    na.omit()
  r <- model_compare_mesa(D = D_impu,predictor = "CombSDB")
  return(r$stats$cat_nri)
})%>%mean()


sapply(1:10,function(s){
  D_impu <- complete(MESA_PSG_nolung$PSGIMPU, action = "long", include = TRUE)%>%
    filter(.imp==s)%>%cbind(MESA_PSG_nolung$OriDat%>%select(idno,HT,Age,Sex,Race,BMI,WHratio,PkYr,SmokingStatus))%>%
    na.omit()
  glm(as.formula(
    "HT ~ CombSDB + Age + as.factor(Sex) + as.factor(Race) + BMI + WHratio + PkYr + as.factor(SmokingStatus)"), 
    family=binomial(link="logit"),data=D_impu)%>%auc(.)
})%>%mean()

sapply(1:10,function(s){
  D_impu <- complete(MESA_HST_nolung$PSGIMPU, action = "long", include = TRUE)%>%
    filter(.imp==s)%>%cbind(MESA_HST_nolung$OriDat%>%select(idno,HT,Age,Sex,Race,BMI,WHratio,PkYr,SmokingStatus))%>%
    na.omit()
  glm(as.formula(
    "HT ~ CombSDB + Age + as.factor(Sex) + as.factor(Race) + BMI + WHratio + PkYr + as.factor(SmokingStatus)"), 
    family=binomial(link="logit"),data=D_impu)%>%auc(.)
})%>%mean()

psg_model_nolung <- with(data = mice::cbind(MESA_PSG_nolung$PSGIMPU,
                                            MESA_PSG_nolung$OriDat[,.(HT,Age,Sex,Race,BMI,WHratio,PkYr,SmokingStatus)]), 
                  exp = glm(as.formula(
                    "HT ~ CombSDB + Age + as.factor(Sex) + as.factor(Race) + BMI + WHratio + PkYr + as.factor(SmokingStatus)"), 
                    family=binomial(link="logit")))%>%
  pool()

hst_model_nolung <- with(data = mice::cbind(MESA_HST_nolung$PSGIMPU,
                                            MESA_HST_nolung$OriDat[,.(HT,Age,Sex,Race,BMI,WHratio,PkYr,SmokingStatus)]), 
                  exp = glm(as.formula(
                    "HT ~ CombSDB + Age + as.factor(Sex) + as.factor(Race) + BMI + WHratio + PkYr + as.factor(SmokingStatus)"), 
                    family=binomial(link="logit")))%>%
  pool()

summary(psg_model_nolung, conf.int = TRUE, exponentiate = TRUE)[2,c(2,7,8)]%>%round(.,2)
summary(hst_model_nolung, conf.int = TRUE, exponentiate = TRUE)[2,c(2,7,8)]%>%round(.,2)
summary(psg_model_nolung)[2,]
summary(hst_model_nolung)[2,]

nri_mesa_psg_nolung <- sapply(1:10,function(s){
  D_impu <- complete(MESA_PSG_nolung$PSGIMPU, action = "long", include = TRUE)%>%
    filter(.imp==s)%>%
    cbind(MESA_PSG_nolung$OriDat%>%select(idno,HT,Age,Sex,Race,BMI,WHratio,PkYr,SmokingStatus))%>%
    inner_join(mesadata%>%select(idno,AHI3),by="idno")%>%
    na.omit()
  r <- model_compare_mesa(D = D_impu,predictor = "CombSDB")
  return(r$stats$cat_nri)
})%>%mean()
nri_mesa_hst_nolung <- sapply(1:10,function(s){
  D_impu <- complete(MESA_HST_nolung$PSGIMPU, action = "long", include = TRUE)%>%
    filter(.imp==s)%>%
    cbind(MESA_HST_nolung$OriDat%>%select(idno,HT,Age,Sex,Race,BMI,WHratio,PkYr,SmokingStatus))%>%
    inner_join(mesadata%>%select(idno,AHI3),by="idno")%>%
    na.omit()
  r <- model_compare_mesa(D = D_impu,predictor = "CombSDB")
  return(r$stats$cat_nri)
})%>%mean()

fev1pp_model <-  with(data = mice::cbind(MESA_PSG$PSGIMPU,mesaD_HST[,.(HT,Age,Sex,Race,BMI,WHR,PkYr,SmokingStatus)]), 
                      exp = glm(as.formula(
                        "HT ~ FEV1pp + Age + as.factor(Sex) + as.factor(Race) + BMI + WHR + PkYr + as.factor(SmokingStatus)"), 
                        family=binomial(link="logit")))%>%
  pool()
fvcpp_model <-  with(data = mice::cbind(MESA_PSG$PSGIMPU,mesaD_HST[,.(HT,Age,Sex,Race,BMI,WHR,PkYr,SmokingStatus)]), 
                      exp = glm(as.formula(
                        "HT ~ FVCpp + Age + as.factor(Sex) + as.factor(Race) + BMI + WHR + PkYr + as.factor(SmokingStatus)"), 
                        family=binomial(link="logit")))%>%
  pool()
summary(fev1pp_model)[2,]
summary(fvcpp_model)[2,]

nri_FEV1pp <- sapply(1:10,function(s){
  D_impu <- complete(MESA_PSG$PSGIMPU, action = "long", include = TRUE)%>%
    filter(.imp==s)%>%
    cbind(MESA_PSG$OriDat%>%select(idno,HT,Age,Sex,Race,BMI,WHratio,PkYr,SmokingStatus))%>%
    inner_join(mesadata%>%select(idno,AHI3),by="idno")%>%
    na.omit()
  r <- model_compare_mesa(D = D_impu,predictor = "FEV1pp")
  return(r$stats$cat_nri)
})%>%mean()
nri_FVCpp <- sapply(1:10,function(s){
  D_impu <- complete(MESA_PSG$PSGIMPU, action = "long", include = TRUE)%>%
    filter(.imp==s)%>%
    cbind(MESA_PSG$OriDat%>%select(idno,HT,Age,Sex,Race,BMI,WHratio,PkYr,SmokingStatus))%>%
    inner_join(mesadata%>%select(idno,AHI3),by="idno")%>%
    na.omit()
  r <- model_compare_mesa(D = D_impu,predictor = "FVCpp")
  return(r$stats$cat_nri)
})%>%mean()

model_coef <- lapply(list(psg_model,hst_model,psg_model_nolung,hst_model_nolung,fev1pp_model,fvcpp_model),function(s)
  summary(s,conf.int = TRUE, exponentiate = F)[2,c(2,7,8)])%>%
  bind_rows()%>%
  set_colnames(c("coef","low","high"))%>%
  mutate(method = c("PSG","HST","PSG_nolung","HST_nolung","FEV1","FVC"))%>%
  rbind(
      AHI_model$stats%>%select(coef,low = coef_low,high = coef_high)%>%mutate(method = "AHI"))
sd_SHHS <- readRDS('~/Documents/HTN_SDB/Result/model_pheno_sd.rds')
model_coef%>%
  inner_join(sd_SHHS,by="method")%>%
  mutate_at(vars(coef,low,high),function(s){exp(s*.$sd)%>%round(.,2)})%>%
  View()
##### Other statistics involved in main text ##### 
PSG_subset <- complete(MESA_PSG$PSGIMPU, action = "long", include = TRUE)%>%
  inner_join(mesaD_PSG%>%
               mutate(.id= 1:n())%>%
               inner_join(mesa6,by="idno")%>%
               filter(HT!=1)%>%
               mutate(mesa6_HT = case_when(htnmed6c==1|htn6c==1~1,TRUE~0))%>%
               select(.id,mesa6_HT,Age,Sex,Race,BMI,WHR,PkYr,SmokingStatus),
             by=".id")%>%as.mids

with(data = PSG_subset, 
     exp = glm(as.formula(
       "mesa6_HT ~ CombSDB + Age + as.factor(Sex) + as.factor(Race) + BMI + WHR + PkYr + as.factor(SmokingStatus)"), 
       family=binomial(link="logit")))%>%
  pool()%>%summary()%>%.[2,]

HST_subset <- complete(MESA_HST$PSGIMPU, action = "long", include = TRUE)%>%
  inner_join(mesaD_HST%>%
               mutate(.id= 1:n())%>%
               inner_join(mesa6,by="idno")%>%
               filter(HT!=1)%>%
               mutate(mesa6_HT = case_when(htnmed6c==1|htn6c==1~1,TRUE~0))%>%
               select(.id,mesa6_HT,Age,Sex,Race,BMI,WHR,PkYr,SmokingStatus),
             by=".id")%>%as.mids

with(data = HST_subset, 
     exp = glm(as.formula(
       "mesa6_HT ~ CombSDB + Age + as.factor(Sex) + as.factor(Race) + BMI + WHR + PkYr + as.factor(SmokingStatus)"), 
       family=binomial(link="logit")))%>%
  pool()%>%summary()%>%.[2,]

PSG_subset_2vcomb <- complete(MESA_PSG$PSGIMPU, action = "long", include = TRUE)%>%
  inner_join(mesaD_PSG%>%
               mutate(.id= 1:n())%>%
               inner_join(mesa6,by="idno")%>%
               mutate(mesa56_HT = case_when(htnmed6c==1|htn6c==1|HT==1~1,TRUE~0))%>%
               select(.id,mesa56_HT,Age,Sex,Race,BMI,WHR,PkYr,SmokingStatus),
             by=".id")%>%as.mids

with(data = PSG_subset_2vcomb, 
     exp = glm(as.formula(
       "mesa56_HT ~ CombSDB + Age + as.factor(Sex) + as.factor(Race) + BMI + WHR + PkYr + as.factor(SmokingStatus)"), 
       family=binomial(link="logit")))%>%
  pool()%>%summary()%>%.[2,]


HST_subset_2vcomb <- complete(MESA_HST$PSGIMPU, action = "long", include = TRUE)%>%
  inner_join(mesaD_HST%>%
               mutate(.id= 1:n())%>%
               inner_join(mesa6,by="idno")%>%
               mutate(mesa56_HT = case_when(htnmed6c==1|htn6c==1|HT==1~1,TRUE~0))%>%
               select(.id,mesa56_HT,Age,Sex,Race,BMI,WHR,PkYr,SmokingStatus),
             by=".id")%>%as.mids

with(data = HST_subset_2vcomb, 
     exp = glm(as.formula(
       "mesa56_HT ~ CombSDB + Age + as.factor(Sex) + as.factor(Race) + BMI + WHR + PkYr + as.factor(SmokingStatus)"), 
       family=binomial(link="logit")))%>%
  pool()%>%summary()%>%.[2,]

mesa_cov_HT <- fread("Data/mesa_sdbtraits_20190827.txt",header=T,data.table = F)
mesa_cov_demgrh <- fread("Data/mesa_ms5527_20190130.csv",header=T,data.table = F)
mesa_cov_demgrh <- mesa_cov_demgrh[,grepALL("idno|bmi|race|age|gender|ppfvc|ppfev1|pkyrs5c|cig5c",colnames(mesa_cov_demgrh),names = T)]%>%
  magrittr::set_colnames(c("idno","BMI","Race","Age","Sex","FVCpp","FEV1pp","PkYr","SmokingStatus"))
mesa_cov_HT <- mesa_cov_HT[,grepALL("idno|waist|hip|med|sbp|dbp",colnames(mesa_cov_HT),names = T)]%>%
  magrittr::set_colnames(c("idno","WAIST","HIP","HTMed","SBP","DBP"))%>%
  mutate(WHratio = WAIST/HIP,
         HT = case_when(SBP>130|DBP>80|HTMed==1~1,TRUE~0))
mesa_cov_lung<- mesa_cov_HT%>%inner_join(mesa_cov_demgrh,by="idno") 
mesa_extended <- readr::read_tsv("Data/MESA_combinedtraits.txt",col_names =T)
mesadata <- read.csv("Data/ruitong-mesa-full-20200423.csv",header = T,stringsAsFactors = F)%>%
  filter(!grepl("Changed the",idno))%>%
  mutate_at(vars(idno),as.integer)%>%
  left_join(mesa_cov_lung,by="idno")%>%
  left_join(mesa_extended,by="idno")%>%as.data.table()
colnames(mesadata)[grep("WHr",colnames(mesadata))] <- "WHR"
glm(as.formula(
  "HT ~ AHI3 + Age + as.factor(Sex) + as.factor(Race) + BMI + WHR + PkYr + as.factor(SmokingStatus)"), 
  family=binomial(link="logit"),data=mesadata)%>%auc(.)
# Correlation between 2 cSP phenotypes
cor(with(MESA_PSG$PSGIMPU,CombSDB)$analyses%>%bind_cols()%>%rowMeans()%>%as.numeric,
    with(MESA_HST$PSGIMPU,CombSDB)$analyses%>%bind_cols()%>%rowMeans()%>%as.numeric())