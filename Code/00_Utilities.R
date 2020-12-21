require(data.table)
require(glmnet)
require(tableone)
require(kableExtra)
require(ggplot2)
require(ggpubr)
require(mice)
require(VIM)
require(lattice)
require(gdata)
require(nlme)
require(naniar)
require(knitr)
require(cluster)
require(factoextra)
require(magrittr)
require(car)
require(stringr)
require(dplyr)
require(nricens)
require(RColorBrewer)
require(ModelMetrics)


#' Logistic Regression Modeling
#'
#' @param formula The formula for logistic regression
#' @param dat The data for fitting logistic regression 
#'
#' @return fitted model object
#' @export
#'
#' @examples 
fit_model_shhs <- function(formula,dat) {
  fit <- glm(as.formula(formula),data=dat,family = "binomial")
  return(fit)
}

#' Process data within SHHS cohort
#'
#' @param exclusionKEY From all SDB traits available, we applied some exclusion criteria quoted in string
#' @param assothres The threshold for association with Incident Hypertension after adjusting requested covariates
#'
#' @return a list of processed analysis dataset, a set of candidate SDB traits and data dictionary 
#' @export
#'
#' @examples
SHHSprocess <- function(exclusionKEY=">=[ ]{0,1}[245]{1}%",assothres=0.2,misthres = 0.3,Drange="complete",lung=T){
  Dfile <- switch(Drange,"complete"="Data/shhs1-dataset-0.15.0.csv",
         "nonEGG"="Data/SHHS1_combined.csv")
  shhs1 <- fread(Dfile,header = T)
  shhs1[,WHratio:=waist/Hip]
  shhs2 <- fread("Data/shhs2-dataset-0.15.0.csv",
                 header = T)[,.(nsrrid,HTNDerv_s2)]
  colnames(shhs2)<-c("SHHSid","IncHT")
  shhs1Cov <- shhs1[,.(nsrrid,age_s1,
                       gender,race,height,
                       HTNDerv_s1,
                       bmi_s1,WHratio,NECK20,CgPkYr,smokstat_s1)]
  colnames(shhs1Cov) <- c("SHHSid","Age","Sex","Race","Height","BaseHT","BMI","WHratio","NeckGrith","PkYr","SmokingStatus")
  ### Data Merging
  if(Drange=="complete"){
    shhs1SDB <- fread("Data/shhs-data-dictionary-0.15.0-variables.csv",header = T)
    repeated_vars <- shhs1SDB%>%filter(grepl("Sleep Architecture",folder))%>%filter(calculation!="" & !grepl("\\*|\\/|\\+|\\-",calculation))%>%pull(id)
    shhs1SDB <- shhs1SDB[grepl("^Measurements/Polysomnography/",folder)&
                           !grepl("Administrative|Quality|Medical Alert",folder)&
                           !grepl("time",type)&
                           !grepl(exclusionKEY,display_name)&!id %in% c("sleep_latency","slp_time",repeated_vars),.(id,folder,display_name)]
    SDBcol <-unlist(sapply(shhs1SDB$id, function(v) grep(paste("^",v,"$",sep=""),colnames(shhs1),value=T,ignore.case=T)))%>%as.character()
    message(paste("We selected between",length(SDBcol)+2,"SDB traits!")) 
    shhsPSG <- shhs1[,setdiff(c(as.character(SDBcol)),c("FEV1","FVC")),with=F]
  }else{
    shhs1SDB <- fread("Data/SHHS1_combined_variables.csv",header = T)
    shhsPSG <- shhs1[,setdiff(colnames(shhs1),c("nsrrid","age_s1",
                                                        "gender","race","height",
                                                        "HTNDerv_s1","FEV1","FVC",
                                                        "bmi_s1","WHratio","NECK20","CgPkYr","smokstat_s1","waist","Hip")),with=F]
    message(paste("We selected between",ncol(shhsPSG)+2,"SDB traits!")) }
  
  if(lung){
    shhsPLM <- shhs1[,c("FEV1","FVC"),with=F]
    shhspp <- cbind(shhs1Cov,shhsPLM)
    ### Update percent predicted FEV1 FVC 
    ppcoef <- get(load("Data/NHANES3.RData"))
    shhspp[,agethres:=ifelse(Age<20,"s","l")]
    shhspp[,Sex:=ifelse(Sex==1,"male","female")]
    shhspp[,Race:=ifelse(Race==1,"white",ifelse(Race==2,"black","other"))]
    shhspp[,Age2:=Age^2][,Height2:=Height^2]
    ### First calculate the predicted FVC and then calculate the percent predicted 
    shhspp$FVCpp <- merge(shhspp,ppcoef[Measure=="FVC",],by=c("Race","agethres","Sex"),all.x=T,sort=F)[
      ,FVCpp:=Intercept+Age*Age_coef+Age2*Age2_coef+Height2*Height2_coef][,.(FVCpp)]
    shhspp[,FVCpp:=FVC/FVCpp]
    shhspp$FEV1pp <- merge(shhspp,ppcoef[Measure=="FEV1",],by=c("Race","agethres","Sex"),all.x=T,sort=F)[
      ,FEV1pp:=Intercept+Age*Age_coef+Age2*Age2_coef+Height2*Height2_coef][,.(FEV1pp)]
    shhspp[,FEV1pp:=FEV1/FEV1pp]
    shhspp[,FEV1FVCr:=FEV1/FVC]
    shhsbase <- cbind(shhspp,shhsPSG)
    shhsbase[,c("Height","agethres","Age2","Height2","FEV1","FVC"):=NULL] 
  }else{
    shhsbase <- cbind(shhs1Cov,shhsPSG)
  }
  shhs <- merge(shhsbase,shhs2,by="SHHSid")
  shhsHP <- shhs[BaseHT==0,]
  shhsHP[,c("Sex","Race","IncHT","SmokingStatus"):=lapply(shhsHP[,.(Sex,Race,IncHT,SmokingStatus)],as.factor)]
  if(lung){
    SDBcol <- if(Drange=="complete"){c("FVCpp","FEV1pp","FEV1FVCr",setdiff(c(as.character(SDBcol)),c("FEV1","FVC")))}else{
      c("FVCpp","FEV1pp","FEV1FVCr",colnames(shhsPSG))
    }
  }else{SDBcol <- colnames(shhsPSG)}
  fitList <- lapply(SDBcol,function(v) 
    fit_model_shhs(paste("IncHT~",v,
                         "+ Age + as.factor(Sex) + as.factor(Race) + BMI + WHratio + NeckGrith + as.factor(SmokingStatus) + PkYr",sep=""),shhsHP))
  singeasso <- sapply(fitList,function(v) coef(summary(v))[2,4])
  ASSOtest <- data.table("Trait"=SDBcol,"singeasso"=singeasso)
  candiSDB <- ASSOtest[singeasso < assothres,Trait]
  message(paste("The candidates for penalized regression is of length",length(candiSDB),"after association threshold!"))
  misscheck <- apply(is.na(shhsHP),2,sum)
  remov <- names(misscheck[misscheck>dim(shhsHP)[1]*misthres])
  candifine <- candiSDB[!candiSDB %in% remov]
  shhs <- na.omit(shhsHP[,c("SHHSid","IncHT","Age","Sex","Race","BMI","WHratio","NeckGrith","PkYr","SmokingStatus",
                          candifine),with=F])
  whole_shhs <- shhsHP[,c("SHHSid","IncHT","Age","Sex","Race","BMI","WHratio","NeckGrith","PkYr","SmokingStatus",
                          candifine),with=F]
  message(paste("The sample size we used for running penalized regression is",dim(shhs)[1],"!"))
  message(paste("The candidates for penalized regression is of length",length(candifine),"after missingness check!"))
  #message(paste("The sample size we used for running penalized regression is",dim(shhsHP)[1]))
  shhs_combinedSDB <- rbind(data.table(id=c("fev1pp","fvcpp"),
                              folder = rep("Measurements/Lung functions!",2),
                              display_name=c("Percent Predicted Forced Expiratory Volume (SHHS1)","Percent Predicted Forced Vital Capacity (SHHS1)")),
                   shhs1SDB)
  return(list(Dat=shhs,candi = candifine,transdat = shhs_combinedSDB,wholeD = whole_shhs))
}

#' Apply Penalized Regression with seed 123, package cv.glmnet
#'
#' @param Xdata The data containing candidate traits and covariates
#' @param Ydata The outcome variable
#' @param Alp The penalized factor, use 1 for lasso regression
#'
#' @return A list of selected traits and achieved AUC
#' @export
#'
#' @examples
SHHSglmnet <- function(Xdata,Ydata,Alp){
  set.seed(123)
  penalCUT <- grep("SmokingStatus2",colnames(Xdata))
  glmmod <- cv.glmnet(x=Xdata,
                      y = Ydata,
                      family ="binomial",
                      intercept=F,
                      alpha=Alp, type.measure="auc",
                      penalty.factor=c(rep(0,penalCUT),rep(1,length((penalCUT+1):dim(Xdata)[2]))))
  tmp_coeffs <- as.matrix(coef(glmmod, s = "lambda.min"))
  coefdt <- data.table("Trait"=rownames(tmp_coeffs),"Effect"= tmp_coeffs[,1])
  sig <- coefdt[Effect!=0,]
  maxauc <- max(glmmod$cvm)
  return(list(COEF = sig,maxAUC = maxauc))}

#' Winsorization Function
#'
#' @param x The numeric vector to be winsorized
#' @param quantile The quantile to be kept, quantile > 0.5 in our case, apply 1-quantile when not
#' @param direction The direction for left and right winsorization
#'
#' @return The winsorized numeric vector 
#' @export
#'
#' @examples
winsorize_quantile <- function(x, quantile, direction){
  stopifnot(is.element(direction, c('left', "right")))
  if (direction == "right"){
    return(pmin(x, quantile(x, quantile)))
  }else{ return(pmax(x, quantile(x, 1-quantile)))}
}
winsorize_examination <- function(D,vec,DR){
  D[,vec,with=F] %>%
    tidyr::gather() %>% 
    ggplot(aes(value)) +
    facet_wrap(~ key, scales = "free") +
    geom_histogram(bins = 40,color="black", fill="white")+
    theme(plot.title = element_text(color='black', hjust = 0.4),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", size = 0.25, fill = NA),
          text = element_text(color='black'),
          panel.grid = element_blank(),
          axis.text = element_text(color='black'),
          strip.background = element_rect(colour="black", fill="gray81", 
                                          size=1, linetype="solid"),legend.position = "none")+
    ggsave(paste0("~/Documents/HTN_SDB/Result/histogram_shh_",DR,dateinfo(),".pdf"),width = 15,height = 15)
}

#' Data pre-processing and penalized regression
#'
#' @param DATA The requested dataset
#' @param candipool The pool of candidate SDB traits
#' @param Alpha The penalized factor, 1 for lasso regression
#' @param sdt Logical flag for whether standardizing the SDB traits or not
#' @param misthres The threshold of missingness
#' @param winthorize Logical flag for whether winsorizing the SDB traits or not
#' @param winRight The character vector for variables to be right winsorized 
#' @param winLeft The character vector for variables to be left winsorized 
#'
#' @return A list of achieved AUC, transformed coefficient set, original SHHS dataset, processed SHHS dataset 
#' @export
#'
#' @examples
penalLogit <- function(D_list,Alpha,sdt,winthorize=F,winRight=NULL,winLeft=NULL,DR){
  set.seed(123)
  DATA <- D_list$Dat
  candifine <- D_list$candi
  shhsANA <- DATA[,-1,with=F]
  if(winthorize){
    if(!is.null(winRight)){
      winr <- shhsANA[,as.data.table(apply(.SD,2,function(v) winsorize_quantile(v,quantile=0.995,direction="right"))),.SDcols=winRight]
      shhsANA <- cbind(shhsANA[,setdiff(colnames(shhsANA),c(winRight)),with=F],winr)}
    if(!is.null(winLeft)){
      winl <- shhsANA[,as.data.table(apply(.SD,2,function(v) winsorize_quantile(v,quantile=0.995,direction="left"))),.SDcols=winLeft]
      shhsANA <- cbind(shhsANA[,setdiff(colnames(shhsANA),c(winLeft)),with=F],winl)
    }
    }
  saveRDS(shhsANA,file=paste0("~/Documents/HTN_SDB/Result/WinsorizedDat_",DR,dateinfo(),".rds"))
  
  if (sdt){sdVEC<- shhsANA[,sapply(.SD, function(v) sd(v)),.SDcols = candifine]
  sdDT <- data.table("Trait"=names(sdVEC),"SD" = sdVEC)
  shhsOSA <- shhsANA[,lapply(.SD, function(v) v/sd(v)),.SDcols = candifine] 
  shhsCOV <- shhsANA[,c("IncHT","Age","Sex","Race","BMI","WHratio","NeckGrith","PkYr","SmokingStatus"),with=F]
  penalDat <- cbind(shhsCOV,shhsOSA)}else{penalDat <- shhsANA}
  Ydat <- as.factor(as.vector(t(penalDat$IncHT)))
  Xdat.model <- model.matrix(~.,penalDat[,-1]) 
  glm.model <- SHHSglmnet(Xdata=Xdat.model,Ydata=Ydat,Alp=Alpha)
  if (sdt){SCALED <- merge(glm.model$COEF,sdDT,by="Trait",all.x=T)
  SCALED$SD[is.na(SCALED$SD)] <- 1
  SCALED[,Effect:=Effect/SD]
  COEFfinal <-SCALED[,.(Trait,Effect)]
  }else{COEFfinal <-glm.model$COEF}
  
  covNAMES <- grep(paste(c("Age","Sex","Race","BMI","WHratio","NeckGrith","PkYr","SmokingStatus"),collapse = "|"),
                   COEFfinal$Trait,value =  T)
  finalOUT <- COEFfinal[!Trait %in% covNAMES,]
  finalOUT$Effect <- as.numeric(finalOUT$Effect)
  finalOUT$niki <-unlist(sapply(finalOUT$Trait, function(v) grep(paste("^",v,"$",sep=""),
                                                                 D_list$transdat$id,value=T,ignore.case=T)))
  FINAL <- merge(finalOUT,D_list$transdat,by.x="niki",by.y="id")
  FINAL <- FINAL[order(-abs(Effect)),]
  message("With this combination of SDB traits, we achieved AUC ",round(glm.model$maxAUC,4), " for predicting hypertension!")
  return(list(coef=FINAL,AUC = glm.model$maxAUC,penalizedPool = penalDat))
}


#' Diagonosis Collection (Not used for final analysis)
#'
#' @param penalcomb The penalized regression return from SHHSpenal()
#'
#' @return A series of plotting and regression testing
#' @export
#'
#' @examples
SDBdiag <- function(penalcomb){
  SDBdat <-na.omit(penalcomb$Dat[,c(penalcomb$coef$Trait,"IncHT","Age",
                                    "Sex","Race","BMI","WHratio","NeckGrith"),with=F])
  
  SDBdat$nSDB <- as.matrix(SDBdat[,penalcomb$coef$Trait,with=F])%*%penalcomb$coef$Effect
  
  p1 <- ggplot(SDBdat, aes(x=nSDB)) + 
    geom_histogram(aes(y=..density..), colour="black", fill="white")+
    geom_density(alpha=.2, fill="#FF6666") 
  
  p2 <- ggplot(SDBdat, aes(x=IncHT, y=nSDB, color=IncHT)) +
    geom_boxplot()
  
  pall <- ggarrange(p1, p2, widths = c(5,5), heights = c(8,8))
  print(pall)
  
  ksTEST <- ks.test(x=SDBdat[IncHT ==1,nSDB], 
                    y=SDBdat[IncHT ==0,nSDB], 
                    alternative = "two.sided")
  print(paste("Two-sample Kolmogorov-Smirnov test has p-value",ksTEST$p))
  
  tTEST <- t.test(x=SDBdat[IncHT ==1,nSDB], 
                  y=SDBdat[IncHT ==0,nSDB], 
                  alternative = "greater")
  
  if (tTEST$p.value<1e-3){
    print(paste("True difference in new SDB trait for people with hypersion versus those without is greater than 0"))}
  
  print(summary(fit_model_shhs("IncHT ~ nSDB + Age + as.factor(Sex) + as.factor(Race)+BMI+ WHratio + NeckGrith",
                               SDBdat))$coef)
}
  
#' Grep sequentially (An extension for self-usage)
#'
#' @param trap The keywords connected by |
#' @param tgt The target character vector
#' @param names Whether to keep the names or only the index
#' @param INV Invert the grep or not 
#'
#' @return The grep output
#' @export
#'
#' @examples
grepALL <- function(trap,tgt,names=F,INV=F){
  vec <- unlist(strsplit(trap,split="|",fixed = T))
  greprlt <- sapply(vec,function(v) grep(v,tgt,value=names,ignore.case = T,invert=INV))
  return(greprlt)
}

#' Process the dataset within MESA accordingly
#'
#' @param penalcomb The penalized result from SHHS cohort
#' @param impuNO The imputation number required
#'
#' @return A list of original MESA dataset, missingness plot and imputed dataset
#' @export
#'
#' @examples
MESAprocess <- function(penalcomb,impuNO){
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
    mutate(FVCpp = FVCpp/100,FEV1pp = FEV1pp/100)%>%as.data.table()
  penalcomb$coef$mesaNAME <-unlist(sapply(penalcomb$coef$Trait, function(v) grep(paste0("^",v,"$"),colnames(mesadata),value=T,ignore.case=T)))
  mesaANA <- mesadata[,c("idno","Age","Sex","Race","BMI","WHratio","HT","PkYr","SmokingStatus",penalcomb$coef$mesaNAME),with=F]
  PSGnames <- grep("FVCpp|FEV1pp",penalcomb$coef$mesaNAME,value=T,invert=T)
  mesaPSG <- mesaANA[!rowSums(is.na(mesaANA[ ,PSGnames,with=F]))==length(PSGnames),]
  mesaPSG <- mesaPSG[!is.na(Race),]
  missdata <- mesaPSG[,colSums(is.na(mesaPSG))!=0,with=F]
  missingpat <- md.pattern(missdata,plot=F)
  ### Perform Imputation over traits to construct 
  PSGimpu <- mesaPSG[,penalcomb$coef$mesaNAME,with=F]
  invisible(capture.output(impurlt <- mice(PSGimpu,seed = 123, m=impuNO)))
  long <- complete(impurlt, action='long', include=TRUE)
  long$CombSDB <- as.numeric(as.matrix(long[,penalcomb$coef$mesaNAME])%*% as.matrix(penalcomb$coef$Effect))
  PSG.imput<- as.mids(long)
  return(list(OriDat=mesaPSG, mispat=missingpat,PSGIMPU=PSG.imput))}

#' Modeling within MESA based imputed dataset
#'
#' @param DAT.PSG The imputed dataset
#' @param DAT.cov The covariates for modeling
#' @param x The exposure
#' @param y The outcome
#' @param modelSTR The modeling string for as.formula()
#' @param randSTR The random covariate string
#' @param causeGEX The logical flag for whether gene expression treated as exposure or outcome
#'
#' @return The list of modeling results for all considered probes
#' @export
#'
#' @examples
midsMODEL <- function(DAT.PSG,DAT.cov,x,y,modelSTR,randSTR=NULL,causeGEX=T){
  message(paste("For Gene Expression Transcript",x,sep=" "))
  covs <- as.character(sapply(strsplit(modelSTR,split=" + ",fixed = T)[[1]],
                              function(v) ifelse(grepl("as.factor",v),
                                                 gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", v, perl=T),v)))
  if(causeGEX){OTHERcov <- DAT.cov[,c(covs,y,randSTR),with=F]}else{OTHERcov <- DAT.cov[,c(covs,x,randSTR),with=F]}
  ANAdat <- mice::cbind(DAT.PSG,OTHERcov)
  lmFML <- paste(paste(y,x,sep=" ~ "),modelSTR,sep=" + ")
  if (is.null(randSTR)){fitMOD <- with(data = ANAdat,exp=lm(formula=as.formula(lmFML)))}else{
    randFML <- as.formula(paste("~1|",randSTR,sep=""))
    fitMOD <- tryCatch(
      {with(data = ANAdat,exp=lme(as.formula(lmFML),random=randFML,na.action = na.omit))},
      warning = function(Warn){
        print(paste("MY_WARNING:  ",Warn))
        fit <- with(data = ANAdat,exp=lm(formula=as.formula(lmFML)))
        return (fit)}, 
      error = function(err) 
      {print(paste("MY_ERROR:  ",err))
        fit = with(data = ANAdat,exp=lm(formula=as.formula(lmFML)))
        return (fit)})}
  suppressWarnings(ModPool <- as.numeric(summary(pool(fitMOD))[2,]))
  names(ModPool)<- c("Estimate","SD","Statistics","df","P.value")
  return(ModPool)}

#' Process midsMODEL() result
#'
#' @param modelLIST The returned model list from midsMODEL()
#' @param probeVEC The probe vector 
#' @param fileNAME The filename for data storage
#'
#' @return
#' @export
#'
#' @examples
modelPROC <- function(modelLIST,probeVEC,fileNAME){
  ModDT <- data.table(do.call(rbind, modelLIST))
  ModDT$Probe <- probeVEC
  ModDT[,fdr:=p.adjust(P.value,method = "BH")]
  save(ModDT,file=fileNAME)}

#' Simply run for getting the date information when running the program
dateinfo<- function(v=Sys.Date()){
  datinfo <- strsplit(as.character(v),split="-")[[1]] 
  datinfo <- paste(substr(datinfo[1],3,4),datinfo[2],datinfo[3],sep="")
  return(datinfo)}

#' Dendrogram Annotation
#'
#' @param n The dendrogram
#' @param clusMember The number of clustered to be considered
#' @param labelColors The labelling colors 
#'
#' @return The annotated dendrogram
#' @export
#'
#' @examples
colLab <- function(n,clusMember,labelColors) {
  if (is.leaf(n)) {
    a <- attributes(n)
    labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
    attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
  }
  n
}

#' Hierarchical Clustering
#'
#' @param k The requested number of clusters
#' @param labelColors The labelling colors of resulted clusters
#' @param data The dataset to be clustered
#' @param jpgname The name for resulted dendrogram
#'
#' @return The way of data being clustered
#' @export
#'
#' @examples
daisyHIER <- function(k,labelColors,data,jpgname){
  DIST <- daisy(data,metric="gower")
  clusters <- hclust(DIST, method="complete")
  clusMember <- cutree(clusters, k)
  clusDendro <- dendrapply(as.dendrogram(clusters), 
                           FUN = colLab,
                           clusMember=clusMember,labelColors=labelColors)
  jpeg(paste('./Result/',jpgname,"_",dateinfo(),".jpg",sep=""),
       width = 4, height = 5, units = 'in', res = 400)
  par(cex=0.4)
  plot(clusDendro)
  n <- nrow(data)
  MidPoint = (clusters$height[n-k] + clusters$height[n-k+1]) / 2
  abline(h = MidPoint, lty=2)
  dev.off()
  return(clusMember)}

#' Format Tableone Object
#'
#' @param tableOne TableOne object
#'
#' @return formatted kable output
#' @export
#'
#' @examples
tableOneformat <- function(tableOne){
  tab1 <- print(tableOne, printToggle = FALSE, noSpaces = T)
  tab1 <- as.data.frame(tab1)
  tab1$p <- as.character(tab1$p)
  rownames(tab1)[1] <- "\\rowcolor{gray!6}Sample Size N"
  tab1$test <- sapply(tab1$p, function(v) if(v==""){" "}else if(as.character(v)=="<0.001"){"***"}else if(
    as.numeric(as.character(v))>=0.001 & as.numeric(as.character(v))<0.01){"**"
  }else if(as.numeric(as.character(v))>=0.01 & as.numeric(as.character(v))<0.05){"*"}else{" "})
  pmsub <- grep("SD",rownames(tab1))
  rownames(tab1) <- sapply(rownames(tab1), function(v){
    cont <- grepl("mean",v)
    if (cont){v <- stringr::str_trim(strsplit(v,split="[(]")[[1]][1],side="right")}
    return(v)})
  tab1[c("0","1","p")] <- sapply(tab1[c("0","1","p")],as.character)
  tab1[,1][pmsub] <- sapply(tab1[,1][pmsub],function(v){
    v <- as.character(v)
    cont <- grepl("[(]",v)
    if(cont){
      mean <- stringr::str_trim(strsplit(v,split="[(]")[[1]][1],side="right")
      sd <- regmatches(v, gregexpr("(?<=\\().*?(?=\\))", v, perl=T))[[1]]
      v <- capture.output(cat(mean,"$\\pm$",sd,sep=""))}
    return(v)})
  tab1[,2][pmsub] <- sapply(tab1[,2][pmsub],function(v){
    v <- as.character(v)
    cont <- grepl("[(]",v)
    if(cont){
      mean <- stringr::str_trim(strsplit(v,split="[(]")[[1]][1],side="right")
      sd <- regmatches(v, gregexpr("(?<=\\().*?(?=\\))", v, perl=T))[[1]]
      v <- capture.output(cat(mean,"$\\pm$",sd,sep=""))}
    return(v)})
  rownames(tab1) <- gsub("%","\\\\%",rownames(tab1))
  rownames(tab1) <- gsub("_","\\\\_",rownames(tab1))
  rownames(tab1)[which(tab1$p=="")[-1]] <- sapply(rownames(tab1)[which(tab1$p=="")[-1]],function(v)
    paste("\\hspace{3mm}",trimws(v,which = "both"),sep=""))
  rownames(tab1)[which(tab1$`0`=="")] <- sapply(rownames(tab1)[which(tab1$`0`=="")],function(v)
    paste("\\rowcolor{gray!6}",trimws(v,which = "both"),sep=""))
  tab1Kable <- kable(tab1, format = "latex",booktabs = T, escape = F,linesep = "")%>%
    kable_styling(latex_options = c("hold_position"))
  return(tab1Kable)}

tableOneCSV <- function(tableOne,savefile){
  table_print <- print(tableOne,
                       exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
  table_print[grepl("SD",rownames(table_print)),1] <- table_print[grepl("SD",rownames(table_print)),1]%>%
    gsub(" \\("," ± ",.)%>%
    gsub("\\)","",.)
  table_print[grepl("SD",rownames(table_print)),2] <- table_print[grepl("SD",rownames(table_print)),2]%>%
    gsub(" \\("," ± ",.)%>%
    gsub("\\)","",.)
  rownames(table_print) <- gsub(" \\(mean \\(SD\\)\\)","",rownames(table_print))
  write.csv(table_print, file =savefile)
}