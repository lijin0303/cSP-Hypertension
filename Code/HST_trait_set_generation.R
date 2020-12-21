require(readr)
require(dplyr)
require(magrittr)
SHHS_data_path <- "~/Downloads/" # Please substitute this with the path to your SHHS1 datasets
setwd(SHHS_data_path)
#===== Pair traits for combination =====
SHHS_vars <- read_csv(paste0("https://sleepdata.org/datasets/shhs/files/m/browser/datasets",
                             "/shhs-data-dictionary-0.15.0-variables.csv"),
                      col_names=T)
SHHS_dat <- read_csv("shhs1-dataset-0.15.0.csv",col_names=T)
# Keep SDB traits from polysomnography  
SHHS_measuresPool <- SHHS_vars%>%
  filter(grepl("^Measurements/Polysomnography/",folder)&
           !grepl("Administrative|Quality|Arousals|Medical Alert|Sleep Architec",folder))
sleepArchi <- SHHS_vars%>%filter(folder=="Measurements/Polysomnography/Sleep Architecture")
# Remove oxygen desaturation 2/4/5 and any traits with arousal 
SHHS_measuresPool%<>%
  filter(!grepl(">=[ ]{0,1}[245]{1}%|arousal",display_name))%>%
  dplyr::select(folder,id,display_name,units,calculation)
Not2combine <- c(SHHS_measuresPool%>%
  filter(grepl("Oxygen Saturation",folder)&(!grepl("oxygen desaturation",display_name)|grepl("Number",display_name)))%>%
  pull(id),
  SHHS_measuresPool%>%
    filter(grepl("Indexes",folder))%>%
    pull(id))
tocombine <- SHHS_measuresPool%>%
  filter(!id %in% Not2combine)
combine_pairs <- tocombine%>%
  mutate(display_name = gsub("supine ","supine,",display_name))%>%
  mutate(remphase = gsub(".*sleep \\((.*)\\),.*","\\1",display_name))%>%
  mutate(merged = gsub(", ",",",gsub(")","",gsub(" \\(.*REM\\)","",display_name))))%>%
  select(-display_name)%>%tidyr::spread(remphase,id)
commonString <- function(s){
  s <- s%>%as.character()
  splitted_common <- Reduce(pmatch,strsplit(s, split = ""))%>%na.omit()%>%as.numeric()
  putBack <- paste0(sapply(splitted_common,function(i) substr(s[1],start=i,stop=i)),collapse = "")
  return(putBack)
}
nameMap <- function(v,ref){
  o <- grep(paste0("^",v,"$") ,ref,ignore.case = T,value = T)%>%as.character()
  return(o)
}
combine_pairs%<>%
  mutate(calculation = case_when(grepl("Minimum",merged)~"pmin",
                                 grepl("Maximum",merged)~"pmax",
                                 grepl("Number of",merged)~"add",
                                grepl("Average",merged)&grepl("per|percent",units)~"*EventTime",
                                grepl("events per hour",units)~"EventsNo/Time",
                                 grepl("Average.*length",merged)~"*EventsNo"),
         UnifiedName = purrr::pmap_chr(list(NREM,REM),~commonString(c(...))),
         NREM_DatCol = purrr::pmap_chr(list(NREM),nameMap,ref=colnames(SHHS_dat)),
         REM_DatCol = purrr::pmap_chr(list(REM),nameMap,ref=colnames(SHHS_dat)))%>%
  arrange(calculation)
#===== Perform calculations =====
min_calc <- pmin(SHHS_dat[,(combine_pairs%>%filter(calculation=="pmin")%>%pull(REM_DatCol))]%>%as.data.frame(),
     SHHS_dat[,(combine_pairs%>%filter(calculation=="pmin")%>%pull(NREM_DatCol))]%>%as.data.frame(),na.rm=T)%>%
  magrittr::set_colnames(combine_pairs%>%filter(calculation=="pmin")%>%pull(UnifiedName))

max_calc <- pmax(SHHS_dat[,(combine_pairs%>%filter(calculation=="pmax")%>%pull(REM_DatCol))]%>%as.data.frame(),
                 SHHS_dat[,(combine_pairs%>%filter(calculation=="pmax")%>%pull(NREM_DatCol))]%>%as.data.frame(),na.rm=T)%>%
  magrittr::set_colnames(combine_pairs%>%filter(calculation=="pmax")%>%pull(UnifiedName))

addition_calc <- (SHHS_dat[,(combine_pairs%>%filter(calculation=="add")%>%pull(REM_DatCol))]%>%
  as.data.frame()%>%mutate_all(~replace(., is.na(.), 0))+
  SHHS_dat[,(combine_pairs%>%filter(calculation=="add")%>%pull(NREM_DatCol))]%>%
  as.data.frame()%>%mutate_all(~replace(., is.na(.), 0)))%>%
  magrittr::set_colnames(combine_pairs%>%filter(calculation=="add")%>%pull(UnifiedName))%>%
  mutate_all(~replace(., .==0, NA))

avgT_Variable <- combine_pairs%>%
  filter(calculation=="*EventsNo")%>%
  mutate(situ = gsub(",","_",gsub("Average | length|>=| oxygen desaturation[s]{0,1}| Apnea|-","",merged)))%>%
  select(situ,NREM_avgT = NREM_DatCol,REM_avgT = REM_DatCol,combName= UnifiedName)%>%arrange(situ)
respiraEvents_Variable <- combine_pairs%>%
  filter(calculation=="add")%>%
  mutate(situ = gsub(",","_",gsub("Number of |>=| oxygen desaturation[s]{0,1}| Apnea|-","",merged)))%>%
  select(situ,NREM_event = NREM_DatCol,REM_event = REM_DatCol)%>%arrange(situ)
avg_calcinfo <- avgT_Variable%>%
  inner_join(respiraEvents_Variable,by="situ")%>%
  mutate(posat = gsub("^[^_]*_([^_]*.*)", "\\1", situ))
SHHS_avglen <- SHHS_dat[,avg_calcinfo%>%tidyr::gather(vari,traitname,-combName,-posat,-situ)%>%pull(traitname)]%>%
  mutate_all(~replace(., is.na(.), 0))
avglen_calc <- ((SHHS_avglen[,avg_calcinfo$NREM_avgT]*SHHS_avglen[,avg_calcinfo$NREM_event]+
                   SHHS_avglen[,avg_calcinfo$REM_avgT]*SHHS_avglen[,avg_calcinfo$REM_event])/
                  (ifelse(SHHS_avglen[,avg_calcinfo$NREM_avgT]==0,0,1)*SHHS_avglen[,avg_calcinfo$NREM_event]+
                     ifelse(SHHS_avglen[,avg_calcinfo$REM_avgT]==0,0,1)*SHHS_avglen[,avg_calcinfo$REM_event]))%>%
  magrittr::set_colnames(avg_calcinfo$combName)%>%
  mutate_all(~replace(., .==0, NA))%>%
  mutate_all(~replace(., .=="NaN", NA))

eventPhour <- combine_pairs%>%
  filter(calculation=="EventsNo/Time")%>%
  mutate(situ = gsub(",","_",gsub(" per hour|>=| oxygen desaturation[s]{0,1}| Apnea|-","",merged)))%>%
  select(situ,NREM_phour = NREM_DatCol,REM_phour = REM_DatCol,combName= UnifiedName)%>%arrange(situ)
eventPhour_calc <- eventPhour%>%
  inner_join(respiraEvents_Variable,by="situ")%>%
  mutate(pose = gsub(".*\\_(.*)\\_.*$", "\\1", situ))%>%
  mutate(REM_time = case_when(pose=="Nonsupine"~"REMepOP",TRUE~"REMepBP"),
         nREM_time = case_when(pose=="Nonsupine"~"NREMepOP",TRUE~"NREMepBP"))
SHHS_eventPhour <- SHHS_dat[,c(eventPhour_calc%>%select(-NREM_phour,-REM_phour,-REM_time,-nREM_time)%>%
                              tidyr::gather(vari,traitname,-combName,-pose,-situ)%>%pull(traitname),
                              "REMepOP","NREMepOP","REMepBP","NREMepBP")]%>%
  mutate_all(~replace(., is.na(.), 0))
eventPhour_calced <- (3600*(SHHS_eventPhour[,eventPhour_calc$NREM_event]+SHHS_eventPhour[,eventPhour_calc$REM_event])/
                  (SHHS_eventPhour[,eventPhour_calc$nREM_time]+SHHS_eventPhour[,eventPhour_calc$REM_time]))%>%
  magrittr::set_colnames(eventPhour_calc$combName)%>%
  mutate_all(~replace(., .=="NaN", NA))%>%
  mutate_all(~replace(., is.infinite(.), NA))

timeall <- lapply(unique(avg_calcinfo$posat),function(s){
  infosub <- avg_calcinfo%>%filter(posat==s)
  REMsub <- (SHHS_avglen[,infosub$REM_avgT]*SHHS_avglen[,infosub$REM_event])%>%rowSums()
  nREMsub <- (SHHS_avglen[,infosub$NREM_avgT]*SHHS_avglen[,infosub$NREM_event])%>%rowSums()
  scomb <- data.frame(nREMsub,REMsub)%>%magrittr::set_colnames(sapply(c("nrem","rem"),function(i) paste(s,i,sep="_")))
  return(scomb)
  })%>%bind_cols()
evenTime_pairs <- combine_pairs%>%filter(calculation=="*EventTime")%>%
  mutate(situ = gsub("-| oxygen desaturation[s]{0,1}|>=","",gsub("^[^,]*,([^,]*.*)", "\\1", merged))%>%gsub(",","_",.))%>%
  mutate(situ_nrem = paste(situ,"nrem",sep="_"),
         situ_rem = paste(situ,"rem",sep="_"))
SHHS_evenTime <- SHHS_dat[,union(evenTime_pairs$NREM_DatCol,evenTime_pairs$REM_DatCol)]%>%
  mutate_all(~replace(., is.na(.), 0))
evenTime_calc <- ((SHHS_evenTime[,evenTime_pairs$NREM_DatCol]*timeall[,evenTime_pairs$situ_nrem]+
                     SHHS_evenTime[,evenTime_pairs$REM_DatCol]*timeall[,evenTime_pairs$situ_rem])/
    (ifelse(SHHS_evenTime[,evenTime_pairs$NREM_DatCol]==0,0,1)*timeall[,evenTime_pairs$situ_nrem]+
       ifelse(SHHS_evenTime[,evenTime_pairs$REM_DatCol]==0,0,1)*timeall[,evenTime_pairs$situ_rem]))%>%
  magrittr::set_colnames(evenTime_pairs$UnifiedName)%>%
  mutate_all(~replace(., .==0, NA))%>%
  mutate_all(~replace(., .=="NaN", NA))
combined_traits <- Reduce(cbind,list(min_calc,max_calc,addition_calc,
                                     avglen_calc,eventPhour_calced,evenTime_calc))

openxlsx::write.xlsx(combine_pairs,"variable_profile.xlsx")
openxlsx::write.xlsx(combined_traits,"combined_traits.xlsx")
