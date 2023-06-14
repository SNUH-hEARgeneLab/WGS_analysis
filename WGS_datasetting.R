library(readxl)
library(dplyr)
library(stringr)
library(reshape2)
library(ComplexHeatmap)
library(ggplot2)
library(colorspace)
library(plotly)
library(tidyr)


#엑셀 원본 불러오기(처음에만 쓰고 최대한 사용 덜하기)
original_data=read_excel("C:\\Users\\SNUH\\Desktop\\WGS_project\\★★★WGS연구_유전 및 표현형 정리_v5.3.xlsx",sheet=1)

####유전자 리스트####
#제외 명단
except <- c('negative','pending','CND','cCMV','B)CND','L)CND')

#유전자 리스트 생성
gene_list <- unique(na.omit(original_data[,15]))
gene_list <- as.data.frame(gene_list %>% filter(!진단유전자 %in% except))
gene_list <- gene_list[,1]

#ploting을 위한 데이터 추출(위에서 지정한 유전자들로만 구성)
data <- original_data[,c(15,27,7:14,29,30,33:40,31,32,41,42,43)] #원본 엑셀의 column이 변하지 않는이상 불변
data <- as.data.frame(data %>% filter(진단유전자 %in% gene_list))

#####1.mutation 분류 및 count####
#mutation을 이어붙인 dataset 생성
data_mut1 <- data[,c(1,3,4,5)]
names(data_mut1) <- c('진단유전자','variant','variant_type_nt','variant_type_aa')
data_mut2 <- data[,c(1,6,7,8)]
names(data_mut2) <- c('진단유전자','variant','variant_type_nt','variant_type_aa')
data_mut <- rbind(data_mut1,data_mut2)
rm(data_mut1,data_mut2)

#Gene Mutation count ranking(Decreasing sorting)
data_mut <- data_mut[!is.na(data_mut$variant),]
variant_rank <- data_mut[,1:2] %>% 
  group_by(`진단유전자`) %>% 
  summarize(v1=n())
variant_rank <- variant_rank[order(variant_rank$v1,decreasing = T),]
variant_rank_list <- as.data.frame(variant_rank)[,1]

#Variant Counting(sorting by Gene count ranking)
variant_count <- data_mut %>% 
  group_by(`진단유전자`,`variant`) %>% 
  summarize(v1=n())
variant_count$진단유전자 <- factor(variant_count$진단유전자,levels=variant_rank_list)
variant_count <- variant_count[order(variant_count$진단유전자),]

names(variant_count) <- c('진단유전자','variant','count')

data_mut2 <- data_mut[,-1]
variant_type <- distinct(merge(variant_count,data_mut2,by='variant'))
variant_type$진단유전자 <- factor(variant_type$진단유전자,levels=variant_rank_list)
variant_type <- variant_type[order(variant_type$진단유전자),]
variant_type <- variant_type %>% 
  relocate(c(variant),.after='진단유전자')
rownames(variant_type)=NULL
rm(data_mut2)
#variant_type 데이터가 가장 잘 정리된 데이터(WGS유전자정보정리-가계제외)

#variant count 가계별
data_mut1 <- data[,c(1,2,3,4,5)]
names(data_mut1) <- c('진단유전자','병록번호','variant','variant_type_nt','variant_type_aa')
data_mut2 <- data[,c(1,2,6,7,8)]
names(data_mut2) <- c('진단유전자','병록번호','variant','variant_type_nt','variant_type_aa')
vrcountdata <- rbind(data_mut1,data_mut2)
rm(data_mut1,data_mut2)

vrcountdata <- vrcountdata[!is.na(vrcountdata$variant),]
vrcountdata <- distinct(vrcountdata)

variant_count2 <- vrcountdata %>% 
  group_by(`진단유전자`,`variant`) %>% 
  summarize(v1=n())

variant_count2$진단유전자 <- factor(variant_count2$진단유전자,levels=variant_rank_list)
variant_count2 <- variant_count2[order(variant_count2$진단유전자),]
names(variant_count2) <- c('진단유전자','variant','count')

variant_type2 <- distinct(merge(variant_count2,data_mut,by='variant'))
variant_type2 <- variant_type2[,-4]
variant_type2 <- variant_type2 %>% 
  relocate(c(variant),.after='진단유전자.x')
names(variant_type2) <- c("진단유전자","variant","count","variant_type_nt","variant_type_aa")
variant_type2$진단유전자 <- factor(variant_type2$진단유전자,levels=variant_rank_list)
variant_type2 <- variant_type2[order(variant_type2$진단유전자),]
rownames(variant_type2)=NULL
#row수가 179가 나와야함
variant_type2 <- variant_type2[-c(18,19,26,27,67,68),]
rownames(variant_type2)=NULL

####2.가계 Count & SNV Mutation Count####
familycount <- data %>% 
  group_by(`진단유전자`) %>% 
  summarize(v1=n())
familycount <- familycount[order(familycount$v1,decreasing = T),]
family_rank_list <- as.data.frame(familycount)[,1]

familylist <- data[,c(1,2,3,6)]
familylist$진단유전자 <- factor(familylist$진단유전자,levels=family_rank_list)
familylist <- familylist[order(familylist$진단유전자),]
rownames(familylist)=NULL

mutation_count <- cbind(sapply(str_split(variant_count$variant,";"),"[",1),variant_count$count)
mutation_count <- as.data.frame(mutation_count)
mutation_count$V1 <- str_sub(mutation_count$V1,-3,-1)
mutation_count$V2 <- as.numeric(mutation_count$V2)
SNVcount <- mutation_count %>% 
  group_by(`V1`) %>% 
  summarize(v1=sum(`V2`))
#리스트 보고 수작업
SNVcount
SNVcount<- SNVcount[c(1:3,5:13),]
rm(mutation_count)

#=====유전자 정보 정리 excel file follow up clear=====#
######oncoPrint data setting#####

#oncodata <- data[,-c(3,6)]
#varianttype 변경 및 통합
onco1 <- data[,c(1,2,4,5,7,8)]

onco1[is.na(onco1)] <- '.'
onco1 <- onco1[,-c(3,5)]

onco1$`variant type_M1_aa` <- gsub(1,"Missense",onco1$`variant type_M1_aa`)
onco1$`variant type_M1_aa` <- gsub(2,"Frameshift",onco1$`variant type_M1_aa`)
onco1$`variant type_M1_aa` <- gsub(3,"Nonsense",onco1$`variant type_M1_aa`)
onco1$`variant type_M1_aa` <- gsub(4,"Inframe_del/dup",onco1$`variant type_M1_aa`)
onco1$`variant type_M1_aa` <- gsub(5,"Splicing",onco1$`variant type_M1_aa`)
onco1$`variant type_M1_aa` <- gsub(6,"SV",onco1$`variant type_M1_aa`)
onco1$`variant type_M1_aa` <- gsub(7,"mtRNA",onco1$`variant type_M1_aa`)

onco1$`variant type_M2_aa` <- gsub(1,"Missense*",onco1$`variant type_M2_aa`)
onco1$`variant type_M2_aa` <- gsub(2,"Frameshift*",onco1$`variant type_M2_aa`)
onco1$`variant type_M2_aa` <- gsub(3,"Nonsense*",onco1$`variant type_M2_aa`)
onco1$`variant type_M2_aa` <- gsub(4,"Inframe_del/dup*",onco1$`variant type_M2_aa`)
onco1$`variant type_M2_aa` <- gsub(5,"Splicing*",onco1$`variant type_M2_aa`)
onco1$`variant type_M2_aa` <- gsub(6,"SV*",onco1$`variant type_M2_aa`)

onco1 <- cbind(onco1,paste(onco1[,3],onco1[,4],sep=";"))
onco1 <- onco1[,c(1,2,5)]
names(onco1) <- c('진단유전자','병록번호','varianttype')
onco1$varianttype <- gsub(".","",onco1$varianttype,fixed = TRUE)

onco2 <- data[,9:25]
onco_fin <- cbind(onco1,onco2)
rm(onco1,onco2)
#initial audio제거
onco_fin <- subset(onco_fin,select=-`initial audio`)

#나이 수정(M으로 표기된)
onco_fin$나이<- gsub("M","000",onco_fin$나이)
onco_fin$나이<- as.numeric(onco_fin$나이)

for(i in 1:nrow(onco_fin)){
  if(onco_fin$나이[i]>48000){
    onco_fin$나이[i] <-  4
  }else if(onco_fin$나이[i]>36000){
    onco_fin$나이[i] <-  3
  }else if(onco_fin$나이[i]>24000){
    onco_fin$나이[i] <-  2
  }else if(onco_fin$나이[i]>12000){
    onco_fin$나이[i] <-  1
  }else if(onco_fin$나이[i]>=1000){
    onco_fin$나이[i] <-  0
  }
}

for(i in 1:nrow(onco_fin)){
  if(as.numeric(onco_fin$나이[i])>=60){
    onco_fin$나이[i] = "60"
  }else if(as.numeric(onco_fin$나이[i])>=50){
    onco_fin$나이[i] = "50" 
  }else if(as.numeric(onco_fin$나이[i])>=40){
    onco_fin$나이[i] = "40" 
  }else if(as.numeric(onco_fin$나이[i])>=30){
    onco_fin$나이[i] = "30"
  }else if(as.numeric(onco_fin$나이[i])>=20){
    onco_fin$나이[i] = "20" 
  }else if(as.numeric(onco_fin$나이[i])>=5){
    onco_fin$나이[i] = "5" 
  }else if(as.numeric(onco_fin$나이[i])>=0){
    onco_fin$나이[i] = "0"
  }
}


onco_fin$진단유전자 <- factor(onco_fin$진단유전자,levels=family_rank_list)
onco_fin <- onco_fin[order(onco_fin$진단유전자),]
#onco_fin <- onco_fin[order(onco_fin$진단유전자),]
rownames(onco_fin)=NULL
#=====onco기본 데이터 세팅 clear=====#
#####!!!oncoPrint ver2-final!!!#####
oncoprintdata2 <- onco_fin[,c(1,2,3,7:18,19)]

oncoprintdata2[is.na(oncoprintdata2)] <- ""

names(oncoprintdata2)[1] <- "Gene"
names(oncoprintdata2)[4] <- "Age_at_Test"


for(i in 1:nrow(oncoprintdata2)){
  if(oncoprintdata2$`hearing loss onset`[i]==1){
    oncoprintdata2$Early_ID[i] = "Early"
    oncoprintdata2$Late_ID[i] = ""
    oncoprintdata2$Adult_onset[i] = ""
  }else if(oncoprintdata2$`hearing loss onset`[i]==2|oncoprintdata2$`hearing loss onset`[i]==3){
    oncoprintdata2$Early_ID[i] = ""
    oncoprintdata2$Late_ID[i] = "Delay"
    oncoprintdata2$Adult_onset[i] = ""
  }else if(oncoprintdata2$`hearing loss onset`[i]==4){
    oncoprintdata2$Early_ID[i] = ""
    oncoprintdata2$Late_ID[i] = ""
    oncoprintdata2$Adult_onset[i] = "Adult"
  }else{
    oncoprintdata2$Early_ID[i] = ""
    oncoprintdata2$Late_ID[i] = ""
    oncoprintdata2$Adult_onset[i] = ""
  }
}

for(i in 1:nrow(oncoprintdata2)){
  if(oncoprintdata2$syndromic_1[i]==2|oncoprintdata2$syndromic_2[i]==2){
    oncoprintdata2$Syndromic[i] = "Syndromic"
  }else{
    oncoprintdata2$Syndromic[i] = ""
  }
}

for(i in 1:nrow(oncoprintdata2)){
  if(oncoprintdata2$`hearing loss type, Rt`[i]==2|oncoprintdata2$`hearing loss type, Lt`[i]==2){
    oncoprintdata2$Mixed[i] = "Mixed"
  }else{
    oncoprintdata2$Mixed[i] = ""
  }
}

for(i in 1:nrow(oncoprintdata2)){
  if(oncoprintdata2$asymmetry[i]==1|oncoprintdata2$asymmetry[i]==2){
    oncoprintdata2$Asymmetric[i] = "Asymmetric"
  }else{
    oncoprintdata2$Asymmetric[i] = ""
  }
}

#Asymmetric과 Mixed는 제외하고 이후 진행
#Severity
for(i in 1:nrow(oncoprintdata2)){
  if(oncoprintdata2$Mixed[i]!="Mixed" && oncoprintdata2$Asymmetric[i]!="Asymmetric"){
    if(oncoprintdata2$`hearing severity, Lt`[i]==1|oncoprintdata2$`hearing severity, Lt`[i]==2|oncoprintdata2$`hearing severity, Rt`[i]==1|oncoprintdata2$`hearing severity, Rt`[i]==2){
      oncoprintdata2$MildMod[i] = "MildMod"
      oncoprintdata2$ModSev[i]=""
      oncoprintdata2$SevProf[i]=""
    }else if(oncoprintdata2$`hearing severity, Lt`[i]==3|oncoprintdata2$`hearing severity, Rt`[i]==3){
      oncoprintdata2$MildMod[i]=""
      oncoprintdata2$ModSev[i] = "ModSev"
      oncoprintdata2$SevProf[i]=""
    }else if(oncoprintdata2$`hearing severity, Lt`[i]==4|oncoprintdata2$`hearing severity, Lt`[i]==5|oncoprintdata2$`hearing severity, Rt`[i]==4|oncoprintdata2$`hearing severity, Rt`[i]==5){
      oncoprintdata2$MildMod[i]=""
      oncoprintdata2$ModSev[i] = ""
      oncoprintdata2$SevProf[i] = "SevProf"
    }else{
      oncoprintdata2$MildMod[i] = ""
      oncoprintdata2$ModSev[i] = ""
      oncoprintdata2$SevProf[i] = ""
    }
  }else{
    oncoprintdata2$MildMod[i] = ""
    oncoprintdata2$ModSev[i] = ""
    oncoprintdata2$SevProf[i] = ""
  }
}

#Configuration
for(i in 1:nrow(oncoprintdata2)){
  if(oncoprintdata2$Mixed[i]!="Mixed" && oncoprintdata2$Asymmetric[i]!="Asymmetric"){
    if(oncoprintdata2$`hearing configuration, Lt`[i]==1|oncoprintdata2$`hearing configuration, Rt`[i]==1){
      oncoprintdata2$Flat[i] = "Flat"
      oncoprintdata2$DownSloping[i] = ""
      oncoprintdata2$Cookie[i] = ""
    }else if(oncoprintdata2$`hearing configuration, Lt`[i]==2|oncoprintdata2$`hearing configuration, Rt`[i]==2|oncoprintdata2$`hearing configuration, Lt`[i]==3|oncoprintdata2$`hearing configuration, Rt`[i]==3){
      oncoprintdata2$Flat[i] = ""
      oncoprintdata2$DownSloping[i] = "DownSlop"
      oncoprintdata2$Cookie[i] = ""
    }else if(oncoprintdata2$`hearing configuration, Lt`[i]==4|oncoprintdata2$`hearing configuration, Rt`[i]==4){
      oncoprintdata2$Flat[i] = ""
      oncoprintdata2$DownSloping[i] = ""
      oncoprintdata2$Cookie[i] = "Cookie"
    }else{
      oncoprintdata2$Flat[i] = ""
      oncoprintdata2$DownSloping[i] = ""
      oncoprintdata2$Cookie[i] = ""
    }
  }else{
    oncoprintdata2$Flat[i] = ""
    oncoprintdata2$DownSloping[i]=""
    oncoprintdata2$Cookie[i]=""
  }
}

#Progressive
for(i in 1:nrow(oncoprintdata2)){
  if(oncoprintdata2$progressive[i]==1){
    oncoprintdata2$Sustan[i] = "Sustan"
    oncoprintdata2$MildNone[i] =""
  }else if(oncoprintdata2$progressive[i]==2){
    oncoprintdata2$Sustan[i] = ""
    oncoprintdata2$MildNone[i] ="Mild"
  }else if(oncoprintdata2$progressive[i]==3){
    oncoprintdata2$Sustan[i]=""
    oncoprintdata2$MildNone[i]="None"
  }else{
    oncoprintdata2$Sustan[i]=""
    oncoprintdata2$MildNone[i]=""
  }
}

#가계3개까지 자르기
oncoprintdata3 <- oncoprintdata2[1:163,]
oncoprintdata3 <- oncoprintdata3[,-16]

data3 <- oncoprintdata3[,c(2,1,16:29)]
data3 <- as.data.frame(t(data3))
names(data3) <- c(1:ncol(data3))
data3 <- data3[-1,]

data4 <- oncoprintdata3[,2:3]
data4 <- as.data.frame(t(data4))
names(data4) <- data4[1,]
data4 <- data4[-1,]
rownames(data4) <- "Variant_Type"

#variant type으로 각각을 표시해주기
for(i in 2:13){
  for(j in 1:ncol(data3)){
    if(data3[i,j]!=""){
      data3[i,j] <- data4[1,j]
    }
    else{
      next
    }
  }
}

col_ord <- colnames(data3)
row_ord <- rownames(data3)

col=c("Missense" = "red", "Frameshift" = "orange","Nonsense"="yellow","Inframe_del/dup"="green","Splicing"="blue","SV"="purple","mtRNA"="black",
      "Missense*" = "red", "Frameshift*" = "orange","Nonsense*"="yellow","Inframe_del/dup*"="green","Splicing*"="blue","SV*"="purple","mtRNA*"="black",
      "GJB2"="lightgoldenrod2","SLC26A4"="lightcyan3","STRC"="lavenderblush2","USH2A"="seagreen4","CDH23"="violetred3","MPZL2"="thistle4",
      "OTOA"="lightyellow","MYO15A"="#63666A","EYA1"="#FFE900","SIX1"="#2AD2C9","WFS1"="#9063CD","ACTG1"="#FCAEBB",
      "COL4A3"="#74D1Ea","LMX1A"="#EA6811","TECTA"="#DDA46F","COL11A1"="#69B3E7","COL1A1"="#E0457B","KCNQ4"="#7C3A2D",
      "MT-TL1"="black","MYO6"="#E4002B","POU4F3"="#0057B8","TMPRSS3"="#ABD156","Sustan"="sienna4","Mild"="#A8dddd","None"="#A8bbbb")

oncoPrint(data3,
          alter_fun = list(
            background = function(x, y, w, h) {
              grid.rect(x,y,w*0.9,h*0.9,
                        gp = gpar(fill = "white", col = "grey"))
            },
            "Missense" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",gp = gpar(fill = col["Missense"], col = NA)),
            "Frameshift" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",gp = gpar(fill = col["Frameshift"], col = NA)),
            "Nonsense" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",gp = gpar(fill = col["Nonsense"], col = NA)),
            "Inframe_del/dup" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",gp = gpar(fill = col["Inframe_del/dup"], col = NA)),
            "SV" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",gp = gpar(fill = col["SV"], col = NA)),
            "Splicing" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",gp = gpar(fill = col["Splicing"], col = NA)),
            "mtRNA" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",gp = gpar(fill = col["mtRNA"], col = NA)),
            "Missense*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="left",gp = gpar(fill = col["Missense*"], col = NA)),
            "Frameshift*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="left",gp = gpar(fill = col["Frameshift*"], col = NA)),
            "Nonsense*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="left",gp = gpar(fill = col["Nonsense*"], col = NA)),
            "Inframe_del/dup*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="left",gp = gpar(fill = col["Inframe_del/dup*"], col = NA)),
            "SV*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="left",gp = gpar(fill = col["SV*"], col = NA)),
            "Splicing*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="left",gp = gpar(fill = col["Splicing*"], col = NA)),
            "mtRNA*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",gp = gpar(fill = col["mtRNA"], col = NA)),
            "GJB2" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["GJB2"], col = NA)),
            "SLC26A4" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["SLC26A4"], col = NA)),
            "STRC" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["STRC"], col = NA)),
            "USH2A" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["USH2A"], col = NA)),
            "CDH23" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["CDH23"], col = NA)),
            "MPZL2" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["MPZL2"], col = NA)),
            "OTOA" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["OTOA"], col = NA)),
            "MYO15A" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["MYO15A"], col = NA)),
            "EYA1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["EYA1"], col = NA)),
            "SIX1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["SIX1"], col = NA)),
            "WFS1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["WFS1"], col = NA)),
            "ACTG1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["ACTG1"], col = NA)),
            "COL4A3" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["COL4A3"], col = NA)),
            "LMX1A" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["LMX1A"], col = NA)),
            "TECTA" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["TECTA"], col = NA)),
            "COL11A1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["COL11A1"], col = NA)),
            "COL1A1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["COL1A1"], col = NA)),
            "KCNQ4" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["KCNQ4"], col = NA)),
            "MT-TL1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["MT-TL1"], col = NA)),
            "MYO6" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["MYO6"], col = NA)),
            "POU4F3" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["POU4F3"], col = NA)),
            "TMPRSS3" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["TMPRSS3"], col = NA)),
            "Sustan" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["Sustan"], col = NA)),
            "Mild" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["Mild"], col = NA)),
            "None" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["None"], col = NA))
          ), col = col,row_names_side = "left",pct_side = "right",show_column_names = FALSE,top_annotation = NULL,right_annotation = NULL,
          alter_fun_is_vectorized = FALSE,column_order = col_ord,row_order=row_ord,show_pct=FALSE)

lgd=Legend(labels=c("GJB2","SLC26A4","STRC","USH2A","CDH23","MPZL2","OTOA","MYO15A","EYA1","SIX1","WFS1","ACTG1",
                    "COL4A3","LMX1A","TECTA","COL11A1","COL1A1","KCNQ4","MT-TL1","MYO6","POU4F3","TMPRSS3",
                    "Missense", "Frameshift","Nonsense","Inframe_del/dup","Splicing","SV","mtRNA",
                    "Sustan","Mild","None"),
           legend_gp=gpar(fill=c("lightgoldenrod2","lightcyan3","lavenderblush2","seagreen4","violetred3","thistle4",
                                 "lightyellow","#63666A","#FFE900","#2AD2C9","#9063CD","#FCAEBB",
                                 "#74D1Ea","#EA6811","#DDA46F","#69B3E7","#E0457B","#7C3A2D",
                                 "black","#E4002B","#0057B8","#ABD156","red","orange","yellow","green","blue","purple","black",
                                 "sienna4","#A8dddd","#A8bbbb")),nrow=1)
dev.new()
draw(lgd)



#####oncoPrint ver1#####
oncoprintdata=dcast(onco_fin,진단유전자~병록번호,value.var = "varianttype")

oncoprintdata <- as.data.frame(oncoprintdata)
row.names(oncoprintdata) <- oncoprintdata[,1]
oncoprintdata <- oncoprintdata[,-1]
oncoprintdata[is.na(oncoprintdata)] <- ""
oncoprintdata <- oncoprintdata[,familylist$병록번호]
names(oncoprintdata)
oncoprintdata <- subset(oncoprintdata,select=-`76573947.1`)
oncoprintdata <- subset(oncoprintdata,select=-`80678315.1`)

data2 <- as.matrix(oncoprintdata)

#유전자 ranker순
col = c("Missense" = "red", "Frameshift" = "orange","Nonsense"="yellow","Inframe_del/dup"="green","Splicing"="blue","SV"="purple","mtRNA"="black",
        "Missense*" = "red", "Frameshift*" = "orange","Nonsense*"="yellow","Inframe_del/dup*"="green","Splicing*"="blue","SV*"="purple","mtRNA*"="black")
col_ord=colnames(data2)

oncoPrint(data2,
          alter_fun = list(
            background = function(x, y, w, h) {
              grid.rect(x,y,w*0.9,h*0.9,
                        gp = gpar(fill = "white", col = "grey"))
            },
            "Missense" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",
                                                        gp = gpar(fill = col["Missense"], col = NA)),
            "Frameshift" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",
                                                          gp = gpar(fill = col["Frameshift"], col = NA)),
            "Nonsense" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",
                                                        gp = gpar(fill = col["Nonsense"], col = NA)),
            "Inframe_del/dup" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",
                                                               gp = gpar(fill = col["Inframe_del/dup"], col = NA)),
            "SV" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",
                                                  gp = gpar(fill = col["SV"], col = NA)),
            "Splicing" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",
                                                        gp = gpar(fill = col["Splicing"], col = NA)),
            "mtRNA" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",
                                                     gp = gpar(fill = col["mtRNA"], col = NA)),
            "Missense*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="left",
                                                         gp = gpar(fill = col["Missense*"], col = NA)),
            "Frameshift*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="left",
                                                           gp = gpar(fill = col["Frameshift*"], col = NA)),
            "Nonsense*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="left",
                                                         gp = gpar(fill = col["Nonsense*"], col = NA)),
            "Inframe_del/dup*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="left",
                                                                gp = gpar(fill = col["Inframe_del/dup*"], col = NA)),
            "SV*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="left",
                                                   gp = gpar(fill = col["SV*"], col = NA)),
            "Splicing*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="left",
                                                         gp = gpar(fill = col["Splicing*"], col = NA)),
            "mtRNA*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",
                                                      gp = gpar(fill = col["mtRNA"], col = NA))
          ), col = col,row_names_side = "left",pct_side = "right",top_annotation = NULL,right_annotation = NULL,column_order = col_ord,
          row_order=family_rank_list,show_column_names = FALSE,show_pct=FALSE
)

lgd=Legend(labels=c("Missense","Frameshift","Nonsense","Inframe_del/dup","Splicing","SV","mtRNA"),
           legend_gp=gpar(fill=c("red","orange","yellow","green","blue","purple","black")),nrow=1)
dev.new()
draw(lgd)

#oncoplot_sub
##한방에 그리기##
oncosubdata <- oncoprintdata2[,c(2,1,17:30)]
oncosubdata <- as.data.frame(t(oncosubdata))
names(oncosubdata) <- oncosubdata[1,]
oncosubdata <- oncosubdata[-1,]
#double primary 제거
oncosubdata <- oncosubdata[,-c(166,154)]

oncosubdata2 <- oncoprintdata2[,2:3]
oncosubdata2 <- as.data.frame(t(oncosubdata2))
names(oncosubdata2) <- oncosubdata2[1,]
oncosubdata2 <- oncosubdata2[-1,]
rownames(oncosubdata2) <- "Variant_Type"
oncosubdata2 <- oncosubdata2[,-c(166,154)]

#variant type으로 각각을 표시해주기
for(i in 2:13){
  for(j in 1:ncol(oncosubdata)){
    if(oncosubdata[i,j]!=""){
      oncosubdata[i,j] <- oncosubdata2[1,j]
    }
    else{
      next
    }
  }
}
rm(oncosubdata2)
oncosubdata <- oncosubdata[-1,]
data2 <- as.data.frame(data2)
oncosubdata <- as.data.frame(oncosubdata)
test <- rbind(data2,oncosubdata)
row_ord <- rownames(test)
col_ord <- names(test)


col=c("Missense" = "red", "Frameshift" = "orange","Nonsense"="yellow","Inframe_del/dup"="green","Splicing"="blue","SV"="purple","mtRNA"="black",
      "Missense*" = "red", "Frameshift*" = "orange","Nonsense*"="yellow","Inframe_del/dup*"="green","Splicing*"="blue","SV*"="purple","mtRNA*"="black",
      "Sustan"="sienna4","Mild"="#A8dddd","None"="#A8bbbb")


oncoPrint(test,
          alter_fun = list(
            background = function(x, y, w, h) {
              grid.rect(x,y,w*0.9,h*0.9,
                        gp = gpar(fill = "white", col = "grey"))
            },
            "Missense" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",gp = gpar(fill = col["Missense"], col = NA)),
            "Frameshift" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",gp = gpar(fill = col["Frameshift"], col = NA)),
            "Nonsense" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",gp = gpar(fill = col["Nonsense"], col = NA)),
            "Inframe_del/dup" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",gp = gpar(fill = col["Inframe_del/dup"], col = NA)),
            "SV" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",gp = gpar(fill = col["SV"], col = NA)),
            "Splicing" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",gp = gpar(fill = col["Splicing"], col = NA)),
            "mtRNA" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",gp = gpar(fill = col["mtRNA"], col = NA)),
            "Missense*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="left",gp = gpar(fill = col["Missense*"], col = NA)),
            "Frameshift*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="left",gp = gpar(fill = col["Frameshift*"], col = NA)),
            "Nonsense*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="left",gp = gpar(fill = col["Nonsense*"], col = NA)),
            "Inframe_del/dup*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="left",gp = gpar(fill = col["Inframe_del/dup*"], col = NA)),
            "SV*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="left",gp = gpar(fill = col["SV*"], col = NA)),
            "Splicing*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="left",gp = gpar(fill = col["Splicing*"], col = NA)),
            "mtRNA*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",gp = gpar(fill = col["mtRNA"], col = NA)),
            "Sustan" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["Sustan"], col = NA)),
            "Mild" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["Mild"], col = NA)),
            "None" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["None"], col = NA))
          ), col = col,row_names_side = "left",pct_side = "right",show_column_names = FALSE,top_annotation = NULL,right_annotation = NULL,
          alter_fun_is_vectorized = FALSE,column_order = col_ord,row_order=row_ord,show_pct=FALSE)


#####figure b,c,d,e,f barplot#####
figuredata <- original_data[,c(15,31,32,34:40,42)]
figuredata <- as.data.frame(figuredata)
figuredata <- figuredata[!is.na(figuredata$진단유전자),]
figuredata[is.na(figuredata)] <- ""
rownames(figuredata)=NULL

#figure B를 위한 데이터셋
figb <- figuredata[figuredata$syndromic_1==2,]
figb <- figuredata[figuredata$syndromic_2==2,]

figbcount<- figb %>% 
  group_by(`진단유전자`) %>% 
  summarize(v1=n())
figbcount <- figbcount[order(figbcount$v1,decreasing = T),]

total <- sum(figbcount$v1)
total2 <- figbcount %>% filter(!진단유전자 %in% except)
total2 <- sum(total2$v1)
total1 <- c('Affected N',total)
total2 <- c('Affected N',total2)
total2 <- as.data.frame(t(total2))
names(total2) <- c('진단유전자','v1')
figbcount <- rbind(total1,total2,figbcount)
figbcount <- as.data.frame(figbcount %>% filter(!진단유전자 %in% except))
names(figbcount) <- c('진단유전자','count')

figbcount$count <- as.numeric(figbcount$count)
figbcount <- figbcount[order(figbcount$count,decreasing = T),]
figbcount <- figbcount[c(1:10),]
rownames(figbcount)=NULL

pal <- c('8'='#FCAEBB','6'='#FCAEBB','4'='#FCAEBB','5'='#FCAEBB','3'='#FCAEBB','2'='#FCAEBB','38'='#FCAEBB','52'='#86D8A4')
ggplot(figbcount)+
  geom_bar(stat="identity",position="identity",mapping=aes(x=진단유전자,y=count,fill=as.factor(count)))+
  scale_x_discrete(limits=c('','','','',rev(figbcount$진단유전자)))+
  theme_classic()+
  scale_fill_manual(values=pal)+
  coord_flip()+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank())

rm(figb,figbcount)

#figure(c)를 위한 데이터셋
figc <- figuredata[figuredata$`hearing loss type, Lt`==2|figuredata$`hearing loss type, Rt`==2,]
rownames(figc)=NULL

figccount<- figc %>% 
  group_by(`진단유전자`) %>% 
  summarize(v1=n())
figccount <- figccount[order(figccount$v1,decreasing = T),]

total <- c('Affected N',sum(figccount$v1))
total2 <- figccount %>% filter(!진단유전자 %in% except)
total2 <- sum(total2$v1)
figccount <- rbind(total,figccount)
figccount$v1 <- as.numeric(figccount$v1)
figccount <- as.data.frame(figccount %>% filter(!진단유전자 %in% except))
figccount <- left_join(figccount,familycount,by='진단유전자')
figccount[1,3] <- figccount[1,2]
figccount[1,2] <- total2
names(figccount) <- c('진단유전자','count','total')
figccount <- figccount[order(figccount$count,figccount$total,decreasing = T),]
#KCNQ4,CRB1제외
rownames(figccount)=NULL
figccount <- figccount[-c(9,11),]
figccount
figccount[1,2] <- 18
figccount[1,3] <- 25

test <- figccount %>% gather(a,n,total,count)
test[test$진단유전자=='Affected N',2][1] <- 'total2'

pal <- c('total'='#86D8A4','count'='#FCAEBB','total2'='#86D8A4')
ggplot(test)+
  geom_bar(stat="identity",position="identity",mapping=aes(x=진단유전자,y=n,fill=a))+
  scale_x_discrete(limits=c('','','','',rev(figccount$진단유전자)))+
  theme_classic()+
  scale_fill_manual(values=pal)+
  coord_flip()+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank())

rm(figc,figccount)

#figure(d)를 위한 데이터셋
figd1 <- figuredata[figuredata$asymmetry==1|figuredata$asymmetry==2,]
rownames(figd1 )=NULL

figd1count<- figd1  %>% 
  group_by(`진단유전자`) %>% 
  summarize(v1=n())
figd1count <- figd1count[order(figd1count$v1,decreasing = T),]

total <- c('Affected N',sum(figd1count$v1))
total2 <- figd1count %>% filter(!진단유전자 %in% except)
total2 <- sum(total2$v1)
figd1count <- rbind(total,figd1count)
figd1count$v1 <- as.numeric(figd1count$v1)
figd1count <- as.data.frame(figd1count %>% filter(!진단유전자 %in% except))
figd1count <- left_join(figd1count,familycount,by='진단유전자')
figd1count[1,3] <- figd1count[1,2]
figd1count[1,2] <- total2
names(figd1count) <- c('진단유전자','count','total')
figd1count <- figd1count[order(figd1count$count,figd1count$total,decreasing = T),]

test <- figd1count %>% gather(a,n,total,count)
test[test$진단유전자=='Affected N',2][1] <- 'total2'

pal <- c('total'='#86D8A4','count'='#FCAEBB','total2'='#86D8A4')
ggplot(test)+
  geom_bar(stat="identity",position="identity",mapping=aes(x=진단유전자,y=n,fill=a))+
  scale_x_discrete(limits=c('','','',rev(figd1count$진단유전자)))+
  theme_classic()+
  scale_fill_manual(values=pal)+
  coord_flip()+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank())

rm(figd1,figd1count)

#figure(e,f)를 위한 데이터셋
figeef <- figuredata[figuredata$`hearing loss type, Lt`!=2 & figuredata$`hearing loss type, Rt`!=2,]
figeef <- figeef[figeef$asymmetry==3,]
rownames(figeef)=NULL

#figure(f)
#f1
figf_1 <- figeef[figeef$`hearing configuration, Rt`==1 | figeef$`hearing configuration, Lt`==1,]
figf1count<- figf_1  %>% 
  group_by(`진단유전자`) %>% 
  summarize(v1=n())
figf1count <- figf1count[order(figf1count$v1,decreasing = T),]

total <- c('Affected N',sum(figf1count$v1))
total2 <- figf1count %>% filter(!진단유전자 %in% except)
total2 <- sum(total2$v1)
figf1count <- rbind(total,figf1count)
figf1count$v1 <- as.numeric(figf1count$v1)
figf1count <- as.data.frame(figf1count %>% filter(!진단유전자 %in% except))
figf1count <- left_join(figf1count,familycount,by='진단유전자')
figf1count[1,3] <- figf1count[1,2]
figf1count[1,2] <- total2
names(figf1count) <- c('진단유전자','count','total')
figf1count <- figf1count[order(figf1count$count,figf1count$total,decreasing = T),]
figf1count2 <- figf1count[figf1count$count>1,]

test <- figf1count2 %>% gather(a,n,total,count)
test[test$진단유전자=='Affected N',2][1] <- 'total2'

pal <- c('total'='#86D8A4','count'='#FCAEBB','total2'='#86D8A4')
ggplot(test)+
  geom_bar(stat="identity",position="identity",mapping=aes(x=진단유전자,y=n,fill=a))+
  scale_x_discrete(limits=c('','',c(rev(figf1count2$진단유전자))))+
  theme_classic()+
  scale_fill_manual(values=pal)+
  coord_flip()+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank())

rm(figf_1,figf1count,figf1count2)


#f2
figf_2 <- figeef[figeef$`hearing configuration, Rt`==2 | figeef$`hearing configuration, Lt`==2 | figeef$`hearing configuration, Rt`==3 | figeef$`hearing configuration, Lt`==3,]
figf2count<- figf_2  %>% 
  group_by(`진단유전자`) %>% 
  summarize(v1=n())
figf2count <- figf2count[order(figf2count$v1,decreasing = T),]

total <- c('Affected N',sum(figf2count$v1))
total2 <- figf2count %>% filter(!진단유전자 %in% except)
total2 <- sum(total2$v1)
figf2count <- rbind(total,figf2count)
figf2count$v1 <- as.numeric(figf2count$v1)
figf2count <- as.data.frame(figf2count %>% filter(!진단유전자 %in% except))
figf2count <- left_join(figf2count,familycount,by='진단유전자')
figf2count[1,3] <- figf2count[1,2]
figf2count[1,2] <- total2
names(figf2count) <- c('진단유전자','count','total')
figf2count <- figf2count[order(figf2count$count,figf2count$total,decreasing = T),]
figf2count2 <- figf2count[figf2count$count>1,]
figf2count2 <- figf2count2[1:12,]

test <- figf2count2 %>% gather(a,n,total,count)
test[test$진단유전자=='Affected N',2][1] <- 'total2'

pal <- c('total'='#86D8A4','count'='#FCAEBB','total2'='#86D8A4')
ggplot(test)+
  geom_bar(stat="identity",position="identity",mapping=aes(x=진단유전자,y=n,fill=a))+
  scale_x_discrete(limits=c(rev(figf2count2$진단유전자)))+
  theme_classic()+
  scale_fill_manual(values=pal)+
  coord_flip()+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank())

rm(figf_2,figf2count,figf2count2)

#f4
figf_4 <- figeef[figeef$`hearing configuration, Rt`==4 | figeef$`hearing configuration, Lt`==4,]
figf4count<- figf_4  %>% 
  group_by(`진단유전자`) %>% 
  summarize(v1=n())
figf4count <- figf4count[order(figf4count$v1,decreasing = T),]

total <- c('Affected N',sum(figf4count$v1))
total2 <- figf4count %>% filter(!진단유전자 %in% except)
total2 <- sum(total2$v1)
figf4count <- rbind(total,figf4count)
figf4count$v1 <- as.numeric(figf4count$v1)
figf4count <- as.data.frame(figf4count %>% filter(!진단유전자 %in% except))
figf4count <- left_join(figf4count,familycount,by='진단유전자')
figf4count[1,3] <- figf4count[1,2]
figf4count[1,2] <- total2
names(figf4count) <- c('진단유전자','count','total')
figf4count <- figf4count[order(figf4count$count,figf4count$total,decreasing = T),]

test <- figf4count %>% gather(a,n,total,count)
test[test$진단유전자=='Affected N',2][1] <- 'total2'

pal <- c('total'='#86D8A4','count'='#FCAEBB','total2'='#86D8A4')
ggplot(test)+
  geom_bar(stat="identity",position="identity",mapping=aes(x=진단유전자,y=n,fill=a))+
  scale_x_discrete(limits=c(rev(figf4count$진단유전자)))+
  theme_classic()+
  scale_fill_manual(values=pal)+
  coord_flip()+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank())

rm(figf_4,figf4count)

#figure(e)
#e1
fige_1 <- figeef[figeef$`hearing severity, Rt`==1 | figeef$`hearing severity, Rt`==2 | figeef$`hearing severity, Lt`==1 | figeef$`hearing severity, Lt`==2,]
fige1count<- fige_1  %>% 
  group_by(`진단유전자`) %>% 
  summarize(v1=n())
fige1count <- fige1count[order(fige1count$v1,decreasing = T),]

total <- c('Affected N',sum(fige1count$v1))
total2 <- fige1count %>% filter(!진단유전자 %in% except)
total2 <- sum(total2$v1)
fige1count <- rbind(total,fige1count)
fige1count$v1 <- as.numeric(fige1count$v1)
fige1count <- as.data.frame(fige1count %>% filter(!진단유전자 %in% except))
fige1count <- left_join(fige1count,familycount,by='진단유전자')
fige1count[1,3] <- fige1count[1,2]
fige1count[1,2] <- total2
names(fige1count) <- c('진단유전자','count','total')
fige1count <- fige1count[order(fige1count$count,fige1count$total,decreasing = T),]
fige1count2 <- fige1count[fige1count$count>1,]

test <- fige1count2 %>% gather(a,n,total,count)
test[test$진단유전자=='Affected N',2][1] <- 'total2'

pal <- c('total'='#86D8A4','count'='#FCAEBB','total2'='#86D8A4')
ggplot(test)+
  geom_bar(stat="identity",position="identity",mapping=aes(x=진단유전자,y=n,fill=a))+
  scale_x_discrete(limits=rev(fige1count2$진단유전자))+
  theme_classic()+
  scale_fill_manual(values=pal)+
  coord_flip()+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank())

rm(fige_1,fige1count,fige1count2)

#e2
fige_2 <- figeef[figeef$`hearing severity, Rt`==3 | figeef$`hearing severity, Lt`==3,]
fige2count<- fige_2  %>% 
  group_by(`진단유전자`) %>% 
  summarize(v1=n())
fige2count <- fige2count[order(fige2count$v1,decreasing = T),]

total <- c('Affected N',sum(fige2count$v1))
total2 <- fige2count %>% filter(!진단유전자 %in% except)
total2 <- sum(total2$v1)
fige2count <- rbind(total,fige2count)
fige2count$v1 <- as.numeric(fige2count$v1)
fige2count <- as.data.frame(fige2count %>% filter(!진단유전자 %in% except))
fige2count <- left_join(fige2count,familycount,by='진단유전자')
fige2count[1,3] <- fige2count[1,2]
fige2count[1,2] <- total2
names(fige2count) <- c('진단유전자','count','total')
fige2count <- fige2count[order(fige2count$count,fige2count$total,decreasing = T),]
fige2count2 <- fige2count[fige2count$count>1,]

test <- fige2count2 %>% gather(a,n,total,count)
test[test$진단유전자=='Affected N',2][1] <- 'total2'

pal <- c('total'='#86D8A4','count'='#FCAEBB','total2'='#86D8A4')
ggplot(test)+
  geom_bar(stat="identity",position="identity",mapping=aes(x=진단유전자,y=n,fill=a))+
  scale_x_discrete(limits=c('','',rev(fige2count2$진단유전자)))+
  theme_classic()+
  scale_fill_manual(values=pal)+
  coord_flip()+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank())

rm(fige_2,fige2count,fige2count2)

#e3
fige_3 <- figeef[figeef$`hearing severity, Rt`==4 | figeef$`hearing severity, Rt`==5 | figeef$`hearing severity, Lt`==4 | figeef$`hearing severity, Lt`==5,]
fige3count<- fige_3  %>% 
  group_by(`진단유전자`) %>% 
  summarize(v1=n())
fige3count <- fige3count[order(fige3count$v1,decreasing = T),]

total <- c('Affected N',sum(fige3count$v1))
total2 <- fige3count %>% filter(!진단유전자 %in% except)
total2 <- sum(total2$v1)
fige3count <- rbind(total,fige3count)
fige3count$v1 <- as.numeric(fige3count$v1)
fige3count <- as.data.frame(fige3count %>% filter(!진단유전자 %in% except))
fige3count <- left_join(fige3count,familycount,by='진단유전자')
fige3count[1,3] <- fige3count[1,2]
fige3count[1,2] <- total2
names(fige3count) <- c('진단유전자','count','total')
fige3count <- fige3count[order(fige3count$count,fige3count$total,decreasing = T),]
fige3count2 <- fige3count[fige3count$count>1,]

test <- fige3count2 %>% gather(a,n,total,count)
test[test$진단유전자=='Affected N',2][1] <- 'total2'

pal <- c('total'='#86D8A4','count'='#FCAEBB','total2'='#86D8A4')
ggplot(test)+
  geom_bar(stat="identity",position="identity",mapping=aes(x=진단유전자,y=n,fill=a))+
  scale_x_discrete(limits=c('','',rev(fige3count2$진단유전자)))+
  theme_classic()+
  scale_fill_manual(values=pal)+
  coord_flip()+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank())

rm(fige_3,fige3count,fige3count2)

#figure(g)
#g1
figg <- figuredata[figuredata$progressive==1,]

figgcount<- figg %>% 
  group_by(`진단유전자`) %>% 
  summarize(v1=n())
figgcount <- figgcount[order(figgcount$v1,decreasing = T),]

total <- sum(figgcount$v1)
total2 <- figgcount %>% filter(!진단유전자 %in% except)
total2 <- sum(total2$v1)
total1 <- c('Affected N',total)
total2 <- c('Affected N',total2)
total2 <- as.data.frame(t(total2))
names(total2) <- c('진단유전자','v1')
figgcount <- rbind(total1,total2,figgcount)
figgcount <- as.data.frame(figgcount %>% filter(!진단유전자 %in% except))
names(figgcount) <- c('진단유전자','count')

figgcount$count <- as.numeric(figgcount$count)
figgcount <- figgcount[order(figgcount$count,decreasing = T),]
figgcount <- figgcount[c(1:16),]
figgcount <- figgcount[-1,]
rownames(figgcount)=NULL

pal <- c('15'='#FCAEBB','6'='#FCAEBB','5'='#FCAEBB','9'='#FCAEBB','3'='#FCAEBB','2'='#FCAEBB','78'='#FCAEBB','83'='#86D8A4')
ggplot(figgcount)+
  geom_bar(stat="identity",position="identity",mapping=aes(x=진단유전자,y=count,fill=as.factor(count)))+
  scale_x_discrete(limits=c(rev(figgcount$진단유전자)))+
  theme_classic()+
  scale_fill_manual(values=pal)+
  coord_flip()+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank())

rm(figg,figgcount)

#g2
figg <- figuredata[figuredata$progressive==2|figuredata$progressive==3,]

figgcount<- figg %>% 
  group_by(`진단유전자`) %>% 
  summarize(v1=n())
figgcount <- figgcount[order(figgcount$v1,decreasing = T),]

total <- sum(figgcount$v1)
total2 <- figgcount %>% filter(!진단유전자 %in% except)
total2 <- sum(total2$v1)
total1 <- c('Affected N',total)
total2 <- c('Affected N',total2)
total2 <- as.data.frame(t(total2))
names(total2) <- c('진단유전자','v1')
figgcount <- rbind(total1,total2,figgcount)
figgcount <- as.data.frame(figgcount %>% filter(!진단유전자 %in% except))
names(figgcount) <- c('진단유전자','count')

figgcount$count <- as.numeric(figgcount$count)
figgcount <- figgcount[order(figgcount$count,decreasing = T),]
figgcount <- figgcount[c(1:19),]
figgcount <- figgcount[-1,]
rownames(figgcount)=NULL

pal <- c('14'='#FCAEBB','13'='#FCAEBB','6'='#FCAEBB','5'='#FCAEBB','8'='#FCAEBB','3'='#FCAEBB','2'='#FCAEBB','81'='#FCAEBB','84'='#86D8A4')
ggplot(figgcount)+
  geom_bar(stat="identity",position="identity",mapping=aes(x=진단유전자,y=count,fill=as.factor(count)))+
  scale_x_discrete(limits=c(rev(figgcount$진단유전자)))+
  theme_classic()+
  scale_fill_manual(values=pal)+
  coord_flip()+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank())

rm(figg,figgcount)

#####stepgraph#####
#graph1
step <- original_data[,c(15,17,18,19,20)]
step <- step[!is.na(step$진단유전자),]
step <- na.omit(step)

`step1` <- as.numeric(count(step[step$step1==1 & step$step2==0 & step$step3==0 & step$step4==0,]))
`step1&2` <- as.numeric(count(step[step$step1==1 & step$step2==1 & step$step3==0 & step$step4==0,]))
`step1&2&3` <- as.numeric(count(step[step$step1==1 & step$step2==1 & step$step3==1 & step$step4==0,]))
`step1&2&3&4` <- as.numeric(count(step[step$step1==1 & step$step2==1 & step$step3==1 & step$step4==1,]))
`step1&2&4` <- as.numeric(count(step[step$step1==1 & step$step2==1 & step$step3==0 & step$step4==1,]))
`step2`<- as.numeric(count(step[step$step1==0 & step$step2==1 & step$step3==0 & step$step4==0,]))
`step2&3` <- as.numeric(count(step[step$step1==0 & step$step2==1 & step$step3==1 & step$step4==0,]))
`step2&3&4` <- as.numeric(count(step[step$step1==0 & step$step2==1 & step$step3==1 & step$step4==1,]))
`step2&4` <- as.numeric(count(step[step$step1==0 & step$step2==1 & step$step3==0 & step$step4==1,]))

sum(`step1`,`step1&2`,`step2`,`step1&2&3`,`step2&3`,`step1&2&3&4`,
    `step1&2&4`,`step2&3&4`,`step2&4`)

stepcount <- as.data.frame(rbind(`step1`,`step1&2`,`step2`,`step1&2&3`,`step2&3`,`step1&2&3&4`,
                                 `step1&2&4`,`step2&3&4`,`step2&4`))
names(stepcount) <- "count"
stepcount <- tibble::rownames_to_column(stepcount,var="step")

pal <- c('58'='#86D8A4','191'='#FCAEBB','28'='#FCAEBB','43'='#D1ACF9','1'='#D1ACF9','13'='#51C6E6',
         '37'='#51C6E6','3'='#51C6E6','4'='#51C6E6')
ggplot(stepcount,aes(x=step,y=count,fill=as.factor(count)))+
  geom_bar(stat="identity",position="identity")+
  scale_x_discrete(limits=rev(stepcount$step))+
  theme_classic()+
  geom_text(aes(x=step,y=count,hjust=-0.1, label=count))+
  coord_flip()+
  scale_fill_manual(values=pal)+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank())

#graph2
step123 <-step[step$step1==1 & step$step2==1 & step$step3==1 & step$step4==0,]
step23 <- step[step$step1==0 & step$step2==1 & step$step3==1 & step$step4==0,]
stepgr2 <- rbind(step123,step23)

stepgr2 <- stepgr2 %>% 
  group_by(`진단유전자`) %>% 
  summarize(v1=n())

total <- sum(stepgr2$v1)
stepgr2 <- stepgr2 %>% filter(!진단유전자 %in% except)
total2 <- sum(stepgr2$v1)
total1 <- c('Affected N',total)
total2 <- c('Affected N',total2)
total2 <- as.data.frame(t(total2))
names(total2) <- c('진단유전자','v1')
gr2count <- rbind(total1,total2,stepgr2)
names(gr2count) <- c('진단유전자','count')
gr2count$count <- as.numeric(gr2count$count)
gr2count <- gr2count[order(gr2count$count,decreasing = T),]
rownames(gr2count)=NULL
rownames(gr2count)


pal <- c('15'='#86D8A4','3'='#86D8A4','2'='#86D8A4','1'='#86D8A4','21'='#FCAEBB','44'='#D1ACF9')
ggplot(gr2count,aes(x=진단유전자,y=count,fill=as.factor(count)))+
  scale_x_discrete(limits=rev(gr2count$진단유전자))+
  geom_bar(stat="identity",position="identity")+
  geom_text(aes(x=진단유전자,y=count,hjust=-0.1, label=count))+
  theme_classic()+
  coord_flip()+
  scale_fill_manual(values=pal)+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank())

#step pie chart
plot_ly(stepcount,labels=~step,values=~count,textposition=ifelse(stepcount$count<10,"outside","inside"),
        textinfo='label+value+percent',
        texttemplate='%{label}<br>n=%{value}<br>%{percent}',pull=c(0,0,0,0,0,0.1,0.1,0.1,0.1),
        marker=list(colors=c("#AAC8FC","#6CD6E4","#D1ACF9","#86D8A4","#FFB3B5","#D5C783","lightyellow",
                             "#93328E","#EA6811","#FFE900","#E0457B","#7C3A2D")),sort=FALSE) %>% 
  add_pie(hole=0.5) %>% 
  layout(xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),showlegend=FALSE)


gr2count <- gr2count[3:6,]
plot_ly(gr2count,labels=~진단유전자,values=~count,textinfo='label+value+percent',
        marker=list(colors=c("#AAC8FC","#6CD6E4","#D1ACF9","#86D8A4","#FFB3B5","#D5C783","lightyellow",
                             "#93328E","#EA6811","#FFE900","#E0457B","#7C3A2D")),sort=FALSE) %>% 
  add_pie(hole=0.5) %>% 
  layout(xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),showlegend=FALSE)




#####classification 적용(supp plot)(2개 이상인 환자들 나누기)#####
a <- oncoprintdata2[c(138:141,152:160),]
a$`functional classification` <- c(7,7,7,7,4,4,4,2,2,2,7,7,7)
oncoprintdata2 <- rbind(oncoprintdata2,a)
rm(a)

for(i in 1:nrow(oncoprintdata2)){
  if(oncoprintdata2$`functional classification` [i]=="3,4"|oncoprintdata2$`functional classification` [i] == "3,7"){
    oncoprintdata2$`functional classification` [i] =  "3"
  }else if(oncoprintdata2$`functional classification`[i]=="1,2"){
    oncoprintdata2$`functional classification`[i]="1"
  }
  else{
    next
  }
}

oncoprintdata2$진단유전자 <- factor(oncoprintdata2$Gene,levels=family_rank_list)
oncoprintdata2 <- oncoprintdata2[order(oncoprintdata2$`functional classification`,oncoprintdata2$Gene),]
rownames(oncoprintdata2)=NULL

oncoraregrid <- oncoprintdata2[,c(2,1,16:30)]
oncoraregrid <- as.data.frame(t(oncoraregrid))
names(oncoraregrid) <- c(1:ncol(oncoraregrid))
oncoraregrid <- oncoraregrid[-1,]
rownames(oncoraregrid)[2] <- c("Classification")

for(i in 2:nrow(oncoraregrid)){
  for(j in 1:ncol(oncoraregrid)){
    if(oncoraregrid[i,j]!=""){
      oncoraregrid[i,j] <- oncoraregrid[2,j]
    }
    else{
      next
    }
  }
}

row_ord <- rownames(oncoraregrid)
col_ord <- colnames(oncoraregrid)

col=c("1"="#951123","2"="#ABEEBA","3"="#FFE760","4"="#A8eeee","5"="#EE55FF","6"="#1a5a54","7"="#63666A","8"="#877DDD",
      "GJB2"="lightgoldenrod2","SLC26A4"="lightcyan3","STRC"="lavenderblush2","USH2A"="seagreen4","CDH23"="violetred3","MPZL2"="thistle4",
      "OTOA"="lightyellow","MYO15A"="#63666A","EYA1"="#FFE900","SIX1"="#2AD2C9","WFS1"="#9063CD","ACTG1"="#FCAEBB",
      "COL4A3"="#74D1Ea","LMX1A"="#EA6811","TECTA"="#DDA46F","COL11A1"="#69B3E7","COL1A1"="#E0457B","KCNQ4"="#7C3A2D",
      "MT-TL1"="black","MYO6"="#E4002B","POU4F3"="#0057B8","TMPRSS3"="#ABD156","COCH"="#EE55FF","COL4A5"="#951123",
      "MITF"="#BBEBBE","MYO7A"="burlywood","OTOGL"="#FFE760","POU3F4"="#0055E7","TMC1"="#9994DD","ANKRD11"="#666666",
      "CASP10"="#FEFAAA","CLCNKA"="#aa11aa","CRB1"="#ab7111","DSPP"="#E2325F","ELMOD3"="#EA5411","ESRRB"="#ffee66",
      "EYA4"="green","GATA3"="blue","GSDME"="darkblue","ILDR1"="purple","KIT"="pink","LOXHD1"="skyblue",
      "MT-RNR1"="grey","MYH9"="white","NLRP3"="#a651be","OTOG"="#1a5a54","PDE6B"="#e4a65a","SLC12A3"="#b1c654",
      "SOX10"="#A8eeee","SSBP1"="#f16aaa","TJP2"="sienna4","SLC12A2"="#abcee5",
      "COQ6"="wheat3", "CTCF"="palegoldenrod", "DNAJC3"="khaki4", "MT-ND5"="honeydew3", "RBM10"="azure", "SMAD4"="lightsalmon3", 
      "SPATA5"="seagreen2")

oncoPrint(oncoraregrid,
          alter_fun = list(
            background = function(x, y, w, h) {
              grid.rect(x,y,w*0.9,h*0.9,
                        gp = gpar(fill = "white", col = "grey"))
            },
            "1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["1"], col = NA)),
            "2" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["2"], col = NA)),
            "3" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["3"], col = NA)),
            "4" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["4"], col = NA)),
            "5" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["5"], col = NA)),
            "6" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["6"], col = NA)),
            "7" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["7"], col = NA)),
            "8" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["8"], col = NA)),
            "GJB2" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["GJB2"], col = NA)),
            "SLC26A4" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["SLC26A4"], col = NA)),
            "STRC" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["STRC"], col = NA)),
            "USH2A" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["USH2A"], col = NA)),
            "CDH23" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["CDH23"], col = NA)),
            "MPZL2" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["MPZL2"], col = NA)),
            "OTOA" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["OTOA"], col = NA)),
            "MYO15A" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["MYO15A"], col = NA)),
            "EYA1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["EYA1"], col = NA)),
            "SIX1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["SIX1"], col = NA)),
            "WFS1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["WFS1"], col = NA)),
            "ACTG1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["ACTG1"], col = NA)),
            "COL4A3" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["COL4A3"], col = NA)),
            "LMX1A" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["LMX1A"], col = NA)),
            "TECTA" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["TECTA"], col = NA)),
            "COL11A1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["COL11A1"], col = NA)),
            "COL1A1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["COL1A1"], col = NA)),
            "KCNQ4" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["KCNQ4"], col = NA)),
            "MT-TL1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["MT-TL1"], col = NA)),
            "MYO6" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["MYO6"], col = NA)),
            "POU4F3" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["POU4F3"], col = NA)),
            "TMPRSS3" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["TMPRSS3"], col = NA)),
            "COCH" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["COCH"], col = NA)),
            "COL4A5" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["COL4A5"], col = NA)),
            "MITF" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["MITF"], col = NA)),
            "MYO7A" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["MYO7A"], col = NA)),
            "OTOGL" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["OTOGL"], col = NA)),
            "POU3F4" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["POU3F4"], col = NA)),
            "TMC1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["TMC1"], col = NA)),
            "ANKRD11" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["ANKRD11"], col = NA)),
            "CASP10" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["CASP10"], col = NA)),
            "CLCNKA" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["CLCNKA"], col = NA)),
            "CRB1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["CRB1"], col = NA)),
            "DSPP" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["DSPP"], col = NA)),
            "ELMOD3" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["ELMOD3"], col = NA)),
            "ESRRB" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["ESRRB"], col = NA)),
            "EYA4" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["EYA4"], col = NA)),
            "GATA3" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["GATA3"], col = NA)),
            "GSDME" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["GSDME"], col = NA)),
            "ILDR1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["ILDR1"], col = NA)),
            "KIT" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["KIT"], col = NA)),
            "LOXHD1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["LOXHD1"], col = NA)),
            "MT-RNR1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["MT-RNR1"], col = NA)),
            "MYH9" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["MYH9"], col = NA)),
            "NLRP3" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["NLRP3"], col = NA)),
            "OTOG" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["OTOG"], col = NA)),
            "PDE6B" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["PDE6B"], col = NA)),
            "SLC12A3" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["SLC12A3"], col = NA)),
            "SOX10" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["SOX10"], col = NA)),
            "SSBP1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["SSBP1"], col = NA)),
            "SLC12A2" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["SLC12A2"], col = NA)),
            "TJP2" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["TJP2"], col = NA)),
            "COQ6" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["COQ6"], col = NA)),
            "CTCF" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["CTCF"], col = NA)),
            "DNAJC3" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["DNAJC3"], col = NA)),
            "MT-ND5" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["MT-ND5"], col = NA)),
            "RBM10" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["RBM10"], col = NA)),
            "SMAD4" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["SMAD4"], col = NA)),
            "SPATA5" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["SPATA5"], col = NA))
          ),
          col = col,row_names_side = "left",pct_side = "right",show_column_names = FALSE,top_annotation = NULL,right_annotation = NULL,
          alter_fun_is_vectorized = FALSE,column_order = col_ord,show_pct=FALSE,row_order = row_ord)

lgd=Legend(labels=c("1","2","3","4","5","6","7","8"),
           legend_gp=gpar(fill=c("#951123","#ABEEBA","#FFE760","#A8eeee","#EE55FF","#1a5a54","#63666A","#877DDD")),nrow=1)
dev.new()
draw(lgd)

#######Full Gene########

oncogene <- oncoprintdata2[,c(2,1)]
oncogene <- as.data.frame(t(oncogene))
names(oncogene) <- c(1:ncol(oncogene))
oncogene <- oncogene[-1,]
rownames(oncogene) <- c("Gene")

col=c("GJB2"="lightgoldenrod2","SLC26A4"="lightcyan3","STRC"="lavenderblush2","USH2A"="seagreen4","CDH23"="violetred3","MPZL2"="thistle4",
      "OTOA"="lightyellow","MYO15A"="#63666A","EYA1"="#FFE900","SIX1"="#2AD2C9","WFS1"="#9063CD","ACTG1"="#FCAEBB",
      "COL4A3"="#74D1Ea","LMX1A"="#EA6811","TECTA"="#DDA46F","COL11A1"="#69B3E7","COL1A1"="#E0457B","KCNQ4"="#7C3A2D",
      "MT-TL1"="black","MYO6"="#E4002B","POU4F3"="#0057B8","TMPRSS3"="#ABD156","COCH"="#EE55FF","COL4A5"="#951123",
      "MITF"="#BBEBBE","MYO7A"="burlywood","OTOGL"="#FFE760","POU3F4"="#0055E7","TMC1"="#9994DD","ANKRD11"="#666666",
      "CASP10"="#FEFAAA","CLCNKA"="#aa11aa","CRB1"="#ab7111","DSPP"="#E2325F","ELMOD3"="#EA5411","ESRRB"="#ffee66",
      "EYA4"="green","GATA3"="blue","GSDME"="darkblue","ILDR1"="purple","KIT"="pink","LOXHD1"="skyblue",
      "MT-RNR1"="grey","MYH9"="white","NLRP3"="#a651be","OTOG"="#1a5a54","PDE6B"="#e4a65a","SLC12A3"="#b1c654",
      "SOX10"="#A8eeee","SSBP1"="#f16aaa","TJP2"="sienna4","SLC12A2"="#abcee5",
      "COQ6"="wheat3", "CTCF"="palegoldenrod", "DNAJC3"="khaki4", "MT-ND5"="honeydew3", "RBM10"="azure", "SMAD4"="lightsalmon3", 
      "SPATA5"="seagreen2")

oncoPrint(oncogene,
          alter_fun = list(
            background = function(x, y, w, h) {
              grid.rect(x,y,w*0.9,h*0.9,
                        gp = gpar(fill = "white", col = "grey"))
            },
            "GJB2" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["GJB2"], col = NA)),
            "SLC26A4" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["SLC26A4"], col = NA)),
            "STRC" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["STRC"], col = NA)),
            "USH2A" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["USH2A"], col = NA)),
            "CDH23" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["CDH23"], col = NA)),
            "MPZL2" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["MPZL2"], col = NA)),
            "OTOA" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["OTOA"], col = NA)),
            "MYO15A" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["MYO15A"], col = NA)),
            "EYA1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["EYA1"], col = NA)),
            "SIX1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["SIX1"], col = NA)),
            "WFS1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["WFS1"], col = NA)),
            "ACTG1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["ACTG1"], col = NA)),
            "COL4A3" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["COL4A3"], col = NA)),
            "LMX1A" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["LMX1A"], col = NA)),
            "TECTA" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["TECTA"], col = NA)),
            "COL11A1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["COL11A1"], col = NA)),
            "COL1A1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["COL1A1"], col = NA)),
            "KCNQ4" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["KCNQ4"], col = NA)),
            "MT-TL1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["MT-TL1"], col = NA)),
            "MYO6" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["MYO6"], col = NA)),
            "POU4F3" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["POU4F3"], col = NA)),
            "TMPRSS3" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["TMPRSS3"], col = NA)),
            "COCH" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["COCH"], col = NA)),
            "COL4A5" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["COL4A5"], col = NA)),
            "MITF" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["MITF"], col = NA)),
            "MYO7A" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["MYO7A"], col = NA)),
            "OTOGL" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["OTOGL"], col = NA)),
            "POU3F4" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["POU3F4"], col = NA)),
            "TMC1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["TMC1"], col = NA)),
            "ANKRD11" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["ANKRD11"], col = NA)),
            "CASP10" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["CASP10"], col = NA)),
            "CLCNKA" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["CLCNKA"], col = NA)),
            "CRB1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["CRB1"], col = NA)),
            "DSPP" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["DSPP"], col = NA)),
            "ELMOD3" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["ELMOD3"], col = NA)),
            "ESRRB" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["ESRRB"], col = NA)),
            "EYA4" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["EYA4"], col = NA)),
            "GATA3" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["GATA3"], col = NA)),
            "GSDME" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["GSDME"], col = NA)),
            "ILDR1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["ILDR1"], col = NA)),
            "KIT" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["KIT"], col = NA)),
            "LOXHD1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["LOXHD1"], col = NA)),
            "MT-RNR1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["MT-RNR1"], col = NA)),
            "MYH9" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["MYH9"], col = NA)),
            "NLRP3" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["NLRP3"], col = NA)),
            "OTOG" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["OTOG"], col = NA)),
            "PDE6B" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["PDE6B"], col = NA)),
            "SLC12A3" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["SLC12A3"], col = NA)),
            "SOX10" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["SOX10"], col = NA)),
            "SSBP1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["SSBP1"], col = NA)),
            "SLC12A2" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["SLC12A2"], col = NA)),
            "TJP2" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["TJP2"], col = NA)),
            "COQ6" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["COQ6"], col = NA)),
            "CTCF" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["CTCF"], col = NA)),
            "DNAJC3" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["DNAJC3"], col = NA)),
            "MT-ND5" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["MT-ND5"], col = NA)),
            "RBM10" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["RBM10"], col = NA)),
            "SMAD4" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["SMAD4"], col = NA)),
            "SPATA5" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["SPATA5"], col = NA))
          ),
          col = col,row_names_side = "left",pct_side = "right",show_column_names = FALSE,top_annotation = NULL,right_annotation = NULL,
          alter_fun_is_vectorized = FALSE,show_pct=FALSE)

lgd=Legend(labels=c("GJB2","SLC26A4","STRC","USH2A","CDH23","MPZL2","OTOA","MYO15A","EYA1","SIX1","WFS1","ACTG1",
                    "COL4A3","LMX1A","TECTA","COL11A1","COL1A1","KCNQ4","MT-TL1","MYO6","POU4F3","TMPRSS3","COCH","COL4A5",
                    "MITF","MYO7A","OTOGL","POU3F4","TMC1","ANKRD11","CASP10","CLCNKA","CRB1","DSPP","ELMOD3","ESRRB",
                    "EYA4","GATA3","GSDME","ILDR1","KIT","LOXHD1","MT-RNR1","MYH9","NLRP3","OTOG","PDE6B","SLC12A3",
                    "SOX10","SSBP1","TJP2","SLC12A2","COQ6", "CTCF", "DNAJC3", "MT-ND5", "RBM10", "SMAD4","SPATA5"),
           legend_gp=gpar(fill=c("lightgoldenrod2","lightcyan3","lavenderblush2","seagreen4","violetred3","thistle4",
                                 "lightyellow","#63666A","#FFE900","#2AD2C9","#9063CD","#FCAEBB",
                                 "#74D1Ea","#EA6811","#DDA46F","#69B3E7","#E0457B","#7C3A2D",
                                 "black","#E4002B","#0057B8","#ABD156","#EE55FF","#951123",
                                 "#BBEBBE","burlywood","#FFE760","#0055E7","#9994DD","#666666",
                                 "#FEFAAA","#aa11aa","#ab7111","#E2325F","#EA5411","#ffee66",
                                 "green","blue","darkblue","purple","pink","skyblue",
                                 "grey","white","#a651be","#1a5a54","#e4a65a","#b1c654",
                                 "#A8eeee","#f16aaa","sienna4","#abcee5","wheat3","palegoldenrod",
                                 "khaki4","honeydew3","azure","lightsalmon3","seagreen2")),nrow=3)
dev.new()
draw(lgd)

#####Fig2 Bar&Pie#####
#예전A
diagno <- as.data.frame(data$진단유전자)
diagno1 <- as.data.frame(original_data %>% filter(진단유전자 %in% gene_list))
diagno2 <- as.data.frame(original_data %>% filter(진단유전자 %in% except))
a <- c("Genetic",nrow(diagno1))
b <- c("Nongenetic",nrow(diagno2))
diagno <- as.data.frame(rbind(a,b))
names(diagno) <- c("Etiology","count")
diagno[1,2] <- 206

plot_ly(diagno,labels=~`Etiology`,values=~count,textposition='none',
        marker=list(colors=c('#FFB3B5','#86D8A4')),sort=FALSE) %>% 
  add_pie(hole=0.5) %>% 
  layout(xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         showlegend=FALSE)

lgd1=Legend(labels=c("Diagnosed","Non Diagnosed"),
            legend_gp=gpar(fill=c('#FFB3B5','#86D8A4')))
dev.new()
draw(lgd1)

#A
#1
names <- familycount$진단유전자

ggplot(data=familycount,mapping=aes(x=reorder(진단유전자,`v1`),y=v1))+
  scale_x_discrete(limits=rev(names))+
  scale_y_continuous(limits=c(0,25),expand=c(0,0))+
  geom_bar(mapping=aes(fill=진단유전자),stat="identity",fill="#42E1a5")+
  coord_flip()+
  geom_label(aes(label=v1),nudge_y=0.6,nudge_x=-0.5,size=4,data=numberdata)+
  theme_classic()+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank())

#2
data_fill <- variant_type[,c(1,3,4,5)]

data_vr1 <- data_fill %>% 
  group_by(진단유전자,variant_type_nt) %>% 
  summarize(v1=sum(count)) 

data_vr2 <- data_fill%>% 
  group_by(진단유전자,variant_type_aa) %>% 
  summarize(v2=sum(count))

names(data_vr1)[2] <- c("varianttype1")
names(data_vr2)[2] <- c("varianttype2")

data_vr1$진단유전자 <- factor(data_vr1$진단유전자, levels = names)
data_vr1<- data_vr1[order(data_vr1$진단유전자, decreasing = FALSE), ]
data_vr1$진단유전자 <- as.character(data_vr1$진단유전자)
data_vr1$varianttype1 <- as.integer(data_vr1$varianttype1)
data_vr1$v1 <- as.integer(data_vr1$v1)

data_vr2$진단유전자 <- factor(data_vr2$진단유전자, levels = names)
data_vr2<- data_vr2[order(data_vr2$진단유전자, decreasing = FALSE), ]
data_vr2$진단유전자 <- as.character(data_vr2$진단유전자)

#proportion bar plot
ggplot(data_vr1,aes(x=진단유전자,y=v1,fill=factor(varianttype1)))+
  scale_x_discrete(limits=rev(names))+
  scale_y_continuous(limits=c(0,43),expand=c(0,0))+
  geom_bar(stat="identity")+
  scale_fill_discrete_qualitative(palette = 'Set 3')+
  theme_classic()+
  theme(legend.position=c(0.75,0.2),legend.title=element_blank(),
        axis.title.y=element_blank(),axis.title.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+  
  coord_flip()


ggplot(data_vr2,aes(x=진단유전자,y=v2,fill=factor(varianttype2)))+
  scale_x_discrete(limits=rev(names))+
  scale_y_continuous(limits=c(0,46),expand=c(0,0))+
  geom_bar(stat="identity")+
  scale_fill_discrete_qualitative(palette = 'Set 3')+
  theme_classic()+
  theme(legend.position=c(0.65,0.2),legend.title=element_blank(),
        axis.title.y=element_blank(),axis.title.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+  
  coord_flip()

lgd1=Legend(labels=c("Missense","Frameshift","Nonsense","Inframe del/dup","Splicing","SV","Mitochondria"),
            legend_gp=gpar(fill=c('#FFB3B5','#D5C783','#86D8A4','#6CD6E4','#AAC8FC','#D1ACF9','#F5B5F0')))
dev.new()
draw(lgd1)



#B
variant_type[is.na(variant_type)] <- ""

#sum(count) or n() 으로 바꿔주면 됨/ n()은 MF3, sum은 Supp
data_vr1 <- variant_type %>% 
  group_by(`variant_type_nt`) %>% 
  summarize(v1=n())
names(data_vr1)[1] <- c("varianttype1")

plot_ly(data_vr1,labels=~varianttype1,values=~v1,textposition='none',
        marker=list(colors=c('#FFB3B5','#D5C783','#86D8A4','#6CD6E4')),sort=FALSE) %>% 
  add_pie(hole=0.5) %>% 
  layout(xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         showlegend=FALSE)

lgd1=Legend(labels=c("SNVs","Short Indel","SVs","Mitochondria"),
            legend_gp=gpar(fill=c("#FFB3B5","#D5C783","#86D8A4","#6CD6E4")),nrow=1)
dev.new()
draw(lgd1)

data_vr2 <- variant_type %>% 
  group_by(`variant_type_aa`) %>% 
  summarize(v1=n())
names(data_vr2)[1] <- c("varianttype2")


plot_ly(data_vr2,labels=~varianttype2,values=~v1,textposition='none',
        marker=list(colors=c('#FFB3B5','#D5C783','#86D8A4','#6CD6E4','#AAC8FC','#D1ACF9','#154689')),sort=FALSE) %>% 
  add_pie(hole=0.5) %>% 
  layout(xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         showlegend=FALSE)

lgd1=Legend(labels=c("Missense","Frameshift","Nonsense","Inframe del/dup","Splicing","SV","Mitochondria"),
            legend_gp=gpar(fill=c('#FFB3B5','#D5C783','#86D8A4','#6CD6E4','#AAC8FC','#D1ACF9','#154689')),nrow=2)
dev.new()
draw(lgd1)

#width=450 이후 자르기
#SNVcount(최종본에 사용 X)
names(SNVcount) <- c('Mutation','Count')

plot_ly(SNVcount,labels=~Mutation,values=~Count,textposition='none',
        marker=list(colors=c('#FFB3B5','#D5C783','#86D8A4','#6CD6E4','#AAC8FC','#D1ACF9',"#154689","#BDE121","#FFA111","#FAF1F1","#459786","#7C3A2D")),sort=FALSE) %>% 
  add_pie(hole=0.5)%>%
  layout(xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         showlegend=FALSE)

lgd1=Legend(labels=c("A>C","A>G","A>T","C>A","C>G","C>T","G>A","G>C","G>T","T>A","T>C","T>G"),
            legend_gp=gpar(fill=c('#FFB3B5','#D5C783','#86D8A4','#6CD6E4','#AAC8FC','#D1ACF9',"#154689","#BDE121","#FFA111","#FAF1F1","#459786","#7C3A2D")),nrow=2)
dev.new()
draw(lgd1)

##inheritance(C)
data_in <- data %>% 
  group_by(`inheritance pattern`) %>% 
  summarize(v1=n())

plot_ly(data_in,labels=~`inheritance pattern`,values=~v1,textposition='none',
        marker=list(colors=c('#FFB3B5','#D5C783','#86D8A4','#6CD6E4','#AAC8FC')),sort=FALSE) %>% 
  add_pie(hole=0.5) %>% 
  layout(xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         showlegend=FALSE)

lgd1=Legend(labels=c("Autosomal Recessive","Autosomal Dominant","De Novo ","X-Linked","Non-Mendelian"),
            legend_gp=gpar(fill=c('#FFB3B5','#D5C783','#86D8A4','#6CD6E4','#AAC8FC')),nrow=2)
dev.new()
draw(lgd1)

######progressive(잠시)########
a <- oncoprintdata2[,c(1,13)]

for(i in 1:nrow(a)){
  if(a$progressive[i]==1){
    a$Sustan[i] = "Sustan"
  }else{
    a$Sustan[i] = ""
  }
}

for(i in 1:nrow(a)){
  if(a$progressive[i]==2){
    a$Mild[i] = "Mild"
  }else{
    a$Mild[i] = ""
  }
}

for(i in 1:nrow(a)){
  if(a$progressive[i]==3){
    a$Non[i] = "Non"
  }else{
    a$Non[i] = ""
  }
}

for(i in 1:nrow(a)){
  if(a$progressive[i]==1|a$progressive[i]==2){
    a$Sus_Mild[i] = "Sus_Mild"
  }else{
    a$Sus_Mild[i] = ""
  }
}

a2 <- a[,c(1,3,4,5,6,5)]
a2 <- as.data.frame(t(a2))
names(a2) <- a2[1,]
a2 <- a2[-1,]

col_ord <- colnames(a2)
row_ord <- rownames(a2)

col=c("Sustan"="#AAC8FC","Mild"="#86D8A4","Non"="lightyellow",
      "Sus_Mild"="#63666A")


oncoPrint(a2,alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x,y,w*0.9,h*0.9,
              gp = gpar(fill = "white", col = "grey"))
  },
  "Sustan" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                            gp = gpar(fill = col["Sustan"], col = NA)),
  "Mild" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                          gp = gpar(fill = col["Mild"], col = NA)),
  "Non" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                         gp = gpar(fill = col["Non"], col = NA)),
  "Sus_Mild" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                              gp = gpar(fill = col["Sus_Mild"], col = NA))),
  col = col,row_names_side = "left",pct_side = "right",show_column_names = FALSE,top_annotation = NULL,right_annotation = NULL,
  alter_fun_is_vectorized = FALSE,column_order = col_ord,row_order=row_ord,show_pct=FALSE)

test <- a2

for(i in 1:nrow(test)){
  for(j in 1:ncol(test)){
    if(test[i,j]!=""){
      test[i,j] <- data4[1,j]
    }
    else{
      next
    }
  }
}


col_ord <- colnames(a2)
row_ord <- rownames(a2)

col=c("Missense" = "red", "Frameshift" = "orange","Nonsense"="yellow","Inframe_del/dup"="green","Splicing"="blue","SV"="purple","mtRNA"="black",
      "Missense*" = "red", "Frameshift*" = "orange","Nonsense*"="yellow","Inframe_del/dup*"="green","Splicing*"="blue","SV*"="purple","mtRNA*"="black")

oncoPrint(a2,
          alter_fun = list(
            background = function(x, y, w, h) {
              grid.rect(x,y,w*0.9,h*0.9,
                        gp = gpar(fill = "white", col = "grey"))
            },
            "Missense" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",
                                                        gp = gpar(fill = col["Missense"], col = NA)),
            "Frameshift" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",
                                                          gp = gpar(fill = col["Frameshift"], col = NA)),
            "Nonsense" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",
                                                        gp = gpar(fill = col["Nonsense"], col = NA)),
            "Inframe_del/dup" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",
                                                               gp = gpar(fill = col["Inframe_del/dup"], col = NA)),
            "SV" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",
                                                  gp = gpar(fill = col["SV"], col = NA)),
            "Splicing" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",
                                                        gp = gpar(fill = col["Splicing"], col = NA)),
            "mtRNA" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",
                                                     gp = gpar(fill = col["mtRNA"], col = NA)),
            "Missense*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="left",
                                                         gp = gpar(fill = col["Missense*"], col = NA)),
            "Frameshift*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="left",
                                                           gp = gpar(fill = col["Frameshift*"], col = NA)),
            "Nonsense*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="left",
                                                         gp = gpar(fill = col["Nonsense*"], col = NA)),
            "Inframe_del/dup*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="left",
                                                                gp = gpar(fill = col["Inframe_del/dup*"], col = NA)),
            "SV*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="left",
                                                   gp = gpar(fill = col["SV*"], col = NA)),
            "Splicing*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="left",
                                                         gp = gpar(fill = col["Splicing*"], col = NA)),
            "mtRNA*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",
                                                      gp = gpar(fill = col["mtRNA"], col = NA)),
            "0" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                 gp = gpar(fill = col["0"], col = NA)),
            "5" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                 gp = gpar(fill = col["5"], col = NA)),
            "20" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                  gp = gpar(fill = col["20"], col = NA)),
            "30" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                  gp = gpar(fill = col["30"], col = NA)),
            "40" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                  gp = gpar(fill = col["40"], col = NA)),
            "50" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                  gp = gpar(fill = col["50"], col = NA)),
            "60" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                  gp = gpar(fill = col["60"], col = NA))
          ), col = col,row_names_side = "left",pct_side = "right",show_column_names = FALSE,top_annotation = NULL,right_annotation = NULL,
          alter_fun_is_vectorized = FALSE,column_order = col_ord,row_order=row_ord,show_pct=FALSE)

#############
row_ord <- rownames(test2)
col_ord <- colnames(test2)

oncoPrint(test2,
          alter_fun = list(
            background = function(x, y, w, h) {
              grid.rect(x,y,w*0.9,h*0.9,
                        gp = gpar(fill = "white", col = "grey"))
            },
            "1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["1"], col = NA)),
            "3" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["3"], col = NA)),
            "4" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["4"], col = NA)),
            "5" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["5"], col = NA)),
            "6" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["6"], col = NA)),
            "7" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["7"], col = NA))
          ),
          col = col,row_names_side = "left",pct_side = "right",show_column_names = FALSE,top_annotation = NULL,right_annotation = NULL,
          alter_fun_is_vectorized = FALSE,column_order = col_ord,row_order = row_ord,show_pct=FALSE)
#######denovo data########
denovo <- data[data$`inheritance pattern`==3,]
denovo <- denovo[,c(3,4,5)]
denovo <- na.omit(denovo)
rownames(denovo)=NULL
a <- c("c.2175_2178delCAAA;p.Asn725LysfsTer23","2","2")
denovo <- rbind(denovo,a)

names(denovo) <- c("variant","variant_type_nt","variant_type_aa")

denovo_vr1 <- denovo %>% 
  group_by(`variant_type_nt`) %>% 
  summarize(v1=n())
names(denovo_vr1)[1] <- c("varianttype1")

plot_ly(denovo_vr1,labels=~varianttype1,values=~v1,textposition='none',
        marker=list(colors=c('#FFB3B5','#D5C783','#86D8A4','#6CD6E4')),sort=FALSE) %>% 
  add_pie(hole=0.5) %>% 
  layout(xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         showlegend=FALSE)

lgd1=Legend(labels=c("SNVs","Short Indel","SVs"),
            legend_gp=gpar(fill=c("#FFB3B5","#D5C783","#86D8A4")),nrow=1)
dev.new()
draw(lgd1)


####
denovo_vr2 <- denovo %>% 
  group_by(`variant_type_aa`) %>% 
  summarize(v1=n())
names(denovo_vr2)[1] <- c("varianttype2")


plot_ly(denovo_vr2,labels=~varianttype2,values=~v1,textposition='none',
        marker=list(colors=c('#FFB3B5','#D5C783','#86D8A4','#6CD6E4','#AAC8FC','#D1ACF9','#154689')),sort=FALSE) %>% 
  add_pie(hole=0.5) %>% 
  layout(xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         showlegend=FALSE)

lgd1=Legend(labels=c("Missense","Frameshift","Nonsense","Inframe del/dup","Splicing","SV"),
            legend_gp=gpar(fill=c('#FFB3B5','#D5C783','#86D8A4','#6CD6E4','#AAC8FC','#D1ACF9')),nrow=2)
dev.new()
draw(lgd1)

denovo_mutation_count <- cbind(sapply(str_split(denovo$variant,";"),"[",1),1)
denovo_mutation_count <- as.data.frame(denovo_mutation_count)
denovo_mutation_count$V1 <- str_sub(denovo_mutation_count$V1,-3,-1)
denovo_mutation_count$V2 <- as.numeric(denovo_mutation_count$V2)
SNVcount <- denovo_mutation_count %>% 
  group_by(`V1`) %>% 
  summarize(v1=sum(`V2`))
SNVcount
SNVcount<- SNVcount[c(1,3:6),]
rm(denovo_mutation_count)

names(SNVcount) <- c('Mutation','Count')
plot_ly(SNVcount,labels=~Mutation,values=~Count,textposition='none',
        marker=list(colors=c('#FFB3B5','#D5C783','#86D8A4','#6CD6E4','#AAC8FC')),sort=FALSE) %>% 
  add_pie(hole=0.5)%>%
  layout(xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         showlegend=FALSE)

lgd1=Legend(labels=c("A>G","C>A","C>T","G>A","T>C"),
            legend_gp=gpar(fill=c('#FFB3B5','#D5C783','#86D8A4','#6CD6E4','#AAC8FC')),nrow=1)
dev.new()
draw(lgd1)

#####onco_full######
oncoprintdata4 <- oncoprintdata2
oncoprintdata4 <- oncoprintdata4[,-16]

data3 <- oncoprintdata4[,c(2,1,16:29)]
data3 <- as.data.frame(t(data3))
names(data3) <- c(1:ncol(data3))
data3 <- data3[-1,]

data4 <- oncoprintdata4[,2:3]
data4 <- as.data.frame(t(data4))
names(data4) <- data4[1,]
data4 <- data4[-1,]
rownames(data4) <- "Variant_Type"

#variant type으로 각각을 표시해주기
for(i in 2:13){
  for(j in 1:ncol(data3)){
    if(data3[i,j]!=""){
      data3[i,j] <- data4[1,j]
    }
    else{
      next
    }
  }
}

data3 <- data3[-1,]

col_ord <- colnames(data3)
row_ord <- rownames(data3)

col=c("Missense" = "red", "Frameshift" = "orange","Nonsense"="yellow","Inframe_del/dup"="green","Splicing"="blue","SV"="purple","mtRNA"="black",
      "Missense*" = "red", "Frameshift*" = "orange","Nonsense*"="yellow","Inframe_del/dup*"="green","Splicing*"="blue","SV*"="purple","mtRNA*"="black",
      "GJB2"="lightgoldenrod2","SLC26A4"="lightcyan3","STRC"="lavenderblush2","USH2A"="seagreen4","CDH23"="violetred3","MPZL2"="thistle4",
      "OTOA"="lightyellow","MYO15A"="#63666A","EYA1"="#FFE900","SIX1"="#2AD2C9","WFS1"="#9063CD","ACTG1"="#FCAEBB",
      "COL4A3"="#74D1Ea","LMX1A"="#EA6811","TECTA"="#DDA46F","COL11A1"="#69B3E7","COL1A1"="#E0457B","KCNQ4"="#7C3A2D",
      "MT-TL1"="black","MYO6"="#E4002B","POU4F3"="#0057B8","TMPRSS3"="#ABD156","Sustan"="sienna4","Mild"="#A8dddd","None"="#A8bbbb")

oncoPrint(data3,
          alter_fun = list(
            background = function(x, y, w, h) {
              grid.rect(x,y,w*0.9,h*0.9,
                        gp = gpar(fill = "white", col = "grey"))
            },
            "Missense" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",gp = gpar(fill = col["Missense"], col = NA)),
            "Frameshift" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",gp = gpar(fill = col["Frameshift"], col = NA)),
            "Nonsense" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",gp = gpar(fill = col["Nonsense"], col = NA)),
            "Inframe_del/dup" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",gp = gpar(fill = col["Inframe_del/dup"], col = NA)),
            "SV" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",gp = gpar(fill = col["SV"], col = NA)),
            "Splicing" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",gp = gpar(fill = col["Splicing"], col = NA)),
            "mtRNA" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",gp = gpar(fill = col["mtRNA"], col = NA)),
            "Missense*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="left",gp = gpar(fill = col["Missense*"], col = NA)),
            "Frameshift*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="left",gp = gpar(fill = col["Frameshift*"], col = NA)),
            "Nonsense*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="left",gp = gpar(fill = col["Nonsense*"], col = NA)),
            "Inframe_del/dup*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="left",gp = gpar(fill = col["Inframe_del/dup*"], col = NA)),
            "SV*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="left",gp = gpar(fill = col["SV*"], col = NA)),
            "Splicing*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="left",gp = gpar(fill = col["Splicing*"], col = NA)),
            "mtRNA*" = function(x, y, w, h) grid.rect(x, y, w*0.45, h*0.9,just="right",gp = gpar(fill = col["mtRNA"], col = NA)),
            "GJB2" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["GJB2"], col = NA)),
            "SLC26A4" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["SLC26A4"], col = NA)),
            "STRC" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["STRC"], col = NA)),
            "USH2A" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["USH2A"], col = NA)),
            "CDH23" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["CDH23"], col = NA)),
            "MPZL2" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["MPZL2"], col = NA)),
            "OTOA" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["OTOA"], col = NA)),
            "MYO15A" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["MYO15A"], col = NA)),
            "EYA1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["EYA1"], col = NA)),
            "SIX1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["SIX1"], col = NA)),
            "WFS1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["WFS1"], col = NA)),
            "ACTG1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["ACTG1"], col = NA)),
            "COL4A3" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["COL4A3"], col = NA)),
            "LMX1A" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["LMX1A"], col = NA)),
            "TECTA" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["TECTA"], col = NA)),
            "COL11A1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["COL11A1"], col = NA)),
            "COL1A1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["COL1A1"], col = NA)),
            "KCNQ4" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["KCNQ4"], col = NA)),
            "MT-TL1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["MT-TL1"], col = NA)),
            "MYO6" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["MYO6"], col = NA)),
            "POU4F3" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["POU4F3"], col = NA)),
            "TMPRSS3" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["TMPRSS3"], col = NA)),
            "Sustan" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["Sustan"], col = NA)),
            "Mild" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["Mild"], col = NA)),
            "None" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9,gp = gpar(fill = col["None"], col = NA))
          ), col = col,row_names_side = "left",pct_side = "right",show_column_names = FALSE,top_annotation = NULL,right_annotation = NULL,
          alter_fun_is_vectorized = FALSE,column_order = col_ord,row_order=row_ord,show_pct=FALSE)

lgd=Legend(labels=c("Missense", "Frameshift","Nonsense","Inframe_del/dup","Splicing","SV","mtRNA",
                    "Sustan","Mild","None"),
           legend_gp=gpar(fill=c("red","orange","yellow","green","blue","purple","black",
                                 "sienna4","#A8dddd","#A8bbbb")),nrow=1)
dev.new()
draw(lgd)
#####잡다 카운트######
step[step$step1==1 & step$step2==0 & step$step3==0 & step$step4==0,] %>%   
  group_by(`진단유전자`) %>% 
  summarize(v1=n())

x <- step[step$step1==1 & step$step2==1 & step$step3==0 & step$step4==0,] %>% 
  group_by(`진단유전자`) %>%
  summarize(v1=n())

step[step$step1==0 & step$step2==1 & step$step3==0 & step$step4==0,] %>% 
  group_by(`진단유전자`) %>%
  summarize(v1=n())

step[step$step1==1 & step$step2==1 & step$step3==1 & step$step4==0,] %>% 
  group_by(`진단유전자`) %>%
  summarize(v1=n())

step[step$step4==1,]%>%
  group_by(`진단유전자`) %>%
  summarize(v1=n())

count(step[step$step4==1,])


