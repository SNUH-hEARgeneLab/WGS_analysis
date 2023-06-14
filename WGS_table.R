install_github("dustinfife/fifer")
install.packages("moonBook")
library(readxl)
library(dplyr)
library(tidyr)
library(moonBook)
library(rcompanion)
library(devtools)
library(fifer)

tabledata <- read.csv("C:\\Users\\SNUH\\Desktop\\WGS_project\\data.csv",header = TRUE,fileEncoding = "euc-kr")
tabledata <- tabledata[!is.na(tabledata$진단유전자),]
tabledata <- tabledata[tabledata$진단유전자!="",]
tabledata[is.na(tabledata)] <- ""
rownames(tabledata)=NULL

#Genetic=1, NonGenetic=0
a <- tabledata$진단유전자 %in% except
for(i in 1:nrow(tabledata)){
  if(a[i]==TRUE){
    tabledata$Diag[i] = "0"
  }
  else{
    tabledata$Diag[i] = "1"
  }
}

#나이수정
tabledata$나이<- gsub("M","000",tabledata$나이)
tabledata$나이<- as.numeric(tabledata$나이)

for(i in 1:nrow(tabledata)){
  if(tabledata$나이[i]>48000){
    tabledata$나이[i] <-  4
  }else if(tabledata$나이[i]>36000){
    tabledata$나이[i] <-  3
  }else if(tabledata$나이[i]>24000){
    tabledata$나이[i] <-  2
  }else if(tabledata$나이[i]>12000){
    tabledata$나이[i] <-  1
  }else if(tabledata$나이[i]>=1000){
    tabledata$나이[i] <-  0
  }
}

supptabledata <- tabledata[,c(29:32,34:41,45,22,44,16,24,43)]
tabledata <- tabledata[,c(29:32,34:41,45,13,44)]
####Main table 1####
#1.Gender / Age / Lossonset / Syndromic (data=test1)
for(i in 1:nrow(tabledata)){
  if(tabledata$inheritance.pattern[i]==2|tabledata$inheritance.pattern[i]==4|tabledata$inheritance.pattern[i]==5){
    tabledata$fam[i]="1"
  }else{
    tabledata$fam[i]="0"
  }
}


test1 <- tabledata[,c(1:4,12,13,15)]
test1 <- as.data.frame(test1)

a <- mytable(Diag~.,data=test1,method=3)
a
summary(a)

median(test1[test1$Diag==1,]$나이)
max(test1[test1$Diag==1,]$나이)
min(test1[test1$Diag==1,]$나이)
median(test1[test1$Diag==0,]$나이)
max(test1[test1$Diag==0,]$나이)
min(test1[test1$Diag==0,]$나이)


#Hearing Loss Type~Configuration,Asym
#한쪽 귀라도 Mixed인 경우 Mixed로 Count

for(i in 1:nrow(tabledata)){
  if(tabledata$hearing.loss.type..Rt[i]==2|tabledata$hearing.loss.type..Lt[i]==2){
    tabledata$Mixed[i] = "1"
  }else{
    tabledata$Mixed[i] = "0"
  }
}

a <- mytable(Diag~Mixed,data=tabledata,method=3)
a

#Mixed인 사람들은 제외하고 진행(data=test2)
test2 <- tabledata[tabledata$Mixed!=1,]

#Asymmetry
a <- mytable(Diag~asymmetry,data=test2,method=3)
a

#Hearing Loss Severity 
test3 <- test2
test2 <- test2[test2$asymmetry==3,]

for(i in 1:nrow(test2)){
  if(test2$hearing.severity..Rt.[i]==1&test2$hearing.severity..Lt.[i]==1){
    test2$MildMod[i] = "1"
    test2$Sever[i]="1"
  }else if(test2$hearing.severity..Rt.[i]==2&test2$hearing.severity..Lt.[i]==2){
    test2$MildMod[i] = "1"
    test2$Sever[i]="1"
  }else if(test2$hearing.severity..Rt.[i]==1&test2$hearing.severity..Lt.[i]==2){
    test2$MildMod[i] = "1"
    test2$Sever[i]="1"
  }else if(test2$hearing.severity..Rt.[i]==2&test2$hearing.severity..Lt.[i]==1){
    test2$MildMod[i] = "1"
    test2$Sever[i]="1"
  }else{
    test2$MildMod[i] = "0"
    test2$Sever[i]="0"
  }
}

for(i in 1:nrow(test2)){
  if(test2$hearing.severity..Rt.[i]==3&test2$hearing.severity..Lt.[i]==3){
    test2$ModSev[i] = "1"
    test2$Sever[i]="2"
  }else{
    test2$ModSev[i] = "0"
  }
}

for(i in 1:nrow(test2)){
  if(test2$hearing.severity..Rt.[i]==4&test2$hearing.severity..Lt.[i]==4){
    test2$SevPro[i] = "1"
    test2$Sever[i]="3"
  }else if(test2$hearing.severity..Rt.[i]==5&test2$hearing.severity..Lt.[i]==5){
    test2$SevPro[i] = "1"
    test2$Sever[i]="3"
  }else if(test2$hearing.severity..Rt.[i]==4&test2$hearing.severity..Lt.[i]==5){
    test2$SevPro[i] = "1"
    test2$Sever[i]="3"
  }else if(test2$hearing.severity..Rt.[i]==5&test2$hearing.severity..Lt.[i]==4){
    test2$SevPro[i] = "1"
    test2$Sever[i]="3"
  }else{
    test2$SevPro[i] = "0"
  }
}

a <- mytable(Diag~MildMod+ModSev+SevPro,data=test2,method=3)
a <- mytable(Diag~Sever,data=test2,method=3)
a

#Configuration
for(i in 1:nrow(test2)){
  if(test2$hearing.configuration..Lt.[i]==1&test2$hearing.configuration..Rt.[i]==1){
    test2$Flat[i] = "1"
    test2$Config[i] = "1"
  }else{
    test2$Flat[i] = "0"
    test2$Config[i]="0"
  }
}

for(i in 1:nrow(test2)){
  if(test2$hearing.configuration..Lt.[i]==2&test2$hearing.configuration..Rt.[i]==2){
    test2$DownSlop[i] = "1"
    test2$Config[i] = "2"
  }else if(test2$hearing.configuration..Lt.[i]==3&test2$hearing.configuration..Rt.[i]==3){
    test2$DownSlop[i] = "1"
    test2$Config[i] = "2"
  }else if(test2$hearing.configuration..Lt.[i]==2&test2$hearing.configuration..Rt.[i]==3){
    test2$DownSlop[i] = "1"
    test2$Config[i] = "2"
  }else if(test2$hearing.configuration..Lt.[i]==3&test2$hearing.configuration..Rt.[i]==2){
    test2$DownSlop[i] = "1"
    test2$Config[i] = "2"
  }else{
    test2$DownSlop[i] = "0"
  }
}

for(i in 1:nrow(test2)){
  if(test2$hearing.configuration..Lt.[i]==4&test2$hearing.configuration..Rt.[i]==4){
    test2$Cookie[i] = "1"
    test2$Config[i] = "3"
  }else{
    test2$Cookie[i] = "0"
  }
}

for(i in 1:nrow(test2)){
  if(test2$hearing.configuration..Lt.[i]==5&test2$hearing.configuration..Rt.[i]==5){
    test2$UpSlop[i] = "1"
    test2$Config[i] = "4"
  }else{
    test2$UpSlop[i] = "0"
  }
}

a <- mytable(Diag~Flat+DownSlop+Cookie+UpSlop,data=test2,method=3)
a <- mytable(Diag~Config,data=test2,method=3)

a

#####Adjusted OR#####
for(i in 1:nrow(test1)){
  if(test1$hearing.loss.onset[i] == "1"){
    test1$EI[i] = "1"
  }
  else{
    test1$EI[i] = "0"
  }
}

for(i in 1:nrow(test1)){
  if(test1$hearing.loss.onset[i] == "4"){
    test1$AO[i] = "1"
  }
  else{
    test1$AO[i] = "0"
  }
}

for(i in 1:nrow(test1)){
  if(test1$syndromic_1[i] == "2"|test1$syndromic_2[i]=="2"){
    test1$SYN[i] = "1"
  }
  else{
    test1$SYN[i] = "0"
  }
}

for(i in 1:nrow(test1)){
  if(test1$syndromic_1[i] == "1"){
    test1$SYN2[i] = "0"
  }
  else{
    test1$SYN2[i] = "1"
  }
}

test1$Diag <- as.numeric(test1$Diag)
test1 <- cbind(test1,tabledata$Mixed)

model <- glm(Diag~EI+AO+Family.history+SYN2,family='binomial'(link='logit'),data=test1)
summary(model)
exp(cbind(OR=coef(model), confint(model)))

#Asym~Config
test2$Diag <- as.numeric(test2$Diag)
test3$Diag <- as.numeric(test3$Diag)

for(i in 1:nrow(test3)){
  if(test3$asymmetry[i] == "1"|test3$asymmetry[i]=="2"){
    test3$ASYM[i] = "1"
  }
  else{
    test3$ASYM[i] = "0"
  }
}

model <- glm(Diag~ASYM,family='binomial'(link='logit'),data=test3)
summary(model)
exp(cbind(OR=coef(model), confint(model)))

model <- glm(Diag~UpSlop,family='binomial'(link='logit'),data=test2)
summary(model)
exp(cbind(OR=coef(model), confint(model)))

#####Supp table 2#####
supptabledata <- supptabledata[supptabledata$WGS!="",]
rownames(supptabledata)=NULL

#1.Gender / Age / Loss onset / Syndromic (data=test1)
test1 <- supptabledata[,c(1:4,12,13,15:17)]
test1 <- as.data.frame(test1)

for(i in 1:nrow(test1)){
  if(test1$WGS.type[i]=="1"|test1$WGS.type[i]=="2"){
    test1$WGS.type[i]="1"
  }else{
    test1$WGS.type[i]="3"
  }
}

a <- mytable(Diag~.,data=test1,method=3)
a
summary(a)

test1.1 <- test1

for(i in 1:nrow(test1.1)){
  if(test1.1$hearing.loss.onset[i]==1){
    test1.1$EI[i]="1"
    test1.1$LI[i]="0"
    test1.1$LD[i]="0"
    test1.1$AO[i]="0"
  }else if(test1.1$hearing.loss.onset[i]==2){
    test1.1$EI[i]="0"
    test1.1$LI[i]="1"
    test1.1$LD[i]="0"
    test1.1$AO[i]="0"
  }else if(test1.1$hearing.loss.onset[i]==3){
    test1.1$EI[i]="0"
    test1.1$LI[i]="0"
    test1.1$LD[i]="1"
    test1.1$AO[i]="0"
  }else{
    test1.1$EI[i]="0"
    test1.1$LI[i]="0"
    test1.1$LD[i]="0"
    test1.1$AO[i]="1"
  }
}

for(i in 1:nrow(test1.1)){
  if(test1.1$carrier[i]==0){
    test1.1$DR[i]="1"
    test1.1$CD[i]="0"
    test1.1$CND[i]="0"
  }else if(test1.1$carrier[i]==1){
    test1.1$DR[i]="0"
    test1.1$CD[i]="1"
    test1.1$CND[i]="0"
  }else{
    test1.1$DR[i]="0"
    test1.1$CD[i]="0"
    test1.1$CND[i]="1"
  }
}

for(i in 1:nrow(test1.1)){
  if(test1.1$WGS.type[i]==1){
    test1.1$SD[i]="1"
    test1.1$TRIO[i]="0"
  }else{
    test1.1$SD[i]="0"
    test1.1$TRIO[i]="1"
  }
}

for(i in 1:nrow(test1.1)){
  if(test1.1$syndromic_1[i]==2|test1.1$syndromic_2[i]==2){
    test1.1$SYN[i]="1"
  }else{
    test1.1$SYN[i]="0"
  }
}

a <- mytable(Diag~.,data=test1.1,method=3)
a
summary(a)

test1.1$Diag <- as.numeric(test1.1$Diag)
model <- glm(Diag~EI+syndromic_1+DR+TRIO,family='binomial'(link='logit'),data=test1.1)
summary(model)
exp(cbind(OR=coef(model), confint(model)))

#Hearing Loss Type~Configuration,Asym
#한쪽 귀라도 Mixed인 경우 Mixed로 Count

for(i in 1:nrow(supptabledata)){
  if(supptabledata$hearing.loss.type..Rt[i]==2|supptabledata$hearing.loss.type..Lt[i]==2){
    supptabledata$Mixed[i] = "1"
  }else{
    supptabledata$Mixed[i] = "0"
  }
}

a <- mytable(Diag~Mixed,data=supptabledata,method=3)
a
summary(a)

#Mixed인 사람들은 제외하고 진행(data=test2)
test2 <- supptabledata[supptabledata$Mixed!=1,]

#Asymmetry
a <- mytable(Diag~asymmetry,data=test2,method=3)
a

#Hearing Loss Severity 
for(i in 1:nrow(test2)){
  if(test2$hearing.severity..Rt.[i]==1&test2$hearing.severity..Lt.[i]==1){
    test2$MildMod[i] = "1"
    test2$severity[i]="1"
  }else if(test2$hearing.severity..Rt.[i]==2&test2$hearing.severity..Lt.[i]==2){
    test2$MildMod[i] = "1"
    test2$severity[i]="1"
  }else if(test2$hearing.severity..Rt.[i]==1&test2$hearing.severity..Lt.[i]==2){
    test2$MildMod[i] = "1"
    test2$severity[i]="1"
  }else if(test2$hearing.severity..Rt.[i]==2&test2$hearing.severity..Lt.[i]==1){
    test2$MildMod[i] = "1"
    test2$severity[i]="1"
  }else{
    test2$MildMod[i] = "0"
    test2$severity[i]="0"
  }
}

for(i in 1:nrow(test2)){
  if(test2$hearing.severity..Rt.[i]==3&test2$hearing.severity..Lt.[i]==3){
    test2$ModSev[i] = "1"
    test2$severity[i]="2"
  }else{
    test2$ModSev[i] = "0"
  }
}

for(i in 1:nrow(test2)){
  if(test2$hearing.severity..Rt.[i]==4&test2$hearing.severity..Lt.[i]==4){
    test2$SevPro[i] = "1"
    test2$severity[i]="3"
  }else if(test2$hearing.severity..Rt.[i]==5&test2$hearing.severity..Lt.[i]==5){
    test2$SevPro[i] = "1"
    test2$severity[i]="3"
  }else if(test2$hearing.severity..Rt.[i]==4&test2$hearing.severity..Lt.[i]==5){
    test2$SevPro[i] = "1"
    test2$severity[i]="3"
  }else if(test2$hearing.severity..Rt.[i]==5&test2$hearing.severity..Lt.[i]==4){
    test2$SevPro[i] = "1"
    test2$severity[i]="3"
  }else{
    test2$SevPro[i] = "0"
  }
}

test3 <- test2[test2$asymmetry==3,]
a <- mytable(Diag~MildMod+ModSev+SevPro,data=test3,method=3)
a
summary(a)

a <- mytable(Diag~severity,data=test3,method=3)
a
summary(a)


#Configuration

for(i in 1:nrow(test2)){
  if(test2$hearing.configuration..Lt.[i]==1&test2$hearing.configuration..Rt.[i]==1){
    test2$Flat[i] = "1"
    test2$config[i] = "1"
  }else{
    test2$Flat[i] = "0"
    test2$config[i] = "0"
  }
}

for(i in 1:nrow(test2)){
  if(test2$hearing.configuration..Lt.[i]==2&test2$hearing.configuration..Rt.[i]==2){
    test2$DownSlop[i] = "1"
    test2$config[i] = "2"
  }else if(test2$hearing.configuration..Lt.[i]==3&test2$hearing.configuration..Rt.[i]==3){
    test2$DownSlop[i] = "1"
    test2$config[i] = "2"
  }else if(test2$hearing.configuration..Lt.[i]==2&test2$hearing.configuration..Rt.[i]==3){
    test2$DownSlop[i] = "1"
    test2$config[i] = "2"
  }else if(test2$hearing.configuration..Lt.[i]==3&test2$hearing.configuration..Rt.[i]==2){
    test2$DownSlop[i] = "1"
    test2$config[i] = "2"
  }else{
    test2$DownSlop[i] = "0"
  }
}

for(i in 1:nrow(test2)){
  if(test2$hearing.configuration..Lt.[i]==4&test2$hearing.configuration..Rt.[i]==4){
    test2$Cookie[i] = "1"
    test2$config[i] = "3"
  }else{
    test2$Cookie[i] = "0"
  }
}

for(i in 1:nrow(test2)){
  if(test2$hearing.configuration..Lt.[i]==5&test2$hearing.configuration..Rt.[i]==5){
    test2$UpSlop[i] = "1"
    test2$config[i] = "4"
  }else{
    test2$UpSlop[i] = "0"
  }
}

test3 <- test2[test2$asymmetry==3,]
a <- mytable(Diag~Flat+DownSlop+Cookie+UpSlop,data=test3,method=3)
a
summary(a)

a <- mytable(Diag~config,data=test3,method=3)
a
summary(a)

#####Classification Statistics/double primary 포함#####
supptabledata2 <- read.csv("C:\\Users\\SNUH\\Desktop\\WGS_project\\testdata2.csv",header = TRUE,fileEncoding = "euc-kr")

supptabledata2 <- supptabledata2[!is.na(supptabledata2$진단유전자),]
supptabledata2 <- supptabledata2[supptabledata2$진단유전자!="",]
supptabledata2[is.na(supptabledata2)] <- ""
rownames(supptabledata2)=NULL

#Genetic=1, NonGenetic=0
a <- supptabledata2$진단유전자 %in% except
for(i in 1:nrow(supptabledata2)){
  if(a[i]==TRUE){
    supptabledata2$Diag[i] = "0"
  }
  else{
    supptabledata2$Diag[i] = "1"
  }
}

#나이수정
supptabledata2$나이<- gsub("M","000",supptabledata2$나이)
supptabledata2$나이<- as.numeric(supptabledata2$나이)

for(i in 1:nrow(supptabledata2)){
  if(supptabledata2$나이[i]>48000){
    supptabledata2$나이[i] <-  4
  }else if(supptabledata2$나이[i]>36000){
    supptabledata2$나이[i] <-  3
  }else if(supptabledata2$나이[i]>24000){
    supptabledata2$나이[i] <-  2
  }else if(supptabledata2$나이[i]>12000){
    supptabledata2$나이[i] <-  1
  }else if(supptabledata2$나이[i]>=1000){
    supptabledata2$나이[i] <-  0
  }
}

supptabledata2 <- supptabledata2[,c(29:32,34:41,45,22,44,16,24,42,43,13)]
supptabledata2 <- supptabledata2[supptabledata2$Diag==1,]
supptabledata2 <- supptabledata2[order(supptabledata2$functional.classification),]
rownames(supptabledata2) <- NULL

a <- supptabledata2[c(59:61,86:95),]
a$functional.classification <- c(2,2,2,4,4,4,7,7,7,7,7,7,7)
supptabledata2 <- rbind(supptabledata2,a)
rm(a)

for(i in 1:nrow(supptabledata2)){
  if(supptabledata2$functional.classification [i]=="3,4"|supptabledata2$functional.classification [i] == "3,7"){
    supptabledata2$functional.classification [i] =  "3"
  }else if(supptabledata2$functional.classification[i]=="1,2"){
    supptabledata2$functional.classification[i]="1"
  }
  else{
    next
  }
}
rownames(supptabledata2) <- NULL

for(i in 1:nrow(supptabledata2)){
  if(supptabledata2$inheritance.pattern[i]==2|supptabledata2$inheritance.pattern[i]==4|supptabledata2$inheritance.pattern[i]==5){
    supptabledata2$fam[i]="1"
  }else{
    supptabledata2$fam[i]="0"
  }
}

names(supptabledata2)[1] <- "Sex"
names(supptabledata2)[2] <- "Age_at_Test"


for(i in 1:nrow(supptabledata2)){
  if(supptabledata2$hearing.loss.onset[i]==1){
    supptabledata2$Early_ID[i] = "1"
    supptabledata2$Late_ID[i] = "0"
    supptabledata2$Adult_onset[i] = "0"
  }else if(supptabledata2$hearing.loss.onset[i]==2|supptabledata2$hearing.loss.onset[i]==3){
    supptabledata2$Early_ID[i] = "0"
    supptabledata2$Late_ID[i] = "1"
    supptabledata2$Adult_onset[i] = "0"
  }else if(supptabledata2$hearing.loss.onset[i]==4){
    supptabledata2$Early_ID[i] = "0"
    supptabledata2$Late_ID[i] = "0"
    supptabledata2$Adult_onset[i] = "1"
  }else{
    supptabledata2$Early_ID[i] = "0"
    supptabledata2$Late_ID[i] = "0"
    supptabledata2$Adult_onset[i] = "0"
  }
}

for(i in 1:nrow(supptabledata2)){
  if(supptabledata2$syndromic_1[i]==2|supptabledata2$syndromic_2[i]==2){
    supptabledata2$Syndromic[i] = "1"
  }else{
    supptabledata2$Syndromic[i] = "0"
  }
}

for(i in 1:nrow(supptabledata2)){
  if(supptabledata2$hearing.loss.type..Rt[i]==2|supptabledata2$hearing.loss.type..Lt[i]==2){
    supptabledata2$Mixed[i] = "1"
  }else{
    supptabledata2$Mixed[i] = "0"
  }
}

for(i in 1:nrow(supptabledata2)){
  if(supptabledata2$asymmetry[i]==1|supptabledata2$asymmetry[i]==2){
    supptabledata2$Asymmetric[i] = "1"
  }else{
    supptabledata2$Asymmetric[i] = "0"
  }
}

#Asymmetric과 Mixed는 제외하고 이후 진행
#Severity
for(i in 1:nrow(supptabledata2)){
  if(supptabledata2$Mixed[i]!="1" && supptabledata2$Asymmetric[i]!="1"){
    if(supptabledata2$hearing.severity..Lt.[i]==1|supptabledata2$hearing.severity..Lt.[i]==2|supptabledata2$hearing.severity..Rt.[i]==1|supptabledata2$hearing.severity..Rt.[i]==2){
      supptabledata2$MildMod[i] = "1"
      supptabledata2$ModSev[i]="0"
      supptabledata2$SevProf[i]="0"
    }else if(supptabledata2$hearing.severity..Lt.[i]==3|supptabledata2$hearing.severity..Rt.[i]==3){
      supptabledata2$MildMod[i]="0"
      supptabledata2$ModSev[i] = "1"
      supptabledata2$SevProf[i]="0"
    }else if(supptabledata2$hearing.severity..Lt.[i]==4|supptabledata2$hearing.severity..Lt.[i]==5|supptabledata2$hearing.severity..Rt.[i]==4|supptabledata2$hearing.severity..Rt.[i]==5){
      supptabledata2$MildMod[i]="0"
      supptabledata2$ModSev[i] = "0"
      supptabledata2$SevProf[i] = "1"
    }else{
      supptabledata2$MildMod[i] = "0"
      supptabledata2$ModSev[i] = "0"
      supptabledata2$SevProf[i] = "0"
    }
  }else{
    supptabledata2$MildMod[i] = "0"
    supptabledata2$ModSev[i] = "0"
    supptabledata2$SevProf[i] = "0"
  }
}

#Configuration
for(i in 1:nrow(supptabledata2)){
  if(supptabledata2$Mixed[i]!="1" && supptabledata2$Asymmetric[i]!="1"){
    if(supptabledata2$hearing.configuration..Lt.[i]==1|supptabledata2$hearing.configuration..Rt.[i]==1){
      supptabledata2$Flat[i] = "1"
      supptabledata2$DownSloping[i] = "0"
      supptabledata2$Cookie[i] = "0"
    }else if(supptabledata2$hearing.configuration..Lt.[i]==2|supptabledata2$hearing.configuration..Rt.[i]==2|supptabledata2$hearing.configuration..Lt.[i]==3|supptabledata2$hearing.configuration..Rt.[i]==3){
      supptabledata2$Flat[i] = "0"
      supptabledata2$DownSloping[i] = "1"
      supptabledata2$Cookie[i] = "0"
    }else if(supptabledata2$hearing.configuration..Lt.[i]==4|supptabledata2$hearing.configuration..Rt.[i]==4){
      supptabledata2$Flat[i] = "0"
      supptabledata2$DownSloping[i] = "0"
      supptabledata2$Cookie[i] = "1"
    }else{
      supptabledata2$Flat[i] = "0"
      supptabledata2$DownSloping[i] = "0"
      supptabledata2$Cookie[i] = "0"
    }
  }else{
    supptabledata2$Flat[i] = "0"
    supptabledata2$DownSloping[i]="0"
    supptabledata2$Cookie[i]="0"
  }
}

#Progressive
for(i in 1:nrow(supptabledata2)){
  if(supptabledata2$progressive[i]==1){
    supptabledata2$Sustan[i] = "1"
    supptabledata2$MildNone[i] ="0"
  }else if(supptabledata2$progressive[i]==2){
    supptabledata2$Sustan[i] = "0"
    supptabledata2$MildNone[i] ="1"
  }else if(supptabledata2$progressive[i]==3){
    supptabledata2$Sustan[i]="0"
    supptabledata2$MildNone[i]="1"
  }else{
    supptabledata2$Sustan[i]="0"
    supptabledata2$MildNone[i]="0"
  }
}

classifi <- supptabledata2[supptabledata2$functional.classification!=8,c(3,19,21:35)]
a <- mytable(functional.classification~.,data=classifi,method=3)
a
summary(a)


classifi2 <- read.csv("C:\\Users\\SNUH\\Desktop\\WGS_project\\2.csv",header = TRUE,fileEncoding = "euc-kr")
rownames(classifi2) <- classifi2[,1]
classifi2 <- t(classifi2[,-1])

chisq.post.hoc(classifi2[,11:12],simulate.p.value=TRUE)
chisq.test(classifi2[,11:12])



#####classification 조합용 bar graph#####
classifi3 <- read.csv("C:\\Users\\SNUH\\Desktop\\WGS_project\\2.csv",header = FALSE,fileEncoding = "euc-kr")
classifi3 <- as.data.frame(t(classifi3))
classifi3[1,1] <- "group"
names(classifi3) <- classifi3[1,]
classifi3 <- classifi3[-1,]

#syndromic(w/o mimics)
syndromic <- classifi3[,c(1,3,32)]
test <- syndromic %>% gather(a,n,Total,`synd(w/omimic)2`)
test$n <- as.numeric(test$n)

pal <- c('Total'='#86D8A4','synd(w/omimic)2'='#FCAEBB')
ggplot(test)+
  geom_bar(stat="identity",position="identity",mapping=aes(x=group,y=n,fill=a))+
  theme_classic()+
  scale_fill_manual(values=pal)+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank())

for(i in 1:7){
  test[i,3] <- test[i,3]-test[i+7,3]
}
ggplot(test)+
  geom_bar(stat="identity",position=position_fill(reverse=TRUE),mapping=aes(x=group,y=n,fill=a))+
  theme_classic()+
  scale_fill_manual(values=pal)+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank())

#syndromic(w mimics)
syndromic <- classifi3[,c(1,10,32)]
test <- syndromic %>% gather(a,n,Total,`Syndromic1`)
test$n <- as.numeric(test$n)

pal <- c('Total'='#86D8A4','Syndromic1'='#FCAEBB')
ggplot(test)+
  geom_bar(stat="identity",position="identity",mapping=aes(x=group,y=n,fill=a))+
  theme_classic()+
  scale_fill_manual(values=pal)+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank())

for(i in 1:7){
  test[i,3] <- test[i,3]-test[i+7,3]
}
test <- test[c(4,5,6,7,11,12,13,14),]
ggplot(test)+
  geom_bar(stat="identity",position=position_fill(reverse=TRUE),mapping=aes(x=group,y=n,fill=a))+
  theme_classic()+
  scale_fill_manual(values=pal)+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank())

#mixed
syndromic <- classifi3[,c(1,12,32)]
test <- syndromic %>% gather(a,n,Total,`Mixed1`)
test$n <- as.numeric(test$n)

pal <- c('Total'='#86D8A4','Mixed1'='#FCAEBB')
ggplot(test)+
  geom_bar(stat="identity",position="identity",mapping=aes(x=group,y=n,fill=a))+
  theme_classic()+
  scale_fill_manual(values=pal)+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank())

for(i in 1:7){
  test[i,3] <- test[i,3]-test[i+7,3]
}
test <- test[c(1,7,8,14),]
ggplot(test)+
  geom_bar(stat="identity",position=position_fill(reverse=TRUE),mapping=aes(x=group,y=n,fill=a))+
  theme_classic()+
  scale_fill_manual(values=pal)+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank())

#flat
syndromic <- classifi3[,c(1,22,32)]
test <- syndromic %>% gather(a,n,Total,`Flat1`)
test$n <- as.numeric(test$n)

pal <- c('Total'='#86D8A4','Flat1'='#FCAEBB')
ggplot(test)+
  geom_bar(stat="identity",position="identity",mapping=aes(x=group,y=n,fill=a))+
  theme_classic()+
  scale_fill_manual(values=pal)+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank())

for(i in 1:7){
  test[i,3] <- test[i,3]-test[i+7,3]
}

test <- test[c(1,3,4,8,10,11),]
ggplot(test)+
  geom_bar(stat="identity",position=position_fill(reverse=TRUE),mapping=aes(x=group,y=n,fill=a))+
  theme_classic()+
  scale_fill_manual(values=pal)+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank())

#downsloping
syndromic <- classifi3[,c(1,24,32)]
test <- syndromic %>% gather(a,n,Total,`DownSloping1`)
test$n <- as.numeric(test$n)

pal <- c('Total'='#86D8A4','DownSloping1'='#FCAEBB')
ggplot(test)+
  geom_bar(stat="identity",position="identity",mapping=aes(x=group,y=n,fill=a))+
  theme_classic()+
  scale_fill_manual(values=pal)+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank())

for(i in 1:7){
  test[i,3] <- test[i,3]-test[i+7,3]
}

test <- test[c(1,6,7,8,13,14),]
ggplot(test)+
  geom_bar(stat="identity",position=position_fill(reverse=TRUE),mapping=aes(x=group,y=n,fill=a))+
  theme_classic()+
  scale_fill_manual(values=pal)+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank())

#####supp table 1####
tabledata

rbind(table(tabledata$성별),round(prop.table(table(tabledata$성별))*100,2))
median(tabledata$나이)
min(tabledata$나이)
max(tabledata$나이)
rbind(table(tabledata$hearing.loss.onset),round(prop.table(table(tabledata$hearing.loss.onset))*100,2))
rbind(table(tabledata$syndromic_1),round(prop.table(table(tabledata$syndromic_1))*100,2))
rbind(table(tabledata$Family.history),round(prop.table(table(tabledata$Family.history))*100,2))
rbind(table(tabledata$Mixed),round(prop.table(table(tabledata$Mixed))*100,2))

tabledata2 <- tabledata[tabledata$Mixed!=1,]
rbind(table(tabledata2$asymmetry),round(prop.table(table(tabledata2$asymmetry))*100,2))

#asymmetry
for(i in 1:nrow(tabledata2)){
  if(tabledata2$asymmetry[i]==1|tabledata2$asymmetry[i]==2){
    tabledata2$Asymmetric[i] = "1"
  }else{
    tabledata2$Asymmetric[i] = "0"
  }
}

#Severity
for(i in 1:nrow(tabledata2)){
  if(tabledata2$Mixed[i]!="1" && tabledata2$Asymmetric[i]!="1"){
    if(tabledata2$hearing.severity..Lt.[i]==1|tabledata2$hearing.severity..Lt.[i]==2|tabledata2$hearing.severity..Rt.[i]==1|tabledata2$hearing.severity..Rt.[i]==2){
      tabledata2$MildMod[i] = "1"
      tabledata2$ModSev[i]="0"
      tabledata2$SevProf[i]="0"
      tabledata2$SevNorm[i] = "0"
    }else if(tabledata2$hearing.severity..Lt.[i]==3|tabledata2$hearing.severity..Rt.[i]==3){
      tabledata2$MildMod[i]="0"
      tabledata2$ModSev[i] = "1"
      tabledata2$SevProf[i]="0"
      tabledata2$SevNorm[i] = "0"
    }else if(tabledata2$hearing.severity..Lt.[i]==4|tabledata2$hearing.severity..Lt.[i]==5|tabledata2$hearing.severity..Rt.[i]==4|tabledata2$hearing.severity..Rt.[i]==5){
      tabledata2$MildMod[i]="0"
      tabledata2$ModSev[i] = "0"
      tabledata2$SevProf[i] = "1"
      tabledata2$SevNorm[i] = "0"
    }else if(tabledata2$hearing.severity..Lt.[i]==0|tabledata2$hearing.severity..Rt.[i]==0){
      tabledata2$MildMod[i]="0"
      tabledata2$ModSev[i] = "0"
      tabledata2$SevProf[i] = "0"     
      tabledata2$SevNorm[i] = "1"
    }else{
      tabledata2$MildMod[i] = "0"
      tabledata2$ModSev[i] = "0"
      tabledata2$SevProf[i] = "0"
      tabledata2$SevNorm[i] = "0"
    }
  }else{
    tabledata2$MildMod[i] = "0"
    tabledata2$ModSev[i] = "0"
    tabledata2$SevProf[i] = "0"
    tabledata2$SevNorm[i] = "0"
  }
}

#Configuration
for(i in 1:nrow(tabledata2)){
  if(tabledata2$Mixed[i]!="1" && tabledata2$Asymmetric[i]!="1"){
    if(tabledata2$hearing.configuration..Lt.[i]==1|tabledata2$hearing.configuration..Rt.[i]==1){
      tabledata2$Flat[i] = "1"
      tabledata2$DownSloping[i] = "0"
      tabledata2$Cookie[i] = "0"
      tabledata2$Confother[i] = "0"
    }else if(tabledata2$hearing.configuration..Lt.[i]==2|tabledata2$hearing.configuration..Rt.[i]==2|tabledata2$hearing.configuration..Lt.[i]==3|tabledata2$hearing.configuration..Rt.[i]==3){
      tabledata2$Flat[i] = "0"
      tabledata2$DownSloping[i] = "1"
      tabledata2$Cookie[i] = "0"
      tabledata2$Confother[i] = "0"
    }else if(tabledata2$hearing.configuration..Lt.[i]==4|tabledata2$hearing.configuration..Rt.[i]==4){
      tabledata2$Flat[i] = "0"
      tabledata2$DownSloping[i] = "0"
      tabledata2$Cookie[i] = "1"
      tabledata2$Confother[i] = "0"
    }else if(tabledata2$hearing.configuration..Lt.[i]==5|tabledata2$hearing.configuration..Rt.[i]==5|tabledata2$hearing.configuration..Lt.[i]==6|tabledata2$hearing.configuration..Rt.[i]==6){
      tabledata2$Flat[i] = "0"
      tabledata2$DownSloping[i] = "0"
      tabledata2$Cookie[i] = "0"    
      tabledata2$Confother[i] = "1"
    }else{
      tabledata2$Flat[i] = "0"
      tabledata2$DownSloping[i] = "0"
      tabledata2$Cookie[i] = "0"
      tabledata2$Confother[i] = "0"
    }
  }else{
    tabledata2$Flat[i] = "0"
    tabledata2$DownSloping[i]="0"
    tabledata2$Cookie[i]="0"
    tabledata2$Confother[i] = "0"
  }
}

rbind(table(tabledata2$MildMod),round(prop.table(table(tabledata2$MildMod))*100,2))
rbind(table(tabledata2$ModSev),round(prop.table(table(tabledata2$ModSev))*100,2))
rbind(table(tabledata2$SevProf),round(prop.table(table(tabledata2$SevProf))*100,2))
rbind(table(tabledata2$SevNorm),round(prop.table(table(tabledata2$SevNorm))*100,2))


rbind(table(tabledata2$Flat),round(prop.table(table(tabledata2$Flat))*100,2))
rbind(table(tabledata2$DownSloping),round(prop.table(table(tabledata2$DownSloping))*100,2))
rbind(table(tabledata2$Cookie),round(prop.table(table(tabledata2$Cookie))*100,2))
rbind(table(tabledata2$Confother),round(prop.table(table(tabledata2$Confother))*100,2))

step <- read.csv("C:\\Users\\SNUH\\Desktop\\WGS_project\\testdata.csv",header = TRUE,fileEncoding = "euc-kr")
step <- step[!is.na(step$진단유전자),]
step <- step[step$진단유전자!="",]
step[is.na(step)] <- ""
step <- step[,17:20]

rbind(table(step$step1),round(prop.table(table(step$step1))*100,2))
rbind(table(step$step2),round(prop.table(table(step$step2))*100,2))
rbind(table(step$step3),round(prop.table(table(step$step3))*100,2))
rbind(table(step$step4),round(prop.table(table(step$step4))*100,2))
##########
install.packages("epitools")
library(epitools)

data <- matrix(c(2,7,167,133),nrow=2,ncol=2,byrow=TRUE)
oddsratio(data)




