#cd /ddn/gs1/home/mackja/UKBB/
#/ddn/gs1/home/mackja/miniconda3/envs/R411/bin/R

#/ddn/gs1/home/mackja/miniconda3/envs/R/bin/R

library(dplyr) 
library(data.table)
setDTthreads(40)

#library("devtools")
library(SeqArray)
library(SNPRelate)
library(GWASTools)
library(GENESIS)
library(stringr)
library(Matrix)

IDs<-fread("/ddn/gs1/shared/ukbiobank/data_versions/data_v1/dataformat_txt/eid.txt")

#Adding in ancestral data from pan-UKBB
Anc<-fread("/ddn/gs1/shared/ukbiobank/housejs_requests/ancestral_assignments_mackja/Files for retman/all_pops_non_eur_pruned_within_pop_pc_covs.tsv")
Anc<-Anc[,c("s", "pop")]

#For bridging from returned dataset
Bridge<-fread("/ddn/gs1/shared/ukbiobank/housejs_requests/ancestral_assignments_mackja/ukb57849bridge31063.txt")
Anc1<-merge(Bridge, Anc, by.x = "V2", by.y = "s", all.y = T)

IDs_Anc<-merge(IDs, Anc1, by.x = "eid", by.y = "V1", all.x = T) 

#called UKBB1
load("/ddn/gs1/home/mackja/UKBB/UKBB_Phenotype.RData")

##Run "Create_UKBB_PhenotypeFile.R" to get full set for code 


#### Load ICD9/10 code selection (modified from Excel table)
#### Note: We will rely on phecode mapping to define controls - normal delivery
#### Controls: 650	females	Normal delivery	pregnancy complications	O800,O801,O808,O809,650,6509	650

input_case9 <- fread("/ddn/gs1/home/mackja/UKBB/Selected_ICD9_Case_UKB.txt")
#input_excl9 <- fread("/ddn/gs1/home/mackja/UKBB/Selected_ICD9_Excl_UKB.txt")
input_control9 <- fread("/ddn/gs1/home/mackja/UKBB/Selected_ICD9_Control_UKB.txt")
input_excl9_control<-fread("/ddn/gs1/home/mackja/UKBB/Selected_ICD9_Excl_Control_UKB.txt")

input_case10 <- fread("/ddn/gs1/home/mackja/UKBB/Selected_ICD10_Case_UKB.txt")
#input_excl10 <- fread("/ddn/gs1/home/mackja/UKBB/Selected_ICD10_Excl_UKB.txt")
input_control10 <- fread("/ddn/gs1/home/mackja/UKBB/Selected_ICD10_Control_UKB.txt")
input_excl10_control<-fread("/ddn/gs1/home/mackja/UKBB/Selected_ICD10_Excl_Control_UKB.txt")

input_self_report<-fread("/ddn/gs1/home/mackja/UKBB/Selected_Report_Case_UKB.txt")
input_self_report_excl<-fread("/ddn/gs1/home/mackja/UKBB/Selected_Report_Excl_UKB.txt")


#  limit to selected codes
selectedICD9_case <- as.character(input_case9[selectable == "Y",`coding`])
#selectedICD9_excl <- as.character(input_excl9[selectable == "Y",`coding`])
selectedICD9_control <- as.character(input_control9[selectable == "Y",`coding`])
selectedICD9_excl_control <- as.character(input_excl9_control[selectable == "Y",`coding`])

selectedICD10_case <- input_case10[selectable == "Y",`coding`]
#selectedICD10_excl <- input_excl10[selectable == "Y",`coding`]
selectedICD10_control <- input_control10[selectable == "Y",`coding`]
selectedICD10_excl_control <- input_excl10_control[selectable == "Y",`coding`]

selected_self_report <- input_self_report[selectable == "Y",`coding`]
selected_self_report_excl <- input_self_report_excl[selectable == "Y", `coding`]

#grepl code

self_report_cols<-cbind(IDs,Self_report)

colcheck<-function(x) {UKBB[grepl(paste("^",selected_self_report,collapse="|",sep=""),x),eid]}
self_report_case<-sapply(self_report_cols,colcheck)
self_report_case_full<-unique(unlist(self_report_case, recursive = TRUE, use.names = FALSE))
length(self_report_case_full)
#[1] 2162 - GH/PE


colcheck<-function(x) {UKBB[grepl(paste("^",selected_self_report_excl,collapse="|",sep=""),x),eid]}
self_report_excl<-sapply(self_report_cols,colcheck)
self_report_excl_full<-unique(unlist(self_report_excl, recursive = TRUE, use.names = FALSE))
length(self_report_excl_full)
#[1] 2547 = GH/PE and Gestational Diabetes

#ICD 10

#Some testing
#selectedICD10_case<-c("O600","O601", "O602", "O603")
#selectedICD10_case<-c("O13") for only GH cases

icd10_cols<-cbind(IDs,Main_ICD10,Sec_ICD10,Inpat_ICD10,Death1_ICD10,Death2_ICD10)

#selectedICD10_case<-c("O13","O16")

selectedICD10_case<-c("O244")
colcheck<-function(x) {UKBB[grepl(paste("^",selectedICD10_case,collapse="|",sep=""),x),eid]}
icd10_case<-sapply(icd10_cols,colcheck)
icd10_case_full<-unique(unlist(icd10_case, recursive = TRUE, use.names = FALSE))
length(icd10_case_full)
#[1] 1544 - HDP

#Just GH-602 - 
#1279 - GH + Unspecified hypertension (Gestational Hypertension) - O13 + O16

##Define Preeclampsia as subtype
#icd10_cols<-cbind(IDs,Main_ICD10,Sec_ICD10,Inpat_ICD10,Death1_ICD10,Death2_ICD10)

selectedICD10_case<-c("O14")
colcheck<-function(x) {UKBB[grepl(paste("^",selectedICD10_case,collapse="|",sep=""),x),eid]}
icd10_case_pe<-sapply(icd10_cols,colcheck)
icd10_case_pe<-unique(unlist(icd10_case_pe, recursive = TRUE, use.names = FALSE))
length(icd10_case_pe)
#[1] 428


#colcheck<-function(x) {UKBB[grepl(paste("^",selectedICD10_excl,collapse="|",sep=""),x),eid]}
#icd10_excl<-sapply(icd10_cols,colcheck)
#icd10_excl_full<-unique(unlist(icd10_excl, recursive = TRUE, use.names = FALSE))
#length(icd10_excl_full)
#[1] 23

#Missed Z33 and Z34 previously for defining pregnancy

#selectedICD10_control<-c("Z349")
colcheck<-function(x) {UKBB[grepl(paste("^",selectedICD10_control,collapse="|",sep=""),x),eid]}
icd10_control<-sapply(icd10_cols,colcheck)
icd10_control_full<-unique(unlist(icd10_control, recursive = TRUE, use.names = FALSE))
length(icd10_control_full)
#[1] 15278
#[2] REDUCED TO NORMAL DELIVERY: 2643
#[3] Normal delivery + other delivery methods: 3547

#selectedICD10_excl_control<-c("Z373")
colcheck<-function(x) {UKBB[grepl(paste("^",selectedICD10_excl_control,collapse="|",sep=""),x),eid]}
icd10_excl_control<-sapply(icd10_cols,colcheck)
icd10_excl_control_full<-unique(unlist(icd10_excl_control, recursive = TRUE, use.names = FALSE))
length(icd10_excl_control_full)
#[1] 6969


#ICD 9

icd9_cols<-cbind(IDs,Main_ICD9,Sec_ICD9,Inpat_ICD9)

#selectedICD9_case<-c("6427")

selectedICD9_case<-c("6488")
colcheck<-function(x) {UKBB[grepl(paste("^",selectedICD9_case,collapse="|",sep=""),x),eid]}
icd9_case<-sapply(icd9_cols,colcheck)
icd9_case_full<-unique(unlist(icd9_case, recursive = TRUE, use.names = FALSE))
length(icd9_case_full)
#[1] 1 for PE


##7/5/2022: Adding a preeclampsia only and gestational hypertension only case group

#colcheck<-function(x) {UKBB[grepl(paste("^",selectedICD9_excl,collapse="|",sep=""),x),eid]}
#icd9_excl<-sapply(icd9_cols,colcheck)
#icd9_excl_full<-unique(unlist(icd9_excl, recursive = TRUE, use.names = FALSE))
#length(icd9_excl_full)
#[1] 0

colcheck<-function(x) {UKBB[grepl(paste("^",selectedICD9_control,collapse="|",sep=""),x),eid]}
icd9_control<-sapply(icd9_cols,colcheck)
icd9_control_full<-unique(unlist(icd9_control, recursive = TRUE, use.names = FALSE))
length(icd9_control_full)
#[1] 0

colcheck<-function(x) {UKBB[grepl(paste("^",selectedICD9_excl_control,collapse="|",sep=""),x),eid]}
icd9_excl_control<-sapply(icd9_cols,colcheck)
icd9_excl_control_full<-unique(unlist(icd9_excl_control, recursive = TRUE, use.names = FALSE))
length(icd9_excl_control_full)
#[1] 266

#Recall the files created above: ID list for each data type and group based on criteria provided in text files
#self_report_case_full,icd10_case_full,icd10_excl_full,icd10_control_full,icd10_excl_control_full
#icd9_case_full,icd9_excl_full,icd9_control_full,icd9_excl_control_full


ehrdata_extracted_self_report_case <- UKBB[eid %in% c(self_report_case_full),]
ehrdata_extracted_self_report_excl <- UKBB[eid %in% c(self_report_excl_full),]
ehrdata_extracted_case_pe <- UKBB[eid %in% c(icd10_case_pe),]
ehrdata_extracted_case <- UKBB[eid %in% c(icd9_case_full,icd10_case_full)] 
#ehrdata_extracted_excl <- UKBB[eid %in% c(icd10_excl_full,icd9_excl_full),]
ehrdata_extracted_control <- UKBB[eid %in% c(icd10_control_full,icd9_control_full),]
ehrdata_extracted_excl_control <- UKBB[eid %in% c(icd10_excl_control_full,icd9_excl_control_full),]


#### Define five subsets
self_report_cases <- unique(ehrdata_extracted_self_report_case$eid)  #n = 2162 [GDM:406 ]
self_report_excl <- unique(ehrdata_extracted_self_report_excl$eid) #n=2547
cases <- unique(ehrdata_extracted_case$eid)  #n = 1279 [GDM:259], All HDP: 1544
cases_pe <- unique(ehrdata_extracted_case_pe$eid) #n = 428
#excl <- unique(ehrdata_extracted_excl$eid)    #n = 23
controls <- unique(ehrdata_extracted_control$eid) #n = 15278, NOW 3547, F-2643
excl_controls <- unique(ehrdata_extracted_excl_control$eid) #n = 7740 (before: 7180)


#### define final controls by excluding individuals that are cases or listed in the excl or in self report:
controls <- controls[which((!controls %in% cases) & (!controls %in% excl_controls) & (!controls %in% self_report_cases) & (!controls %in% self_report_excl))] #8732

ehrdata_extracted_self_report_case$group_self<- 1
self_report_cases_full<-ehrdata_extracted_self_report_case[!duplicated(ehrdata_extracted_self_report_case$eid), ]

self_report_excl_full<-ehrdata_extracted_self_report_excl[!duplicated(ehrdata_extracted_self_report_excl$eid), ]

ehrdata_extracted_case$group_icd<- 1
cases_full <- ehrdata_extracted_case[!duplicated(ehrdata_extracted_case$eid), ] 

ehrdata_extracted_case_pe$group_icd<- 1
cases_pe <- ehrdata_extracted_case_pe[!duplicated(ehrdata_extracted_case_pe$eid), ] 

#ehrdata_extracted_excl$group_icd <- -1
#excl_full <- ehrdata_extracted_excl[!duplicated(ehrdata_extracted_excl$eid), ]

ehrdata_extracted_control$group_icd<- 0
ehrdata_extracted_control$group_self<- 0
controls_full <- ehrdata_extracted_control[!duplicated(ehrdata_extracted_control$eid), ]

ehrdata_extracted_excl_control$group_icd <- -1
excl_control_full <- ehrdata_extracted_excl_control[!duplicated(ehrdata_extracted_excl_control$eid), ]

controls_full0 <- controls_full[!(controls_full$eid %in% self_report_excl_full$eid),] #14914
controls_full1 <- controls_full0[!(controls_full0$eid %in% excl_control_full$eid),] 
#8732
controls_full2 <- controls_full1[!(controls_full1$eid %in% self_report_cases_full$eid),] #n = 8214

#nrow(controls_full2)

#cases_full1 <- cases_full[!(cases_full$eid %in% excl_full$eid),]  #n = 1604
#self_report_cases_full1 <- self_report_cases_full[!(self_report_cases_full$eid %in% excl_full$eid),] #n = 2161

#SUMMARY
##EXCLUSION CRITERIA FOR CASES: ANY MULTIPLE PREGNANCY
##EXCLUSION CRITERIA FOR CONTROLS: ANY MULTIPLE PREGNANCY AND/OR PREGNANCY-RELATED HYPERTENSION


#THUS, THERE ARE 4,343 PARTICIPANTS WHO REPORTED HYPERTENSION IN PREGNANCY

#Merging process!
Cases_Controls<- as.data.frame(rbind(cases_full, controls_full2,fill=TRUE))

Cases_Controls<- as.data.frame(rbind(cases_pe, controls_full2,fill=TRUE))

Cases_Controls<- as.data.frame(rbind(self_report_cases_full, controls_full2,fill=TRUE))

#FINAL QC
#Cases_Controls<-Cases_Controls[,-c(749)]
sum(duplicated(Cases_Controls$eid))
nrow(Cases_Controls) #3919 [GDM:2638]
Cases_Controls<-dplyr::filter(Cases_Controls, Cases_Controls$`22001-0.0`== 0 )
nrow(Cases_Controls) #3795,XX participants marked as male removed [GDM:2566]


#Group variable: Cases = 1, Controls = 0
#N = 3800

table(Cases_Controls$group_icd)

#HDP
#0    1
#2315 1480


#GH

#0    1
#2315 1225

#GH-ONLY

#0    1
#2315 582

#PE
#0    1
#2315  410


#GDM
#   0    1
#2315 216

table(Cases_Controls$group_self)
#HDP
#0    1
#2315 2068



#Let's include ethnic categories
#Cases/Controls by Race/Ethnicity: https://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=1001
table(Cases_Controls$'21000', Cases_Controls$group_icd)

Cases_Controls$Ethnic_Cat[Cases_Controls$'21000' =="1" | Cases_Controls$'21000' =="1001" | Cases_Controls$'21000' =="1002" | Cases_Controls$'21000' =="1003"] = "White"
Cases_Controls$Ethnic_Cat[Cases_Controls$'21000' =="2" | Cases_Controls$'21000' =="2001" | Cases_Controls$'21000' =="2002" | Cases_Controls$'21000' =="2003"| Cases_Controls$'21000' =="2004" | Cases_Controls$'21000' =="6" | Cases_Controls$'21000' =="3" | Cases_Controls$'21000' =="3001" | Cases_Controls$'21000' =="3002" | Cases_Controls$'21000' =="3003" | Cases_Controls$'21000' =="3004" | Cases_Controls$'21000' =="5" | is.na(Cases_Controls$'21000')| Cases_Controls$'21000' =="-3" | Cases_Controls$'21000' =="-1"  ] = "Mixed/Other"
Cases_Controls$Ethnic_Cat[Cases_Controls$'21000' =="4" | Cases_Controls$'21000' =="4001" | Cases_Controls$'21000' =="4002" | Cases_Controls$'21000' =="4003"] = "Black"



#SUBSET FILE TO KEEP RELEVANT COLUMNS

#"34-0.0" (YOB)     "31-0.0"(SEX)      "22000-0.0"(Batch) "21001-0.0"(BMI)
#"group_icd"   "group_self"  "Ethnic_Cat"



#Cases_Controls1<-Cases_Controls[,c("eid","34-0.0","31-0.0", "22000-0.0", "21001-0.0", "group_icd", "group_self","Ethnic_Cat")]

Cases_Controls1<-Cases_Controls[,c("eid","21022-0.0","22001-0.0", "22000-0.0", "21001-0.0", "group_icd", "group_self","Ethnic_Cat")]

##Merge in ancestry information
Cases_Controls2<-merge(Cases_Controls1, IDs_Anc, by="eid", all.x = T)

Cases_Controls2 <- Cases_Controls2[!duplicated(Cases_Controls2$eid), ]

Cases_Controls2<-Cases_Controls2[Cases_Controls2$eid %in% UKBB1$eid,]

table(Cases_Controls2$group_icd)

#HDP
#0    1
#2311 1480

#GH Only
#0    1
#2311  582

#GH + unspec
#0    1
#2311 1225

#PE only
#2311  410

#GDM
#0    1
#2311  216



#SAVE PHENO FILE

fwrite(Cases_Controls2,"/ddn/gs1/home/mackja/UKBB/APOs/GestHyp_Preeclampsia/HDP.txt",sep="\t")


#########

##Let's get the IDs to filter the PLINK fileset

CC1<-fread("/ddn/gs1/home/mackja/UKBB/APOs/GestHyp_Preeclampsia/HDP.txt")

CC1_ids<-CC1[,c("eid","eid")]

names(CC1_ids)<-NULL

fwrite(CC1_ids, "/ddn/gs1/home/mackja/UKBB/APOs/GestHyp_Preeclampsia/HDP_IDs.txt",sep="\t")
          
#################################################################################
#Let's calculate PCs for UKBB study

#DOUBLY FIRST: FILTER KINSHIP
##Threshold of 0.0884

#Ensure that the IDs match between pheno file and geno file before PC calc, remove related individuals

IDs<-fread("/ddn/gs1/home/mackja/UKBB/APOs/GestHyp_Preeclampsia/HDP_IDs.txt")
rel<-fread("/ddn/gs1/home/mackja/UKBB/KING_GRM/UKBB_KINGunrelated_toberemoved.txt")
IDs<-IDs[!IDs$V1 %in% rel$V1][[1]]


#snpgdsBED2GDS(bed.fn = "/data/mackja/UKBB/all_chromosomes/ukb_cal_all_v2.bed", 
#              bim.fn = "/data/mackja/UKBB/all_chromosomes/ukb_cal_all_v2.bim", 
#              fam.fn = "/data/mackja/UKBB/all_chromosomes/ukb_cal_all_v2.fam", 
#              out.gdsfn = "UKBB_full.gds")

#Reading in GDS data + LD pruning


gds_ukbb <- snpgdsOpen("/ddn/gs1/home/mackja/UKBB/KING_GRM/UKBB_full.gds")

#SNPs are kept if MAF>5%, nonmonomorphic, call rate (95%) , limiting to autosomal SNPs, LD calculated by correlation coefficient
#Sliding window of bp: 1,000,000

#LD PRUNING TO CONSIDER FOR NON-EURO POPS: On the contrary, if you're dealing with a population with massive haplotype diversity (such as sub-Saharan Africans) and you would assume that your study population is not large enough to be a true representative of all haplotypes in the population you may want to use a more relaxed r^2 threshold for pruning.

#snpset <- snpgdsLDpruning(gds, method="corr", slide.max.bp=5e6, 
#                          ld.threshold=sqrt(0.1), verbose=FALSE, maf=.05,missing.rate=.05)

#pruned <- unlist(snpset, use.names=FALSE)
#length(pruned)  #110845

#Checkpoint: Saved the snpset file--can load in here:
#save(snpset, file="LDPruned_SNPs.RData")

load("/ddn/gs1/home/mackja/UKBB/KING_GRM/LDPruned_SNPs.RData")
pruned <- unlist(snpset, use.names=FALSE)
#length(pruned) #110845


#Pairwise Measures of Ancestry Divergence

KINGmat<-snpgdsIBDKING(gds_ukbb,sample.id=IDs,num.thread=100)
snpgdsClose(gds_ukbb)

KINGmat1<-kingToMatrix(KINGmat)

geno_ukbb <- GdsGenotypeReader(filename = "/ddn/gs1/home/mackja/UKBB/KING_GRM/UKBB_full.gds")
genoData <- GenotypeData(geno_ukbb)

pcair <- pcair(genoData,kinobj = KINGmat1,divobj=KINGmat1, snp.include = pruned,sample.include=IDs,num.cores = 100)



#####MERGING PCs WITH PHENO FILE####
PCs1<-as.data.frame(pcair$vectors)

PCs1 <- cbind(rownames(PCs1),PCs1)
rownames(PCs1) <- NULL
colnames(PCs1) <- c(names(PCs1)) #to not write all the column names

colnames(PCs1)[1] <- "eid" 
names(PCs1)

PC_11<-PCs1[,1:11]
PC_11$eid<-as.integer(PC_11$eid)

#Bringing in the Phenotype file

CC1<-fread("/ddn/gs1/home/mackja/UKBB/APOs/GestHyp_Preeclampsia/HDP.txt")
CC1<-CC1[,-c("V2")]
CC1$eid.1<-CC1$eid

PC_pheno_final1<-merge(CC1, PC_11, by = "eid")

PC_pheno_final1<-rename(PC_pheno_final1, `#FID`=eid, IID=eid.1, Age=`21022-0.0`, Batch=`22000-0.0`, BMI=`21001-0.0`)

PC_pheno_final1<-PC_pheno_final1[,c("#FID","IID","group_icd","group_self","Age","BMI","Ethnic_Cat","Batch","pop","22001-0.0","V1","V2","V3","V4","V5","V6","V7","V8","V9","V10")]

PC_pheno_final1$pop[PC_pheno_final1$pop==""] = "NONE"

fwrite(PC_pheno_final1, "/ddn/gs1/home/mackja/UKBB/APOs/GDM/GDM_PCs.txt", sep="\t", quote = F , na=NA)

showfile.gds(closeall=TRUE)
######################################################################################################################################
#"34-0.0" (YOB)     "31-0.0"(SEX)      "22000-0.0"(Batch) "21001-0.0"(BMI)
#"group_icd"   "group_self"  "Ethnic_Cat"


library(data.table)
library(tableone)

#Create tableone
pheno<-fread("/ddn/gs1/home/mackja/UKBB/APOs/GestHyp_Preeclampsia/HDP_PCs.txt")
IDs<-fread("/ddn/gs1/home/mackja/UKBB/GWAS/PLINK_HDP.group_icd.glm.firth.id")

pheno<-pheno[pheno$IID %in% IDs$IID]

Vars<-c("Age","BMI","Ethnic_Cat")

Cat_Var<-c("Ethnic_Cat")

Cont_Var<-c("Age","BMI")

print(CreateTableOne
      (vars=Vars, factorVars = Cat_Var,
        strata="group_icd",data=pheno, includeNA= T, testExact = fisher.test, addOverall = T)) 


print(CreateContTable(vars=Cont_Var,strata="group_icd",data=pheno, addOverall = T
),nonnormal=Cont_Var)






