#  install.packages("devtools")
library("insect")   ##    devtools::install_github("shaunpwilkinson/insect") 
library(ape)
library('taxize') #  install.packages("taxize")
library(ShortRead) # packageVersion("ShortRead")  #  install.packages("ShortRead")
library(seqinr)  #  install.packages("seqinr")
library(dada2); packageVersion("dada2") # devtools::install_github("benjjneb/dada2") 
library(dplyr)
library(tidyr)

#################################################################################
#################################################################################
#################################################################################
#################################################################################
#  Mar 2021 done 
# arctic surface 
#  - rs surface euka and 18smini , Dammam , 18stoek, 18smini , euka02 02 and 03 - with ortega
#    with sarah - arctic surface , euka02, euka03   :  arctic core euka02 euka03 18smini 18s_stoeck
# with sarah - atacama 18smini euka02

###### !!!!!!!!!!!!!!!!!! set universal varibles   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
###    change everything with till next section
# enter study name for folder - examples  "extraction_test18" , "RS_surface", "Dammam", "Arctic_core" ,
##          "Atacama"    Arctic_surface
study_name<-"RS_surface"  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# enter primer name,  "co1" , "euka02" , "euka03"   , "rbclmini"  ,  "18smini"  ,  "vert" , "18s_stoeck" 
primer_name<-"18smini"

## if replicate- mulitiple same primer used in same study
## should be no if 1st replicate or no replicate
### if yes will add number to outputs
replicate<-"no"
#rep_label<-"rep23"

# set computer to get correct directory,  "win"  or "mac"
computer<-"mac"  

dir<-"C:/Users/nathangeraldi/Dropbox/"
## location to put filtered reads , one up from fastqfile below!!!!!!!!!!!!!!!!!!!!!
filtfile<-"eDNA_db/Arctic_surface/AS_2_euka02_CLL7R"

## location of miseq fastq files - within dir below
#fastqfile<-"eDNA_db/other_miseq/Extraction_test2_aug2018/primer_cut_vert_exttest2"
#fastqfile<-"eDNA_db/Nojood/noj_18suni_B5W44/noj_18suni_B5W44_cut"  #B5W44 APM90 
#
#fastqfile<-"eDNA_db/Red_sea_surface/RS_euka02_CLBR2/RS_euka02_CLBR2_cut"


#fastqfile<-"eDNA_db/Red_sea_surface/RS_18smini_BGDBR/RS_18smini_cut"
#fastqfile<-"eDNA_db/Red_sea_surface/RSS_CO1/RS_CO1_cut"
#fastqfile<-"eDNA_db/Red_sea_surface/RSS_vert/RS_vert_cut"
#fastqfile<-"eDNA_db/Red_sea_surface/RS_rbclmini_BF9VN/RS_rbclmini_cut"
#
#fastqfile<-"eDNA_db/Arctic_surface/AS_rbclmini_BR3D3/AS_rbclmini_cut"
#fastqfile<-"eDNA_db/Arctic_surface/AS_euka02_BR345/AS_euka02_cut"
#fastqfile<-"eDNA_db/Arctic_surface/AS-2_CO1_CWC3J/AS-2_CO1_cut"
fastqfile<-"eDNA_db/Arctic_surface/AS_2_euka02_CLL7R/AS_2_euka02_CLL7R_cut"

#fastqfile<-"eDNA_db/Dammam_cores/DAM_Euka02_BFB79/DAM_Euka02_cut"
#fastqfile<-"eDNA_db/Dammam_cores/Dammam_core_18SV4/DAM_18sv4_cut"
#fastqfile<-"eDNA_db/Dammam_cores/Dammam_core_CO1/primer_cut_co1_dammam"
#fastqfile<-"eDNA_db/Dammam_cores/DC_18smini/DAM_18smini_cut"
#fastqfile<-"eDNA_db/Dammam_cores/DC_vert/DAM_vert_cut"
#fastqfile<-"eDNA_db/Dammam_cores/DAM_rbclmini_BTWPY/DAM_rbcl_cut"
#fastqfile<-"eDNA_db/Dammam_cores/Dammam_core_Euka02_rep23_J74L5/Dam_euka02_rep23_J74L5_cut"
#
#fastqfile<-"eDNA_db/Arctic_core/Arctic_core_18SV4/AC_18sV4_cut"
#fastqfile<-"eDNA_db/Arctic_core/AC_euka02_BJ34K/AC_euka02_cut"
#fastqfile<-"eDNA_db/Arctic_core/AC_vert/AC_vert_cut"
#fastqfile<-"eDNA_db/Arctic_core/AC_18smini-BTWR2/AC_18smini_cut"
#fastqfile<-"eDNA_db/Arctic_core/AC_rbclmini-BTKHB/AC_rbcl_cut"
#fastqfile<-"eDNA_db/Arctic_core/AC-2_CO1_CWDR5/AC-2_CO1_cut"
#fastqfile<-"eDNA_db/Arctic_core/AC-2_18SV4_CW3MM/AC-2_18SV4_cut"
#fastqfile<-"eDNA_db/Arctic_core/AC-2_euka02_CWBMM/AC-2_euka02_CWBMM_cut"
#
#fastqfile<-"eDNA_db/Atacama/Atacama_CO1_C6RMT/Atacama_CO1_cut_check"
#fastqfile<-"eDNA_db/Atacama/Atacama_vert_CJPDF/Atacama_vert_cut"
#fastqfile<-"eDNA_db/Atacama/Atacama_euka02_C6RNK/Atacama_euka02_cut"

#  Location for summary table --- need to create folder !!!!!
out_file<-"Documents/KAUST/eDNA/R/pipe_summary"

# sample data  !! can skip - not used til post process
#sample<-openxlsx::read.xlsx("C:/Users/geraldn/Dropbox/Documents/KAUST/eDNA/Samples_Data/Extraction test aug18/EXPERIMENT DESIGN 2018 PCR TEMPLATE.xlsx",
 #                               startRow=1, sheet=1, check.names=T)

#################################################################################
#################################################################################
#################  get primer specific information
# euka02
if (primer_name=="euka02" | primer_name=="euka03") {
trunclength_FandT<-c(105,105)
final_trim<-seq(100,150)# peak length 108-110
insect_ref<-"euka02_learn.rds"
#dada_ref<-"SILVA_138_trimmed_euka02_dada2_names.fasta"
   dada_ref<-"SILVA_138_trimmed_euka02_dada2_names_with_ortega_macroalgae.fasta"
#  dada_ref<-"SILVA_138_trimmed_euka02_dada2_names_with_ortega_bachmann_macroalgae.fasta"
}

# co1
if (primer_name=="co1") {
trunclength_FandT<-c(260,180)
final_trim<-seq(280,355) # peak length 313
insect_ref<-"CO1_marine_from_insect.rds"
dada_ref<-"MIDORI_LONGEST_GB239_CO1_trimmed_co1_dada2_names.fasta"
}
# rbclmini
if (primer_name=="rbclmini") {
  trunclength_FandT<-c(180,180)
  #  trunclength_FandT<-c(170,160)
  final_trim<-seq(180,240) # peak length 300
  insect_ref<-"minirbcl_learn.rds"
  dada_ref<-"ncbi_rbcl_trimmed_rbclmini_ncbi_names_rbclmini_dada2_names.fasta"
  }
# 18smini
if (primer_name=="18smini") {
  trunclength_FandT<-c(150,150)
  final_trim<-seq(150,200)# peak length 165
  insect_ref<-"18smini_learn.rds"
  dada_ref<-"SILVA_138_trimmed_18smini_dada2_names.fasta"
  }
# vert
if (primer_name=="vert") {
  trunclength_FandT<-c(120,110)
  final_trim<-seq(115,140)# peak length 123
  insect_ref<-"12s_ruiz_learn.rds"
  dada_ref<-"ncbi_vert_trimmed_ncbi_names_vert_dada2_names.fasta"
  #dada_ref<-"ncbi_12s_euk_only_dada2_mock_virtualPCR.fasta"
  }
#  18s_stoeck
if (primer_name=="18s_stoeck") {
  trunclength_FandT<-c(280,200)
  final_trim<-seq(331,431) # peak 381,384
  insect_ref<-"18s_stoeck_learn.rds"
  dada_ref<-"SILVA_138_trimmed_18s_stoeck_dada2_names.fasta"
  }
#################  get computer specific information
if (computer=="mac") {
  dir<-"/Users/nathangeraldi/Dropbox/"
  multithreadTF<-TRUE}
### kaust windows
if (computer=="win") {
  dir<-"C:/Users/geraldn/Dropbox/"
  multithreadTF<-FALSE}
########################## load or reopen workspace
setwd(paste(dir,out_file,sep=""))
#     load(workspname)
#     save.image(workspname)
##    setting file names -- should need to change

###  set project folder location
workspname<-paste(out_file,"/" , study_name,"__", primer_name,".rdata", sep="")
out_name<-paste(study_name,primer_name,sep="__")

###  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##  if project has more than one miseq per primer pair, use next
# this is old need to just add one number to letter above last, becuase need to sort alphabetically - does not work with adding number !!!!!
# out_name<-paste(study_name,primer_name,"2",sep="__")
#if (replicate=="yes")  {
#workspname<-paste(out_file,"/",study_name,"__","euka03",".rdata", sep="")
#out_name<-paste(study_name,"euka03",sep="__")
#}


### file for filtered fastq
filter_folder<-paste("filtered",out_name)  #

###path to fastq files
path<-paste(dir,fastqfile, sep="")

## path to filtered files
filt_path<-paste(dir,filtfile, sep="")

###path_out for summaries
path_out<-paste(dir,out_file, sep="")

##  set working directory baed on path   
setwd(path)


################################################################################################################
################################################################################################################

################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
## get list of files - for ncbi deposit

#  path<-"/Users/geraldn/Dropbox/eDNA_db/Dammam_cores/Dammam_core_18SV4/Lane1/version_01"
fns <- list.files(path)  
 mes<-data.frame(fns) %>% 
    mutate(num=as.numeric(gsub("_.*$", "", fns))) %>% 
 #  arrange(num) %>% 
   filter(grepl("R2",fns))

 
 
 ################################################################################################################
 ################################################################################################################
 
 
 
################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
##        begin pipline     ######
## create list of files
#  path<-"/Users/geraldn/Dropbox/eDNA_db/Dammam_cores/DAM_Euka02_BFB79/Lane1/version_01"
fns <- list.files(path)            
##########    filter and trim   #################################################################
fnFs <- sort(list.files(path, pattern="_R1_cut.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_cut.fastq", full.names = TRUE))
## extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)  # use 1,  clean up miseq double names, use later on
## plot errors   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#    
#    plotQualityProfile(fnRs[10])
#################         begin filtering and trimming
# create empty folers
filt_path2 <- file.path(filt_path, filter_folder) # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path2, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path2, paste0(sample.names, "_R_filt.fastq.gz"))
# filter
#   if need to change - trunclength_FandT<-c(210,180)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=trunclength_FandT,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=T,
                     compress=TRUE, multithread=multithreadTF) #
#  maxEE- default 2,2 , uses Q score, EE = sum(10^(-Q/10)), sum(10^(-5/10))
# could cut primers here - trimLeft = 0, trimRight = 0,  trim-left-f 17 and trim-left-r 21
#    head(out) 
### learn error rates
# if sample has no reads pass filter
seq_end<-length(filtFs)  # c(1:21,23:seq_end)
errF <- learnErrors(filtFs[1:seq_end], multithread=TRUE)
errR <- learnErrors(filtRs[1:seq_end], multithread=TRUE)
###  plot
plotErrors(errF, nominalQ=TRUE)
###  dereplicate #####################################
## need to correct length if one smample does not have sequences - look at out
#  length()
#  seq_end<-length(filtFs)-1
seq_end<-length(filtFs)
derepFs <- derepFastq(filtFs[1:seq_end], verbose=TRUE)
derepRs <- derepFastq(filtRs[1:seq_end], verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names[1:seq_end]
names(derepRs) <- sample.names[1:seq_end]  # c(1:21,23:seq_end)
##  sample inference- read correction
dadaFs <- dada(derepFs, err=errF, multithread=T, pool=TRUE)#, pool=TRUE
dadaRs <- dada(derepRs, err=errR, multithread=T, pool=TRUE)
###  dadaFs[[1]]
#####  merge  #############################################
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
seqtab <- makeSequenceTable(mergers)
dim(seqtab)  # head(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
# reduce to more specific length
final_trim
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% final_trim]
dim(seqtab2)
table(nchar(getSequences(seqtab2)))
## remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
# see percen removed
sum(seqtab.nochim)/sum(seqtab)
###################################################################
##### produce summary table
# sum(getUniques(out))     sum(getUniques(out))
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
hist(log(colSums(seqtab.nochim)))  # number of reads per unique sequence
###############################################################################
###################  save and read
###############################################################################
write.table(track, paste(path_out,"/", out_name ,"_summary.csv",sep=""),row.names=T, sep=",")
saveRDS(seqtab.nochim, paste(path_out,"/",out_name,"_seqtab.rds",sep=""))



################################################################################################################
################################################################################################################
####
####
####      Clear environment (need to rerun 1-153) or   Close and re-open R      ####################
###
#####
################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################
########     assign taxonomy begin             #################################################################
################################################################################################################
###
###  set up paths and folders if starting here
####################  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  #######################
#  ############# !!!!!!    run lines 1 to 172     !!!!!!!     ###
####
##########################################################################################
##########            get seqtab
setwd(path_out)
#   out_name<-"Arctic_surface__euka02"
seqtab <- readRDS(paste(path_out,"/",out_name,"_seqtab.rds",sep=""))  ## this is pooled nochim  see line 141

#      mes<-colSums(seqtab)    # mes[1:100]

##########################################################################################################
######   for insect      !!!!!!!!!!!!!!!!!
#################################################################################################
##   get reference dendogram
insect_ref_path<-"C:/Users/geraldn/Dropbox/eDNA_db/reference_data/insect_learn"  # setwd(insect_ref_path)
ff <- list.files(insect_ref_path)   #  list.files
ref<-readRDS(paste(insect_ref_path,"/",insect_ref,sep=""))  # ref<-readRDS("18smini_learn.rds")

#################################################
# alter seqtab table to fit insect
x <- char2dna(colnames(seqtab))
## name the sequences sequentially
names(x) <- paste0("ASV", seq_along(x))    # head(x)
#############################################
##  assign taxonomy by insect
# default threshold= 0.9, decay= TRUE
#  !!!! specify levels same as dada2 ????
tax <- classify(x, ref, threshold = 0.5, decay = TRUE, cores = "autodetect", ping = 1,
                ranks = c("superkingdom","kingdom", "phylum", "class", "order", "family", "genus",
                          "species"))
#   to do ?? to match dada2 rename columns\
#    ???tax<-tax 
## combine to take quick loook    names(tax)
tax2 <- tax %>%
    group_by(taxID,superkingdom,kingdom,phylum,class,order,family,genus,species,taxon) %>%
    summarize(mean=mean(score), n=n())
  
#
setwd(path_out)
write.table(tax, paste(out_name,"_taxass_insect.csv", sep=""), sep=",")
#
###############################################################################
###############################################################################
####  assign using dada2   !!!!!!!!!!!!!!!!!!!!!!!
#############################################################
##   get reference fasta and fix names
path_ref<-paste(dir,"eDNA_db/reference_data/dada_final", sep="")
# reference file name
# dada_ref<-"SILVA_138_trimmed_euka02_dada2_names_with_ortega_bachmann_macroalgae.fasta"
ref<-dada_ref
#
list.files(path_ref)
setwd(path_ref)  #
#  check library
#      lib<- read.fasta(file = ref, as.string = TRUE, strip.desc=TRUE)  # 
#fast1<-getAnnot(lib)##     head(lib)  ## psych::describe(getLength(lib))
##
set.seed(100) # Initialize random number generator for reproducibility
##   try  tryRC=TRUE  # reverse of ?  also taxLevels=c()  default minBoot=50
taxass50<- assignTaxonomy(seqtab, ref, multithread = T, verbose = T, minBoot=50,
                          taxLevels=c("Superkingdom","Kingdom","Phylum", "Class", "Order", "Family", "Genus","Species"))  
taxass70<- assignTaxonomy(seqtab, ref, multithread = T, verbose = T, minBoot=70,
                          taxLevels=c("Superkingdom","Kingdom","Phylum", "Class", "Order", "Family", "Genus","Species")) 
taxass90<- assignTaxonomy(seqtab, ref, multithread = T, verbose = T, minBoot=90,
                          taxLevels=c("Superkingdom","Kingdom","Phylum", "Class", "Order", "Family", "Genus","Species"))
#unname(head(taxass))
setwd(path_out)
#  !!!!  need to check out put !!!!!!!!!!!!!!!!!
#  out_name<-"Arctic_surface__euka02_with_ortega_and_Bachm"
write.table(taxass50, paste(out_name,"_taxass50.csv", sep=""), sep=",",row.names=T, col.names = NA)
write.table(taxass70, paste(out_name,"_taxass70.csv", sep=""), sep=",",row.names=T, col.names = NA)
write.table(taxass90, paste(out_name,"_taxass90.csv", sep=""), sep=",",row.names=T, col.names = NA)

taxall<-data.frame(taxass70)
mess<-taxall %>% 
  # dplyr::filter(Phylum=="Chordata") #%>% 
  dplyr::filter(complete.cases(Species)) %>% 
  dplyr::group_by_all() %>% 
  summarise(n())

###############################################################################
###############################################################################

###############################################################################
###############################################################################

###############################################################################
###############################################################################
####    mess ------ old code left behind, do not use
###############################################################################
setwd(paste(dir,"Documents/KAUST/eDNA/R/export/dada",sep=""))
save.image("dada_kapser_tax.rdata")
# Close R, Re-open R
setwd(paste(dir,"Documents/KAUST/eDNA/R/export/dada",sep=""))
load("dada.rdata")


############################## put together



taxonlist<-as.data.frame(taxseq)
mess<- taxonlist %>% 
  group_by(Class) %>% 
  summarise(n())



setwd(paste(dir,"eDNA_db/reference_data/dada_final", sep=""))
#  run lines 1-85 for set up

taxseq<- seqinr::read.fasta(file =ref, as.string = TRUE, strip.desc=TRUE)  # # psych::describe(getLength(taxseq))
fast1<-getAnnot(taxseq)##  get only anotation data  
head(fast1)  # head(mess)
#fastg<-gsub(" .*$","",fast1)  ## keep only before first " "   head(fastg)  #  use for NCBI
fastg<-gsub("\\..*$","",fast1)  ## keep only before first "."  ### use for silva
seq1 <- getSequence(taxseq) # head(seq1)
