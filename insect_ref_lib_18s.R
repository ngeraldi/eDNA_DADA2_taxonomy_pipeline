########
# script to make reference library for insect pipeline
#  for 18s primers  
# based on Siva SSU 132 database

# 3 main steps
#  1-insilica primer cuts
#  2- Prep libraries with dada2 lineage and fro insect dendogram format
#  3- make ref dendogram for insect

# LAST RUN  ---    March 29  2019

##     devtools::install_github("shaunpwilkinson/insect") 
library("insect")

# set computer to get correct directory,  "win"  or "mac"
computer<-"win"      #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## ref library - within dir below  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fastqfile<-"eDNA_db/reference_data/silva"
refname<-"SILVA_132_SSURef_Nr99_tax_silva_trunc_DNA.fasta"
# refname<-"SILVA_132_SSURef_Nr99_tax_silva_trunc_DNA_dada_names.fasta"  # add dada names after trim !!

#################  get computer specific information
if (computer=="mac") {
  dir<-"/Users/geraldn/Dropbox/"
  multithreadTF<-TRUE}
### kaust windows
if (computer=="win") {
  dir<-"C:/Users/geraldn/Dropbox/"
  multithreadTF<-FALSE}
###path to fastq files
path<-paste(dir,fastqfile, sep="")
##  set working directory baed on path
setwd(path)
##########################################################################################
##########################################################################################
##### 1-  perform in silica test of primers
##########################################################################################
###  get reference library
x<-readFASTA(file=paste(dir,fastqfile,"/",refname, sep=""))
##########################################################################################
### set primers  !!!!!!!!!!!!!!!!!^^^^^^^^^^^^^^^^^!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# enter primer name,  "co1" , "euka02"  , "rbclmini"  ,  "18smini"  ,  "vert" , "18s_stoeck" 
primer_name<-"18s_stoeck" 
#################  get primer specific information
# euka02
if (primer_name=="euka02") {
  trunclength_FandT<-c(105,105)
  final_trim<-seq(105,150)# peak length 108-110
  maxlen<-200
  minlen<-90
  insect_ref<-"euka02_learn.rds"
  dada_ref<-"SILVA_132_SSURef_Nr99_tax_silva_trunc_DNA_trimmed_euka02_dada_names.fasta"
  PRIMER1='TTTGTCTGSTTAATTSCG'
  PRIMER2='CACAGACCTGTTATTGC' 
}
# 18smini
if (primer_name=="18smini") {
  trunclength_FandT<-c(150,150)
  final_trim<-seq(150,190)# peak length 165
  maxlen<-200
  minlen<-120
  insect_ref<-"18smini_learn.rds"
  dada_ref<-"silva_18smini_bothprim_trimmed_dada_names.fasta"
  PRIMER1="CCCTGCCHTTTGTACACAC"
  PRIMER2="CCTTCYGCAGGTTCACCTAC"
}

#  18s_stoeck
if (primer_name=="18s_stoeck") {
  trunclength_FandT<-c(280,200)
  maxlen<-500
  minlen<-300
  final_trim<-seq(331,431) # peak 381,384
  insect_ref<-"18suniV4_learn.rds"
  dada_ref<-"silva_18suniV4_bothprim_dada_names.fasta"
  PRIMER1="CCAGCASCYGCGGTAATTCC"   
  PRIMER2="ACTTTCGTTCTTGATYRA"  
}

##   18s V9 - Tara    72839  Amaral-Zettler et al. 2009
#PRIMER1="TTGTACACACCGCCC"
#PRIMER2="CCTTCYGCAGGTTCACCTAC"
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
################################################
# run in silica and trim primers
trimmedf<-virtualPCR(x, up = PRIMER1, down = PRIMER2, rcdown = TRUE, trimprimers = TRUE,
           minfsc = 50, minrsc = 50, minamplen = minlen, maxamplen = maxlen,
           maxNs = 0.02, partialbind = TRUE, cores = "autodetect", quiet = FALSE)
# set working directory     head(trimmedf)
setwd(path)

## export to fasta
fname<-paste("SILVA_132_SSURef_Nr99_tax_silva_trunc_DNA_trimmed_",primer_name,"_silva_names.fasta",sep="")
writeFASTA(trimmedf, file=fname)

################################################  export lineage data
#  data.table::fwrite(tax2,file=paste(path,"/","Silva132_dada_and_silva_lineage.csv", sep=""),row.names=F,sep=",")
##########################################################################################
##########################################################################################
##### 2 -  get dada2 and  insect::learn format   
##########################################################################################
############################################################################
# usual start
library(ape)
library('taxize')
library(ShortRead); packageVersion("ShortRead")
library(seqinr)
library("taxonomizr")
library("dplyr")
#
# with silva 18s
dir<-"C:/Users/geraldn/Dropbox/"
path<-paste(dir,"eDNA_db/reference_data/silva", sep="")
file<-list.files(path)
#  run lines 1-85 for set up
fname<-paste("SILVA_132_SSURef_Nr99_tax_silva_trunc_DNA_trimmed_",primer_name,"_silva_names.fasta",sep="")

taxseq<- seqinr::read.fasta(file =paste(path,"/",fname, sep=""), as.string = TRUE, strip.desc=TRUE)  # # psych::describe(getLength(taxseq))
fast1<-getAnnot(taxseq)##  get only anotation data  
head(fast1)  # head(mess)
#fastg<-gsub(" .*$","",fast1)  ## keep only before first " "   head(fastg)  #  use for NCBI
fastg<-gsub("\\..*$","",fast1)  ## keep only before first "."  ### use for silva
seq1 <- getSequence(taxseq) 
setwd("/Users/geraldn/Dropbox/eDNA_db/reference_data/NCBI_tax_taxonomizr")   # getwd()
taxon_id<-accessionToTaxa(fastg,"accessionTaxa.sql",version='base') # head(taxon_id)
tax2<-data.frame(cbind(fastg,taxon_id))
names(tax2)[1]<-"accession"
names(tax2)[2]<-"taxID"    # head(df)
tax2<-dplyr::mutate_all(tax2,as.character)
tax<-getTaxonomy(tax2$taxID, "accessionTaxa.sql", desiredTaxa = c("superkingdom",
        "phylum", "class", "order", "family", "genus", "species"), mc.cores = 1,debug = FALSE)
tax<-data.frame(tax,stringsAsFactors = F)
############  get lineages
tax3<- tax2 %>%     # head(tax3)
  mutate(num=as.numeric(row.names(tax2))) %>%  
  bind_cols(tax) %>%   
  tidyr::unite(dada2_lineage, superkingdom:species,sep=";",remove=F) %>%
  tidyr::unite(insect_lineage, accession,taxID, sep="|",remove=F) %>% 
  filter(complete.cases(accession)) %>% 
  filter(complete.cases(taxID)) %>% 
  mutate(taxID=as.numeric(taxID)) %>% 
  filter(taxID>0)  # 0 is unknown - need to remove
### check taxIDs against insect ncbi database
db<-insect::taxonomy(db = "NCBI", synonyms = F) # download ncbi taxonomy database directly (about 2 min) or can save and use below
## remove taxIDs not in db    head(tax3$num)
tax3 <-tax3 %>%
    inner_join(db)
# export fasta
  setwd(path)     # getwd()
write.fasta(sequences = seq1[tax3$num], names = tax3$insect_lineage, file.out = "silva_18suniV4_bothprim_insect_names.fasta", open = "w")  #
# and 
write.fasta(sequences = seq1[tax3$num], names = tax3$dada2_lineage, file.out = "silva_18suniV4_bothprim_dada_names.fasta", open = "w")  #

##########################################################################################
##########################################################################################
##### 3- make insect learned tree
##########################################################################################
############################################################################
#  steps after using virtualPCR
# also need to use "taxonomy" to download ncbi taxonomy
# 1- get_lineage - need taxaIDS
# 2- learn   -  fasta  names- AF296345|8017 (accession|taxID) (|salmon)  and sequences
# use #2 to classify

library("insect")
# path<-
setwd(path)  
x<-insect::readFASTA(file ="silva_18suniV4_bothprim_insect_names.fasta")
db<-insect::taxonomy(db = "NCBI", synonyms = F) # download ncbi taxonomy database directly (about 2 min) or can save and use below
# head(db)

insect_learn<- learn(x, db, model = NULL, refine = "Viterbi", iterations = 50,
      nstart = 20, minK = 2, maxK = 2, minscore = 0.9, probs = 0.5,
      retry = TRUE, resize = TRUE, maxsize = max(sapply(x, length)),
      recursive = TRUE, cores = "autodetect", quiet = F)

# Found 36315 unique sequences for Euka02 Nov3, 2018
# Found 24347 unique sequences for 18smini Nov12, 2018
# Found 54595 unique sequences for 18suniV4 Jan29, 2019


setwd("C:/Users/geraldn/Dropbox/eDNA_db/reference_data/insect_learn")
saveRDS(insect_learn,"18suniV4_learn.rds")  ## then readRDS

