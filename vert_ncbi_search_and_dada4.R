######################################################################################################
######################################################################################################


#  this code steps through downloading 12S sequences from NCBI 
## then....  run code to clean names and cross reference with worrms - create_ref_lib_18s 
# it primarily uses rentrez, taxonomizer and insect packages
##  following steps should be followed
## find search using ncbi website and copy exact search terms then.....
# 1- get sequences
#  2- virtual PCR
# 3- switch to other script - create_ref_lib_18s 

##     last run Nov 2020   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

######################################################################################################
######################################################################################################
library(ape)
library('taxize')
library(ShortRead); packageVersion("ShortRead")
library(seqinr)
library(dplyr)
library(tidyr)
library("bold")
library("taxonomizr")
library("rentrez")  #  devtools::install_github("ropensci/rentrez")

############################################################################
# usual start - should not need now, only once
#######   get data   ##############################################
### kaust PC
dir<-"C:/Users/geraldn/Dropbox/"
#   mac
dir<-"/Users/nathangeraldi/Dropbox/"
######################################################################################################
######################################################################################################
###  try download ncbi
####################################################################################
###   vert
path<-paste(dir,"eDNA_db/reference_data/ncbi_vert", sep="")
setwd(path)   # getwd()
#### for rbcl database
###############################################################################################
#####    search ncbi for sequences  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
### run in ncbi website search to make sure search terms are correct and get number of records-numr
######################################################################################################
######################################################################################################
#  download seqeunces
############################
##     !!!!   WILL ADD to EXISTING FILE, so delete old file or change name to old one   !!!!!!
#############################
##### used below  !!!!!!!!!
sch<-"(Vertebrata[Organism] OR Vertebrata[Organism] OR vertebrata[All Fields]) AND (Metazoa[Organism] OR animals[All Fields]) AND 12s[All Fields]) AND 00000000001[SLEN] : 00000100000[SLEN]"
numr<-179159

## set api key ########################################################
rentrez::set_entrez_key("27ecad18cf3950adb97a95bb779e8a510209")
Sys.getenv("ENTREZ_KEY")
#
x<-rentrez::entrez_search(db="nuccore", term = sch, retmax=numr, use_history=TRUE)  # length(x)
#  x$ids #gives you the NCBI ids!
####use web history
x$web_history
####summ <- entrez_summary(db="nuccore", web_history=x$web_history) #doesnt' work
#seq(1,numr,100)
stp<-1000
stp2<-stp-1
for( seq_start in seq(1,numr,stp)){
  recs <- entrez_fetch(db="nuccore", web_history=x$web_history,
                       rettype="fasta", retmax=stp, retstart=seq_start)
  cat(recs, file="12s_vert_1_100k_nov2020.fasta", append=TRUE)  #  ! set name !!!!! or will overwrite !!!!!!!!!!!!!!!!!!!!
  cat(seq_start+stp2, "sequences downloaded\r")
  # print status to console
  print(seq_start)
  print(Sys.time()) 
}
# head(recs)


###############################################################################################################
####################################################################################################################
################     TRIM   !!!!!!!!!!!!!!!!!!!!! 
# set computer to get correct directory,  "win"  or "mac"

###  get reference library
dir<-"/Users/nathangeraldi/Dropbox/"
fastqfile<-"eDNA_db/reference_data/ncbi_vert/"
file_name="12s_vert_1_100k_nov2020.fasta"
#list.files()
x<- readFASTA(file =paste0(dir,fastqfile,file_name)  )  # # psych::describe(getLength(taxseq))

####   vert
PRIMER1<-"ACTGGGATTAGATACCCC"
PRIMER2<-"TAGAACAGGCTCCTCTAG"


################################################
# run in silica   from 24564 seq to 18075
trimmedf<-insect::virtualPCR(x, up = PRIMER1, down = PRIMER2, rcdown = TRUE, trimprimers = T,
                             minfsc = 50, minrsc = 50, minamplen = 50, maxamplen = 400,
                             maxNs = 0.02, partialbind = TRUE, cores = "autodetect", quiet = FALSE)

## export to fasta
writeFASTA(trimmedf, file= paste0(dir,fastqfile,"ncbi_vert_trimmed_ncbi_names.fasta")  )


# for nov 2020 or 179158 seq, 111331 were retained after trim







##########################################################################################
##########################################################################################
#####   !!!!!!!!!!!!!!!!!!!!!!1 STOP  all below is old  !!!!!!!!!!!!!!!!!!!!!!!!!
#### go to create_ref_lib_18s to clean names and format reference databases
##########################################################################################
############################################################################

#####   !!!!!!!!!!!!!!!!!!!!!!1 STOP  all below is old  !!!!!!!!!!!!!!!!!!!!!!!!!

######################################################################################################
######################################################################################################



######################################################################################################
######################################################################################################




####################################################################################################################
################ use taxonomizr to get ncbi lineage  
##### load fasta file and get lineage
#list.files()
fast<- read.fasta(file ="12s_vert_1_100k_apr_19.fasta", as.string = TRUE, strip.desc=TRUE)  # # psych::describe(getLength(taxseq))
dat<-getAnnot(fast)##  head(dat)  length(dat)
seq <- getSequence(fast)  #  head(seq)
dat1<-gsub("\\..*$","",dat)  ## keep only before first "."      head(dat1)
dat2<-gsub(" .*$","",dat) ## keep only before first " "      head(dat2)
#
setwd("/Users/geraldn/Dropbox/eDNA_db/reference_data/NCBI_tax_taxonomizr")   # getwd()
#
ids<-accessionToTaxa(dat1,"accessionTaxa.sql", version='base')  ## ids[1:20]    dat1[1:100]
#  check number of NA, alwasy a few
summary(is.na(ids)) 
range(which(is.na(ids)))
dat3<-data.frame(cbind(dat1,dat2,ids))  # names(dat3)
names(dat3)[1]<-"accession"
names(dat3)[2]<-"accession2"
names(dat3)[3]<-"taxaID"
# ids2<-accessionToTaxa(dat2[1:100],"accessionTaxa.sql")  #
####   from entrex get id
taxx<-getTaxonomy(dat3$taxaID,"accessionTaxa.sql",desiredTaxa = c("superkingdom","kingdom",
              "phylum", "class", "order", "family", "genus", "species"), mc.cores = 1,debug = FALSE)
tax2<-data.frame(taxx,stringsAsFactors = F)  ## head(tax2)   names(tax2)
tax3<-cbind(dat3,tax2)
summary(is.na(tax3$superkingdom)) 
summary(tax3$phylum=="Chordata") #   unique(tax3$phylum)   unique(tax3$superkingdom)
## export accession numbers and lineage if you want????
tax2<-dplyr::mutate_all(tax3,as.character)
dir<-"C:/Users/geraldn/Dropbox/"
path1<-paste(dir,"eDNA_db/reference_data/ncbi_vert", sep="")
setwd(path1) # getwd() 
data.table::fwrite(tax2,file="ncbi_12s_lineage_2019_april.csv",row.names=F,sep=",")


##############################################################################  \
#####################  limit to euk and clean
dir<-"C:/Users/geraldn/Dropbox/"
path1<-paste(dir,"eDNA_db/reference_data/ncbi_vert", sep="")
setwd(path1)
# if starting here import   getwd()
tax2<-data.table::fread(file="ncbi_12s_lineage_2019_april.csv",header = T, sep=",")
all<- seqinr::read.fasta(file =paste(path1,"/","12s_vert_1_100k_apr_19.fasta", sep=""), as.string = TRUE, strip.desc=TRUE) 
seq<-seqinr::getSequence(all)##   
anot<-seqinr::getAnnot(all) # head(anot)

############  make data base of  lineages  and filter bad entries
#Set un-wanted terms
unwanteds <- c("fungal","endophyte","symbiont","uncultured", 
               "unclassified", "unidentified", "metagenome",
               "predicted", "environmental", "unknown")#
unwanteds <- paste0(unwanteds, collapse = "|")

###  prep for dada2 naming     names(tax2)
tax3<- tax2 %>%
  mutate(num=as.numeric(row.names(tax2))) %>% 
  dplyr::rename(taxID=taxaID)  %>% 
  mutate_all(na_if,"")  %>% 
  #bind_cols(tax) %>%   
  tidyr::unite(dada2_lineage, superkingdom:species,sep=";",remove=F) %>%
  tidyr::unite(insect_lineage, accession,taxID, sep="|",remove=F) %>% 
  filter(complete.cases(accession)) %>% 
  filter(complete.cases(taxID)) %>% 
  mutate(taxID=as.numeric(taxID)) %>% 
  filter(taxID>0) %>%  # 0 is unknown - need to remove
  filter_at(vars(genus,species), any_vars(complete.cases(.)))    %>%  # no species or genus
  filter_at(vars(phylum:genus), any_vars(complete.cases(.)))    %>%  # no phylum through genus
  filter(!grepl(unwanteds, dada2_lineage, fixed = FALSE, ignore.case = TRUE)) %>%  # remove unwanted terms
  filter(!duplicated(accession))   #
  


### check taxIDs against insect ncbi database  ########################
#setwd("C:/Users/geraldn/Dropbox/eDNA_db/reference_data/NCBI_tax_insect")
#db <- readRDS("C:/Users/geraldn/Dropbox/eDNA_db/reference_data/NCBI_tax_insect/db_insect.rds")  
db<-insect::taxonomy(db = "NCBI", synonyms = F) # download ncbi taxonomy database directly (about 2 min) or can save and use below
## remove taxIDs not in db    head(tax3$num)
tax3 <-tax3 %>%
  inner_join(db)
# remove all but Chordata or eukaryotes
tax4<- tax3 %>%        # 
  filter(superkingdom=="Eukaryota" | phylum=="Chordata")
#### export in both insect and dada2 format
# export fasta   fast[1:5]  length(unique(nam$acc))
dir<-"C:/Users/geraldn/Dropbox/"
path1<-paste(dir,"eDNA_db/reference_data/ncbi_vert", sep="")
setwd(path1) # getwd() 
write.fasta(sequences = seq[tax4$num], names = tax4$insect_lineage, file.out = "ncbi_12s_euk_only_insect.fasta", open = "w")  #
write.fasta(sequences = seq[tax4$num], names = tax4$dada2_lineage, file.out = "ncbi_12s_euk_only_dada2.fasta", open = "w")  
###############################################################################
########################################################################################################



##########################################################################################
##########################################################################################
##### 3- make insect reference databases
##########################################################################################
############################################################################
#  steps after using virtualPCR
# also need to use "taxonomy" to download ncbi taxonomy
# 1- learn   -  fasta  names- AF296345|8017 (accession|taxID) (|salmon)  and sequences
# use #2 to classify

library("insect")
dir<-"C:/Users/geraldn/Dropbox/"
path1<-paste(dir,"eDNA_db/reference_data/ncbi_vert", sep="")
setwd(path1)  
x<-insect::readFASTA(file ="ncbi_12s_euk_only_insect_virtualPCR.fasta")  # head(x)
db<-insect::taxonomy(db = "NCBI", synonyms = F) # download ncbi taxonomy database directly (about 2 min) or can save and use below
# head(db)

vert_learn<- learn(x, db, model = NULL, refine = "Viterbi", iterations = 50,
                     nstart = 20, minK = 2, maxK = 2, minscore = 0.9, probs = 0.5,
                     retry = TRUE, resize = TRUE, maxsize = max(sapply(x, length)),
                     recursive = TRUE, cores = "autodetect", quiet = F)

# Found 3730 unique sequences for 12s April 21, 2019


setwd("C:/Users/geraldn/Dropbox/eDNA_db/reference_data/insect_learn")
saveRDS(vert_learn,"12s_ruiz_learn.rds")  ## then readRDS







