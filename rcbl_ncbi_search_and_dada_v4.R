######################################################################################################
######################################################################################################

#  this code steps through downloading rbcl sequences from NCBI 
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


###    1- get sequences
######################################################################################################
######################################################################################################
#######   get data   ##############################################
### kaust PC
dir<-"C:/Users/geraldn/Dropbox/"
#   mac
dir<-"/Users/nathangeraldi/Dropbox/"

######################################################################################################
######################################################################################################
######################################################################################################
#########################'#############################################################################
#####        BOLD  -- don't use !!!!!!  ##########################################################################
#  if useing Bold??  dont because not good id number or taxonomy and not many 1912 total  1523 RBCL
path<-paste(dir,"eDNA_db/reference_data/bold", sep="")
setwd(path)   # getwd()
dat<-read.table("iBOL_phase_6.50_plants.tsv", sep="\t", header=TRUE, quote = "")  #  names(dat)
mes<- dat   %>%    # 
  filter(grepl("rbcLa", processid)) %>% 
  filter(grepl("Alismatal", order_reg)) 
### can not find right id to match to taxonomy????
id<-7985671
mes3<-bold_tax_id(id, dataTypes = "basic", includeTree = FALSE,
            response = FALSE)
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
### use EMBL plant database    -- don't use  !!!!!!
##  very bad format, nothing on web to import or convert in R
path<-paste(dir,"eDNA_db/reference_data/embl/plants", sep="")
setwd(path)   # getwd()

# start loop for all files
dat<-read.table("rel_std_pln_17_r134.dat", header=F, quote = "", skip=3, n=10)  #  names(dat)
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
##############################
###  Download from  ncbi   !!!   best  !!!!!!
####################################################################################
###   rcbl
path<-paste(dir,"eDNA_db/reference_data/ncbi_rbcl", sep="")
setwd(path)   # getwd()
#### for rbcl database
###############################################################################################
#####    search ncbi for sequences  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

##### used below  !!!!!!!!!
sch<-"(rbcl[All Fields] AND ((Embryophyta[Organism] OR Plants[All Fields] OR Chlorophyta[Organism] OR Phaeophyceae[Organism] OR Rhodophyta[Organism]) AND 00000000001[SLEN] : 00000010000[SLEN]) NOT (Metazoa[Organism] OR animals[All Fields]))"
# copy and paste above into ncbi and put in numr
numr<- 256994 # on Nov 2020, before - Mar 2018 -232193

## set api key ########################################################
rentrez::set_entrez_key("27ecad18cf3950adb97a95bb779e8a510209")
Sys.getenv("ENTREZ_KEY")
#
##    get sequence ids
x<-rentrez::entrez_search(db="nuccore", term = sch, retmax=numr, use_history=TRUE)  # length(x)
#  x$ids #gives you the NCBI ids!
####use web history
x$web_history

######################################################################################################
######################################################################################################
#  download seqeunces
############################
##     !!!!   WILL ADD to EXISTING FILE, so delete old file or change name to old one   !!!!!!
#############################

####   get sequences - downloads to  getwd() ,  need to use loop with ncbi or times out
#seq(1,numr,100)
stp<-1000
stp2<-stp-1
for( seq_start in seq(1,numr,stp)){
  recs <- entrez_fetch(db="nuccore", web_history=x$web_history,
                       rettype="fasta", retmax=stp, retstart=seq_start)
  cat(recs, file="rbcl_embry_plants_1_10k_nov2020.fasta", append=TRUE)  #  ! set name !!!!! or will overwrite
  cat(seq_start+stp2, "sequences downloaded\r")
  # print status to console
  print(seq_start)
  print(Sys.time()) 
}
# head(recs)

######################################################################################################
######################################################################################################
##########################################################################################
########################## 2- virtualPCR   !!!!!!!!!!!!!!!!!!!!!!!!!!!!

###  get reference library
dir<-"/Users/nathangeraldi/Dropbox/"
path1<-paste(dir,"eDNA_db/reference_data/ncbi_rbcl", sep="")
setwd(path1) 
#list.files()
x<- readFASTA(file ="rbcl_embry_plants_1_10k_nov2020.fasta")  # # psych::describe(getLength(taxseq))
####    minirbcl    head(x)  summary(x)
PRIMER1="GTTGGATTCAAAGCTGGTGTTA"
PRIMER2="CVGTCCAMACAGTWGTCCATGT"
maxlen<-500
minlen<-100

################################################
# run in silica and trim primers
trimmedf<-virtualPCR(x, up = PRIMER1, down = PRIMER2, rcdown = TRUE, trimprimers = TRUE,
                     minfsc = 50, minrsc = 50, minamplen = minlen, maxamplen = maxlen,
                     maxNs = 0.02, partialbind = TRUE, cores = "autodetect", quiet = FALSE)

## export to fasta
# set working directory
setwd(path1)   #  getwd()

## export to fasta
writeFASTA(trimmedf, file="ncbi_rbcl_trimmed_ncbi_names.fasta")



##########################################################################################
##########################################################################################
#####   !!!!!!!!!!!!!!!!!!!!!!1 STOP  all below is old  !!!!!!!!!!!!!!!!!!!!!!!!!
#### go to create_ref_lib_18s to clean names and format reference databases
##########################################################################################
############################################################################



######################################################################################################
######################################################################################################



######################################################################################################
######################################################################################################




######################################################################################################
######################################################################################################
####################################################################################
##########################################################################################
##### 3 - Clean sequences for good taxonomy then get dada2 and  insect::learn format   
##########################################################################################
############################################################################
# usual start
library(ape)
library('taxize')
library(ShortRead); packageVersion("ShortRead")
library(seqinr)
library("taxonomizr")
library("dplyr")

##### reload fasta file and get lineage
dir<-"C:/Users/geraldn/Dropbox/"
path1<-paste(dir,"eDNA_db/reference_data/ncbi_rbcl", sep="")
setwd(path1) 
#list.files()
fast<- read.fasta(file ="ncbi_rbcl_virtualPCR_rbclmini_ncbi_names.fasta", as.string = TRUE, strip.desc=TRUE)  # # psych::describe(getLength(taxseq))
fast1<-getAnnot(fast)##  head(dat)  length(dat)
seq1 <- getSequence(fast)  #  head(seq1)    psych::describe(getLength(fast))
#  head(fast)
fastg<-gsub("\\..*$","",fast1)  ## keep only before first "."      head(dat1) # use for ncbi
#fastg<-gsub(" .*$","",fast1)  ## keep only before first " "   head(fastg)  #  use for silva
########
setwd("/Users/geraldn/Dropbox/eDNA_db/reference_data/NCBI_tax_taxonomizr")
taxon_id<-accessionToTaxa(fastg,"accessionTaxa.sql",version='base') # head(taxon_id)
tax2<-data.frame(cbind(fastg,taxon_id))
names(tax2)[1]<-"accession"
names(tax2)[2]<-"taxID"    # head(df)
tax2<-dplyr::mutate_all(tax2,as.character)
tax<-getTaxonomy(tax2$taxID, "accessionTaxa.sql", desiredTaxa = c("superkingdom","kingdom",
                                                                  "phylum", "class", "order", "family", "genus", "species"), mc.cores = 1,debug = FALSE)
tax<-data.frame(tax,stringsAsFactors = F)
############  make data base of  lineages  and filter bad entries
#Set un-wanted terms
unwanteds <- c("fungal","endophyte","symbiont","uncultured", 
               "unclassified", "unidentified", "metagenome",
               "predicted", "environmental", "unknown")#
unwanteds <- paste0(unwanteds, collapse = "|")

tax3<- tax2 %>%     # head(tax3)    names(tax3)
  mutate(num=as.numeric(row.names(tax2))) %>%  
  bind_cols(tax) %>%   
  tidyr::unite(dada2_lineage, superkingdom:species,sep=";",remove=F) %>%
  tidyr::unite(insect_lineage, accession,taxID, sep="|",remove=F) %>% 
  filter(complete.cases(accession)) %>% 
  filter(complete.cases(taxID)) %>% 
  mutate(taxID=as.numeric(taxID)) %>% 
  filter(taxID>0) %>%  # 0 is unknown - need to remove
  filter_at(vars(genus,species), any_vars(complete.cases(.)))    %>%  # no species or genus
  filter_at(vars(phylum:genus), any_vars(complete.cases(.)))    %>%  # no phylum through genus
  filter(!grepl(unwanteds, dada2_lineage, fixed = FALSE, ignore.case = TRUE)) %>%  # remove unwanted terms
  filter(!duplicated(accession))   # remove duplicated accession numbers - silva should not have these, what -  must after decimal ?????????????

# sanity check
mes<-tax3 %>% 
  filter(grepl("labrax", dada2_lineage, fixed = FALSE, ignore.case = TRUE))
### check taxIDs against insect ncbi database  ########################
#setwd("C:/Users/geraldn/Dropbox/eDNA_db/reference_data/NCBI_tax_insect")
#db <- readRDS("C:/Users/geraldn/Dropbox/eDNA_db/reference_data/NCBI_tax_insect/db_insect.rds")  
db<-insect::taxonomy(db = "NCBI", synonyms = F) # download ncbi taxonomy database directly (about 2 min) or can save and use below
## remove taxIDs not in db    head(tax3$num)
tax3 <-tax3 %>%
  inner_join(db)
###   re-make orignal DNAbin to match filter sequences
#    head(x2)  names(x2)
seq2<-seq1[tax3$num]
nnnn<-as.character(tax3$accession) # head(nnnn)  mes<-nnnn[duplicated(nnnn)]
names(seq2) <- c(nnnn)

## purge suspect sequences
z<-tax3$taxID
tdb <- prune_taxonomy(db, taxIDs = z)
set.seed(999)
#  rm(db)  rm(fast1)
zz <- purge(seq2, tdb, cores=18) # none to remove


# export fasta
setwd(path1)     # getwd()
primer_name<-"rbclmini"
fname_out_insect<-paste("ncbi_rbcl_trimmed_",primer_name,"_insect_names.fasta",sep="")
fname_out_dada<-paste("ncbi_rbcl_trimmed_",primer_name,"_dada2_names.fasta",sep="")
write.fasta(sequences = seq1[tax3$num], names = tax3$insect_lineage, file.out = fname_out_insect, open = "w")  #
# and 
write.fasta(sequences = seq1[tax3$num], names = tax3$dada2_lineage, file.out = fname_out_dada, open = "w")  #


##########################################################################################
##########################################################################################
##### 4- make insect reference databases
##########################################################################################
############################################################################

library("insect")
dir<-"C:/Users/geraldn/Dropbox/"
path1<-paste(dir,"eDNA_db/reference_data/ncbi_rbcl", sep="")
setwd(path1)  
fname_out_insect<-"ncbi_rbcl_trimmed_rbclmini_insect_names.fasta"
x<-insect::readFASTA(file = fname_out_insect)
db<-insect::taxonomy(db = "NCBI", synonyms = F) # download ncbi taxonomy database directly (about 2 min) or can save and use below
# head(db)

rbcl_learn<- learn(x, db, model = NULL, refine = "Viterbi", iterations = 50,
                     nstart = 20, minK = 2, maxK = 2, minscore = 0.9, probs = 0.5,
                     retry = TRUE, resize = TRUE, maxsize = max(sapply(x, length)),
                     recursive = TRUE, cores = "autodetect", quiet = F)

# Found 12669 unique sequences for rbcl april 18th, 2019

setwd("C:/Users/geraldn/Dropbox/eDNA_db/reference_data/insect_learn")
saveRDS(rbcl_learn,"minirbcl_learn.rds")  ## then readRDS






##  if want some additional info.
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
##############################################################################################################
### import different fasta and make csv with yes for minirbcl

dir<-"C:/Users/geraldn/Dropbox/"
path<-paste(dir,"eDNA_db/reference_data/ncbi_rbcl", sep="")
setwd(path)
all<- seqinr::read.fasta(file =paste(path,"/","ncbi_rbcl_all_dada_names.fasta", sep=""), as.string = TRUE, strip.desc=TRUE) 
ins_trim<- seqinr::read.fasta(file =paste(path,"/","ncbi_rbcl_virtualPCR_rbclmini_insect_names.fasta", sep=""), as.string = TRUE, strip.desc=TRUE) 
#   get annote
all<-seqinr::getAnnot(all)##      head(all2)
ins_trim2<-seqinr::getAnnot(ins_trim)##     head(ins_trim3)
ins_trim3<-data.frame(do.call(rbind,ins_trim2))
ins_trim3$mm<-"yes"
names(ins_trim3)[1]<-"lineage"

#
#all<-gsub("\\..*$","",all)  # get accesssion number before first period
#ins_trim<-gsub("\\..*$","",ins_trim)  #  head(ins_trim)
#t18seuka02<-gsub("\\..*$","",t18seuka02)  #
### get csv with tax ids     getwd()  list.files()
dir<-"C:/Users/geraldn/Dropbox/"
path1<-paste(dir,"eDNA_db/reference_data/ncbi_rbcl", sep="")
setwd(path1) # getwd() 
x<-data.table::fread(file="ncbi_rbcl_lineage.csv",header = T,sep=",")


x2<-x %>%
  unite(lineage_insect, accession:taxsID,sep="|",remove=F)  %>%
  dplyr::mutate(rbclmini_amp="no") %>%
  dplyr::mutate(rbclmini_amp=replace(rbclmini_amp,lineage_insect %in% ins_trim3$lineage, "yes" ))   %>%
  filter(lineage_insect !=";;;;;;")

mes<-x2[x2$rbclmini_amp=="yes",]
################################################  export lineage data
data.table::fwrite(x2,file=paste(path,"/","ncbi_rbcl_dada_lineage_rbclmini_insilico.csv", sep=""),row.names=F,sep=",")








