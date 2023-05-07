



########
# script to make reference library for insect pipeline
#  for 18s primers  
# based on Siva SSU 132 database

# 5 main steps
# 1-  update SILVA and conver to DNA
# 2-  update taxonomizer to get ncbi naming
#  3-insilica primer cuts
#  4- Prep libraries with dada2 lineage and fro insect dendogram format
#  4.5 - correct taxonomy with WORMS
#  5- make ref dendogram for insect

############################################################################
########################################
## update SILVA  !!!!!!
#   https://www.arb-silva.de/download/archive/
# get most recent release - go to export
#  something like -- SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta.gz
## last updated Nov 2020


############################################################################
####    for taxonomizer - quickest way to get ncbi taxonomy data  #######################################
#  devtools::install_github("sherrillmix/taxonomizr")
#  setwd("/Users/nathangeraldi/Dropbox/eDNA_db/reference_data/NCBI_tax_taxonomizr")
#  taxonomizr::prepareDatabase('accessionTaxa.sql')  # only if want to update data base, need to delete old, 
# LAST RUN  ---    March 29  2019

##     devtools::install_github("shaunpwilkinson/insect") 
library("insect")

# set computer to get correct directory,  "win"  or "mac"
computer<-"mac"      #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## ref library - within dir below  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fastqfile<-"eDNA_db/reference_data/silva"
#  location of taxonomizer
taxonz_file<-"eDNA_db/reference_data/NCBI_tax_taxonomizr" # location of taxonomizr within dir

# refname<-"SILVA_132_SSURef_Nr99_tax_silva_trunc_DNA_dada_names.fasta"  # add dada names after trim !!

#################  get computer specific information
if (computer=="mac") {
  dir<-"/Users/nathangeraldi/Dropbox/"
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
###################  need to  convert from rna to dna - only first time using new SILVA   ###################################################
library(Biostrings)
##   devtools::install_github("GuillemSalazar/FastaUtils")
library(FastaUtils)

#  setwd(path);    RNA2DNA(file ="SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta",out="SILVA_138_SSURef_Nr99_tax_silva_trunc_DNA.fasta")
##########################################################################################
##########################################################################################

##########################################################################################
##### 1-  perform in silica test of primers

### set primers  !!!!!!!!!!!!!!!!!^^^^^^^^^^^^^^^^^!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# enter primer name,  "co1" , "euka02"  , "rbclmini"  ,  "18smini"  ,  "vert" , "18s_stoeck" 
primer_name<-"18s_stoeck" 
#################  get primer specific information
# euka02
if (primer_name=="euka02") {
  trunclength_FandT<-c(105,105)
  final_trim<-seq(70,200)# peak length 108-110
  maxlen<-200
  minlen<-70
  insect_ref<-"euka02_learn.rds"
  dada_ref<-"SILVA_138_SSURef_Nr99_tax_silva_trunc_DNA_trimmed_euka02_dada_names.fasta"
  PRIMER1='TTTGTCTGSTTAATTSCG'
  PRIMER2='CACAGACCTGTTATTGC' 
}
# 18smini
if (primer_name=="18smini") {
  trunclength_FandT<-c(150,150)
  final_trim<-seq(100,220)# peak length 165
  maxlen<-220
  minlen<-100
  insect_ref<-"18smini_learn.rds"
  dada_ref<-"SILVA_138_SSURef_Nr99_tax_silva_trunc_DNA_trimmed_18smini_dada_names.fasta"
  PRIMER1="CCCTGCCHTTTGTACACAC"
  PRIMER2="CCTTCYGCAGGTTCACCTAC"
}

#  18s_stoeck
if (primer_name=="18s_stoeck") {
  trunclength_FandT<-c(280,200)
  maxlen<-500
  minlen<-300
  final_trim<-seq(300,500) # peak 381,384
  insect_ref<-"18suniV4_learn.rds"
  dada_ref<-"SILVA_138_SSURef_Nr99_tax_silva_trunc_DNA_trimmed_18sUniv_dada_names"
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




##########################################################################################
##########################################################################################
##########################################################################################
##########################   begin reference library creations
## STEP 1 - virtualPCR
########################## 

refname<-"SILVA_138_SSURef_Nr99_tax_silva_trunc_DNA.fasta"
###  get reference library
x<-readFASTA(file=paste(dir,fastqfile,"/",refname, sep=""))
##########################################################################################

################################################
# run in silica and trim primers
trimmedf<-virtualPCR(x, up = PRIMER1, down = PRIMER2, rcdown = TRUE, trimprimers = TRUE,
           minfsc = 50, minrsc = 50, minamplen = minlen, maxamplen = maxlen,
           maxNs = 0.02, partialbind = TRUE, cores = "autodetect", quiet = FALSE)
# set working directory     head(trimmedf)
setwd(path)

## export to fasta
fname<-paste("SILVA_138_SSURef_Nr99_tax_silva_trunc_DNA_trimmed_",primer_name,"_silva_names.fasta",sep="")
writeFASTA(trimmedf, file=fname)

# nov 2020, euka -52252 non ambiguous sequences ,18smini - 62741
################################################  export lineage data
#  data.table::fwrite(tax2,file=paste(path,"/","Silva132_dada_and_silva_lineage.csv", sep=""),row.names=F,sep=",")

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


#  !! next verstion add any additional sequences here - move up from 2.5 !!!!!!!!!!!!

##########################################################################################
##########################################################################################
##### STEP  2 -  get worrms names and then ..  dada2 and  insect::learn format   
##########################################################################################
############################################################################
# usual start
library(ape)
library('taxize')
library(ShortRead); packageVersion("ShortRead")
library(seqinr)
library("taxonomizr")
library("tidyverse")
# devtools::install_github("ropensci/worrms")  # don't need to load but will use
#

############################################################################
##  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
#  run lines 37-113 for set up for 18S or next sections for other primers
##  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
#  prepareDatabase('accessionTaxa.sql')  # only if want to update data base, need to delete old 


fname<-paste("SILVA_138_SSURef_Nr99_tax_silva_trunc_DNA_trimmed_",primer_name,"_silva_names.fasta",sep="")
####   import fast, get headers and get taxonomy for ref data format
taxseq<- seqinr::read.fasta(file =paste(path,"/",fname, sep=""), as.string = TRUE, strip.desc=TRUE)  # # psych::describe(getLength(taxseq))


############################################################################
############################################################################
############################################################################
############# for other primers start here, skip if ran lines 37-113

dir<-"/Users/nathangeraldi/Dropbox/"
##
primer_name<-"co1"  #   "co1"  , "rbclmini"  ,   "vert" 
ref_file<-"eDNA_db/reference_data/Midori_meta_COI/"
f_name_file<-"MIDORI_LONGEST_GB239_CO1_trimmed"
##
primer_name<-"rbclmini"   #   "co1"  , "rbclmini"  ,   "vert" 
ref_file<-"eDNA_db/reference_data/ncbi_rbcl/"
f_name_file<-"ncbi_rbcl_trimmed_rbclmini_ncbi_names"
##
primer_name<-"vert"   #  
ref_file<-"eDNA_db/reference_data/ncbi_vert/"
f_name_file<-"ncbi_vert_trimmed_ncbi_names"

###   run for all non-silva  !!!
taxseq<- read.fasta(file = paste0(dir,ref_file,f_name_file,".fasta") , as.string = TRUE, strip.desc=TRUE)  # # psych::describe(getLength(taxseq))
taxonz_file<-"eDNA_db/reference_data/NCBI_tax_taxonomizr" 


############################################################################
############################################################################
############################################################################
####### now start running for all primers
##

setwd(paste(dir,taxonz_file, sep="") )
fast1<-getAnnot(taxseq)##  get only anotation data  
head(fast1)  # head(mess)
#fastg<-gsub(" .*$","",fast1)  ## keep only before first " "   head(fastg)  #  
fastg<-gsub("\\..*$","",fast1)  ## keep only before first "."  ### use for silva and ncbi
seq1 <- getSequence(taxseq) 

# get taxa id and taxonomy using taxonomizr
taxon_id<-accessionToTaxa(fastg,"accessionTaxa.sql",version='base') # head(taxon_id)
tax2<-data.frame(cbind(fastg,taxon_id))
names(tax2)[1]<-"accession"
names(tax2)[2]<-"taxID"    # head(df)
tax2<-dplyr::mutate_all(tax2,as.character)
tax<-getTaxonomy(tax2$taxID, "accessionTaxa.sql", desiredTaxa = c("superkingdom","kingdom",
        "phylum", "class", "order", "family", "genus", "species"), mc.cores = 1,debug = FALSE)
tax<-data.frame(tax,stringsAsFactors = F)


######################################################################################################
############  make data base of  lineages  and filter incomplete entries and clean up taxa names
#Set un-wanted terms
unwanteds <- c("fungal","endophyte","symbiont","uncultured", 
                "unclassified", "unidentified", "metagenome",
               "predicted", "environmental", "unknown", "unassigned")#
unwanteds <- paste0(unwanteds, collapse = "|")

tax2.2<- tax2 %>%     # head(tax2.2)    names(tax3)
  mutate(num=as.numeric(row.names(tax2))) %>%  
  bind_cols(tax) %>%   
  filter(complete.cases(accession)) %>% 
  filter(complete.cases(taxID)) %>% 
  mutate(taxID=as.numeric(taxID)) %>% 
  filter(taxID>0) %>%  # 0 is unknown - need to remove
  tidyr::unite(lineage, superkingdom:species,sep=";",remove=F) %>%
  filter(!grepl(unwanteds, lineage, fixed = FALSE, ignore.case = TRUE)) %>%  # remove unwanted terms
  filter(!duplicated(accession)) %>%   # remove duplicated accession numbers - silva should not have these, what -  must after decimal ?????????????
  mutate(species_old=species) %>% 
  mutate(species=gsub("(sp. ).*", "\\1", species)  ) %>%  # clean names remove everything after sp. in speceis 
  ## remove following symbols - note - perhaps change cf to sp., but most have sp.
  mutate(species=gsub("[", "", species, fixed = TRUE)  ) %>% 
  mutate(species=gsub("]", "", species, fixed = TRUE)  ) %>% 
  mutate(species=gsub(".cf ", "", species, fixed = TRUE)  ) %>% 
  mutate(species=gsub("cf. ", "", species, fixed = TRUE)  ) %>% 
  mutate(species=gsub("gen. ", "", species, fixed = TRUE)  ) %>% 
  mutate(species=gsub(".aff ", "", species, fixed = TRUE)  ) %>% 
  mutate(species=gsub("aff. ", "", species, fixed = TRUE)  ) %>% 
  mutate(species=gsub("n. ", "", species, fixed = TRUE)  ) %>% 
  mutate(species=gsub("nr. ", "", species, fixed = TRUE)  ) %>% 
  ## remove everything after and including 2nd space
  mutate(species = stringr::word(species, 1 ,2)) %>% 
  # remove white space, just in case
  mutate(species=trimws(species)) %>% 
  mutate(class=trimws(class))
  

  #################################################################################################################
#################################################################################################################

  ##  get worms names and lineage !!!!!!!!!!!!!!!!!
  ### test
#    xx<-tax2.2$species[1:5] # set species name column
#    sp_info <- worrms::wm_records_names(name = xx, fuzzy = F, marine_only = F) # get worrms info
#    worm_info <- bind_rows(sp_info)    #  make dataframe,   names(my_sp_taxo)   xx[10800:10900]


#################################################################################################################
#################################################################################################################
#     load if previous or run lines 288- 336
worm_info2 <- read_csv(paste0(path,"_",primer_name,"_worrms_data.csv"))
#################################################################################################################
#################################################################################################################


#####   begin loop to get basic worrms info -- if didn't load above
# set pattern
xx<-unique(tax2.2$species) # set species name column
z<-100  ## batch size
first<-1    

stt<-seq(first,length(xx),by=z)  # head(stt)   
stt2<-seq(first,length(xx),by=z*100)  # for taking nap, or may get kicked off worrms
stt2<-stt2[-1] # remove first

for (i in stt) {
  
  tryCatch({
  
  if (i== last(stt)){
    j<-length(xx) } else{j<-i+(z-1)}
  
  g<- xx[i:j]
  gg <- worrms::wm_records_names(name = g, fuzzy = F, marine_only = F) # get worrms info
  ggg<-bind_rows(gg) 
  
  if (i==1){
    worm_info<-ggg}
  else{
    worm_info<-rbind(worm_info,ggg)}
  print(j)  
  print(Sys.time())
  # make rest
  if (i %in% stt2){
    print("taking a quick namp")
    Sys.sleep(300)
  }
  
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}


## tidy
#   names(worm_info)     names(worm_info2)
worm_info2 <- worm_info %>% 
  select(AphiaID, scientificname, valid_AphiaID, valid_name, kingdom, phylum, class, order, family, genus, isMarine) %>% 
  rename(AphiaID_old = AphiaID, ncbi_species = scientificname, AphiaID = valid_AphiaID,  worrms_species=valid_name ) %>% 
  rename_at( vars(kingdom:genus), function(x) paste0("worrms_", x)  ) %>% 
  distinct(ncbi_species, .keep_all = TRUE)

#####
# export data if need to run again
# and  save data collected
write_csv(worm_info2, paste0(path,"_",primer_name,"_worrms_data.csv"))


############################################################################3
#   begin for all data - if loaded or run above code   !!!!!!!
 #################################################################################################################
#################################################################################################################

##   add wormms data, make best lineage - use worms names then fill in with ncbi
tax2.5<- tax2.2 %>%   # names(tax2.2)   names(tax2.5)
  rename_at( vars(kingdom:species), function(x) paste0("ncbi_", x)  )  %>% 
  left_join(worm_info2) %>% 
  # mutate(kingdom=worrms_kingdom,phylum=worrms_phylum, class=worrms_class,order=worrms_order,family=worrms_family,genus=worrms_genus,species=worrms_species) %>% 
  mutate(kingdom=ifelse(complete.cases(worrms_kingdom), worrms_kingdom, ncbi_kingdom) )   %>% 
  mutate(phylum=ifelse(complete.cases(worrms_phylum), worrms_phylum, ncbi_phylum) )   %>% 
  mutate(class=ifelse(complete.cases(worrms_class), worrms_class, ncbi_class) )   %>% 
  mutate(order=ifelse(complete.cases(worrms_order), worrms_order, ncbi_order) )   %>% 
  mutate(family=ifelse(complete.cases(worrms_family), worrms_family, ncbi_family ))   %>% 
  mutate(genus=ifelse(complete.cases(worrms_genus), worrms_genus, ncbi_genus) )   %>% 
  mutate(species=ifelse(complete.cases(worrms_species), worrms_species, ncbi_species) )  
  
 
####  correction lineage  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ncbi uses metazoa but change to worrms anamalia
# names(tax2.5)
tax2.6<- tax2.5 %>% 
  mutate( kingdom=if_else (kingdom=="Metazoa",  "Animalia" , kingdom )  ) %>% 
  
  mutate( class=if_else (class=="Dinoflagellata incertae sedis",  "Dinophyceae" , class )  ) %>% 
  mutate( order=if_else (order=="Gymnodiniaceae",  "Gymnodiniales " , order )  ) %>% 
  

  mutate( phylum=if_else (class=="Dinophyceae",  "Myzozoa" , phylum )  ) %>% 
  mutate( kingdom=if_else (class=="Dinophyceae",  "Chromista" , kingdom )  ) %>% 
  
  
  mutate( phylum=if_else (class=="Dictyochophyceae",  "Ochrophyta" , phylum )  ) %>% 
  mutate( kingdom=if_else (class=="Dictyochophyceae",  "Chromista" , kingdom )  ) %>% 

  mutate( phylum=if_else (class=="Ichthyosporea",  "Choanozoa" , phylum )  ) %>% 
  mutate( kingdom=if_else (class=="Ichthyosporea",  "Protozoa" , kingdom )  ) %>% 
  
  mutate( phylum=if_else (class=="Chrysophyceae",  "Ochrophyta" , phylum )  ) %>% 
  mutate( kingdom=if_else (class=="Chrysophyceae",  "Chromista" , kingdom )  ) %>%
  
  mutate( phylum=if_else (class=="Pelagophyceae",  "Ochrophyta" , phylum )  ) %>% 
  mutate( kingdom=if_else (class=="Pelagophyceae",  "Chromista" , kingdom )  ) %>% 
  
  mutate( phylum=if_else (class=="Phaeophyceae",  "Ochrophyta" , phylum )  ) %>% 
  mutate( kingdom=if_else (class=="Phaeophyceae",  "Chromista" , kingdom )  ) %>% 
  
  mutate( kingdom=if_else (phylum=="Rhodophyta",  "Plantae" , kingdom )  ) %>% # Rhodophyta
  mutate( kingdom=if_else (kingdom=="Viridiplantae",  "Plantae" , kingdom )  ) %>%
  
  mutate( class=if_else (order=="Thraustochytrida",  "Labyrinthulea" , class )  ) %>% 
  mutate( phylum=if_else (order=="Thraustochytrida",  "Bigyra" , phylum )  ) %>% 
  mutate( kingdom=if_else (order=="Thraustochytrida",  "Chromista" , kingdom )  ) %>% 
  mutate( class=if_else (order=="Labyrinthulida",  "Labyrinthulea" , class )  ) %>% 
  mutate( phylum=if_else (order=="Labyrinthulida",  "Bigyra" , phylum )  ) %>% 
  mutate( kingdom=if_else (order=="Labyrinthulida",  "Chromista" , kingdom )  ) %>% 
  mutate( class=if_else (order=="Bicosoecida",  "" , class )  ) %>% 
  mutate( phylum=if_else (order=="Bicosoecida",  "Bigyra" , phylum )  ) %>% 
  mutate( kingdom=if_else (order=="Bicosoecida",  "Chromista" , kingdom )  ) %>% 

mutate( class=if_else (order=="Chaetonotida",  "" , class )  ) %>% 
  mutate( phylum=if_else (order=="Chaetonotida",  "Gastrotricha" , phylum )  ) %>% 
  mutate( kingdom=if_else (order=="Chaetonotida",  "Animalia" , kingdom )  ) %>% 

mutate( kingdom=if_else (phylum=="Apicomplexa",  "Chromista" , kingdom )  ) %>% 

  mutate( class=if_else (class=="Coscinodiscophyceae",  "Bacillariophyceae" , class )  ) %>% 
mutate( phylum=if_else (phylum=="Bacillariophyta",  "Bacillariophyceae" , phylum )  ) %>% 
  mutate( phylum=if_else (phylum=="Bacillariophyceae",  "Ochrophyta" , phylum )  ) %>% 
  mutate( kingdom=if_else (phylum=="Ochrophyta",  "Chromista" , kingdom )  ) %>%

mutate( class=if_else (order=="Phaeocystales",  "Prymnesiophyceae" , class )  ) %>% 
  mutate( phylum=if_else (class=="Prymnesiophyceae",  "Haptophyta" , phylum )  ) %>% 
  
  mutate( kingdom=if_else (phylum=="Ciliophora",  "Chromista" , kingdom )  ) %>% 
  
  mutate( phylum=if_else (class=="Oomycota",  "Oomycota" , phylum )  ) %>% 
  mutate( kingdom=if_else (class=="Oomycota",  "Chromista" , kingdom )  ) %>% 
  
  mutate( phylum=if_else (class=="Bigyra",  "Bigyra" , phylum )  ) %>% 
  mutate( kingdom=if_else (class=="Bigyra",  "Chromista" , kingdom )  ) %>% 
  
mutate( kingdom=if_else (phylum=="Haptista",  "Chromista" , kingdom )  ) %>% 
  
mutate( phylum=if_else (class=="Cryptophyceae",  "Cryptophyta" , phylum )  ) %>% 
  mutate( kingdom=if_else (class=="Cryptophyceae",  "Chromista" , kingdom )  ) %>% 

  mutate( kingdom=if_else (phylum=="Euglenozoa",  "Protozoa" , kingdom )  )  %>% 
mutate( kingdom=if_else (phylum=="Firmicutes",  "Bacteria" , kingdom )  )  %>% 
mutate( kingdom=if_else (phylum=="Proteobacteria",  "Bacteria" , kingdom )  ) %>% 

  mutate( kingdom=if_else (phylum=="Evosea",  "Protozoa" , kingdom )  ) %>% 
  mutate( kingdom=if_else (phylum=="Euglenozoa",  "Protozoa" , kingdom )  ) %>% 
  
mutate( phylum=if_else (order=="Macrodasyida",  "Gastrotricha" , phylum )  ) %>% 
  mutate( kingdom=if_else (order=="Macrodasyida",  "Animalia" , kingdom )  ) %>% 
  
  mutate( phylum=if_else (class=="Acantharea",  "Radiozoa" , phylum )  ) %>% 
  mutate( kingdom=if_else (class=="Acantharea",  "Chromista" , kingdom )  ) %>% 
  mutate( phylum=if_else (class=="Acantharai",  "Radiozoa" , phylum )  ) %>% 
  mutate( kingdom=if_else (class=="Acantharia",  "Chromista" , kingdom )  ) %>% 
  
mutate( class=if_else (class=="Choanoflagellata",  "Choanoflagellatea" , class )  ) %>% 
mutate( phylum=if_else (class=="Choanoflagellatea",  "Choanozoa" , phylum )  ) %>% 
  mutate( kingdom=if_else (class=="Choanoflagellatea",  "Protozoa" , kingdom )  ) %>% 
  
mutate( class=if_else (order=="Thaumatomonadida",  "Imbricatea" , class )  ) %>% 
  mutate( phylum=if_else (class=="Imbricatea",  "Cercozoa" , phylum )  ) %>% 
  
  mutate( class=if_else (class=="Polycystinea",  "Polycystina" , class )  ) %>% 
  mutate( phylum=if_else (class=="Polycystina",  "Radiozoa" , phylum )  ) %>% 
  mutate( kingdom=if_else (class=="Polycystina",  "Chromista" , kingdom )  ) %>%  
  
  mutate( class=if_else (order=="Leptomyxida",  "Tubulinea" , class )  ) %>% 
  mutate( phylum=if_else (class=="Tubulinea",  "Amoebozoa" , phylum )  ) %>% 
  mutate( kingdom=if_else (class=="Tubulinea",  "Protozoa" , kingdom )  ) %>%  


mutate( kingdom=if_else (phylum=="Discosea",  "Protozoa" , kingdom )  )  %>% 
mutate( kingdom=if_else (phylum=="Tubulinea",  "Protozoa" , kingdom )  )  %>% 
mutate( kingdom=if_else (phylum=="Endomyxa",  "Chromista" , kingdom )  )  %>% 
mutate( kingdom=if_else (phylum=="Cercozoa",  "Chromista" , kingdom )  )  
  
mes<- tax2.6 %>% 
  filter(phylum=="Rhodophyta")

  ####    get data2 and insect lineage and clean up, remove taxa with limited taxonomic information
tax3<- tax2.6 %>% 
  filter_at(vars(family,species), any_vars(complete.cases(.)))    %>%  # no species ,genus, and family
  filter_at(vars(phylum:genus), any_vars(complete.cases(.)))  %>%   # no phylum through genus
  filter_at(vars(phylum:class), any_vars(complete.cases(.)))  %>%   # no phylum and class  !!!!!!!
  tidyr::unite(dada2_lineage, superkingdom,kingdom:species,sep=";",remove=F) %>%
  tidyr::unite(insect_lineage, accession,taxID, sep="|",remove=F) 



# sanity check
mes<-tax3 %>% 
  filter(grepl("labrax", dada2_lineage, fixed = FALSE, ignore.case = TRUE))

#  check taxanomy is in insect data and correct - only for insect
### check taxIDs against insect ncbi database  ########################
db<-insect::taxonomy(db = "NCBI", synonyms = F) # download ncbi taxonomy database directly (about 2 min) or can save and use below
## remove taxIDs not in db or will cause error for insect training   head(tax3$num)
tax4 <-tax3 %>%
    inner_join(db)
###   re-make orignal DNAbin to match filter sequences
#    head(x2)  names(x2)
## remove next 3 lines????
seq2<-seq1[tax4$num]
nnnn<-as.character(tax4$accession) # head(nnnn)  mes<-nnnn[duplicated(nnnn)]
names(seq2) <- c(nnnn)

## purge suspect sequences
#z<-tax3$taxID
#tdb <- prune_taxonomy(db, taxIDs = z)
#set.seed(999)
#zz <- purge(seq2, tdb, cores=16) # none to remove

########################################################################################
########################################################################################
########################################################################################
# export fasta
########################################################################################

########################################################################################
## for 18S primers
  setwd(path)     # getwd()
  fname_out_insect<-paste("SILVA_138_trimmed_",primer_name,"_insect_names.fasta",sep="")
  fname_out_dada<-paste("SILVA_138_trimmed_",primer_name,"_dada2_names.fasta",sep="")
write.fasta(sequences = seq1[tax4$num], names = tax4$insect_lineage, file.out = fname_out_insect, open = "w")  #
# and 
write.fasta(sequences = seq1[tax3$num], names = tax3$dada2_lineage, file.out = fname_out_dada, open = "w")  #
# and  save data collected
write_csv(tax3, paste("SILVA_138_trimmed_",primer_name,"_ncbi_and_worrms_data.csv"))

########################################################################################
############# for other primers

setwd(paste(dir,ref_file, sep="") )

fname_out_insect<-paste(f_name_file,"_",primer_name,"_insect_names.fasta",sep="")
fname_out_dada<-paste(f_name_file,"_",primer_name,"_dada2_names.fasta",sep="")
write.fasta(sequences = seq1[tax4$num], names = tax4$insect_lineage, file.out = fname_out_insect, open = "w")  #
# and 
write.fasta(sequences = seq1[tax3$num], names = tax4$dada2_lineage, file.out = fname_out_dada, open = "w")  #
# and  save data collected
write_csv(tax3, paste(f_name_file,"_",primer_name,"_ncbi_and_worrms_data.csv"))




##########################################################################################
##########################################################################################
##### 2.7- add additional barcodes !!
##########################################################################################
############################################################################

## euka02 dada2 - with aleja's and sarahs
dir<-"/Users/nathangeraldi/Dropbox/"
A <- insect::readFASTA(file =paste0(dir,"eDNA_db/reference_data/silva/SILVA_138_trimmed_euka02_dada2_names.fasta") )
B <- insect::readFASTA(file =paste0(dir,"eDNA_db/reference_data/Macrophyte_additions/Ortega_Barcoding_Sequences_dada2.fasta") )
B2 <- insect::readFASTA(file =paste0(dir,"eDNA_db/reference_data/Macrophyte_additions/macroalgae_sequence_arctic_Sarah_Oct2020.fasta") )

# mes<-getAnnot(B2)##

C <- insect::join(A,B)
C2 <- insect::join(A,B,B2)
insect::writeFASTA(C, file=paste0(dir,"eDNA_db/reference_data/silva/SILVA_138_trimmed_euka02_dada2_names_with_ortega_macroalgae.fasta") )
insect::writeFASTA(C2, file=paste0(dir,"eDNA_db/reference_data/silva/SILVA_138_trimmed_euka02_dada2_names_with_ortega_bachmann_macroalgae.fasta") )




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
fname_out_insect<-paste("SILVA_138_trimmed_",primer_name,"_insect_names.fasta",sep="")
x<-insect::readFASTA(file =fname_out_insect)
# only load if didn't above 
#    db<-insect::taxonomy(db = "NCBI", synonyms = F) # download ncbi taxonomy database directly (about 2 min) or can save and use below
# head(db)

insect_learn<- learn(x, db, model = NULL, refine = "Viterbi", iterations = 50,
      nstart = 20, minK = 2, maxK = 2, minscore = 0.9, probs = 0.5,
      retry = TRUE, resize = TRUE, maxsize = max(sapply(x, length)),
      recursive = TRUE, cores = "autodetect", quiet = F)

# Found 25130 (was 36315) unique sequences for Euka02 mare 29 2019 adn same nov, 2020
# Found 10107 (was 24347) unique sequences for 18smini MAr29, 2019
# Found 35653  (was 54959) unique sequences for 18s_stoeck MAr29, 2019


setwd(paste(dir,"eDNA_db/reference_data/insect_learn", sep="") )
saveRDS(insect_learn, paste(primer_name,"_learn.rds", sep=""))  ## then readRDS

#    setwd("C:/Users/geraldn/Dropbox/eDNA_db/reference_data/NCBI_tax_insect")
#    saveRDS(db, "db_insect.rds")  

