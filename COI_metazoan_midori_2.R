library(ape)
library('taxize')
library(ShortRead); packageVersion("ShortRead")
#library(ggplot2); packageVersion("ggplot2")
library(seqinr)
#library(dada2); packageVersion("dada2")
###set directory for specific computer Dropbox
library(dplyr)
library(tidyr)
library(taxonomizr)
#
## this code trims CO1 reference seqeunces from midori for with leray co1 primers. 
## then run code to clean names and cross reference with worrms - create_ref_lib_18s 
#########################  download data
##  go to http://reference-midori.info/download.php#
# in latest genbank release , rdp,  longest then download CO1 file and put in relevant folder

#  last updated Nov 2020 with genbank 239

########################## 1- virtualPCR   !!!!!!!!!!!!!!!!!!!!!!!!!!!!

###  get reference library
dir<-"/Users/nathangeraldi/Dropbox/"
path1<-paste(dir,"eDNA_db/reference_data/Midori_meta_COI", sep="")
setwd(path1) 
#list.files()
x<- readFASTA(file ="MIDORI_LONGEST_GB239_CO1_RDP 2.fasta")  # # psych::describe(getLength(taxseq))
####    minirbcl    head(x)  summary(x)
PRIMER1="GGWACWGGWTGAACWGTWTAYCCYCC"
PRIMER2="TANACYTCNGGRTGNCCRAARAAYCA"
maxlen<-450
minlen<-150

################################################################################################
################################################################################################
# run in silica and trim primers
trimmedf<-virtualPCR(x, up = PRIMER1, down = PRIMER2, rcdown = TRUE, trimprimers = TRUE,
                     minfsc = 50, minrsc = 50, minamplen = minlen, maxamplen = maxlen,
                     maxNs = 0.02, partialbind = TRUE, cores = "autodetect", quiet = FALSE)

## export to fasta
# set working directory
setwd(path1)   #  getwd()

## export to fasta
writeFASTA(trimmedf, file="MIDORI_LONGEST_GB239_CO1_trimmed.fasta")

################################################################################################
## as of Nov 2020, go to create_ref_lib_18s to get correct format for dada2 and insect
##
## the rest is old and can be removed

################################################################################################
################################################################################################
################################################################################################
################################################################################################





## the following is old - do not run
################################################################################################
#######################  2 - fix/get taxonomy
#  get reference fasta and fix names
setwd(paste(dir,"eDNA_db/reference_data/Midori_meta_COI", sep="") )
   taxseq<- read.fasta(file ="MIDORI_LONGEST_GB239_CO1_trimmed.fasta", as.string = TRUE, strip.desc=TRUE)  # # psych::describe(getLength(taxseq))
   fast1<-getAnnot(taxseq)##get only anotation data  head(fast1)
   seq1 <- getSequence(taxseq)
   fastg<-gsub("\\..*$","",fast1)  ## keep only before first "."  ### head(fastg)
   setwd(paste(dir,"eDNA_db/reference_data/NCBI_tax_taxonomizr", sep="") )
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
   

   
   
   # export fasta

   fname_out_insect<-paste("SILVA_132_trimmed_",primer_name,"_insect_names.fasta",sep="")
   fname_out_dada<-paste("SILVA_132_trimmed_",primer_name,"_dada2_names.fasta",sep="")
   write.fasta(sequences = seq1[tax3$num], names = tax3$insect_lineage, file.out = fname_out_insect, open = "w")  #
   # and 
   
   setwd("/Users/geraldn/Dropbox/eDNA_db/reference_data/Midori_meta_COI")
   write.fasta(sequences = seq1[tax3$num], names = tax3$dada2_lineage, file.out = "MIDORI_LONGEST_20180221_COI_trimmed_dada.fasta", open = "w")  #








