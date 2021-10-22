library(ohchibi)


set.seed(130816)
######Set the work directory
setwd("/Users/isaisalasgonzalez/Documents/pi/scripts/")

######Open the csv files
Tab <- read.table(file = "../rawdata/Rawcounts_all.csv", header = T, sep = ",",
                  quote = "", comment.char = "")

Tax <- read.table(file = "../rawdata/taxonomy_all.csv", header = T, sep = ",",
                  quote = "", comment.char = "")

ASV <- read.table(file = "../rawdata/Tax_all.csv", header = T, sep = ",",
                  quote = "", comment.char = "")

Map <- read.csv(file = "../rawdata/Metadata_16S_phosphate.csv", header = T, sep =",")

######Add the ASV code in the taxonmy
mertax <- ASV[,c(6,7)]
Tax_merge <- cbind(mertax, Tax)

######Create the DADA_Id factor
Tab$DADA2_Header <- sapply(strsplit(Tab$Samples, "_"), `[`, 1)

###### Prepare the different object
row.names(Tab) <- Tab$DADA2_Header
Tab <- Tab[,1:28761]

###### Transform factors in metadata
colnames(Map)[c(2,7,9,11)] <- c("DADA_Header","Pots","Genotype","Fraction")
Map$Lines <- as.factor(Map$Lines)

Map$Genotype <- as.factor(Map$Genotype)
Map$Fraction <- Map$Fraction %>%
  gsub(pattern = "Leaf",replacement = "Shoot") %>%
  factor(levels = c("Soil","Root","Shoot"))
Map$Tray <- as.factor(Map$Tray)
Map$Group <- as.factor(Map$Group)
Map$DADA_Header <- paste0("R", Map$Rep, Map$DADA_Header)
Map$Rep <- paste0("Rep",Map$Rep) %>% factor
row.names(Map) <- Map$DADA_Header


#Create dataset for the phosphate transporter project
usable_samples <- intersect(Map$DADA_Header,rownames(Tab)) 
Map <- match(usable_samples,Map$DADA_Header) %>%
  Map[.,] %>% droplevels
Tab <- match(usable_samples,rownames(Tab)) %>%
  Tab[.,] %>% droplevels
#Remove Tab column with only 0
Tab <- Tab[,-which(colSums(Tab) ==0)]

###### Modify the structure of the ASV file to contain the proper qiime format needed by AMOR
Tax_merge <- Tax_merge %>%
  tidyr::unite(tax, Kingdom, Phylum, Class, Order, Family, Genus, sep = ";", remove = F)

Tax_merge <- Tax_merge[,1:3]

temp_new<-NULL
for(element in Tax_merge$tax){
  temp_vec<-unlist(strsplit(x = element,split = ";"))
  qiime_format<-paste("Root; k__",temp_vec[1],"; p__",temp_vec[2],"; c__",temp_vec[3],"; o__",temp_vec[4],"; f__",temp_vec[5],"; g__",temp_vec[6],sep="")
  temp_new<-c(temp_new,qiime_format)
}
Tax_merge$Taxonomy <- factor(temp_new)

Tax_merge <- Tax_merge[match(colnames(Tab),Tax_merge$ASV),]
row.names(Tax_merge) <- Tax_merge$ASV

#Create the dataset object for RawCounts
Dat <- create_dataset(Tab = Tab %>% t,Map = Map, Tax = Tax_merge)

#Create contaminants to remove later
contam_otus <- c(grep(pattern = "chloroplast", as.character(Dat$Tax$Taxonomy),ignore.case = TRUE),
                 grep(pattern = "mitochondri", as.character(Dat$Tax$Taxonomy),ignore.case = TRUE),
                 grep(pattern = "oomycete", as.character(Dat$Tax$Taxonomy),ignore.case = TRUE),
                 grep(pattern = "k__NA", as.character(Dat$Tax$Taxonomy),ignore.case = TRUE),
                 grep(pattern = "p__NA", as.character(Dat$Tax$Taxonomy),ignore.case = TRUE),
                 grep(pattern = "k__Eukaryota", as.character(Dat$Tax$Taxonomy),ignore.case = TRUE),
                 grep(pattern = "k__Archaea", as.character(Dat$Tax$Taxonomy),ignore.case = TRUE)
)
contam_otus <- row.names(Dat$Tax)[contam_otus]

# Filter Dataset
Dat_filter <- remove_taxons(Dat = Dat, taxons = contam_otus)
Dat_filter <- clean(Dat = Dat_filter,verbose = TRUE)
contam_otus <- Dat$Tax[ contam_otus, ]
Dat <- Dat_filter

#Compute the Depth per Sample
Dat$Map$Depth <- colSums(Dat$Tab)

#Compute the Log Depth per Sample
Dat$Map$LDepth <-log(colSums(Dat$Tab))

##Plot usable reads
Dat <- clean(Dat)

#Remove samples with less than 500 reads
Dat_raw <- remove_samples(Dat,row.names(Dat$Map)[ Dat$Map$Depth < 100])
Dat_raw <- clean(Dat_raw)

#Create taxonomy data.frame
mdf <- Dat_raw$Tax
df_tax <- mdf$Taxonomy %>% as.character %>%
  strsplit(split  = "\\;") %>%
  unlist %>%  gsub(pattern = "[a-z]__",replacement = "",) %>%
  gsub(pattern = " ",replacement = "") %>%
  matrix(data = .,ncol = 7,byrow = T) %>% as.data.frame
colnames(df_tax) <-c("Root","Kingdom","Phylum",
                     "Class","Order","Family","Genus")
df_tax <- cbind(mdf$ASV,df_tax)
colnames(df_tax)[1] <- "ASV"
df_tax$ASV <- df_tax$ASV %>% factor
df_tax$Root <- df_tax$Root %>% factor
df_tax$Kingdom <- df_tax$Kingdom %>% factor
df_tax$Phylum <- df_tax$Phylum %>% factor
df_tax$Class <- df_tax$Class %>% factor
df_tax$Order <- df_tax$Order %>% factor
df_tax$Family <- df_tax$Family %>% factor
df_tax$Genus <- df_tax$Genus %>% factor
row.names(df_tax) <- df_tax$ASV

#Rarefaction
Dat_rar <- rarefaction.Dataset(x = Dat_raw,sample =100)
Dat_rar <- clean(Dat_rar)

###Relative abundance 
Dat_ra <- Dat_rar
#Here we are scaling each column (sample) by the total number of reads in that sample
Dat_ra$Tab<-scale(x = Dat_ra$Tab,center = F,scale = colSums(Dat_ra$Tab))

##Save the raw , rarefied and relative abundance datasets into a global structure
Dat_amplicon <- list(RawCounts=Dat_raw,Rarefied=Dat_rar,RelativeAbundance=Dat_ra, df_tax = df_tax)
saveRDS(object = Dat_amplicon,file = "../cleandata/dat_asv_phosphate.RDS")

rm(list=ls())
gc()
