library(ohchibi)
library(emmeans)
library(Rmisc)
library(palettesPM)
library(multcomp)
library(limma)

set.seed(130816)

######Define color for the genotypes
paleta_gen <- c("#A6CEE3", "#60FDF4", 
                "#B565F0","#1F74F0",
                "#E5C494", "#8CD17D",
                "#FDB462", "#F1CE63", 
                "#499894", "#86BCB6", "#E15759",
                "#E03558","#79706E", "#BAB0AC", "#D37295", 
                "#FABFD2", "#DAA520")
names(paleta_gen) <- c("Col-0", "WS", 
                       "bik1", "Si81psl1",
                       "phf1", "phr1",
                       "phr1/phl1", "pht1,1;1,4/pPHT1,4:PHT1,4 GFP # 2-1",
                       "pht1,1;1,4/pPHT1,4:PHT1,4 GFP # 3-1",
                       "pht1,1;1,4/pPHT1,4:PHT1,4;S258A - GFP #1-2",
                       "pht1,1;1,4/pPHT1,4:PHT1,4;S258A - GFP #2-3",
                       "pht1,1", "pht1,1;1,4", "pht1,4", "pht1,4/S1-L11_138643",
                       "pht1,4/SAIL_1225_G08", "Bulk soil")


palette_variance <- paletteer_d("dutchmasters::pearl_earring",11)[c(1:4,7)] %>%
  c(.,"white")
names(palette_variance) <- c("Fraction","Genotype","Phosphate",
                             "Fraction:Phosphate","Fraction:Genotype","Residual")



distfun <- function(x,method) vegan::vegdist(x = x, method = "bray")

Dat <- readRDS(file = "../cleandata/dat_asv_Pht1-BIK.RDS")

Dat_sub <- Dat$Rarefied 

#Permanova
Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~   Genotype + Rep ,
                      data = Dat_sub$Map,
                      #strata = Dat_sub$Map$Rep,
                      permutations = 9999)
mypermanova



mcap <- oh.cap(Tab = Dat_sub$Tab ,Map = Dat_sub$Map ,
               formula = "Genotype + Condition(Rep)",
               distfun = distfun,perms = 9999)


#Summary plot
Map_cap <- mcap$Map_cap
df_1 <- summarySE(data = Map_cap,measurevar = "CAP1",groupvars = "Genotype")
df_2 <- summarySE(data = Map_cap,measurevar = "CAP2",groupvars = "Genotype")
df_3 <- summarySE(data = Map_cap,measurevar = "CAP3",groupvars = "Genotype")
df_4 <- summarySE(data = Map_cap,measurevar = "CAP4",groupvars = "Genotype")


merged <- merge(df_1[,c("Genotype","CAP1","ci")],df_2[,c("Genotype","CAP2","ci")],by = "Genotype") 
colnames(merged)[c(3,5)] <- c("ci.CAP1","ci.CAP2")


xlab_id <- paste0("CAP1 (",mcap$variance_explained_axis[1],"%)")
ylab_id <- paste0("CAP2 (",mcap$variance_explained_axis[2],"%)")
annotation_id <- paste0("R2 = ",round(mypermanova$aov.tab$R2[1],2),"\n","p = ",
                        format.pval(mypermanova$aov.tab$`Pr(>F)`[1]))
merged$Background <- match(merged$Genotype,Dat_sub$Map$Genotype) %>%
  Dat_sub$Map$Background[.]

merged$Genotype <- merged$Genotype %>%
  gsub(pattern = "Si81",replacement = "bik1") %>%
  factor(levels = c("Col-0","bik1","WS","pht1,1;1,4"))
mcap$Map_cap$Genotype <- mcap$Map_cap$Genotype %>% 
  gsub(pattern = "Si81",replacement = "bik1") %>%
  factor(levels = c("Col-0","bik1","WS","pht1,1;1,4"))

p1 <- ggplot(data = merged,aes(CAP1,CAP2)) +
  geom_point(data = mcap$Map_cap,mapping = aes(CAP1,CAP2,color = Genotype),alpha = 0.5,size = 2) +
  geom_linerange(mapping = aes(xmin = CAP1 - ci.CAP1,xmax = CAP1 + ci.CAP1),alpha = 0.1) +
  geom_linerange(mapping = aes(ymin = CAP2 - ci.CAP2,ymax = CAP2 + ci.CAP2),alpha = 0.1) +
  geom_point(aes(color = Genotype),size = 6) +
  scale_x_continuous(limits = c(-1.3,2.2),oob = rescale_none) +
  scale_y_continuous(limits = c(-1.3,1.3),oob = rescale_none) +
  theme_ohchibi() +
  scale_shape_manual(values = c(16,15),na.value = 16) +
  ggtitle(label = "Root fraction") +
  scale_color_manual(values = paleta_gen) +
  xlab(label = xlab_id) +
  ylab(label = ylab_id) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
  ) +
  annotate(geom = "text",x = 2,y = 1,label = annotation_id)

#### Phylogram #####
Dat_rab <- Dat$RelativeAbundance 
#Remove batch effect
Tab <- Dat_rab$Tab
Tab_batch <- removeBatchEffect(x = Tab,Dat_rab$Map$Rep)

Dat_batch <- create_dataset(Tab = Tab_batch,Map = Dat_rab$Map,Tax = Dat_rab$Tax)
#Dat_batch <- Dat_rab

#Collapse taxonomy
Dat_phyla <- Dat_batch %>%
  collapse_by_taxonomy.Dataset(Dat = .,level = 3)

#Split Tab and Map to create another structure
Tab  <- Dat_phyla$Tab
Map <- Dat_phyla$Map

rownames(Tab) <- Tab %>% rownames %>%
  strsplit(split = "\\;") %>% unlist %>%
  grep(pattern = "p__",value = T) %>% 
  gsub(pattern = "p__",replacement = "") %>%
  gsub(pattern = " ",replacement = "")


rownames(Tab) <- rownames(Tab) %>%
  gsub(pattern = "Acidobacteriota",replacement = "Acidobacteria") %>%
  gsub(pattern = "Actinobacteriota",replacement = "Actinobacteria") %>%
  gsub(pattern = "Bacteroidota",replacement = "Bacteroidetes") %>%
  gsub(pattern = "Gemmatimonadota",replacement = "Gemmatimonadetes") %>%
  gsub(pattern = "Verrucomicrobiota",replacement = "Verrucomicrobia") 


mphyla <- palettesPM::pm.names.phyla()
paleta <- palettesPM::pm.colors.phyla()


Tab_pr <- Tab[which((rownames(Tab) %in% mphyla )),]
Other <- Tab[which(!(rownames(Tab) %in% mphyla )),] %>%
  colSums
Tab <- rbind(Tab_pr,Other)

#Create new phyla dataset
Dat_phyla <- create_dataset(Tab = Tab,Map = Map)

Dat_phyla$Map$FractionByBackground <- paste0(Dat_phyla$Map$Fraction," ",Dat_phyla$Map$Background) %>%
  factor(levels = c("Soil Bulk soil","Soil Col-0","Soil WS","Root Col-0","Root WS","Shoot Col-0","Shoot WS"))

res <- chibi.phylogram(Tab = Dat_phyla$Tab,Map = Dat_phyla$Map,
                       facet_formula = "FractionByBackground+Genotype",
                       size_ticks_x = 0,size_strip_text = 35,size_axis_text = 25,
                       legend_proportion_size = 4,size_legend_text = 30,
                       size_axis_title = 0,font_family = "Arial",size_title_text = 35)

mdata <- res$p_raw$data


mdata$Genotype <- mdata$Genotype %>%
  gsub(pattern = "Si81",replacement = "bik1") %>%
  factor(levels = c("Col-0","bik1","WS","pht1,1;1,4"))
mdata$FractionByBackground <- mdata$FractionByBackground %>%
  gsub(pattern = "Root Col-0",replacement = "Col-0 background") %>%
  gsub(pattern = "Root WS",replacement = "WS background")

p2 <- ggplot(data = mdata,aes(Genotype,Abundance,fill = Taxon)) +
  geom_bar(stat = "identity", position = "fill", width = 1) +
  palettesPM::scale_fill_phyla(name = "Phylum") +
  facet_grid(.~FractionByBackground,space = "free",scales = "free") +
  scale_y_continuous(breaks = seq(0,1,0.1),limits = c(0,1),expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme_ohchibi() +
  xlab(label = "Genotype") +
  ylab(label = "Relative abundance (%)") +
  theme(
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1)
  )



#### Model

glm_model <- function(Dat_model = Dat, ref = "Col-0"){
  #Structure to hold the results
  Dat_model$Map$Genotype <- Dat_model$Map$Genotype %>% relevel(ref = ref)
  Res_pv <- NULL
  #Loop over each ESV
  for(i in 1:nrow(Dat_model$Tab)){
    cat("Working on",i,"\n")
    esv <- rownames(Dat_model$Tab)[i]
    y <- Dat_model$Tab[i,]
    Map <- Dat_model$Map
    Map$y <- y
    Map$LDepth <- log(Map$Depth)
    
    #Main effects model###
    #Add the error catch function
    mod_nbinom <- tryCatch(MASS::glm.nb(formula = y~ Genotype  +  offset(LDepth) + Rep,
                                        data = Map),error = function(e) list(call = NA))
    mod_poisson<- tryCatch(glm(formula = y~ Genotype   +  offset(LDepth) + Rep,
                               data = Map,family = poisson(link = "log")),error = function(e) list(call = NA))
    mod_poisson_zinf <- tryCatch(zeroinfl(formula = y~ Genotype  +  offset(LDepth)  | 1 + Rep,
                                          data = Map,dist = "poisson",link = "log"),error = function(e) list(call = NA))
    
    mod_nbinom_zinf <- tryCatch(zeroinfl(formula = y~ Genotype +  offset(LDepth) | 1 + Rep,
                                         data = Map,dist = "negbin",link = "log"),error = function(e) list(call = NA))
    #Create dataframe of AIC
    models <- c("mod_nbinom","mod_poisson","mod_poisson_zinf","mod_nbinom_zinf")
    s_models <- NULL
    for(model in models){
      num <- all(!((is.na(as.character(get(model)$call)))))
      if(num == TRUE){
        temp <- data.frame(model=model,AIC=AIC(get(model)))
        s_models <- rbind(s_models,temp)
      }
    }
    if(is.null(s_models)){
      cat("Did not enter\n")
    }else{
      s_models <- s_models[with(s_models,order(AIC)),]
      
      #Best model
      m1 <- get(as.character(s_models[1,1]))
      #Num genotype
      m1_em <- emmeans(m1,specs = "Genotype",adjust="None")
      df_pvals_genotype <- pairs(m1_em,reverse = T,adjust = "None") %>% as.data.frame
      #Subset only the ones that contain the reference value
      df_pvals_genotype <- df_pvals_genotype$contrast %>% grep(pattern =ref ,value = F) %>%
        df_pvals_genotype[.,] %>% droplevels
      df_pvals_genotype$Factor <- rep("Genotype",nrow(df_pvals_genotype))
      
      #Attach
      #Pvals
      res <- rbind(df_pvals_genotype)
      res$Id <- rep(esv,nrow(res))
      res$Model <- rep(s_models[1,1],nrow(res))
      res$contrast <- res$contrast %>% gsub(pattern = " ",replacement = "|")
      res$Genotype <- res$contrast %>% gsub(pattern = "\\|.*",replacement = "")
      Res_pv <- rbind(Res_pv,res)
      
    }
    
  }
  return(Res_pv)
}

Dat_sub <- Dat$RawCounts  %>%
  subset.Dataset(Genotype == "Col-0"  | Genotype == "Si81" ) %>%
  subset.Dataset(Fraction == "Root",drop = T,clean = T) 



total <- Dat_sub$Tab %>% sum
Dat_filter <- measurable_taxa.Dataset(Dat = Dat_sub,
                                      min_samples_otu =2,min_reads_otu =  )
Dat_filter <- clean(Dat_filter)
(Dat_filter$Tab %>% sum)/total

res_otu <- glm_model(Dat_model = Dat_filter,ref = "Col-0")
res_col <- res_otu



Dat_sub <- Dat$RawCounts  %>%
  subset.Dataset(Genotype == "WS"  | Genotype == "pht1,1;1,4" ) %>%
  subset.Dataset(Fraction == "Root",drop = T,clean = T) 



total <- Dat_sub$Tab %>% sum
Dat_filter <- measurable_taxa.Dataset(Dat = Dat_sub,
                                      min_samples_otu =2,min_reads_otu =  )
Dat_filter <- clean(Dat_filter)
(Dat_filter$Tab %>% sum)/total

res_otu <- glm_model(Dat_model = Dat_filter,ref = "WS")
res_ws <- res_otu

res_col$padj <- p.adjust(p = res_col$p.value,method = "fdr")
res_ws$padj <- p.adjust(p = res_ws$p.value,method = "fdr")

df_glm <- rbind(res_col,res_ws)


#Create heatmap of enrichment
df_tax <- Dat$df_tax


#Define threshold for enrichment
threshold <- 0.05
hclust_method <- "ward.D"

#### ASV ####
df_temp <- df_glm %>%
  subset(padj < threshold) %>%
  droplevels 

mids <- df_temp$Id %>% unique

df_temp <- which(df_glm$Id %in% mids) %>%
  df_glm[.,]


Tab <- acast(data = df_temp,formula =Genotype~Id,
             value.var = "estimate",fill = 0 )
mclust_id <- hclust(d = as.dist(1-cor(Tab)),method = hclust_method)
order_ids <- mclust_id$order %>% mclust_id$labels[.]

df_temp$Phylum <- match(df_temp$Id,df_tax$ASV) %>%
  df_tax$Phylum[.] %>% as.character

df_temp$Genus <- match(df_temp$Id,df_tax$ASV) %>%
  df_tax$Genus[.] %>% as.character

df_temp$Nom <- paste0(df_temp$Genus," ",df_temp$Id)
order_nom <- with(df_temp,order(Id)) %>%
  df_temp$Nom[.] %>% unique

df_temp$Nom <- df_temp$Nom %>% factor(levels = order_nom)
df_temp$Nom <- df_temp$Nom %>% 
  gsub(pattern = "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium",replacement = "Rhiobium")

#Order according to hitmap
df_temp$Phylum<- df_temp$Phylum %>%
  gsub(pattern = "Acidobacteriota",replacement = "Acidobacteria") %>%
  gsub(pattern = "Actinobacteriota",replacement = "Actinobacteria") %>%
  gsub(pattern = "Bacteroidota",replacement = "Bacteroidetes") %>%
  gsub(pattern = "Gemmatimonadota",replacement = "Gemmatimonadetes") %>%
  gsub(pattern = "Verrucomicrobiota",replacement = "Verrucomicrobia") 

mphyla <- palettesPM::pm.names.phyla()

df_temp$Phylum[which(!(df_temp$Phylum %in% mphyla))] <- "Other"

df_temp$Phylum[df_temp$Phylum %>% grep(pattern = "Patescibacteria")] <- "Other"
df_temp$Phylum[df_temp$Phylum %>% grep(pattern = "Gemmatimonadetes")] <- "Other"
df_temp$Phylum[df_temp$Phylum %>% grep(pattern = "Cyanobacteria")] <- "Other"

df_temp$Phylum <- df_temp$Phylum %>%
  factor(levels = mphyla)

df_temp$Background <- "Col0 background"
df_temp$Background[df_temp$contrast %>% grep(pattern = "WS")] <- "WS backgroud"
df_temp$Background <- df_temp$Background %>% factor

df_temp$Significance <- NA
df_temp$Significance[which(df_temp$padj < 0.05)] <- "q < 0.05"

df_temp$Genotype <- df_temp$Genotype %>% gsub(pattern = "Si81",replacement = "bik1")
p3 <- ggplot(data = df_temp ,aes(Nom,Genotype)) +
  geom_raster(aes(fill = estimate)) +
  geom_tile(aes(color = Significance), fill = "#00000000", 
            size = 0.3, width = 0.98, height = 0.98) +
  facet_grid(Background~Phylum,space = "free",scales = "free") +
  scale_color_manual(values = "black", na.value = "#00000000")+
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-5,5),oob = squish,na.value = "grey89",name = "Estimate\nvs Col-0|WS") +
  theme_ohchibi(size_panel_border = 0.1)  +
  scale_y_discrete(expand = c(0,0)) +
  theme(
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size = 2),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.ticks.x  = element_line(size = 0.1),
    axis.ticks.y = element_line(size = 0.1),
  ) +
  xlab(label = "ASV Id")

compositiona <- egg::ggarrange(p1,p2,nrow =1,widths = c(1,0.2))
oh.save.pdf(p = compositiona,outname = "pi_sub_cap_phylogram.pdf",
            outdir = "../figures/",width = 10,height = 5)

oh.save.pdf(p = p3,outname = "pi_sub_model.pdf",
            outdir = "../figures/",width = 10,height = 5)

rm(list=ls())
dev.off()
gc()
