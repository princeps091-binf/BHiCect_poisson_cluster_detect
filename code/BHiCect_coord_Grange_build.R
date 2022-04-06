# Produce Coordinate tables for BHiCect clusters
library(tidyverse)
library(GenomicRanges)
library(furrr)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#--------------------------------------------------------------------
res_file<-"~/Documents/multires_bhicect/data/GM12878/spec_res/"
out_folder<-"./data/GRanges/BHiCect_Grange/GM12878/"

chr_set<-unlist(lapply(strsplit(grep("chr",list.files(res_file),value=T),split="_"),'[',1))
for (chromo in chr_set){
  print(chromo)
  load(paste0(res_file,chromo,"_spec_res.Rda"))
  chr_cl_tbl<-tibble(chr=chromo,cl=names(chr_spec_res$cl_member),bins=chr_spec_res$cl_member)
  chr_cl_tbl<-chr_cl_tbl %>% mutate(res=map_chr(cl,function(x){
    return(strsplit(x,split="_")[[1]][1])
  }))
  plan(multisession, workers = 3)
  
  chr_cl_tbl<-chr_cl_tbl %>% mutate(GRange=future_pmap(list(chr,bins,res),function(chr,bins,res){
    return(GRanges(seqnames=chr,
                   ranges = IRanges(start=as.numeric(bins),
                                    end=as.numeric(bins)+res_num[res]-1
                   )))
    
    
  }))
  
  save(chr_cl_tbl,file=paste0(out_folder,chromo,"_BHiCect_cl.Rda"))
  
}

