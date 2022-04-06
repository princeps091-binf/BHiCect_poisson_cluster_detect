library(tidyverse)
library(GenomicRanges)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set

#-----------------------------------------
## Utility function
build_GRange_fn<-function(CAGE_file){
  cage_coord_tbl<-get(base::load(CAGE_file))
  tmp_obj<-names(mget(base::load(CAGE_file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  cage_coord_tbl<-cage_coord_tbl %>% filter(!(is.na(start)))
  cage_chr_Grange<-GRanges(seqnames=cage_coord_tbl$chr,
                           ranges = IRanges(start=as.numeric(cage_coord_tbl$start),
                                            end=as.numeric(cage_coord_tbl$end)
                           ))
  return(cage_chr_Grange)
}

#-----------------------------------------

CAGE_file<-"./data/CAGE_union_coord_H1_tbl.Rda"
out_file<-"./data/GRanges/CAGE_union_H1_Grange.Rda"

cage_chr_Grange<-build_GRange_fn(CAGE_file)
save(cage_chr_Grange,file=out_file)
