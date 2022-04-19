library(GenomicRanges)
library(tidyverse)
library(data.tree)
library(igraph)
#--------------------------
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#--------------------------
obj_in_fn<-function(file){
  obj_out<-get(load(file))
  tmp_obj<-names(mget(load(file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(obj_out)
}
#-------------------------
dagger_file<-"./data/pval_tbl/DAGGER/HMEC_poisson_DAGGER_01.Rda"
mres_dagger_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/DAGGER_tbl/HMEC_union_dagger_tbl.Rda"
emp_pval_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/pval_tbl/CAGE_union_HMEC_pval_tbl.Rda"
#-------------------------

dagger_tbl<-obj_in_fn(dagger_file)
mres_dagger_tbl<-obj_in_fn(mres_dagger_file)
emp_pval_tbl<-obj_in_fn(emp_pval_file)
dagger_tbl %>%
  dplyr::rename(cl=node) %>% 
  left_join(.,emp_pval_tbl) %>%
  mutate(res=fct_relevel(res,res_set)) %>% 
  ggplot(.,aes(-log10(pois.pval),-log10(emp.pval)))+
  geom_point()+
  facet_wrap(res~.,scales="free")+
  xlim(c(0,10))
  