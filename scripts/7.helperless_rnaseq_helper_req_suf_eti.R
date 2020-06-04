library(ohchibi)



setwd('/home/isai/Documents/results/helperless//scripts/')
set.seed(130816)



#Try complext heatmap implementation
#Solution here  to add number of interesections
#https://support.bioconductor.org/p/118557/



#Red ojbect
mlist <- readRDS(file = "../cleandata/list_contrasts.RDS")

########## 4 hours ###############

### Analysis for RPS4

#Up#
r <- intersect(mlist$Up[["col_s4_4h_up"]],mlist$Up[["nrg1_s4_4h_up"]])

adr1_sufficient_4h_s4_up <- r[ ! ( r %in% mlist$Up[["hpk_s4_4h_up"]] ) ]


adr1_requiring_4h_s4_up <- adr1_sufficient_4h_s4_up[!(adr1_sufficient_4h_s4_up %in% mlist$Up[["adr1_s4_4h_up"]])] 


nrg1_adr1_redundant_4h_s4 <- adr1_sufficient_4h_s4_up[(adr1_sufficient_4h_s4_up %in% mlist$Up[["adr1_s4_4h_up"]])] 

#The contrary

s <- intersect(mlist$Up[["col_s4_4h_up"]],mlist$Up[["adr1_s4_4h_up"]]) 

nrg1_sufficient_4h_s4_up <- s[!(s %in% mlist$Up[["hpk_s4_4h_up"]])] 

nrg1_requiring_4h_s4_up <- nrg1_sufficient_4h_s4_up[!(nrg1_sufficient_4h_s4_up %in% mlist$Up[["nrg1_s4_4h_up"]])] 

helper_dependent_4h_s4_up <- mlist$Up[["col_s4_4h_up"]][!(mlist$Up[["col_s4_4h_up"]] %in% mlist$Up[["hpk_s4_4h_up"]])]

col_eti_4h_s4_up  =  mlist$Up[["col_s4_4h_up"]]


#Down#
### Analysis for RPS4
r <- intersect(mlist$Down[["col_s4_4h_down"]],mlist$Down[["nrg1_s4_4h_down"]])

adr1_sufficient_4h_s4_down <- r[ ! ( r %in% mlist$Down[["hpk_s4_4h_down"]] ) ]


adr1_requiring_4h_s4_down <- adr1_sufficient_4h_s4_down[!(adr1_sufficient_4h_s4_down %in% mlist$Down[["adr1_s4_4h_down"]])] 


nrg1_adr1_redundant_4h_s4 <- adr1_sufficient_4h_s4_down[(adr1_sufficient_4h_s4_down %in% mlist$Down[["adr1_s4_4h_down"]])] 

#The contrary

s <- intersect(mlist$Down[["col_s4_4h_down"]],mlist$Down[["adr1_s4_4h_down"]]) 

nrg1_sufficient_4h_s4_down <- s[!(s %in% mlist$Down[["hpk_s4_4h_down"]])] 

nrg1_requiring_4h_s4_down <- nrg1_sufficient_4h_s4_down[!(nrg1_sufficient_4h_s4_down %in% mlist$Down[["nrg1_s4_4h_down"]])] 

helper_dependent_4h_s4_down <- mlist$Down[["col_s4_4h_down"]][!(mlist$Down[["col_s4_4h_down"]] %in% mlist$Down[["hpk_s4_4h_down"]])]

col_eti_4h_s4_down  =  mlist$Down[["col_s4_4h_down"]]

### Analysis for RPt2

#Up#
r <- intersect(mlist$Up[["col_t2_4h_up"]],mlist$Up[["nrg1_t2_4h_up"]])

adr1_sufficient_4h_t2_up <- r[ ! ( r %in% mlist$Up[["hpk_t2_4h_up"]] ) ]


adr1_requiring_4h_t2_up <- adr1_sufficient_4h_t2_up[!(adr1_sufficient_4h_t2_up %in% mlist$Up[["adr1_t2_4h_up"]])] 


nrg1_adr1_redundant_4h_t2 <- adr1_sufficient_4h_t2_up[(adr1_sufficient_4h_t2_up %in% mlist$Up[["adr1_t2_4h_up"]])] 

#The contrary

s <- intersect(mlist$Up[["col_t2_4h_up"]],mlist$Up[["adr1_t2_4h_up"]]) 

nrg1_sufficient_4h_t2_up <- s[!(s %in% mlist$Up[["hpk_t2_4h_up"]])] 

nrg1_requiring_4h_t2_up <- nrg1_sufficient_4h_t2_up[!(nrg1_sufficient_4h_t2_up %in% mlist$Up[["nrg1_t2_4h_up"]])] 

helper_dependent_4h_t2_up <- mlist$Up[["col_t2_4h_up"]][!(mlist$Up[["col_t2_4h_up"]] %in% mlist$Up[["hpk_t2_4h_up"]])]

col_eti_4h_t2_up  =  mlist$Up[["col_t2_4h_up"]]

#Down#
### Analysis for RPt2
r <- intersect(mlist$Down[["col_t2_4h_down"]],mlist$Down[["nrg1_t2_4h_down"]])

adr1_sufficient_4h_t2_down <- r[ ! ( r %in% mlist$Down[["hpk_t2_4h_down"]] ) ]


adr1_requiring_4h_t2_down <- adr1_sufficient_4h_t2_down[!(adr1_sufficient_4h_t2_down %in% mlist$Down[["adr1_t2_4h_down"]])] 


nrg1_adr1_redundant_4h_t2 <- adr1_sufficient_4h_t2_down[(adr1_sufficient_4h_t2_down %in% mlist$Down[["adr1_t2_4h_down"]])] 

#The contrary

s <- intersect(mlist$Down[["col_t2_4h_down"]],mlist$Down[["adr1_t2_4h_down"]]) 

nrg1_sufficient_4h_t2_down <- s[!(s %in% mlist$Down[["hpk_t2_4h_down"]])] 

nrg1_requiring_4h_t2_down <- nrg1_sufficient_4h_t2_down[!(nrg1_sufficient_4h_t2_down %in% mlist$Down[["nrg1_t2_4h_down"]])] 

helper_dependent_4h_t2_down <- mlist$Down[["col_t2_4h_down"]][!(mlist$Down[["col_t2_4h_down"]] %in% mlist$Down[["hpk_t2_4h_down"]])]

col_eti_4h_t2_down  =  mlist$Down[["col_t2_4h_down"]]


### Analysis for RPm1

#Up#
r <- intersect(mlist$Up[["col_m1_4h_up"]],mlist$Up[["nrg1_m1_4h_up"]])

adr1_sufficient_4h_m1_up <- r[ ! ( r %in% mlist$Up[["hpk_m1_4h_up"]] ) ]


adr1_requiring_4h_m1_up <- adr1_sufficient_4h_m1_up[!(adr1_sufficient_4h_m1_up %in% mlist$Up[["adr1_m1_4h_up"]])] 


nrg1_adr1_redundant_4h_m1 <- adr1_sufficient_4h_m1_up[(adr1_sufficient_4h_m1_up %in% mlist$Up[["adr1_m1_4h_up"]])] 

#The contrary

s <- intersect(mlist$Up[["col_m1_4h_up"]],mlist$Up[["adr1_m1_4h_up"]]) 

nrg1_sufficient_4h_m1_up <- s[!(s %in% mlist$Up[["hpk_m1_4h_up"]])] 

nrg1_requiring_4h_m1_up <- nrg1_sufficient_4h_m1_up[!(nrg1_sufficient_4h_m1_up %in% mlist$Up[["nrg1_m1_4h_up"]])] 

helper_dependent_4h_m1_up <- mlist$Up[["col_m1_4h_up"]][!(mlist$Up[["col_m1_4h_up"]] %in% mlist$Up[["hpk_m1_4h_up"]])]

col_eti_4h_m1_up  =  mlist$Up[["col_m1_4h_up"]]

#Down#
### Analysis for RPm1
r <- intersect(mlist$Down[["col_m1_4h_down"]],mlist$Down[["nrg1_m1_4h_down"]])

adr1_sufficient_4h_m1_down <- r[ ! ( r %in% mlist$Down[["hpk_m1_4h_down"]] ) ]


adr1_requiring_4h_m1_down <- adr1_sufficient_4h_m1_down[!(adr1_sufficient_4h_m1_down %in% mlist$Down[["adr1_m1_4h_down"]])] 


nrg1_adr1_redundant_4h_m1 <- adr1_sufficient_4h_m1_down[(adr1_sufficient_4h_m1_down %in% mlist$Down[["adr1_m1_4h_down"]])] 

#The contrary

s <- intersect(mlist$Down[["col_m1_4h_down"]],mlist$Down[["adr1_m1_4h_down"]]) 

nrg1_sufficient_4h_m1_down <- s[!(s %in% mlist$Down[["hpk_m1_4h_down"]])] 

nrg1_requiring_4h_m1_down <- nrg1_sufficient_4h_m1_down[!(nrg1_sufficient_4h_m1_down %in% mlist$Down[["nrg1_m1_4h_down"]])] 

helper_dependent_4h_m1_down <- mlist$Down[["col_m1_4h_down"]][!(mlist$Down[["col_m1_4h_down"]] %in% mlist$Down[["hpk_m1_4h_down"]])]

col_eti_4h_m1_down  =  mlist$Down[["col_m1_4h_down"]]


########## 8 hours ###############
### Analysis for RPS4

#Up#
r <- intersect(mlist$Up[["col_s4_8h_up"]],mlist$Up[["nrg1_s4_8h_up"]])

adr1_sufficient_8h_s4_up <- r[ ! ( r %in% mlist$Up[["hpk_s4_8h_up"]] ) ]


adr1_requiring_8h_s4_up <- adr1_sufficient_8h_s4_up[!(adr1_sufficient_8h_s4_up %in% mlist$Up[["adr1_s4_8h_up"]])] 


nrg1_adr1_redundant_8h_s4 <- adr1_sufficient_8h_s4_up[(adr1_sufficient_8h_s4_up %in% mlist$Up[["adr1_s4_8h_up"]])] 

#The contrary

s <- intersect(mlist$Up[["col_s4_8h_up"]],mlist$Up[["adr1_s4_8h_up"]]) 

nrg1_sufficient_8h_s4_up <- s[!(s %in% mlist$Up[["hpk_s4_8h_up"]])] 

nrg1_requiring_8h_s4_up <- nrg1_sufficient_8h_s4_up[!(nrg1_sufficient_8h_s4_up %in% mlist$Up[["nrg1_s4_8h_up"]])] 

helper_dependent_8h_s4_up <- mlist$Up[["col_s4_8h_up"]][!(mlist$Up[["col_s4_8h_up"]] %in% mlist$Up[["hpk_s4_8h_up"]])]

col_eti_8h_s4_up  =  mlist$Up[["col_s4_8h_up"]]

#Down#
### Analysis for RPS4
r <- intersect(mlist$Down[["col_s4_8h_down"]],mlist$Down[["nrg1_s4_8h_down"]])

adr1_sufficient_8h_s4_down <- r[ ! ( r %in% mlist$Down[["hpk_s4_8h_down"]] ) ]


adr1_requiring_8h_s4_down <- adr1_sufficient_8h_s4_down[!(adr1_sufficient_8h_s4_down %in% mlist$Down[["adr1_s4_8h_down"]])] 


nrg1_adr1_redundant_8h_s4 <- adr1_sufficient_8h_s4_down[(adr1_sufficient_8h_s4_down %in% mlist$Down[["adr1_s4_8h_down"]])] 

#The contrary

s <- intersect(mlist$Down[["col_s4_8h_down"]],mlist$Down[["adr1_s4_8h_down"]]) 

nrg1_sufficient_8h_s4_down <- s[!(s %in% mlist$Down[["hpk_s4_8h_down"]])] 

nrg1_requiring_8h_s4_down <- nrg1_sufficient_8h_s4_down[!(nrg1_sufficient_8h_s4_down %in% mlist$Down[["nrg1_s4_8h_down"]])] 

helper_dependent_8h_s4_down <- mlist$Down[["col_s4_8h_down"]][!(mlist$Down[["col_s4_8h_down"]] %in% mlist$Down[["hpk_s4_8h_down"]])]

col_eti_8h_s4_down  =  mlist$Down[["col_s4_8h_down"]]

### Analysis for RPt2

#Up#
r <- intersect(mlist$Up[["col_t2_8h_up"]],mlist$Up[["nrg1_t2_8h_up"]])

adr1_sufficient_8h_t2_up <- r[ ! ( r %in% mlist$Up[["hpk_t2_8h_up"]] ) ]


adr1_requiring_8h_t2_up <- adr1_sufficient_8h_t2_up[!(adr1_sufficient_8h_t2_up %in% mlist$Up[["adr1_t2_8h_up"]])] 


nrg1_adr1_redundant_8h_t2 <- adr1_sufficient_8h_t2_up[(adr1_sufficient_8h_t2_up %in% mlist$Up[["adr1_t2_8h_up"]])] 

#The contrary

s <- intersect(mlist$Up[["col_t2_8h_up"]],mlist$Up[["adr1_t2_8h_up"]]) 

nrg1_sufficient_8h_t2_up <- s[!(s %in% mlist$Up[["hpk_t2_8h_up"]])] 

nrg1_requiring_8h_t2_up <- nrg1_sufficient_8h_t2_up[!(nrg1_sufficient_8h_t2_up %in% mlist$Up[["nrg1_t2_8h_up"]])] 

helper_dependent_8h_t2_up <- mlist$Up[["col_t2_8h_up"]][!(mlist$Up[["col_t2_8h_up"]] %in% mlist$Up[["hpk_t2_8h_up"]])]

col_eti_8h_t2_up  =  mlist$Up[["col_t2_8h_up"]]


#Down#
### Analysis for RPt2
r <- intersect(mlist$Down[["col_t2_8h_down"]],mlist$Down[["nrg1_t2_8h_down"]])

adr1_sufficient_8h_t2_down <- r[ ! ( r %in% mlist$Down[["hpk_t2_8h_down"]] ) ]


adr1_requiring_8h_t2_down <- adr1_sufficient_8h_t2_down[!(adr1_sufficient_8h_t2_down %in% mlist$Down[["adr1_t2_8h_down"]])] 


nrg1_adr1_redundant_8h_t2 <- adr1_sufficient_8h_t2_down[(adr1_sufficient_8h_t2_down %in% mlist$Down[["adr1_t2_8h_down"]])] 

#The contrary

s <- intersect(mlist$Down[["col_t2_8h_down"]],mlist$Down[["adr1_t2_8h_down"]]) 

nrg1_sufficient_8h_t2_down <- s[!(s %in% mlist$Down[["hpk_t2_8h_down"]])] 

nrg1_requiring_8h_t2_down <- nrg1_sufficient_8h_t2_down[!(nrg1_sufficient_8h_t2_down %in% mlist$Down[["nrg1_t2_8h_down"]])] 

helper_dependent_8h_t2_down <- mlist$Down[["col_t2_8h_down"]][!(mlist$Down[["col_t2_8h_down"]] %in% mlist$Down[["hpk_t2_8h_down"]])]

col_eti_8h_t2_down  =  mlist$Down[["col_t2_8h_down"]]


### Analysis for RPm1

#Up#
r <- intersect(mlist$Up[["col_m1_8h_up"]],mlist$Up[["nrg1_m1_8h_up"]])

adr1_sufficient_8h_m1_up <- r[ ! ( r %in% mlist$Up[["hpk_m1_8h_up"]] ) ]


adr1_requiring_8h_m1_up <- adr1_sufficient_8h_m1_up[!(adr1_sufficient_8h_m1_up %in% mlist$Up[["adr1_m1_8h_up"]])] 


nrg1_adr1_redundant_8h_m1 <- adr1_sufficient_8h_m1_up[(adr1_sufficient_8h_m1_up %in% mlist$Up[["adr1_m1_8h_up"]])] 

#The contrary

s <- intersect(mlist$Up[["col_m1_8h_up"]],mlist$Up[["adr1_m1_8h_up"]]) 

nrg1_sufficient_8h_m1_up <- s[!(s %in% mlist$Up[["hpk_m1_8h_up"]])] 

nrg1_requiring_8h_m1_up <- nrg1_sufficient_8h_m1_up[!(nrg1_sufficient_8h_m1_up %in% mlist$Up[["nrg1_m1_8h_up"]])] 

helper_dependent_8h_m1_up <- mlist$Up[["col_m1_8h_up"]][!(mlist$Up[["col_m1_8h_up"]] %in% mlist$Up[["hpk_m1_8h_up"]])]

col_eti_8h_m1_up  =  mlist$Up[["col_m1_8h_up"]]

#Down#
### Analysis for RPm1
r <- intersect(mlist$Down[["col_m1_8h_down"]],mlist$Down[["nrg1_m1_8h_down"]])

adr1_sufficient_8h_m1_down <- r[ ! ( r %in% mlist$Down[["hpk_m1_8h_down"]] ) ]


adr1_requiring_8h_m1_down <- adr1_sufficient_8h_m1_down[!(adr1_sufficient_8h_m1_down %in% mlist$Down[["adr1_m1_8h_down"]])] 


nrg1_adr1_redundant_8h_m1 <- adr1_sufficient_8h_m1_down[(adr1_sufficient_8h_m1_down %in% mlist$Down[["adr1_m1_8h_down"]])] 

#The contrary

s <- intersect(mlist$Down[["col_m1_8h_down"]],mlist$Down[["adr1_m1_8h_down"]]) 

nrg1_sufficient_8h_m1_down <- s[!(s %in% mlist$Down[["hpk_m1_8h_down"]])] 

nrg1_requiring_8h_m1_down <- nrg1_sufficient_8h_m1_down[!(nrg1_sufficient_8h_m1_down %in% mlist$Down[["nrg1_m1_8h_down"]])] 

helper_dependent_8h_m1_down <- mlist$Down[["col_m1_8h_down"]][!(mlist$Down[["col_m1_8h_down"]] %in% mlist$Down[["hpk_m1_8h_down"]])]

col_eti_8h_m1_down  =  mlist$Down[["col_m1_8h_down"]]


################ We have all the comparison lists ################
rm(r)
rm(s)
rm(mlist)


save.image(file = "../cleandata/helperless_contrasts_results.RData")
