################################################################################
##########             CONSORT Flow Diagram                           ##########
################################################################################


# Load libraries
library('diagram')
library('ggplot2')
library('tidyr')
library('dplyr')

# Generate the consort diagram
consort_diag<-function(consort_tbl){
  tbl<-data.frame(CNT_TYPE=c("Initial",
                             "Has_at_least_2_SCr",
                             "Initial_GFR_below_15",
                             "KRT_within_48hr",
                             "Burn_patients",
                             "Pre_ESRD",
                             "Pre_KRT",
                             "No_Baseline", 
                             "nonAKI",
                             "Total",
                             "AKI1",
                             "AKI2",
                             "AKI3"),
                  label_txt=c("Inpatient visit with LOS >= 2\nand of age >= 18",
                              "Has at least 2 SCr record",
                              "Excluded: Initial eGFR below 15",
                              "Excluded: KRT with 48 hours since \nadmission",
                              "Excluded: Burn patients",
                              "Excluded: Pre-existance of \nkidney failure",
                              "Excluded: Pre-existance of \ndialysis and kidney transplantation",
                              "Excluded: No SCr baseline CKD patients",  #
                              "Excluded: Non-AKI encounters",  #
                              "Total eligible encounters",
                              "AKI1 at onset",
                              "AKI2 at onset",
                              "AKI3 at onset"),
                  stringsAsFactors=F) %>%
    left_join(consort_tbl, by="CNT_TYPE") %>%
    replace_na(list(ENC_CNT=0)) %>%
    mutate(cnt_ref=ifelse(CNT_TYPE %in% c("Initial","Has_at_least_2_SCr","Total"),ENC_CNT,NA)) %>%
    fill(cnt_ref,.direction="down") %>%
    mutate(cnt_ref=ifelse(CNT_TYPE=="Has_at_least_2_SCr",lag(cnt_ref,n=1L),cnt_ref)) %>%
    mutate(ENC_PROP=round(ENC_CNT/cnt_ref,4)) %>%
    mutate(label_val=paste0("(",ENC_CNT,",",ENC_PROP*100,"%)")) %>%
    mutate(label=paste(label_txt,"\n",label_val)) %>%
    mutate(node_id=c(2,5,7,9,10,12,13,15,16,20,22,23,24))
  
  #prepare canvas
  par(mfrow=c(1,1))
  par(mar=c(0,0,0,0))
  openplotmat()

  elpos<-coordinates(rep(3,8))
  fromto<-matrix(ncol=2,byrow=T,
                 c(2,5,
                   5,8,
                   8,7,
                   8,9,
                   8,11,
                   11,10,
                   11,12,
                   11,14,
                   14,13,
                   14,15,
                   14,17,
                   17,16,
                   17,20,
                   20,19,
                   20,21,
                   19,22,
                   20,23,
                   21,24
                 ))
  ##draw arrows
  arrpos <- matrix(ncol = 2, nrow = nrow(fromto))
  for (i in 1:nrow(fromto)){
    arrpos[i, ] <- straightarrow (to = elpos[fromto[i, 2], ],
                                  from = elpos[fromto[i, 1], ],
                                  lwd = 1, arr.pos = 0.6, arr.length = 0.3)
  }
  
  ##draw nodes
  for(i in 1:nrow(tbl)){
    textrect(elpos[tbl$node_id[i],],
             radx=0.15,
             rady=0.05,
             lab=tbl$label[i],
             font=4,
             cex=0.7,
             shadow.size=0)
  }
}



consort_tbl <- data.frame(
  CNT_TYPE = c("Initial", 
               "Has_at_least_2_SCr", 
               "Initial_GFR_below_15", 
               "KRT_within_48hr", 
               "Burn_patients", 
               "Pre_ESRD", 
               "Pre_KRT", 
               "No_Baseline", 
               "nonAKI",
               "Total",
               "AKI1",
               "AKI2",
               "AKI3"),
  ENC_CNT = c(991972, 
              761620, 
              40828, 
              13827, 
              954, 
              51430, 
              57333, 
              2069,
              509862, 
              130699, 
              101057, 
              22120,  
              7522   
  )
)

consort_diag(consort_tbl)