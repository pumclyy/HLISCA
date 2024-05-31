##############################################################################
# Script information                                                      
# Title: MPRA
# Author: Erping Long
# Date: 2022-09-07
# Description: None
##############################################################################

library(tidyverse)
library(ggrepel)
library(viridis)
library(gplots)

#creating input data

#reading in table of RNA and DNA input files
transfections = read_delim("RNA_DNA_input_files_all.txt", delim = '\t', col_names = T, col_types = cols())

MPRA_data = data_frame()

##reading in design file
variantinfo <- read_csv("../Oligo_Library_w_ALL_hg19coords_Cell_Type_20220829.csv",
						            col_names = T, col_types = 'cccccccc')	

##reading data

for (i in 1:nrow(transfections)){
	rna_data <- read_delim(as.character(transfections[i,'RNA']), 
                                          col_names = T, 
                                          delim = '\t', 
                                          col_types = 'c-d-d') %>% 
                                          mutate(Tag = str_trunc(str_to_upper(Tag), 12, 'right', ellipsis = ""))
	rna_data <- variantinfo %>% left_join(rna_data, by = 'Tag')

	dna_data <- read_delim(as.character(transfections[i,'DNA']), 
                                          col_names = T, 
                                          delim = '\t', 
                                          col_types = 'c-d-d') %>% 
                                          mutate(Tag = str_trunc(str_to_upper(Tag), 12, 'right', ellipsis = ""))

	rna_data <- left_join(rna_data, dna_data, by = 'Tag', suffix = c("_RNA", "_DNA")) %>% mutate(TPM_RNA=TPM_RNA + 1) %>% mutate(TPM_DNA = TPM_DNA + 1) %>% filter(TPM_DNA > 1.6) %>% mutate(Ratio = TPM_RNA/TPM_DNA)  %>% mutate(Transfection = paste0('T', i))
	
	if(grepl('H520', transfections[i, 'RNA'] )){
		rna_data <- rna_data %>% mutate(Cell_Type = 'H520')
	} else if(grepl('A549', transfections[i, 'RNA'])){
		rna_data <- rna_data %>% mutate(Cell_Type = 'A549')
	}

            if(grepl('DMSO', transfections[i, 'RNA'] )){
		rna_data <- rna_data %>% mutate(Exposure = 'DMSO')
	} else if(grepl('BaP', transfections[i, 'RNA'])){
		rna_data <- rna_data %>% mutate(Exposure = 'BaP')
	}

	
	MPRA_data <- bind_rows(MPRA_data, rna_data)
}



Kaitest <- function(mydata=data,outlier.p.cut=0.05){
  suppressMessages(require(MASS))
  mydata$Ratio<-log2(mydata$Ratio+0.01)
  mydata.2 <- mydata
  
  
  if(test_type == 'Cell_type'){
	  model.1 <-lm(Ratio ~ Type + Strand + Cell_Type + Exposure, data=mydata.2)
	  model.2 <- lm(Ratio ~ Type + Strand + Cell_Type + Exposure + Type:Cell_Type, data=mydata.2)
  } else if(test_type == 'Exposure'){
	  model.1 <-lm(Ratio ~ Type + Strand + Cell_Type + Exposure, data=mydata.2)
	  model.2 <- lm(Ratio ~ Type + Strand + Cell_Type + Exposure + Type:Exposure, data=mydata.2)
  }
  
  ## Test the effect of Cell Type
  test.result.1<-anova(model.1, model.2, test="Chisq")
  aa<-dim(test.result.1)
  p.val <- exp(pchisq((test.result.1$`Sum of Sq`[2])*(test.result.1$Res.Df[2])/(test.result.1$RSS[2]),test.result.1$Df[2],lower.tail = FALSE,log.p = TRUE))
  return(p.val)
}

test.fun<-function(model.full, model.reduced, var.type="HC3")
{
  v.1<-vcovHC(model.full, type=var.type)
  v.2<-vcovHC(model.reduced, type=var.type)
  if (dim(v.1)[1] > dim(v.2)[1])
  {
    v.full<-v.1
    v.red<-v.2
    beta.full<-model.full$coefficient

  }
  if (dim(v.1)[1] < dim(v.2)[1])
  {
    v.full<-v.2
    v.red<-v.1
    beta.full<-model.reduced$coefficient  
    
  }
  if (dim(v.1)[1] == dim(v.2)[1])
  {
    stop('the two models are not nested')
  }
  par.full<-row.names(v.full)
  par.red<-row.names(v.red)
  
  par.test<-par.full[!par.full%in% par.red]
  test.df<-length(par.test)
  test.beta<-beta.full[par.test]
  test.v <-v.full[par.test, par.test]
  test.beta<-matrix(test.beta, byrow=T, ncol=1, nrow=test.df)
  test.v<-as.matrix(test.v)
  test.stat<-t(test.beta) %*% solve(test.v) %*% test.beta
  test.stat<-as.numeric(test.stat)
  pval<-pchisq(test.stat, df=test.df,lower.tail = FALSE)
 
  pval
 }

Kaitest_robust <- function(mydata=data,variant=snp, outlier.p.cut=0.05)
{
  suppressMessages(require(MASS))
  suppressMessages(require("sandwich"))
  my.outcome<-log2(mydata$Ratio+0.01)
  mydata<-data.frame(my.outcome, mydata)
  mydata.2 <- mydata
  
  ## test the cell type effect
  
  if(test_type == 'Cell_type'){
	  model.0 <-lm(Ratio ~ Type + Strand + Cell_Type + Exposure, data=mydata.2)
	  model.1 <- lm(Ratio ~ Type + Strand + Cell_Type + Exposure + Type:Cell_Type, data=mydata.2)
  } else if(test_type == 'Exposure'){
	  model.0 <-lm(Ratio ~ Type + Strand + Cell_Type + Exposure, data=mydata.2)
	  model.1 <- lm(Ratio ~ Type + Strand + Cell_Type + Exposure + Type:Exposure, data=mydata.2)
  }

  
  p.val.1 <- NA
  p.val.2 <- NA
  
  ## Test the effect of cell type using the standard wald test
  p.val.1<-test.fun(model.0, model.1, "const")
  
  ## Test the effect of cell Type using the robust test
  p.val.2<-test.fun(model.0, model.1, "HC3")
  
  return(c(p.val.1,p.val.2))
  
}


# testing cell type effect
result_kai_all1 <- data_frame(Variant=character(),pvalue_annovar=numeric(),pvalue_robust_const=numeric(),pvalue_robust_hc3=numeric())
test_type <- 'Cell_type'

tmpdata <- MPRA_data %>% filter(Type %in% c('ref', 'alt')) %>% dplyr::select(ID,Variant,REF,ALT,Type,Strand,Type_Strand,Transfection,Ratio,Cell_Type,Exposure) %>% unique()

for(snp in unique(tmpdata$Variant)){
  mydata <- tmpdata %>% filter(Variant==snp)
  
  tmpout1 <-NA
  tmpout2 <- c(NA,NA)
  
  try({ 
  tmpout1 <- Kaitest(mydata)
  tmpout2 <- Kaitest_robust(mydata, snp)
  })
  
 if (is.na(tmpout2[2])){print(snp)}

  result_kai_all1=bind_rows(
    result_kai_all1,
    data_frame(Variant=snp, pvalue_annovar=tmpout1,pvalue_robust_const=tmpout2[1],pvalue_robust_hc3=tmpout2[2])
  )}


result_kai_all1 <- result_kai_all1 %>%  mutate(FDR_annovar=p.adjust(pvalue_annovar,method = "BH"),FDR_robust_const=p.adjust(pvalue_robust_const,method = "BH"),FDR_robust_hc3=p.adjust(pvalue_robust_hc3,method = "BH")) 



## testing BaP effect

result_kai_all2 <- data_frame(Variant=character(),pvalue_annovar=numeric(),pvalue_robust_const=numeric(),pvalue_robust_hc3=numeric())
test_type <- 'Exposure'

tmpdata <- MPRA_data %>% filter(Type %in% c('ref', 'alt')) %>% dplyr::select(ID,Variant,REF,ALT,Type,Strand,Type_Strand,Transfection,Ratio,Cell_Type,Exposure) %>% unique()

for(snp in unique(tmpdata$Variant)){
  mydata <- tmpdata %>% filter(Variant==snp)
 
  tmpout1 <-NA
  tmpout2 <- c(NA,NA)
  
  try({ 
  tmpout1 <- Kaitest(mydata)
  tmpout2 <- Kaitest_robust(mydata, snp)
  })

  if (is.na(tmpout2[2])){print(snp)}

  result_kai_all2=bind_rows(
    result_kai_all2,
    data_frame(Variant=snp, pvalue_annovar=tmpout1,pvalue_robust_const=tmpout2[1],pvalue_robust_hc3=tmpout2[2])
  )}

result_kai_all2 <- result_kai_all2 %>%  mutate(FDR_annovar=p.adjust(pvalue_annovar,method = "BH"),FDR_robust_const=p.adjust(pvalue_robust_const,method = "BH"),FDR_robust_hc3=p.adjust(pvalue_robust_hc3,method = "BH"))


#######result output

write.csv(result_kai_all1,"./CellType_All.csv")
write.csv(result_kai_all2,"./BaP_All.csv")





