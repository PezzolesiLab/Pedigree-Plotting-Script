library(tidyverse)
library(dplyr)
library(kinship2)
library(knitr)
library(reticulate)

setwd('Y:/IRB_00096551/UPDB_Pedigrees_10_2021/FullPedigree')

# Load pedigree file
families <- read_csv('David_Marcus_ICDRenal_ALL_T_20211011.txt')

# Make pedigree object
pedList <- pedigree(famid = families$kindredID, id = families$egoUPDBID, dad = families$paUPDBID, mom = families$maUPDBID, sex = families$sex, status = families$deceased)
#write_csv(families, "~/Documents/Projects/eGFR/EDW_IHC_RFD_Pedigrees/20190717DavidMarcusM5Pedigree/mappingDiabetics/kinship2_ready_pedigree_dm.csv")

#TODO: It appears we are missing some parents in our pedigrees, may need to 
# correct for that...
# 
# reference drawing_full_pedigrees.Rmd from box for fixing the missing parent's

