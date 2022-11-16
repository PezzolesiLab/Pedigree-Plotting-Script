library(devtools)
library(GENLIB)
#install_url("https://cran.r-project.org/src/contrib/Archive/GENLIB/GENLIB_1.0.6.tar.gz")
#install.packages("D:/Users/u1123911/Downloads/GENLIB_1.0.8.tar.gz", repos = NULL, type="source")
library(statnet)
library(networkR)
library(dplyr)
library(gap)
library(kinship2)
library(data.tree)
library(treemap)
library(networkD3)
library(pedigree)
library(ggenealogy)
library(stringr)
library(stringi)
library(tidyverse)

setwd("/Users/bradyneeley/Work")

Pedigree_File_Original <- read_tsv("c5orf42_ped.tsv", col_types = "nnnncnnnccnnnncnccnn")

####################START DRAWING PEDIGREES

###### Multiple Pedigrees
Families_list <- unique(Pedigree_File_Original$kindredID)

# I think these next three lines go unused.
#updbids <- Pedigree_File$UPDB_ID[!is.na(Pedigree_File$UPDB_ID)]
#distids <- Pedigree_File$DIST_ID[!is.na(Pedigree_File$DIST_ID)]
#personids <- Pedigree_File$egoUPDBID[!is.na(Pedigree_File$egoUPDBID)]

affectionTable <- select(Pedigree_File_Original, kindredID, egoUPDBID, paUPDBID, maUPDBID, sex, affection, Slope, esrd, dmType)

affectionTable$sex <- ifelse(affectionTable$sex == "M", 1, 2)

#Leave out Scotts "dummy" data from kinship2
affectionTable$paUPDBID[is.na(affectionTable$paUPDBID)] <- 0
affectionTable$maUPDBID[is.na(affectionTable$maUPDBID)] <- 0
affectionTable$egoUPDBID[affectionTable$egoUPDBID < 0] <- 1
sum(is.na(Pedigree_File_Original$paUPDBID)) + sum(is.na(Pedigree_File_Original$maUPDBID))
affectionTable$paUPDBID[affectionTable$paUPDBID < 0] <- 1
sum(affectionTable$paUPDBID == 0) + sum(affectionTable$maUPDBID == 0)

#colnames(Families_list) <- c("ped_number","DIST_IDs","count")
#Families_list$DIST_IDs <- gsub(";", ",", Families_list$DIST_IDs)
#row.names(Families_list) <- Families_list$ped_number

for (ped_number in Families_list) {    
# myFunc <- function(pedNumber) {
  # The 2 lines below did not seem useful...
  #print(ped_number)
  #ped_number <- 26174
  
  ### Make an affection column for Rapid Decliners called affectionRD
  Pedigree_File <- affectionTable %>% filter(kindredID == ped_number)
  # if NA then 0, if ESRD then 1
  Pedigree_File$esrd <- ifelse(is.na(Pedigree_File$esrd), 0, 1)
  Pedigree_File$slope <- Pedigree_File$Slope
  RD <- Pedigree_File$Slope <= -5
  Pedigree_File$Slope <- RD
  Pedigree_File$Slope[is.na(Pedigree_File$Slope)] <- 0
  # Rapid decliners added as 1s to the rapid_decliner column
  Pedigree_File$rapid_decliner<- ifelse((Pedigree_File$Slope == TRUE | Pedigree_File$esrd == 1), '990000', '999999')
  Pedigree_File$rapid_decliner[is.na(Pedigree_File$rapid_decliner)] <- '0'
  
  ### Format the dmType column into a new column for Cranefoot called dmTypeNum
  Pedigree_File$dmTypeNum <- Pedigree_File$dmType
  Pedigree_File$dmTypeNum[Pedigree_File$dmType == "T1D"] <- "41"
  Pedigree_File$dmTypeNum[Pedigree_File$dmType == "T2D"] <- "43"
  Pedigree_File$dmTypeNum[Pedigree_File$dmType == "T3D"] <- "43"
  Pedigree_File$dmTypeNum[Pedigree_File$dmType == "IMPUTED_T1D"] <- "21"
  Pedigree_File$dmTypeNum[Pedigree_File$dmType == "IMPUTED_T2D"] <- "23"
  # TODO: Change the columns here that are NA vals
  
  
  # If 3+ rows, grabs cols 2-7 and 11, and name them...
  if (nrow(Pedigree_File) > 3) {
    #Pedigre_to_genlib <- Pedigree_File[ , c(2,3,4,5,6,7,11)]
    #colnames(Pedigre_to_genlib) <- c("ind","father","mother","sex","deceased","hypoglycemia","T1D")
    
    #  Delete any duplicate rows for an individual
    Pedigree_File <- distinct(Pedigree_File)
    #Pedigre_to_genlib <- distinct(Pedigre_to_genlib, ind, .keep_all = TRUE)
 
    Pedigre_to_genlib <- Pedigree_File   
    Pedigre_to_genlib <- group_by(Pedigre_to_genlib,egoUPDBID)
    
    #Not doing anything in 3135s case except make the table with these colnames...
    Pedigre_to_genlib_summary <- summarize(Pedigre_to_genlib, father = max(paUPDBID),
                                           mother = max(maUPDBID),sex = median(sex),
                                           Total = n())
    
    #validate_trio_consistency(Pedigre_to_genlib_summary$ind,Pedigre_to_genlib_summary$father,
    #                          Pedigre_to_genlib_summary$mother,Pedigre_to_genlib_summary$sex)
    
    Missing_parents <- as.data.frame(Pedigre_to_genlib_summary[ , 1:4])
    colnames(Missing_parents)[1] = "ind"
    Missing_parents$father[is.na(Missing_parents$father)] <- 0
    Missing_parents$mother[is.na(Missing_parents$mother)] <- 0
    Missing_parents$ind <- abs(Missing_parents$ind)
    
    # This generates an object for genlib and uses autocomplete to add father/mother ids that
    # aren't present in the individual id column to the individual column and marks them as not
    # having parents...
    Pedigree_to_genlib <- gen.genealogy(Missing_parents, autoComplete = TRUE)
    
    ###PROBANDS LIST
    #Probands = as.character(c(1049451,1086888,1259537,1337787,1505327,156127,2644237,2660265,2751457,2761151,315159,336820,47986,625687,796196,842366,915818,965663))
    
    # make this list the list of hypoglycemics (Jose)
    #Probands <- Pedigree_File %>% filter(hypoglycemia == 1 | T1D == 1)
    #Pedigree_Branchs <- gen.branching(Pedigre_to_genlib, pro = abs(Probands$egoUPDBID))
    # TODO: Learn more about probands to see if it would help the tree to have them listed,
    # just doesnt seem super relevant here...
    Pedigree_Branchs <- gen.branching(Pedigree_to_genlib)
    
    finalPed <- gen.genout(Pedigree_Branchs)
    finalPed <- merge(finalPed, Pedigree_File[ , c("egoUPDBID", "rapid_decliner", "dmTypeNum", "dmType", "slope")], by.x = "ind", by.y = "egoUPDBID")
    #finalPed <- merge(finalPed, Pedigree_File[ , c("egoUPDBID", "rapid_decliner", "dmTypeNum")], by.x = "ind", by.y = "egoUPDBID", all.x = TRUE)
    finalPed$dmTypeNum[is.na(finalPed$dmTypeNum)] <- '0'
    #finalPed$hypoglycemia[is.na(finalPed$hypoglycemia)] <- 0
    
    ######################################################################################
    ####----------------------------------CRANEFOOT------------------------------------###
    ######################################################################################
    
    ####Cranefoot Pedigree
    finalPed$family <- as.character(ped_number)
    finalPed$sex <- ifelse(finalPed$sex == 1, "M", "F")
    
    #### START TO GENERATE PEDIGREE FILES
    
    ###Create cranefoot FILES for all families
    #SGLT2_Families_grp <- group_by(SGLT2_PEDIGREE,FAMILY)
    #SGLT2_Families_summary <- summarize(SGLT2_Families_grp, Members = n())
    #SGLT2_Families_List <- SGLT2_Families_summary$FAMILY
    
    Family_Number <- as.character(ped_number)
    
    Config_File <- c(paste0("PedigreeFile       ", Family_Number, "_Family.txt"),
                     paste0("PedigreeName       ", Family_Number, "_Family"),
                     "SubgraphVariable   family",
                     "NameVariable       ind",
                     "FatherVariable     father",
                     "MotherVariable     mother",
                     "GenderVariable     sex",
                     "TextVariable       ind",
                     "TextVariable       dmType",
                     "TextVariable       slope",
                     "PatternVariable    dmTypeNum",
                     "ColorVariable      rapid_decliner",
                     "ColorInfo          Rapid_Decliner       990000",
                     "PatternInfo        T1D       41",
                     "PatternInfo        T2D       43",
                     "PatternInfo        T3D       43",
                     "PatternInfo        IMPUTED_T1D       21",
                     "PatternInfo        IMPUTED_T2D       23",
                     "PageOrientation    landscape",
                     "PageSize    auto")
    
    writeLines(Config_File, paste0("config.txt"))
    write.table(finalPed,file=paste0(Family_Number,"_Family.txt"),
                quote = FALSE, sep ="\t", row.names = FALSE, col.names = TRUE, na = "")
    
    # make sure cranefoot executable is in the same directory in the D: drive that I'm running this in
    command = "./cranefoot config.txt"
    
    system(command, intern = FALSE,
           ignore.stdout = FALSE, ignore.stderr = FALSE,
           wait = TRUE, input = NULL, timeout = 0)
    
    #CraneFoot_File <- as.data.frame(readLines(paste0(Family_Number,"_Family.ps")))
    CraneFoot_File <- readLines(paste0(Family_Number,"_Family.ps"))
    
    count = 0
    for (i in 1:length(CraneFoot_File)) {
      if (CraneFoot_File[i]=="showpage") { #print(i)
        if(count==0) {
          CraneFoot_File[i] <- "% Line Removed for Printing"
        }
        count=count+1
      }
    }
    writeLines(CraneFoot_File, paste0(Family_Number,"_Family.ps"))
    
    # STOP HERE FOR NOW...
    
    
    #####Convertion DIST_ID to PERSON_ID
    Family_Dist_to_Person_ID <- Probands[,c("DIST_ID","PERSON_ID")]
    Family_Dist_to_Person_ID$DIST_ID <- gsub("DIST_","",Family_Dist_to_Person_ID$DIST_ID)
    Family_Dist_to_Person_ID$PERSON_ID <- paste0("P",Family_Dist_to_Person_ID$PERSON_ID)
    
    write.table(Family_Dist_to_Person_ID, file = paste0("Family_",ped_number,".New_Names"), sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    ####DM proband
    Probands_filtered <- filter(Pedigree_Exomiser,is.na(XLC)==FALSE)
    Probands_filtered$Patient_ID <- paste0("P",Probands_filtered$Patient_ID)
    
    DM_probands <- filter(Probands_filtered,DIABETES==2)
    ESRD_probands <- filter(Probands_filtered,ESRD==2)
    RAPID_DECLINERS_probands <- filter(Probands_filtered,Rapid_Decliner==2)
    
    write.table(DM_probands[2], file = paste0("Fam_",ped_number,".DM"), sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(ESRD_probands[2], file = paste0("Fam_",ped_number,".ESRD"), sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(RAPID_DECLINERS_probands[2], file = paste0("Fam_",ped_number,".RD"), sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
  } else {
    write.table(ped_number, file = paste0("Fam_",ped_number,".NO_PEDIGREE"), sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }###END pedigree no exist else
}###END Pedigree loop  

############################################################################################################
######BAMS LIST

Families_list <-read.csv("Families_list.csv", sep=",")
colnames(Families_list) <- c("ped_number","DIST_IDs","count")
Families_list$DIST_IDs <- gsub(";", ",", Families_list$DIST_IDs)
row.names(Families_list) <- Families_list$ped_number

Selected_Pedigrees <- c(12885, 15681, 1681, 1685, 19787, 22764, 2395, 29921, 30214, 3683, 4067, 47560, 538383, 540746, 553250, 556638, 562325, 564456, 565563, 568085, 575055, 575279, 582371, 5827, 586518, 589261, 594728, 597177, 760212, 770662, 91500)

Selected_Families <- Families_list[Families_list$ped_number %in% Selected_Pedigrees,]

for (ped_number in Selected_Families$ped_number) {    
  Probands <- as.data.frame(unlist(str_split(Selected_Families[as.character(ped_number),"DIST_IDs"],",")))
  colnames(Probands) <- "Probands_list"
  
  Probands$Probands_list <-paste0("reheaded_bams/",as.character(Probands$Probands_list),".bam")
  
  write.table(Probands, file = paste0("Fam_",ped_number,".Bams_list"), sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
}###END for bams  