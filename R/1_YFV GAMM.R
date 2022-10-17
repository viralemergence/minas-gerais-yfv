
# YFV GAMM ####

library(ggregplot); library(tidyverse); library(readxl); library(magrittr); library(cowplot)
library(colorspace); library(chron); library(fs); library(splines); library(coda)
library(rnaturalearth); library(rnaturalearthdata); library(patchwork); library(fs)
library(MASS)

dir_create("Figures")

source("R/0_YFV Import.R")

FamilyList <- c("binomial", "binomial", "gaussian", "gaussian", "binomial", "gaussian")

names(FamilyList) <- Resps

IMList <- list()

ClashList <- list(c("Year", "YearMonth"),
                  c("Month", "YearMonth"))

# Covar <- c("Year", "Month", "Taxon")
Covar <- c("Taxon")

AddCovar <- c("Area of sampling", "NHP carcass preservation status",
              "Year", "Month", "YearMonth",
              "LocalDensity",
              "Population", "GDR_gross_domestic_product") %>% 
  str_replace_all(" ", "_")

# Resps <- Resps[1:4]
Resps <- Resps[2:3]

BAMList <- list()

r <- 1

for(r in r:length(Resps)){
  
  print(Resps[r])
  
  TestDF <- 
    YFV %>% 
    mutate(YearMonth = paste(Year, Month, sep = "_")) %>% 
    dplyr::select(all_of(c(Covar, Resps[r], AddCovar[1:5])), 
                  X, Y, Date) %>% 
    na.omit %>% 
    filter(Year %in% 2017:2018) %>%
    # mutate_at(c("Population", "GDR_gross_domestic_product"), ~log(.x)) %>% 
    mutate_at(c("Month", "Year"), as.factor) %>% 
    mutate_at("Date", as.numeric) %>% 
    # mutate_at("Date", ~.x - min(.x)) %>% 
    droplevels
  
  TestDF %<>% mutate_at(vars(contains("LocalDensity")), ~kader:::cuberoot(.x))
  
  # HostInclude <- 
  #   TestDF %>% 
  #   group_by(Taxon) %>% 
  #   count() %>% 
  #   filter(n > 5, !Taxon == "NA") %>% 
  #   pull(Taxon)
  # 
  # TestDF %<>% filter(Taxon %in% HostInclude)
  
  TestDF %<>% 
    mutate_at("Taxon", ~str_split(.x, "_") %>% map_chr(1)) %>% 
    filter(!Taxon %in% c("NA", "Sapajus", "Cebidae", "Pan"))
  
  TestDF$Response <- TestDF[,Resps[r]]
  
  TestDF %>% nrow %>% print
  
  AddCovar2 <- 
    c(
      glue::glue("s(Date, bs = '{c('ad', 'cr', 'cs', 'tp', 'ts')}')")[1],
      "s(X) + s(Y) + ti(X, Y)"
      # "LocalDensity",
      # "s(LocalDensity)"
    )
  
  # ClashList <- list(AddCovar2[str_detect(AddCovar2, "Date")],
  #                   AddCovar2[str_detect(AddCovar2, "X")],
  #                   AddCovar2[str_detect(AddCovar2, "Density")])
  
  # Maybe try an spde effect?
  
  BAM1 <- BAMModelAdd(Data = TestDF,
                      Response = Resps[r],
                      Explanatory = AddCovar[1:5] %>% 
                        setdiff(c("Year", "Month", "YearMonth")) %>% 
                        c("Taxon") %>% 
                        c(AddCovar2[1:2]), 
                      # Add = AddCovar2[3:4],
                      # Clashes = ClashList,
                      Family = FamilyList[[Resps[r]]],
                      AllModels = F
  )
  
  BAMList[[Resps[r]]] <- BAM1
  
}

# Deviance Contributions ####

DevList <- list()

r <- 1

for(r in r:length(Resps)){
  
  print(Resps[r])
  
  TestDF <- 
    YFV %>% 
    mutate(YearMonth = paste(Year, Month, sep = "_")) %>% 
    dplyr::select(all_of(c(Covar, Resps[r], AddCovar[1:5])), 
                  X, Y, Date) %>% 
    na.omit %>% 
    filter(Year %in% 2017:2018) %>%
    # mutate_at(c("Population", "GDR_gross_domestic_product"), ~log(.x)) %>% 
    mutate_at(c("Month", "Year"), as.factor) %>% 
    mutate_at("Date", as.numeric) %>% 
    # mutate_at("Date", ~.x - min(.x)) %>% 
    droplevels
  
  TestDF %<>% mutate_at(vars(contains("LocalDensity")), ~kader:::cuberoot(.x))
  
  # HostInclude <- 
  #   TestDF %>% 
  #   group_by(Taxon) %>% 
  #   count() %>% 
  #   filter(n > 5, !Taxon == "NA") %>% 
  #   pull(Taxon)
  # 
  # TestDF %<>% filter(Taxon %in% HostInclude)
  
  TestDF %<>% 
    mutate_at("Taxon", ~str_split(.x, "_") %>% map_chr(1)) %>% 
    filter(!Taxon %in% c("NA", "Sapajus", "Cebidae", "Pan"))
  
  TestDF$Response <- TestDF[,Resps[r]]
  
  TestDF %>% nrow %>% print
  
  AddCovar2 <- 
    c(
      glue::glue("s(Date, bs = '{c('ad', 'cr', 'cs', 'tp', 'ts')}')")[1],
      "s(X) + s(Y) + ti(X, Y)"
      # "LocalDensity",
      # "s(LocalDensity)"
    )
  
  ClashList <- list(AddCovar2[str_detect(AddCovar2, "Date")],
                    AddCovar2[str_detect(AddCovar2, "X")],
                    AddCovar2[str_detect(AddCovar2, "Density")])
  
  # Maybe try an spde effect?
  
  Expl <- 
    AddCovar[1:5] %>% 
    setdiff(c("Year", "Month", "YearMonth")) %>% 
    c("Taxon") %>% 
    c(AddCovar2[1:2])
  
  for(q in Expl){
    
    BAM1 <- BAMModelAdd(Data = TestDF,
                        Response = Resps[r],
                        Explanatory = Expl %>% setdiff(q), 
                        # Add = AddCovar2[3:4],
                        Clashes = ClashList,
                        Family = FamilyList[[Resps[r]]],
                        AllModels = F
    )
    
    DevList[[Resps[r]]][[q]] <- BAM1$FinalModel
    
  }
  
}

BAMList[[1]]$FinalModel %>% deviance
DevList[[1]] %>% map(deviance)




