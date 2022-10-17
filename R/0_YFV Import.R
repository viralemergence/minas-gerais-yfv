
# 0_YFV Import ####

rm(list = ls())

library(ggregplot); library(tidyverse); library(readxl); library(magrittr); library(cowplot)
library(colorspace); library(chron); library(fs); library(patchwork)
library(rnaturalearth); library(rnaturalearthdata)

world <- ne_countries(scale = "medium", returnclass = "sf")

dir_create("Output Files")

theme_set(theme_cowplot() + theme(strip.background = element_rect(fill = "white")))

# Importing and exploring ####

Path <- "Data" %>% list.files(full.names = T) %>% last

FileList <- 
  Path %>% 
  last %>% 
  excel_sheets() %>% 
  purrr::set_names() %>% 
  map(read_excel, path = Path) %>% 
  map(data.frame) %>% 
  map(~.x %>% rename_all(function(a) str_replace_all(a, "[.]", "_")))

FileList[[1]] %>% 
  GregHeader ->
  YFV

YFV %<>% dplyr::select(1:`adress complement`)

YFV %<>% rename_all(~.x %>% str_trim %>% str_replace_all(" ", "_") %>% CamelConvert)

YFV %<>% rename(PCR_YFV = `PCR_YFV_(liver_or_brain)`)

Resps <- 
  YFV %>% dplyr::select(contains("YFV")) %>% 
  names #%>% 
# setdiff("BRAIN_YFV_PCR") %>% 
# setdiff("PCR_YFV_(liver_or_brain)")

Resps

# RespRename <- Resps
# 
# names(RespRename) <- c("PCR_Overall", "PCR_LIVER", "qPCR_LIVER", "YFV_Per_Gram", 
#                        "PCR_BRAIN", "qPCR_BRAIN")

RespRename <- c("PCR_Overall", "PCR_LIVER", "qPCR_LIVER", "YFV_Per_Gram", 
                "PCR_BRAIN", "qPCR_BRAIN")# %>% 
# paste0("^", ., "$")

names(RespRename) <- paste0("^", Resps, "$")

YFV %<>% rename_all(~str_replace_all(.x, RespRename))

names(YFV)[names(YFV) == "ESTIMATED_yfv_RNA_COPIES_PER_GRAM_OF_LIVER_(10E_(-0.2986_x_Cq_+_13.5094))"] <-
  "YFV_Per_Gram"

Resps <-
  YFV %>% dplyr::select(PCR_Overall:qPCR_BRAIN) %>% names

YFV %<>% mutate_at(Resps[c(1, 2, 5)], ~.x %>% str_replace_all(c("POS" = "1", "NEG" = "0")) %>% as.numeric)

YFV %<>% mutate_at(Resps, ~.x %>% as.character %>% as.numeric)

YFV %<>% 
  mutate_at(c("Longitude", "Latitude"), 
            ~as.numeric(as.character(.x))) %>% # %>% round(8)) %>% 
  rename(X = Longitude, 
         Y = Latitude)

YFV %<>% 
  mutate_at("Species,_genus_or_family", ~str_replace_all(.x, " ", "_")) %>% 
  rename(Taxon = `Species,_genus_or_family`)

YFV %<>% dplyr::select(-Year)

# Sorting out dates ####

SplitDates1 <- 
  YFV %>% 
  mutate_at("Date", ~.x %>% as.character %>% as.numeric) %>% 
  filter(!is.na(Date)) %>% 
  pull(Date) %>% 
  month.day.year(c(month = 1, 
                   day = 1,
                   year = 1900)) %>% 
  bind_cols() %>% 
  rename_all(CamelConvert)

SplitDates2 <- 
  YFV %>% separate(Date, into = c("Month", "Day", "Year"), sep = "/") %>% 
  mutate_at(c("Month", "Day", "Year"), 
            ~as.numeric(as.character(.x))) %>% 
  filter(Month < 13)

YFV <- 
  YFV %>% 
  mutate_at("Date", ~.x %>% as.character %>% as.numeric) %>% 
  filter(!is.na(Date)) %>% 
  bind_cols(SplitDates1) %>% 
  bind_rows(SplitDates2)

YFV$Month %>% qplot

YFV %<>% 
  mutate(Date = lubridate::ymd(paste(Year, Month, Day, sep = "-"))) %>% 
  arrange(Date)

YFV %<>% mutate(YearMonth = paste(Year, Month, sep = "_"))

# Attaching other data ####

FileList$`Brazilian cities pop GRD` %>% 
  dplyr::select(Municipality, GDR_gross_domestic_product) %>% 
  left_join(YFV, ., by = "Municipality") ->
  YFV

YFV %<>% 
  mutate_at(c("Population", "GDR_gross_domestic_product"), 
            ~as.numeric(as.character(.x))) %>% 
  mutate_at(c("Area_of_sampling", "NHP_carcass_preservation_status"), 
            ~ifelse(.x == "NA", NA, .x)) %>% 
  mutate_at(c("Area_of_sampling"), 
            ~str_replace_all(.x, "urban-rural", "urban_rural") %>% 
              str_replace_all(" ", "_") %>% 
              CamelConvert %>% 
              factor(levels = c("Rural", "Urban_rural_interface", "Urban")))

# Adding density ####
# https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-rev11

# "Data" %>% dir_ls(regexp = "zip$") %>% unzip(exdir = "Data")

DensityRaster <- "Data" %>% dir_ls(regexp = "tif$") %>% raster::raster()

YFV$LocalDensity <- DensityRaster %>% raster::extract(YFV[,c("X", "Y")])

YFV %>% ggplot(aes(X, Y, colour = log(LocalDensity))) + geom_point()

# Importing epizootics ####

Path <- "Data" %>% list.files(full.names = T) %>% first

FileList <- 
  Path %>% 
  last %>% 
  excel_sheets() %>% 
  purrr::set_names() %>% 
  map(read_excel, path = Path) %>% 
  map(data.frame) %>% 
  map(~.x %>% rename_all(function(a) str_replace_all(a, "[.]", "_")))

FileList[[1]] %>% 
  ggplot(aes(long, lat)) + 
  geom_point(aes(colour = total_number_of_animals)) + coord_sf() +
  labs(colour = "Animals")

SummedEpizootics <- FileList[[1]] %>% 
  group_by(date_of_occurrence) %>% 
  summarise_at("total_number_of_animals", sum) %>% 
  filter(date_of_occurrence > "2015-01-01") %>% 
  rename(Date = date_of_occurrence) %>% 
  mutate_at("Date", ~lubridate::ymd(.x))

FullDates <- 
  (lubridate::ymd("2017-01-01")):(lubridate::ymd("2019-01-01")) %>% 
  month.day.year(c(month = 1, 
                   day = 1,
                   year = 1970)) %>% 
  bind_cols() %>% 
  mutate(Date = paste(year, month, day, sep = "-") %>% lubridate::ymd()) %>% 
  dplyr::select(Date) %>% left_join(SummedEpizootics) %>% 
  mutate_at("total_number_of_animals", 
            ~replace_na(.x, 0))


ggplot() + 
  geom_point(data = FullDates, 
             aes(Date, total_number_of_animals/max(total_number_of_animals))) + 
  geom_smooth(data = FullDates, 
              aes(Date, total_number_of_animals/max(total_number_of_animals)))

if(0){
  
  FitList %>% filter(Area_of_sampling == "rural", 
                     Taxon == "Callithrix", 
                     NHP_carcass_preservation_status == "good", 
                     # Date == 0, 
                     Y == last(Y), X == last(X)) %>% 
    mutate(Facet = "Samples") %>% 
    ggplot(aes(Date, Fit)) + geom_line(colour = AlberColours[[1]]) +
    geom_vline(data = YearDF, aes(xintercept = Date), #aes(x = c(0, 365), y = c(0, 0)), 
               lty = 2, alpha = 0.6, colour = "grey") + 
    geom_ribbon(alpha = 0.3, fill = AlberColours[[1]], 
                aes(ymin = Lower, ymax = Upper)) +
    lims(x = c(min(TestDF$Date) - 10, NA)) +
    labs(y = "YFV Prevalence", y = "Day") +
    # theme(axis.title.y = element_text(vjust = -25)) +
    geom_point(data = TestDF %>% mutate(Facet = "Samples"), aes(Date, Response), 
               colour = AlberColours[[1]],
               alpha = 0.3) +
    geom_text(data = LabelDF %>% mutate(Facet = "Epizootics"), 
              aes(Date, Fit*30, label = Label, hjust = HJust), 
              # hjust = 0,
              alpha = 0.3) +
    # geom_point(data = YearDF, #aes(x = c(0, 365), y = c(0, 0)), 
    #            size = 3, colour = AlberColours[[3]], shape = 2) +
    scale_fill_continuous_sequential(palette = AlberPalettes[[2]])  + 
    geom_point(data = FullDates %>% mutate(Facet = "Epizootics"), alpha = 0.1,
               aes(as.numeric(Date),
                   total_number_of_animals), #/max(total_number_of_animals)), 
               colour = AlberColours[[2]]) + 
    geom_smooth(data = FullDates %>% mutate(Facet = "Epizootics"), 
                colour = AlberColours[[2]], 
                fill = AlberColours[[2]], 
                alpha = 0.3,
                aes(as.numeric(Date), 
                    total_number_of_animals)) + #/max(total_number_of_animals))) +
    facet_grid(Facet ~ ., scales = "free") +
    labs(y = "")
  
}

# Turning qPCR into non-nonsense ####

YFV %<>% mutate_at(Resps[str_detect(Resps, "qPCR")], ~-.x)

YFV %<>% 
  mutate(OldDate = Date) %>% 
  mutate_at("Date", ~as.numeric(.x - min(.x)) + 11)

# Human Cases ####

Humans <- "Data" %>% dir_ls(regex = "HUMAN") %>% last %>% read_xls %>% as.data.frame

Humans %<>% rename_all(~.x %>% str_trim %>% tolower %>% str_replace_all(" ", "_") %>% CamelConvert)

HumanCases <- 
  Humans %>% #head
  group_by(`Yf_symptons_onset`) %>% 
  rename(OnsetDate = Yf_symptons_onset) %>% 
  # summarise_at("total_number_of_animals", sum) %>% 
  count %>% rename(Count = n) %>% 
  filter(OnsetDate > "2015-01-01") %>% arrange(OnsetDate) %>% 
  # rename(Date = OnsetDate) %>% 
  mutate_at("OnsetDate", ~lubridate::ymd(.x))

HumanDates <- 
  (lubridate::ymd("2017-01-01")):(lubridate::ymd("2019-01-01")) %>% 
  month.day.year(c(month = 1, 
                   day = 1,
                   year = 1970)) %>% 
  bind_cols() %>% 
  mutate(Date = paste(year, month, day, sep = "-") %>% lubridate::ymd()) %>% 
  dplyr::select(Date) %>% left_join(HumanCases, by = c("Date" = "OnsetDate")) %>% 
  mutate_at("Count", 
            ~replace_na(.x, 0))

HumanDates %>% 
  ggplot(aes(Date, Count)) + #geom_line(colour = AlberColours[[1]]) +
  geom_vline(data = YearDF, aes(xintercept = Date), #aes(x = c(0, 365), y = c(0, 0)), 
             lty = 2, alpha = 0.6, colour = "grey") + 
  # geom_ribbon(alpha = 0.3, fill = AlberColours[[1]], 
  #             aes(ymin = Lower, ymax = Upper)) +
  lims(x = c(min(TestDF$Date) - 10, NA)) +
  labs(y = "YFV Infection", y = "Day") +
  # theme(axis.title.y = element_text(vjust = -25)) +
  # geom_point(data = TestDF %>% mutate(Facet = "Samples"), aes(Date, Response), 
  #            colour = AlberColours[[1]],
  #            alpha = 0.3) +
  geom_text(data = LabelDF %>% mutate(Facet = "Human YFV Cases"), 
            aes(Date, Fit*30, label = Label, hjust = HJust), 
            # hjust = 0,
            alpha = 0.3) +
  geom_point(data = HumanDates %>% mutate(Facet = "Human YFV Cases"), alpha = 0.1,
             aes(as.numeric(as.factor(Date)),
                 Count), #/max(total_number_of_animals)), 
             colour = AlberColours[[5]]) + 
  geom_smooth(data = HumanDates %>% 
                mutate(Facet = "Human YFV Cases"), 
              colour = AlberColours[[5]], 
              fill = AlberColours[[5]], 
              alpha = 0.3,
              method = "loess", 
              # method.args = list(family = "nb"),
              aes(as.numeric(as.factor(Date)), 
                  Count)) + #/max(total_number_of_animals))) +
  facet_grid(Facet ~ ., scales = "free") +
  labs(y = "")
