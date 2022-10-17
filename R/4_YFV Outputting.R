
# 4_YFV GAMM Outputting ####

source("R/0_YFV Import.R")

source("R/1_YFV GAMM.R")

# Sampling map? ####

(Map <- 
   YFV %>% 
   ggplot(aes(X, Y)) + 
   geom_point(colour = AlberColours[[3]], alpha = 0.3) +
   # geom_point(aes(colour = as.factor(PCR_LIVER)), alpha = 0.3) +
   geom_sf(data = world %>% filter(name == "Brazil"),
           inherit.aes = F, 
           # lty = 3,
           fill = NA, colour = "black", size = 1) +
   coord_sf() +
   # facet_wrap(~Year) +
   labs(x = NULL, y = NULL))

ggsave("Figures/MapFigure.jpeg", 
       units = "mm", 
       height = 120, width = 120, 
       dpi = 600)

YFV %>% 
  ggplot(aes(X, Y)) + 
  geom_point(colour = AlberColours[[3]], alpha = 0.3) +
  # geom_point(aes(colour = as.factor(PCR_LIVER)), alpha = 0.3) +
  geom_sf(data = world %>% filter(name == "Brazil"),
          inherit.aes = F, 
          # lty = 3,
          fill = NA, colour = "black", size = 1) +
  coord_sf() +
  facet_wrap(~Year) +
  labs(x = NULL, y = NULL)

ggsave("Figures/MapFigure2.jpeg", 
       units = "mm", 
       height = 120, width = 120, 
       dpi = 600)

# Model smooths ####

r = 1

RespLabels <- c("YFV Infection", "YFV Intensity")

names(RespLabels) <- Resps

names(AlberColours)[c(1, 3)] <- c(Resps)

TimePlotList <- SpacePlotList <- EffectDFList <- list()

for(r in r:length(Resps)){
  
  print(Resps[r])
  
  BAM1 <- BAMList[[Resps[r]]]
  
  TestDF <- BAM1$Data
  
  if(FamilyList[[Resps[r]]] == "gaussian"){
    
    TestDF %<>% mutate_at("Response", ~c(scale(.x)))
    
  }
  
  Coefs <- coef(BAM1$FinalModel)
  VC <- vcov(BAM1$FinalModel)
  
  sim <- mvrnorm(1000, mu = Coefs, Sigma = VC)
  
  MCMCEffects <- sim %>% as.data.frame %>%
    map(~.x %>% as.mcmc %>% HPDinterval)
  
  EffectsDF <- 
    summary(BAM1$FinalModel)[1:4] %>% 
    map(~.x[1:length(summary(BAM1$FinalModel)[[1]])]) %>% 
    bind_cols() %>% 
    mutate(Name = summary(BAM1$FinalModel)[[1]] %>% attr("names")) %>%
    rename(Estimate = p.coeff, P = p.pv, Var = `p.t`) %>%
    # mutate_at("se", ~.x * Var) %>% 
    mutate(RawLower = Estimate - se, 
           RawUpper = Estimate + se) %>%
    as.data.frame
  
  EffectsDF <- 
    MCMCEffects[1:nrow(EffectsDF)] %>% 
    map(~tibble(Lower = .x[1], 
                Upper = .x[2])) %>% 
    bind_rows(.id = "Name") %>% 
    left_join(EffectsDF, .)
  
  EffectsDF %<>% bind_rows(
    
    data.frame(Name = c("Area_of_samplingRural", "NHP_carcass_preservation_statusBad", "TaxonAlouatta"),
               Estimate = c(0, 0, 0))
    
  ) %>% #arrange(desc(Name)) %>% 
    mutate_at("Name", ~factor(.x, levels = (unique(.x))))
  
  EffectsDF %<>% 
    mutate(Sig = ifelse(P < 0.05, "*", ""),
           SigLocation = Upper + 0.2)
  
  EffectDFList[[Resps[r]]] <- EffectsDF
  
  # (FixedEffects <-
  #     EffectsDF %>% 
  #     ggplot(aes(Name, Estimate)) + 
  #     geom_hline(yintercept = 0, lty = 2, alpha = 0.4) +
  #     geom_errorbar(aes(ymin = Lower, ymax = Upper),
  #                   width = 0.2) +
  #     geom_point(size = 4) + #scale_x_reverse() +
  #     geom_point(colour = AlberColours[[1]], size = 2.5) + #scale_x_reverse() +
  #     coord_flip() +
  #     geom_text(aes(y = SigLocation, label = Sig),
  #               colour = AlberColours[[1]]) +
  #     scale_x_discrete(limits = c(glue::glue("Taxon{c('Callithrix', 'Callicebus', 'Alouatta')}"),
  #                                 "Gap1",
  #                                 glue::glue("NHP_carcass_preservation_status{c('good', 'intermediate', 'Bad')}"), 
  #                                 "Gap2",
  #                                 rev(glue::glue("Area_of_sampling{c('Rural', 'Urban_rural_interface', 'Urban')}")),
  #                                 "Gap3",
  #                                 # "LocalDensity",
  #                                 # "Gap4",
  #                                 "(Intercept)"),
  #                      labels = c("Taxon: Callithrix",
  #                                 "Taxon: Callicebus",
  #                                 "Taxon: Alouatta", 
  #                                 "",
  #                                 "Carcass: Good",
  #                                 "Carcass: Intermediate", 
  #                                 "Carcass: Bad",
  #                                 "",
  #                                 "Urban",
  #                                 "Urban-Rural Interface", 
  #                                 "Rural", 
  #                                 "",
  #                                 # "Population Density",
  #                                 # "",
  #                                 "Intercept")))
  
  # Getting model outputs ####
  
  XRange <- seq(from = min(TestDF$X) - 1,
                to = max(TestDF$X) + 1,
                length = 50)  %>% 
    c(mean(TestDF$X))
  
  YRange <- seq(from = min(TestDF$Y) - 1,
                to = max(TestDF$Y) + 1,
                length = 50)  %>% 
    c(mean(TestDF$Y))
  
  DateRange <- seq(from = min(TestDF$Date),
                   to = max(TestDF$Date),
                   length = 50)  %>% 
    c(mean(TestDF$Date))
  
  # DensityRange <- seq(from = min(BAM1$Data$LocalDensity),
  #                     to = max(BAM1$Data$LocalDensity),
  #                     length = 50)  %>% 
  # c(mean(BAM1$Data$LocalDensity))
  
  FitList <- expand.grid(
    # HRO = HRORange,
    X = XRange,
    Y = YRange,
    # X = last(XRange),
    # Y = last(YRange),
    Date = DateRange,
    # LocalDensity = DensityRange,
    Area_of_sampling = unique(TestDF$Area_of_sampling),
    NHP_carcass_preservation_status = unique(TestDF$NHP_carcass_preservation_status),
    Taxon = unique(TestDF$Taxon)
    
  )
  
  # FitList %>% filter()
  
  FitPredictions  <- predict.gam(BAM1$FinalModel, 
                                 newdata = FitList, 
                                 se.fit = T)
  
  if(FamilyList[[Resps[r]]] == "binomial"){
    
    FitList[,c("Fit","Lower", "Upper")] <- 
      logistic(with(FitPredictions, cbind(fit, fit - se.fit, fit + se.fit)))
    
  }else{
    
    FitList[,c("Fit","Lower", "Upper")] <-
      with(FitPredictions, cbind(fit, fit - se.fit, fit + se.fit))
    
  }
  
  YearDF <- data.frame(Date = c(min(TestDF$Date) - 9, 
                                min(TestDF$Date) - 9 + 365, 
                                min(TestDF$Date) - 9 + 365 + 365), Fit = 0)
  
  LabelDF <- data.frame(Date = c(min(TestDF$Date), 
                                 min(TestDF$Date) + 365, 
                                 min(TestDF$Date) + 365 + 365), 
                        Fit = 0.935,
                        Label = c("2017", "2018", "2019"),
                        HJust = c(0, 0, 1.25))
  
  Jitter <- 0.25
  
  (TimePlotList[[Resps[r]]] <-
      FitList %>% filter(Area_of_sampling == "Rural", 
                         Taxon == "Callithrix", 
                         NHP_carcass_preservation_status == "good", 
                         # Date == 0, 
                         Y == last(Y), X == last(X)) %>% 
      mutate(Facet = RespLabels[[Resps[r]]]) %>% 
      ggplot(aes(Date, Fit)) + 
      geom_line(colour = AlberColours[[Resps[r]]]) +
      geom_vline(data = YearDF, aes(xintercept = Date), #aes(x = c(0, 365), y = c(0, 0)), 
                 lty = 2, alpha = 0.6, colour = "grey") + 
      geom_ribbon(alpha = 0.3, fill = AlberColours[[Resps[r]]], 
                  aes(ymin = Lower, ymax = Upper)) +
      lims(x = c(min(TestDF$Date) - 10, NA)) +
      labs(y = RespLabels[[Resps[r]]], y = "Day") +
      # theme(axis.title.y = element_text(vjust = -25)) +
      geom_point(data = TestDF %>% mutate(Facet = RespLabels[[Resps[r]]]), 
                 aes(Date, Response), 
                 colour = AlberColours[[Resps[r]]],
                 alpha = 0.3) +
      # geom_text(data = LabelDF %>% mutate(Facet = "Epizootics"), 
      #           aes(Date, Fit*30, label = Label, hjust = HJust), 
      #           # hjust = 0,
      #           alpha = 0.3) +
      # geom_point(data = YearDF, #aes(x = c(0, 365), y = c(0, 0)), 
      #            size = 3, colour = AlberColours[[3]], shape = 2) +
      scale_fill_continuous_sequential(palette = AlberPalettes[[2]])  + 
      # geom_point(data = FullDates %>% mutate(Facet = "Epizootics"), alpha = 0.1,
      #            aes(as.numeric(Date),
      #                total_number_of_animals), #/max(total_number_of_animals)), 
      #            colour = AlberColours[[2]]) + 
      # geom_smooth(data = FullDates %>% mutate(Facet = "Epizootics"), 
      #             colour = AlberColours[[2]], 
      #             fill = AlberColours[[2]], 
      #             alpha = 0.3,
      #             aes(as.numeric(Date), 
      #                 total_number_of_animals)) + #/max(total_number_of_animals))) +
      facet_grid(Facet ~ ., scales = "free") +
      labs(y = ""))
  # plot_annotation(tag_levels = "A") +
  
  ggsave(glue::glue("Figures/OutputFigure1{Resps[r]}.jpeg"), 
         units = "mm",
         # height = 200, 
         height = 120, 
         width = 180)
  
  Jitter <- 0.1
  
  (SpacePlotList[[Resps[r]]] <- FitList %>% 
      mutate_at("Fit", ~ifelse(.x > 1.5 | .x < -1.5, NA, .x)) %>% 
      # filter(Fit < 1.5, Fit > -1.5) %>% 
      filter(Area_of_sampling == "Rural", 
             Taxon == "Callithrix", 
             NHP_carcass_preservation_status == "good", 
             Date == min(Date), !Y == last(Y), !X == last(X)) %>% 
      ggplot(aes(X, Y, fill = Fit)) + 
      geom_tile(alpha = 1) +
      # geom_contour(aes(z = Fit), colour = "white", alpha = 0.8) + 
      # metR::geom_text_contour(aes(z = Fit), colour = "white", size = 2.5, hjust = 0.5, vjust = 1.1, check_overlap = T) +
      geom_point(data = TestDF, inherit.aes = F, 
                 # shape = 2, #colour = "white",
                 aes(X, Y), alpha = 0.35, 
                 position = position_jitter(w = Jitter, h = Jitter)) + #0.6) +
      scale_fill_continuous_sequential(palette = AlberPalettes[[3]], 
                                       limits = c(NA, NA)) +
      geom_sf(data = world %>% filter(name == "Brazil"),
              inherit.aes = F, 
              # lty = 3,
              fill = NA, colour = "black", size = 1.5) +
      coord_sf(xlim = c(-52, -40), 
               ylim = c(-24, -14)) +
      labs(fill = RespLabels[[Resps[r]]], x = NULL, y = NULL)) 
  
  # (SpacePlotList[[Resps[r]]] /
  #     (FixedEffects +
  #        # scale_y_reverse() +
  #        labs(x = NULL))) + 
  #   
  #   plot_layout(heights = c(1, 1)) 
  # 
  # ggsave(glue::glue("Figures/OutputFigure2{Resps[r]}.jpeg"), 
  #        units = "mm", 
  #        height = 200, width = 180)
  
}

# Effects ####

FullEffectDF <- 
  EffectDFList %>% bind_rows(.id = "Response") %>% 
  mutate_at("Response", ~str_split(.x, "_") %>% map_chr(1)) %>% 
  mutate_at("Response", ~str_replace_all(.x, c("qPCR" = "Intensity",
                                               "^PCR$" = "Infection")))

FullEffectDF %>% 
  ggplot(aes(Name, Estimate, colour = Response, group = Response)) + 
  geom_hline(yintercept = 0, lty = 2, alpha = 0.4) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                position = position_dodge(w = 0.5),
                width = 0.2) +
  geom_point(size = 2, colour = "black",
             position = position_dodge(w = 0.5)) + #scale_x_reverse() +
  geom_point(size = 1.5,
             position = position_dodge(w = 0.5)) + #scale_x_reverse() +
  coord_flip() +
  geom_text(#data = FullEffectDF %>% filter(Sig == "*"),
    aes(y = SigLocation, label = Sig),
    # shape = "*",
    size = 4,
    vjust = 0.8,
    show.legend = F,
    position = position_dodge(w = 0.5)) +
  scale_x_discrete(limits = c(glue::glue("Taxon{c('Callithrix', 'Callicebus', 'Alouatta')}"),
                              "Gap1",
                              glue::glue("NHP_carcass_preservation_status{c('good', 'intermediate', 'Bad')}"),
                              "Gap2",
                              rev(glue::glue("Area_of_sampling{c('Rural', 'Urban_rural_interface', 'Urban')}")),
                              "Gap3",
                              # "LocalDensity",
                              # "Gap4",
                              "(Intercept)"),
                   labels = c("Taxon: Callithrix",
                              "Taxon: Callicebus",
                              "Taxon: Alouatta",
                              "",
                              "Carcass: Good",
                              "Carcass: Intermediate",
                              "Carcass: Bad",
                              "",
                              "Urban",
                              "Urban-Rural Interface",
                              "Rural",
                              "",
                              # "Population Density",
                              # "",
                              "Intercept")) +
  scale_colour_manual(
    values = c(AlberColours[[1]],
               AlberColours[[2]],
               AlberColours[[3]],
               AlberColours[[4]])
  ) +
  labs(x = NULL) +
  theme(legend.position = c(0.1, 0.9))

ggsave("Figures/EffectsFigure.jpeg", 
       units = "mm", 
       height = 120, width = 150, 
       dpi = 300)

# Space ####

(Map |
   (SpacePlotList[[1]] + 
      scale_y_continuous(breaks = -c(8:15)*2) +
      scale_fill_continuous_sequential(palette = AlberPalettes[[1]]))|
   SpacePlotList[[2]] + 
   scale_y_continuous(breaks = -c(8:15)*2)) +
  plot_layout(guides = "collect", ncol = 3) +
  plot_annotation(tag_levels = "A")

ggsave("Figures/SpaceFigure.jpeg", 
       units = "mm", 
       height = 120, width = 320, 
       dpi = 300)

# Time ####

(Epizootics <-
  FitList %>% filter(Area_of_sampling == "rural", 
                     Taxon == "Callithrix", 
                     NHP_carcass_preservation_status == "good", 
                     # Date == 0, 
                     Y == last(Y), X == last(X)) %>% 
  mutate(Facet = "Samples") %>% 
  ggplot(aes(Date, Fit)) + geom_line(colour = AlberColours[[1]]) +
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
  geom_text(data = LabelDF %>% mutate(Facet = "Epizootic Cases"), 
            aes(Date, Fit*30, label = Label, hjust = HJust), 
            # hjust = 0,
            alpha = 0.3) +
  # geom_point(data = YearDF, #aes(x = c(0, 365), y = c(0, 0)), 
  #            size = 3, colour = AlberColours[[3]], shape = 2) +
  scale_fill_continuous_sequential(palette = AlberPalettes[[2]])  + 
  geom_point(data = FullDates %>% mutate(Facet = "Epizootic Cases"), alpha = 0.1,
             aes(as.numeric(as.factor(Date)),
                 total_number_of_animals), #/max(total_number_of_animals)), 
             colour = AlberColours[[2]]) + 
  geom_smooth(data = FullDates %>% 
                mutate(Facet = "Epizootic Cases"), 
              colour = AlberColours[[2]], 
              fill = AlberColours[[2]], 
              alpha = 0.3,
              method = "gam", formula = y ~ s(x, bs = 'cs'), 
              method.args = list(family = "nb"),
              aes(as.numeric(as.factor(Date)), 
                  total_number_of_animals)) + #/max(total_number_of_animals))) +
  facet_grid(Facet ~ ., scales = "free") +
  labs(y = ""))

(HumanCases <-
  FitList %>% filter(Area_of_sampling == "rural", 
                     Taxon == "Callithrix", 
                     NHP_carcass_preservation_status == "good", 
                     # Date == 0, 
                     Y == last(Y), X == last(X)) %>% 
  mutate(Facet = "Samples") %>% 
  ggplot(aes(Date, Fit)) + geom_line(colour = AlberColours[[1]]) +
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
  # geom_text(data = LabelDF %>% mutate(Facet = "Epizootic Cases"), 
  #           aes(Date, Fit*30, label = Label, hjust = HJust), 
  #           # hjust = 0,
  #           alpha = 0.3) +
  # # geom_point(data = YearDF, #aes(x = c(0, 365), y = c(0, 0)), 
  # #            size = 3, colour = AlberColours[[3]], shape = 2) +
  # scale_fill_continuous_sequential(palette = AlberPalettes[[2]])  + 
  # geom_point(data = FullDates %>% mutate(Facet = "Epizootic Cases"), alpha = 0.1,
  #            aes(as.numeric(as.factor(Date)),
  #                total_number_of_animals), #/max(total_number_of_animals)), 
  #            colour = AlberColours[[2]]) + 
  # geom_smooth(data = FullDates %>% 
  #               mutate(Facet = "Epizootic Cases"), 
  #             colour = AlberColours[[2]], 
  #             fill = AlberColours[[2]], 
  #             alpha = 0.3,
  #             aes(as.numeric(as.factor(Date)), 
  #                 total_number_of_animals)) + #/max(total_number_of_animals))) +
  facet_grid(Facet ~ ., scales = "free") +
  labs(y = "") +
  geom_text(data = LabelDF %>% mutate(Facet = " Human YFV Cases"), 
            aes(Date, Fit*30, label = Label, hjust = HJust), 
            # hjust = 0,
            alpha = 0.3) +
  geom_point(data = HumanDates %>% mutate(Facet = " Human YFV Cases"), alpha = 0.1,
             aes(as.numeric(as.factor(Date)),
                 Count), #/max(total_number_of_animals)), 
             colour = AlberColours[[5]]) + 
  geom_smooth(data = HumanDates %>% 
                filter(Date < "2018-06-01") %>% 
                mutate(Facet = " Human YFV Cases"), 
              colour = AlberColours[[5]], 
              fill = AlberColours[[5]], 
              alpha = 0.3,
              # method = "loess", 
              # method.args = list(family = "nb"),
              method = "gam", formula = y ~ s(x, bs = 'cs'), 
              method.args = list(family = "nb"),
              aes(as.numeric(as.factor(Date)), 
                  Count)) +
    lims(y = c(NA, 30)))

(Epizootics/TimePlotList[[1]]/TimePlotList[[2]]/HumanCases) +
  plot_annotation(tag_levels = "A")

ggsave("Figures/TimeFigure.jpeg", units = "mm", height = 200, width = 150, dpi = 300)

# Results ####

Model <- BAMList[[1]]$FinalModel

Table <- 
  BAMList %>% 
  map("FinalModel") %>% 
  map(function(Model){
    
    DF <- 
      Model %>% 
      summary %>% 
      extract(c(1, 2, 4)) %>% 
      map(~.x[1:7]) %>% 
      bind_cols() %>% 
      rename(Estimate = 1, SE = 2, P = 3) %>% 
      mutate_all(~round(.x, 3)) %>% 
      rownames_to_column("Effect") %>% 
      mutate(Effect = 
               Model %>% 
               summary %>% 
               extract(c(1)) %>% 
               unlist %>% 
               names %>% str_remove("^p.coeff.")) %>% 
      data.frame %>% 
      mutate(Sig = ifelse(P < 0.001, "***", ifelse(P < 0.05, "*", ""))) %>% 
      mutate_at("P", ~ifelse(.x < 0.001, "<0.001", .x)) %>% 
      mutate_at("P", ~paste0(P, Sig)) %>% 
      dplyr::select(-Sig)
    
    DF2 <- 
      data.frame(
        
        P = Model %>% summary %>% extract("s.pv"),
        Effect = Model %>% summary %>% extract2("chi.sq") %>% names
        
      )
    
    DF %>% bind_rows(DF2)
    
  })

NameReplace <-  c("Taxon: Callithrix",
                  "Taxon: Callicebus",
                  "Taxon: Alouatta",
                  "Carcass: Good",
                  "Carcass: Intermediate",
                  "Carcass: Bad",
                  "Urban",
                  "Urban-Rural Interface",
                  "Rural",
                  "Intercept")

names(NameReplace) <- c(glue::glue("Taxon{c('Callithrix', 'Callicebus', 'Alouatta')}"),
                        glue::glue("NHP_carcass_preservation_status{c('good', 'intermediate', 'Bad')}"),
                        rev(glue::glue("Area_of_sampling{c('Rural', 'Urban_rural_interface', 'Urban')}")),
                        "(Intercept)")

Table[[1]]$Effect %<>% str_replace_all(NameReplace)
Table[[2]]$Effect %<>% str_replace_all(NameReplace)

Table %>% write.csv("EstimateTable.csv", row.names = F)

# Table ####

BAMList %>% 
  map("FinalModel") %>% 
  map(summary)

# Deviance estimates ####

DataList <- BAMList %>% map("Data")

RandomPredictionList <- DevianceList <- list()

RealPredictions <- InterceptPredictions <- list()

RepCovar <- list("Area_of_sampling", "NHP_carcass_preservation_status", "Taxon", 
                 "Date", c("X", "Y"))

y <- Resps[1]

Iterations <- 100

for(y in Resps){
  
  print(y)
  
  RandomPredictionList[[y]] <- DevianceList[[y]] <- list()
  
  RealOutcomes <- DataList[[y]][,y]
  
  RealPredictions[[y]] <- predict.bam(BAMList[[y]]$FinalModel, 
                                      newdata = DataList[[y]])# %>% logistic
  
  InterceptPredictions[[y]] <- rep(mean(RealPredictions[[y]]), nrow(DataList[[y]]))
  
  x <- RepCovar[[1]]
  
  for(x in RepCovar){
    
    print(x)
    
    for(i in 1:Iterations){
      
      print(i)
      
      PredDF <- DataList[[y]]
      
      PredDF[,x] <- PredDF %>% slice(sample(1:n())) %>% 
        dplyr::select(all_of(x))
      # apply(2, mean)
      
      Predictions <- predict.bam(BAMList[[y]]$FinalModel, 
                                 newdata = PredDF)
      
      RandomPredictionList[[y]][[paste(x, collapse = "_")]][[i]] <- Predictions
      
    }
  }
}

MeanSD <- function(a){
  
  list(Mean = mean(a), 
       SD = sd(a), 
       N = length(a), 
       SE = sd(a)/(length(a)^0.5))
  
}

y <- Resps[1]

PCRDeviances <- 
  RandomPredictionList[[y]] %>% 
  map(function(a){
    
    map(a, function(b){
      
      cor(b, RealPredictions[[y]], method = "spearman")
      
    }) %>% unlist
    
  }) %>% 
  bind_cols() %>% 
  map(MeanSD) %>% 
  map("Mean")

cor(DataList[[y]][,y], RealPredictions[[y]], method = "spearman")

PCRDeviances %<>% 
  unlist %>% add(-1) %>% prop.table %>% 
  multiply_by(cor(DataList[[y]][,y], RealPredictions[[y]], method = "spearman")) %>% 
  as.data.frame %>% 
  rownames_to_column() %>% 
  rename(Var = 1, Estimate = 2) %>% 
  mutate(Response = "Infection")

y <- Resps[2]

qPCRDeviances <- 
  RandomPredictionList[[y]] %>% 
  map(function(a){
    
    map(a, function(b){
      
      cor(b, RealPredictions[[y]], method = "spearman")
      
    }) %>% unlist
    
  }) %>% 
  bind_cols() %>% 
  map(MeanSD) %>% 
  map("Mean")

cor(DataList[[y]][,y], RealPredictions[[y]], method = "spearman")

qPCRDeviances %<>% 
  unlist %>% add(-1) %>% prop.table %>% 
  multiply_by(cor(DataList[[y]][,y], RealPredictions[[y]], method = "spearman")) %>% 
  as.data.frame %>% 
  rownames_to_column() %>% 
  rename(Var = 1, Estimate = 2) %>% 
  mutate(Response = "Intensity")

(DevianceDF <- 
    PCRDeviances %>% 
    bind_rows(qPCRDeviances) %>% 
    dplyr::select(Response, Var, R = Estimate)) %>% 
  pivot_wider(names_from = "Response", values_from = "R")

DevianceDF %>% 
  # write.csv("DevianceDF.csv", row.names = F) %>% 
  write.csv("Table_Y.csv", row.names = F)

PrintEffectDF <- 
  BAMList %>% map("FinalModel") %>% map(summary) %>% map("s.table") %>% 
  map(c(as.data.frame, rownames_to_column)) %>% 
  bind_rows(.id = "Response") %>% 
  mutate_at("Response", ~str_replace_all(.x, c("^PCR_LIVER" = "Infection",
                                               "^qPCR_LIVER" = "Intensity"))) %>% 
  dplyr::select(Response, Variable = rowname, P = `p-value`, edf, Chi_or_F = 4) %>% 
  bind_rows(FullEffectDF %>% 
              dplyr::select(Response, Variable = Name, 
                            Estimate, Lower95 = Lower, Upper95 = Upper, 
                            P, Sig) %>% 
              mutate_at("Variable", as.character) %>% 
              arrange(Variable), .) %>% 
  arrange(Response) %>%  
  mutate(Sig = ifelse(P < 0.05, "*", "")) %>% 
  mutate_if(is.numeric, ~round(.x, 3)) %>% 
  mutate_all(as.character) %>% 
  mutate_all(~replace_na(.x, "")) %>% 
  mutate_at("P", ~str_replace_all(.x, "^0$", "<0.001"))

NameReplace <-  c("Taxon: Callithrix",
                  "Taxon: Callicebus",
                  "Taxon: Alouatta",
                  "Carcass: Good",
                  "Carcass: Intermediate",
                  "Carcass: Bad",
                  "Urban",
                  "Urban-Rural Interface",
                  "Rural",
                  "Intercept", 
                  "Date",
                  "Longitude", 
                  "Latitude", 
                  "Long:Lat")

names(NameReplace) <- c(glue::glue("Taxon{c('Callithrix', 'Callicebus', 'Alouatta')}"),
                        glue::glue("NHP_carcass_preservation_status{c('good', 'intermediate', 'Bad')}"),
                        rev(glue::glue("Area_of_sampling{c('Rural', 'Urban_rural_interface', 'Urban')}")),
                        "[(]Intercept[)]", 
                        "s[(]Date[)]", 
                        "s[(]X[)]", "s[(]Y[)]", "ti[(]X,Y[)]")

PrintEffectDF %>% 
  mutate_at("Variable", ~str_replace_all(.x, NameReplace)) %>% 
  dplyr::select(Response, Variable, 
                Estimate, Lower95, Upper95, 
                Sig, edf, Chi_or_F, P) %>% 
  # write.csv("EffectsDF.csv", row.names = F) %>% 
  write.csv("Table_X.csv", row.names = F)

