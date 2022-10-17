
# 5_YFV Results ####

# pseudo-R2

y <- Resps[1]
cor(DataList[[y]][,y], RealPredictions[[y]], method = "spearman")

y <- Resps[2]
cor(DataList[[y]][,y], RealPredictions[[y]], method = "spearman")

# The temporal smooth was both highly significant (P_) 

PrintEffectDF %>% filter(str_detect(Var, "Date"))

# and accounted for a large proportion of the model deviance (R_). 

DevianceDF %>% filter(str_detect(Var, "Date"))

# There were strong spatial patterns of prevalence and intensity which were highly significant (P_) 

PrintEffectDF %>% filter(str_detect(Var, "X"))

PrintEffectDF %>% filter(str_detect(Var, "Y"))

# and accounted for large amounts of deviance (R_). 

DevianceDF %>% filter(str_detect(Var, "X"))

DevianceDF %>% filter(str_detect(Var, "Y"))

# Samples taken from the urban-rural interface had lower intensity (P_), 

PrintEffectDF %>% filter(str_detect(Var, "Area"))

# but not lower prevalence (P>0.05), 
# than those taken from rural areas. Samples from fully urban areas had both lower intensity (P_) 

# and prevalence (P_). 

# The genus Callithrix (marmosets) had lower prevalence and intensity than those in the genus Alouatta 
# (howler monkeys; P_ and P_ respectively) 

PrintEffectDF %>% filter(str_detect(Var, "Taxon"))

# and Callicebus (titis; P_ and P_ respectively).


