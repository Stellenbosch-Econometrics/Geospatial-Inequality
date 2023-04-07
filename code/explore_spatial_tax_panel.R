###########################
# Explore Spatial Tax Panel
###########################

library(fastverse)

STP3_files <- list.files("data/spatial_tax_panel/Spatial_Tax_Panel_v3/Data", full.names = TRUE)
names(STP3_files) <- flast(strsplit(STP3_files, "/", fixed = TRUE)) %>% sub(pat = ".csv", rep = "")
STP3_Municipal <- STP3_files[names(STP3_files) %ilike% "Gini|MedianIncome" & startsWith(names(STP3_files), "Municipal")] %>% lapply(fread)
STP3_H7 <- STP3_files[names(STP3_files) %ilike% "Gini|MedianIncome" & startsWith(names(STP3_files), "hex7")] %>% lapply(fread)

# FTE = The number of full time equivalent (FTE) employees.
# gini = The gini coefficient calculated based on the weighted incomes of all employees within an aggregation.

# Spatial Tax Panel: Municipal
MUN <- fread("data/spatial_tax_panel/Spatial_Tax_Panel_v3/Data/Municipal_MedianIncome.csv") %>% fselect(-FTE) %>% 
       merge(fread("data/spatial_tax_panel/Spatial_Tax_Panel_v3/Data/Municipal_Gini.csv"), by = .c(CAT_B, TaxYear)) %>% 
       fmutate(FTE = as.double(FTE))

descr(MUN)
hist(MUN$gini, breaks = 100)

MUN_wide <- MUN %>% dcast(CAT_B ~ TaxYear, value.var = .c(MedianIncome, FTE, gini)) %T>% 
  fwrite("data/spatial_tax_panel/Spatial_Tax_Panel_v3/Data/Municipal_MedianIncome_Gini_wide.csv")

# Spatial Tax Panel: Hexagons
H7 <- fread("data/spatial_tax_panel/Spatial_Tax_Panel_v3/Data/hex7_MedianIncome.csv") %>% fselect(-FTE, -CAT_B) %>% 
      merge(fread("data/spatial_tax_panel/Spatial_Tax_Panel_v3/Data/hex7_Gini.csv"), by = .c(hex7, TaxYear)) %>% 
      fmutate(FTE = as.double(FTE))

descr(H7)
hist(H7$gini, breaks = 100)

H7_wide <- H7 %>% dcast(hex7 ~ TaxYear, value.var = .c(MedianIncome, FTE, gini)) %T>% 
  fwrite("data/spatial_tax_panel/Spatial_Tax_Panel_v3/Data/hex7_MedianIncome_Gini_wide.csv")
