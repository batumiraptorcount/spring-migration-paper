##################
### Read packages 
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("dplyr","tidyr","data.table","ggplot2","lubridate","maptools","RColorBrewer","tidyverse",
              "ggridges","cowplot","openxlsx","showtext","sysfonts","gridExtra")
ipak(packages)

################################################################################
### PREPARE DATA BOTH SEASONS ###
################################################################################

## READ FILTERED DATA ##
spring.nr <- read.csv("data - prepared/BRC_Spring_Non-raptor-data_FILTERED.csv", header=TRUE, stringsAsFactors=FALSE)
autumn.nr <- read.csv("data - prepared/BRC_Autumn_Non-raptor-data_FILTERED.csv", header=TRUE, stringsAsFactors=FALSE)

### ADD COLUMN DEFINING SEASONS ###
season <- c('spring') #spring
spring.nr$season <- season
season <- c('autumn') #autumn
autumn.nr$season <- season

### SUBSET RECENT AUTUMN YEARS ###
autumn.nr.recent <- select(filter(autumn.nr, year %in% c(2018, 2019, 2021)),c(X:season)) #subset daily totals

### NEW DATAFRAME SPRING AND AUTUMN COMBINED ###
comb_non_raptors <- rbind(autumn.nr.recent, spring.nr)

### CALCULATE CUMULATIVE ANNUAL STATISTICS ###
comb_non_raptors_fin <- comb_non_raptors %>% 
  group_by(species,year,season) %>% 
  summarise(yr.tot = sum(nr))

### ADD VALUES FOR MISSING YEARS ###
#define season & year combinations
required_years <- tibble(
  season = c("spring", "spring", "spring", "autumn", "autumn", "autumn"),
  year = c(2019, 2020, 2022, 2018, 2019, 2021)
)

#define species names
species_list <- tibble(species = unique(comb_non_raptors_fin$species))

#combine species with possible season/year
expanded_data <- species_list %>%
  expand_grid(required_years)

#combine annual totals and add zeros for missing years
comb_non_raptors_fin_complete <- expanded_data %>%
  full_join(comb_non_raptors_fin, by = c("species", "season", "year")) %>%
  replace_na(list(yr.tot = 0))

### CALCULATE MEAN SEASONAL STATISTICS ###
comb_non_raptors_fin_seasonal <- comb_non_raptors_fin_complete %>%
  group_by(species,season) %>%
  summarise(mean.yrtot  = round(mean(yr.tot[yr.tot > 0])), 
            sd.yrtot    = round(sd(yr.tot[yr.tot > 0])),
            n.yrs       = length(year[yr.tot > 0])
            )

### Save table S1
write.xlsx(comb_non_raptors_fin_seasonal, "table.S1.nonraptors.xlsx")





