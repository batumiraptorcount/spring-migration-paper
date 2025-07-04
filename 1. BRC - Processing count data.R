##################
### Read packages 
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("data.table","ggplot2","lubridate","maptools","RColorBrewer","tidyverse",
              "ggridges","cowplot","openxlsx","showtext","sysfonts","gridExtra")
ipak(packages)

################################################################################
                     ### PREPARE DATA BOTH SEASONS ###
################################################################################
## READ FILTERED DATA ##
spring.daily <- read.csv("data - prepared/Output - spring/Daily totals - FILTERED/BRC_Spring_Daily-totals_FILTERED.csv", header=TRUE, stringsAsFactors=FALSE)
spring.annual <- read.csv("data - prepared/Output - spring/Annual totals - FILTERED/BRC_Spring_Annual-totals_FILTERED.csv", header=TRUE, stringsAsFactors=FALSE)
autumn.daily <- read.csv("data - prepared/Output - autumn/Daily totals - FILTERED/BRC_Autumn_Daily-totals_FILTERED.csv", header=TRUE, stringsAsFactors=FALSE)
autumn.annual <- read.csv("data - prepared/Output - autumn/Annual totals - FILTERED/BRC_Autumn_Annual-totals_FILTERED.csv", header=TRUE, stringsAsFactors=FALSE)

### ADD COLUMN DEFINING SEASONS ###
season <- c('spring') #spring
spring.daily$season <- season
spring.annual$season <- season   
season <- c('autumn') #autumn
autumn.daily$season <- season
autumn.annual$season <- season

### SUBSET RECENT AUTUMN YEARS ###
autumn.daily.recent <- select(filter(autumn.daily, year %in% c(2018, 2019, 2021)),c(X:season)) #subset daily totals
autumn.annual.recent <- select(filter(autumn.annual, year %in% c(2018, 2019, 2021)),c(X:season)) #subset annual totals

### NEW DATAFRAME SPRING AND AUTUMN COMBINED ###
comb_daily <- rbind(autumn.daily.recent, spring.daily) #combine daily totals
comb_annual <- rbind(autumn.annual.recent, spring.annual) #combine annual totals

### MERGE DAILY AND ANNUAL TOTALS ### 
comb_fin <- merge(comb_daily,comb_annual[,c("englishname","year","yrtot","season")],all.x=TRUE)

### Set yrtot to 0 if yrtot = NA ###
comb_fin$yrtot <- ifelse(is.na(comb_fin$yrtot)==TRUE,0,comb_fin$yrtot)

### ORDER DATA BASED ON SPECIES AND DATE ###
comb_fin <- comb_fin[order(comb_fin$englishname,comb_fin$date),]

### CALCULATE CUMULATIVE ANNUAL STATISTICS ###
comb_fin <- comb_fin %>% 
  group_by(englishname,year,season) %>% 
  mutate(cumtot = cumsum(ifelse(is.na(daytot)==TRUE,0,daytot)),
         cumprop = cumtot/yrtot,
         yr.daymax = max(daytot,na.rm=TRUE))%>%
  ungroup()

################################################################################
           ### DETECT AND FILTER EARLY SPECIES && LATE DAYS SPRING ###
################################################################################

### Calculate proportion of birds passing before 21st March in 2020/2022
percent.early.passage <- comb_fin %>%
  filter(year>2019 & season=="spring") %>%
  group_by(englishname,latinname,year) %>%
  reframe(prop.early  = (cumtot[yday==80]/yrtot)*100) %>% 
  distinct()

### Calculate mean proportion and add indication of early species
mean.percent.early.passage <- percent.early.passage %>%
  group_by(englishname,latinname) %>%
  reframe(mean.prop.early = mean(prop.early), #calculate mean proportion
          early.species = ifelse(mean.prop.early > 3, TRUE, FALSE)) %>% #indicate if mean is above 3%
  distinct()

### Save list of early species as calculated above
early.species <- as.character(mean.percent.early.passage$englishname[mean.percent.early.passage$early.species == TRUE])

### Filter out 2019 data for early species 
comb_fin <- comb_fin %>%
  filter(!(season == "spring" & englishname %in% early.species & year == 2019))

comb_daily <- comb_daily %>%
  filter(!(season == "spring" & englishname %in% early.species & year == 2019))

### Filter out 2019 data after 26th May 
comb_fin <- comb_fin %>%
  filter(!(season == "spring" & year == 2019 & yday > 146))
  
comb_daily <- comb_daily %>%
  filter(!(season == "spring" & year == 2019 & yday > 146))

### SAVE FILTERED BASE DATAFRAMES ###
saveRDS(comb_daily, "comb_daily.RDS") # DAILY totals both seasons combined
saveRDS(comb_annual, "comb_annual.RDS") # ANNUAL totals both seasons combined
saveRDS(comb_fin, "comb_fin.RDS") # daily & annual & seasons combined

################################################################################
            ### Set various species selections - use where needed ###
################################################################################

non.raptors <- c("Stork spec","Black stork", "White stork") #add others if needed

uid <- c("Common/Lesser kestrel","buzzard spec","Stork spec", "hawk spec",
         "Falcon spec","Ringtail harriers", "eagle spec", "raptor spec", 
         "Hobby/Red-footed falcon", "harrier spec", "raptor spec (medium)",
         "honey buzzard spec") #use for unidentified species

test.species <- c("Montagu's harrier","Pallid harrier","Marsh harrier",
                  "Black kite","Honey buzzard","Booted eagle","Short-toed eagle",
                  "Lesser spotted eagle","Steppe buzzard","Ringtail harriers", 
                  "eagle spec", "raptor spec (medium)") #use for table 2 & 3 (and tests)

covered.species.spring <- c("Black kite","Steppe buzzard","Eurasian sparrowhawk",
                            "Honey buzzard","Marsh harrier","Levant sparrowhawk",
                            "Booted eagle","Short-toed eagle","Lesser spotted eagle",
                            "Greater spotted eagle","Hen harrier","Pallid harrier",
                            "Imperial eagle","Steppe eagle","Montagu's harrier") #species covered in spring

covered.species.autumn <- c("Montagu's harrier","Pallid harrier","Marsh harrier",
                            "Black kite","Honey buzzard","Steppe buzzard","Booted eagle",
                            "Short-toed eagle","Lesser spotted eagle","Osprey") #species covered in autumn

covered.species.combined <- c("Montagu's harrier","Pallid harrier","Marsh harrier",
                              "Black kite","Honey buzzard","Booted eagle",
                              "Short-toed eagle","Lesser spotted eagle",
                              "Steppe buzzard") #species covered in both seasons

################################################################################
                        ### SUMMARIZE ALL DATA ###
################################################################################

### READ DATA ###
comb_daily <- readRDS("comb_daily.RDS") # DAILY totals both seasons combined
comb_fin   <- readRDS("comb_fin.RDS")   # daily & annual & seasons combined

                          #####################
                          #### PER SPECIES ####
                          #####################

### Summarize annual statistics per species ###
annual.quantiles <- comb_fin %>%
  group_by(englishname,latinname,speciesgroup,season,year)%>%
  reframe(yr.tot = unique(yrtot),
          yr.daymax = unique(yr.daymax),
          q01 = min(date[cumprop >= .01]),
          q05 = min(date[cumprop >= .05]),
          q25 = min(date[cumprop >= .25]),
          q50 = min(date[cumprop >= .50]),
          q75 = min(date[cumprop >= .75]),
          q95 = min(date[cumprop >= .95]),
          q05.doy = yday(q05),
          q25.doy = yday(q25),
          q50.doy = yday(q50),
          q75.doy = yday(q75),
          q95.doy = yday(q95),
          core.50.duration = as.numeric(difftime(q75,q25,units="days"))+1,
          main.90.duration = as.numeric(difftime(q95,q05,units="days"))+1,
          tail.20.duration = as.numeric(difftime(q95,q75,units="days"))+1,
          peak.date = date[which.max(daytot.fin)],
          peak_percentage = (yr.daymax / yr.tot)*100) %>%
  distinct()

### Summarize seasonal statistics per species ###
seasonal.quantiles <- annual.quantiles %>%
  group_by(englishname,latinname,speciesgroup,season) %>%
  reframe(n.yrs       = length(year[yr.tot > 0]),
          sd.yrtot    = round(sd(yr.tot[yr.tot > 0])),
          mean.yrtot  = round(mean(yr.tot[yr.tot > 0])),
          year.record = round(max(yr.tot,na.omit=TRUE)),
          year.lowest = round(min(yr.tot,na.omit=TRUE)),
          day.record  = round(max(yr.daymax,na.omit=TRUE)),
          mean.daymax = round(mean(yr.daymax,na.omit=TRUE)),
          sd.daymax   = round(sd(yr.daymax[yr.tot > 0])),
          mean.q01    = round(mean(yday(q01[yr.tot > 0]))),
          mean.q05    = round(mean(yday(q05[yr.tot > 0]))),
          mean.q25    = round(mean(yday(q25[yr.tot > 0]))),
          mean.q50    = round(mean(yday(q50[yr.tot > 0]))),
          mean.q75    = round(mean(yday(q75[yr.tot > 0]))),
          mean.q95    = round(mean(yday(q95[yr.tot > 0]))),
          mean.peak   = round(mean(yday(peak.date[yr.tot > 0]),na.omit=TRUE)),
          sd.q05      = round(sd(yday(q05[yr.tot > 0]))),
          sd.q25      = round(sd(yday(q25[yr.tot > 0]))),
          sd.q50      = round(sd(yday(q50[yr.tot > 0]))),
          sd.q75      = round(sd(yday(q75[yr.tot > 0]))),
          sd.q95      = round(sd(yday(q95[yr.tot > 0]))),
          sd.peak     = round(sd(yday(peak.date[yr.tot > 0]))),
          mean.main.90.duration = round(mean(main.90.duration)),
          sd.main.90.duration   = round(sd(main.90.duration[yr.tot > 0])),
          mean.core.50.duration = round(mean(core.50.duration)),
          sd.core.50.duration   = round(sd(core.50.duration[yr.tot > 0])),
          mean.tail.20.duration = round(mean(tail.20.duration)),
          sd.tail.20.duration   = round(sd(tail.20.duration[yr.tot > 0])),
          mean.peak.percentage  = round(mean(peak_percentage)),
          sd.peak.percentage    = round(sd(peak_percentage)),
          min.peak.percentage   = min(peak_percentage),
          max.peak.percentage   = max(peak_percentage),
          dummy.date.q05        = if_else(season == "spring",
                                    as.Date(parse_date_time(x = paste(9999, mean.q05), orders = "yj")),
                                    NA_Date_),
          dummy.date.q50        = if_else(season == "spring",
                                   as.Date(parse_date_time(x = paste(9999, mean.q50), orders = "yj")),
                                   NA_Date_), 
          dummy.date.q95        = if_else(season == "spring",
                                    as.Date(parse_date_time(x = paste(9999, mean.q95), orders = "yj")),
                                    NA_Date_),
          dummy.date.peak       = if_else(season == "spring",
                                   as.Date(parse_date_time(x = paste(9999, mean.peak), orders = "yj")),
                                   NA_Date_)
          ) %>%
  distinct()

### Save seasonal quantiles
saveRDS(seasonal.quantiles, "seasonal.quantiles.RDS") 


                        ##########################
                        #### SPECIES COMBINED ####
                        ##########################

### Summarize annual statistics all species ### - ADD FILTER RAPTORS ONLY
mean.annual.allspecs <- comb_fin %>%
  filter(!(englishname %in% non.raptors)) %>%
  group_by(englishname, latinname, season, year) %>%
  summarise(yr.tot = unique(yrtot), .groups = 'drop') %>%
  group_by(season, year) %>%
  summarise(sum.yr.tot = sum(yr.tot), .groups = 'drop')

### Summarize seasonal statistics all species ###
mean.season.allspecs <- mean.annual.allspecs %>%
  group_by(season) %>%
  summarise(mean.season.tot = mean(sum.yr.tot),
            sd.mean.season.tot = round(sd(sum.yr.tot)), .groups = 'drop')


################################################################################
                          ### CREATE TABLES ###
################################################################################

###########################################
### TABLE 1 - ALL SPECIES - SPRING ONLY ###
###########################################

# Select variables wanted in the final table
table.1 <- seasonal.quantiles %>%  
  filter(season=="spring") %>%
  subset(select = c("englishname","latinname","mean.yrtot","sd.yrtot","dummy.date.q50","sd.q50","dummy.date.peak","sd.peak","mean.daymax","sd.daymax","n.yrs"))

#Write table
write.xlsx(table.1, "Output tables and figures/Tables/table1.xlsx")

################################################################
### TABLE 2 - SEASONAL MEANS - BOTH SEASONS - INCLUDING TESTS###
################################################################

### Select seasonal means
count.avg.spring <- seasonal.quantiles %>%
  filter(englishname %in% test.species, season == "spring") %>%
  select(englishname,latinname,spring.mean.yrtot=mean.yrtot,spring.sd.yrtot=sd.yrtot)

### Select autumn core period
count.avg.autumn <- seasonal.quantiles %>%
  filter(englishname %in% test.species, season == "autumn") %>%
  select(englishname,latinname,autumn.mean.yrtot=mean.yrtot,autumn.sd.yrtot=sd.yrtot)

count.comb.avg <- merge(count.avg.spring,count.avg.autumn,by=c('englishname','latinname'))

write.xlsx(count.comb.avg, "Output tables and figures/Tables/table2.xlsx")

### Perform tests for difference in mean total per species
count_t_test_results <- list()

for(i in test.species) {
  # Filter data 
  count.spring_data <- annual.quantiles %>% filter(englishname == i & season == "spring") %>% select(yr.tot)
  count.autumn_data <- annual.quantiles %>% filter(englishname == i & season == "autumn") %>% select(yr.tot)
  
  # Perform t-test (unpaired)
  count.t_test_result <- t.test(count.spring_data$yr.tot, count.autumn_data$yr.tot, alternative = "two.sided", paired = FALSE)
  
  # Store result 
  count_t_test_results[[i]] <- count.t_test_result
}

print(count_t_test_results)

### test for season totals
spring.data.yrtot <- mean.annual.allspecs %>% filter(season == "spring") %>% select(sum.yr.tot)
autumn_data.yrtot <- mean.annual.allspecs %>% filter(season == "autumn") %>% select(sum.yr.tot)

t.test(spring.data.yrtot$sum.yr.tot, autumn_data.yrtot$sum.yr.tot, alternative = "two.sided", paired = FALSE) # Performing the t-test (unpaired)


################################################################
### TABLE 3 - MAIN 90 - BOTH SEASONS - INCLUDING TESTS ###
################################################################

### Select spring core period
core.spring.90 <- seasonal.quantiles %>%
  filter(season=="spring" & englishname %in% test.species) %>%
  select(englishname,latinname,spring.core.90=mean.main.90.duration,spring.sd.core=sd.main.90.duration) %>%
  distinct()

### Select autumn core period
core.autumn.90 <- seasonal.quantiles %>%
  filter(season=="autumn" & englishname %in% test.species) %>%
  select(englishname,latinname,autumn.core.90=mean.main.90.duration,autumn.sd.core=sd.main.90.duration) %>%
  distinct()

comb.core.periods.90 <- merge(core.spring.90,core.autumn.90,by=c('englishname','latinname'))

write.xlsx(comb.core.periods.90, "Output tables and figures/Tables/table3.xlsx")

### Perform test for difference in 90 duration between seasons per species ###
t_test_results_quantiles <- list()

for(i in test.species) {
  # Filter data 
  spring_data_quantile <- annual.quantiles %>% filter(englishname == i & season=="spring") %>% select(main.90.duration)
  autumn_data_quantile <- annual.quantiles %>% filter(englishname == i & season=="autumn") %>% select(main.90.duration)
  
  # Perform t-test
  t_test_result_quant <- t.test(spring_data_quantile$main.90.duration, autumn_data_quantile$main.90.duration, alternative = "two.sided", paired = FALSE)
  
  # Store result 
  t_test_results_quantiles[[i]] <- t_test_result_quant
}

print(t_test_results_quantiles)

################################################################
### TABLE 4 - CORE 50 - BOTH SEASONS - INCLUDING TESTS ###
################################################################

### Select spring core period
core.spring.50 <- seasonal.quantiles %>%
  filter(season=="spring" & englishname %in% test.species) %>%
  select(englishname,latinname,spring.core.50=mean.core.50.duration,spring.sd.core=sd.core.50.duration) %>%
  distinct()

### Select autumn core period
core.autumn.50 <- seasonal.quantiles %>%
  filter(season=="autumn" & englishname %in% test.species) %>%
  select(englishname,latinname,autumn.core.50=mean.core.50.duration,autumn.sd.core=sd.core.50.duration) %>%
  distinct()

comb.core.periods.50 <- merge(core.spring.50,core.autumn.50,by=c('englishname','latinname'))

write.xlsx(comb.core.periods.50, "Output tables and figures/Tables/table4.xlsx")

### Perform test for difference in 90 duration between seasons per species ###
t_test_results_quantiles <- list()

for(i in test.species) {
  # Filter data 
  spring_data_quantile <- annual.quantiles %>% filter(englishname == i & season=="spring") %>% select(core.50.duration)
  autumn_data_quantile <- annual.quantiles %>% filter(englishname == i & season=="autumn") %>% select(core.50.duration)
  
  # Perform t-test
  t_test_result_quant <- t.test(spring_data_quantile$core.50.duration, autumn_data_quantile$core.50.duration, alternative = "two.sided", paired = FALSE)
  
  # Store result 
  t_test_results_quantiles[[i]] <- t_test_result_quant
}

print(t_test_results_quantiles)


##################################################################
### TABLE 5 - PEAK PERCENTAGE - BOTH SEASONS - INCLUDING TESTS ###
##################################################################

### Select spring peak percentage
percent.spring.peak <- seasonal.quantiles %>%
  filter(season=="spring" & englishname %in% test.species) %>%
  select(englishname,latinname,spring.peak.percent=mean.peak.percentage,spring.sd.peak=sd.peak.percentage) %>%
  distinct()

### Select autumn peak percentage
percent.autumn.peak <- seasonal.quantiles %>%
  filter(season=="autumn" & englishname %in% test.species) %>%
  select(englishname,latinname,autumn.peak.percent=mean.peak.percentage,autumn.sd.peak=sd.peak.percentage) %>%
  distinct()

comb.peak.percentage <- merge(percent.spring.peak,percent.autumn.peak,by=c('englishname','latinname'))

write.xlsx(comb.peak.percentage, "Output tables and figures/Tables/table5.xlsx")

### Perform test for difference in 90 duration between seasons per species ###
t_test_results_quantiles <- list()

for(i in test.species) {
  # Filter data 
  spring_data_quantile <- annual.quantiles %>% filter(englishname == i & season=="spring") %>% select(peak_percentage)
  autumn_data_quantile <- annual.quantiles %>% filter(englishname == i & season=="autumn") %>% select(peak_percentage)
  
  # Perform t-test
  t_test_result_quant <- t.test(spring_data_quantile$peak_percentage, autumn_data_quantile$peak_percentage, alternative = "two.sided", paired = FALSE)
  
  # Store result 
  t_test_results_quantiles[[i]] <- t_test_result_quant
}

print(t_test_results_quantiles)
