##################
### Read packages 
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("patchwork","ggforce","data.table","ggplot2","lubridate","maptools","RColorBrewer","tidyverse",
              "ggforce","ggridges","cowplot","openxlsx","showtext","sysfonts","gridExtra","summarytools","dplyr")
ipak(packages)

################################################################################
                            ### READ DATA ###
################################################################################
comb_daily <- readRDS("comb_daily.RDS") # daily totals 
seasonal.quantiles <- readRDS("seasonal.quantiles.RDS") # seasonal quantiles

################################################################################
        #### VISUALIZING DAILY PASSAGE AND QUANTILES SPRING ####
################################################################################

### Adjust species names - and add early species indicator (*)
name_mapping <- c("Black kite" = "Black Kite", #names to replace
                  "Steppe buzzard" = "Steppe Buzzard",
                  "Eurasian sparrowhawk" = "Eurasian Sparrowhawk",
                  "Honey buzzard" = "Honey Buzzard",
                  "Marsh harrier" = "Marsh Harrier*",
                  "Levant sparrowhawk" = "Levant Sparrowhawk",
                  "Booted eagle" = "Booted Eagle",
                  "Short-toed eagle" = "Short-toed Eagle",
                  "Lesser spotted eagle" = "Lesser Spotted Eagle",
                  "Greater spotted eagle" = "Greater Spotted Eagle",
                  "Hen harrier" = "Hen Harrier*",
                  "Pallid harrier" = "Pallid Harrier",
                  "Imperial eagle" = "Imperial Eagle*",
                  "Steppe eagle" = "Steppe Eagle*",
                  "Montagu's harrier" = "Montagu's Harrier")

### Select only species fully covered in spring
d.fin <- comb_daily %>%
  filter(season=="spring", englishname %in% covered.species.spring) 

### Rename species with capitals
d.fin <- d.fin %>% 
  mutate(englishname = ifelse(englishname %in% names(name_mapping), name_mapping[englishname], englishname))

### Create dummy dates for plotting time series in multiple years on same axis
# Set class to date
d.fin$date <- date(d.fin$date)
# Create daily dummy dates accounting for leap year 2020
d.fin$dummy.date <- as.Date(ifelse(d.fin$year %in% c(2019,2022),as.Date(paste("9999",format(d.fin$date,"%m-%d"),sep="-"),tz='UTC'),as.Date(paste("9999",format(d.fin$date+1,"%m-%d"),sep="-"),tz='UTC')),origin="1970-01-01") 
# Add month value
d.fin$mth <- as.factor(month(d.fin$date))

### Add max. day record per species
d.fin <- d.fin %>%
  group_by(englishname) %>%
  mutate(day.record = max(max(daytot.fin), na.rm = TRUE)) %>%
  ungroup()

### Add 0 values for all missing dates per species per year
d.fin <- d.fin %>%
  group_by(englishname, year) %>%
  complete(dummy.date = {
    current_year <- as.integer(cur_group()$year)
    if (current_year == 2019) {
      seq(as.Date("9999-03-21"), as.Date("9999-05-31"), by = "1 day")
    } else if (current_year %in% c(2020, 2022)) {
      seq(as.Date("9999-03-01"), as.Date("9999-05-26"), by = "1 day")
    }
  }) %>%
  replace_na(list(daytot.fin = 0)) %>%
  mutate(daytot.fin = ifelse(daytot.fin < 0, 0, daytot.fin)) %>%
  ungroup()

### Summarize daily statistics per species
daily.means.fig <- d.fin %>%
  group_by(englishname,dummy.date) %>%
  reframe(daily.mean = round(mean(daytot.fin)),
          daily.sd   = round(sd(daytot.fin)),
          daily.min  = min(daytot.fin),
          daily.max  = max(daytot.fin)) 

### Select desired species for figure
seasonal.quantiles.figure <- seasonal.quantiles %>% 
  filter(season=="spring", englishname %in% covered.species.spring) %>%
  mutate(englishname = ifelse(englishname %in% names(name_mapping), name_mapping[englishname], englishname))

### Adjust order of species according to median passage date
# Extract desired order
seasonal.quantiles.figure <- seasonal.quantiles.figure[order(seasonal.quantiles.figure$mean.q50),]
# Select species names
order.graphs <- unique(seasonal.quantiles.figure$englishname)
# Adjust order in all relevant dataframes
daily.means.fig <- daily.means.fig %>%
  mutate(englishname = factor(englishname, levels = order.graphs)) %>%
  arrange(englishname)
d.fin <- d.fin %>%
  mutate(englishname = factor(englishname, levels = order.graphs)) %>%
  arrange(englishname)
seasonal.quantiles.figure <- seasonal.quantiles.figure %>%
  mutate(englishname = factor(englishname, levels = order.graphs)) %>%
  arrange(englishname)

### PLOT FIGURE (PAGE 1) ###
p.wrap.1 <- ggplot() +
  geom_rect(data=d.fin,aes(xmin=as.Date("9999-03-01")-0.5,xmax=as.Date("9999-03-31")+0.5,ymin=0,ymax=day.record+6),fill="grey80",alpha=.8)+
  geom_rect(data=d.fin,aes(xmin=as.Date("9999-04-01")-0.5,xmax=as.Date("9999-04-30")+0.5,ymin=0,ymax=day.record+6),fill="grey90",alpha=.8)+
  geom_rect(data=d.fin,aes(xmin=as.Date("9999-05-01")-0.5,xmax=as.Date("9999-05-25")+0.5,ymin=0,ymax=day.record+6),fill="grey80",alpha=.8)+
  geom_ribbon(data=daily.means.fig,aes(x=dummy.date, ymin = daily.min, ymax = daily.max), fill='steelblue1',col='steelblue3',alpha=.4,linewidth=.15) +
  geom_line(data=daily.means.fig,aes(x=dummy.date, y=daily.mean), color = "steelblue4", linewidth = 1) +
  geom_vline(data=seasonal.quantiles.figure,aes(xintercept = dummy.date.q05), linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_vline(data=seasonal.quantiles.figure,aes(xintercept = dummy.date.q95), linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_vline(data=seasonal.quantiles.figure,aes(xintercept = dummy.date.q50), color = "black", linewidth = 1) +
  scale_x_date(limits=c(as.Date("9999-03-01")-0.5,as.Date("9999-05-25")+0.5),date_breaks="12 days",date_labels="%d-%m",expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  ylab("Average daily count") +
  xlab("Date") +
  theme_bw() +
  facet_wrap_paginate(~englishname, nrow = 4, ncol = 2, scales = "free_y", page = 1) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.text.x	= element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=24,face='bold',margin = margin(t = 15)),
        axis.title.y = element_text(size=24,face='bold',margin = margin(t = 15)),
        plot.title = element_blank(),
        panel.spacing = unit(1.5, "lines"),
        strip.text = element_text(size = 20, face = "bold"),
        plot.margin = margin(12, 15, 12, 12))

### PLOT FIGURE (PAGE 2) ###
p.wrap.2 <- ggplot() +
  geom_rect(data=d.fin,aes(xmin=as.Date("9999-03-01")-0.5,xmax=as.Date("9999-03-31")+0.5,ymin=0,ymax=day.record+6),fill="grey80",alpha=.8)+
  geom_rect(data=d.fin,aes(xmin=as.Date("9999-04-01")-0.5,xmax=as.Date("9999-04-30")+0.5,ymin=0,ymax=day.record+6),fill="grey90",alpha=.8)+
  geom_rect(data=d.fin,aes(xmin=as.Date("9999-05-01")-0.5,xmax=as.Date("9999-05-25")+0.5,ymin=0,ymax=day.record+6),fill="grey80",alpha=.8)+
  geom_ribbon(data=daily.means.fig,aes(x=dummy.date, ymin = daily.min, ymax = daily.max), fill='steelblue1',col='steelblue3',alpha=.4,linewidth=.15) +
  geom_line(data=daily.means.fig,aes(x=dummy.date, y=daily.mean), color = "steelblue4", linewidth = 1) +
  geom_vline(data=seasonal.quantiles.figure,aes(xintercept = dummy.date.q05), linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_vline(data=seasonal.quantiles.figure,aes(xintercept = dummy.date.q95), linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_vline(data=seasonal.quantiles.figure,aes(xintercept = dummy.date.q50), color = "black", linewidth = 1) +
  scale_x_date(limits=c(as.Date("9999-03-01")-0.5,as.Date("9999-05-25")+0.5),date_breaks="12 days",date_labels="%d-%m",expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  ylab("Average daily count") +
  xlab("Date") +
  theme_bw() +
  facet_wrap_paginate(~englishname, ncol = 2, nrow = 4, scales = "free_y", page = 2) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.text.x	= element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=24,face='bold',margin = margin(t = 15)),
        axis.title.y = element_text(size=24,face='bold',margin = margin(t = 15)),
        plot.title = element_blank(),
        panel.spacing = unit(1.5, "lines"),
        strip.text = element_text(size = 20, face = "bold"),
        plot.margin = margin(12, 15, 12, 12))

### SAVE PLOTS FIGURE ###
# (PAGE 1)
ggsave(plot=p.wrap.1,filename=paste("Output tables and figures/graphs wrapped - p1.png"),dpi=360,width=16.5,height=19)
# (PAGE 2)
ggsave(plot=p.wrap.2,filename=paste("Output tables and figures/graphs wrapped - p2.png"),dpi=360,width=16.5,height=19)

################################################################################
        #### VISUALIZING QUANTILES SPRING VS AUTUMN ####
################################################################################

#############################
## AUTUMN VS SPRING MEDIAN ##
#############################

### Select median dates of species for autumn 
autumn.data.median <- seasonal.quantiles %>%
  filter(englishname %in% covered.species.combined, season == "autumn") %>%
  select(englishname, mean.q50)

### Select median dates of species for spring 
spring.data.median <- seasonal.quantiles %>%
  filter(englishname %in% covered.species.combined, season == "spring") %>%
  select(englishname, mean.q50)

### Merge autumn and spring data
merged.data.median <- inner_join(autumn.data.median, spring.data.median, by = "englishname", suffix = c("_autumn", "_spring"))

### Fit linear model
m.median      <- summary(lm(mean.q50_spring ~ mean.q50_autumn, merged.data.median))
# Save confidence interval of LM to use in plot
limits.median <- as.data.frame(predict(lm(mean.q50_spring ~ mean.q50_autumn, merged.data.median), 
                     se.fit = TRUE, interval = "confidence")$fit)

### Plot median dates spring vs. autumn
p_median_spring_vs_autumn <- ggplot(merged.data.median) +
  geom_smooth(aes(x = mean.q50_autumn, y = mean.q50_spring), method = "lm", color = "black", se = FALSE, linewidth = 1) +
  geom_line(aes(x = mean.q50_autumn,y = limits.median$upr), linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_line(aes(x = mean.q50_autumn,y = limits.median$lwr), linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_point(aes(x = mean.q50_autumn, y = mean.q50_spring, color = englishname),size = 2.5) +
  scale_color_brewer(palette="Paired") +
  theme_bw() +
  labs(x = "Median passage date in autumn (julian day)", y = "Median passage date in spring (julian day)", color = "Species") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

###########################
## AUTUMN VS SPRING PEAK ##
###########################

### Select peak percentages of species for autumn 
autumn.data.peak <- seasonal.quantiles %>%
  filter(englishname %in% covered.species.combined, season == "autumn") %>%
  select(englishname, mean.peak.percentage)

### Select peak percentages of species for spring 
spring.data.peak <- seasonal.quantiles %>%
  filter(englishname %in% covered.species.combined, season == "spring") %>%
  select(englishname, mean.peak.percentage)

### Merge autumn and spring data
merged.data.peak <- inner_join(autumn.data.peak, spring.data.peak, by = "englishname", suffix = c("_autumn", "_spring"))

### Fit linear model
m.peak      <- summary(lm(mean.peak.percentage_spring ~ mean.peak.percentage_autumn, merged.data.peak))
# Save confidence interval of LM to use in plot
limits.peak <- as.data.frame(predict(lm(mean.peak.percentage_spring ~ mean.peak.percentage_autumn, merged.data.peak), 
                se.fit = TRUE, interval = "confidence")$fit)

### Plot peak percentage spring vs. autumn
p_peak_spring_vs_autumn <- ggplot(merged.data.peak, aes(x = mean.peak.percentage_autumn)) +
  geom_ribbon(aes(ymin=min(limits.peak$lwr), ymax=pmin(mean.peak.percentage_autumn,mean.peak.percentage_autumn)), fill = "grey70", alpha=0.5) +
  geom_smooth(aes(y=mean.peak.percentage_spring), method = "lm", color = "black", se = FALSE, linewidth = 1) +
  geom_line(aes(y = limits.peak$upr), linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_line(aes(y = limits.peak$lwr), linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_point(aes(y = mean.peak.percentage_spring, color = englishname),size = 2.5) +
  scale_color_brewer(palette="Paired") +
  theme_bw() +
  labs(x = "Percentage of passage on peak days in autumn", y = "Percentage of passage on peak days in spring", color = "Species") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

##########################
## AUTUMN VS SPRING 90% ##
##########################

### Select main passage (90) period of species for autumn 
autumn.90 <- seasonal.quantiles %>%
  filter(englishname %in% covered.species.combined, season == "autumn") %>%
  select(englishname, mean.main.90.duration)

### Select main passage (90) period of species for spring 
spring.90 <- seasonal.quantiles %>%
  filter(englishname %in% covered.species.combined, season == "spring") %>%
  select(englishname, mean.main.90.duration)

### Merge autumn and spring data
merged.data.90 <- inner_join(autumn.90, spring.90, by = "englishname", suffix = c("_autumn", "_spring"))

### Fit linear model
m.90 <- summary(lm(mean.main.90.duration_spring ~ mean.main.90.duration_autumn, merged.data.90))
# Save confidence interval of LM to use in plot
limits.90 <- as.data.frame(predict(lm(mean.main.90.duration_spring ~ mean.main.90.duration_autumn, merged.data.90), 
                     se.fit = TRUE, interval = "confidence")$fit)

### Plot main passage (90) period spring vs. autumn
p_90_spring_vs_autumn <- ggplot(merged.data.90, aes(x=mean.main.90.duration_autumn)) +
  geom_ribbon(aes(ymin=min(limits.90$lwr), ymax=pmin(mean.main.90.duration_autumn,mean.main.90.duration_autumn)), fill = "grey70", alpha=0.5) +
  #geom_abline(slope=1, intercept = 0, color="grey80", linewidth=0.5) +
  geom_smooth(aes(y=mean.main.90.duration_spring), method = "lm", color = "black", se = FALSE, linewidth = 1) +
  geom_line(aes(y = limits.90$upr), linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_line(aes(y = limits.90$lwr), linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_point(aes(y = mean.main.90.duration_spring, color = englishname),size = 2.5) +
  scale_color_brewer(palette="Paired") +
  coord_cartesian(clip = 'off') + 
  theme_bw() +
  labs(x = "Central 90% migration duration in autumn (days)", y = "Central 90% migration duration in spring (days)", color = "Species") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

#############################
 ## AUTUMN VS SPRING 50% ##
#############################

### Select core passage (50) period of species for autumn 
autumn.50 <- seasonal.quantiles %>%
  filter(englishname %in% covered.species.combined, season == "autumn") %>%
  select(englishname, mean.core.50.duration)

### Select core passage (50) period of species for spring 
spring.50 <- seasonal.quantiles %>%
  filter(englishname %in% covered.species.combined, season == "spring") %>%
  select(englishname, mean.core.50.duration)

### Merge autumn and spring data
merged.data.50 <- inner_join(autumn.50, spring.50, by = "englishname", suffix = c("_autumn", "_spring"))

### Fit linear model
m.50      <- summary(lm(mean.core.50.duration_spring ~ mean.core.50.duration_autumn, merged.data.50))
# Save confidence interval of LM to use in plot
limits.50 <- as.data.frame(predict(lm(mean.core.50.duration_spring ~ mean.core.50.duration_autumn, merged.data.50), 
                   se.fit = TRUE, interval = "confidence")$fit)

### Plot core passage (50) period spring vs. autumn
p_50_spring_vs_autumn <- ggplot(merged.data.50, aes(x=mean.core.50.duration_autumn)) +
  geom_ribbon(aes(ymin=min(limits.50$lwr), ymax=pmin(mean.core.50.duration_autumn,mean.core.50.duration_autumn)), fill = "grey70", alpha=0.5) +
  #geom_abline(slope=1, intercept = 0, color="grey80", linewidth=0.5) +
  geom_smooth(aes(y=mean.core.50.duration_spring), method = "lm", color = "black", se = FALSE, linewidth = 1) +
  geom_line(aes(y = limits.50$upr), linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_line(aes(y = limits.50$lwr), linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_point(aes(y = mean.core.50.duration_spring, color = englishname),size = 2.5) +
  scale_color_brewer(palette="Paired") +
  coord_cartesian(clip = 'off') + 
  theme_bw() +
  labs(x = "Central 50% migration duration in autumn (days)", y = "Central 50% migration duration in spring (days)", color = "Species") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

###########################
### COMBINE & SAVE PLOT ###
###########################

### Create function to extract species legend
get_legend <- function(p) {
  g <- ggplotGrob(p)
  legend <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
  return(legend)
}

### Extract species legend from relevant plot
legend <- get_legend(p_median_spring_vs_autumn + theme(legend.position = "top"))

### Combine all plots, include the legend and add plot labels
combined_plot <- plot_grid(
  legend, 
  plot_grid(
    p_median_spring_vs_autumn, 
    p_peak_spring_vs_autumn, 
    p_90_spring_vs_autumn, 
    p_50_spring_vs_autumn, 
    ncol = 2,
    labels = "AUTO",
    label_x = 0, label_y = 1,
    hjust = -4.5, vjust = 2.5
  ), 
  ncol = 1, 
  rel_heights = c(0.1, 1))

### Save final plot
ggsave(plot = combined_plot, filename=paste0("Output tables and figures/quantiles.png"),dpi=360,width=10,height=8) 


################################################################################
                   #### CALCULATING & VISUALIZING PEAK DAYS ####
################################################################################

### Define species selections (if not done previously in Processing_fin file)
# Unidentified species groups
uid.peak <- c("Common/Lesser kestrel","buzzard spec","Stork spec", "Falcon spec",
              "raptor spec", "Hobby/Red-footed falcon", "harrier spec","honey buzzard spec")

# Covered species
covered.species.spring.peak <- c("Black kite","Steppe buzzard","Eurasian sparrowhawk",
                                 "Honey buzzard","Marsh harrier","Levant sparrowhawk",
                                 "Booted eagle","Short-toed eagle","Lesser spotted eagle",
                                 "Greater spotted eagle","Hen harrier","Pallid harrier",
                                 "Imperial eagle","Steppe eagle","Montagu's harrier", 
                                 "eagle spec", "Ringtail harriers", "raptor spec (medium)", "hawk spec")

### Filter and order species in dataframe (relevant for order of species within final plot)
mean.peak.days.spring <- seasonal.quantiles %>%
  filter(!(englishname %in% uid.peak),englishname %in% covered.species.spring.peak,season=="spring") %>%
  mutate(englishname = factor(englishname, levels = c("hawk spec", "Levant sparrowhawk", 
                                                      "Eurasian sparrowhawk","eagle spec", 
                                                      "Short-toed eagle", "Steppe eagle", 
                                                      "Greater spotted eagle", "Imperial eagle", 
                                                      "Lesser spotted eagle",  "raptor spec (medium)", 
                                                      "Honey buzzard", "Black kite", "Steppe buzzard", 
                                                      "Booted eagle",
                                                      "Ringtail harriers", "Montagu's harrier", 
                                                      "Pallid harrier", "Marsh harrier", 
                                                      "Hen harrier")))

### Manually adjust speciesgroup of Booted Eagle
mean.peak.days.spring <- mean.peak.days.spring %>%
  mutate(speciesgroup = ifelse(englishname == "Booted eagle", "Buzzards & Kites", speciesgroup))

### Plot the peak percentages of species in spring
p_peak_spring <- ggplot(mean.peak.days.spring, aes(x=mean.peak.percentage, y=englishname, color = speciesgroup)) + #ORDER ON SPECIESGROUP AND PHENOLOGY
  geom_point(size=2.25) +
  geom_errorbar(aes(xmin=min.peak.percentage, xmax=max.peak.percentage), width=0.4, linewidth=1) + 
  scale_color_brewer(palette = "Set2") +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) + 
  theme_bw() +
  labs(x = "% seasonal total on peak days") + 
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.margin = margin(10, 15, 10, 12))

### Save final plot
ggsave(plot=p_peak_spring, filename=paste0("Output tables and figures/peak days spring.png"),dpi=360,width=7,height=5)
