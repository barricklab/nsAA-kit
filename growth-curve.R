library(readxl)
library(dplyr)
library(magrittr)
library(tidyr)
library(ggplot2)

theme_set(theme_bw(base_size = 16))
line_thickness = 1.5
# theme_update(legend.background = element_rect(colour =alpha('white',0),fill=alpha('white', 0)), 
#              legend.text = element_text(size = 10), 
#              legend.title = element_text(size=10),
#              legend.key.size = unit(0.5, "cm"),
#              panel.grid.minor = element_blank(), 
#              panel.grid.major=element_blank(),
#              panel.border=element_rect(fill=NA), 
#              legend.key=element_rect(color=NA, fill=NA))

parseIcontrolOutput <- function(data, lookupTable, numMeasurementTypes = 1, numWells = 96) {
  rowsPerMeasure <- numWells + 4
  
  data %T>%
  {
    colnames(.) <- c("well", paste(data[2,][2:ncol(data)], data[3,][2:ncol(data)], sep ="|"))
  } %>%
  {
    measurementLocations <- seq_along(rownames(data))[seq(1, (rowsPerMeasure * numMeasurementTypes), rowsPerMeasure)]
    mutate(., measurementType = rep(data[[1]][measurementLocations], each = rowsPerMeasure)) %>%
    {
      metadataRows <- unlist(lapply(measurementLocations,
                                    function(i) {
                                      seq(i, i+3)
                                    }))
      slice(., -c(metadataRows))
    }
  } %>%
    gather(key, value, -well, -measurementType) %>%
    separate(key, into = c("cycle", "time"), sep ="\\|") %>%
    mutate(time = (as.numeric(time))/3600) %>%
    inner_join(lookupTable, by = "well") %>%
    mutate(well = as.factor(well), 
           measurementType = as.factor(measurementType), 
           cycle = as.factor(cycle), 
           medium = as.factor(medium), 
           strain = factor(strain), 
           replicate = as.factor(replicate),
           supplements = as.factor(supplements)) %>%
    arrange(well)
}
masterTbl <- tbl_df(data.frame())
rawDataTbl <- read_excel(file.path("~/box/c321-growthcurve-f500-101115.xlsx"), col_names = FALSE, skip = 0) %>% slice(1:100)
rawKeyTbl <- read_excel(file.path("~/github/research/raw-data/10-11-15 growth curve key.xlsx"))
masterTbl <- parseIcontrolOutput(rawDataTbl, rawKeyTbl) %>% bind_rows(masterTbl) %>% mutate(value = as.numeric(value))

levels(masterTbl$measurementType) <- c("OD[600]")

masterTbl$strain <- factor(masterTbl$strain,
                      levels = c('MG1655','MJH116','OPT97','OPT99','OPT100','MJH184','OPT103','OPT104','OPT105', 'blank'))
masterTbl$supplements <- factor(masterTbl$supplements,
                           levels = c('Zeo', 'none', 'Zeo Gent IodoY'))
expTbl <- masterTbl %>% filter(strain != 'blank') %>% 
  group_by(time, strain, medium, supplements) %>%
  summarise(n = n(), meanValue = mean(value), sd = sd(value)) %>% 
  mutate(se = sd/sqrt(n), LCI = meanValue+qt(0.025,df=(n-1))*se, UCI = meanValue+qt(0.975, df=(n-1))*se) #%>% 
  #filter(medium == "LB", supplements == "Zeo Gent IodoY")


ggplot(expTbl, aes(x=time, y=meanValue, color=strain)) + 
  geom_line() +
  facet_grid(medium ~ supplements) + #, labeller=label_parsed) +
  ggtitle("24 Hours Growth of OPT Lines and Ancestors by Medium") +
  scale_y_log10(limits = c(0.065,1.0), breaks=c(0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.825)) +
  scale_x_continuous(limits = c(0, 24.5)) +
  ylab(bquote(log[10] ~ OD[600])) +
  xlab("Time (h)") +
  geom_errorbar(aes(x=time, ymin=LCI, ymax=UCI), size=.1, width=.25)

# expTbl2 <- expTbl %>% filter(medium == "M9") %>% group_by(time, strain, supplements)
# 
# ggplot(aes(x=time, y=meanValue)) + 
#   facet_grid(. ~ supplements) +
#   geom_line(size=1.125, color=strain) +
#   geom_errorbar(aes(x=time, ymin=LCI, ymax=UCI), size=.25, width=.1) +
#   scale_y_log10(limits = c(0.065,0.9), breaks=c(0.125, 0.25, 0.375, 0.5, 0.625, 0.75)) +
#   scale_x_continuous(limits = c(0, 17.5)) +
#   ggtitle("17 Hours Growth of OPT Lines and Ancestors in M9") +
#   ylab(bquote(log[10] ~ OD[600])) +
#   xlab("Time (h)")
  
# expTbl %>% filter(supplements != "Zeo", strain != "OPT105") %>% group_by(time, strain) %>%
# ggplot(aes(x=time, y=meanValue, color=strain)) + 
#   geom_line(size=1.125) +
#   geom_errorbar(aes(x=time, ymin=LCI, ymax=UCI), size=.25, width=.1) +
#   scale_y_log10(limits = c(0.065,0.9), breaks=c(0.125, 0.25, 0.375, 0.5, 0.625, 0.75)) +
#   scale_x_continuous(limits = c(0, 17.5)) +
#   ggtitle("17 Hours Growth of OPT Lines and Ancestors in M9") +
#   ylab(bquote(log[10] ~ OD[600])) +
#   xlab("Time (h)")
# 
