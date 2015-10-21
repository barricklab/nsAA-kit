library(readxl)
library(dplyr)
library(magrittr)
library(tidyr)
library(ggplot2)

#USE R VERSION 3.1.3, VERSION 3.2.1 WILL NOT CURRENTLY FUNCTION

filePath <- "~/box/"
dataFile <- "c321-growthcurve-f500-101115.xlsx"
keyFile <- "key-c321-growthcurve-f500-101115.xlsx"

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
rawDataTbl <- read_excel(file.path(paste(filePath, dataFile, sep="")), col_names = FALSE, skip = 0) %>% slice(1:100)
rawKeyTbl <- read_excel(file.path(paste(filePath, keyFile, sep="")))
masterTbl <- parseIcontrolOutput(rawDataTbl, rawKeyTbl) %>% bind_rows(masterTbl) %>% mutate(value = as.numeric(value))

levels(masterTbl$measurementType) <- c("OD[600]")

masterTbl$strain <- factor(masterTbl$strain,
                           levels = c('MG1655','MJH116','OPT97','OPT99','OPT100','MJH184','OPT103','OPT104','OPT105', 'blank'))
masterTbl$supplements <- factor(masterTbl$supplements,
                                levels = c('Zeo', 'none', 'Zeo Gent IodoY'))
expTbl <- masterTbl %>% filter(strain != 'blank') %>% 
  group_by(time, strain, medium, supplements) %>%
  summarise(n = n(), meanValue = mean(value), sd = sd(value)) %>% 
  mutate(se = sd/sqrt(n), LCI = meanValue+qt(0.025,df=(n-1))*se, UCI = meanValue+qt(0.975, df=(n-1))*se) 

ggplot(expTbl, aes(x=time, y=meanValue, color=strain)) + 
  geom_line() +
  facet_grid(medium ~ supplements) +
  ggtitle("24 Hours Growth of OPT Lines and Ancestors by Medium") +
  scale_y_log10(limits = c(0.065,1.0), breaks=c(0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.825)) +
  scale_x_continuous(limits = c(0, 24.5)) +
  ylab(bquote(log[10] ~ OD[600])) +
  xlab("Time (h)") +
  theme(legend.position="right") +
  geom_errorbar(aes(x=time, ymin=LCI, ymax=UCI), size=.1, width=.25)

# levels(masterTbl$measurementType) <- c("OD[600]")
# 
# masterTbl$strain <- factor(masterTbl$strain,
#                        levels = c('MG1655','MJH116','OPT97','OPT99','OPT100','MJH184','OPT103','OPT104','OPT105', 'blank'))
# expTbl <- masterTbl %>% filter(strain != 'blank') %>% group_by(cycle, time, strain) #%>% filter(time > 1 , value > 0.01)
#   #%>% summarise(n = n(), meanValue = mean(value), sd = sd(value)) %>% 
#   #mutate(se = sd/sqrt(n), LCI = meanValue+qt(0.025,df=(n-1))*se, UCI = meanValue+qt(0.975, df=(n-1))*se)
# 
# 
# theme_set(theme_bw(base_size = 16))
# theme_update(panel.border=element_rect(color="black", fill=NA, size=1), legend.key=element_rect(color=NA, fill=NA))
# line_thickness = 0.8
# theme_update(legend.position=c(0.985,0.075), legend.background = element_rect(colour =alpha('white',0),fill=alpha('white', 0)), 
#              legend.text = element_text(size = 10), legend.title = element_text(size=10),
#              legend.key.size = unit(0.5, "cm"),
#              panel.grid.minor=element_blank(), panel.grid.major=element_blank())
# 
# expTbl %>% ggplot(aes(x=time, y=value), measurementType) + 
#   geom_line(aes(color=replicate)) +
#   facet_grid(medium ~ strain, labeller=label_parsed) +
#   #scale_y_continuous(trans=log10_trans()) +
#   scale_x_discrete(breaks=c(0,5,10,15,20,25)) +
#   ggtitle("24 Hours Growth of OPT Lines and Ancestors by Medium") +
#   #ylab(bquote(log[10] ~ OD[600])) +
#   xlab("Time (h)") #+
# #  geom_errorbar(aes(x=time, ymin=LCI, ymax=UCI), color=replicate, size=.05, width=.05)
