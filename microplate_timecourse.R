library(readxl)
library(dplyr)
library(magrittr)
library(tidyr)
library(ggplot2)

theme_set(theme_bw(base_size = 8))
theme_update(panel.border=element_rect(color="black", fill=NA, size=1), legend.key=element_rect(color=NA, fill=NA))
line_thickness = 0.8

parseExperimentData <- function(data, lookupTable, numMeasurementTypes = 3, numWells = 96) {
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
      mutate(well = as.factor(well), measurementType = as.factor(measurementType), cycle = as.factor(cycle), 
             medium = as.factor(medium), clone = as.factor(clone), aaRS = as.factor(aaRS), 
             AA = as.factor(AA), codon = as.factor(codon)) %>%
      arrange(well)
}
masterTbl <- tbl_df(data.frame())
for(i in excel_sheets(file.path("~/data2011.xlsx"))){
  rawDataTbl <- read_excel(file.path("~/data2011.xlsx"), sheet = i, col_names = FALSE, skip = 52) %>% slice(1:300)
  rawKeyTbl <- read_excel(file.path("~/key2011.xlsx"), sheet = i)
  masterTbl <- parseExperimentData(rawDataTbl, rawKeyTbl) %>% bind_rows(masterTbl)
}

masterTbl <- masterTbl %>% filter(cycle != "NA")
masterTbl <- masterTbl %>% mutate(aaRS = as.factor(aaRS))

expTbl <- masterTbl %>% filter(codon != "none")

levels(expTbl$measurementType) <- c("GFS", "OD[600]", "RFS")
expTbl$AA <- relevel(expTbl$AA, "present")
expTbl <- expTbl %>% 
  group_by(aaRS, AA, codon,clone, measurementType, time) %>% summarise(n=n(), meanValue = mean(value))#, sd=sd(value)) %>% 
  #mutate(se=sd/sqrt(n),LCI=meanValue+qt(0.025,df=(n-1))*se,UCI=meanValue+qt(0.975, df=(n-1))*se)

theme_update(axis.title.y = element_blank(), panel.grid.minor=element_blank(), panel.grid.major=element_blank())
ggplot(expTbl) + geom_line(aes(x=time, y=meanValue, color=codon, linetype=AA)) + 
  facet_grid(clone ~ measurementType ~ aaRS, scales="free_y", labeller=label_parsed) + 
  xlab("Time (h)") #+ 
  #geom_errorbar(aes(x=time, ymin=LCI, ymax=UCI, color=codon), size=.1, width=.05) #+ 
  #scale_x_continuous(limits=c(0, 16))
