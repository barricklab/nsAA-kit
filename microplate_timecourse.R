library(readxl)
library(dplyr)
library(magrittr)
library(tidyr)
library(ggplot2)

# Works using R 3.2.2 (Fire Safety) and dplyr 0.4.3, does not work using R 3.1.3 (Smooth Sidewalk) and dplyr 0.4.2
filePath <- "~/box/"
dataFile <- "grouped-data-output.xlsx"
keyFile <- "grouped-key.xlsx"

parseIcontrolOutput <- function(data, lookupTable, numMeasurementTypes = 3, numWells = 96) {
  rowsPerMeasure <- numWells + 4
  data %T>% {
    colnames(.) <- c("well", paste(data[2,][2:ncol(data)],
                                   data[3,][2:ncol(data)], sep ="|"))
  } %>% {
    measurementLocations <- seq_along(rownames(data))[seq(1, (rowsPerMeasure *
                                                                numMeasurementTypes), rowsPerMeasure)]
    mutate(., measurementType = rep(data[[1]][measurementLocations],
                                    each = rowsPerMeasure)) %>% {
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
    mutate(well = as.factor(well), measurementType = as.factor(measurementType),
           cycle = as.factor(cycle),
           medium = as.factor(medium), clone = as.factor(clone), aaRS = as.factor(aaRS),
           AA = as.factor(AA), codon = as.factor(codon), replicate = as.factor(replicate)) %>%
    arrange(well)
}
masterTbl <- tbl_df(data.frame())
for(i in excel_sheets(file.path(paste(filePath, dataFile, sep="")))) {
  rawDataTbl <- read_excel(file.path(paste(filePath, dataFile, sep="")), sheet = i, col_names = FALSE, skip = 0) %>% slice(1:300)
  rawKeyTbl <- read_excel(file.path(paste(filePath, keyFile, sep="")), sheet = i)
  masterTbl <- parseIcontrolOutput(rawDataTbl, rawKeyTbl) %>% mutate(plate = as.factor(i)) %>% bind_rows(masterTbl)
}

masterTbl <- masterTbl %>% mutate(aaRS = as.factor(aaRS)) %>% spread(measurementType, value)

expTbl <- masterTbl %>%
  filter(codon != "blank" & codon != "none") %>%
  group_by(cycle, time, AA, aaRS, codon, replicate) %>%
  select(-medium, -well) %>%
  unite(commonType, time, AA, replicate, sep="|")

blankTbl <- masterTbl %>%
  filter(medium == "LB", clone == "blank", codon == "blank") %>%
  group_by(cycle, time, AA, aaRS, codon, replicate) %>%
  select(-medium, -well) %>%
  unite(commonType, time, AA, replicate, sep="|")

expPlusBlankTbl <- expTbl %>% inner_join(blankTbl, by = "commonType") %>%
  separate(commonType, into = c("time", "AA", "replicate"), sep = "\\|") %>%
  select(-clone.y, -aaRS.y, -codon.y) %>%
  transmute(time = time, AA = AA, replicate = replicate, aaRS = aaRS.x, clone = clone.x,
            codon = codon.x, rawGFS = GFS.x, rawOD600 = OD600.x, rawRFS = RFS.x,
            blankGFS = GFS.y, blankOD600 = OD600.y, blankRFS = RFS.y) %>%
  mutate(netGFS = rawGFS-blankGFS, netOD600=rawOD600-blankOD600, netRFS=rawRFS-blankRFS,
         normGFS = netGFS/netOD600, normRFS = netRFS/netOD600, ratio = normRFS/normGFS) %>%
  gather("measurementType", "value", -time, -clone, -aaRS, -codon, -AA, -replicate) %>%
  mutate(time=as.numeric(time), AA=as.factor(AA), replicate=as.factor(replicate)) %>%
  group_by(time, aaRS, measurementType, codon)

rawBlankNetTbl <- expTbl %>% inner_join(blankTbl, by = "commonType") %>%
  separate(commonType, into = c("time", "AA", "replicate"), sep = "\\|") %>%
  select(-clone.y, -aaRS.y, -codon.y) %>%
  transmute(time = time, AA = AA, replicate = replicate, aaRS = aaRS.x, clone = clone.x,
            codon = codon.x, rawGFS = GFS.x, rawOD600 = OD600.x, rawRFS = RFS.x,
            blankGFS = GFS.y, blankOD600 = OD600.y, blankRFS = RFS.y) %>%
  mutate(netGFS = rawGFS-blankGFS, netOD600=rawOD600-blankOD600, netRFS=rawRFS-blankRFS) %>%
  gather("measurementType", "value", -time, -clone, -aaRS, -codon, -AA, -replicate) %>%
  mutate(time=as.numeric(time), AA=as.factor(AA), replicate=as.factor(replicate)) %>%
  group_by(time, aaRS, measurementType, codon)

presentRawBlankNetTbl <- rawBlankNetTbl %>% filter(AA != "absent") %>% mutate(AA=as.factor("present"))
absentRawBlankNetTbl <- rawBlankNetTbl %>% filter(AA == "absent") %>% mutate(AA=as.factor("absent"))

plotRawBlankNetTbl <- presentRawBlankNetTbl %>% full_join(absentRawBlankNetTbl) %>% mutate(AA=as.factor(AA)) %>%
  group_by(time, aaRS, measurementType, codon, AA)
plotRawBlankNetTbl$AA <- factor(plotRawBlankNetTbl$AA, levels=c("present","absent"))
ggplot(plotRawBlankNetTbl, aes(time, value, group=interaction(replicate, AA, codon, clone), color=codon, linetype=AA)) + 
   geom_line() +
   facet_grid(measurementType ~ aaRS, scales = "free_y")

normGFSTbl <- expTbl %>% inner_join(blankTbl, by = "commonType") %>%
  separate(commonType, into = c("time", "AA", "replicate"), sep = "\\|") %>%
  select(-clone.y, -aaRS.y, -codon.y) %>%
  transmute(time = time, AA = AA, replicate = replicate, aaRS = aaRS.x, clone = clone.x,
            codon = codon.x, rawGFS = GFS.x, rawOD600 = OD600.x, rawRFS = RFS.x,
            blankGFS = GFS.y, blankOD600 = OD600.y, blankRFS = RFS.y) %>%
  mutate(netGFS = rawGFS-blankGFS, netOD600=rawOD600-blankOD600, netRFS=rawRFS-blankRFS) %>% 
  transmute(time = time, AA = AA, replicate = replicate, aaRS = aaRS, clone = clone, codon = codon,
            normGFS = netGFS/netOD600) %>%
  gather("measurementType", "value", -time, -clone, -aaRS, -codon, -AA, -replicate) %>%
  mutate(time=as.numeric(time), AA=as.factor(AA), replicate=as.factor(replicate)) %>%
  group_by(time, aaRS, measurementType, codon)

presentNormGFSTbl <- normGFSTbl %>% filter(AA != "absent") %>% mutate(AA=as.factor("present"))
absentNormGFSTbl <- normGFSTbl %>% filter(AA == "absent") %>% mutate(AA=as.factor("absent"))

plotNormGFSTbl <- presentNormGFSTbl %>% full_join(absentNormGFSTbl) %>% mutate(AA=as.factor(AA)) %>%
  group_by(time, aaRS, measurementType, codon, AA)
plotNormGFSTbl$AA <- factor(plotNormGFSTbl$AA, levels=c("present","absent"))
ggplot(plotNormGFSTbl, aes(time, value, group=interaction(replicate, AA, codon, clone), color=codon, linetype=AA)) + 
  geom_line() +
  coord_cartesian(ylim=c(0, 75000)) +
  facet_grid(measurementType ~ aaRS)

normRFSTbl <- expTbl %>% inner_join(blankTbl, by = "commonType") %>%
  separate(commonType, into = c("time", "AA", "replicate"), sep = "\\|") %>%
  select(-clone.y, -aaRS.y, -codon.y) %>%
  transmute(time = time, AA = AA, replicate = replicate, aaRS = aaRS.x, clone = clone.x,
            codon = codon.x, rawGFS = GFS.x, rawOD600 = OD600.x, rawRFS = RFS.x,
            blankGFS = GFS.y, blankOD600 = OD600.y, blankRFS = RFS.y) %>%
  mutate(netGFS = rawGFS-blankGFS, netOD600=rawOD600-blankOD600, netRFS=rawRFS-blankRFS) %>% 
  transmute(time = time, AA = AA, replicate = replicate, aaRS = aaRS, clone = clone, codon = codon, normRFS = netRFS/netOD600) %>%
  gather("measurementType", "value", -time, -clone, -aaRS, -codon, -AA, -replicate) %>%
  mutate(time=as.numeric(time), AA=as.factor(AA), replicate=as.factor(replicate)) %>%
  group_by(time, aaRS, measurementType, codon)

presentNormRFSTbl <- normRFSTbl %>% filter(AA != "absent") %>% mutate(AA=as.factor("present"))
absentNormRFSTbl <- normRFSTbl %>% filter(AA == "absent") %>% mutate(AA=as.factor("absent"))

plotNormRFSTbl <- presentNormRFSTbl %>% full_join(absentNormRFSTbl) %>% mutate(AA=as.factor(AA)) %>%
  group_by(time, aaRS, measurementType, codon, AA)
plotNormRFSTbl$AA <- factor(plotNormRFSTbl$AA, levels=c("present","absent"))
ggplot(plotNormRFSTbl, aes(time, value, group=interaction(replicate, AA, codon, clone), color=codon, linetype=AA)) + 
  geom_line() +
  coord_cartesian(ylim=c(0, 12500)) +
  facet_grid(measurementType ~ aaRS)

amberTbl <- expPlusBlankTbl %>%
  filter(codon == "amber") %>%
  group_by(time, aaRS, AA, clone, codon, replicate) %>%
  unite(commonCodon, time, aaRS, AA, clone, replicate, measurementType, sep="|")

tyrosineTbl <- expPlusBlankTbl %>%
  filter(codon == "tyrosine") %>%
  group_by(time, aaRS, AA, clone, codon, replicate) %>%
  unite(commonCodon, time, aaRS, AA, clone, replicate, measurementType, sep="|")

amberTyrosineTbl <- amberTbl %>% inner_join(tyrosineTbl, by = "commonCodon") %>%
  separate(commonCodon, into = c("time", "aaRS", "AA", "clone", "replicate", "measurementType"), sep = "\\|") %>%
  select(-codon.x, -codon.y) %>%
  transmute(time = time, AA = AA, aaRS = aaRS, clone = clone,
    replicate = replicate, amberRatio = value.x, tyrosineRatio = value.y)

absentTbl <- amberTyrosineTbl %>%
  filter(AA == "absent") %>%
  group_by(time, aaRS, clone, replicate) %>%
  unite(commonAA, time, aaRS, clone, replicate, sep="|")

presentTbl <- amberTyrosineTbl %>%
  filter(AA != "absent") %>%
  group_by(time, aaRS, clone, replicate) %>%
  unite(commonAA, time, aaRS, clone, replicate, sep="|")

absentPresentTbl <- absentTbl %>% inner_join(presentTbl, by = "commonAA") %>%
  separate(commonAA, into = c("time", "aaRS", "clone", "replicate"), sep = "\\|") %>%
  select(-AA.x, -AA.y) %>%
  transmute(time = time, aaRS = aaRS, clone = clone,  replicate = replicate,
    absentAmberRatio = amberRatio.x, absentTyrosineRatio =  tyrosineRatio.x,
    presentAmberRatio = amberRatio.y, presentTyrosineRatio = tyrosineRatio.y) %>%
  mutate(decodingEfficiencyAbsent = absentAmberRatio/absentTyrosineRatio,
    decodingEfficiencyPresent = presentAmberRatio/presentTyrosineRatio,
    misincorporationValue = decodingEfficiencyAbsent/decodingEfficiencyPresent)

finalTbl <- absentPresentTbl %>%
  mutate(aaRS = as.factor(aaRS), time = as.numeric(time)) %>%
  group_by(time, aaRS) %>%
  summarise(n=n(), meanAbsentAmberRatio = mean(absentAmberRatio),
    meanAbsentTyrosineRatio = mean(absentTyrosineRatio),
    meanPresentAmberRatio = mean(presentAmberRatio),
    meanPresentTyrosineRatio = mean(presentTyrosineRatio),
    meanDecodingEfficiencyAbsent = mean(decodingEfficiencyAbsent),
    meanDecodingEfficiencyPresent = mean(decodingEfficiencyPresent),
    meanMisincorporationValue = mean(misincorporationValue)) %>%
  gather("statisticType", "statistic", 4:10)
mAARtbl <- finalTbl %>% filter(statisticType=="meanAbsentAmberRatio")
ggplot(mAARtbl) + geom_line(aes(time, statistic, group=aaRS)) +
  coord_cartesian(ylim=c(0,12000)) +
  scale_x_continuous(limits=c(0,15), breaks=seq(0,15,5)) +
  facet_grid(statisticType ~ aaRS)
mATRtbl <- finalTbl %>% filter(statisticType=="meanAbsentTyrosineRatio")
ggplot(mATRtbl) + geom_line(aes(time, statistic, group=aaRS)) +
  coord_cartesian(ylim=c(0,12000)) +
  scale_x_continuous(limits=c(0,15), breaks=seq(0,15,5)) +
  facet_grid(statisticType ~ aaRS)
mPARtbl <- finalTbl %>% filter(statisticType=="meanPresentAmberRatio")
ggplot(mPARtbl) + geom_line(aes(time, statistic, group=aaRS)) +
  coord_cartesian(ylim=c(0,12000)) +
  scale_x_continuous(limits=c(0,15), breaks=seq(0,15,5)) +
  facet_grid(statisticType ~ aaRS) +
  geom_hline(yintercept=1, linetype=2)
mPTRtbl <- finalTbl %>% filter(statisticType=="meanPresentTyrosineRatio")
ggplot(mPTRtbl) + geom_line(aes(time, statistic, group=aaRS)) +
  coord_cartesian(ylim=c(0,12000)) +
  scale_x_continuous(limits=c(0,15), breaks=seq(0,15,5)) +
  facet_grid(statisticType ~ aaRS) +
  geom_hline(yintercept=1, linetype=2)
mDEAtbl <- finalTbl %>% filter(statisticType=="meanDecodingEfficiencyAbsent")
ggplot(mDEAtbl) + geom_line(aes(time, statistic, group=aaRS)) +
  coord_cartesian(ylim=c(0,1.5)) +
  scale_x_continuous(limits=c(0,15), breaks=seq(0,15,5)) +
  facet_grid(statisticType ~ aaRS) +
  geom_hline(yintercept=1, linetype=2)
mDEPtbl <- finalTbl %>% filter(statisticType=="meanDecodingEfficiencyPresent")
ggplot(mDEPtbl) + geom_line(aes(time, statistic, group=aaRS)) +
  coord_cartesian(ylim=c(0,1.5)) +
  scale_x_continuous(limits=c(0,15), breaks=seq(0,15,5)) +
  facet_grid(statisticType ~ aaRS) +
  geom_hline(yintercept=1, linetype=2)
mMVtbl <- finalTbl %>% filter(statisticType=="meanMisincorporationValue")
ggplot(mMVtbl) + geom_line(aes(time, statistic, group=aaRS)) +
  coord_cartesian(ylim=c(0,5)) +
  scale_x_continuous(limits=c(0,15), breaks=seq(0,15,5)) +
  facet_grid(statisticType ~ aaRS) +
  geom_hline(yintercept=1, linetype=2)
ggplot(finalTbl) + geom_line(aes(time, statistic, group=aaRS)) +
  scale_x_continuous(limits=c(0,15), breaks=seq(0,15,5)) +
  facet_grid(statisticType ~ aaRS, scales="free_y") +
  geom_hline(yintercept=1, linetype=2)