library(readxl)
library(dplyr)
library(magrittr)
library(tidyr)
library(ggplot2)
library(gridExtra)
require(cowplot)

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
    mutate(well = as.factor(well),
           plate_column=as.factor(gsub("\\d", "", well, perl=T)),
           plate_row=as.factor(gsub("\\D", "", well, perl=T)),
           measurementType = as.factor(measurementType),
           cycle = as.factor(cycle),
           medium = as.factor(medium), clone = as.factor(clone), aaRS = as.factor(aaRS),
           AA = as.factor(AA), codon = as.factor(codon), replicate = as.factor(replicate)) %>%
    arrange(well)


}

analyze_microplate <- function(dataFile, keyFile, output_base_name) {

  ## Set for debugging
  #dataFile = "input_examples/3-aminotyrosine-032615.xlsx"
  #keyFile = "input_examples/key-3-aminotyrosine-032615.xlsx"
  #output_base_name = "output/3-aminotyrosine-032615"

  ## Global settings for graphs
  theme_set(theme_bw(base_size = 24))
  theme_update(panel.border=element_rect(color="black", fill=NA, size=1), legend.key=element_rect(color=NA, fill=NA))

  plot.width = 6
  plot.height = 4

  masterTbl <- tbl_df(data.frame())
  for(i in excel_sheets(dataFile)) {
    rawDataTbl <- read_excel(dataFile, sheet = i, col_names = FALSE, skip = 0) %>% slice(1:300)
    rawKeyTbl <- read_excel(keyFile, sheet = i)
    masterTbl <- parseIcontrolOutput(rawDataTbl, rawKeyTbl) %>% mutate(plate = as.factor(i)) %>% bind_rows(masterTbl)
  }

  ## fix order of columns
  masterTbl$plate_row = factor(masterTbl$plate_row,  levels = order(levels(masterTbl$plate_row)))

  masterTbl <- masterTbl %>% mutate(aaRS = as.factor(aaRS)) %>% spread(measurementType, value)
  #write.csv(masterTbl, "masterTbl.csv")

  ## Now make graphs of raw data per well!
  # ggplot(masterTbl, aes(time, GFS, color=codon, linetype=AA)) +
  #   geom_line() +
  #   coord_cartesian() +
  #   facet_grid(plate_column ~ plate_row, scales="fixed")
  # ggsave(paste0(output_base_name, ".raw.plate.GFS.pdf"), width=2*plot.width, height=2*plot.height)
  #
  # ggplot(masterTbl, aes(time, RFS, color=codon, linetype=AA)) +
  #   geom_line() +
  #   coord_cartesian() +
  #   facet_grid(plate_column ~ plate_row, scales="fixed")
  # ggsave(paste0(output_base_name, ".raw.plate.RFS.pdf"), width=2*plot.width, height=2*plot.height)
  #
  # ggplot(masterTbl, aes(time, OD600, color=codon, linetype=AA)) +
  #   geom_line() +
  #   coord_cartesian() +
  #   facet_grid(plate_column ~ plate_row, scales="fixed")
  # ggsave(paste0(output_base_name, ".raw.plate.OD600.pdf"), width=12, height=8)

  expTbl <- masterTbl %>%
    filter(codon != "blank" & codon != "none") %>%
    group_by(cycle, time, AA, aaRS, codon, replicate) %>%
    select(-medium, -well) %>%
    unite(commonType, time, AA, sep="|")
  #write.csv(expTbl, "expTbl.csv")

  blankTbl <- masterTbl %>%
    filter(medium == "LB", clone == "blank", codon == "blank") %>%
    group_by(cycle, time, AA, aaRS, codon) %>%
    select(-medium, -well) %>%
    summarise(blankOD600 = mean(OD600), blankRFS = mean(RFS), blankGFS = mean(GFS)) %>%
    # , n = n(), sd = sd(value)) %>%
    #mutate(se = sd/sqrt(n), LCI = blankValue+qt(0.025,df=(n-1))*se, UCI = meanValue+qt(0.975, df=(n-1))*se) %>%
    unite(commonType, time, AA, sep="|")

  #write.csv(blankTbl, "blankTbl.csv")


  ## These tables are just for the graphs

  expBlankMasterTbl <- expTbl %>% inner_join(blankTbl, by = "commonType") %>%
    separate(commonType, into = c("time", "AA"), sep = "\\|") %>%
  #  select(-aaRS.y, -codon.y, -cycle.y, plate) %>%
    transmute(time = time, AA = AA, replicate = replicate, aaRS = aaRS.x, clone = clone, codon = codon.x,
              rawGFS = GFS, rawOD600 = OD600, rawRFS = RFS, blankOD600 = blankOD600, blankRFS = blankRFS, blankGFS = blankGFS)

  expPlusBlankTbl <- expBlankMasterTbl %>% mutate(netGFS = rawGFS-blankGFS, netOD600=rawOD600-blankOD600, netRFS=rawRFS-blankRFS,
           normGFS = netGFS/netOD600, normRFS = netRFS/netOD600, ratio = normRFS/normGFS) %>%
    gather("measurementType", "value", -time, -clone, -aaRS, -codon, -AA, -replicate) %>%
    mutate(time=as.numeric(time), AA=as.factor(AA), replicate=as.factor(replicate)) %>%
    group_by(time, aaRS, measurementType, codon)


  rawBlankNetTbl <- expBlankMasterTbl %>%
    mutate(netGFS = rawGFS-blankGFS, netOD600=rawOD600-blankOD600, netRFS=rawRFS-blankRFS) %>%
    gather("measurementType", "value", -time, -clone, -aaRS, -codon, -AA, -replicate) %>%
    mutate(time=as.numeric(time), AA=as.factor(AA), replicate=as.factor(replicate)) %>%
    group_by(time, aaRS, measurementType, codon)

  ## This is a stopgap fix to get the linetype assignment to show present as solid/absent as dashed while plotting
  presentRawBlankNetTbl <- rawBlankNetTbl %>% filter(AA == "present") %>% mutate(AA=as.factor("present"))
  absentRawBlankNetTbl <- rawBlankNetTbl %>% filter(AA == "absent") %>% mutate(AA=as.factor("absent"))

  plotRawBlankNetTbl <- presentRawBlankNetTbl %>% full_join(absentRawBlankNetTbl) %>% mutate(AA=as.factor(AA)) %>%
    group_by(time, aaRS, measurementType, codon, AA)
  plotRawBlankNetTbl$AA <- factor(plotRawBlankNetTbl$AA, levels=c("present","absent"))

  ggplot(plotRawBlankNetTbl, aes(time, value, group=interaction(replicate, AA, codon, clone), color=codon, linetype=AA)) +
    geom_line() +
    coord_cartesian() +
    facet_grid(measurementType ~ aaRS, scales = "free_y")
  ggsave(paste0(output_base_name, ".all.signals.pdf"), width=plot.width, height=7*plot.height)

  ## Continue with main tables

  ## Find the exponential time window

  netOD600Tbl <- expBlankMasterTbl %>%
    mutate(netOD600=rawOD600-blankOD600) %>%
    mutate(time=as.numeric(time), AA=as.factor(AA), replicate=as.factor(replicate)) %>%
    group_by(time, aaRS, codon)

  # minExpOD600 = 0.3
  # maxExpOD600 = 0.4
  # maxFilteredExpOD600Tbl = netOD600Tbl %>% filter(netOD600 >= maxExpOD600)
  # maxExpTime = min(maxFilteredExpOD600Tbl$time)
  # minFilteredExpOD600Tbl = netOD600Tbl %>% filter(netOD600 <= minExpOD600)
  # minExpTime = max(minFilteredExpOD600Tbl$time)

  p1 = ggplot(netOD600Tbl, aes(time, netOD600, group=interaction(replicate, AA, codon, clone), color=codon, linetype=AA)) +
    geom_line() +
    coord_cartesian(ylim=c(0, 0.8), xlim=c(0,15)) #+
    # ggsave(paste0(output_base_name, ".net_OD600.pdf"), width=plot.width, height=plot.height)

  ## Return to calculations of ratios

  normGFSTbl <- expBlankMasterTbl %>%
    mutate(netGFS = rawGFS-blankGFS, netOD600=rawOD600-blankOD600, netRFS=rawRFS-blankRFS) %>%
    transmute(time = time, AA = AA, replicate = replicate, aaRS = aaRS, clone = clone, codon = codon,
              normGFS = netGFS/netOD600) %>%
    gather("measurementType", "value", -time, -clone, -aaRS, -codon, -AA, -replicate) %>%
    mutate(time=as.numeric(time), AA=as.factor(AA), replicate=as.factor(replicate)) %>%
    group_by(time, aaRS, measurementType, codon)

  p3 = ggplot(normGFSTbl, aes(time, value, group=interaction(replicate, AA, codon, clone), color=codon, linetype=AA)) +
    geom_line() +
    coord_cartesian(ylim=c(-5000, 45000), xlim=c(0,15)) #+
  # ggsave(paste0(output_base_name, ".normalized_GFP.pdf"), width=plot.width, height=plot.height)

  normRFSTbl <- expBlankMasterTbl %>%
    mutate(netGFS = rawGFS-blankGFS, netOD600=rawOD600-blankOD600, netRFS=rawRFS-blankRFS) %>%
    transmute(time = time, AA = AA, replicate = replicate, aaRS = aaRS, clone = clone, codon = codon, normRFS = netRFS/netOD600) %>%
    gather("measurementType", "value", -time, -clone, -aaRS, -codon, -AA, -replicate) %>%
    mutate(time=as.numeric(time), AA=as.factor(AA), replicate=as.factor(replicate)) %>%
    group_by(time, aaRS, measurementType, codon)

  p2 = ggplot(normRFSTbl, aes(time, value, group=interaction(replicate, AA, codon, clone), color=codon, linetype=AA)) +
    geom_line() +
    coord_cartesian(ylim=c(-2000, 8000), xlim=c(0,15)) #+
  # ggsave(paste0(output_base_name, ".normalized_RFP.pdf"), width=plot.width, height=plot.height)


  ## Construct final tables of normalized RFP and GFP signals

  amberNormRFSTbl <- normRFSTbl %>%
    filter(codon == "amber") %>%
    group_by(time, aaRS, AA, clone, codon, replicate) %>%
    unite(commonCodon, time, aaRS, AA, clone, replicate, sep="|")
  #write.csv(amberTbl, "amberTbl.csv")

  amberNormGFSTbl <- normGFSTbl %>%
    filter(codon == "amber") %>%
    group_by(time, aaRS, AA, clone, codon, replicate) %>%
    unite(commonCodon, time, aaRS, AA, clone, replicate, sep="|")
  #write.csv(amberTbl, "amberTbl.csv")

  tyrosineNormRFSTbl <- normRFSTbl %>%
    filter(codon == "tyrosine") %>%
    group_by(time, aaRS, AA, clone, codon, replicate) %>%
    unite(commonCodon, time, aaRS, AA, clone, replicate, sep="|")
  #write.csv(amberTbl, "amberTbl.csv")

  tyrosineNormGFSTbl <- normGFSTbl %>%
    filter(codon == "tyrosine") %>%
    group_by(time, aaRS, AA, clone, codon, replicate) %>%
    unite(commonCodon, time, aaRS, AA, clone, replicate, sep="|")
  #write.csv(amberTbl, "amberTbl.csv")

  ## Construct tables of RFP/GFP ratios

  amberGFSToRFSRatioTbl <- amberNormGFSTbl %>%
    left_join(amberNormRFSTbl, by = "commonCodon") %>%
    mutate(value = value.x / value.y) %>%
    select(commonCodon, value)

  tyrosineGFSToRFSRatioTbl <- tyrosineNormGFSTbl %>%
    left_join(tyrosineNormRFSTbl, by = "commonCodon") %>%
    mutate(value = value.x / value.y) %>%
    select(commonCodon, value)

  ## Construct table of Amber to Tyrosine ratio of RFP/GFP ratios

  amberToTyrosineRatioGFSToRFSRatioTbl <- amberGFSToRFSRatioTbl %>%
    left_join(tyrosineGFSToRFSRatioTbl, by = "commonCodon") %>%
    mutate(value = value.x / value.y) %>%
    select(commonCodon, value) %>%
    separate(commonCodon, into = c("time", "aaRS", "AA", "clone", "replicate"), sep = "\\|") %>%
    mutate(time=as.numeric(time), AA=as.factor(AA), clone=as.factor(clone), replicate=as.factor(replicate))

  p4 = ggplot(amberToTyrosineRatioGFSToRFSRatioTbl, aes(time, value, group=interaction(replicate, AA, clone), linetype=AA, color=replicate)) +
    geom_line() +
    coord_cartesian(ylim=c(-0.2, 1.8), xlim=c(0,15)) +
    geom_hline(aes(yintercept=1), linetype="dashed")

  # ggsave(paste0(output_base_name, ".RRE.pdf"), width=plot.width, height=plot.height)

  ## Construct tables of RFP/GFP ratios

  absentAmberToTyrosineRatioGFSToRFSRatioTbl <- amberToTyrosineRatioGFSToRFSRatioTbl %>%
    filter(AA == "absent") %>%
    group_by(time, aaRS, clone, replicate) %>%
    unite(commonAA, time, aaRS, clone, replicate, sep="|")

  presentAmberToTyrosineRatioGFSToRFSRatioTbl <- amberToTyrosineRatioGFSToRFSRatioTbl %>%
    filter(AA != "absent") %>%
    group_by(time, aaRS, clone, replicate) %>%
    unite(commonAA, time, aaRS, clone, replicate, sep="|")

  absentToPresentRatioAmberToTyrosineRatioGFSToRFSRatioTbl <- absentAmberToTyrosineRatioGFSToRFSRatioTbl %>%
    left_join(presentAmberToTyrosineRatioGFSToRFSRatioTbl, by = "commonAA") %>%
    mutate(value = value.x / value.y) %>%
    select(commonAA, value) %>%
    separate(commonAA, into = c("time", "aaRS", "clone", "replicate"), sep = "\\|") %>%
    mutate(time=as.numeric(time), clone=as.factor(clone), replicate=as.factor(replicate))

  p5 = ggplot(absentToPresentRatioAmberToTyrosineRatioGFSToRFSRatioTbl, aes(time, value, group=interaction(replicate, clone), color=clone)) +
    geom_line() +
    coord_cartesian(ylim=c(-0.2, 1.8), xlim=c(0,15)) +
    geom_hline(aes(yintercept=1), linetype="dashed")

  # ggsave(paste0(output_base_name, ".R3C.pdf"), width=plot.width, height=plot.height)

  ## Plotting with 95% confidence interval

  x <- absentToPresentRatioAmberToTyrosineRatioGFSToRFSRatioTbl %>%
    group_by(time, clone) %>%
    summarize( n=n(), mean=mean(value),sd=sd(value), var=var(value)  ) %>%
    mutate (CI.size = sd * qt(0.975,n-1)/sqrt(n-1)) %>%
    mutate(CI95.L = mean-CI.size, CI95.U = mean+CI.size )

  ggplot(subset(x), aes(time, mean)) +
    geom_line(aes(color=clone)) +
    #geom_ribbon(aes(ymin=CI95.L, ymax=CI95.U, fill=clone), alpha=0.3) +
    geom_errorbar(aes(ymin=CI95.L, ymax=CI95.U, color=clone), width=0, alpha=0.5) +
    coord_cartesian(ylim=c(-0.2, 1.8), xlim=c(0,15)) +
    geom_hline(aes(yintercept=1), linetype="dashed")

  ggsave(paste0(output_base_name, ".R3C.CI.pdf"), width=plot.width, height=plot.height)

  plot_grid(p1, p2, p3, p4, p5, align='vh', ncol=1)
  ggsave(paste0(output_base_name, ".final.pdf"), width=plot.width, height=plot.height*5)
}
