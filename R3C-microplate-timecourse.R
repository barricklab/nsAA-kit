require(readxl)
require(magrittr)
require(tidyr)
require(dplyr)
require(ggplot2)
require(gridExtra)
require(cowplot)

parseIcontrolOutput <- function(data, lookupTable, numMeasurementTypes = 3, numWells = 96) {
  
  rowsPerMeasure <- numWells + 4
  data = data[1:(rowsPerMeasure*numMeasurementTypes),]
  data = data %T>% {
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
    inner_join(lookupTable, by = "well")
  
  
  data = data %>% mutate(well = as.factor(well))
  data = data %>% mutate(plate_column=as.factor(gsub("\\d", "", well, perl=T)))
  data = data %>% mutate(plate_row=as.factor(gsub("\\D", "", well, perl=T)))
  data = data %>% mutate(measurementType = as.factor(measurementType))
  data = data %>% mutate(cycle = as.factor(cycle))
  
  if ( "well" %in% names(data) ) {
    data = data %>% mutate(well = as.factor(well))
  }
  
  if ( "medium" %in% names(data) ) {
    data = data %>% mutate(medium = as.factor(medium))
  }
  
  if ( "clone" %in% names(data) ) {
    data = data %>% mutate(clone = as.factor(clone))
  }
  
  if ( "aaRS" %in% names(data) ) {
    data = data %>% mutate(aaRS = as.factor(aaRS))
  }
  
  if ( "AA" %in% names(data) ) {
    data = data %>% mutate(AA = as.factor(AA))
  }
  
  if ( "codon" %in% names(data) ) {
    data = data %>% mutate(codon = as.factor(codon))
  }
  
  if ( "replicate" %in% names(data) ) {
    data = data %>% mutate(replicate = as.factor(replicate))
  }

  return(data)
}

average_over_OD_window <- function(output_base_name, OD.range)
{
  ## For debug
  #output_base_name = "output/3-aminotyrosine-032615"
  return.averages = c();
  
  minExpOD600 = OD.range[1];
  maxExpOD600 = OD.range[2];
  
  signalTbl = read.csv(paste0(output_base_name, ".signal.csv"))
  RRETbl = read.csv(paste0(output_base_name, ".RRE.csv"))
  R3CTbl = read.csv(paste0(output_base_name, ".R3C.csv"))
  
  clone.list = levels(as.factor(signalTbl$clone))
  for (this.clone in clone.list) {
    this.clone.data = signalTbl %>% filter(clone == this.clone)
    replicate.list = levels(as.factor(this.clone.data$replicate))
    
    #now figure out the value to keep
    for (this.replicate in replicate.list) {
      
      this.replicate.data = this.clone.data %>% filter(replicate == this.replicate)
      
      ##average OD600 over all 4 curves used in calculating R3C value for this replicate!!
      avgOD600Tbl = this.replicate.data %>% group_by(time) %>% summarize(avgOD600 = mean(netOD600))
      
      maxFilteredExpOD600Tbl = avgOD600Tbl %>% filter(avgOD600 >= maxExpOD600)
      maxExpTime = min(maxFilteredExpOD600Tbl$time)
      
      minFilteredExpOD600Tbl = avgOD600Tbl %>% filter(avgOD600 <= minExpOD600)
      minExpTime = max(minFilteredExpOD600Tbl$time)
      
      ##now
      
      R3C.window = R3CTbl %>% filter (replicate == this.replicate) %>% filter(clone == this.clone) %>% filter(time >= minExpTime) %>% filter(time <= maxExpTime)
      this.R3C = mean(R3C.window$R3C)
      
      RRE.window = RRETbl %>% filter (replicate == this.replicate) %>% filter(clone == this.clone) %>% filter(time >= minExpTime) %>% filter(time <= maxExpTime)
      this.RRE.present = mean(RRE.window$RRE.present)
      
      new.row = data.frame(aaRS = as.character(R3C.window$aaRS[1]), clone=this.clone, replicate=this.replicate, minTime = minExpTime, maxTime = maxExpTime, R3C=this.R3C, RRE.present=this.RRE.present)
      return.averages = return.averages %>% bind_rows(new.row)
      
      ##add row to output table
    }
  } 
  
  return(return.averages);
}

summarize_microplates <- function(summaryOutputFileName, outputFileList, OD.range) {
  
  #debug
  #OD.range = c(0.3, 0.4)
  #summaryOutputFileName = "output/summary_test"
  #outputFileList = c("output/3-aminotyrosine-032615", "output/3-iodotyrosine-032615")
  
  all.averages = c();
  for (dataFile in outputFileList) {
    cat(dataFile, "\n")
    this.averages = average_over_OD_window(dataFile, OD.range)
    all.averages = all.averages %>% bind_rows(this.averages)
    
#     if (length(levels(as.factor(this.averages$clone))) > 1) {
#       per_clone_model = lm(R3C~clone, this.averages)
#       overall_model = lm(R3C~1, this.averages)
#       a = anova(overall_model, per_clone_model,  test="LRT")
#       cat("p-value for LRT that clones have different mean: ", a["Pr(>Chi)"][[1]][2], "\n")
#     }
  }
  
  all.averages = all.averages %>% mutate(aaRS = as.factor(aaRS), clone = as.factor(clone), replicate = as.factor(replicate))
  
  ## Global settings for graphs
  theme_set(theme_bw(base_size = 12))
  theme_update(panel.border=element_rect(color="black", fill=NA, size=1), legend.key=element_rect(color=NA, fill=NA))
  
  plot.width = 6
  plot.height = 4
  
  ## Graphs with per-clone bars
  all.averages.RRE.present.ci =  all.averages %>%
                  group_by(aaRS, clone) %>%
                  summarize(n=n(), mean=mean(RRE.present),sd=sd(RRE.present),R3C.var=var(RRE.present)  ) %>%
                  mutate (ci = sd * qt(0.975,n-1)/sqrt(n-1))
  
  ggplot(all.averages.RRE.present.ci, aes(x=aaRS, y=mean, fill=clone)) + 
    geom_bar(position=position_dodge(), stat="identity") + 
    geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci),
                  width=0.2,                    # Width of the error bars
                  position=position_dodge(0.9))
  
  ggsave(paste0(summaryOutputFileName, ".per-clone.RRE.present.pdf"), width=plot.width, height=plot.height)
  
  all.averages.R3C.ci =  all.averages %>%
    group_by(aaRS, clone) %>%
    summarize(n=n(), mean=mean(R3C),sd=sd(R3C),R3C.var=var(R3C)  ) %>%
    mutate (ci = sd * qt(0.975,n-1)/sqrt(n-1))
  
  ggplot(all.averages.R3C.ci, aes(x=aaRS, y=mean, fill=clone)) + 
    geom_bar(position=position_dodge(), stat="identity") + 
    geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci),
                  width=0.2,                    # Width of the error bars
                  position=position_dodge(0.9))
  
  ggsave(paste0(summaryOutputFileName, ".per-clone.R3C.pdf"), width=plot.width, height=plot.height)
  
  ## Graphs with per-clone bars
  all.averages.RRE.present.ci =  all.averages %>%
    group_by(aaRS) %>%
    summarize(n=n(), mean=mean(RRE.present),sd=sd(RRE.present),R3C.var=var(RRE.present)  ) %>%
    mutate (ci = sd * qt(0.975,n-1)/sqrt(n-1))
  
  ggplot(all.averages.RRE.present.ci, aes(x=aaRS, y=mean)) + 
    geom_bar(position=position_dodge(), stat="identity") + 
    geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci),
                  width=0.2,                    # Width of the error bars
                  position=position_dodge(0.9))
  
  ggsave(paste0(summaryOutputFileName, ".per-aaRS.RRE.present.pdf"), width=plot.width, height=plot.height)
  
  all.averages.R3C.ci =  all.averages %>%
    group_by(aaRS) %>%
    summarize(n=n(), mean=mean(R3C),sd=sd(R3C),R3C.var=var(R3C)  ) %>%
    mutate (ci = sd * qt(0.975,n-1)/sqrt(n-1))
  
  ggplot(all.averages.R3C.ci, aes(x=aaRS, y=mean)) + 
    geom_bar(position=position_dodge(), stat="identity") + 
    geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci),
                  width=0.2,                    # Width of the error bars
                  position=position_dodge(0.9))
  
  ggsave(paste0(summaryOutputFileName, ".per-aaRS.R3C.pdf"), width=plot.width, height=plot.height)
  
  write.csv(all.averages, paste0(summaryOutputFileName, ".csv"))
}

## OD range is the range to perform the calculation of R3C across
analyze_microplate <- function(dataFile, keyFile, output_base_name) {

  ## Set for debugging
  dataFile = "input-examples/tyrosine-081415.xlsx"
  keyFile = "input-examples/key-tyrosine-081415.xlsx"
  output_base_name = "output/tyrosine-081415"

  ## Global settings for graphs
  theme_set(theme_bw(base_size = 24))
  theme_update(panel.border=element_rect(color="black", fill=NA, size=1), legend.key=element_rect(color=NA, fill=NA))

  plot.width = 6
  plot.height = 4

  rawDataTbl <- read_excel(dataFile, col_names = FALSE, skip = 0) %>% slice(1:300)
  rawKeyTbl <- read_excel(keyFile)
  masterTbl <- parseIcontrolOutput(rawDataTbl, rawKeyTbl)

  ## fix order of columns
  masterTbl$plate_row = factor(masterTbl$plate_row,  levels = order(levels(masterTbl$plate_row)))

  #Set so that plus amino acid ends up as solid line
  AA_levels = c("present", "absent", "none")
  
  masterTbl <- masterTbl %>% mutate(aaRS = as.factor(aaRS), AA = factor(AA, levels=AA_levels)) %>% spread(measurementType, value)
  #write.csv(masterTbl, "masterTbl.csv")

  make.raw = F
  if (make.raw) {
  ## Now make graphs of raw data per well!
   ggplot(masterTbl, aes(time, GFS, color=codon, linetype=AA)) +
     geom_line() +
     coord_cartesian() +
     facet_grid(plate_column ~ plate_row, scales="fixed")
   ggsave(paste0(output_base_name, ".raw.plate.GFS.pdf"), width=2*plot.width, height=2*plot.height)
  #
   ggplot(masterTbl, aes(time, RFS, color=codon, linetype=AA)) +
     geom_line() +
     coord_cartesian() +
     facet_grid(plate_column ~ plate_row, scales="fixed")
   ggsave(paste0(output_base_name, ".raw.plate.RFS.pdf"), width=2*plot.width, height=2*plot.height)
  #
   ggplot(masterTbl, aes(time, OD600, color=codon, linetype=AA)) +
     geom_line() +
     coord_cartesian() +
     facet_grid(plate_column ~ plate_row, scales="fixed")
   ggsave(paste0(output_base_name, ".raw.plate.OD600.pdf"), width=12, height=8)
  }
  
  expTbl <- masterTbl %>%
      filter(codon != "blank" & codon != "none") %>%
      group_by(cycle, time, AA, aaRS, codon, replicate) %>%
      ungroup() %>%
      select(-medium, -well) %>%
      unite(commonType, time, AA, sep="|")
  #write.csv(expTbl, "expTbl.csv")
  
  blankTbl <- masterTbl %>%
    filter(medium == "LB", clone == "blank", codon == "blank") %>%
    group_by(cycle, time, AA, aaRS, codon) %>%
    select(-medium, -well) %>%
    summarise(blankOD600 = mean(OD600), blankRFS = mean(RFS), blankGFS = mean(GFS)) %>%
    ungroup() %>%
    # , n = n(), sd = sd(value)) %>%
    #mutate(se = sd/sqrt(n), LCI = blankValue+qt(0.025,df=(n-1))*se, UCI = meanValue+qt(0.975, df=(n-1))*se) %>%
    unite(commonType, time, AA, sep="|")
  
    #write.csv(blankTbl, "blankTbl.csv")
    
    
    ## These tables are just for the graphs
    
  expBlankMasterTbl <- expTbl %>% inner_join(blankTbl, by = "commonType") %>%
      separate(commonType, into = c("time", "AA"), sep = "\\|") %>%
      transmute(time = time, AA = factor(AA, levels=AA_levels), replicate = replicate, aaRS = aaRS.x, clone = clone, codon = codon.x,
          rawGFS = GFS, rawOD600 = OD600, rawRFS = RFS, blankOD600 = blankOD600, blankRFS = blankRFS, blankGFS = blankGFS)
  
    expPlusBlankTbl <- expBlankMasterTbl %>% mutate(netGFS = rawGFS-blankGFS, netOD600=rawOD600-blankOD600, netRFS=rawRFS-blankRFS,
                  normGFS = netGFS/netOD600, normRFS = netRFS/netOD600,   ratio = normRFS/normGFS) %>%
        gather("measurementType", "value", -time, -clone, -aaRS, -codon, -AA, -replicate) %>%
        mutate(time=as.numeric(time), AA=factor(AA, levels=AA_levels), replicate=as.factor(replicate)) %>%
        group_by(time, aaRS, measurementType, codon)
    
  
  rawBlankNetTbl <- expBlankMasterTbl %>%
    mutate(netGFS = rawGFS-blankGFS, netOD600=rawOD600-blankOD600, netRFS=rawRFS-blankRFS) %>%
    gather("measurementType", "value", -time, -clone, -aaRS, -codon, -AA, -replicate) %>%
    mutate(time=as.numeric(time), AA=factor(AA, levels=AA_levels), replicate=as.factor(replicate)) %>%
    group_by(time, aaRS, measurementType, codon)

  ## This is a stopgap fix to get the linetype assignment to show present as solid/absent as dashed while plotting
  presentRawBlankNetTbl <- rawBlankNetTbl %>% filter(AA == "present") %>% mutate(AA=as.factor("present"))
  absentRawBlankNetTbl <- rawBlankNetTbl %>% filter(AA == "absent") %>% mutate(AA=as.factor("absent"))

  plotRawBlankNetTbl <- presentRawBlankNetTbl %>% full_join(absentRawBlankNetTbl) %>% mutate(AA=factor(AA, levels=AA_levels)) %>%
    group_by(time, aaRS, measurementType, codon, AA)
  plotRawBlankNetTbl$AA <- factor(plotRawBlankNetTbl$AA, levels=c("present","absent"))

  ggplot(plotRawBlankNetTbl, aes(time, value, group=interaction(replicate, AA, codon, clone), color=codon, linetype=AA)) +
    geom_line() +
    coord_cartesian() +
    facet_grid(measurementType ~ aaRS, scales = "free_y")
  ggsave(paste0(output_base_name, ".all.signals.pdf"), width=plot.width, height=7*plot.height)

  ## Return to calculations of ratios

  netOD600Tbl <- expBlankMasterTbl %>%
    mutate(netOD600=rawOD600-blankOD600) %>%
    mutate(time=as.numeric(time), AA=factor(AA, levels=AA_levels), replicate=as.factor(replicate)) %>%
    group_by(time, aaRS, codon)
  
  p1 = ggplot(netOD600Tbl, aes(time, netOD600, group=interaction(replicate, AA, codon, clone), color=codon, linetype=AA)) +
    geom_line() +
    coord_cartesian(ylim=c(0, 0.8), xlim=c(0,15)) #+
  # ggsave(paste0(output_base_name, ".net_OD600.pdf"), width=plot.width, height=plot.height)
  
  
  normGFSTbl <- expBlankMasterTbl %>%
    mutate(netGFS = rawGFS-blankGFS, netOD600=rawOD600-blankOD600, netRFS=rawRFS-blankRFS) %>%
    transmute(time = time, AA = AA, replicate = replicate, aaRS = aaRS, clone = clone, codon = codon,
              normGFS = netGFS/netOD600) %>%
    gather("measurementType", "value", -time, -clone, -aaRS, -codon, -AA, -replicate) %>%
    mutate(time=as.numeric(time), AA=factor(AA, levels=AA_levels), replicate=as.factor(replicate)) %>%
    group_by(time, aaRS, measurementType, codon)

  p3 = ggplot(normGFSTbl, aes(time, value, group=interaction(replicate, AA, codon, clone), color=codon, linetype=AA)) +
    geom_line() +
    coord_cartesian(ylim=c(-5000, 45000), xlim=c(0,15)) #+
  # ggsave(paste0(output_base_name, ".normalized_GFP.pdf"), width=plot.width, height=plot.height)

  normRFSTbl <- expBlankMasterTbl %>%
    mutate(netGFS = rawGFS-blankGFS, netOD600=rawOD600-blankOD600, netRFS=rawRFS-blankRFS) %>%
    transmute(time = time, AA = AA, replicate = replicate, aaRS = aaRS, clone = clone, codon = codon, normRFS = netRFS/netOD600) %>%
    gather("measurementType", "value", -time, -clone, -aaRS, -codon, -AA, -replicate) %>%
    mutate(time=as.numeric(time), AA=factor(AA, levels=AA_levels), replicate=as.factor(replicate)) %>%
    group_by(time, aaRS, measurementType, codon)

  p2 = ggplot(normRFSTbl, aes(time, value, group=interaction(replicate, AA, codon, clone), color=codon, linetype=AA)) +
    geom_line() +
    coord_cartesian(ylim=c(-2000, 8000), xlim=c(0,15)) #+
  # ggsave(paste0(output_base_name, ".normalized_RFP.pdf"), width=plot.width, height=plot.height)


  ## Construct final tables of normalized RFP and GFP signals

  amberNormRFSTbl <- normRFSTbl %>%
    filter(codon == "amber") %>%
    ungroup() %>%
    unite(commonCodon, time, aaRS, AA, clone, replicate, sep="|")
  #write.csv(amberTbl, "amberTbl.csv")

  amberNormGFSTbl <- normGFSTbl %>%
    filter(codon == "amber") %>%
    ungroup() %>%
    unite(commonCodon, time, aaRS, AA, clone, replicate, sep="|")
  #write.csv(amberTbl, "amberTbl.csv")

  tyrosineNormRFSTbl <- normRFSTbl %>%
    filter(codon == "tyrosine") %>%
    ungroup() %>%
    unite(commonCodon, time, aaRS, AA, clone, replicate, sep="|")
  #write.csv(amberTbl, "amberTbl.csv")

  tyrosineNormGFSTbl <- normGFSTbl %>%
    filter(codon == "tyrosine") %>%
    ungroup() %>%
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
    mutate(time=as.numeric(time), AA=factor(AA, levels=AA_levels), clone=as.factor(clone), replicate=as.factor(replicate))

  p4 = ggplot(amberToTyrosineRatioGFSToRFSRatioTbl, aes(time, value, group=interaction(replicate, AA, clone), linetype=AA, color=replicate)) +
    geom_line() +
    coord_cartesian(ylim=c(-0.2, 1.8), xlim=c(0,15)) +
    geom_hline(aes(yintercept=1), linetype="dashed")

  # ggsave(paste0(output_base_name, ".RRE.pdf"), width=plot.width, height=plot.height)

  ## Construct tables of RFP/GFP ratios

  absentAmberToTyrosineRatioGFSToRFSRatioTbl <- amberToTyrosineRatioGFSToRFSRatioTbl %>%
    filter(AA == "absent") %>%
    unite(commonAA, time, aaRS, clone, replicate, sep="|")

  presentAmberToTyrosineRatioGFSToRFSRatioTbl <- amberToTyrosineRatioGFSToRFSRatioTbl %>%
    filter(AA != "absent") %>%
    unite(commonAA, time, aaRS, clone, replicate, sep="|")

  absentToPresentRatioAmberToTyrosineRatioGFSToRFSRatioTbl <- absentAmberToTyrosineRatioGFSToRFSRatioTbl %>%
    left_join(presentAmberToTyrosineRatioGFSToRFSRatioTbl, by = "commonAA") %>%
    mutate(value = value.x / value.y) %>%
    select(commonAA, value) %>%
    separate(commonAA, into = c("time", "aaRS", "clone", "replicate"), sep = "\\|") %>%
    mutate(time=as.numeric(time), clone=as.factor(clone), replicate=as.factor(replicate))

  if (length(levels(absentToPresentRatioAmberToTyrosineRatioGFSToRFSRatioTbl$clone))>1) {
    
  p5 = ggplot(absentToPresentRatioAmberToTyrosineRatioGFSToRFSRatioTbl, aes(time, value, group=interaction(replicate, clone), color=clone)) +
    geom_line() +
    coord_cartesian(ylim=c(-0.2, 1.8), xlim=c(0,15)) +
    geom_hline(aes(yintercept=1), linetype="dashed")
  } else {
    p5 = ggplot(absentToPresentRatioAmberToTyrosineRatioGFSToRFSRatioTbl, aes(time, value, group=interaction(replicate, clone), color=replicate)) +
      geom_line() +
      coord_cartesian(ylim=c(-0.2, 1.8), xlim=c(0,15)) +
      geom_hline(aes(yintercept=1), linetype="dashed")
  }

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

  
  ## Write processed data files
  
  signalTbl = expBlankMasterTbl %>%
              mutate(netOD600=rawOD600-blankOD600) %>%
              mutate(netGFS=rawGFS-blankGFS) %>%
              mutate(netRFS=rawRFS-blankRFS) %>%
              mutate(normGFS = netGFS/netOD600) %>%
              mutate(normRFS = netRFS/netOD600)
  write.csv(signalTbl, paste0(output_base_name, ".signal.csv"))
  
  RRETbl = amberToTyrosineRatioGFSToRFSRatioTbl %>%
           rename(RRE=value)
  RRETbl.present = RRETbl %>% filter(AA=="present") %>% rename(RRE.present = RRE)
  RRETbl.absent = RRETbl %>% filter(AA=="absent") %>% rename(RRE.absent = RRE)
  RRETbl = RRETbl.absent
  RRETbl$RRE.present = RRETbl.present$RRE.present
  write.csv(RRETbl, paste0(output_base_name, ".RRE.csv"))
  
  R3CTbl = absentToPresentRatioAmberToTyrosineRatioGFSToRFSRatioTbl %>%
           rename(R3C=value)
  write.csv(R3CTbl, paste0(output_base_name, ".R3C.csv"))
}

