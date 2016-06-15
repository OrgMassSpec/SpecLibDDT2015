# Generate MSP text file:

# x <- MetaDataKelthane[, c("CompoundName", "SpectrumFilename")]
# names(x) <- c("compound", "filename")
# x <- x[!duplicated(x), ]
# 
# # Remove K4 and K5 spectra since these are not recorded in the metadata (incorrect assignment later removed).
# y <- droplevels(SpecData[SpecData$filename %in% c("K1","K2","K3","K6","K7","K8","K9","K10","K11"), ])
# WriteMspFile(y, x, filename = "SpecLibKelthane2015.txt", 
#              comment = "Newly Identified DDT-Related Compounds Accumulating in Southern California Bottlenose Dolphins, Mackintosh et al., 2015")



LibraryReportKelthane <- function(spectra = SpecData,
  metadata = MetaDataKelthane,
  pdfFile = "SpecLibKelthane2015.pdf",
  pdfTitle = "SpecLibKelthane2015 Library",
  xMin = 40) {
  
  pdf(file = pdfFile, width = 10.5, height = 8, title = pdfTitle, paper = "usr")
  
  
  #---------- Title Page ----------
  
  grid.newpage()
  
  pushViewport(viewport(layout = grid.layout(nrow = 4, ncol = 1, heights = unit(c(0.1, 0.2, 0.5, 0.2), "npc"))))
  
  pushViewport(viewport(layout.pos.row = 1))
  
  grid.text("SpecLibKelthane2015 Mass Spectral Library", y = 0.5, gp = gpar(cex = 1.25))
  grid.lines(x = unit(c(0,1), "npc"), y = unit(c(0,0), "npc"))
  
  popViewport()
  
  pushViewport(viewport(layout.pos.row = 2))
  
  grid.text("Kelthane-Dicofol Technical Mixture (500 ppm)", y = 0.5, gp = gpar(cex = 2))
  
  popViewport()
  
  pushViewport(viewport(layout.pos.row = 3))
  
  grid.text("Newly Identified DDT-Related Compounds Accumulating in\nSouthern California Bottlenose Dolphins", y = 0.9, gp = gpar(cex = 1.25))
  grid.text("Susan A. Mackintosh, Nathan G. Dodder, Nellie J. Shaul, Lihini I. Aluwihare, Keith A. Maruya,\nSusan J. Chivers, Kerri Danil, David W. Weller, Eunha Hoh", y = 0.6, gp = gpar(cex = 1.25))
  grid.text("Web Reference: http://OrgMassSpec.github.io", y = 0.3, gp = gpar(cex = 1.25))
  
  popViewport()
  
  pushViewport(viewport(layout.pos.row = 4))
  
  session.info <- sessionInfo()
  
  grid.text(paste("Prepared:", Sys.time()), y = 0.8)
  grid.text(paste("SpecLibDDT version", session.info$otherPkgs$SpecLibDolphin2014$Version), y = 0.65)
  grid.text(paste("OrgMassSpecR version", session.info$otherPkgs$OrgMassSpecR$Version), y = 0.5)
  grid.text(session.info$R.version$version.string, y = 0.35)
  
  popViewport()
  
  # Prepare metadata
  
  uniqueCompounds <- as.character(unique(metadata$CompoundName))
  
  # DrawSpectrum <- function(currentFilename) {
  
  for(i in 1:length(uniqueCompounds)) {
    
    currentMetadata <- metadata[metadata$CompoundName == uniqueCompounds[i], ][1, ]
    
    # currentMetadata <- metadata[metadata$CompoundName == currentCompound, ]
    
    currentSpectrum <- spectra[spectra$filename == as.character(currentMetadata$SpectrumFilename), ]
    
    isomerMetadata <- metadata[metadata$CompoundName == uniqueCompounds[i], 
      c("Isomer", "DetectedinPacificDolphin", 
        "TypicallyMonitored", "RetentionTime1D", "RetentionTime2D")]
    
    names(isomerMetadata) <- c("Isomer", "InDolphin", "Monitored", "RT1D", "RT2D")
    
    row.names(isomerMetadata) <- seq(1:nrow(isomerMetadata))
    
    message("Making spectrum for ", currentMetadata$Compound)
    
    page <- (1:length(uniqueCompounds))[i] + 1
    
    grid.newpage()
    
    # Write information at top of page
    
    pushViewport(viewport(layout = grid.layout(nrow = 5, 
      ncol = 1, 
      heights = unit(c(0.05, 0.075, 0.45, 0.375, 0.05), "npc"))))
    
    pushViewport(viewport(layout.pos.row = 1))
    
    grid.text(paste("Name:", currentMetadata$CompoundName), 
      x = 0, y = 0.5, just = c("left", "center"), gp = gpar(cex = 1.25))
    
    grid.text(paste("Technical Mixture:", currentMetadata$TechnicalMixture),
      x = 0.67, y = 0.5, just = c("left", "center"), gp = gpar(cex = 1.25))
    
    grid.lines(x = unit(c(0, 1), "npc"), y = unit(c(0,0), "npc"))
    
    popViewport()
    
    
    #---------- Compound Info ----------
    
    pushViewport(viewport(layout.pos.row = 2))
    
    grid.text(paste("Compound Class: ", currentMetadata$Class), x = 0, y = 0.75, hjust = 0) 
    grid.text(paste("Instrument:", currentMetadata$Instrument), x = 0, y = 0.375, hjust = 0)    
    
    if(!is.na(currentMetadata$Formula)) {
      grid.text(paste("Elemental Formula:", currentMetadata$Formula), x = 0.67, y = 0.75, hjust = 0)
    } else {
      grid.text(paste("Elemental Formula:"), x = 0.67, y = 0.75, hjust = 0)
    }
    
    if(!is.na(currentMetadata$Comment)) {
      grid.text(paste("Comment:", currentMetadata$Comment), x = 0, y = 0, hjust = 0)
    } else {
      grid.text(paste("Comment:"), x = 0, y = 0, hjust = 0)
    }
    
    popViewport()
    
    
    #---------- Draw spectrum ----------
    
    pushViewport(viewport(layout.pos.row = 3))
    
    currentSpectrum$percent.intensity <- with(currentSpectrum, intensity / max(intensity) * 100)
    
    # Calculate molecular weight to set x-axis upper limit
    
    if(!is.na(currentMetadata$Formula)) {
      
      mw <- MolecularWeight(formula = ListFormula(currentMetadata$Formula))
      xMax <- mw + (mw * 0.03)
      
    } else {
      
      m <- max(currentSpectrum$mz) 
      xMax <- m + (m * 0.03)
      
    }
    
    plot.data <-currentSpectrum[currentSpectrum$mz >= xMin & currentSpectrum$mz <= xMax, ]
    
    pushViewport(plotViewport(c(3.75, 3.5, 1.5, 1)))
    pushViewport(dataViewport(xscale = c(xMin, xMax),
      yscale = c(0, 110)))
    
    grid.rect()
    p.ticks <- pretty(plot.data$mz, n = 10)
    x.ticks <- p.ticks[p.ticks >= xMin & p.ticks <= xMax]
    grid.xaxis(at = x.ticks)
    grid.yaxis(at = c(0, 25, 50, 75, 100))
    
    grid.segments(plot.data$mz,
      plot.data$percent.intensity,
      plot.data$mz,
      rep(0, length(plot.data$intensity)),
      default.units = "native",
      gp = gpar(lwd = 0.75))
    
    ## print m/z values in plot
    
    display.values <- plot.data$mz[plot.data$display == TRUE]
    if (length(display.values) > 0) {
      grid.text(display.values,
        x = display.values,
        y = plot.data$percent.intensity[plot.data$display == TRUE] + 5,
        default.units = "native",
        gp = gpar(col = "blue"))
    }
    
    grid.text("intensity (%)", x = unit(-3.2, "lines"), rot = 90)
    grid.text("m/z", y = unit(-2.5, "lines"))
    
    # Automatic display of values every 25 units
    
    plotSegments <- cut(plot.data$mz, breaks = xMax/25)
    
    plot.data <- cbind(plot.data, plotSegments)
    
    for(i in 1:length(unique(plotSegments))) {
      
      segmentTmp <- plot.data[plot.data$plotSegments == unique(plotSegments)[i], ]
      
      maxPeak <- segmentTmp[segmentTmp$percent.intensity == max(segmentTmp$percent.intensity, na.rm = TRUE), ]
      
      if(maxPeak$percent.intensity[1] >= 5) {
        
        grid.text(maxPeak$mz,
          x = maxPeak$mz,
          y = maxPeak$percent.intensity + 5,
          default.units = "native",
          gp = gpar(col = "blue"))
      }
      
    }
    
    popViewport(3)
    
    # Define area below spectrum
    
    pushViewport(viewport(layout.pos.row = 4))
    
    pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 5, 
      widths = unit(c(0.05, 0.55, 0.05, 0.3, 0.05), "npc"))))
    
    
    #---------- Isomer Table ----------
    
    pushViewport(viewport(layout.pos.col = 2))
    
    grid.rect(x = unit(0, "npc"), height = unit(1, "npc"), width = unit(1, "npc"), gp = gpar(col = "black"), just = "left")
    
    grid.text(paste(capture.output(isomerMetadata), collapse = "\n"), gp = gpar(cex = 0.75,
      fontfamily = "mono"), vjust = "bottom")
    
    grid.text("Isomer Information", x = unit(0.02, "npc"), rot = 90)
    
    popViewport(1)
    
    
    #---------- Write fragment ion identifications ----------
    
    pushViewport(viewport(layout.pos.col = 4))
    
    grid.rect(x = unit(0, "npc"), height = unit(1, "npc"), width = unit(1, "npc"), gp = gpar(col = "black"), just = "left")
    
    grid.text("m/z [Fragment]", x = unit(0.075, "npc"), y = unit(0.9, "npc"), gp = gpar(col = "blue"), just = "left")
    
    grid.lines(x = unit(c(0, 1), "npc"), y = unit(0.8, "npc"), gp = gpar(col = "black"))
    
    if(!is.na(currentMetadata$FragmentIdentification)) {
      
      fragmentText <- gsub(pattern = "; ", replacement = "\n", x = currentMetadata$FragmentIdentification, fixed = TRUE)
      
      grid.text(fragmentText, x = unit(0.075, "npc"), y = unit(0.7, "npc"), gp = gpar(col = "black"), just = c("left", "top"))
      
    }
    
    popViewport(3) 
    
    # Write filename in corner
    
    pushViewport(viewport(layout.pos.row = 5))
    
    grid.text(paste("Filename: ", currentMetadata$SpectrumFilename, ", Page: ", page, sep = ""),
      x = 1,
      just = c("right", "top"),
      gp = gpar(col = "dark grey"))
    
    grid.text("InDolphin = In Pacific Dolphin (Shaul, et al., 2015); Monitored = Typically Monitored; RT1D and RT2D = 1st and 2nd Dimension Retention Times",
      x = 0,
      just = c("left", "top"),
      gp = gpar(col = "dark grey", cex = 0.75))
    
  }
  
  #   sapply(metadata$SpectrumFilename, 
  #          DrawSpectrum)
  
  graphics.off()
  
}

