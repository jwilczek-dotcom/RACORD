rm(list=ls(all=TRUE))

#..............................................................................................
# R A C O R D  (Recherche Automatique de la Ceramique par ORDinateur) -  Server
#..............................................................................................
# Last update: 2021/01/28

# SET PATH TO APPLICATION FOLDER:
WD <- c("C:/RACORD/")
dir_profiles <<- c("C:/RACORD/profiles/")

# ..............................................................................................
# Notes:
#   RTC:  no weights are used in calculation of values which are shown as titles in figures plotted in the bottom
# ..............................................................................................
# DO NOT ALTER
# libraries
library(shiny)
library(sp)       # surf, rechant
library(MASS)     # Bezier
library(kmlShape) # DouglasPeuckerNbPoints


# inits
def.par <- par(no.readonly = TRUE)
setwd(WD)
source(paste(WD,"functions.r",sep=""))
set.seed(1234)
options(shiny.maxRequestSize=200*1024^2)
dur <<- 5         # duration of the message in Shiny
S.nb <<- 0
FileName <<-NULL
dataset <<- NULL
sources <<- NULL
targets <<- NULL
standardize_radius <<- F
pts.nb <<- 200
P.TARGET.O <<- NULL
P.SOURCE.O <<- NULL
TARG <<- NULL
SOUR <<- NULL
bla <<- NULL
PtsNb <<- 0
Cex <<- 0.5
shoow <<- T
sS <<- NULL
sT <<- NULL
SD <<- 10
SelNormalisation_RUNNED <<- ""
SelSampling_RUNNED <<- ""
PtsNb_RUNNED <<- 0
PtsDist_RUNNED <<- 0
PtsRechant4Breaking <<- 1000


#..............................................................................................
shinyServer(function(input, output, session) {
  
  dataInputDatasetTargets <- reactive({
    inFile <- input$fileTargets
    if (is.null(inFile)) return(NULL)
    FileName <<- inFile$name
    FileExt <- substr(FileName, nchar(FileName)-2, nchar(FileName))
    if (FileExt!="txt") { showNotification(paste("Please load TXT file"), duration = dur); return() }
    file <- inFile$datapath
    filen <- gsub("[\\]", "/", file)
    filen <- paste(substr(filen, 1, nchar(filen)-1), inFile$name, sep="")
    file.rename(file, filen)
    targets <- read.table(filen, header=T, sep=";", dec = ",", stringsAsFactors=F)
    indiv <<- as.matrix(targets$ID)
    return(targets)
  })
  
  dataInputDatasetSources <- reactive({
    inFile <- input$fileSources
    if (is.null(inFile)) return(NULL)
    FileName <<- inFile$name
    FileExt <- substr(FileName, nchar(FileName)-2, nchar(FileName))
    if (FileExt!="txt") { showNotification(paste("Please load TXT file"), duration = dur); return() }
    file <- inFile$datapath
    filen <- gsub("[\\]", "/", file)
    filen <- paste(substr(filen, 1, nchar(filen)-1), inFile$name, sep="")
    file.rename(file, filen)
    sources <- read.table(filen, header=T, sep=";", dec = ",", stringsAsFactors=F)
    indiv.sources <<- as.matrix(sources$ID)
    return(sources)
  })
  
  observeEvent(input$fileTargets, {
    targets <<- dataInputDatasetTargets()
    output$dataTableTargets <- DT::renderDataTable(
      DT::datatable(data <- targets, rownames=F))
    P.TARGET.O <<- funLoadTargetProfiles()
  })
  
  observeEvent(input$fileSources, {
    sources <<- dataInputDatasetSources()
    output$dataTableSources <- DT::renderDataTable(
      DT::datatable(data <- sources, rownames=F)
    )
    P.SOURCE.O <<- funLoadSourceProfiles()
  })
  
  output$dataTableTargets = DT::renderDataTable(targets, server = FALSE)
  
  output$dataTableSources = DT::renderDataTable(sources, server = FALSE)
  
  funLoadTargetProfiles <- function(){
    P.TARGET.O <<- vector("list", length(indiv))
    proc = 1/length(P.TARGET.O)
    withProgress(message = 'Loading targets', value = 0, {
      for (i in 1:length(indiv)) {
        ID <- indiv[i]
        xy <- read.table(paste(dir_profiles,ID,".txt", sep=""), sep=";", stringsAsFactors=F)
        xy <- cbind(xy[,1],-xy[,2])
        P.TARGET.O[[i]] <- profil.xy(xy, plots=F)
        setProgress(proc*i)
      }
    })
    
    return(P.TARGET.O)
  }
  
  funLoadSourceProfiles <- function(){
    P.SOURCE.O <<- vector("list", length(indiv.sources))
    for (i in 1:length(indiv.sources)) {
      ID <- indiv.sources[i]
      xy <- read.table(paste(dir_profiles,ID,".txt", sep=""), sep=";", stringsAsFactors=F)
      xy <- cbind(xy[,1],-xy[,2])
      P.SOURCE.O[[i]] <- profil.xy(xy, plots=F)
    }
    return(P.SOURCE.O)
  }
  
  
  funSelNormalisationAndSampling <- function() {}    ### NOT USED, SERVES FOR CODE SEPARATION
  
  observeEvent(c(input$selSampling,input$selNormalisation),{
    if(input$selSampling=="Sampling_equidistant"){
      if(input$selNormalisation=="Normalisation_radius"){
        updateSliderInput(session, 'sldPtsDistance', label='Distance between points', min=0.05, max=1, step=0.05, value = 0.05)
        updateSliderInput(session, 'sldICP_rigid_z', label='z', min=0.05, max=1, step=0.05, value = 0.5)
        updateSliderInput(session, 'sldICP_R_Attribution', label=NULL, value=c(-1,1), min=-10, max=10, step=0.1)
      }
      else if(input$selNormalisation=="Normalisation_rim"){
        updateSliderInput(session, 'sldPtsDistance', label='Distance between points', min=0.005, max=0.3, step=0.01, value = 0.01)
        updateSliderInput(session, 'sldICP_rigid_z', label='z (in cm)', min=0.005, max=1, step=0.01, value = 0.1)
        updateSliderInput(session, 'sldICP_R_Attribution', label=NULL, value=c(-0.01,0.01), min=-1, max=1, step=0.01)
      }
    }
  })
  
  
  plotALL <- function(){}    ### NOT USED, SERVES FOR CODE SEPARATION
  
  observeEvent(c(input$dataTableTargets_rows_selected,
                 input$dataTableSources_rows_selected,
                 input$selNormalisation,
                 input$selSampling,
                 input$sldNbOfPoints,
                 input$sldPtsDistance,
                 input$selShow,
                 
                 input$cbxAdvancedDCT,
                 input$sldDCT_nbHarm,
                 
                 input$cbxAdvancedRDP,
                 input$sldRDP_nbPts,
                 
                 input$cbxAdvancedICP_Rigid,
                 input$sldICP_rigid_z,
                 
                 input$cbxAdvancedICP,
                 input$sldICP_r,
                 input$sldICP_theta,
                 input$sldICP_size
  ), {
    
    sT <<- input$dataTableTargets_rows_selected
    sS <<- input$dataTableSources_rows_selected
    
    SelShow <<- input$selShow
    SelMethod <<- input$selMethod
    
    if (SelShow=='selShowAnalyse') {
      showAnalyse()
    }
    
    else if (SelShow=='selShowResults') {
      # Test if source is selected
      if (is.null(sS)) { }
      else if (length(sS)>1) {
        showNotification(paste("Select just one source"), duration = dur)
      }
      else if (length(sS)==1) {
        # Test if target is selected
        if (is.null(sT)) {
          output$plotTemp1 <- renderPlot({})
          output$plotTemp2 <- renderPlot({})
          output$plotTemp3 <- renderPlot({})
          output$plotTemp4 <- renderPlot({})
        }
        # if more targets are selected show them all
        else if (length(sT)>1) {
          showAnalyse()
          output$plotTemp3 <- renderPlot({})
          output$plotTemp4 <- renderPlot({})
        }
        # if only one target is selected not, show details
        else if (length(sT)==1) {
          if (SelMethod=='') { }
          else if (SelMethod=='selBEZ') { plotResult_BEZ() }
          else if (SelMethod=='selDCT') { plotResult_DCT() }
          else if (SelMethod=='selRDP') { plotResult_RDP() }
          else if (SelMethod=='selRTC') { plotResult_RTC() }
          else if (SelMethod=='selICP_Rim') { plotResult_ICP_Rim() }
          else if (SelMethod=='selICP_Rigid') { plotResult_ICP_Rigid() }
          else if (SelMethod=='selICP') { plotResult_ICP() }
        }
      }
    }
  })
  
  
  observeEvent(input$sldBreaks, {
    MaxVal = floor(100/input$sldBreaks)
    updateSliderInput(session, 'sldMinPerc', label='Minimum % size of one fragment', min=5, max=MaxVal, step=1, value = floor(MaxVal*2/3))
  })
  
  plotResult_BEZ <- function() {
    if (is.null(P.SOURCE.O[[sS]]$res_BEZ)) { notifyNotAnalysed() }
    else {
      x_d <<- P.SOURCE.O[[sS]]$res_BEZ[sT]
      res_log <<- P.SOURCE.O[[sS]]$res_BEZ_log
      
      x_SelNormalisation <<- res_log$SelNormalisation
      x_SelSampling      <<- res_log$SelSampling
      x_PtsNb            <<- as.numeric(res_log$PtsNb)
      x_PtsDist          <<- as.numeric(res_log$PtsDist)

      # 0) Preparations
      SOURCE <- vector("list", length(sS))
      for(i in 1:length(SOURCE)){ SOURCE[[i]] <- P.SOURCE.O[[sS]] }
      bla <- funFun(LIST = SOURCE,
                    selNormalisation = x_SelNormalisation,
                    selSampling = x_SelSampling,
                    Radius = sources$Diameter[[sS]]/2,
                    PtsNb = x_PtsNb,
                    PtsDist = x_PtsDist)
      SOURCE <- funFunDct(LIST = bla, PtsDist = x_PtsDist)
      
      TARGET <- vector("list", length(sT))
      for(i in 1:length(TARGET)){ TARGET[[i]] <- P.TARGET.O[[sT]] }
      bla <- funFun(LIST = TARGET,
                    selNormalisation = x_SelNormalisation,
                    selSampling = x_SelSampling,
                    Radius = targets$Diameter[[sT]]/2,
                    PtsNb = x_PtsNb,
                    PtsDist = x_PtsDist)
      TARGET <- funFunDct(LIST = bla, PtsDist = x_PtsDist)
      
      P.source <- SOURCE[[1]]
      P.target <- TARGET[[1]]
      
      if ( min(P.source$s)<min(P.target$s) && max(P.source$s)>max(P.target$s)) {
        showNotification("No way, target is not long enough...", duration = dur)
      }
      else {
        s_min <- which(P.target$s==min(P.source$s))
        s_max <- which(P.target$s==max(P.source$s))
        sel <- s_min:s_max
        
        xl <- min(c(P.target$profil[,1],P.source$profil[,1]))
        yl <- min(c(P.target$profil[,2],P.source$profil[,2]))
        
        # plot the plot
        output$plotTemp1 <- renderPlot({
          plot(0, asp=1, type="n", xlab="r", ylab="z", xlim=c(xl,0), ylim=c(yl,0), main=paste("d = ", round(x_d,2)))
          polygon(P.target$profil, col=rgb(0,0,0,0.2), border=rgb(0,0,0,0.2))
          polygon(P.source$profil, col=rgb(1,0,0,0.2), border=rgb(1,0,0,0.2))
          lines(P.target$extint_seg[sel,], col=rgb(0,0,0,0.2), lwd=2)
          lines(P.source$extint_seg, col=rgb(1,0,0,0.2), lwd=2)
        })
        output$plotTemp2 <- renderPlot({ })
        output$plotTemp3 <- renderPlot({ })
        output$plotTemp4 <- renderPlot({
          compare.bez(P.source, P.target, method="rmsd", plots=T)
        })
      }
    }
  }
  
  plotResult_DCT <- function() {
    if (is.null(P.SOURCE.O[[sS]]$res_DCT)) { notifyNotAnalysed() }
    else {
      x_d <<- P.SOURCE.O[[sS]]$res_DCT[sT]
      res_log <<- P.SOURCE.O[[sS]]$res_DCT_log
      
      x_SelNormalisation <<- res_log$SelNormalisation
      x_SelSampling      <<- res_log$SelSampling
      x_PtsNb            <<- as.numeric(res_log$PtsNb)
      x_PtsDist          <<- as.numeric(res_log$PtsDist)
      x_SldDCT_nbHarm    <<- as.numeric(res_log$SldDCT_nbHarm)
      
      # 0) Preparations
      SOURCE <- vector("list", length(sS))
      for(i in 1:length(SOURCE)){ SOURCE[[i]] <- P.SOURCE.O[[sS]] }
      bla <- funFun(LIST = SOURCE,
                    selNormalisation = x_SelNormalisation,
                    selSampling = x_SelSampling,
                    Radius = sources$Diameter[[sS]]/2,
                    PtsNb = x_PtsNb,
                    PtsDist = x_PtsDist)
      SOURCE <- funFunDct(LIST = bla, PtsDist = x_PtsDist)
      
      TARGET <- vector("list", length(sT))
      for(i in 1:length(TARGET)){ TARGET[[i]] <- P.TARGET.O[[sT]] }
      bla <- funFun(LIST = TARGET,
                    selNormalisation = x_SelNormalisation,
                    selSampling = x_SelSampling,
                    Radius = targets$Diameter[[sT]]/2,
                    PtsNb = x_PtsNb,
                    PtsDist = x_PtsDist)
      TARGET <- funFunDct(LIST = bla, PtsDist = x_PtsDist)
      
      P.source <- SOURCE[[1]]
      P.target <- TARGET[[1]]
      
      if ( min(P.source$s)<min(P.target$s) && max(P.source$s)>max(P.target$s)) {
        showNotification("No way, target is not long enough...", duration = dur)
      }
      else {
        s_min <- which(P.target$s==min(P.source$s))
        s_max <- which(P.target$s==max(P.source$s))
        sel <- s_min:s_max
        
        xl <- min(c(P.target$profil[,1],P.source$profil[,1]))
        yl <- min(c(P.target$profil[,2],P.source$profil[,2]))
        
        # plot the plot
        output$plotTemp1 <- renderPlot({
          plot(0, asp=1, type="n", xlab="r", ylab="z", xlim=c(xl,0), ylim=c(yl,0), main=paste("d = ", round(x_d,2)))
          polygon(P.target$profil, col=rgb(0,0,0,0.2), border=rgb(0,0,0,0.2))
          polygon(P.source$profil, col=rgb(1,0,0,0.2), border=rgb(1,0,0,0.2))
          lines(P.target$extint_seg[sel,], col=rgb(0,0,0,0.2), lwd=2)
          lines(P.source$extint_seg, col=rgb(1,0,0,0.2), lwd=2)
        })
        output$plotTemp2 <- renderPlot({ })
        output$plotTemp3 <- renderPlot({ })
        output$plotTemp4 <- renderPlot({
          compare.dct(P.source, P.target, nharm=x_SldDCT_nbHarm, method="rmsd", plots=T)
        })
      }
    }
  }
  
  plotResult_RDP <- function() {
    if (is.null(P.SOURCE.O[[sS]]$res_RDP)) { notifyNotAnalysed() }
    else {
      x_d <<- P.SOURCE.O[[sS]]$res_RDP[sT]
      res_log <<- P.SOURCE.O[[sS]]$res_RDP_log
      
      x_SelNormalisation <<- res_log$SelNormalisation
      x_SelSampling      <<- res_log$SelSampling
      x_PtsNb            <<- as.numeric(res_log$PtsNb)
      x_PtsDist          <<- as.numeric(res_log$PtsDist)
      x_SldRDP_nbPts     <<- res_log$SldRDP_nbPts
      
      # 0) Preparations
      SOURCE <- vector("list", length(sS))
      for(i in 1:length(SOURCE)){ SOURCE[[i]] <- P.SOURCE.O[[sS]] }
      bla <- funFun(LIST = SOURCE,
                    selNormalisation = x_SelNormalisation,
                    selSampling = x_SelSampling,
                    Radius = sources$Diameter[[sS]]/2,
                    PtsNb = x_PtsNb,
                    PtsDist = x_PtsDist)
      SOURCE <- funFunDct(LIST = bla, PtsDist = x_PtsDist)
      
      TARGET <- vector("list", length(sT))
      for(i in 1:length(TARGET)){ TARGET[[i]] <- P.TARGET.O[[sT]] }
      bla <- funFun(LIST = TARGET,
                    selNormalisation = x_SelNormalisation,
                    selSampling = x_SelSampling,
                    Radius = targets$Diameter[[sT]]/2,
                    PtsNb = x_PtsNb,
                    PtsDist = x_PtsDist)
      TARGET <- funFunDct(LIST = bla, PtsDist = x_PtsDist)
      
      P.source <- SOURCE[[1]]
      P.target <- TARGET[[1]]
      
      if ( min(P.source$s)<min(P.target$s) && max(P.source$s)>max(P.target$s)) {
        showNotification("No way, target is not long enough...", duration = dur)
      }
      else {
        s_min <- which(P.target$s==min(P.source$s))
        s_max <- which(P.target$s==max(P.source$s))
        sel <- s_min:s_max
        
        xl <- min(c(P.target$profil[,1],P.source$profil[,1]))
        yl <- min(c(P.target$profil[,2],P.source$profil[,2]))
        
        # plot the plot
        output$plotTemp1 <- renderPlot({
          plot(0, asp=1, type="n", xlab="r", ylab="z", xlim=c(xl,0), ylim=c(yl,0), main=paste("d = ", round(x_d,2)))
          polygon(P.target$profil, col=rgb(0,0,0,0.2), border=rgb(0,0,0,0.2))
          polygon(P.source$profil, col=rgb(1,0,0,0.2), border=rgb(1,0,0,0.2))
          lines(P.target$extint_seg[sel,], col=rgb(0,0,0,0.2), lwd=2)
          lines(P.source$extint_seg, col=rgb(1,0,0,0.2), lwd=2)
        })
        output$plotTemp2 <- renderPlot({ })
        output$plotTemp3 <- renderPlot({ })
        output$plotTemp4 <- renderPlot({
          compare.rdp(P.source, P.target, method="rmsd", nbPoints=x_SldRDP_nbPts, plots=T)
        })
      }
    }
  }
  
  plotResult_RTC <- function() {
    if (is.null(P.SOURCE.O[[sS]]$res_RTC)) { notifyNotAnalysed() }
    else {
      x_d <<- P.SOURCE.O[[sS]]$res_RTC[sT]
      res_log <<- P.SOURCE.O[[sS]]$res_RTC_log
      
      x_SelNormalisation <<- res_log$SelNormalisation
      x_SelSampling      <<- res_log$SelSampling
      x_PtsNb            <<- as.numeric(res_log$PtsNb)
      x_PtsDist          <<- as.numeric(res_log$PtsDist)
      x_w                <<- as.numeric(res_log$w)
      # x_SldRTC_W         <<- as.numeric(res_log$SldRTC_W)

      # 0) Preparations
      SOURCE <- vector("list", length(sS))
      for(i in 1:length(SOURCE)){ SOURCE[[i]] <- P.SOURCE.O[[sS]] }
      bla <- funFun(LIST = SOURCE,
                    selNormalisation = x_SelNormalisation,
                    selSampling = x_SelSampling,
                    Radius = sources$Diameter[[sS]]/2,
                    PtsNb = x_PtsNb,
                    PtsDist = x_PtsDist)
      SOURCE <- funFunRtc(LIST = bla, PtsDist = x_PtsDist)

      TARGET <- vector("list", length(sT))
      for(i in 1:length(TARGET)){ TARGET[[i]] <- P.TARGET.O[[sT]] }
      bla <- funFun(LIST = TARGET,
                    selNormalisation = x_SelNormalisation,
                    selSampling = x_SelSampling,
                    Radius = targets$Diameter[[sT]]/2,
                    PtsNb = x_PtsNb,
                    PtsDist = x_PtsDist)
      TARGET <- funFunRtc(LIST = bla, PtsDist = x_PtsDist)
      
      P.source <- SOURCE[[1]]
      P.target <- TARGET[[1]]
      
      if ( min(P.source$s)<min(P.target$s) && max(P.source$s)>max(P.target$s)) {
        showNotification("No way, target is not long enough...", duration = dur)
      }
      else {
        s_min <- which(P.target$s==min(P.source$s))
        s_max <- which(P.target$s==max(P.source$s))
        sel <- s_min:s_max
        
        xl <- min(c(P.target$profil[,1],P.source$profil[,1]))
        yl <- min(c(P.target$profil[,2],P.source$profil[,2]))
        
        # plot the plot
        output$plotTemp1 <- renderPlot({
          plot(0, asp=1, type="n", xlab="r", ylab="z", xlim=c(xl,0), ylim=c(yl,0), main=paste("d = ", round(x_d,2)))
          polygon(P.target$profil, col=rgb(0,0,0,0.2), border=rgb(0,0,0,0.2))
          polygon(P.source$profil, col=rgb(1,0,0,0.2), border=rgb(1,0,0,0.2))
          lines(P.target$extint_seg[sel,], col=rgb(0,0,0,0.2), lwd=2)
          lines(P.source$extint_seg, col=rgb(1,0,0,0.2), lwd=2)
        })
        output$plotTemp2 <- renderPlot({
          compare.rtc(P.source, P.target, w=x_w, selfun="Ts", method="rmsd.RTC", plots=T)  
        })
        output$plotTemp3 <- renderPlot({
          compare.rtc(P.source, P.target, w=x_w, selfun="Rs", method="rmsd.RTC", plots=T)
        })
        output$plotTemp4 <- renderPlot({
          compare.rtc(P.source, P.target, w=x_w, selfun="Ks", method="rmsd.RTC", plots=T)
        })
      }
    }
  }
  
  
  
  
  plotResult_ICP_Rim <- function() {
    if (is.null(P.SOURCE.O[[sS]]$res_ICP_rim)) { notifyNotAnalysed() }
    else {
      x_d <<- P.SOURCE.O[[sS]]$res_ICP_rim[sT]
      res_log <<- P.SOURCE.O[[sS]]$res_ICP_rim_log

      x_par <<- c(0,0,1,0)
      
      x_SelNormalisation <<- res_log$SelNormalisation
      x_SelSampling      <<- res_log$SelSampling
      x_PtsNb            <<- as.numeric(res_log$PtsNb)
      x_PtsDist          <<- as.numeric(res_log$PtsDist)
      
      # 0) Preparations
      SOURCE <- vector("list", length(sS))
      for(i in 1:length(SOURCE)){ SOURCE[[i]] <- P.SOURCE.O[[sS]] }
      bla <- funFun(LIST = SOURCE,
                    selNormalisation = x_SelNormalisation,
                    selSampling = x_SelSampling,
                    Radius = sources$Diameter[[sS]]/2,
                    PtsNb = x_PtsNb,
                    PtsDist = x_PtsDist)
      if (SelSampling=="Sampling_equidistant") { for(i in 1:length(bla)) { bla[[i]]$extern <- bla[[i]]$extern_seg; bla[[i]]$intern <- bla[[i]]$intern_seg } }
      SOURCE <- bla
      
      TARGET <- vector("list", length(sT))
      for(i in 1:length(TARGET)){ TARGET[[i]] <- P.TARGET.O[[sT]] }
      bla <- funFun(LIST = TARGET,
                    selNormalisation = x_SelNormalisation,
                    selSampling = x_SelSampling,
                    Radius = targets$Diameter[[sT]]/2,
                    PtsNb = x_PtsNb,
                    PtsDist = x_PtsDist)
      if (SelSampling=="Sampling_equidistant") { for(i in 1:length(bla)) { bla[[i]]$extern <- bla[[i]]$extern_seg; bla[[i]]$intern <- bla[[i]]$intern_seg } }
      TARGET <- bla
      
      P.source <<- SOURCE[[1]]
      P.target <<- TARGET[[1]]
      
      # Source (adjustement with the rim)
      adj <- max(P.source$profil[,2])
      P.source$profil[,2] <- P.source$profil[,2]-adj
      P.source$extern[,2] <- P.source$extern[,2]-adj
      P.source$intern[,2] <- P.source$intern[,2]-adj
      
      # plot the plot
      output$plotTemp1 <- renderPlot({
        icp.p.f(P.source=P.source, P.target = P.target, is.rim=F, par = x_par, projection = F, minfun = "rmsd",
                 method=c("rtsz"), EI="EI", verbose = F, plots=T, plots.full=F)
      })
      output$plotTemp2 <- renderPlot({ })
      output$plotTemp3 <- renderPlot({ })
      output$plotTemp4 <- renderPlot({
        icp.p.all(P.source=P.source, P.target=P.target, is.rim=F, par = x_par, projection = F, minfun = "rmsd",
                 method=c("rtsz"), EI="EI", verbose = F, plots=T, plots.full=T)
      })
    }
  }
  
  plotResult_ICP_Rigid <- function() {
    if (is.null(P.SOURCE.O[[sS]]$res_ICP_rigid)) {  notifyNotAnalysed() }
    else {
      x_d <<- P.SOURCE.O[[sS]]$res_ICP_rigid[sT,5]
      res_log <<- P.SOURCE.O[[sS]]$res_ICP_rigid_log
      
      x_par <<- P.SOURCE.O[[sS]]$res_ICP_rigid[sT,1:4]
      
      x_SelNormalisation <<- res_log$SelNormalisation
      x_SelSampling      <<- res_log$SelSampling
      x_PtsNb            <<- as.numeric(res_log$PtsNb)
      x_PtsDist          <<- as.numeric(res_log$PtsDist)
      
      # 0) Preparations
      SOURCE <- vector("list", length(sS))
      for(i in 1:length(SOURCE)){ SOURCE[[i]] <- P.SOURCE.O[[sS]] }
      bla <- funFun(LIST = SOURCE,
                    selNormalisation = x_SelNormalisation,
                    selSampling = x_SelSampling,
                    Radius = sources$Diameter[[sS]]/2,
                    PtsNb = x_PtsNb,
                    PtsDist = x_PtsDist)
      if (SelSampling=="Sampling_equidistant") { for(i in 1:length(bla)) { bla[[i]]$extern <- bla[[i]]$extern_seg; bla[[i]]$intern <- bla[[i]]$intern_seg } }
      SOURCE <- bla
      
      TARGET <- vector("list", length(sT))
      for(i in 1:length(TARGET)){ TARGET[[i]] <- P.TARGET.O[[sT]] }
      bla <- funFun(LIST = TARGET,
                    selNormalisation = x_SelNormalisation,
                    selSampling = x_SelSampling,
                    Radius = targets$Diameter[[sT]]/2,
                    PtsNb = x_PtsNb,
                    PtsDist = x_PtsDist)
      if (SelSampling=="Sampling_equidistant") { for(i in 1:length(bla)) { bla[[i]]$extern <- bla[[i]]$extern_seg; bla[[i]]$intern <- bla[[i]]$intern_seg } }
      TARGET <- bla
      
      P.source <<- SOURCE[[1]]
      P.target <<- TARGET[[1]]
      
      # Source (adjustement with the rim)
      adj <- max(P.source$profil[,2])
      P.source$profil[,2] <- P.source$profil[,2]-adj
      P.source$extern[,2] <- P.source$extern[,2]-adj
      P.source$intern[,2] <- P.source$intern[,2]-adj
      
      # plot the plot
      output$plotTemp1 <- renderPlot({
        icp.p.f(P.source=P.source, P.target = P.target, is.rim=F, par = x_par, projection = F, minfun = "rmsd",
                 method=c("rtsz"), EI="EI", verbose = F, plots=T, plots.full=F)
      })
      output$plotTemp2 <- renderPlot({ })
      output$plotTemp3 <- renderPlot({ })
      output$plotTemp4 <- renderPlot({
        icp.p.all(P.source=P.source, P.target = P.target, is.rim=F, par = x_par, projection = F, minfun = "rmsd",
                 method=c("rtsz"), EI="EI", verbose = F, plots=T, plots.full=T)
      })
    }
  }
  
  plotResult_ICP <- function() {
    if (is.null(P.SOURCE.O[[sS]]$res_ICP)) { notifyNotAnalysed() }
    else {
      x_d <<- P.SOURCE.O[[sS]]$res_ICP[sT,5]
      res_log <<- P.SOURCE.O[[sS]]$res_ICP_log
      
      x_par <<- P.SOURCE.O[[sS]]$res_ICP[sT,1:4]
      
      x_SelNormalisation <<- res_log$SelNormalisation
      x_SelSampling      <<- res_log$SelSampling
      x_PtsNb            <<- as.numeric(res_log$PtsNb)
      x_PtsDist          <<- as.numeric(res_log$PtsDist)
      
      # 0) Preparations
      SOURCE <- vector("list", length(sS))
      for(i in 1:length(SOURCE)){ SOURCE[[i]] <- P.SOURCE.O[[sS]] }
      bla <- funFun(LIST = SOURCE,
                    selNormalisation = x_SelNormalisation,
                    selSampling = x_SelSampling,
                    Radius = sources$Diameter[[sS]]/2,
                    PtsNb = x_PtsNb,
                    PtsDist = x_PtsDist)
      if (SelSampling=="Sampling_equidistant") { for(i in 1:length(bla)) { bla[[i]]$extern <- bla[[i]]$extern_seg; bla[[i]]$intern <- bla[[i]]$intern_seg } }
      SOURCE <- bla
      
      TARGET <- vector("list", length(sT))
      for(i in 1:length(TARGET)){ TARGET[[i]] <- P.TARGET.O[[sT]] }
      bla <- funFun(LIST = TARGET,
                    selNormalisation = x_SelNormalisation,
                    selSampling = x_SelSampling,
                    Radius = targets$Diameter[[sT]]/2,
                    PtsNb = x_PtsNb,
                    PtsDist = x_PtsDist)
      if (SelSampling=="Sampling_equidistant") { for(i in 1:length(bla)) { bla[[i]]$extern <- bla[[i]]$extern_seg; bla[[i]]$intern <- bla[[i]]$intern_seg } }
      TARGET <- bla
      
      P.source <<- SOURCE[[1]]
      P.target <<- TARGET[[1]]
      
      # Source (adjustement with the rim)
      adj <- max(P.source$profil[,2])
      P.source$profil[,2] <- P.source$profil[,2]-adj
      P.source$extern[,2] <- P.source$extern[,2]-adj
      P.source$intern[,2] <- P.source$intern[,2]-adj
      
      # plot the plot
      output$plotTemp1 <- renderPlot({
        icp.p.f(P.source=P.source, P.target = P.target, is.rim=F, par = x_par, projection = F, minfun = "rmsd",
                 method=c("rtsz"), EI="EI", verbose = F, plots=T, plots.full=T)
      })
      output$plotTemp2 <- renderPlot({ })
      output$plotTemp3 <- renderPlot({ })
      output$plotTemp4 <- renderPlot({
        icp.p.all(P.source=P.source, P.target = P.target, is.rim=F, par = x_par, projection = F, minfun = "rmsd",
                 method=c("rtsz"), EI="EI", verbose = F, plots=T, plots.full=T)
      })
    }
  }
  
  
  funRun <- function() {}    ### NOT USED, SERVES FOR CODE SEPARATION
  
  observeEvent(input$btnRun, {
    # test whatever somethin is selected
    if(is.null(P.TARGET.O)) { showNotification(paste("No targets loaded. Please load targets"), duration = dur); return() }
    if(is.null(sS)) { showNotification(paste("No fragment selected. Please select source fragment"), duration = dur); return() }
    showNotification(paste("Matching started"), duration = dur)
    
    sS <<- input$dataTableSources_rows_selected
    
    # Parameters overall
    SelNormalisation <<- input$selNormalisation
    SelSampling <<- input$selSampling
    PtsNb <<- input$sldNbOfPoints
    PtsDist <<- input$sldPtsDistance
    
    SelMethod <<- input$selMethod
    
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (SelMethod=='selBEZ' || SelMethod=='selDCT' || SelMethod=='selRDP' || SelMethod=='selRTC' || SelMethod == 'selselALL') {
      SelSampling <<- "Sampling_equidistant"
    }
    
    # Parameters selBEZ

    # Parameters selDCT
    SldDCT_nbHarm <<- input$sldDCT_nbHarm
    
    # Parameters selRDP
    SldRDP_nbPts <<- input$sldRDP_nbPts
    
    # Parameters selRTC
    w <<- 0       # weigth of points along outline (TO DO)
    SldRTC_W <<- input$sldRTC_W
    
    # Parameters selICP_Rim
    
    # Parameters selICP_Rigid
    SldICP_rigid_z <<- input$sldICP_rigid_z
    
    # Parameters selICP
    SldICP_r <<- input$sldICP_r
    SldICP_theta <<- deg2rad(input$sldICP_theta)
    SldICP_size <<- input$sldICP_size
    
    SldICP_MaxIter <<- input$sldICP_MaxIter
    SldICP_Temp <<- input$sldICP_Temp
    SldICP_Tmax <<- input$sldICP_Tmax
    
    # This part is the same for all methods
    SOURCE_RUN <<- funFun(LIST = P.SOURCE.O,
                          selNormalisation = SelNormalisation,
                          selSampling = SelSampling,
                          Radius = sources$Diameter/2,
                          PtsNb = PtsNb,
                          PtsDist = PtsDist)
    
    # targets do not need to be recalculated if nothing has changed (sources may change when new break pottery is made)
    if (SelNormalisation==SelNormalisation_RUNNED && SelSampling==SelSampling_RUNNED  && PtsNb==PtsNb_RUNNED && PtsDist==PtsDist_RUNNED) {
      showNotification(paste("No need to recalculate targets"), duration = dur)
    }
    else {
      TARGET_RUN <<- funFun(LIST = P.TARGET.O,
                            selNormalisation = SelNormalisation,
                            selSampling = SelSampling,
                            Radius = targets$Diameter/2,
                            PtsNb = PtsNb,
                            PtsDist = PtsDist)
    }
    
    # Run calculations
    if (SelMethod=='selBEZ' || SelMethod=='selALL') { run_BEZ() }
    if (SelMethod=='selDCT' || SelMethod=='selALL') { run_DCT() }
    if (SelMethod=='selRDP' || SelMethod=='selALL') { run_RDP() }
    if (SelMethod=='selRTC' || SelMethod=='selALL') { run_RTC() }
    if (SelMethod=='selICP_Rim' || SelMethod=='selALL') { run_ICP_Rim() }
    if (SelMethod=='selICP_Rigid') { run_ICP_Rigid() }
    if (SelMethod=='selICP') { run_ICP() }

    # Update counters
    SelNormalisation_RUNNED <<- SelNormalisation
    SelSampling_RUNNED <<- SelSampling
    PtsNb_RUNNED <<- PtsNb
    PtsDist_RUNNED <<- PtsDist
    
    # Get information about what have been analysed
    BEZ <- rep("",length(P.SOURCE.O))
    DCT <- rep("",length(P.SOURCE.O))
    RDP <- rep("",length(P.SOURCE.O))
    RTC <- rep("",length(P.SOURCE.O))
    ICP_rim <- rep("",length(P.SOURCE.O))
    ICP_rigid <- rep("",length(P.SOURCE.O))
    ICP <- rep("",length(P.SOURCE.O))
    
    for (i in 1:length(P.SOURCE.O)) {
      if (!is.null(P.SOURCE.O[[i]]$res_BEZ)) { BEZ[i] <- "O"  }
      if (!is.null(P.SOURCE.O[[i]]$res_DCT)) { DCT[i] <- "O"  }
      if (!is.null(P.SOURCE.O[[i]]$res_RDP)) { RDP[i] <- "O"  }
      if (!is.null(P.SOURCE.O[[i]]$res_RTC)) { RTC[i] <- "O"  }
      if (!is.null(P.SOURCE.O[[i]]$res_ICP_rim)) { ICP_rim[i] <- "O"  }
      if (!is.null(P.SOURCE.O[[i]]$res_ICP_rigid)) { ICP_rigid[i] <- "O"  }
      if (!is.null(P.SOURCE.O[[i]]$res_ICP)) { ICP[i] <- "O"  }
    }
    sources <<- sources[,1:2]
    sources <<- cbind(sources, BEZ, DCT, RDP, RTC, ICP_rim, ICP_rigid, ICP)
    output$dataTableSources = DT::renderDataTable(sources, server = FALSE)
    showNotification(paste("COMPARISON FINISHED"), duration = dur)
  })
  
  run_BEZ <- function() {
    # Preparations
    SOURCE <- SOURCE_RUN
    TARGET <- TARGET_RUN
    
    SOURCE <- funFunDct(LIST = SOURCE, PtsDist = PtsDist)
    TARGET <- funFunDct(LIST = TARGET, PtsDist = PtsDist)
    
    # Calculations
    i=1
    for (i in 1:length(sS)) {
      P.source <- SOURCE[[sS[i]]]
      res <- rep(1000, length(TARGET))
      proc <- 1/length(TARGET)
      withProgress(message = paste('Matching profile - ', i, '/', length(sS)), value = 0, {
        for (j in 1:length(TARGET)) {
          P.target <- TARGET[[j]]
          if ( min(P.source$s)>=min(P.target$s) && max(P.source$s)<=max(P.target$s)) {            # check if the source curve fits to target curve
            res[j] <- compare.bez(P.source, P.target, method="rmsd", plots=F)
          }
          else{  }
          setProgress(proc*j)
        }
      })
      P.SOURCE.O[[sS[i]]]$res_BEZ <<- res
      P.SOURCE.O[[sS[i]]]$res_BEZ_log <<- list(SelNormalisation=SelNormalisation, SelSampling=SelSampling, PtsNb=PtsNb, PtsDist=PtsDist)
    }
    showNotification(paste("Matching finished"), duration = dur)
  }
  
  run_DCT <- function() {
    # Preparations
    SOURCE <- SOURCE_RUN
    TARGET <- TARGET_RUN
    
    SOURCE <- funFunDct(LIST = SOURCE, PtsDist = PtsDist)
    TARGET <- funFunDct(LIST = TARGET, PtsDist = PtsDist)
    
    # Calculations
    for (i in 1:length(sS)) {
      P.source <- SOURCE[[sS[i]]]
      res <- rep(1000, length(TARGET))
      proc <- 1/length(TARGET)
      withProgress(message = paste('Matching profile - ', i, '/', length(sS)), value = 0, {
        for (j in 1:length(TARGET)) {
          P.target <- TARGET[[j]]
          if ( min(P.source$s)>=min(P.target$s) && max(P.source$s)<=max(P.target$s)) {            # check if the source curve fits to target curve
            res[j] <- compare.dct(P.source, P.target, nharm=SldDCT_nbHarm, method="rmsd", plots=F)
          }
          else{  }
          setProgress(proc*j)
        }
      })
      P.SOURCE.O[[sS[i]]]$res_DCT <<- res
      P.SOURCE.O[[sS[i]]]$res_DCT_log <<- list(SelNormalisation=SelNormalisation, SelSampling=SelSampling, PtsNb=PtsNb, PtsDist=PtsDist, SldDCT_nbHarm=SldDCT_nbHarm)
    }
    showNotification(paste("Matching finished"), duration = dur)
  }

  run_RDP <- function() {
    # Preparations
    SOURCE <- SOURCE_RUN
    TARGET <- TARGET_RUN
    
    SOURCE <- funFunDct(LIST = SOURCE, PtsDist = PtsDist)
    TARGET <- funFunDct(LIST = TARGET, PtsDist = PtsDist)
    
    # Calculations
    for (i in 1:length(sS)) {
      P.source <- SOURCE[[sS[i]]]
      res <- rep(1000, length(TARGET))
      proc <- 1/length(TARGET)
      withProgress(message = paste('Matching profile - ', i, '/', length(sS)), value = 0, {
        for (j in 1:length(TARGET)) {
          P.target <- TARGET[[j]]
          res[j] <- compare.rdp(P.source, P.target, method="rmsd", nbPoints=SldRDP_nbPts, plots=F)
          setProgress(proc*j)
        }
      })
      P.SOURCE.O[[sS[i]]]$res_RDP <<- res
      P.SOURCE.O[[sS[i]]]$res_RDP_log <<- list(SelNormalisation=SelNormalisation, SelSampling=SelSampling, PtsNb=PtsNb, PtsDist=PtsDist, SldRDP_nbPts=SldRDP_nbPts)
    }
    showNotification(paste("Matching finished"), duration = dur)
  }
  
  run_RTC <- function() {
    # Preparations
    SOURCE <- SOURCE_RUN
    TARGET <- TARGET_RUN
    
    SOURCE <- funFunRtc(LIST = SOURCE, PtsDist = PtsDist)
    TARGET <- funFunRtc(LIST = TARGET, PtsDist = PtsDist)
    
    # 1) Define parameters
    if (SldRTC_W[1]==0.33 && (SldRTC_W[2]==0.66)) {
      W <<- c(1,1,1)/3
    } else {
      WR <- SldRTC_W[1]
      WT <- SldRTC_W[2]-SldRTC_W[1]
      WK <- 1-WR-WT
      W <<- c(WR, WT, WK)
    }
    
    # Calculations
    for (i in 1:length(sS)) {
      P.source <- SOURCE[[sS[i]]]
      res <- rep(1000, length(TARGET))
      proc <- 1/length(TARGET)
      withProgress(message = paste('Matching profile - ', i, '/', length(sS)), value = 0, {
        for (j in 1:length(TARGET)) {
          P.target <- TARGET[[j]]
          if ( min(P.source$s)>=min(P.target$s) && max(P.source$s)<=max(P.target$s)) {            # check if the source curve fits to target curve
            res[j] <- compare.rtc.full(P.source, P.target, W=W, plots=F) 
          }
          else{  }
          setProgress(proc*j)
        }
      })
      P.SOURCE.O[[sS[i]]]$res_RTC <<- res
      P.SOURCE.O[[sS[i]]]$res_RTC_log <<- list(SelNormalisation=SelNormalisation, SelSampling=SelSampling, PtsNb=PtsNb, PtsDist=PtsDist, w=w, SldRTC_W=SldRTC_W)
    }
    showNotification(paste("Matching finished"), duration = dur)
  }

  run_ICP_Rim <- function() {
    # Preparations
    SOURCE <- SOURCE_RUN
    TARGET <- TARGET_RUN
    
    if (SelSampling=="Sampling_equidistant") { for(i in 1:length(SOURCE)) { SOURCE[[i]]$extern <- SOURCE[[i]]$extern_seg; SOURCE[[i]]$intern <- SOURCE[[i]]$intern_seg } }
    if (SelSampling=="Sampling_equidistant") { for(i in 1:length(TARGET)) { TARGET[[i]]$extern <- TARGET[[i]]$extern_seg; TARGET[[i]]$intern <- TARGET[[i]]$intern_seg } }
    SOURCE <- SOURCE
    TARGET <- TARGET
    
    # Calculations
    for (i in 1:length(sS)) {
      P.source <- SOURCE[[sS[i]]]
      
      # Source (adjustement with the rim)
      adj <- max(P.source$profil[,2])
      P.source$profil[,2] <- P.source$profil[,2]-adj
      P.source$extern[,2] <- P.source$extern[,2]-adj
      P.source$intern[,2] <- P.source$intern[,2]-adj
      
      res <- rep(1000, length(TARGET))
      proc <- 1/length(TARGET)
      withProgress(message = paste('Matching profile - ', i, '/', length(sS)), value = 0, {
        j=1
        for (j in 1:length(TARGET)) {
          P.target <- TARGET[[j]]
          par <- c(0,0,1,0)
          res[j] <- icp.p.f(P.source=P.source, P.target = P.target, is.rim=F, par = par, projection = F, minfun = "rmsd",
                             method=c("rtsz"), EI="EI", verbose = T, plots=F)
          setProgress(proc*j)
        }
      })
      P.SOURCE.O[[sS[i]]]$res_ICP_rim <<- res
      P.SOURCE.O[[sS[i]]]$res_ICP_rim_log <<- list(SelNormalisation=SelNormalisation, SelSampling=SelSampling, PtsNb=PtsNb, PtsDist=PtsDist)
      showNotification(paste("Matching finished"), duration = dur)
    }
  }
  
  run_ICP_Rigid <- function() {
    # Preparations
    SOURCE <- SOURCE_RUN
    TARGET <- TARGET_RUN
    
    if (SelSampling=="Sampling_equidistant") { for(i in 1:length(SOURCE)) { SOURCE[[i]]$extern <- SOURCE[[i]]$extern_seg; SOURCE[[i]]$intern <- SOURCE[[i]]$intern_seg } }
    if (SelSampling=="Sampling_equidistant") { for(i in 1:length(TARGET)) { TARGET[[i]]$extern <- TARGET[[i]]$extern_seg; TARGET[[i]]$intern <- TARGET[[i]]$intern_seg } }
    SOURCE <- SOURCE
    TARGET <- TARGET
    
    for (i in 1:length(sS)) {
      P.source <- SOURCE[[sS[i]]]

      # Source (adjustement with the rim)
      adj <- max(P.source$profil[,2])
      P.source$profil[,2] <- P.source$profil[,2]-adj
      P.source$extern[,2] <- P.source$extern[,2]-adj
      P.source$intern[,2] <- P.source$intern[,2]-adj
      
      RES <- matrix(numeric(0),length(TARGET),5)
      colnames(RES) <- c("r", "theta", "size", "z", "result")
      proc <- 1/length(TARGET)
      withProgress(message = paste('Matching profile - ', i, '/', length(sS)), value = 0, {
        for (j in 1:length(TARGET)) {
          P.target <- TARGET[[j]]
          
          # Definition of the Z searching parameter (method=c("rtsz"))
          Z <- c(min(P.target$profil[,2])-min(P.source$profil[,2]),0)

          ZZ <- c(seq(0, abs(Z[1]), by=SldICP_rigid_z), abs(Z[1]))
          ZZ <- -ZZ
          RES_TEMP <- matrix(numeric(0),length(ZZ),5)
          colnames(RES_TEMP) <- c("r", "theta", "size", "z", "result")
          for(zi in 1:length(ZZ)) {
            par = c(0,0,1,ZZ[zi])
            RES_TEMP[zi,1:4] <- par
            RES_TEMP[zi,5] <- icp.p.f(P.source=P.source, P.target = P.target, is.rim=F, par = par, projection = F, minfun = "rmsd",
                                       method=c("rtsz"), EI="EI", verbose = F, plots=F)
          }
          RES[j,] <- RES_TEMP[which.min(RES_TEMP[,5]),]
          setProgress(proc*j)
        }
      })
      P.SOURCE.O[[sS[i]]]$res_ICP_rigid <<- RES
      P.SOURCE.O[[sS[i]]]$res_ICP_rigid_log <<- list(SelNormalisation=SelNormalisation, SelSampling=SelSampling, PtsNb=PtsNb, PtsDist=PtsDist)
      showNotification(paste("Matching finished"), duration = dur)
    }
  }
  
  run_ICP <- function() {
    # Preparations
    SOURCE <- SOURCE_RUN
    TARGET <- TARGET_RUN
    
    if (SelSampling=="Sampling_equidistant") { for(i in 1:length(SOURCE)) { SOURCE[[i]]$extern <- SOURCE[[i]]$extern_seg; SOURCE[[i]]$intern <- SOURCE[[i]]$intern_seg } }
    if (SelSampling=="Sampling_equidistant") { for(i in 1:length(TARGET)) { TARGET[[i]]$extern <- TARGET[[i]]$extern_seg; TARGET[[i]]$intern <- TARGET[[i]]$intern_seg } }
    SOURCE <- SOURCE
    TARGET <- TARGET
    
    # Calculations
    for (i in 1:length(sS)) {
      P.source <- SOURCE[[sS[i]]]
      
      # Source (adjustement with the rim)
      adj <- max(P.source$profil[,2])
      P.source$profil[,2] <- P.source$profil[,2]-adj
      P.source$extern[,2] <- P.source$extern[,2]-adj
      P.source$intern[,2] <- P.source$intern[,2]-adj
      
      RES <- matrix(numeric(0),length(TARGET),5)
      colnames(RES) <- c("r", "theta", "size", "z", "result")
      
      proc <- 1/length(TARGET)
      withProgress(message = paste('Matching profile - ', i, '/', length(sS)), value = 0, {
        for (j in 1:length(TARGET)) {
          P.target <- TARGET[[j]]
          
          # Definition of the Z searching parameter (method=c("rtsz"))
          Z <- c(min(P.target$profil[,2])-min(P.source$profil[,2]),0)
          
          # Definition of limits (R,Theta,Size,Z)
          low.lim <- c(SldICP_r[1], SldICP_theta[1], SldICP_size[1], Z[1])
          up.lim  <- c(SldICP_r[2], SldICP_theta[2], SldICP_size[2], Z[2])
          
          # Optimisation - SANN
          C <- list(maxit=SldICP_MaxIter, temp=SldICP_Temp, tmax=SldICP_Tmax, trace=TRUE, REPORT=10000)
          D <- 4
          s <- c(0,0,1,0)
          par <- c(0,0,1,0)
          
          hchange <- function (par, lower, upper, dist, round=TRUE,...) {
            step <- dist(D,...) # slight step
            if (round) step=round(step)
            par1 <- par+step
            return(ifelse(par1<lower,lower,ifelse(par1>upper,upper,par1)))
          }

          bchange <- function(par) { D=length(par); hchange(par, lower=low.lim, upper=up.lim, rnorm, mean=0, sd=SD, round = F) }
          minfun <- function(par) { icp.p.f(P.source, P.target, is.rim=F, par=par, projection=F, EI="EI", minfun="rmsd", method=c("rtsz"), plots=F, plots.full=F, verbose=F) }
          solution <- optim(s, minfun, gr=bchange, method="SANN", control=C)
          RES[j,] <- c(solution$par,solution$value)
          setProgress(proc*j)
        }
      })
      P.SOURCE.O[[sS[i]]]$res_ICP <<- RES
      P.SOURCE.O[[sS[i]]]$res_ICP_log <<- list(SelNormalisation=SelNormalisation, SelSampling=SelSampling, PtsNb=PtsNb, PtsDist=PtsDist)
      showNotification(paste("Matching finished"), duration = dur)
    }
  }
  
  
  funSelShow <- function(){}    ### NOT USED, SERVES FOR CODE SEPARATION
  
  observeEvent(c(input$selShow), {
    SelShow <<- input$selShow
    
    if (SelShow=='selShowAnalyse') {
      output$dataTableSources = DT::renderDataTable(sources, server = FALSE)
      output$dataTableTargets = DT::renderDataTable(targets, server = FALSE)
      
      output$plotTemp2 <- renderPlot({ plot.blank() })
      output$plotTemp3 <- renderPlot({ plot.blank() })
      output$plotTemp4 <- renderPlot({ plot.blank() })
    }
    else if (SelShow=='selShowResults') {
      showResults()
      
      output$plotTemp2 <- renderPlot({ plot.blank() })
      output$plotTemp3 <- renderPlot({ plot.blank() })
      output$plotTemp4 <- renderPlot({ plot.blank() })
    }
  })
  
  showAnalyse <- function() {
    # Parameters
    SelNormalisation <<- input$selNormalisation
    SelSampling <<- input$selSampling
    PtsNb <<- input$sldNbOfPoints
    PtsDist <<- input$sldPtsDistance
    
    # Parameters selDCT
    CbxAdvancedDCT <<- input$cbxAdvancedDCT
    SldDCT_nbHarm <<- input$sldDCT_nbHarm

    # Parameters selRDP
    CbxAdvancedRDP <<- input$cbxAdvancedRDP
    SldRDP_nbPts <<- input$sldRDP_nbPts

    # Parameters selICP_Rigid
    CbxAdvancedICP_Rigid <<- input$cbxAdvancedICP_Rigid
    SldICP_rigid_z <<- input$sldICP_rigid_z

    # Parameters selICP
    CbxAdvancedICP <<- input$cbxAdvancedICP
    SldICP_r <<- input$sldICP_r
    SldICP_theta <<- deg2rad(input$sldICP_theta)
    SldICP_size <<- input$sldICP_size

    if (is.null(sT) && is.null(sS)) {
      output$plotTemp1<- renderPlot({ plot.blank() })
    }
    else {
      # Prepare sources
      if (!is.null('P.SOURCE.O') && !is.null(sS) ) {
        SOUR <<- vector("list", length(sS))
        for(i in 1:length(SOUR)){ SOUR[[i]] <<- P.SOURCE.O[[sS[i]]] }
        SOUR <<- funFun(LIST=SOUR,
                        selNormalisation=SelNormalisation,
                        selSampling=SelSampling,
                        Radius=sources$Diameter[sS]/2,
                        PtsNb=PtsNb,
                        PtsDist=PtsDist)
      }
      
      # Prepare targets
      if (!is.null('P.TARGET.O') && !is.null(sT) ) {
        TARG <<- vector("list", length(sT))
        for(i in 1:length(TARG)){ TARG[[i]] <<- P.TARGET.O[[sT[i]]] }
        TARG <<- funFun(LIST=TARG,
                        selNormalisation=SelNormalisation,
                        selSampling=SelSampling,
                        Radius=targets$Diameter[sT]/2,
                        PtsNb=PtsNb,
                        PtsDist=PtsDist)
      }
      
      # Calculate limits of the plot
      xl = 0; yl = 0
      if(!is.null(sT)) {
        for(i in 1:length(TARG)) {
          xln <- min(TARG[[i]]$profil[,1])
          yln <- min(TARG[[i]]$profil[,2])
          if(xln<xl) { xl=xln }
          if(yln<yl) { yl=yln }
        }
      }
      if(!is.null(sS)) {
        for(i in 1:length(SOUR)) {
          xln <- min(SOUR[[i]]$profil[,1])
          yln <- min(SOUR[[i]]$profil[,2])
          if(xln<xl) { xl=xln }
          if(yln<yl) { yl=yln }
        }
      }
      
      # plot the plot
      output$plotTemp1<- renderPlot({
        plot(0, asp=1, type="n", xlab="r", ylab="z", xlim=c(xl,0), ylim=c(yl,0))
        
        if(!is.null(sT) ) {
          for (i in 1:length(TARG)) {
            polygon(TARG[[i]]$profil, col="lightgrey", border = "darkgrey", lwd=1)
            if (SelSampling=='Sampling_segments') {
              points(TARG[[i]]$extern, col="blue", cex=Cex)
              points(TARG[[i]]$intern, col="red", cex=Cex)
            }
            else if (SelSampling=='Sampling_equidistant') {
              points(TARG[[i]]$extern_seg, col="blue", cex=Cex)
              points(TARG[[i]]$intern_seg, col="red", cex=Cex)
            }
          }
        }
        
        if(!is.null(sS) ) {
          for (i in 1:length(SOUR)) {
            polygon(SOUR[[i]]$profil, col="darkgrey", border = "black", lwd=1)
            if (SelSampling=='Sampling_segments') {
              points(SOUR[[i]]$extern, col="blue", cex=Cex, pch=19)
              points(SOUR[[i]]$intern, col="red", cex=Cex, pch=19)
              
              if (CbxAdvancedDCT==T) {
                if(SldDCT_nbHarm>dim(SOUR[[i]]$profil)[1]/2) {
                  SldRDP_nbPts <- floor(dim(SOUR[[i]]$profil)[1]/2)
                  showNotification(paste("Maximal number of DCT harmonics for a given fragment is", SldRDP_nbPts), duration = dur)
                }
                points(idct(Sk=dct(SOUR[[i]]$profil), harm=SldDCT_nbHarm, N=dim(SOUR[[i]]$profil)[1]), col="green", cex=Cex*1.5, pch=19)
                lines(idct(Sk=dct(SOUR[[i]]$profil), harm=SldDCT_nbHarm, N=dim(SOUR[[i]]$profil)[1]), col="green", lwd=2)
              }
              
              if (CbxAdvancedRDP==T) {
                if (SldRDP_nbPts>dim(SOUR[[i]]$profil)[1]) {
                  SldRDP_nbPts <- dim(SOUR[[i]]$profil)[1]
                  showNotification(paste("Maximal number of points of RDP for a given fragment is", SldRDP_nbPts), duration = dur)
                }
                points(DouglasPeuckerNbPoints(SOUR[[i]]$profil[,1], SOUR[[i]]$profil[,2], SldRDP_nbPts), col="green", cex=Cex*1.5, pch=19)
                lines(DouglasPeuckerNbPoints(SOUR[[i]]$profil[,1], SOUR[[i]]$profil[,2], SldRDP_nbPts), col="green", lwd=2)
              }

            }
            else if (SelSampling=='Sampling_equidistant') {
              points(SOUR[[i]]$extern_seg, col="blue", cex=Cex, pch=19)
              points(SOUR[[i]]$intern_seg, col="red", cex=Cex, pch=19)
              
              if (CbxAdvancedDCT==T) {
                if(SldDCT_nbHarm>dim(SOUR[[i]]$extint_seg)[1]/2) {
                  SldRDP_nbPts <- floor(dim(SOUR[[i]]$extint_seg)[1]/2)
                  showNotification(paste("Maximal number of DCT harmonics for a given fragment is", SldRDP_nbPts), duration = dur)
                }
                points(idct(Sk=dct(SOUR[[i]]$extint_seg), harm=SldDCT_nbHarm, N=dim(SOUR[[i]]$extint_seg)[1]), col="green", cex=Cex*1.5, pch=19)
                lines(idct(Sk=dct(SOUR[[i]]$extint_seg), harm=SldDCT_nbHarm, N=dim(SOUR[[i]]$extint_seg)[1]), col="green", lwd=2)
              }
              
              if (CbxAdvancedRDP==T) {
                if (SldRDP_nbPts>dim(SOUR[[i]]$extint_seg)[1]) {
                  SldRDP_nbPts <- dim(SOUR[[i]]$extint_seg)[1]
                  showNotification(paste("Maximal number of points of RDP for a given fragment is", SldRDP_nbPts), duration = dur)
                }
                points(DouglasPeuckerNbPoints(SOUR[[i]]$extint_seg[,1], SOUR[[i]]$extint_seg[,2], SldRDP_nbPts), col="green", cex=Cex*1.5, pch=19)
                lines(DouglasPeuckerNbPoints(SOUR[[i]]$extint_seg[,1], SOUR[[i]]$extint_seg[,2], SldRDP_nbPts), col="green", lwd=2)
              }
              
            }
          }
        }
        
        if (CbxAdvancedICP_Rigid==T) {
          z <- seq(from=0, to=abs(yl), by = SldICP_rigid_z)
          z <- -z
          z <- cbind(0,z)
          points(z, cex=Cex, pch=3)
        }
        
        if (CbxAdvancedICP==T) {
          if (length(sS)==1 && length(sT)==1) {
            P.source <<- SOUR[[1]]
            
            # Source (adjustement with the rim)
            adj <- max(P.source$profil[,2])
            P.source$profil[,2] <- P.source$profil[,2]-adj
            P.source$extern[,2] <- P.source$extern[,2]-adj
            P.source$intern[,2] <- P.source$intern[,2]-adj
            
            P.target <<- TARG[[1]]
            
            # To make it simple
            R <- SldICP_r
            Theta <- SldICP_theta
            Size <- SldICP_size
            
            # Definition of the Z searching parameter (method=c("rtsz"))
            Z <- c(min(P.target$profil[,2])-min(P.source$profil[,2]),0)
            
            # Definition of limits (R,Theta,Size,Z)
            low.lim <- c(R[1], Theta[1], Size[1], Z[1])
            up.lim  <- c(R[2], Theta[2], Size[2], Z[2])

            # pPlot equilibre
            profil.plot(profil.transform(P.source, par=c(0,0,1,Z[1])), add=T)
            profil.plot(profil.transform(P.source, par=c(0,0,1,Z[2])), add=T)
            
            profil.plot(profil.transform(P.source, par=c(R[1], Theta[1], Size[1], Z[1])), add=T)
            profil.plot(profil.transform(P.source, par=c(R[1], Theta[1], Size[1], Z[2])), add=T)
            profil.plot(profil.transform(P.source, par=c(R[1], Theta[1], Size[2], Z[2])), add=T)
            profil.plot(profil.transform(P.source, par=c(R[1], Theta[2], Size[2], Z[2])), add=T)
            profil.plot(profil.transform(P.source, par=c(R[2], Theta[2], Size[2], Z[2])), add=T)
            
            profil.plot(profil.transform(P.source, par=c(R[1], Theta[2], Size[1], Z[1])), add=T)
            profil.plot(profil.transform(P.source, par=c(R[1], Theta[2], Size[2], Z[1])), add=T)
            
            profil.plot(profil.transform(P.source, par=c(R[1], Theta[2], Size[1], Z[2])), add=T)
            profil.plot(profil.transform(P.source, par=c(R[2], Theta[1], Size[2], Z[1])), add=T)
            
            profil.plot(profil.transform(P.source, par=c(R[2], Theta[2], Size[2], Z[1])), add=T)
            profil.plot(profil.transform(P.source, par=c(R[2], Theta[1], Size[1], Z[1])), add=T)
            
            profil.plot(profil.transform(P.source, par=c(R[2], Theta[1], Size[1], Z[2])), add=T)
          }
        }
      })
    }
    return(NULL)
  }
  
  showResults <- function() {
    sS <<- input$dataTableSources_rows_selected
    if (is.null(sS)) { }
    else if (length(sS)>1) {
      # updateSelectInput(session, 'selShow', label=NULL, c("Show Results" = "", "Normal" = "selShowAnalyse", "Results" = "selShowResults"), selected="selShowAnalyse")
      showNotification(paste("Select just one source"), duration = dur)
    }
    
    else {
      if (input$selMethod=="selBEZ") {
        if (is.null(P.SOURCE.O[[sS]]$res_BEZ)) { notifyNotAnalysed() }
        else { res <- P.SOURCE.O[[sS]]$res_BEZ }
      }
      
      else if (input$selMethod=="selDCT") {
        if (is.null(P.SOURCE.O[[sS]]$res_DCT)) { notifyNotAnalysed() }
        else { res <- P.SOURCE.O[[sS]]$res_DCT }
      }

      else if (input$selMethod=="selRDP") {
        if (is.null(P.SOURCE.O[[sS]]$res_RDP)) { notifyNotAnalysed() }
        else { res <- P.SOURCE.O[[sS]]$res_RDP }
      }

      else if (input$selMethod=="selRTC") {
        if (is.null(P.SOURCE.O[[sS]]$res_RTC)) { notifyNotAnalysed() }
        else { res <- P.SOURCE.O[[sS]]$res_RTC }
      }
      
      else if (input$selMethod=="selICP_Rim") {
        if (is.null(P.SOURCE.O[[sS]]$res_ICP_rim)) { notifyNotAnalysed() }
        else { res <- P.SOURCE.O[[sS]]$res_ICP_rim }
      }
      
      else if (input$selMethod=="selICP_Rigid") {
        if (is.null(P.SOURCE.O[[sS]]$res_ICP_rigid)) { notifyNotAnalysed() }
        else { res <- P.SOURCE.O[[sS]]$res_ICP_rigid[,5] }
      }
      
      else if (input$selMethod=="selICP") {
        if (is.null(P.SOURCE.O[[sS]]$res_ICP)) { notifyNotAnalysed() }
        else { res <- P.SOURCE.O[[sS]]$res_ICP[,5] }
      }
      
      if(exists("res")) {
        Value = res
        results <- cbind(targets, Value)
        output$dataTableTargets = DT::renderDataTable(results, server = FALSE)
      }
      else {
        output$dataTableTargets = DT::renderDataTable(targets, server = FALSE)
      }
    }
  }
  
  observeEvent(input$btnBreakPottery, {
    if(!is.null(sS)) {
      if (length(sS)>1) {
        showNotification(paste("Select just one source"), duration = dur)
        return()
      }
      sS <<- sS
      SldBreaks <<- input$sldBreaks
      SldMinPerc <<-input$sldMinPerc
      NumPtsSetSeed <<-input$numPtsSetSeed
      CbxOnlyRim <<- input$cbxOnlyRim
      MinNbOfPts <<- PtsRechant4Breaking/100*SldMinPerc
      
      P.SOURCE.O1 <<- P.SOURCE.O
      n = length(P.SOURCE.O)
      
      profil <<- P.SOURCE.O[[sS]]
      profil <<- profil.rechant(profil, pts=PtsRechant4Breaking)
      
      NEW <<- profil.break(profil = profil, breaks = SldBreaks, only_rim = CbxOnlyRim, min.pts = MinNbOfPts, seed=NumPtsSetSeed, plots=F)
      N <- length(NEW)
      
      P.SOURCE.O <<- vector('list', n+N)
      for(i in 1:n) { P.SOURCE.O[[i]] <<- P.SOURCE.O1[[i]] }
      for(i in 1:N) { P.SOURCE.O[[n+i]] <<- NEW[[i]] }
      
      # Update Diameters
      for(i in 1:N) {
        if (i==1) {
          temp <<- sources[sS,]
          sources <<- rbind(sources,temp)
          
          # Prepare ratios for other sherds
          r <- sources[sS,]$Diameter/2
          prof <- P.SOURCE.O[[sS]]
          # Get the first index of extern and intern
          temp = which.max(c  (prof$extern[1,2],prof$intern[1,2]))
          if (temp==1) { px <- abs(prof$extern[1,1]) } else { px <- abs(prof$intern[1,1]) }
          rat <- r/px
          
        } else {
          prof <- NEW[[i]]
          temp = which.max(c  (prof$extern[1,2],prof$intern[1,2]))
          if (temp==1) { px <- abs(prof$extern[1,1]) } else { px <- abs(prof$intern[1,1]) }
          temp <- sources[sS,]
          temp$Diameter <- px*rat   * 2
          sources <<- rbind(sources,temp)
        }
      }
      P.SOURCE.O <<- P.SOURCE.O
      sources <<- sources
    }
    output$dataTableSources = DT::renderDataTable(sources, server = FALSE)
  })
  
  observeEvent(input$selMethod, {
    if (input$selMethod=='selRTC') {
      updateSelectInput(session, 'selSampling', label=NULL,
                        c("Select sampling" = "",
                          "Segments" = "Sampling_segments",
                          "Equidistant" = "Sampling_equidistant"),
                        selected="Sampling_equidistant")
    }
    if (input$selShow=='selShowResults') {
      showResults()
      output$plotTemp2 <- renderPlot({ plot.blank() })
      output$plotTemp3 <- renderPlot({ plot.blank() })
      output$plotTemp4 <- renderPlot({ plot.blank() })
    }
  })
  
  notifyNotAnalysed <- function() {
    showNotification(paste("Selected individual was not analysed by selected method"), duration = dur)
  }
  
})
