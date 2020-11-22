getTICs <- function(files, method='BPC'){
  path <- readfiles(files)
  tics <- lapply(path,function(filename){
    splitname <- strsplit(filename,"\\.")[[1]]
    if(tolower(splitname[length(splitname)]) == "cdf")
    {
      msobj <- openMSfile(filename,backend="netCDF")
    }else{
      msobj <- openMSfile(filename)
    }
    peakInfo <- peaks(msobj)
    headerInfo <- header(msobj)
    whMS1 <- which(headerInfo$msLevel==1)
    peakInfo <- peakInfo[whMS1]
    peakInfo <- lapply(peakInfo, function(spectrum) {
      keep <- spectrum[,2] > 1e-6
      output <- as.data.frame(spectrum[keep,,drop = FALSE])
      colnames(output) <- c('mz','intensity')
      return(output)
    })
    tic <- sapply(peakInfo,function(s){
      if (method=='BPC'){return(max(s$intensity))}else{return(sum(s$intensity))}
    })
    rt <- round(headerInfo$retentionTime[whMS1],3)
    return(cbind(rt,tic))
  })
  return(list(path=path,tics=tics))
}

viewTICs <- function(tics){

  ui <- fluidPage(
    titlePanel("TIC Viewer"),
    sidebarLayout(
      sidebarPanel(
        sliderInput("inds",
                    "index of sample:",
                    min = 1,
                    max = length(tics$path),
                    value = 1,
                    step=1)
      ),
      mainPanel(
        h5(textOutput("caption")),
        plotlyOutput("Plot")
      )
    )
  )

  server <- function(input, output) {
    output$caption <- renderText({
      tics$path[[input$inds]]
    })
    output$Plot <- renderPlotly({
      tic.i <- tics$tics[[input$inds]]
      plot_ly(x=tic.i[,1], y=tic.i[,2], type = 'scatter', mode = 'lines') %>%
        layout(xaxis = list(title = 'Retention Time (s)'),
               yaxis = list (title = 'Intensity'))
    })
  }
  shinyApp(ui = ui, server = server)
}

getMS <- function(filename){

  splitname <- strsplit(filename,"\\.")[[1]]
  if(tolower(splitname[length(splitname)]) == "cdf")
  {
    msobj <- openMSfile(filename,backend="netCDF")
  }else{
    msobj <- openMSfile(filename)
  }
  peakInfo <- peaks(msobj)
  headerInfo <- header(msobj)
  whMS1 <- which(headerInfo$msLevel==1)
  peakInfo <- peakInfo[whMS1]
  peakInfo <- lapply(peakInfo, function(spectrum) {
    keep <- spectrum[,2] > 1e-6
    output <- as.data.frame(spectrum[keep,,drop = FALSE])
    colnames(output) <- c('mz','intensity')
    return(output)
  })

  rt <- round(headerInfo$retentionTime[whMS1],3)

  return(list(rt=rt,MS=peakInfo))
}

stem <- function(x,y,pch=5,linecol=1,col='blue',cex.lab=1.2,cex.axis=1.3,font=2,...){
  if (missing(y)){
    y = x
    x = 1:length(x)}
  plot(x,y,pch=pch,type='n',col=col,cex.lab=cex.lab,cex.axis=cex.axis,font=font,xlab="m/z",ylab="Intensity")
  for (i in 1:length(x)){
    lines(c(x[i],x[i]), c(0,y[i]),pch=pch,col=col,cex.lab=cex.lab,cex.axis=cex.axis,font=font)
  }
  lines(c(x[1]-2,x[length(x)]+2), c(0,0),pch=pch,col=col,cex.lab=cex.lab,cex.axis=cex.axis,font=font)
}

viewMS <- function(MS){

  ui <- fluidPage(
    titlePanel("MS Viewer"),
    sidebarLayout(
      sidebarPanel(
        sliderInput("inds",
                    "index of scan:",
                    min = 1,
                    max = length(MS$rt),
                    value = 1,
                    step=1)
      ),
      mainPanel(
        h5(textOutput("caption")),
        plotOutput("Plot")
      )
    )
  )

  server <- function(input, output) {
    output$caption <- renderText({
      paste('MS of retention time: ',MS$rt[[input$inds]], ' s')
    })
    output$Plot <- renderPlot({
      MS.i <- MS$MS[[input$inds]]
      stem(MS.i[,1], MS.i[,2])
    })
  }
  shinyApp(ui = ui, server = server)
}

viewPICs <- function(pics){

  ui <- fluidPage(
    titlePanel("PIC Viewer"),
    sidebarLayout(
      sidebarPanel(
        sliderInput("inds",
                    "index of pic:",
                    min = 1,
                    max = length(pics$pics),
                    value = 1,
                    step=1)
      ),
      mainPanel(
        plotlyOutput("Plot"),
        plotlyOutput("Plot1"),
        h4 (textOutput("evals"))
      )
    )
  )

  server <- function(input, output) {
    output$Plot <- renderPlotly({
      pic.i <- pics$pics[[input$inds]]
      peak.i <- pics$peaks[[input$inds]]
      plot_ly(x=pics$scantime[pic.i[,1]], y=pic.i[,2], type = 'scatter', mode = 'lines', showlegend= FALSE) %>%
        layout(xaxis = list(tick0 = 0, title = 'Retention Time (s)'),
               yaxis = list(title = 'Intensity')) %>%
        add_markers(x=pics$scantime[pic.i[peak.i$peakIndex,1]], y=pic.i[peak.i$peakIndex,2])
    })
    output$Plot1 <- renderPlotly({
      pic.i <- pics$pics[[input$inds]]
      pic.i <- pic.i[!is.na(pic.i[,3]),]
      plot_ly(x=pics$scantime[pic.i[,1]], y=pic.i[,3], color=pic.i[,2] ,type = 'scatter') %>%
        layout(xaxis = list(tick0 = 0, title = 'Retention Time (s)'),
               yaxis = list(title = 'M/Z'))
    })
    output$evals <- renderText({
      pic.i <- pics$pics[[input$inds]]
      gf <- gaussfit(pic.i)
      sp <- getSharpness(pic.i)
      paste('Evaluation result: ', 'gaussian fitness ', gf, ' sharpness ', sp)
    })
  }
  shinyApp(ui = ui, server = server)
}

viewGroups <- function(groups){

  ui <- fluidPage(
    titlePanel("groups Viewer"),
    sidebarLayout(
      sidebarPanel(
        sliderInput("inds",
                    "index of group:",
                    min = 1,
                    max = max(groups$peakmat[,'group']),
                    value = 1,
                    step=1)
      ),
      mainPanel(
        h4 (textOutput("info"),
        plotlyOutput("Plot")
        )
      )
    )
  )

  server <- function(input, output) {
    output$Plot <- renderPlotly({
      p <- plot_ly()
      candidates <- groups$peakmat[groups$peakmat[,'group']==input$inds,]
      for (i in 1:nrow(candidates)){
        pic <- groups$picset[[candidates[i,'sample']]]$pics[[candidates[i,'index']]]
        p <- add_trace(p, x=pic[,1], y=pic[,2], type = 'scatter', mode = 'lines')
      }
      p
    })

    output$info <- renderText({
      candidates <- groups$peakmat[groups$peakmat[,'group']==input$inds,]
      mz_s <- min(candidates[,'mzmin'])
      mz_e <- max(candidates[,'mzmax'])
      paste('m/z from ', mz_s, 'to', mz_e)
    })
  }
  shinyApp(ui = ui, server = server)
}

viewAlign <- function(groups_raw, groups_align){

  peakmat_raw <- as.data.table(groups_raw$peakmat)
  setkey(peakmat_raw, group)
  peakmat_align <- as.data.table(groups_align$peakmat)
  setkey(peakmat_align, group)

  ui <- fluidPage(
    titlePanel("groups Viewer"),
    sidebarLayout(
      sidebarPanel(
        sliderInput("inds",
                    "index of group:",
                    min = 1,
                    max = max(groups_raw$peakmat[,'group']),
                    value = 1,
                    step=1)
      ),
      mainPanel(
        h4 (textOutput("info"),
            plotlyOutput("Plot"),
            textOutput("mcc1"),
            plotlyOutput("Plot1"),
            textOutput("mcc2")
        )
      )
    )
  )

  server <- function(input, output) {
    candidates_raw <- reactive({
      peakmat_raw[.(input$inds)]
    })

    candidates_align <- reactive({
      peakmat_align[.(input$inds)]
    })

    output$Plot <- renderPlotly({
      p <- plot_ly()%>%
        layout(showlegend = FALSE)
      candidates <- candidates_raw()
      for (i in 1:nrow(candidates)){
        a <- as.numeric(candidates[i,'sample'])
        b <- as.numeric(candidates[i,'index'])
        pic <- groups_raw$picset[[a]]$pics[[b]]
        p <- add_trace(p, x=pic[,1], y=pic[,2], type = 'scatter', mode = 'lines')
      }
      p
    })

    output$Plot1 <- renderPlotly({
      p <- plot_ly()%>%
        layout(showlegend = FALSE)
      candidates <- candidates_align()
      for (i in 1:nrow(candidates)){
        a <- as.numeric(candidates[i,'sample'])
        b <- as.numeric(candidates[i,'index'])
        pic <- groups_align$picset[[a]]$pics[[b]]
        p <- add_trace(p, x=pic[,1], y=pic[,2], type = 'scatter', mode = 'lines')
      }
      p
    })

    output$mcc1 <- renderText({
      picset <- groups_raw$picset
      gpi <- candidates_raw()
      sams <- unlist(gpi[,'sample'])
      inds <- unlist(gpi[,'index'])
      apics <- lapply(1:nrow(gpi),function(s){
        picset[[sams[s]]]$pics[[inds[s]]]
      })
      paste('The mcc of raw pics: ', round(.cal_mcc(apics),4))
    })

    output$mcc2 <- renderText({
      picset <- groups_align$picset
      gpi <- candidates_align()
      sams <- unlist(gpi[,'sample'])
      inds <- unlist(gpi[,'index'])
      apics <- lapply(1:nrow(gpi),function(s){
        picset[[sams[s]]]$pics[[inds[s]]]
      })
      paste('The mcc of aligned pics: ', round(.cal_mcc(apics),4))
    })

    output$info <- renderText({
      candidates <- candidates_raw()
      mz_s <- min(candidates[,'mzmin'])
      mz_e <- max(candidates[,'mzmax'])
      paste('m/z from ', mz_s, 'to', mz_e)
    })
  }
  shinyApp(ui = ui, server = server)
}

viewPseudospecturm <- function(groups){

  ui <- fluidPage(
    titlePanel("Pseudospecturm Viewer"),
    sidebarLayout(
      sidebarPanel(
        sliderInput("inds",
                    "index of scan:",
                    min = 1,
                    max = max(groups$group.info[,'cluster']),
                    value = 1,
                    step=1)
      ),
      mainPanel(
        h4(textOutput("caption")),
        plotOutput("Plot"),
        tableOutput("View")
      )
    )
  )

  server <- function(input, output) {
    output$caption <- renderText({
      paste('Pseudospecturm of cluster: ', input$inds)
    })

    output$Plot <- renderPlot({
      spectrumA <- getPseudospecturm(groups, input$inds)
      stem(spectrumA[,'mz'], spectrumA[,'intensity'])
    })

    output$View <- renderTable({
      spectrumA <- getPseudospecturm(groups, input$inds)
      as.data.frame(spectrumA)
    })
  }

  shinyApp(ui = ui, server = server)
}
