analyst.RF <- function(labels, data){
  
  labels <- as.factor(labels)
  data <- as.matrix(data)
  
  ui <- navbarPage("Random Forest",
                   tabPanel(
                     "Over View",
                     sidebarLayout(
                       sidebarPanel(
                         numericInput("trees", "Number of trees to grow:", 100),
                         numericInput("nodes", "Maximum number of nodes:", 5),
                         selectInput("normalization", "Normalization:", choices=c('None','Sum','Interior Label')),
                         selectInput("scaling", "Scale:", choices=c('None', 'Center','Scale','Auto-scaling')),
                         numericInput("interior", "Index of Interior Label (only for normalization with interior babel):", 1)
                       ),
                       mainPanel(
                         h4 ("Random Forest classification"),
                         tableOutput("cmx"),
                         h4 ("Error Rates"),
                         plotlyOutput("mod.plot")
                       )
                     )
                   ),
                   tabPanel(
                     "MDP plot",
                     sidebarLayout(
                       sidebarPanel(
                         numericInput("mds.num", "Number of dimensions for the scaling coordinates:", 2)
                       ),
                       mainPanel(
                         h4 ("Multi-dimensional Scaling Plot"),
                         plotOutput("mds.plot")
                       )
                     )
                   ),
                   tabPanel(
                     "Variable Importance",
                     sidebarLayout(
                       sidebarPanel(
                         numericInput("vip.num", "How many variables to show ?:", 10)
                       ),
                       mainPanel(
                         h4 ("Variable Importance"),
                         plotOutput("imp.plot"),
                         tableOutput("imp.table")
                       )
                     )
                   ),
                   tabPanel(
                     "Trees",
                     sidebarLayout(
                       sidebarPanel(
                         numericInput("tree.num", "which tree to extract ?", 1)
                       ),
                       mainPanel(
                         tableOutput("tree.print")
                       )
                     )
                   )
  )
  
  server <- function(input, output) {
    model <- reactive({
      ColNames <- paste('mz:', data['mz',], 'rt:', data['rt',])
      data <- data[-(1:6),]
      colnames(data) <- ColNames
      
      if (input$normalization=='Sum'){
        sums <- rowSums(data)
        data <- t(t(data)/sums)
      }else if (input$normalization=='Interior Label'){
        interior <- data[,input$interior]
        data <- t(t(data)/interior)
      }
      
      if (input$scaling == 'Center'){
        data <- scale(data, center=TRUE, scale=FALSE)
      } else if (input$scaling == 'Scale'){
        data <- scale(data, center=FALSE, scale=TRUE)
      } else if (input$scaling == 'Auto-scaling'){
        data <- scale(data, center=TRUE, scale=TRUE)
      }
      
      randomForest(data, labels, ntree=input$trees, maxnodes=input$nodes,
                   importance=TRUE, proximity=TRUE)
    })
    
    output$cmx <- renderTable({
      as.data.frame(model()$confusion)
    })
    
    output$mod.plot <-  renderPlotly({
      points <- as.data.frame(plot(model()))
      p <- plot_ly() %>%
        layout(xaxis = list(title = 'Trees'),
               yaxis = list(title = 'Error'))
      for (i in 1:ncol(points)){
        p <- add_trace(p, x = 1:nrow(points), y = points[,i], name = colnames(points)[i], type = 'scatter', mode = 'lines')
      }
      p
    })
    
    output$mds.plot <- renderPlot({
      MDSplot(model(), labels, k=input$mds.num)$points
    })
    
    output$imp.plot <- renderPlot({
      varImpPlot(model(), n.var = min(input$vip.num, nrow(model()$importance)), main=NULL)
    })
    
    output$imp.table <- renderTable({
      imp.tab <- importance(model())
      as.data.frame(imp.tab[1:min(nrow(imp.tab), input$vip.num),])
    })
    
    output$tree.print <-  renderTable({
      getTree(model(), k=input$tree.num, labelVar=TRUE)
    })
    
  }
  
  shinyApp(ui = ui, server = server)
}

analyst.OPLS <- function(labels, data){
  labels <- as.vector(labels)
  data <- as.matrix(data)

  ui <- fluidPage(
    titlePanel("PLS or OPLS classification"),
    sidebarLayout(
      sidebarPanel(
        numericInput("predI", "Number of components:", 2),
        numericInput("orthoI", "Number of orthogonal components (for OPLS only):", 1),
        selectInput("method", "Method:", choices=c('PLS','OPLS')),
        selectInput("normalization", "Normalization:", choices=c('Sum','Interior Label')),
        selectInput("scaling", "Scale:", choices=c('none','center','scale','auto-scaling')),
        numericInput("interior", "Index of Interior Label (only for normalization with interior babel):", 1)
      ),
      mainPanel(
        h4 ("Overview"),
        tableOutput("overview.plot"),
        h4 ("Score Plot"),
        plotlyOutput("score.plot"),
        h4 ("Variable Importance"),
        plotOutput("imp.plot")
      )
    )
  )

  server <- function(input, output) {
    model <- reactive({
      ColNames <- paste('mz:', data['mz',], 'rt:', data['rt',])
      data <- data[-(1:6),]
      colnames(data) <- ColNames

      if (input$normalization=='Sum'){
        sums <- rowSums(data)
        data <- t(t(data)/sums)
      }else if (input$normalization=='Interior Label'){
        interior <- data[,input$interior]
        data <- t(t(data)/interior)
      }

      if (input$scaling == 'center'){
        data <- scale(data, center=TRUE, scale=FALSE)
      } else if (input$scaling == 'scale'){
        data <- scale(data, center=FALSE, scale=TRUE)
      } else if (input$scaling == 'auto-scaling'){
        data <- scale(data, center=TRUE, scale=TRUE)
      }

      if (input$method == 'PLS') {
        f <- opls(data, labels, predI=input$predI)
      } else if (input$method == 'OPLS'){
        f <- opls(data, labels, predI=input$predI, orthoI=input$orthoI)
      }
      f
    })

    output$overview.plot <- renderTable({
      res <- cbind(rownames(model()@modelDF), model()@modelDF)
      colnames(res)[1] <- ' '
      as.data.frame(res)
    })


    output$score.plot <- renderPlotly({
      if (input$method == 'PLS'){
        scores <- model()@scoreMN
        p <- plot_ly(x=scores[,1], y=scores[,2], color=labels) %>%
          layout(xaxis = list(title = paste('Component 1 (', round(model()@modelDF$R2X[1]*100, 2), '%)')),
                 yaxis = list(title = paste('Component 2 (', round(model()@modelDF$R2X[2]*100, 2), '%)')))
      } else if (input$method == 'OPLS'){
        comp1 <- model()@scoreMN[,1]
        comp2 <- model()@orthoScoreMN[,1]
        p <- plot_ly(x=comp1, y=comp2, color=labels) %>%
          layout(xaxis = list(title = paste('T-score')),
                 yaxis = list(title = paste('Orthogonal T-score')))
      }
      p
    })

    output$imp.plot <- renderPlot({
      imp <- getVipVn(model())
      if (is.null(names(imp))){names(imp) <- 1:length(imp)}
      tops <- sort(imp, decreasing=T)
      tops <- tops[1:min(20,length(tops))]
      tops <- sort(tops)
      dotchart(tops,names(tops))
    })
  }

  shinyApp(ui = ui, server = server)
}
