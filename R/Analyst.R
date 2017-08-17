analyst.RF <- function(labels, data){
  library(shiny)
  library(plotly)
  library(randomForest)

  labels <- as.factor(labels)
  data <- as.matrix(data)

  ui <- fluidPage(
    titlePanel("Random Forest"),
    sidebarLayout(
      sidebarPanel(
        numericInput("trees", "Number of trees to grow:", 100),
        numericInput("nodes", "Maximum number of nodes:", 5),
        selectInput("normalization", "Normalization:", choices=c('None','Sum','Interior Label')),
        selectInput("scaling", "Scale:", choices=c('center','scale','auto-scaling')),
        numericInput("interior", "Index of Interior Label (only for normalization with interior babel):", 1)
      ),
      mainPanel(
        h4 ("Confusion Matrix"),
        tableOutput("cmx"),
        h4 ("Multi-dimensional Scaling Plot"),
        plotlyOutput("mds.plot"),
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

      randomForest(data, labels, ntree=input$trees, maxnodes=input$nodes,
                   importance=TRUE, proximity=TRUE)
    })

    output$cmx <- renderTable({
      as.data.frame(model()$confusion)
    })

    output$mds.plot <- renderPlotly({
      points <- MDSplot(model(), labels, k=2)$points
      p <- plot_ly(x=points[,1], y=points[,2], color=labels) %>%
        layout(xaxis = list(title = 'Dim 1'),
               yaxis = list(title = 'Dim 2'))
      p
    })

    output$imp.plot <- renderPlot({
      varImpPlot(model(), main=NULL)
    })
  }

  shinyApp(ui = ui, server = server)
}

analyst.OPLS <- function(labels, data){
  library(shiny)
  library(plotly)
  library(ropls)
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
