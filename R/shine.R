
##################################################
shinierInput <- function(FUN, id, namez, ...){
  inputs = character(length(namez))
  for (i in 1:length(namez))
  {
    inputs[i] = as.character(FUN(paste0(id, namez[i]), ...))
  }
  inputs
}

# returns c(nrow,ncol) based on npic to minimise whitespace in plot display
# based on the assumption of slightly more available width than height
matrix_layout = Vectorize(
  function(npic){
  # larger than 3,3 seems not supported
    if(npic == 1) return(c(1,1))
    if(npic == 2) return(c(1,2))
    if(npic <= 4) return(c(2,2))
    if(npic <= 6) return(c(2,3))
    return(c(3,3))
})


###########################################################  
#' Interactive test-item analysis
#'
#' Open a shiny application for interactive item-test
#' analysis on the database
#'
#'
#' @param db A handle to the database, i.e. the output of \code{create_new_project}
#' or \code{open_project}
#' @return
#' An object that represents the application. Printing the object or passing it to \code{shiny::runApp} will run the app.
iTIA <- function(db) {
  
  message('processing data - this can take some time depending on the size of your data.')
  respData = get_resp_data(db, extra_columns='response', summarised=FALSE)
  
  theTIA = tia_tables(respData, type='averaged')
  
  ourItems = theTIA$itemStats$item_id
  ourBooklets = theTIA$testStats$booklet_id
  
  # do some rounding and aesthetic renaming
  theTIA$itemStats = theTIA$itemStats %>%
    mutate(pvalue = round(.data$pvalue,3), rit = round(.data$rit,3), rir = round(.data$rir,3), 
           meanScore = round(.data$meanScore,2), sdScore = round(.data$sdScore,2)) %>%
    rename(Rit = 'rit', Rir = 'rir', Pvalue = 'pvalue')
  
  theTIA$testStats = theTIA$testStats %>%
    mutate(alpha = round(.data$alpha,3), meanP = round(.data$meanP,3), meanRit = round(.data$meanRit,3), meanRir = round(.data$meanRir,3))
  
  
  sidebar = dashboardSidebar(
    sidebarMenu(
      menuItem("Booklets", tabName = "booklets"),
      menuItem("Items", tabName = "items")
    )
  )
  
  body = dashboardBody(
    tabItems(
      tabItem(tabName = "booklets", dataTableOutput("booklets"),
              lapply(ourBooklets,
                     function(i){
                       bsModal(paste0("myModal", i),
                               # "Local independence (left), Interaction model (right)",
                               "Models", 
                               paste0("btn", i),
                               size = "large",
                               plotOutput(paste0("plot", i)))
                     })
      ),
      tabItem(tabName = "items", dataTableOutput("items"),
              lapply(ourItems, function(i){
                bsModal(paste0("myModal", i),
                        "Distractor plot",
                        paste0("btn", i),
                        size = "large",
                        plotOutput(paste0("plot", i)))
              })
      )
    )
  )
  
  # Put them together into a dashboardPage
  ui = dashboardPage(
    dashboardHeader(title = "Interactive Test-Item Analysis"),
    sidebar,
    body
  )
  
  server <- function(input, output, session) {
    tia = reactive({
      theTIA
    })
    
    lapply(ourItems, function(i){
      lo = matrix_layout(nrow(filter(respData$design, .data$item_id==i)))
      # renderPlot automatically delays execution
      output[[paste0("plot", i)]] = renderPlot(distractor_plot(respData,i,nr=lo[1],nc=lo[2]))
      observeEvent(input[[paste0("btn", i)]], {
        toggleModal(session, paste0("myModal", i), "open")
      })
    })
    
    lapply(ourBooklets, function(i){
      # we delay execution of fit_inter until the moment the plot is requested
      delayedAssign('mo', fit_inter(respData %>% filter(.data$booklet_id==i)))
      #mo = fit_inter(respData %>% filter(.data$booklet_id==i))
      output[[paste0("plot", i)]] = renderPlot(
        plot(mo, overlay=TRUE, nc=2)
      )
      observeEvent(input[[paste0("btn", i)]], {
        toggleModal(session, paste0("myModal", i), "open")
      })
    })
    
    output$booklets = DT::renderDataTable({
      B=tia()[[2]]
      Plots = shinierInput(actionButton, "btn", ourBooklets, label = "Show")
      B = cbind(B, Plots)
      B
    }, extensions = 'Buttons', options=list(dom = 'Bfrtip',
                                            buttons = c('copy', 'csv', 'excel', 'pdf', 'print', 'pageLength'),
                                            orderClasses = TRUE,
                                            preDrawCallback = DT::JS("function() {
                                                                     Shiny.unbindAll(this.api().table().node()); }"),
                                            drawCallback = DT::JS("function() {
                                                                  Shiny.bindAll(this.api().table().node()); } "),
                                            lengthMenu = list(c(10, 25, -1), c('10', '25', 'All')),
                                            autoWidth = TRUE, scrollX=TRUE), escape=FALSE)
    
    
    output$items = DT::renderDataTable({
      A=tia()[[1]]
      Plots <- shinierInput(actionButton, "btn", ourItems, label = "Show")
      A = cbind(A,Plots)
      A
    }, extensions = 'Buttons', options = list(dom = 'Bfrtip',
                                              buttons = c('copy', 'csv', 'excel', 'pdf', 'print', 'pageLength'),
                                              orderClasses = TRUE,
                                              pageLength = 25,
                                              preDrawCallback = DT::JS("function() {
                                                                       Shiny.unbindAll(this.api().table().node()); }"),
                                              drawCallback = DT::JS("function() {
                                                                    Shiny.bindAll(this.api().table().node()); } "),
                                              lengthMenu = list(c(10, 25, -1), c('10', '25', 'All')),
                                              autoWidth = TRUE,
                                              scrollX=TRUE
    ), escape=FALSE)
  } # end of server
  shinyApp(ui, server)
}






#############################
#' Interactive model display
#'
#' Opens up a shiny application with item statistics and interactive
#' plots for the Rasch and Interaction models
#'
#'
#' @param db A handle to the database, i.e. the output of \code{create_new_project}
#' or \code{open_project}
#' @param booklet booklet_id of the booklet that will be shown
#' @return
#' An object that represents the application. Printing the object or passing it to \code{shiny::runApp} will run the app.
iModels <- function(db, booklet) {
  
  message('processing data - this can take some time depending on the size of your data')
  
  respData = get_resp_data(db, qtpredicate = quote(booklet_id == booklet), summarised=FALSE)
  models = fit_inter(respData)
  
  tia = tia_tables(respData, type='raw')$itemStats 
  sidebar = dashboardSidebar(
    checkboxInput("show", "Show data", value = FALSE, width = NULL),
    checkboxInput("summate", "As scores", value = TRUE, width = NULL)
  )
  
  body = dashboardBody(
    dataTableOutput("items"),
    lapply(tia$item_id,
           function(i){
             bsModal(paste0("myModal", i),
                     "Models",
                     paste0("btn", i),
                     size = "large",
                     plotOutput(paste0("plot", i)))
           })
  )
  # Put them together into a dashboardPage
  ui = dashboardPage(
    dashboardHeader(title = "Interactive Models"),
    sidebar,
    body
  )
  
  server <- function(input, output, session) {
    output$items = DT::renderDataTable({
      A=tia[tia$booklet_id==booklet,c("item_id","meanScore","maxScore","pvalue","rit","rir","n")]
      A$pvalue=round(A$pvalue,3)
      A$rit=round(A$rit,3)
      A$rir=round(A$rir,3)
      A$meanScore=round(A$meanScore,2)
      Plots = shinierInput(actionButton,  "btn", tia$item_id, label = "Show")
      A = cbind(A,Plots)
      A
    },extensions = 'Buttons', options = list(
      dom = 'Bfrtip',
      buttons = c('copy', 'csv', 'excel', 'pdf', 'print', 'pageLength'),
      orderClasses = TRUE,
      pageLength = 25,
      preDrawCallback = DT::JS("function() {
                               Shiny.unbindAll(this.api().table().node()); }"),
      drawCallback = DT::JS("function() {
                            Shiny.bindAll(this.api().table().node()); } "),
      lengthMenu = list(c(10, 25, -1), c('10', '25', 'All')),
      autoWidth = TRUE,
      scrollX=TRUE
      ), escape=FALSE)
    lapply(tia$item_id, function(i){
      output[[paste0("plot", i)]] =
        renderPlot(plot(models, item=i,
                        summate=input$summate,
                        overlay=FALSE,
                        nc=1, nr=1, curtains=10,
                        show.observed=input$show
        ))
      observeEvent(input[[paste0("btn", i)]], {
        toggleModal(session, paste0("myModal", i), "open")
      })
    })
}
  shinyApp(ui, server)
}
