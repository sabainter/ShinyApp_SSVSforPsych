options(digits=9)
options(warn=-1) #Suppress Warnings
options(shiny.maxRequestSize=20*1024^2) #max file size = 20 MB
require(foreign)
library(shiny)
library(data.table)
library(ggplot2)
library(dplyr)
library(formattable)
library(shinyalert)

# Comment

# Thomas' comment

# toKeep function--------------

toKeep <- function(x) is.numeric(x) | is.factor(x)

# SSVS Function--------------

SSVS<-function(y,x,xp,
                             runs=20000,burn=5000,update=1000,
                             a1=0.01,b1=0.01,prec.beta=0.1,inprob=0.5){
  
  n  <- length(y)
  np <- nrow(xp)
  p  <- ncol(x)
  
  #initial values:
  
  int   <- mean(y)
  beta  <- rep(0,p)
  alpha <- rep(0,p)
  delta <- rep(0,p)
  taue  <- 1/var(y)
  
  #keep track of stuff:
  
  keep.beta           <- matrix(0,runs,p)
  colnames(keep.beta) <- colnames(x)
  keep.int<-keep.taue <- rep(0,runs)
  keep.yp             <- matrix(0,runs,np)
  
  #LET'S ROLL:
  for(i in 1:runs){
    
    taue  <- rgamma(1,n/2+a1,sum((y-int-x%*%beta)^2)/2+b1)
    int   <- rnorm(1,mean(y-x%*%beta),1/sqrt(n*taue))
    
    #update alpha
    z     <- x%*%diag(delta)
    V     <- solve(taue*t(z)%*%z+prec.beta*diag(p))
    M     <- taue*t(z)%*%(y-int)
    alpha <- V%*%M+t(chol(V))%*%rnorm(p)
    beta  <- alpha*delta
    
    #update inclusion indicators: 
    r <- y-int-x%*%beta
    for(j in 1:p){
      r         <- r+x[,j]*beta[j]
      log.p.in  <- log(inprob)-0.5*taue*sum((r-x[,j]*alpha[j])^2)
      log.p.out <- log(1-inprob)-0.5*taue*sum(r^2)
      diff      <- log.p.in-log.p.out
      diff      <- ifelse(diff>10,10,diff)
      p.in      <- exp(diff)/(1+exp(diff))
      delta[j]  <- rbinom(1,1,p.in)
      beta[j]   <- delta[j]*alpha[j]
      r         <- r-x[,j]*beta[j]
    }
    
    #Make predictions:
    yp <- rnorm(np,int+xp%*%beta,1/sqrt(taue))
    
    #Store the output:
    keep.beta[i,] <- beta
    keep.int[i]   <- int
    keep.taue[i]  <- taue
    keep.yp[i,]   <- yp
    
    #if(i%%update==0){
    #  plot(beta,main=paste("Iteration",i))
    #  abline(0,0)
    #}
    
    # Incrementally track the progress in the Shiny App
    incProgress(1/runs)

  }
  
  list(beta = keep.beta[burn:runs,],
       int  = keep.int[burn:runs],
       taue = keep.taue[burn:runs],
       pred = keep.yp[burn:runs,])}

# User interface-------------

ui <- shinyUI(
  
  fluidPage(
    #useShinyalert()  # Instantiate shiny alert. LIKELY WILL NOT USE,
    tags$style(type="text/css",
               ".shiny-output-error { visibility: hidden; }",
               ".shiny-output-error:before { visibility: hidden; }",
               "h1, h3 {
               text-align:center;
               }"
 
    ),
    
    titlePanel(h1("SSVSforPsych",
                  h3("An interactive web tool for performing stochastic search variable selection"))),
  
  sidebarPanel(
    
    #tags$strong("How to cite"),
    
    #p("This web tool may be cited in APA style in the following manner:"),
    
    #p("Bainter, S. A., McCauley, T. G., Wager, T., & Losin, E. A. R. (2020). Improving practices for selecting a subset of important predictors in psychology: An application to predicting pain. Advances in Methods in Psychological Science, XX(X), XXXX-XXXX. https://doi.org/10.1177/2515245919885617"),
    
    #tags$strong("The purpose of SSVS, and how to use this tool"),
    
    #p("The overall goal of SSVS is to provide information about the relative importance of predictors, accounting for uncertainty in which other predictors are included in the model. SSVS samples thousands of regression models in order to characterize the model uncertainty regarding both the predictor set and the regression parameters. The models are selected using a sampling process that is designed to select among “good” models, that is, models with high probability."),
      
    #p("After sampling, rather than selecting the best model according to a specified criterion (e.g., the best Akaike’s or Bayesian information criterion or the highest model R2), researchers can examine the proportion of times each predictor was selected, which provides information about which predictors reliably predict the outcome, accounting for uncertainty in the other predictors in the model. Please see Bainter, McCauley, Wager, and Losin (2020) for more details."),
    
    #p("The key quantity obtained by SSVS is the marginal inclusion probability (MIP), which is the proportion of times each predictor was included in the sampled models. Predictors with higher MIPs are consistent predictors of the dependent variable, accounting for uncertainty in the other variables included in the models."),
    
    tags$strong("Pre-requisites and constraints"),
    
    p("Unfortunately, SSVS cannot handle variables with a categorical data structure. In addition, SSVS requires predictor variables to be standardized before analysis: because the priors have a fixed scale, predictors on different scales will differentially influence results. Because of this, SSVSforPsych automatically standardizes all predictors selected for analysis. Finally, SSVS cannot analyze data with missing values, so please upload data that has complete cases."),
    
    #htmlOutput("missingNote"),
    
    fileInput('file1', 
              'Choose File (max size = 30 MB)',
              accept=c('text/csv',
                       'text/comma-separated-values,text/plain',
                       '.csv',
                       '.xls',
                       '.xlsx',
                       '.sav')),
  
    # Prior inclusion probablility
    tags$strong("Prior inclusion probablility"),
    p("Please select the prior inclusion probablility value. Please consider the number of variables that you have when selecting the prior. For example, if you include 75
predictors and set the prior inclusion probability at 0.50 (as Bainter et al. (2020) do in their empirical analysis), this implies the prior belief that 37.5 predictors belong in the true model. The value selected will influence the magnitude of the MIPs, but the relative pattern of MIPS should remain fairly consistent."),
    uiOutput("ui.PriorProb"),
    tags$hr(),
  
    # Burn-in iterations
    tags$strong("Burn-in iterations"),
    p("Please select the number of burn-in iterations. Burn-ins are the number of discarded warmup iterations used to achieve  Markov chain Monte Carlo (MCMC) estimation convergence. You may increase the number of burn-in iterations if you are having convergence issues."),
    uiOutput("ui.BurnIn"),
    tags$hr(),
    
    # Iterations
    tags$strong("Total number of iterations (inclusive of burn-in)"),
    p("Please select the total number of iterations. The number of iterations indicates the number of models sampled by SSVS. We recommend 10,000 iterations as a reasonable number of iterations that balances precision and computational efficiency."),
    uiOutput("ui.Runs"),
    tags$hr(),
    
    # Fixed or random start values
    tags$strong("Fixed or random start values"),
    p("Please indicate whether you want to use fixed or random start values in your analysis. Fixed values will result in a perfectly reproducible solution produced in each run. Random values will result in slighly variant, but highly similar solutions produced in each run."),
    uiOutput("ui.randomFixed"),
    tags$hr(),
    
    # Dependent variable
    tags$strong("Dependent variable"),
    p("Please select the dependent variable in your analysis."),
    uiOutput("ui.Dependent"),
    tags$hr(),
    
    # Predictor variables
    tags$strong("Predictor variables"),
    p("Please select the predictor variables in your analysis."),
    uiOutput("ui.Predictors"),
    tags$hr(),
    
    # Go button
    uiOutput("ui.go"),
    tags$hr(),
    
    # Acknowledgements
    tags$strong("Acknowledgements"),
    p("Thank you to Dr. Brian Reich of NCSU (https://www4.stat.ncsu.edu/~reich/) for posting R code for SSVS, which was used to help build this app."),
    
    # Contact us
    tags$strong("Contact Us"),
    p("If you encounter any problems, please contact us at ssvsforpsych@gmail.com.")
    
   # htmlOutput("citation")
   
  ),
  
  mainPanel(
    
    # Output: Data file ----
    tabsetPanel(id = "inTabset",
                tabPanel("Background",
                         h4("The purpose of SSVS, and how to use this tool:"),
                         
                         p("The overall goal of SSVS is to provide information about the relative importance of predictors, accounting for uncertainty in which other predictors are included in the model. SSVS samples thousands of regression models in order to characterize the model uncertainty regarding both the predictor set and the regression parameters. The models are selected using a sampling process that is designed to select among “good” models, that is, models with high probability."),
                         
                         p("After sampling, rather than selecting the best model according to a specified criterion (e.g., the best Akaike’s or Bayesian information criterion or the highest model R2), researchers can examine the proportion of times each predictor was selected, which provides information about which predictors reliably predict the outcome, accounting for uncertainty in the other predictors in the model. Please see Bainter, McCauley, Wager, and Losin (2020) for more details."),
                         
                         p("The key quantity obtained by SSVS is the marginal inclusion probability (MIP), which is the proportion of times each predictor was included in the sampled models. Predictors with higher MIPs are consistent predictors of the dependent variable, accounting for uncertainty in the other variables included in the models."),
                         h4("How to cite:"),
                         
                         p("This web tool may be cited in APA style in the following manner:"),
                         
                         p("Bainter, S. A., McCauley, T. G., Wager, T., & Losin, E. A. R. (2020). Improving practices for selecting a subset of important predictors in psychology: An application to predicting pain. Advances in Methods in Psychological Science, XX(X), XXXX-XXXX. https://doi.org/10.1177/2515245919885617"),
                         
                ),
                tabPanel("Data descriptives",
                         tableOutput("varnames")),
                         #formattableOutput("colors"), # Table with color formatting (red for factors) but I don't think we need anymore
                tabPanel("SSVS analysis results",
                         h4("Status:"),
                         verbatimTextOutput('contents'),
                         downloadButton('download', 'Download table of results'),
                         tableOutput("resultsDF")),
                tabPanel("SSVS plots",
                         downloadButton('downloadPlot.final', 'Download plot'),
                         plotOutput("resultsPlot"))
    ) # tabsetPanel end
  ) # mainPanel end
)) # UI end

server<-function(input, output, session) {
  values<-reactiveValues()
  
  # Upload the data---------------------------------------------------
  
  df<-reactive({
    if (!is.null(input$file1)){
      if(tools::file_ext(input$file1) == "csv"){
        read.csv(input$file1$datapath, header = TRUE)
      }
      else if (tools::file_ext(input$file1) == "xlsx"){
        file.rename(input$file1$datapath,
                    paste(input$file1$datapath, ".xlsx", sep=""))
        as.data.frame(readxl::read_excel(paste(input$file1$datapath, ".xlsx", sep=""), 1,col_names = TRUE))
      }
      else if (tools::file_ext(input$file1) == "xls"){
        file.rename(input$file1$datapath,
                    paste(input$file1$datapath, ".xls", sep=""))
        as.data.frame(readxl::read_xls(paste(input$file1$datapath, ".xls", sep=""), 1,col_names = TRUE))
      }
      else if (tools::file_ext(input$file1) == "sav"){
        foreign::read.spss(input$file1$datapath, to.data.frame=TRUE, use.value.labels = F)
      }
    }
  })
  
  #Only allow users to include data that is numeric------------------------------

  #Produce a table of data descriptives------------------------------
  
  output$varnames <- renderTable({
    if (is.null(input$file1))
      return()
    desc <- as.data.frame(psych::describe(df()))[,c("mean",
                                                    "sd",
                                                    "min",
                                                    "max",
                                                    "skew")]
    varClass <- as.vector(sapply(df(),class))
    varmissing <- sapply(df(), function(x) sum(is.na(x)))
    
    varlist <- cbind(colnames(df()),
                     desc,
                     varClass,
                     varmissing)
    colnames(varlist) <- c("Variable",
                           colnames(desc),
                           "Object_type",
                           "NA values")
    varlist
  })
  
  # Below is the code for the color formatted table
  # output$colors <- renderFormattable({
  #   if (is.null(input$file1))
  #     return()
  #   desc <- as.data.frame(psych::describe(df()))[,c("mean",
  #                                                   "sd",
  #                                                   "min",
  #                                                   "max",
  #                                                   "skew")]
  #   desc[,] <- round(desc[,],2)
  #   varClass <- as.vector(sapply(df(),class))
  #   varlist <- cbind(colnames(df()),
  #                    desc,
  #                    varClass)
  #   colnames(varlist) <- c("Variable",
  #                          colnames(desc),
  #                          "Object_type")
  #   varlist <- as.data.frame(varlist)
  #   row.names(varlist) <- c()
  #   formattable(varlist,list(
  #   Variable = formatter(.tag = "span", width = '1px'),
  #   Object_type = formatter("span", 
  #                       style = x ~ style(
  #                           color = ifelse(x == 'factor','red','black')))
  #   ))
  # })
  # 
  
  # Status bar
  output$contents <- renderText({
    varClass <- as.vector(sapply(df(),class))
    if (any(varClass == 'factor')){
      return("Warning: Factor variables are not allowed!")
    }
    return("Waiting for analysis to begin")
  })
  
  # Select the predictor variables------------------------------
  
  output$ui.Predictors <- renderUI ({
    if (is.null(input$file1))
      return()
  
  # Includes an option for select all so that all variables can be simultaneously selected
  choices = c("Select all", colnames(df())) 
  observe({
    if("Select all" %in% input$Predictors)
      selected_choices = choices[-1] # Choose all of the choices except for "Select all"
    else
      selected_choices = input$Predictors
    updateSelectInput(session, "Predictors", selected = selected_choices)
  
  })
  selectInput(inputId = "Predictors",
              label = "",
              choices = choices,
              selected = "",
              multiple = TRUE,
              selectize = TRUE)
  
 
  })
  

  # Select the dependent variable------------------------------

  output$ui.Dependent <- renderUI ({
    if (is.null(input$file1))
      return()

  selectInput(inputId = "Dependent",
              label = "",
              choices = colnames(df()),
              selected = "")

  })
  
  # Select the prior probablity value------------------------------
  
  output$ui.PriorProb <- renderUI ({
    
  selectInput(inputId = "PriorValue", label = "",
                 choices = c("0.1", "0.2", "0.3", "0.4","0.5", "0.6", "0.7", "0.8", "0.9"),
                 selected = "0.5")
    
  })
  
  # Select the burn-in value------------------------------
  
  output$ui.BurnIn <- renderUI ({
    
    selectInput(inputId = "BurnInValue", label = "",
                choices = c("1000", "2000", "3000", "4000","5000", "6000", "7000", "8000", "9000", "10000"),
                selected = "1000")
    
  })
  
  # Select the runs value------------------------------
  
  output$ui.Runs <- renderUI ({
    
    selectInput(inputId = "RunsValue", label = "",
                choices = c("10000", "11000", "12000", "13000", "14000", "15000", "16000", "17000", "18000","19000", "20000"),
                selected = "10000")
    
  })
  
  # Select whether a fixed or random start value will be used------------------------------
  
  output$ui.randomFixed <- renderUI ({
    radioButtons(inputId="randomFixed", label = NULL, inline = TRUE,
                 c("Fixed" = "fix", "Random" = "rnd"),
                 selected = "Fixed")
  })
  
  # Create an object that returns NAs------------------------------
  
  # TESTING TESTING TESTING TESTING TESTING
  output$ui.navals<-renderUI ({
    if (is.null(input$file1))
      return()
    
    textInput(inputId="naval",
              label="Missing (NA) values code:",
              value = "NA",
              width = "100%")  
  })
  
  # Create a button that runs the analysis------------------------------
  
  output$ui.go <- renderUI ({
    if (is.null(input$file1))
      return()
    actionButton("go","Run Analysis", width = "100%", icon = icon("fa fa-thumbs-up"))
  })
  
  # Create a citation for the results------------------------------
  
  # output$citation<-renderUI({
  #   str1 <- paste0("Acknowledgements:")
  #   str2 <- paste0("Thank you to Dr. Brian Reich of NCSU (https://www4.stat.ncsu.edu/~reich/) for posting R code for SSVS, which was used to help build this app. ")
  #   str3 <- paste0("")
  #   str4 <- paste0("Preprint available at https://psyarxiv.com/j8t7s/")
  #   str5 <- paste0("")
  #   str6 <- paste0("If you encounter any problems, please contact us at ssvsforpsych@gmail.com")
  #   
  #   HTML(paste(str1, str2, str3, str4, str5, str6, sep = '<br/><br/>'))
  # })
  
  # output$missingNote<-renderUI({
  #   str1 <- paste0("Please upload data that does not have any missing values, as the SSVS application is not equipped to handle missing data.")
  #   
  #  # HTML(paste(str1,sep = '<br/><br/>'))
  # })
    
  # Analyses code -----------------------------------------------------------

  observeEvent(input$go, { #Once the "go" button is hit, InterActive looks at all the ui input and runs the model.
   
  dftrue <- df() #retain an original copy; REPLACED THE ORIGINAL CODE: dftrue<-plotelements$data<-df()

  # TESTING TESTING TESTING TESTING TESTING
  # Replace all missing values in dataset with NA
  reactive ({
    dftrue[dftrue==input$naval] <- NA
  })
    
  mydata <- dftrue #make a copy of the data for variable rescaling and analyses

  # Create a data.frame of the data??? Not really sure what the utility of this function is.
  mydata <- as.data.frame(data.matrix(mydata))
  
  # TESTING TESTING TESTING TESTING TESTING
  # Only include data that is numeric or factor
  mydata %>% select_if(toKeep)
  
  # define predictor variables
  preds <- input$Predictors
  
  # define outcome variable
  dependent <- input$Dependent

  # Scale the predictors 
  mydata[,preds] <- scale(mydata[,preds],  scale = FALSE)
  

  # Create matrices for the independent and dependent variables to use in SSVS.
  
  # predictor variables
  values$predsSSVS <- as.matrix(mydata[,preds])
  
  # outcome variable
  values$dependentSSVS <- as.matrix(mydata[,dependent])
  
  # Pop up alert using shinyalert (LIKELY WILL NOT USE)
  # if (sum(is.na(values$dependentSSVS))>0){
  #   shinyalert(paste(c("Oops!", sum(is.na(values$dependentSSVS)),"missing values found. Would you like to use only complete cases? This will reduce your sample size to",length(complete.cases(values$dependentSSVS)),"from",length(values$dependentSSVS)))
  #                       ,type = "error")
  #   values$dependentSSVS <- complete.cases(values$dependentSSVS)
  # }
  
  # Pop up window using base Shiny (PLACEHOLDER CODE)
  dataModal <- function(failed = FALSE) {
    modalDialog(
      textInput("dataset", "Choose data set",
                placeholder = 'Try "mtcars" or "abc"'
      ),
      span('(Try the name of a valid data object like "mtcars", ',
           'then a name of a non-existent object like "abc")'),
      if (failed)
        div(tags$b("Invalid name of data object", style = "color: red;")),
      
      footer = tagList(
        modalButton("Cancel"),
        actionButton("ok", "OK")
      )
    )
  }
  
  # Activate pop-up window if dependent variable has missing values
  if (sum(is.na(values$dependentSSVS))>0){
  showModal(dataModal())
  }
  
  # randomFixed set.seed
  if (input$randomFixed == "fix"){ 
    set.seed(0820)
  }
  
  ##### Running the SSVS function
  
  # Run the SSVS analysis -----------------------------------------------------------
  
  # Setup the analysis
  # set.seed(0820)
  n <- nrow(values$predsSSVS)
  p <- ncol(values$predsSSVS)
  xp     <- matrix(0,25,as.numeric(p))
  xp[,1] <- seq(-3,3,length=25)
  values$xp <- xp

  
  # Create a table of the results -----------------------
  output$resultsDF <- renderTable({
    if (is.null(input$file1))
      return()
    
    # I may need to add a button here for is.null(input$go)
    withProgress(message = "Please wait, SSVS analysis is running", value=0.1, {
      
      values$ssvsResults <- SSVS(y = values$dependentSSVS, x = values$predsSSVS, xp = values$xp, runs = as.numeric(input$RunsValue), burn = as.numeric(input$BurnInValue), inprob = as.numeric(input$PriorValue))
      
      ssvsResults<-values$ssvsResults
      
      # Feed the non-reactive results into a table, and make that table a reactive value
      resultsTable<-as.data.frame(apply(ssvsResults$beta!=0,2,mean))
      resultsTable$var<-rownames(resultsTable)
      resultsTable$DV<-as.character(input$Dependent)
      names(resultsTable)<-c("Inclusion_probability","Variable_name","Dependent_variable")
      resultsTable<-resultsTable[order(-resultsTable$Inclusion_probability),]
      
      # Save the results as a reactive value
      values$resultsTable<-resultsTable
       }
    )
    
    # Print the results in a table
    values$resultsTable
    
  })
  
  # Updating the status bar

  #values$dependentSSVS
  
  output$contents <- renderText({
    varClass <- as.vector(sapply(df(),class))
    #return(str(values$dependentSSVS))
    # if (class(values$dependentSSVS == 'factor')){
    #   return("Error. Please convert factor variables to perform the analysis.")
    # }
     if (sum(is.na(values$dependentSSVS))>0){
       return(paste(c("Error.",sum(is.na(values$dependentSSVS)),"missing values found. Please remove to perform the analysis")))
     }
     if (ncol(values$predsSSVS) <= 1){
       return("Error. Please select at least two predictors")
     }
     return(paste(c("Calculation complete.",input$RunsValue, "MCMC iterations run, results for",as.numeric(input$RunsValue)
                    -as.numeric(input$BurnInValue), "iterations post-warmup shown below for",length(complete.cases(values$dependentSSVS)),"complete cases.")))
  })

  # Download the results table -----------------------
  
  output$download <- downloadHandler(
    filename = function(){"SSVS results.csv"}, 
    content = function(x){
      write.csv(values$resultsTable, x)
    })
  
  

  
  
  # Create a plot of the results -----------------------
  output$resultsPlot <- renderPlot({
    
    # Read in the results that are saved as reactive values
    ssvsResults<-values$ssvsResults
    
    # We'll need to recreate a dataframe of the results here
    plotDF<-as.data.frame(apply(ssvsResults$beta!=0,2,mean))
    plotDF$var<-rownames(plotDF)
    plotDF$DV<-as.character(input$Dependent)
    names(plotDF)<-c("Inclusion_probability","Variable_name","Dependent_variable")
    plotDF<-plotDF[order(-plotDF$Inclusion_probability),]
    
    plotDF$threshold<-ifelse(plotDF$Inclusion_probability>.5, 1, 0)
    plotDF$threshold<-as.factor(plotDF$threshold)
    levels(plotDF$threshold) <- c('< .50', '> .50')
    
    temp.Plot<-ggplot(data=plotDF) +
      geom_point(aes(x = plotDF$Variable_name, 
                     y = plotDF$Inclusion_probability, 
                     shape = plotDF$threshold), 
                 size = 2) +
      labs(y = "Inclusion Probability", 
           x = "Predictor variables",
           shape = "MIP threshold", 
           title = paste("Inclusion Probability for", as.character(input$Dependent))) +
      scale_y_continuous(limits = c(0,1.1), breaks = c(0, .25, .5, .75, 1)) +
      theme_classic() +
      geom_vline(xintercept = nrow(plotDF)+.5, linetype = 1, size = .5, alpha = .2) +
      geom_hline(yintercept = .5, linetype = 2) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            panel.spacing = unit(0, "lines"), 
            strip.background = element_blank(),
            strip.placement = "outside")
    
    values$temp.Plot <- temp.Plot
    
    values$temp.Plot
    
  })
  
  # Download the plot of the results -----------------------
  
  output$downloadPlot.final <- downloadHandler(
    filename = 'SSVSplot.png',
    content = function(file) {
      ggsave(file,values$temp.Plot, width = 7, height = 4, units = "in")
    })
  
  }) # Make an observation for the "go" button
  
 } # Server end

shinyApp(ui = ui, server = server)