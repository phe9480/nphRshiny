library(shiny)
library(ggplot2)
library(DT)
library(nphRshiny)
library(tidyverse)
library(plotly)
library(mvtnorm)
library(gsDesign)
library(parallel)  
library(foreach)  
library(doParallel)
library(nphsim)

# Define server logic required to draw a histogram
function(input, output, session) {
  
  showError <- function(message) {
    # Ensure the message is a single string
    if (length(message) > 1) {
      message <- paste(message, collapse=", ")  # Collapses multiple values into one string
    }
    shiny::showNotification(as.character(message), type = "error")
  }
  
  # for test parameter
  output$inFunControl <- renderUI({
    inputs <- list()
    if (input$testFunControl == "F"){
      inputs[['conMethod']] <- radioButtons("conMethod", "Please choose one method to input. ",choices = list("Median","Survival"))
    }
    
    if (input$testFunControl == "A") {
      message <- "Please provide the function formula for the control arm. "
      inputs[['A']] <- HTML(paste("<p>", message, "</p>"))
      inputs[['H0']] <- textInput("custH0", "Hazard function of control")
      inputs[['H1']] <- textInput("custS0", "Survival function of control")
    }
    
    if (input$testFunControl == "B") {
      message2 <- "Please select the number of pieces and input the corresponding piecewise rate for control arm"
      inputs[['A']] <- HTML(paste("<p>", message2, "</p>"))
      inputs[['piece']] <- selectInput("con", "Number of pieces for control", choices = 1:5)
    }
    if (input$testFunControl == "C") {
      message3 <- "Please input the corresponding weibull parameters for control arm"
      inputs[['A']] <- HTML(paste("<p>", message3, "</p>"))
      inputs[['shape1_c']] <- textInput("shape1_c", "Shape")
      inputs[['scale1_c']] <- textInput("scale1_c", "Scale")
      
    }
    if (input$testFunControl == "D") {
      message4 <- "Please input the corresponding mixture cure rate of exponential for control arm"
      inputs[['A']] <- HTML(paste("<p>", message4, "</p>"))
      inputs[['p2_c']] <- textInput("p2_c", "Cure rate")
      inputs[['lam2_c']] <- textInput("lam2_c", "Hazard rate") 
    }
    if (input$testFunControl == "E") {
      message5 <- "Please input the corresponding mixture cure rate of weibull for control arm"
      inputs[['A']] <- HTML(paste("<p>", message5, "</p>"))
      inputs[['p3_c']] <- textInput("p3_c", "Cure rate")
      inputs[['shape3_c']] <- textInput("shape3_c", "Shape")
      inputs[['scale3_c']] <- textInput("scale3_c", "Scale") 
      
    }
    
    do.call(tagList, inputs)
  })
  
  output$inFunExperiment <- renderUI({
    inputs <- list()
    if (input$testFunExperiment == "F"){
      
      inputs[['expMethod']] <- radioButtons("expMethod", "Please choose one method to input. ",choices = list("Median","Survival"))
    }
    
    if (input$testFunExperiment == "A") {
      message <- "Please provide the function formula for the experimental arm." 
      inputs[['A']] <- HTML(paste("<p>", message, "</p>"))
      inputs[['H0']] <- textInput("custH1", "Hazard function of experiment")
      inputs[['H1']] <- textInput("custS1", "Survival function of experiment")
    }
    
    if (input$testFunExperiment == "B") {
      message2 <- "Please select the number of pieces and input the corresponding piecewise rate for experimental arm"
      inputs[['A']] <- HTML(paste("<p>", message2, "</p>"))
      inputs[['piece']] <- selectInput("exp", "Number of pieces for experiment", choices = 1:5)
    }
    if (input$testFunExperiment == "C") {
      message3 <- "Please input the corresponding weibull parameters for control arm"
      inputs[['A']] <- HTML(paste("<p>", message3, "</p>"))
      inputs[['shape1_e']] <- textInput("shape1_e", "Shape")
      inputs[['scale1_e']] <- textInput("scale1_e", "Scale")
      
    }
    if (input$testFunExperiment == "D") {
      message4 <- "Please input the corresponding mixture cure rate of exponential for control arm"
      inputs[['A']] <- HTML(paste("<p>", message4, "</p>"))
      inputs[['p2_e']] <- textInput("p2_e", "Cure rate")
      inputs[['lam2_e']] <- textInput("lam2_e", "Hazard rate") 
    }
    if (input$testFunExperiment == "E") {
      message5 <- "Please input the corresponding mixture cure rate of weibull for control arm"
      inputs[['A']] <- HTML(paste("<p>", message5, "</p>"))
      inputs[['p3_e']] <- textInput("p3_e", "Cure rate")
      inputs[['shape3_e']] <- textInput("shape3_e", "Shape")
      inputs[['scale3_e']] <- textInput("scale3_e", "Scale") 
      
    }
    
    do.call(tagList, inputs)
  })
  
  observeEvent(input$con, {
    req(input$testFunControl == "B")
    output$moreInputsControl <- renderUI({
      numCon <- as.numeric(input$con)
      controlInputs <- lapply(1:numCon, function(i) {
        fluidRow(
          column(6, textInput(paste("tc", i, sep = "_"), "Starting Time (months)", value = if (i == 1) 0 else NA)),
          column(6, textInput(paste("hc", i, sep = "_"), "Hazard Rate"))
        )
      })
      do.call(tagList, controlInputs)
    })
  })
  
  observeEvent(input$exp, {
    req(input$testFunExperiment == "B")
    output$moreInputsExperiment <- renderUI({
      numExp <- as.numeric(input$exp)
      experimentInputs <- lapply(1:numExp, function(j) {
        fluidRow(
          column(6, textInput(paste("te", j, sep = "_"), "Starting Time (months)", value = if (j == 1) 0 else NA)),
          column(6, textInput(paste("he", j, sep = "_"), "Hazard Rate"))
        )
      })
      do.call(tagList, experimentInputs)
    })
  })
  
  observeEvent(input$conMethod, {
    req(input$testFunControl == "F")
    output$medConOut <- renderUI({
      if (input$conMethod == "Median"){
        textInput("medC","Median of Control")
      } else{
        fluidRow(
          column(6,
                 textInput("surC","Survival")),
          column(6,
                 textInput("tCon","At (month)"))
        )
        
      }
      
    })
  })
  
  observeEvent(input$expMethod, {
    req(input$testFunExperiment == "F")
    output$medExpOut <- renderUI({
      if (input$expMethod == "Median"){
        textInput("medE","Median of Experimental")
      } else{
        fluidRow(
          column(6,
                 textInput("surE","Survival")),
          column(6,
                 textInput("tExp","At (month)"))
        )
      }
      
    })
  })
  
  output$paramInput <- renderUI({
    if (input$spFun == "Hwang-Shih-DeCani") {
      list(
        textInput("hsd", "Gamma Parameter", placeholder = "Enter value (e.g., -3 or -4)"),
        h6("* Usually choose -3 or -4. When choosing -4, similar to Lan-DeMets O'BrienFleming")
      )
    } else if (input$spFun == 'Haybittle-Peto') {
      textInput("hp", "Interim Type I error (one-sided)", placeholder = "Choosing from Range: (0, 0.01)")
    } else if (input$spFun == 'Bespoke') {
      num_looks <- as.numeric(input$look)
      inputs <- lapply(1:num_looks, function(i) {
        if (i == num_looks) {
          textInput(paste("IA", i, sep = ""), "FA", value = as.numeric(input$t1e))
        } else {
          textInput(paste("IA", i, sep = ""), sprintf("IA%d", i), value = "")
        }
      })
      list(
        h5("Cumulative type I error (one-sided) up to each analysis"),
        do.call(tagList, inputs)
      )
     
    } else if (input$spFun %in% c("Lan-DeMets O'BrienFleming", "Lan-DeMets Pocock")) {
      return(NULL)  # No input needed for these stopping functions
    }
  })
  
  

# for accrual tab
  
  output$inACC <- renderUI({
    inputs <- list()
    
    if ("A" %in% input$accFun) {
      inputs[['wt']] <- list(
        h6("Cumulative recruitment distribution follows the form of (t/A)^w. When w=1, it's uniform enrollment."),
        textInput("wt", "Weight w"),
        textInput('dur',
                  'Enrollment duration A (month)')
      )
        
    }
    
    if ("B" %in% input$accFun) {
      inputs[['cuAcc']] <- list(
        div(
          h6("The enrollment distribution should be a proper cumulative distributed function (CDF) with a domain at (0, Inf)."),
          div(
            h6("Example 1:", style = "font-weight: bold;"),
            "Uniform enrollment in 24 months has CDF: t/24 * as.numeric(t<=24) + as.numeric(t>24);",
            style = "margin-bottom: 10px;"
          ),
          div(
            h6("Example 2:", style = "font-weight: bold;"),
            "Two-piece enrollment for a total of 300 patients with 10 patients per month for the first 6 months and 20 patients per month for the remaining 12 months.",
            "It has the CDF of (10*t/300)*as.numeric(t<=6) + ((60 + 20*(t-6))/300)*as.numeric(t>6 & t<=18) + as.numeric(t>18);",
            style = "margin-bottom: 10px;"
          )
        ),
        textInput("cuAcc","Customized recruitment distribution")
      )
        
    }
    
    do.call(tagList, inputs)
  })
  
  output$stOut <- renderUI({
    req(input$look)
    ui_elements <- lapply(1:as.numeric(input$look), function(i) {
      choice_id <- paste("stChoice", i, sep = "")
      
      # Radio buttons for selecting test type
      rb <- radioButtons(choice_id, 
                         label = sprintf("Analysis %d", i),
                         choices = c("logrank" = "logrank", "maximum of multiple tests(maxcombo)" = "maxcombo"))
      
      # Conditional UI elements for 'maxcombo'
      max_combo_ui <- conditionalPanel(
        condition = sprintf("input['%s'] == 'maxcombo'", choice_id),
        checkboxInput(paste("logrank", i), "Log rank"),
        fluidRow(
          column(2,checkboxInput(paste("fh1", i), "Fleming-Harrington test")),
          column(3,textInput(paste("rho1", i), "rho", "")),
          column(3,textInput(paste("gamma1", i), "gamma", "")),
          column(3,textInput(paste("stau1", i), "S(tau)", value='0'))
        ),
        fluidRow(
          column(2,checkboxInput(paste("fh2", i), "Fleming-Harrington test")),
          column(3,textInput(paste("rho2", i), "", "")),
          column(3,textInput(paste("gamma2", i), "", "")),
          column(3,textInput(paste("stau2", i), "", value='0'))
        ),
        fluidRow(
          column(2,checkboxInput(paste("fh3", i), "Fleming-Harrington test")),
          column(3,textInput(paste("rho3", i), "", "")),
          column(3,textInput(paste("gamma3", i), "", "")),
          column(3,textInput(paste("stau3", i), "", value='0'))
        ),
        fluidRow(
          column(2,checkboxInput(paste("custB", i), "Customized")),
          column(5,textInput(paste("custT", i), "Customized Function", "")),
          column(5,h6("* For customized function, please use 'a**b' to represent 'a^b' for the calculation of a to the power of b"))
        ),
        h6("*S(tau) is the truncation of survival rate at time tau. When S(tau) = 0, there is no truncation. When S(tau) = 0.5, the Fleming-Harrington test is truncated at median, i.e., for any survival time T > tau (median), the weight is calculated for survival rate of 0.5. This truncation method limits the weight up to a threshold S(tau).")
        
      )
      
      tagList(rb, max_combo_ui)
    })
    
    
    
    do.call(tagList, ui_elements)
  })
  
  wss <- reactive({
    req(input$look)
    ws_res <- list()
    for (i in 1:as.numeric(input$look)){
      if(input[[paste("stChoice", i, sep = "")]]=="logrank"){
        ws_cur <- list(function(s){1})
      }
      else {
        index <- 1
        ws_cur <- list()
        fw <- function(s, rho, gamma, s.tau){
          s.til = apply(cbind(s, s.tau), MARGIN=1,FUN=max)
          w = s.til^rho*(1-s.til)^gamma
          return(w)
        }
        
        if(input[[paste("logrank", i)]]){
          ws_cur[[index]] <- function(s){1}
          index <- index+1
          #showError(index)
          
        } 
        if (input[[paste("fh1", i)]]){
          rho <- as.numeric(input[[paste("rho1",i)]])
          gamma <- as.numeric(input[[paste("gamma1",i)]])
          stau <- as.numeric(input[[paste("stau1",i)]])
          ws_cur[[index]] <- function(s){fw(s,gamma=gamma,rho=rho,s.tau=stau)}
          index <- index+1
          #showError(ws_cur)
          #showError(index)
        } 
        if (input[[paste("fh2", i)]]){
          rho <- as.numeric(input[[paste("rho2",i)]])
          gamma <- as.numeric(input[[paste("gamma2",i)]])
          stau <- as.numeric(input[[paste("stau2",i)]])
          ws_cur[[index]] <- function(s){fw(s,rho=rho,gamma=gamma,s.tau=stau)}
          index <- index+1
          #showError(ws_cur)
        } 
        if (input[[paste("fh3", i)]]){
          rho <- as.numeric(input[[paste("rho3",i)]])
          gamma <- as.numeric(input[[paste("gamma3",i)]])
          stau <- as.numeric(input[[paste("stau3",i)]])
          ws_cur[[index]] <- function(s){fw(s,rho=rho,gamma=gamma,s.tau=stau)}
          index <- index+1
          #showError(ws_cur)
        } 
        if (input[[paste("custB", i)]]){
          ftext <- paste("function(s){",input[[paste("custT",i)]],"}")
          ws_cur[[index]] <- eval(parse(text=ftext))
          #showError(ws_cur)
        }
      }
      ws_res[[i]] <- ws_cur
      #showError(ws_res)
    }
    return(ws_res)
  })
  
  output$fDef <- renderText({
    if (!is.null(wss) && exists("wss")) {
      func <- wss()  
      if (!is.null(func)) {
        capture.output(func) 
      } else {
        "Function is being initialized..."  
      }
    } else {
      "Awaiting input..."  
    }
  })

  
# for wlr.power()
  
  F.entry <- reactive({
    req(input$accFun)  # Ensure these inputs are available
    dur <- as.numeric(input$dur)
    if (input$accFun == "A" && !is.null(input$wt)) {
      req(input$wt, input$dur)
      wt <- as.numeric(input$wt)  
      
      funcText <- sprintf("function(t) { (( t / %s) ^ %s) * as.numeric(t <= %s) + as.numeric(t > %s)}",
                          dur,wt, dur, dur)
      return(eval(parse(text=funcText)))
    } else if (input$accFun == "B" ) {
      req(input$cuAcc)
      funcText <- paste("function(t){",input$cuAcc,"}")
      
      return(eval(parse(text=funcText)))
    }
  })
  
  h0 <- reactive({
    req(input$testFunControl)  # Ensure testFunControl is selected
    if(input$testFunControl == "F"){
      if(input$conMethod == "Median"){
        funcText <- paste("function(t){log(2)/",input$medC,"}")
        return(eval(parse(text=funcText)))
      } else {
        funcText <- paste("function(t){-log(",input$surC,")/",input$tCon,"}")
        return(eval(parse(text=funcText)))
      }
    }else if (input$testFunControl == "A") {
      req(input$custH0)  
      funcText <- paste("function(t) {", input$custH0, "}")
      return(eval(parse(text=funcText)))
    } else if (input$testFunControl == "B") {
      req(input$con)  
      hc_values <- sapply(1:input$con, function(i) {
        hc <- input[[paste0("hc_", i)]]
        tc <- as.numeric(input[[paste0("tc_", i)]])
        tc_next <- if (i < input$con) {
          as.numeric(input[[paste0("tc_", i + 1)]])
        } else {
          "+Inf"  # Assume open interval for the last segment
        }
        if (is.na(hc) || is.na(tc) || is.na(tc_next)) {
          return(NULL)  
        } else {
          if (i == input$con) {
            return(paste(hc, "*as.numeric(t >=", tc, ")"))
          } else {
            return(paste(hc, "*as.numeric(t >=", tc, " & t < ", tc_next, ")"))
          }
        }
      })
      
      funcText <- paste("function(t) {", paste(hc_values, collapse=" + "), "}")
      return(eval(parse(text = funcText)))
    } else if (input$testFunControl == "C"){
      shape = eval(str2lang(input$shape1_c))
      scale = eval(str2lang(input$scale1_c))
      wh <- function(t) {
        (shape / scale) * (t / scale)^(shape - 1)
      }
      return(wh)
    } else if (input$testFunControl == "D"){
      lam = eval(str2lang(input$lam2_c))
      p = eval(str2lang(input$p2_c))
      ss0 <- function(t){
        return(exp(-lam*t))
      }
      s <- function(t){
        p + (1-p)*ss0(t)
      }
      f = function(t){dmcr(t,p=p,alpha=lam)}
      h = function(t){f(t)/s(t)}
      return(h)
    } else if (input$testFunControl == "E"){
      p = eval(str2lang(input$p3_c))
      shape = eval(str2lang(input$shape3_c))
      scale = eval(str2lang(input$scale3_c))
      alpha = scale^(-shape)
      ss0 <- function(t){
        exp(-alpha*t^shape)
      }
      s <- function(t){
        p + (1-p)*ss0(t)
      }
      f = function(t){dmcr(t,p=p,alpha=alpha, gamma=shape)}
      h = function(t){f(t)/s(t)}
      return(h)
      
    }
    
  })
  
  h1 <- reactive({
    req(input$distEXP)
    if(input$distEXP == "A"){
      hratio = eval(str2lang(input$ratioExp))
      hh1 <- function(t) {
        h0_t <- h0()(t)
        h0_t*hratio
      }
      return(hh1)
    } else {
      req(input$testFunExperiment)  
      if(input$testFunExperiment == "F"){
        if(input$expMethod == "Median"){
          funcText <- paste("function(t){log(2)/",input$medE,"}")
          return(eval(parse(text=funcText)))
        } else {
          funcText <- paste("function(t){-log(",input$surE,")/",input$tExp,"}")
          return(eval(parse(text=funcText)))
        }
        
      } else if (input$testFunExperiment == "A") {
        req(input$custH1)  
        funcText <- paste("function(t) {", input$custH1, "}")
        return(eval(parse(text=funcText)))
      } else if (input$testFunExperiment == "B") {
        req(input$exp)  
        he_values <- sapply(1:input$exp, function(i) {
          he <- input[[paste0("he_", i)]]
          te <- as.numeric(input[[paste0("te_", i)]])
          te_next <- if (i < input$exp) {
            as.numeric(input[[paste0("te_", i + 1)]])
          } else {
            "+Inf"  # Assume open interval for the last segment
          }
          if (is.na(he) || is.na(te) || is.na(te_next)) {
            return(NULL)  
          } else {
            if (i == input$exp) {
              return(paste(he, "*as.numeric(t >=", te, ")"))
            } else {
              return(paste(he, "*as.numeric(t >=", te, " & t < ", te_next, ")"))
            }
          }
        })
        
        funcText <- paste("function(t) {", paste(he_values, collapse=" + "), "}")
        return(eval(parse(text = funcText)))
      } else if (input$testFunExperiment == "C"){
        shape = eval(str2lang(input$shape1_e))
        scale = eval(str2lang(input$scale1_e))
        wh <- function(t) {
          (shape / scale) * (t / scale)^(shape - 1)
        }
        return(wh)
      } else if (input$testFunExperiment == "D"){
        lam = eval(str2lang(input$lam2_e))
        p = eval(str2lang(input$p2_e))
        ss0 <- function(t){
          return(exp(-lam*t))
        }
        s <- function(t){
          p + (1-p)*ss0(t)
        }
        f = function(t){dmcr(t,p=p,alpha=lam)}
        h = function(t){f(t)/s(t)}
        return(h)
      } else if (input$testFunExperiment == "E"){
        p = eval(str2lang(input$p3_e))
        shape = eval(str2lang(input$shape3_e))
        scale = eval(str2lang(input$scale3_e))
        alpha = scale^(-shape)
        ss0 <- function(t){
          exp(-alpha*t^shape)
        }
        s <- function(t){
          p + (1-p)*ss0(t)
        }
        f = function(t){dmcr(t,p=p,alpha=alpha, gamma=shape)}
        h = function(t){f(t)/s(t)}
        return(h)
      }
    }
    
  })
  
  logHR <- reactive({
    
    logHr <- function(t) {
      # Evaluate h1 and h0 at each time point t
      h1_t <- h1()(t)
      h0_t <- h0()(t)
      
      # Calculate the log hazard ratio
      log(h1_t / h0_t)
    }
    return(logHr)
  })
  
  s0 <- reactive({
    req(input$testFunControl)  
    if (input$testFunControl == "F") {
      
      if(input$conMethod == "Median"){
        funcText <- paste("function(t) {exp(-log(2)/", input$medC, "*t)}")
        return(eval(parse(text=funcText)))
      } else {
        funcText <- paste("function(t){exp(log(",input$surC,")/",input$tCon,"*t)}")
        return(eval(parse(text=funcText)))
      }
    } else if (input$testFunControl == "A") {
      req(input$custS0)
      funcText <- paste("function(t) {", input$custS0, "}")
      return(eval(parse(text=funcText)))
    } else if (input$testFunControl == "B") {
      req(input$con) 
      prev <- list()
      s <- list()
      
      for (i in 1:input$con) {
        hc <- input[[paste0("hc_", i)]]
        tc <- as.numeric(input[[paste0("tc_", i)]])
        tc_next <- if (i == input$con) {
          "+Inf" 
        } else {
          as.numeric(input[[paste0("tc_", i + 1)]])
        }
        if (i != input$con){
          pre <- ifelse(i==1,0,prev[[i-1]])
          if (i == 1){
            s[[i]] <- paste("(",pre,"+",hc,"*t","-",hc,"*",tc,")", "*", "as.numeric(t >=", tc, "& t <=", tc_next, ")")
          } else{
            s[[i]] <- paste("(",pre,"+",hc,"*t","-",hc,"*",tc,")", "*", "as.numeric(t >", tc, "& t <=", tc_next, ")")
          }
          
          func <- eval(parse(text= paste("function(t){",s[[i]],"}")))
          prev[[i]] <- as.numeric(func(tc_next))
        }
        else {
          pre <- ifelse(i==1,0,prev[[i-1]])
          if (i == 1){
            s[[i]] <- paste("(",pre,"+",hc,"*t","-",hc,"*",tc,")", "*", "as.numeric(t >=", tc,")")
          } else {
            s[[i]] <- paste("(",pre,"+",hc,"*t","-",hc,"*",tc,")", "*", "as.numeric(t >", tc,")")
          }
          
        }
        
      }
      H0 <- paste(s,collapse = "+")
      
      funcText <- paste("function(t){exp(-(",H0,"))}")
      return(eval(parse(text=funcText)))
      
    } else if (input$testFunControl == "C"){
      shape = eval(str2lang(input$shape1_c))
      scale = eval(str2lang(input$scale1_c))
      w <- function(t) {
        exp(-(t / scale)^shape)
      }
      return(w)
    } else if (input$testFunControl == "D"){
      lam = eval(str2lang(input$lam2_c))
      p = eval(str2lang(input$p2_c))
      ss0 <- function(t){
        return(exp(-lam*t))
      }
      s <- function(t){
        p + (1-p)*ss0(t)
      }
      return(s)
    } else if (input$testFunControl == "E"){
      p = eval(str2lang(input$p3_c))
      shape = eval(str2lang(input$shape3_c))
      scale = eval(str2lang(input$scale3_c))
      alpha = scale^(-shape)
      ss0 <- function(t){
        exp(-alpha*t^shape)
      }
      s <- function(t){
        p + (1-p)*ss0(t)
      }
      return(s)
    }
  })
  
  s1 <- reactive({
    req(input$distEXP)
    if(input$distEXP == "A"){
      hratio = eval(str2lang(input$ratioExp))
      ss1 <- function(t) {
        s0_t <- s0()(t)
        s0_t^hratio
      }
      return(ss1)
    } else{
      req(input$testFunExperiment)  
      if (input$testFunExperiment == "F") {
        if(input$expMethod == "Median"){
          funcText <- paste("function(t) {exp(-log(2)/", input$medE, "*t)}")
          return(eval(parse(text=funcText)))
        } else {
          funcText <- paste("function(t){exp(log(",input$surE,")/",input$tExp,"*t)}")
          return(eval(parse(text=funcText)))
        }
        
      } else if (input$testFunExperiment == "A") {
        req(input$custS1)
        funcText <- paste("function(t) {", input$custS1, "}")
        return(eval(parse(text=funcText)))
      } else if (input$testFunExperiment == "B") {
        req(input$exp) 
        prev <- list()
        s <- list()
        
        for (i in 1:input$exp) {
          he <- input[[paste0("he_", i)]]
          te <- as.numeric(input[[paste0("te_", i)]])
          te_next <- if (i == input$exp) {
            "+Inf" 
          } else {
            as.numeric(input[[paste0("te_", i + 1)]])
          }
          if (i != input$exp){
            pre <- ifelse(i==1,0,prev[[i-1]])
            if (i == 1){
              s[[i]] <- paste("(",pre,"+",he,"*t","-",he,"*",te,")", "*", "as.numeric(t >=", te, "& t <=", te_next, ")")
            } else {
              s[[i]] <- paste("(",pre,"+",he,"*t","-",he,"*",te,")", "*", "as.numeric(t >", te, "& t <=", te_next, ")")
            }
            
            func <- eval(parse(text= paste("function(t){",s[[i]],"}")))
            prev[[i]] <- as.numeric(func(te_next))
          }
          else {
            pre <- ifelse(i==1,0,prev[[i-1]])
            if (i == 1){
              s[[i]] <- paste("(",pre,"+",he,"*t","-",he,"*",te,")", "*", "as.numeric(t >=", te,")")
            } else{
              s[[i]] <- paste("(",pre,"+",he,"*t","-",he,"*",te,")", "*", "as.numeric(t >", te,")")
            }
            
          }
          
        }
        H1 <- paste(s,collapse = "+")
        
        funcText <- paste("function(t){exp(-(",H1,"))}")
        return(eval(parse(text=funcText)))
      } else if (input$testFunExperiment == "C"){
        shape = eval(str2lang(input$shape1_e))
        scale = eval(str2lang(input$scale1_e))
        w <- function(t) {
          exp(-(t / scale)^shape)
        }
        return(w)
      } else if (input$testFunExperiment == "D"){
        lam = eval(str2lang(input$lam2_e))
        p = eval(str2lang(input$p2_e))
        ss0 <- function(t){
          return(exp(-lam*t))
        }
        s <- function(t){
          p + (1-p)*ss0(t)
        }
        return(s)
      } else if (input$testFunExperiment == "E"){
        p = eval(str2lang(input$p3_e))
        shape = eval(str2lang(input$shape3_e))
        scale = eval(str2lang(input$scale3_e))
        alpha = scale^(-shape)
        ss0 <- function(t){
          exp(-alpha*t^shape)
        }
        s <- function(t){
          p + (1-p)*ss0(t)
        }
        return(s)
      }
    }
    
  })
  
  output$functionDef <- renderText({
    # Check if F.entry is defined and not null
    if (!is.null(F.entryS) && exists("F.entryS")) {
      func <- F.entryS()  # Call the reactive expression to get the current function
      if (!is.null(func)) {
        capture.output(func) 
      } else {
        "Function is being initialized..."  
      }
    } else {
      "Awaiting input..."  
    }
  })
  
  output$funcDef0 <- renderText({
    if (!is.null(s0) && exists("s0")) {
      func <- s0()  
      if (!is.null(func)) {
        capture.output(func) 
      } else {
        "Function is being initialized..."  
      }
    } else {
      "Awaiting input..."  
    }
  })

  
  observeEvent(input$plotBtn, {
    output$hz <- renderPlot({
      
      plot_S(
        S = list(h0(), h1()),
        Tmax = 50,
        leg= list(x = 30, y = 1, txt = c("Control Arm", "Experimental Arm")),
        param = list(xlab = "Time Since First Subject Randomized (mo)", ylab = "Hazard", main
                     = "Hazard Curve Per Study Design")
      )
      
    })
  })
  
  
  
  
  showPlot0 <- reactiveVal(FALSE)
  
  observeEvent(input$plotBtn0, {
    showPlot0(TRUE)  # Set the reactive value to TRUE when button is clicked
  })
  
  output$sur <- renderPlot({
    req(showPlot0())
    plot_S(
      S = list(s0(), s1()),
      Tmax = 50
    )
    
  })
  
  output$showPlot0 <- reactive({ showPlot0() })
  outputOptions(output, "showPlot0", suspendWhenHidden = FALSE)
  
  data_ac <- reactive({
    req(F.entry())  
    x_values <- 1:50
    y_values <- sapply(x_values, F.entry())
    data.frame(x = x_values, y = y_values)
  })
  
  showPlot2 <- reactiveVal(FALSE)
  
  observeEvent(input$plotBtn2, {
    showPlot2(TRUE)  # Set the reactive value to TRUE when button is clicked
  })
  
  output$accCurve <- renderPlot({
    req(showPlot2(),data_ac())
    ggplot(data_ac(), aes(x = x, y = y)) +
      geom_line(lwd=3,color = 'seagreen3') +
      labs(title = "Recruitment Distribution Over Time",
           x = "Time (month)", y = "Enrollment") +
      theme_minimal()
    
    
  })
  
  output$showPlot2 <- reactive({ showPlot2() })
  outputOptions(output, "showPlot2", suspendWhenHidden = FALSE)
  
  data_dpc <- reactive({
    req(G0())  
    x_values <- 1:50
    y_values <- sapply(x_values, G0())
    data.frame(x = x_values, y = y_values)
  })
  
  data_dpe <- reactive({
    req(G1())  
    x_values <- 1:50
    y_values <- sapply(x_values, G1())
    data.frame(x = x_values, y = y_values)
  })
  
  
  showPlot3 <- reactiveVal(FALSE)
  
  observeEvent(input$plotBtn3, {
    showPlot3(TRUE)  # Set the reactive value to TRUE when button is clicked
  })
  
  output$dpCurve1 <- renderPlot({
    req(showPlot3(),data_dpc())
    ggplot(data_dpc(), aes(x = x, y = y)) +
      geom_line(lwd=3,color = 'seagreen3') +
      labs(title = "Dropoff Distribution over Time for Control",
           x = "Time (month)", y = "Dropoff") +
      theme_minimal()
  })
  
  output$dpCurve2 <- renderPlot({
    req(showPlot3(),data_dpe())
    ggplot(data_dpe(), aes(x = x, y = y)) +
      geom_line(lwd=3,color='blue3') +
      labs(title = "Dropoff Distribution over Time for Experimental",
           x = "Time (month)", y = "Dropoff") +
      theme_minimal()
  })
  output$showPlot3 <- reactive({ showPlot3() })
  outputOptions(output, "showPlot3", suspendWhenHidden = FALSE)
  
  G0 <- reactive({
    req(input$conDP)
    co <- log(1-as.numeric(input$conDP)/12)
    funcText <- paste("function(t){1-exp(", co,"*t)}")
    return(eval(parse(text=funcText)))
  })
  
  G1 <- reactive({
    req(input$expDP) 
    ex <- log(1-as.numeric(input$expDP)/12)
    funcText <- paste("function(t){1-exp(", ex,"*t)}")
    return(eval(parse(text=funcText)))
  })
  
  G <- reactive({
    req(input$conDP, input$expDP, input$ratio)
    r = as.numeric(input$ratio)
    Gt = function(t){
      g1_t <- G1()(t)
      g0_t <- G0()(t)
      (1-r)*g0_t+r*g1_t
    }
    return(Gt)
  })
  
  
  output$dcoInput <- renderUI({
    req(input$look) 
    num_looks <- as.numeric(input$look)
    input_list <- lapply(1:num_looks, function(i) {
      numericInput(paste("dco_", i), label = sprintf("DCO for Analysis %d", i),value=20)
    })
    
    do.call(tagList, input_list)
  })

  
  p_res <- eventReactive(input$prun, {
    req(input$look, input$ratio, input$n,input$spFun)  # Ensure all inputs are ready
    dco <- c()
    for (i in 1:as.numeric(input$look)) {
      dco_input <- input[[paste("dco_", i)]]
      if (is.null(dco_input) || dco_input == "") {
        showError("All DCO values must be numeric and not empty")
        return(NULL)  # Stops the execution and won't proceed to calculation
      }
      if (is.na(as.numeric(dco_input))) {
        showError("Invalid numeric value for DCO")
        return(NULL)
      }
      dco <- c(dco, as.numeric(dco_input))
    }
    
    overall.alpha <- as.numeric(input$t1e)
    param <- NULL
    p1 <- NULL
    cum <- NULL
    
    if(input$spFun == "Hwang-Shih-DeCani"){
      sf <- "HSD"
      param <- as.numeric(input$hsd)
    } else if (input$spFun == "Lan-DeMets O'BrienFleming"){
      sf <- "LDOF"
    } else if (input$spFun == "Lan-DeMets Pocock"){
      sf <- "LDPK"
    } else if (input$spFun == "Haybittle-Peto"){
      sf <- "Haybittle-Peto"
      p1 <- as.numeric(input$hp)
    } else if (input$spFun == "Bespoke"){
      req(input$look)
      sf <- "Bespoke"
      cum <- c()
      for (i in 1:as.numeric(input$look)){
        curr <- input[[paste("IA", i, sep = "")]]
        cum <- c(cum, as.numeric(curr))
      }
    }
    
    
    ws <- list()
    for (i in 1:as.numeric(input$look)) {
      ws[[i]] <- list(function(s) { 1 } )
    }
    
    # Assuming `wlr.power.maxcombo` is a predefined function in your environment
    power <- wlr.power.maxcombo(
      r = as.numeric(input$ratio),
      n = as.numeric(input$n),
      DCO = dco,
      overall.alpha = overall.alpha,
      sf = sf,
      param = param,
      p1 = p1,
      cum.alpha = cum,
      h0 = h0(),  
      S0 = s0(),
      h1 = h1(),
      S1 = s1(),
      f.ws = wss(),
      Lambda = F.entry(),  
      G0 = G0(),
      G1 = G1(),
      mu.method = input$mu,
      cov.method = input$cov
    )
    round(data.frame(power$design),6)
  })
  
  output$powerResults <- renderDataTable({
    req(p_res())
    datatable(
      p_res(),
      extensions = 'Buttons',  # Enable buttons extension
      options = list(
        dom = 'Bfrtip',  # Positioning of table elements: buttons, filtering input, table info, pagination
        buttons = list('colvis'),  # Add a column visibility button
        scrollX = TRUE,
        pageLength = 10,
        autoWidth = TRUE
      )
    )
    
  })
  
  output$powerPlot1 <- renderPlotly({
    req(p_res())  # Ensure data is ready
    p=ggplot(p_res(), aes(x = DCO, y = cum.power)) +  # make sure column names are all lowercase or as named in the dataframe
      geom_line(lty=2,lwd=2, color = 'seagreen3') +
      geom_point(size=8, color = 'seagreen3')+
      labs(x = "Time since 1st subject randomized at each analysis", y = "Power") +
      theme_minimal()
    ggplotly(p)
  })
  
  output$powerPlot2 <- renderPlotly({
    req(p_res())  # Ensure data is ready
    p = ggplot(p_res(), aes(x = targetEvents, y = cum.power)) +  # Corrected to 'targetEvents' if that's the actual column name
      geom_line(lty=2,lwd=2, color = 'seagreen3') +
      geom_point(size=8, color = 'seagreen3')+
      labs(x = "Number of events at each analyis", y = "Power") +
      theme_minimal()
    ggplotly(p)
  })
  
  observeEvent(input$eve, {
    
    output$eplot <- renderPlot({
      req(input$n,input$tmax,input$ratio)
      e <<- plot_events(
        n=as.numeric(input$n),
        Tmax=as.numeric(input$tmax),
        r=as.numeric(input$ratio),
        h0=h0(),S0=s0(),
        h1=h1(),S1=s1(),
        Lambda=F.entry(),
        G0=G0(),G1=G1()
      )
    })
  })
  
  observeEvent(input$fdco,{
    req(input$targetEvents, input$n, input$ratio)
    output$dco <- renderText({
      paste("DCO:", fDCO(events = input$targetEvents, 
                                   r = as.numeric(input$ratio), 
                                   h0=h0(),S0=s0(),
                                   h1=h1(),S1=s1(),
                                   Lambda=F.entry(), 
                                   n = as.numeric(input$n)))
    })
  })
  
  observeEvent(input$eve,{
    req(input$n,input$tmax,input$ratio)
    output$eTable <- renderDataTable({
      data.frame(round(e,4))
    })
  })
  
  
  observeEvent(input$bButton, {
    
    output$bPlot <- renderPlotly({
      req(input$bu,input$cov,p_res())
      if(input$bu == "Z-scale"){
        p=ggplot(p_res(), aes(x = DCO, y = bd)) +  
          geom_line(lty=2,lwd=2,color='seagreen3') +
          geom_point(size=8,color='seagreen3')+
          labs(x = "Time since 1st subject randomized at each analysis", y = "Boundary (Z-scale)") +
          theme_minimal()
        ggplotly(p)
        
      } else if(input$bu == "p-value(one-sided)"){
        p=ggplot(p_res(), aes(x = DCO, y = p)) +  
          geom_line(lty=2,lwd=2,color='seagreen3') +
          geom_point(size=8,color='seagreen3')+
          labs(x = "Time since 1st subject randomized at each analysis", y = "Boundary (p-value, one-sided)") +
          theme_minimal()
        ggplotly(p)
        
      } else if(input$bu == "Hazard ratio"){
        if(input$cov == "H0"){
          p=ggplot(p_res(), aes(x = DCO, y = CV.HR.H0)) +  
            geom_line(lty=2,lwd=2, color='seagreen3') +
            geom_point(size=8, color='seagreen3')+
            labs(x = "Time since 1st subject randomized at each analysis", y = "Boundary (Hazard ratio)") +
            theme_minimal()
          ggplotly(p)
        } else {
          p=ggplot(p_res(), aes(x = DCO, y = CV.HR.H1)) +  
            geom_line(lty=2,lwd=2,color='seagreen3') +
            geom_point(size=8,color='seagreen3')+
            labs(x = "Time since 1st subject randomized at each analysis", y = "Boundary (Hazard ratio)") +
            theme_minimal()
          ggplotly(p)
        }
      }
      
    })
  })
  
  output$aaPlot <- renderPlotly({
    req(p_res())
    p=ggplot(p_res(), aes(x = DCO, y = cum.alpha)) +  
      geom_line(lty=2,lwd=2,color='seagreen3') +
      geom_point(size=8,color='seagreen3')+
      labs(x = "Time since 1st subject randomized at each analysis", y = "Cumulative Alpha Spending") +
      theme_minimal()
    ggplotly(p)
  })
  
  observeEvent(input$ahrButton, {
    
    output$ahrPlot <- renderPlot({
      req(input$ahr,input$tmaxA)
      plot_AHR(
        n = as.numeric(input$n),
        Tmax = as.numeric(input$tmaxA),
        r = as.numeric(input$ratio),
        h0 = h0(),  
        S0 = s0(),
        h1 = h1(),
        S1 = s1(),
        Lambda = F.entry(),
        method = input$ahr,
        G = G()
      )
      
      
    })
  }) 
  
  
  observeEvent(input$alpha, {
    
    output$as <- renderPlot({
      req(input$spFun)
      oa <- as.numeric(input$t1e)
      param <- NULL
      p1 <- NULL
      cum <- NULL
      
      if(input$spFun == "Hwang-Shih-DeCani"){
        req(input$hsd)
        sf <- "HSD"
        param <- as.numeric(input$hsd)
      } else if (input$spFun == "Lan-DeMets O'BrienFleming"){
        sf <- "LDOF"
      } else if (input$spFun == "Lan-DeMets Pocock"){
        sf <- "LDPK"
      } else if (input$spFun == "Haybittle-Peto"){
        req(input$hp)
        sf <- "Haybittle-Peto"
        p1 <- as.numeric(input$hp)
      } else if (input$spFun == "Bespoke"){
        req(input$look)
        sf <- "Bespoke"
        cum <- c()
        for (i in 1:as.numeric(input$look)){
          curr <- input[[paste("IA", i, sep = "")]]
          cum <- c(cum, as.numeric(curr))
        }
      }
      
      t = seq(0.01, 1, by=0.01)
      a = rep(NA, length(t))
      for (i in 1:length(t)){
        a[i] = f.alpha(overall.alpha = oa, 
                       sf=sf, timing = c(t[i], 1), p1 = p1,
                       cum.alpha=cum,param=param)[1]
      }
      plot(t, a, type="l", xlab="Information Time", ylab = "IA Alpha Spending", 
           col = "seagreen3", lwd = 5)
    })
  })
  
  

### SIMULATION
  output$paramInputS <- renderUI({
    if (input$spFunS == "Hwang-Shih-DeCani") {
      list(
        textInput("hsdS", "Gamma Parameter", placeholder = "Enter value (e.g., -3 or -4)"),
        h6("* Usually choose -3 or -4. When choosing -4, similar to Lan-DeMets O'BrienFleming")
      )
    } else if (input$spFunS == 'Haybittle-Peto') {
      textInput("hpS", "Interim Type I error (one-sided)", placeholder = "Choosing from Range: (0, 0.01)")
    } else if (input$spFunS == 'Bespoke') {
      num_looks <- as.numeric(input$lookS)
      inputs <- lapply(1:num_looks, function(i) {
        if (i == num_looks) {
          textInput(paste("IAS", i, sep = ""), "FA", value = as.numeric(input$t1eS))
        } else {
          textInput(paste("IAS", i, sep = ""), sprintf("IA%d", i), value = "")
        }
      })
      list(
        h5("Cumulative type I error (one-sided) up to each analysis"),
        do.call(tagList, inputs)
      )
      
    } else if (input$spFunS %in% c("Lan-DeMets O'BrienFleming", "Lan-DeMets Pocock")) {
      return(NULL)  # No input needed for these stopping functions
    }
  })
  
  output$inACCS <- renderUI({
    inputs <- list()
    
    if ("A" %in% input$accFunS) {
      inputs[['wt']] <- list(
        h6("Cumulative recruitment distribution follows the form of (t/A)^w. When w=1, it's uniform enrollment."),
        textInput("wtS", "Weight w"),
        textInput('durS',
                  'Enrollment duration A (month)')
      )
      
    }
    
    if ("B" %in% input$accFunS) {
      inputs[['cuAcc']] <- list(
        div(
          h6("The enrollment distribution should be a proper cumulative distributed function (CDF) with a domain at (0, Inf)."),
          div(
            h6("Example 1:", style = "font-weight: bold;"),
            "Uniform enrollment in 24 months has CDF: t/24 * as.numeric(t<=24) + as.numeric(t>24);",
            style = "margin-bottom: 10px;"
          ),
          div(
            h6("Example 2:", style = "font-weight: bold;"),
            "Two-piece enrollment for a total of 300 patients with 10 patients per month for the first 6 months and 20 patients per month for the remaining 12 months.",
            "It has the CDF of (10*t/300)*as.numeric(t<=6) + ((60 + 20*(t-6))/300)*as.numeric(t>6 & t<=18) + as.numeric(t>18);",
            style = "margin-bottom: 10px;"
          )
        ),
        textInput("cuAccS","Customized recruitment distribution")
      )
      
    }
    
    do.call(tagList, inputs)
  })
  
  
  output$stOutS <- renderUI({
    req(input$lookS)
    ui_elements <- lapply(1:as.numeric(input$lookS), function(i) {
      choice_id <- paste("stChoiceS", i, sep = "")
      
      # Radio buttons for selecting test type
      rb <- radioButtons(choice_id, 
                         label = sprintf("Analysis %d", i),
                         choices = c("logrank" = "logrank", "maximum of multiple tests(maxcombo)" = "maxcombo"))
      
      # Conditional UI elements for 'maxcombo'
      max_combo_ui <- conditionalPanel(
        condition = sprintf("input['%s'] == 'maxcombo'", choice_id),
        checkboxInput(paste("logrankS", i), "Log rank"),
        fluidRow(
          column(2,checkboxInput(paste("fh1S", i), "Fleming-Harrington test")),
          column(3,textInput(paste("rho1S", i), "rho", "")),
          column(3,textInput(paste("gamma1S", i), "gamma", "")),
          column(3,textInput(paste("stau1S", i), "S(tau)", value='0'))
        ),
        fluidRow(
          column(2,checkboxInput(paste("fh2S", i), "Fleming-Harrington test")),
          column(3,textInput(paste("rho2S", i), "", "")),
          column(3,textInput(paste("gamma2S", i), "", "")),
          column(3,textInput(paste("stau2S", i), "", value='0'))
        ),
        fluidRow(
          column(2,checkboxInput(paste("fh3S", i), "Fleming-Harrington test")),
          column(3,textInput(paste("rho3S", i), "", "")),
          column(3,textInput(paste("gamma3S", i), "", "")),
          column(3,textInput(paste("stau3S", i), "", value='0'))
        ),
        fluidRow(
          column(2,checkboxInput(paste("custBS", i), "Customized")),
          column(5,textInput(paste("custTS", i), "Customized Function", "")),
          column(5,h6("* For customized function, please use 'a**b' to represent 'a^b' for the calculation of a to the power of b"))
        ),
        h6("*S(tau) is the truncation of survival rate at time tau. When S(tau) = 0, there is no truncation. When S(tau) = 0.5, the Fleming-Harrington test is truncated at median, i.e., for any survival time T > tau (median), the weight is calculated for survival rate of 0.5. This truncation method limits the weight up to a threshold S(tau).")
        
      )
      
      tagList(rb, max_combo_ui)
    })
    
    do.call(tagList, ui_elements)
  })
  
  wss_s <- reactive({
    req(input$lookS)
    ws_ress <- list()
    for (i in 1:as.numeric(input$lookS)){
      if(input[[paste("stChoiceS", i, sep = "")]]=="logrank"){
        ws_curs <- list(function(s){1})
      }
      else {
        index <- 1
        ws_curs <- list()
        fw <- function(s, rho, gamma, s.tau){
          s.til = apply(cbind(s, s.tau), MARGIN=1,FUN=max)
          w = s.til^rho*(1-s.til)^gamma
          return(w)
        }
        
        if(input[[paste("logrankS", i)]]){
          ws_curs[[index]] <- function(s){1}
          index <- index+1
          #showError(index)
          
        } 
        if (input[[paste("fh1S", i)]]){
          rho <- as.numeric(input[[paste("rho1S",i)]])
          gamma <- as.numeric(input[[paste("gamma1S",i)]])
          stau <- as.numeric(input[[paste("stau1S",i)]])
          ws_curs[[index]] <- function(s){fw(s,gamma=gamma,rho=rho,s.tau=stau)}
          index <- index+1
          #showError(ws_cur)
          #showError(index)
        } 
        if (input[[paste("fh2S", i)]]){
          rho <- as.numeric(input[[paste("rho2S",i)]])
          gamma <- as.numeric(input[[paste("gamma2S",i)]])
          stau <- as.numeric(input[[paste("stau2S",i)]])
          ws_curs[[index]] <- function(s){fw(s,rho=rho,gamma=gamma,s.tau=stau)}
          index <- index+1
          #showError(ws_cur)
        } 
        if (input[[paste("fh3S", i)]]){
          rho <- as.numeric(input[[paste("rho3S",i)]])
          gamma <- as.numeric(input[[paste("gamma3S",i)]])
          stau <- as.numeric(input[[paste("stau3S",i)]])
          ws_curs[[index]] <- function(s){fw(s,rho=rho,gamma=gamma,s.tau=stau)}
          index <- index+1
          #showError(ws_cur)
        } 
        if (input[[paste("custBS", i)]]){
          ftext <- paste("function(s){",input[[paste("custTS",i)]],"}")
          ws_curs[[index]] <- eval(parse(text=ftext))
          #showError(ws_cur)
        }
      }
      ws_ress[[i]] <- ws_curs
      #showError(ws_res)
    }
    ws_result = list()
    ws_result[[1]] = ws_ress
    return(ws_result)
  })
  
  output$fDefS <- renderText({
    if (!is.null(wss_s) && exists("wss_s")) {
      func <- wss_s()  
      if (!is.null(func)) {
        capture.output(func) 
      } else {
        "Function is being initialized..."  
      }
    } else {
      "Awaiting input..."  
    }
  })
  
  output$eveInputS <- renderUI({
    req(input$lookS) 
    num_looks <- as.numeric(input$lookS)
    input_list <- lapply(1:num_looks, function(i) {
      numericInput(paste("eveS_", i), label = sprintf("Target events for Analysis %d", i),value=0)
    })
    
    do.call(tagList, input_list)
  })
  
  
  output$inFunC <- renderUI({
    inputs <- list()
    if (input$testFunC == "F"){
      inputs[['conMethodS']] <- radioButtons("conMethodS", "Please choose one method to input. ",choices = list("Median","Survival"))
    }
    
    if (input$testFunC == "A") {
      message <- HTML(paste(
        "Please provide the function formula for the control arm. ",
        sep = ""
      ))
      inputs[['A']] <- HTML(paste("<p>", message, "</p>"))
      inputs[['H1']] <- textInput("custS0S", "Survival function of control")
    }
    
    if (input$testFunC == "B") {
      message2 <- "Please select the number of pieces and input the corresponding piecewise rate for control arm"
      inputs[['A']] <- HTML(paste("<p>", message2, "</p>"))
      inputs[['piece']] <- selectInput("conS", "Number of pieces for control", choices = 1:5)
    }
    
    if (input$testFunC == "C") {
      message3 <- "Please input the corresponding weibull parameters for control arm"
      inputs[['A']] <- HTML(paste("<p>", message3, "</p>"))
      inputs[['shape1_cS']] <- textInput("shape1_cS", "Shape")
      inputs[['scale1_cS']] <- textInput("scale1_cS", "Scale")
      
    }
    if (input$testFunC == "D") {
      message4 <- "Please input the corresponding mixture cure rate of exponential for control arm"
      inputs[['A']] <- HTML(paste("<p>", message4, "</p>"))
      inputs[['p2_cS']] <- textInput("p2_cS", "Cure rate")
      inputs[['lam2_cS']] <- textInput("lam2_cS", "Hazard rate") 
    }
    if (input$testFunC == "E") {
      message5 <- "Please input the corresponding mixture cure rate of weibull for control arm"
      inputs[['A']] <- HTML(paste("<p>", message5, "</p>"))
      inputs[['p3_cS']] <- textInput("p3_cS", "Cure rate")
      inputs[['shape3_cS']] <- textInput("shape3_cS", "Shape")
      inputs[['scale3_cS']] <- textInput("scale3_cS", "Scale") 
      
    }
    
    do.call(tagList, inputs)
  })
  
  output$inFunE <- renderUI({
    inputs <- list()
    if (input$testFunE == "F"){
     inputs[['expMethodS']] <- radioButtons("expMethodS", "Please choose one method to input. ",choices = list("Median","Survival"))
    }
    if (input$testFunE == "A") {
      message <- "Please provide the function formula for the experimental arm."
      inputs[['A']] <- HTML(paste("<p>", message, "</p>"))
      inputs[['H1']] <- textInput("custS1S", "Survival function of experiment")
    }
    
    if (input$testFunE == "B") {
      message2 <- "Please select the number of pieces and input the corresponding piecewise rate for experimental arm"
      inputs[['A']] <- HTML(paste("<p>", message2, "</p>"))
      inputs[['piece']] <- selectInput("expS", "Number of pieces for experiment", choices = 1:5)
    }
    
    if (input$testFunE == "C") {
      message3 <- "Please input the corresponding weibull parameters for experimental arm"
      inputs[['A']] <- HTML(paste("<p>", message3, "</p>"))
      inputs[['shape1_eS']] <- textInput("shape1_eS", "Shape")
      inputs[['scale1_eS']] <- textInput("scale1_eS", "Scale")
      
    }
    if (input$testFunE == "D") {
      message4 <- "Please input the corresponding mixture cure rate of exponential for experimental arm"
      inputs[['A']] <- HTML(paste("<p>", message4, "</p>"))
      inputs[['p2_eS']] <- textInput("p2_eS", "Cure rate")
      inputs[['lam2_eS']] <- textInput("lam2_eS", "Hazard rate") 
    }
    if (input$testFunE == "E") {
      message5 <- "Please input the corresponding mixture cure rate of weibull for experimental arm"
      inputs[['A']] <- HTML(paste("<p>", message5, "</p>"))
      inputs[['p3_eS']] <- textInput("p3_eS", "Cure rate")
      inputs[['shape3_eS']] <- textInput("shape3_eS", "Shape")
      inputs[['scale3_eS']] <- textInput("scale3_eS", "Scale") 
      
    }
    
    do.call(tagList, inputs)
  })
  
  observeEvent(input$conS, {
    req(input$testFunC == "B")
    output$moreInputsC <- renderUI({
      numCon <- as.numeric(input$conS)
      controlInputs <- lapply(1:numCon, function(i) {
        fluidRow(
          column(6, textInput(paste("tcS", i, sep = "_"), "Starting Time (months)", value = if (i == 1) 0 else NA)),
          column(6, textInput(paste("hcS", i, sep = "_"), "Hazard Rate"))
        )
      })
      do.call(tagList, controlInputs)
    })
  })
  
  observeEvent(input$expS, {
    req(input$testFunE == "B")
    output$moreInputsE <- renderUI({
      numExp <- as.numeric(input$expS)
      experimentInputs <- lapply(1:numExp, function(j) {
        fluidRow(
          column(6, textInput(paste("teS", j, sep = "_"), "Starting Time (months)", value = if (j == 1) 0 else NA)),
          column(6, textInput(paste("heS", j, sep = "_"), "Hazard Rate"))
        )
      })
      do.call(tagList, experimentInputs)
    })
  })
  
  observeEvent(input$conMethodS, {
    req(input$testFunC == "F")
    output$medConOutS <- renderUI({
      if (input$conMethodS == "Median"){
        textInput("medCS","Median of Control")
      } else{
        fluidRow(
          column(6,
                 textInput("surCS","Survival")),
          column(6,
                 textInput("tConS","At (month)"))
        )
        
      }
      
    })
  })
  
  observeEvent(input$expMethodS, {
    req(input$testFunE == "F")
    output$medExpOutS <- renderUI({
      if (input$expMethodS == "Median"){
        textInput("medES","Median of Experimental")
      } else{
        fluidRow(
          column(6,
                 textInput("surES","Survival")),
          column(6,
                 textInput("tExpS","At (month)"))
        )
        
      }
      
    })
  })
  
  cuts0 <- reactive({
    req(input$conS)
    cut <- c()
    for (i in 1:as.numeric(input$conS)){
      cut <- c(cut,as.numeric(input[[paste("tcS", i, sep = "_")]]))
    }
    return(cut)
  })
  
  cuts1 <- reactive({
    req(input$expS)
    cut <- c()
    for (i in 1:as.numeric(input$expS)){
      cut <- c(cut,as.numeric(input[[paste("teS", i, sep = "_")]]))
    }
    return(cut)
  })
  
  F.entryS <- reactive({
    req(input$accFunS)  # Ensure these inputs are available
    dur <- as.numeric(input$durS)
    if (input$accFunS == "A" && !is.null(input$wtS)) {
      req(input$wtS, input$durS)
      wt <- as.numeric(input$wtS)  
      
      funcText <- sprintf("function(t) { (( t / %s) ^ %s) * as.numeric(t <= %s) + as.numeric(t > %s)}",
                          dur,wt, dur, dur)
      return(eval(parse(text=funcText)))
    } else if (input$accFunS == "B" ) {
      req(input$cuAccS)
      funcText <- paste("function(t){",input$cuAccS,"}")
      return(eval(parse(text=funcText)))
    }
  })
  
  h0S <- reactive({
    req(input$testFunC)  # Ensure testFunControl is selected
    
    if (input$testFunC == "B") {
      req(input$conS)  
      hc_values <- c()
      for (j in 1:as.numeric(input$conS)) {
        a <- eval(str2lang(input[[paste("hcS", j, sep = "_")]]))
        hc_values[j] <- as.numeric(a)
      }
      return(hc_values)
    }
  })
  
  h1S <- reactive({
    req(input$testFunE)  # Ensure testFunControl is selected
    
    if (input$testFunE == "B") {
      req(input$expS)  
      he_values <- c()
      for (i in 1:as.numeric(input$expS)) {
        a <- eval(str2lang(input[[paste("heS", i, sep = "_")]]))
        he_values[i] <- as.numeric(a)
        }
      return(he_values)
    }
  })
  
  s0S <- reactive({
    req(input$testFunC)  
    
    if (input$testFunC == "A") {
      req(input$custS0S)
      funcText <- paste("function(t) {", input$custS0S, "}")
      return(eval(parse(text=funcText)))
    } 
  })
  
  s1S <- reactive({
    req(input$testFunE)  
    
    if (input$testFunE == "A") {
      req(input$custS1S)
      funcText <- paste("function(t) {", input$custS1S, "}")
      return(eval(parse(text=funcText)))
    } 
  })
  output$funcDef <- renderText({
    if (!is.null(h1S) && exists("h1S")) {
      func <- h1S()  
      if (!is.null(func)) {
        capture.output(func) 
      } else {
        "Function is being initialized..."  
      }
    } else {
      "Awaiting input..."  
    }
  })
  
  p_ress <- eventReactive(input$prunS, {
    req(input$lookS, input$ratioS, input$nS,input$spFunS)  # Ensure all inputs are ready
    eve <- c()
    for (i in 1:as.numeric(input$lookS)) {
      eve_input <- input[[paste("eveS_", i)]]
      if (is.null(eve_input) || eve_input == "") {
        showError("All DCO values must be numeric and not empty")
        return(NULL)  # Stops the execution and won't proceed to calculation
      }
      if (is.na(as.numeric(eve_input))) {
        showError("Invalid numeric value for DCO")
        return(NULL)
      }
      eve <- c(eve, as.numeric(eve_input))
    }
    
    overall.alpha <- as.numeric(input$t1eS)
    param <- NULL
    p1 <- NULL
    cum <- NULL
    
    if(input$spFunS == "Hwang-Shih-DeCani"){
      sf <- "HSD"
      param <- as.numeric(input$hsdS)
    } else if (input$spFunS == "Lan-DeMets O'BrienFleming"){
      sf <- "LDOF"
    } else if (input$spFunS == "Lan-DeMets Pocock"){
      sf <- "LDPK"
    } else if (input$spFunS == "Haybittle-Peto"){
      sf <- "Haybittle-Peto"
      p1 <- as.numeric(input$hpS)
    } else if (input$spFunS == "Bespoke"){
      req(input$lookS)
      sf <- "Bespoke"
      cum <- c()
      for (i in 1:as.numeric(input$lookS)){
        curr <- input[[paste("IAS", i, sep = "")]]
        cum <- c(cum, as.numeric(curr))
      }
    }
    dist0 = NULL
    lam0 = NULL
    shape0 = NULL
    scale0 = NULL
    p10 = NULL
    S0 = NULL
    cuts0 = NULL
    dist1 = NULL
    lam1 = NULL
    shape1 = NULL
    scale1 = NULL
    p11 = NULL
    S1 = NULL
    cuts1 = NULL
    HR = NULL
    
    if(input$testFunC=="A"){
      dist0 = "customized"
      S0 = s0S()
      
    } else if(input$testFunC=="B"){
      dist0 = "piecewise exponential"
      lam0 = h0S()
      cuts0 = cuts0()
    } else if (input$testFunC=="C"){
      dist0 = "weibull"
      shape0 = as.numeric(eval(str2lang(input$shape1_cS)))
      scale0 = as.numeric(eval(str2lang(input$scale1_cS)))
    } else if(input$testFunC=="D"){
      dist0 = "mixture cure rate of exponential"
      p10 = as.numeric(eval(str2lang(input$p2_cS)))
      lam0 = as.numeric(eval(str2lang(input$lam2_cS)))
    } else if(input$testFunC=="E"){
      dist0 = "mixture cure rate of weibull"
      p10 = as.numeric(eval(str2lang(input$p3_cS)))
      scale0 = as.numeric(eval(str2lang(input$scale3_cS)))
      shape0 = as.numeric(eval(str2lang(input$shape3_cS)))
    } else if(input$testFunC=="F"){
      dist0 = "exponential"
      cuts0 = NULL
      if(input$conMethodS == 'Median'){
        lam0 = log(2)/as.numeric(eval(str2lang(input$medCS)))
      } else {
        SR = as.numeric(eval(str2lang(input$surCS)))
        ATO = as.numeric(eval(str2lang(input$tConS)))
        lam0 = -log(as.numeric(SR))/ATO
      }
      
    }
    
    if (input$distEXPS == 'A'){
      dist1 = "Proportional Hazards"
      HR = as.numeric(input$ratioExpS)
    } else{
      if(input$testFunE=="A"){
        dist1 = "customized"
        S1 = s1S()
      } else if(input$testFunE=="B"){
        dist1 = "piecewise exponential"
        lam1 = h1S()
        cuts1 = cuts1()
      } else if (input$testFunE=="C"){
        dist1 = "weibull"
        shape1 = as.numeric(eval(str2lang(input$shape1_eS)))
        scale1 = as.numeric(eval(str2lang(input$scale1_eS)))
      } else if(input$testFunE=="D"){
        dist1 = "mixture cure rate of exponential"
        p11 = as.numeric(eval(str2lang(input$p2_eS)))
        lam1 = as.numeric(eval(str2lang(input$lam2_eS)))
      } else if(input$testFunE=="E"){
        dist1 = "mixture cure rate of weibull"
        p11 = as.numeric(eval(str2lang(input$p3_eS)))
        scale1 = as.numeric(eval(str2lang(input$scale3_eS)))
        shape1 = as.numeric(eval(str2lang(input$shape3_eS)))
      } else if(input$testFunE=="F"){
        dist1 = "exponential"
        cuts1 = NULL
        if(input$expMethodS == 'Median'){
          lam1 = log(2)/as.numeric(eval(str2lang(input$medES)))
        } else {
          SR = as.numeric(eval(str2lang(input$surES)))
          ATO = as.numeric(eval(str2lang(input$tExpS)))
          lam1 = -log(as.numeric(SR))/ATO
        }
      }}
      parall <- if_else(input$parS == 'T', TRUE, FALSE)
      
      pw_exp0 = (input$testFunC == 'F' | input$testFunC == 'B')
      pw_exp1 = (input$testFunE == 'F' | input$testFunE == 'B' | input$distEXPS == 'A')
      
      if (pw_exp0 & pw_exp1){
        if(input$distEXPS == 'A'){
          lam1 = lam0 * HR
          cuts1 = cuts0
        }
        powerS <- simulation.nphDesign.pwexp(
          nSim = as.numeric(input$nSim),
          r = as.numeric(input$ratioS),
          N = as.numeric(input$nS),
          A = as.numeric(input$durS),
          w = as.numeric(input$wtS),
          targetEvents = eve,
          overall.alpha = overall.alpha,
          sf = sf,
          param = param,
          p1 = p1,
          cum.alpha = cum, 
          lam0 = lam0, cuts0 = cuts0,
          lam1 = lam1, cuts1 = cuts1,
          fws.options = wss_s(),
          Lambda = F.entryS(),  
          drop0 = as.numeric(input$conDPS)/12,
          drop1 = as.numeric(input$expDPS)/12,
          logrank = input$lr,
          parallel = parall,
          n.cores = as.numeric(input$coreS),
          seed = as.numeric(input$seed)
        )
          
      } else {
        powerS <- simulation.nphDesign(
          nSim = as.numeric(input$nSim),
          r = as.numeric(input$ratioS),
          n = as.numeric(input$nS),
          targetEvents = eve,
          overall.alpha = overall.alpha,
          sf = sf,
          param = param,
          p1 = p1,
          cum.alpha = cum, dist0 = dist0,
          lam0 = lam0, shape0 = shape0, scale0 = scale0,
          p10 = p10, S0 = S0, cuts0 = cuts0,
          dist1 = dist1, HR = HR,
          lam1 = lam1, shape1 = shape1, scale1 = scale1,
          p11 = p11, S1 = S1, cuts1 = cuts1,
          fws.options = wss_s(),
          Lambda = F.entryS(),  
          drop0 = as.numeric(input$conDPS)/12,
          drop1 = as.numeric(input$expDPS)/12,
          logrank = input$lr,
          parallel = parall,
          n.cores = as.numeric(input$coreS),
          seed = as.numeric(input$seed)
        )
      }
      
      
      
      if (input$lr == 'N') {
        data.frame(
          Analysis = 1:as.numeric(input$lookS),
          TargetEvents = eve,
          Cumulative.Power = c(powerS$cum.pow),
          row.names = NULL
        )
      } else if (input$lr == 'Y') {
        data.frame(
          Analysis = 1:as.numeric(input$lookS),
          TargetEvents = eve,
          Cumulative.Power = c(powerS$cum.pow),
          LogRank.Cumulative.Power = c(powerS$lr.cum.pow),
          row.names = NULL
        )
      }
    
    
  })
  
  output$power <- renderDataTable({
    req(p_ress())
    datatable(
      p_ress(),
      extensions = 'Buttons',  
      options = list(
        dom = 'Bfrtip',  
        buttons = list('colvis'), 
        scrollX = TRUE,
        pageLength = 10,
        rownames = FALSE
      )
    )
    
  })
  
  output$sim1 <- renderPlotly({
    req(p_ress())
    if(input$lr == 'N'){
      p = ggplot(p_ress(),aes(x=TargetEvents, y=Cumulative.Power)) +
        geom_line(lty=2,lwd=3,color='seagreen3')+
        geom_point(size=8,color='seagreen3')+
        labs(x='Target Events')+
        theme_minimal()
      ggplotly(p)
    } else{
      long_data <- p_ress() %>%
        pivot_longer(
          cols = c(Cumulative.Power, LogRank.Cumulative.Power),
          names_to = "PowerType",
          values_to = "Power"
        )
      
      p = ggplot(long_data, aes(x = TargetEvents, y = Power, color = PowerType, shape = PowerType, linetype = PowerType)) +
        geom_line(lwd = 3) +  
        geom_point(size = 8) +  
        labs(x = 'Target Events', y = 'Power') +
        theme_minimal() +
        scale_color_manual(values = c('seagreen3', 'blue3')) +
        scale_shape_manual(values = c(19, 17)) +
        scale_linetype_manual(values = c("dotdash", "dashed"))
      ggplotly(p)
    }
    
  })
  
}



