#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(DT)
library(nphRshiny)
library(tidyverse)
 
shinyUI(navbarPage(title = 'Demo',
                   tabPanel('Design',
                            h3('Here is the test example for design part'),
                        
                            # Create tabs 
                            tabsetPanel(
                              tabPanel("Type I error",
                                       fluidRow(
                                  
                                         column(12,
                                                selectInput("look", 
                                                            "Number of analyses that require alpha spending",
                                                            choices = c(1:5)),
                                                
                                                textInput("t1e",
                                                          "Overall Type I error (one-sided), eg. 0.025 for single primary"),
                                                fluidRow(
                                                  column(6,
                                                         selectInput("spFun",
                                                                     "Alpha Spending Function",
                                                                     choices=c("Lan-DeMets O'BrienFleming","Lan-DeMets Pocock","Hwang-Shih-DeCani",'Haybittle-Peto',"Bespoke"))
                                                  ),
                                                  column(6,uiOutput("paramInput"))),
                                                
                                                actionButton("alpha", "Display IA Alpha Spending Curve for Exploration"),
                                                plotOutput("as")
                                                )) 
                                       ),
                              tabPanel("Distributions", 
                                       fluidRow(
                                         column(6,
                                                radioButtons("testFunControl", "Choose a function for Control",
                                                             choices = list("Piecewise exponential" = "B","Weibull" = "C","Mixture cure rate of exponential"="D",
                                                                            "Mixture cure rate of weibull"='E',"Customized distribution" = "A")),
                                                uiOutput("inFunControl"),
                                                conditionalPanel(
                                                  condition = "input.testFunControl === 'B'",
                                                  uiOutput("moreInputsControl")
                                                )),
                                         column(6,
                                                radioButtons("testFunExperiment", "Choose a function for Experiment",
                                                             choices = list("Piecewise exponential" = "B","Weibull" = "C","Mixture cure rate of exponential"="D",
                                                                            "Mixture cure rate of weibull"='E',"Customized distribution" = "A")),
                                                uiOutput("inFunExperiment"),
                                                conditionalPanel(
                                                  condition = "input.testFunExperiment === 'B'",
                                                  uiOutput("moreInputsExperiment")
                                                ))
                                       ),
                                       h6("To customize function: Example 1: if follows an exponential distribution with a median of 12 months, then the survival function is exp(-log(2)/12 * t).
                                           Example 2: if follows a piecewise exponential function with hazards log(2)/12 before 6 months and log(2)/24 after 6 months, which is log(2)/12*as.numeric(t<=6) + log(2)/24*as.numeric(t>6).
                                          The survival function is exp(-log(2)/12*t)*as.numeric(t <= 6) + exp(-log(2)/12 * 6 - log(2)/24 * t + log(2)/24 *6)*as.numeric(t>6)"),
                                       
                                       verbatimTextOutput("funcDef0"),
                                       actionButton("plotBtn0", "Display Survival Curves"),
                                       plotOutput("sur")
                                       
                              ),
               
                              tabPanel("Accrual / Dropouts",
                                       #fluidRow(column(12,
                                                       #selectInput('saf',
                                                                  # 'Subjects are followed',
                                                                  # choices = c('Until end of study','For fixed period')))),
                                       fluidRow(column(6,
                                                       h4('Accrual Info'),
                                                       radioButtons("accFun", "Choose a method for accrual",
                                                                    choices = list("Cumulative recruitment distribution" = "A", "Customized recruitment distribution " = "B")),
                                                       
                                                       uiOutput("inACC"),
                                                       actionButton("plotBtn2", "Display Accural Curve"),
                                                       conditionalPanel(
                                                         condition = "output.showPlot2",
                                                         plotOutput("accCurve"))
                                                       ),
                                                column(6,
                                                       h4('Dropout Rate'),
                                                       h6("* For example, 3% dropoff in 12 months, input 0.03"),
                                                       textInput('conDP','Percentage of dropout in 12 months for control'),
                                                       textInput('expDP','Percentage of dropout in 12 months for experiment'),
                                                       
                                                       actionButton("plotBtn3", "Display Dropout Curve"),
                                                       conditionalPanel(
                                                         condition = "output.showPlot3",
                                                         plotOutput("dpCurve1"),
                                                         plotOutput("dpCurve2"))
                                                       )
                                                )#,
                                       #fluidRow(column(6,
                                       #                verbatimTextOutput("functionDef")
                                       #                ))
                                       ),
                              tabPanel("Statistical Test",
                                       uiOutput("stOut"),
                                       verbatimTextOutput("fDef")
                                       ),
                              tabPanel("Power",
                                       sidebarPanel(
                                         textInput("ratio",
                                                   "Allocation ratio (experiment/control)"),
                                         h6("* For example, allocation ratio = 1 for 1:1 randomization. Allocation ratio = 2 for 2:1 randomization."),
                                         textInput("n",
                                                   "Sample size"),
                                         uiOutput("dcoInput"),
                                         selectInput("mu", "Choose the method for mean of weighted logrank test",
                                                     choices = c("Schoenfeld","H1")),
                                         selectInput("cov", "Choose the method for covariance matrix calculation",
                                                     choices = c("H0", "H1", "H1.LA")),
                                         actionButton("prun", "Calculate Power")
                                          
                                       ),
                                       mainPanel(
                                         dataTableOutput("powerResults"),
                                         plotOutput("powerPlot1"),
                                         plotOutput("powerPlot2")
                                       )
                                       
                                       ),
                              tabPanel("Events",
                                       fluidRow(
                                         column(4,
                                                numericInput("tmax","Maximum time",value=50),
                                                actionButton("eve","Calculate events"))
                                       ),
                                       fluidRow(
                                         column(6,
                                                dataTableOutput("eTable")
                                                ),
                                         column(6,
                                                plotOutput("eplot"))
                                       )
                                
                              ),
                              tabPanel("Boundary",
                                       fluidRow(
                                         selectInput(
                                           "bu",
                                           "Select the boundary to plot",
                                           choices = c("Z-scale","p-value(one-sided)","Hazard ratio")
                                         ),
                                         actionButton("bButton","Display Boundary Plot")
                                       ),
                                       fluidRow(
                                         column(9,
                                                plotOutput("bPlot"))
                                       )
                                       
                                       ),
                              tabPanel("Actual Alpha Spending",
                                       column(9,
                                              plotOutput("aaPlot"))
                                       )
                            )),
                   tabPanel('Simulation',
                            h3("Here's the simulation."),
                            tabsetPanel(
                              tabPanel("Type I error",
                                       fluidRow(
                                         
                                         column(12,
                                                selectInput("lookS", 
                                                            "Number of analyses that require alpha spending",
                                                            choices = c(1:5)),
                                                
                                                textInput("t1eS",
                                                          "Overall Type I error (one-sided), eg. 0.025 for single primary"),
                                                fluidRow(
                                                  column(6,
                                                         selectInput("spFunS",
                                                                     "Alpha Spending Function",
                                                                     choices=c("Lan-DeMets O'BrienFleming","Lan-DeMets Pocock","Hwang-Shih-DeCani",'Haybittle-Peto',"Bespoke"))
                                                  ),
                                                  column(6,uiOutput("paramInputS")))#,
                                                
                                               # actionButton("alpha", "Display IA Alpha Spending Curve for Exploration"),
                                               # plotOutput("as")
                                         )) 
                              ),
                              tabPanel("Distributions", 
                                       fluidRow(
                                         column(6,
                                                radioButtons("testFunC", "Choose a function for Control",
                                                             choices = list("Piecewise exponential" = "B","Weibull" = "C","Mixture cure rate of exponential"="D",
                                                                            "Mixture cure rate of weibull"='E',"Customized distribution" = "A")),
                                                uiOutput("inFunC"),
                                                conditionalPanel(
                                                  condition = "input.testFunC === 'B'",
                                                  uiOutput("moreInputsC")
                                                )),
                                         column(6,
                                                radioButtons("testFunE", "Choose a function for Experiment",
                                                             choices = list("Piecewise exponential " = "B","Weibull" = "C","Mixture cure rate of exponential"="D",
                                                                            "Mixture cure rate of weibull"='E',"Customized distribution" = "A")),
                                                uiOutput("inFunE"),
                                                conditionalPanel(
                                                  condition = "input.testFunE === 'B'",
                                                  uiOutput("moreInputsE")
                                                ))
                                       ),
                                       h6("To customize function: Example 1: if follows an exponential distribution with a median of 12 months, then the survival function is exp(-log(2)/12 * t).
                                           Example 2: if follows a piecewise exponential function with hazards log(2)/12 before 6 months and log(2)/24 after 6 months, which is log(2)/12*as.numeric(t<=6) + log(2)/24*as.numeric(t>6).
                                          The survival function is exp(-log(2)/12*t)*as.numeric(t <= 6) + exp(-log(2)/12 * 6 - log(2)/24 * t + log(2)/24 *6)*as.numeric(t>6)"),
                                       verbatimTextOutput("funcDef")
                                       
                              ),
                              
                              tabPanel("Accrual / Dropouts",
                                       fluidRow(column(6,
                                                       h4('Accrual Info'),
                                                       radioButtons("accFunS", "Choose a method for accrual",
                                                                    choices = list("Cumulative recruitment distribution" = "A", "Customized recruitment distribution " = "B")),
                                                       
                                                       uiOutput("inACCS")
                                       ),
                                       column(6,
                                              h4('Dropout Rate'),
                                              h6("* For example, 3% dropoff in 12 months, input 0.03"),
                                              textInput('conDPS','Percentage of dropout in 12 months for control'),
                                              textInput('expDPS','Percentage of dropout in 12 months for experiment')
                                            
                                       )
                                       ),
                                       fluidRow(column(6,
                                                       verbatimTextOutput("functionDef")
                                                       ))
                              ),
                              tabPanel("Statistical Test",
                                       uiOutput("stOutS"),
                                       verbatimTextOutput("fDefS")
                              ),
                              tabPanel("Simulation",
                                       sidebarPanel(
                                         textInput("nSim","Number of trials"),
                                         textInput("ratioS",
                                                   "Allocation ratio (experiment/control)"),
                                         h6("* For example, allocation ratio = 1 for 1:1 randomization. Allocation ratio = 2 for 2:1 randomization."),
                                         textInput("nS",
                                                   "Sample size"),
                                         uiOutput("eveInputS"),
                                         selectInput("lr","Log-rank test requested",choices = c('N','Y')),
                                         actionButton("prunS", "Start simulation")
                                         
                                       ),
                                       mainPanel(
                                         dataTableOutput("power")
                                       )
                                       
                              )
                            )
                   
                   ))

)