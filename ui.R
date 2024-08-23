

library(shiny)
library(DT)
library(nphRshiny)
library(tidyverse)
library(mvtnorm)
library(gsDesign)
library(plotly)
library(ggplot2)
library(parallel)  
library(foreach)   
library(doParallel)

 
shinyUI(navbarPage(title = tags$div(tags$img(src="https://i.postimg.cc/QCn8wk28/Trial-Nexus-Logo.png", height="35px"), 
                                     "TrialNexus"),
                   tabPanel(
                     'Home',
                     tags$div(
                       img(src = "https://i.postimg.cc/QCn8wk28/Trial-Nexus-Logo.png", height = "400px", width = "60%"),
                       style = "text-align: center; margin: 0 auto;"
                     ),
                     
                     h3('Introduction'),
                    
                     p("With significant scientific breakthroughs, targeted cancer drugs have entered clinical trials with complex mechanisms of action and unique features. Traditional statistical distributions, such as the simple exponential distribution for time-to-event endpoints, may struggle to adequately model the expected survival curves for these innovative cancer therapies."),
                     
                     p("To address this challenge, TrialNexus was created through a collaborative effort among industry colleagues. Currently, TrialNexus is in development, with ongoing updates provided by volunteer developers. The tool is not intended for formal protocol or regulatory submission purposes, and the developers assume no responsibility for its accuracy. If you’re interested in volunteering and contributing to the development, please contact us at ",
                       a("oncotrialdesign@gmail.com", href = "mailto:oncotrialdesign@gmail.com"), "."),
                     
                     p(" "), 
                     p("TrialNexus offers the following features:"),
                     
                     tags$ul(
                       tags$li("Flexible Distribution: Users can select from various distributions for each study arm, such as piecewise exponential, mixture cure rate, or customized survival distributions. The experimental arm also includes a convenient 'Proportional Hazards' option that can be applied alongside any survival distribution in the control arm."),
                       tags$li("Flexible Enrollment Curve: Users have the option to utilize common non-uniform pattern functions or create a customized enrollment curve to better suit their study's needs."),
                       tags$li("Flexible Statistical Testing at Each Analysis: At each analysis point, users can specify the type of statistical test to be applied, including the logrank test, the maximum of multiple weighted logrank tests, or even a customized weighted logrank test.")
                     ),
                     
                     p("These features are seamlessly integrated into both the study design and simulations modules. Additionally, TrialNexus includes a comprehensive collection of common statistical graphics used in study design.")
                   )
                   ,
                   tabPanel('Design',
                            h3('Group Sequential Design With Complex Survival Distribution'),
                        
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
                                                radioButtons("testFunControl", "Distribution for control arm",
                                                             choices = list("Exponential"="F","Piecewise exponential" = "B","Weibull" = "C","Mixture cure rate of exponential"="D",
                                                                            "Mixture cure rate of weibull"='E',"Customized distribution" = "A")),
                                                uiOutput("inFunControl"),
                                                conditionalPanel(
                                                  condition = "input.testFunControl === 'B'",
                                                  uiOutput("moreInputsControl")
                                                ),
                                                conditionalPanel(
                                                  condition = "input.testFunControl === 'F'",
                                                  uiOutput("medConOut")
                                                )),
                                         column(6,
                                                radioButtons("distEXP", "Distribution for experimental arm",choices = list("Proportional Hazards"="A","General"="B")),
                                                conditionalPanel(
                                                  condition = "input.distEXP === 'A'",
                                                  textInput("ratioExp","Hazard ratio")
                                                ),
                                                conditionalPanel(
                                                  condition = "input.distEXP === 'B'",
                                                  radioButtons("testFunExperiment", "General",
                                                               choices = list("Exponential"="F","Piecewise exponential" = "B","Weibull" = "C","Mixture cure rate of exponential"="D",
                                                                              "Mixture cure rate of weibull"='E',"Customized distribution" = "A")),
                                                  uiOutput("inFunExperiment"),
                                                  conditionalPanel(
                                                    condition = "input.testFunExperiment === 'B'",
                                                    uiOutput("moreInputsExperiment")
                                                  ),
                                                  conditionalPanel(
                                                    condition = "input.testFunExperiment === 'F'",
                                                    uiOutput("medExpOut")
                                                  )
                                                )
                                                )
                                       ),
                                       h6("To customize function: Example 1: if follows an exponential distribution with a median of 12 months, then the survival function is exp(-log(2)/12 * t).
                                           Example 2: if follows a piecewise exponential function with hazards log(2)/12 before 6 months and log(2)/24 after 6 months, which is log(2)/12*as.numeric(t<=6) + log(2)/24*as.numeric(t>6).
                                          The survival function is exp(-log(2)/12*t)*as.numeric(t <= 6) + exp(-log(2)/12 * 6 - log(2)/24 * t + log(2)/24 *6)*as.numeric(t>6)"),
                                       
                                       #verbatimTextOutput("funcDef0"),
                                       fluidRow(
                                         column(6,
                                                actionButton("plotBtn0", "Display Survival Curves"),
                                                plotOutput("sur")),
                                         column(6,
                                                actionButton("plotBtn", "Display Hazard Curves"),
                                                plotOutput("hz"))
                                       )
                                       
                                       
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
                                                       textInput('conDP','Control Arm: Probability of drop out in 12 months'),
                                                       textInput('expDP','Experiment Arm: Probability of drop out in 12 months'),
                                                       
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
                                       uiOutput("stOut")#,
                                       #verbatimTextOutput("fDef")
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
                                         plotlyOutput("powerPlot1"),
                                         plotlyOutput("powerPlot2")
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
                                                plotlyOutput("bPlot"))
                                       )
                                       
                                       ),
                              tabPanel("Actual Alpha Spending",
                                       column(9,
                                              plotlyOutput("aaPlot"))
                                       ),
                              tabPanel("Average Hazard Ratio",
                                       fluidRow(
                                         numericInput("tmaxA","Maximum time",value=50),
                                         selectInput(
                                           "ahr",
                                           "Select the method",
                                           choices = c("Geometric Schoenfeld","Kalbfleisch and Prentice")
                                         ),
                                         actionButton("ahrButton","Display Average Hazard Ratio Plot")
                                       ),
                                       fluidRow(
                                         column(9,
                                                plotOutput("ahrPlot"))
                                       )
                                       
                              )
                            )),
                   tabPanel('Simulation',
                            h3("Simulation For Group Sequential Trial With Complex Survival Distribution"),
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
                                                             choices = list("Exponential"="F","Piecewise exponential" = "B","Weibull" = "C","Mixture cure rate of exponential"="D",
                                                                            "Mixture cure rate of weibull"='E',"Customized distribution" = "A")),
                                                uiOutput("inFunC"),
                                                conditionalPanel(
                                                  condition = "input.testFunC === 'B'",
                                                  uiOutput("moreInputsC")
                                                ),
                                                conditionalPanel(
                                                  condition = "input.testFunC === 'F'",
                                                  uiOutput("medConOutS")
                                                )),
                                         column(6,
                                                radioButtons("distEXPS", "Distribution for experimental arm",choices = list("Proportional Hazards"="A","General"="B")),
                                                conditionalPanel(
                                                  condition = "input.distEXPS === 'A'",
                                                  textInput("ratioExpS","Hazard ratio")
                                                ),
                                                conditionalPanel(
                                                  condition = "input.distEXPS === 'B'",
                                                  radioButtons("testFunE", "Choose a function for Experiment",
                                                               choices = list("Exponential"="F","Piecewise exponential " = "B","Weibull" = "C","Mixture cure rate of exponential"="D",
                                                                              "Mixture cure rate of weibull"='E',"Customized distribution" = "A")),
                                                  uiOutput("inFunE"),
                                                  conditionalPanel(
                                                    condition = "input.testFunE === 'B'",
                                                    uiOutput("moreInputsE")
                                                  ),
                                                  conditionalPanel(
                                                    condition = "input.testFunE === 'F'",
                                                    uiOutput("medExpOutS")
                                                  ))
                                                )
                                               
                                       ),
                                       h6("To customize function: Example 1: if follows an exponential distribution with a median of 12 months, then the survival function is exp(-log(2)/12 * t).
                                           Example 2: if follows a piecewise exponential function with hazards log(2)/12 before 6 months and log(2)/24 after 6 months, which is log(2)/12*as.numeric(t<=6) + log(2)/24*as.numeric(t>6).
                                          The survival function is exp(-log(2)/12*t)*as.numeric(t <= 6) + exp(-log(2)/12 * 6 - log(2)/24 * t + log(2)/24 *6)*as.numeric(t>6)")#,
                                       #verbatimTextOutput("funcDef")
                                       
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
                                       )
                              ),
                              tabPanel("Statistical Test",
                                       uiOutput("stOutS")#,
                                       #verbatimTextOutput("fDefS")
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
                                         selectInput("parS","Parallel requested",choices = c('F','T')),
                                         selectInput("coreS","Number of cores",choices = c(1:4)),
                                         numericInput("seed","Seed",value=2022),
                                         actionButton("prunS", "Start simulation")
                                         
                                       ),
                                       mainPanel(
                                         dataTableOutput("power"),
                                         plotlyOutput('sim1')
                                       )
                                       
                              )
                            )
                   
                   ))

)