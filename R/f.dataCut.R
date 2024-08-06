#' Data Cut Utility Function to Determine the Datasets According to Target Events
#' 
#' Determine the dataset for each analysis according to specified target event.
#' 
#' @param  data Input dataset
#' @param  targetEvents A vector of targetEvents for planned analyses
#' @param  DCO A vector of data cutoffs for planned analyses. Either targetEvents or DCO is required.
#' 
#' @return An object with the components:
#' \describe{
#' \item{[[1]]}{Dataset for First Analysis}
#' \item{[[2]]}{Dataset for Second Analysis}
#' \item{[[L]]}{Dataset for Lth Analysis}
#' }
#' 
#' @export 
#' 
f.dataCut = function(data, targetEvents = 397, DCO = NULL) {
  data0 = data
  data0.order <- data0[order(as.numeric(data0$calendarTime)), ] #order by calendar time
  data.event <- data0.order[data0.order$cnsr == 0,] #Events Only
  
  data.event$event.seq <- seq.int(nrow(data.event)) #event sequence number
  if(is.null(DCO)){
    #Data cutoff in calendar time added to the original dataframe as a variable
    data0$calendarCutoff = as.numeric(data.event$calendarTime[data.event$event.seq == targetEvents])
  } else {
    data0$calendarCutoff = DCO
  }
  data0$survTimeCut = ifelse(as.numeric(data0$calendarTime) <= as.numeric(data0$calendarCutoff), as.numeric(data0$survTime), as.numeric(data0$calendarCutoff) - as.numeric(data0$enterTime))
  data0$cnsrCut = ifelse(as.numeric(data0$calendarTime) <= as.numeric(data0$calendarCutoff), data0$cnsr, 1)
  
  return(data0)
}
