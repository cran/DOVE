#' Toy Dataset For Illustration
#'
#' This data set is provided for the purposes of illustrating the use of
#' the software. It was simulated under a priority-tier dependent crossover 
#' design.
#' 
#' 
#' @usage data(doveData)
#'
#' @format doveData is a data.frame containing 40,000 participants The
#'   data.frame contains 7 columns, 
#'   \describe{
#'   \item{entry.time}{The entry time in days}
#'   \item{event.time}{The event time in days}
#'   \item{event.status}{The event indicator (1=event; 0=censored)}
#'   \item{vaccine.time}{The time of vaccination in days; 
#'                       NA if not vaccinated}
#'   \item{vaccine.status}{The indicator of vaccination 
#'                         (1 = vaccinated; 0 = not vaccinated)}
#'   \item{priority}{A composite baseline risk score taking values 1-5}
#'   \item{sex}{A binary indicator of sex (male/female)}
#'   }
#'
#' @name doveData
#' @keywords datasets
NULL
