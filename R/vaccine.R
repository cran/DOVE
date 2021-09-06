#' Specify Vaccination Variables and Entry Time
#'
#' This function is used in the model statements of dove() and dove2() to specify
#'  the vaccination time, vaccination status, and entry time.
#'
#' For participants that were vaccinated, times must obey
#'   entry_time \eqn{\le} vaccination_time. If a case is
#'   found to violate this relationship, its entry_time is set to NA
#'   and it is removed from the analysis.
#'
#' @param entry_time The variable for the time when 
#'   the participant enters the trial. Entry times must be
#'   non-negative and complete.
#'
#' @param vaccination_status The variable indicating the vaccination 
#'   status: 1 = vaccinated; 0 = not vaccinated. Vaccination status
#'   must be binary, integer (or be able to be cast as integer without
#'   loss of information), and complete.
#'
#' @param vaccination_time The variable for the time when 
#'   vaccination takes place. Vaccination time must be non-negative for
#'   vaccinated participants and can be
#'   any non-negative value, NA or Inf if a participant was not 
#'   vaccinated during the trial.
#'
#' @returns This function is intended to be used only in the model statements 
#'  of dove() and dove2(). The result, a matrix, is used internally.
#'    
#' @name vaccine
#' @rdname vaccine
#' @export
#
# Note for developer: any vaccination times provided as NA are set to Inf.
# all non-vaccinated participants have vaccination times set to Inf.
# entry times are set to NA to indicate that time relationships are not
# satisfied. the returned matrix is therefore not necessarily complete 
# (though NA values will only be present for entry times) and
# any missing values indicate that the data should be excluded from the
# analysis

vaccine <- function(entry_time, vaccination_status, vaccination_time) {

  if (missing(x = entry_time) || 
      missing(x = vaccination_status) ||
      missing(x = vaccination_time)) {
    stop("must provide entry_time, vaccination_status, and vaccination_time", 
         call. = FALSE)
  }

  ### entry time
  
  # must be provided as a numeric vector.
  # must be complete
  # must be non-negative

  entry_time <- .basicTests_noNA(x = entry_time, name = "entry_time")

  ### time of vaccination
  
  # must be provided as a numeric vector.
  # can be incomplete (NA or Inf)
  # must be non-negative

  vaccination_time <- .basicTests_NAorInf(x = vaccination_time, 
                                          name = "vaccination_time")

  ### vaccination status

  # must be integer-like binary 0/1
  # must be complete

  vaccination_status <- .basicTests_noNA_binary(x = vaccination_status, 
                                                name = "vaccination_status")

  # set vaccination time to Inf for all participants not vaccinated
  vaccination_time[vaccination_status == 0L] <- Inf

  # ensure that vaccination times for vaccinated participants are well
  # defined
  if (any(is.infinite(vaccination_time[vaccination_status == 1L]))) {
    stop("encountered NA or Inf vaccination time for vaccinated participant",
         call. = FALSE)
  }

  # inform user if all participants are vaccinated
  if (all(vaccination_status == 1L)) {
    message("all participants received vaccine")
  }

  # must satisfy entry_time <= vaccination_time

  tst <- entry_time > {vaccination_time + 1e-8}

  violate <- sum(tst)

  if (violate > 0L) {
    message(violate, 
            ifelse(test = violate > 1L, 
                   yes = " cases do not ", 
                   no = " case does not "),
            "satisfy required entry_time <= vaccination_time relationship; ",
            ifelse(test = violate > 1L, 
                   yes = "cases removed", 
                   no = "case removed"))
    entry_time[tst] <- NA
  }

  if (length(x = entry_time) != length(x = vaccination_status) ||
      length(x = entry_time) != length(x = vaccination_time)) {
    stop("vaccine() inputs must be of same length", call. = FALSE)
  }
  
  dm <- cbind(entry_time, vaccination_time, vaccination_status)

  cname <- c("entryTime", "vacTime", "vacStatus")
  dimnames(x = dm) <- list(NULL, cname)

  class(x = dm) <- "vaccine"

  return( dm )
}


# test of input vector that must satisfy:
#   must be provided as a numeric vector.
#   must be complete
#   must be non-negative
.basicTests_noNA <- function(x, name) {

  if (anyNA(x = x)) {
    stop(name, " must be complete", call. = FALSE)
  }
  
  if (!is.numeric(x = x)) {
    stop(name, " must be a numeric vector", call. = FALSE)
  }

  if (any(x < 0.0)) {
    stop(name, " must be non-negative", call. = FALSE)
  }

  return( x )
}

# test of input vector that must satisfy:
#   must be provided as a numeric vector.
#   can be incomplete (NA or Inf w/ NA replaced by Inf)
#   must be non-negative
.basicTests_NAorInf <- function(x, name) {

  x[is.na(x = x)] <- Inf

  if (!is.numeric(x = x)) {
    stop(name, " must be a numeric vector", call. = FALSE)
  }

  if (any(x < 0L)) {
    stop(name, " must be non-negative", call. = FALSE)
  }

  return( x )
}

# test of input vector that must satisfy:
#   must be integer-like binary of value 0/1
#   must be complete
.basicTests_noNA_binary <- function(x, name) {

  if (anyNA(x = x)) {
    stop(name, " must be complete", call. = FALSE)
  }
  
  if (is.logical(x = x)) {
    # status provided as T/F - convert to integer
    x <- as.integer(x = x)
  } else if (is.factor(x = x)) {
    # status provided as a factor - convert to integer level ids
    x <- match(x = levels(x = x)[x],
               table = levels(x = x)) - 1L
  } else if (is.numeric(x = x)) {
    if (!is.integer(x = x)) {
      tmp <- as.integer(x = round(x = x, digits = 0L))
      if (!isTRUE(x = all.equal(target = x, current = tmp))) {
        stop(name, " is not integer", call. = FALSE)
      }
      x <- tmp
    }
  } else {
    stop(name, " is not integer", call. = FALSE)
  }

  # if coded as 1,2 rather than 0,1 shift coding
  if (max(x) == 2L) x <- x - 1L

  # ensure that coding is 0,1
  tst <- x %in% c(0L,1L)
  if (any(!tst)) {
    stop(" invalid ", name, " value encountered", call. = FALSE)
  }

  return( x )
}
