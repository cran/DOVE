#' Specify Vaccination Variables and Entry Time
#'
#' This function is used in the model statement of dove() to specify
#'  the vaccination time, vaccination status, and entry time.
#'
#' @param entry_time The variable for the time when 
#'   the participant enters the trial. 
#'
#' @param vaccination_status The variable indicating the vaccination 
#'   status: 1 = vaccinated; 0 = not vaccinated.
#'
#' @param vaccination_time The variable for the time when 
#'   vaccination takes place. 
#'
#' @returns This function is intended to be used only in the model statement 
#'  of dove(). The result, a matrix, is used internally.
#'    
#' @name vaccine
#' @rdname vaccine
#' @export

vaccine <- function(entry_time, vaccination_status, vaccination_time) {

  ### entry time

  # must be provided as a numeric vector. 

  if (missing(x = entry_time)) {
    stop("must provide a entry_time argument", 
         call. = FALSE)
  }

  if (!is.numeric(x = entry_time)) {
    stop ("entry_time is not numeric", call. = FALSE)
  }

  ### time of vaccination

  # must be provided as a numeric vector. 

  if (missing(x = vaccination_time)) {
    stop("must provide a vaccination_time argument", 
         call. = FALSE)
  }

  if (!is.numeric(x = vaccination_time)) {
    stop ("vaccination_time is not numeric", call. = FALSE)
  }

  if (length(x = entry_time) != length(x = vaccination_time)) {
    stop("all inputs must be of same length", call. = FALSE)
  }

  ### vaccination status

  # must be provided as a numeric or logical vector

  if (missing(x = vaccination_status)) {
    stop("must provide a vaccination_status argument", call. = FALSE)
  }
  if (length(x = vaccination_status) != length(x = vaccination_time)) {
    stop("all inputs must be of same length", call. = FALSE)
  }
  if (any(is.na(x = vaccination_status))) {
    stop("vaccination_status cannot contain missing data", call. = FALSE)
  }

  if (is.logical(x = vaccination_status)) {
    # status provided as T/F - convert to integer
    stat <- as.integer(x = vaccination_status)
  } else if (is.factor(x = vaccination_status)) {
    # status provided as a factor - convert to integer level ids
    stat <- match(x = levels(x = vaccination_status)[vaccination_status],
                 table = levels(x = vaccination_status)) - 1L
  } else if (is.numeric(x = vaccination_status)) {
    # status provided as numeric - convert to integer
    stat <- as.integer(x = round(x = vaccination_status, digits = 0L))
  } else {
    stop("invalid vaccination_status input, must be logical or numeric", 
         call. = FALSE)
  }

  # if coded as 1,2 rather than 0,1 shift coding
  if (max(stat, na.rm = TRUE) == 2L) stat <- stat - 1L

  # ensure that coding is 0,1
  temp <- stat == 0L | stat == 1L
  if (any(!temp & !is.na(x = stat))) {
    stop(" invalid vaccination_status value encountered", call. = FALSE)
  }
  stat <- ifelse(test = temp, yes = stat, no = NA)

  vaccination_time[stat == 0L] <- NA

  tst <- {vaccination_time < entry_time} & {stat == 1L}

  if (any(tst, na.rm = TRUE)) {
    message(sum(tst, na.rm = TRUE), 
            " have vaccination_time < entry_time; cases set to NA")
    entry_time[tst] <- NA
  }

  if (all(stat == 1L, na.rm = TRUE)) {
    message("all participants received vaccine")
  }

  dm <- cbind(entry_time, vaccination_time, stat)

  cname <- c("entry_time", "vaccination_time", "vaccination_status")
  dimnames(x = dm) <- list(NULL, cname)

  return( dm )
}
