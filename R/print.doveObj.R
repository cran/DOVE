#-------------------------------------------------------------------------------
# method to print key results to screen
#-------------------------------------------------------------------------------
# method is exported
#-------------------------------------------------------------------------------
# Print Analysis Results
#
# Prints the key results of the analysis.
#
# @param x A doveObj object. The value returned by dove().
#
# @param ... Ignored. 
#
# @name print
# @examples
# data(doveData)
#
# set.seed(1234)
#
# ind <- sample(1:nrow(x = doveData), 2500, FALSE)
#
# # NOTE: This sample size is chosen for example only -- larger data sets
# # should be used.
# # See the vignette for a full analysis of the doveData dataset
#
# res <- dove(formula = Surv(event.time, event.status) ~ priority + sex + 
#                        vaccine(entry.time, vaccine.status, vaccine.time),
#             data = doveData[ind,])
#
# print(x = res)
# @method print doveObj
# @export
print.doveObj <- function(x, ...) {
  attributes(x = x) <- NULL
  x <- unclass(x = x)
  print(x = x, ...)
}
