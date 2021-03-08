### R code from vignette source 'dove-vignette.Rnw'

###################################################
### code chunk number 1: init
###################################################
origOpts <- options()
on.exit(options(origOpts))
options(continue="  ", width=70, prompt=" ")
options(SweaveHooks=list(fig=function() par(mar=c(4.1, 4.1, .3, 1.1))))
pdf.options(pointsize=8) #text in graph about the same as regular text
library(DOVE, quietly=TRUE)


###################################################
### code chunk number 2: approx1 (eval = FALSE)
###################################################
## vaccine(entry_time, vaccination_status, vaccination_time)


###################################################
### code chunk number 3: doveStructure (eval = FALSE)
###################################################
## dove(formula, data, plots = TRUE, timePts = NULL, bandwidth = NULL)


###################################################
### code chunk number 4: data
###################################################
data(doveData)


###################################################
### code chunk number 5: dataShow
###################################################
head(doveData)


###################################################
### code chunk number 6: dataSummary
###################################################
summary(doveData)


###################################################
### code chunk number 7: result
###################################################
result <- dove(formula = Surv(event.time, event.status) ~ priority + sex + 
                 vaccine(entry.time, vaccine.status, vaccine.time), 
               data = doveData)


###################################################
### code chunk number 8: beta
###################################################
result$covariates


###################################################
### code chunk number 9: cumulative
###################################################
head(result$vaccine$efficacy)
tail(result$vaccine$efficacy)


###################################################
### code chunk number 10: period
###################################################
result$vaccine$period_efficacy


