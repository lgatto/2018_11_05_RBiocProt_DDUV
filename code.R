## Environment

suppressPackageStartupMessages(library("MSnbase"))
library("pRoloc")
library("pRolocdata")
library("pRolocGUI")
library("magrittr")

## Loading data into R

f0 <- dir(system.file("extdata", package = "pRolocdata"),
          full.names = TRUE,
          pattern = "Dunkley2006")
basename(f0)

res <- readMSnSet2(f0, ecol = 5:20)
res

## Accessing data and meta-data

data(hyperLOPIT2015)
x <- hyperLOPIT2015
exprs(x)
pData(x)
fData(x)


## Missing data
data(naset)
naplot(naset)

filterNA(naset, pNA = 0.5) %>%
    naplot

naset %>%
    filterNA( pNA = 0.5) %>%
    impute(method = "knn") %>%
    naplot

res <- naset %>%
    filterNA( pNA = 0.5) %>%
    impute(method = "knn")

## Normalisation

par(mfrow = c(1, 2))
boxplot(exprs(res))

res %>%
    normalise(method = "vsn") %>%
    exprs %>%
    boxplot

## Dimensionality reduction

hl <- hyperLOPIT2015
plot2D(hl)


fvarLabels(hl)
plot2D(hl, fcol = "final.assignment")

ptz <- exp(fData(hl)$svm.score) - 1
plot2D(hl, fcol = "final.assignment", cex = ptz)

setStockcol(paste0(getStockcol(), 80))
plot2D(hl, fcol = "final.assignment", cex = ptz)

## Classification

data(tan2009r1)

## params <- svmOptimisation(tan2009r1, fcol = "markers.orig",
##               times = 100, xval = 5,
##               verbose = FALSE)

fn <- dir(system.file("extdata", package = "pRoloc"),
      full.names = TRUE, pattern = "params.rda")
load(fn, verbose = TRUE)
params

plot(params)
levelPlot(params)

svmres <- svmClassification(tan2009r1, fcol = "markers.orig",
                            sigma = 1, cost = 1)
fvarLabels(svmres)

svmres <- svmClassification(tan2009r1, fcol = "markers.orig",
                            assessRes = params)
processingData(svmres)
fvarLabels(svmres)

## Interactive visualisation
library("pRolocGUI")
pRolocGUI::pRolocVis(hl)

browseURL("https://lgatto.shinyapps.io/christoforou2015/")
