## ----setup, echo = FALSE, results = "hide", message = FALSE, warnings = FALSE----
suppressWarnings(RNGversion("3.5.2"))
pkgs <- c("partykit", "coin", "strucchange", "Formula", "survival",
          "sandwich", "party", "TH.data", "knitr")
pkgs <- sapply(pkgs, require, character.only = TRUE)
if (pkgs["party"]) detach(package:party)
set.seed(290875)

## ----fail, results = "asis", echo = FALSE---------------------------
if (any(!pkgs))
{
    cat(paste("Package(s)", paste(names(pkgs)[!pkgs], collapse = ", "), 
        "not available, stop processing.",
        "\\end{document}\n"))
    knitr::knit_exit()
}

## ----party-max, echo = TRUE, results = "hide"-----------------------
ctree_control(teststat = "max")

## ----party-quad, echo = TRUE, results = "hide"----------------------
ctree_control(teststat = "quad")

## ----party-Bonf, echo = TRUE, results = "hide"----------------------
ctree_control(testtype = "Bonferroni")

## ----party-minsplit, echo = TRUE, results = "hide"------------------
ctree_control(minsplit = 20)

## ----party-maxsurr, echo = TRUE, results = "hide"-------------------
ctree_control(maxsurrogate = 3)

## ----party-data, echo = TRUE----------------------------------------
ls <- data.frame(y = gl(3, 50, labels = c("A", "B", "C")), 
                 x1 = rnorm(150) + rep(c(1, 0, 0), c(50, 50, 50)), 
                 x2 = runif(150))

## ----party-formula, echo = TRUE, results = "hide"-------------------
library("partykit")
ctree(y ~ x1 + x2, data = ls)

## ----party-fitted, echo = TRUE--------------------------------------
ct <- ctree(y ~ x1 + x2, data = ls)

## ----party-print, echo = TRUE---------------------------------------
ct

## ----party-plot, echo = TRUE,  fig.height = 5, fig.width = 8--------
plot(ct)

## ----party-nodes, echo = TRUE---------------------------------------
ct[1]

## ----party-nodelist, echo = TRUE------------------------------------
class(ct[1])

## ----party-predict, echo = TRUE-------------------------------------
predict(ct, newdata = ls)

## ----party-treeresponse, echo = TRUE--------------------------------
predict(ct, newdata = ls[c(1, 51, 101),], type = "prob")

## ----party-where, echo = TRUE---------------------------------------
predict(ct, newdata = ls[c(1,51,101),], type = "node")

## ----party-sctest, echo = TRUE--------------------------------------
if (require("strucchange"))
    print(sctest(ct))

## ----treepipit-ctree, echo = TRUE-----------------------------------
data("treepipit", package = "coin")
tptree <- ctree(counts ~ ., data = treepipit)

## ----treepipit-plot, echo = TRUE,  fig.height = 5, fig.width = 8----
plot(tptree, terminal_panel = node_barplot)

## ----treepipit-x, echo = FALSE--------------------------------------
p <- info_node(node_party(tptree))$p.value
n <- table(predict(tptree, type = "node"))

## ----glaucoma-ctree, echo = TRUE------------------------------------
data("GlaucomaM", package = "TH.data")
gtree <- ctree(Class ~ ., data = GlaucomaM)

## ----glaucoma-x, echo = FALSE, results = "hide"---------------------
sp <- split_node(node_party(gtree))$varID

## ----glaucoma-plot, echo = FALSE,  fig.height = 6, fig.width = 10----
plot(gtree)

## ----glaucoma-plot-inner, echo = FALSE,  fig.height = 7, fig.width = 10----
plot(gtree, inner_panel = node_barplot, 
     edge_panel = function(...) invisible(), tnex = 1)

## ----glaucoma-plot2, echo = TRUE, eval = FALSE----------------------
# plot(gtree)

## ----glaucoma-plot-inner-bar, echo = TRUE, eval = FALSE-------------
# plot(gtree, inner_panel = node_barplot,
#      edge_panel = function(...) invisible(), tnex = 1)

## ----glaucoma-prediction, echo = TRUE-------------------------------
table(predict(gtree), GlaucomaM$Class)

## ----glaucoma-classprob, echo = TRUE,  fig.height = 4, fig.width = 5----
prob <- predict(gtree, type = "prob")[,1] + 
                runif(nrow(GlaucomaM), min = -0.01, max = 0.01)
splitvar <- character_split(split_node(node_party(gtree)), 
                            data = data_party(gtree))$name
plot(GlaucomaM[[splitvar]], prob, 
     pch = as.numeric(GlaucomaM$Class), ylab = "Conditional Class Prob.",
     xlab = splitvar)
abline(v = split_node(node_party(gtree))$breaks, lty = 2)
legend(0.15, 0.7, pch = 1:2, legend = levels(GlaucomaM$Class), bty = "n")

## ----GBSGS-ctree, echo = TRUE---------------------------------------
data("GBSG2", package = "TH.data")  
library("survival")
(stree <- ctree(Surv(time, cens) ~ ., data = GBSG2))

## ----GBSG2-plot, echo = TRUE,  fig.height = 6, fig.width = 10-------
plot(stree)

## ----GBSG2-KM, echo = TRUE------------------------------------------
pn <- predict(stree, newdata = GBSG2[1:2,], type = "node")
n <- predict(stree, type = "node")
survfit(Surv(time, cens) ~ 1, data = GBSG2, subset = (n == pn[1]))
survfit(Surv(time, cens) ~ 1, data = GBSG2, subset = (n == pn[2]))

## ----mammo-ctree, echo = TRUE---------------------------------------
data("mammoexp", package = "TH.data")
mtree <- ctree(ME ~ ., data = mammoexp)

## ----mammo-plot, echo = TRUE,  fig.height = 6, fig.width = 13-------
plot(mtree)

## ----spider-ctree, echo = TRUE--------------------------------------
data("HuntingSpiders", package = "partykit")
sptree <- ctree(arct.lute + pard.lugu + zora.spin + pard.nigr +
  pard.pull + aulo.albi + troc.terr + alop.cune + pard.mont + alop.acce +
  alop.fabr + arct.peri ~ herbs + reft + moss + sand + twigs + water,
  data = HuntingSpiders, teststat = "max", minsplit = 5,
  pargs = GenzBretz(abseps = .1, releps = .1))

## ----spider-plot1, echo = TRUE,  fig.height = 6, fig.width = 13-----
plot(sptree, terminal_panel = node_barplot)

## ----spider-plot2, echo = TRUE,  fig.width = 10, fig.height = 22----
plot(sptree)

## ----party-setup, echo = FALSE, results = "hide", message = FALSE, warnings = FALSE----
library("party")
set.seed(290875)

## ----party-airq, echo = TRUE----------------------------------------
data("airquality", package = "datasets")
airq <- subset(airquality, !is.na(Ozone))
(airct_party <- party::ctree(Ozone ~ ., data = airq, 
     controls = party::ctree_control(maxsurrogate = 3)))
mean((airq$Ozone - predict(airct_party))^2)

## ----partykit-airq, echo = TRUE-------------------------------------
(airct_partykit <- partykit::ctree(Ozone ~ ., data = airq, 
     control = partykit::ctree_control(maxsurrogate = 3, 
                                       numsurrogate = TRUE)))
mean((airq$Ozone - predict(airct_partykit))^2)
table(predict(airct_party, type = "node"),
      predict(airct_partykit, type = "node"))
max(abs(predict(airct_party) - predict(airct_partykit)))

## ----party-partykit-airq, echo = TRUE-------------------------------
airct_party@tree$criterion
info_node(node_party(airct_partykit))

## ----party-partykit-iris, echo = TRUE-------------------------------
(irisct_party <- party::ctree(Species ~ .,data = iris))
(irisct_partykit <- partykit::ctree(Species ~ .,data = iris, 
control = partykit::ctree_control(splitstat = "maximum")))
table(predict(irisct_party, type = "node"),
      predict(irisct_partykit, type = "node"))

## ----party-iris, echo = TRUE----------------------------------------
tr_party <- treeresponse(irisct_party, newdata = iris)

## ----partykit-iris, echo = TRUE-------------------------------------
tr_partykit <- predict(irisct_partykit, type = "prob", 
                       newdata = iris)
max(abs(do.call("rbind", tr_party) - tr_partykit))

## ----party-partykit-mammoexp, echo = TRUE---------------------------
### ordinal regression
data("mammoexp", package = "TH.data")
(mammoct_party <- party::ctree(ME ~ ., data = mammoexp))
### estimated class probabilities
tr_party <- treeresponse(mammoct_party, newdata = mammoexp)
(mammoct_partykit <- partykit::ctree(ME ~ ., data = mammoexp))
### estimated class probabilities
tr_partykit <- predict(mammoct_partykit, newdata = mammoexp, type = "prob")
max(abs(do.call("rbind", tr_party) - tr_partykit))

## ----party-partykit-GBSG2, echo = TRUE------------------------------
data("GBSG2", package = "TH.data")
(GBSG2ct_party <- party::ctree(Surv(time, cens) ~ .,data = GBSG2))
(GBSG2ct_partykit <- partykit::ctree(Surv(time, cens) ~ .,data = GBSG2))

## ----party-partykit-KM, echo = TRUE---------------------------------
tr_party <- treeresponse(GBSG2ct_party, newdata = GBSG2)
tr_partykit <- predict(GBSG2ct_partykit, newdata = GBSG2, type = "prob")
all.equal(lapply(tr_party, function(x) unclass(x)[!(names(x) %in% "call")]),
lapply(tr_partykit, function(x) unclass(x)[!(names(x) %in% "call")]),
check.names = FALSE)

## ----nf-alpha, echo = TRUE------------------------------------------
(airct_partykit_1 <- partykit::ctree(Ozone ~ ., data = airq, 
     control = partykit::ctree_control(maxsurrogate = 3, alpha = 0.001,
                                       numsurrogate = FALSE)))
depth(airct_partykit_1)
mean((airq$Ozone - predict(airct_partykit_1))^2)

## ----nf-maxstat, echo = TRUE----------------------------------------
(airct_partykit <- partykit::ctree(Ozone ~ ., data = airq, 
     control = partykit::ctree_control(maxsurrogate = 3, splittest = TRUE,
                                       testtype = "MonteCarlo"))) 

## ----nf-nmax, echo = TRUE-------------------------------------------
(irisct_partykit_1 <- partykit::ctree(Species ~ .,data = iris, 
control = partykit::ctree_control(splitstat = "maximum", nmax = 25)))
table(predict(irisct_partykit), predict(irisct_partykit_1))

## ----nf-multiway, echo = TRUE---------------------------------------
GBSG2$tgrade <- factor(GBSG2$tgrade, ordered = FALSE)
(GBSG2ct_partykit <- partykit::ctree(Surv(time, cens) ~ tgrade,
    data = GBSG2, control = partykit::ctree_control(multiway = TRUE, 
                                                    alpha = .5)))

## ----nf-cluster, echo = TRUE----------------------------------------
airq$month <- factor(airq$Month)
(airct_partykit_3 <- partykit::ctree(Ozone ~ Solar.R + Wind + Temp, data = airq,
cluster = month, control = partykit::ctree_control(maxsurrogate = 3)))
info_node(node_party(airct_partykit_3))
mean((airq$Ozone - predict(airct_partykit_3))^2)

## ----nf-ytrafo-1, echo = TRUE---------------------------------------
### with weight-dependent log-rank scores
### log-rank trafo for observations in this node only (= weights > 0)
h <- function(y, x, start = NULL, weights, offset, estfun = TRUE, object = FALSE, ...) {
    if (is.null(weights)) weights <- rep(1, NROW(y))
    s <- logrank_trafo(y[weights > 0,,drop = FALSE])
    r <- rep(0, length(weights))
    r[weights > 0] <- s
    list(estfun = matrix(as.double(r), ncol = 1), converged = TRUE, unweighted = TRUE)
}
partykit::ctree(Surv(time, cens) ~ ., data = GBSG2, ytrafo = h)

## ----nf-ytrafo-2, echo = TRUE---------------------------------------
### normal varying intercept / varying coefficient model (aka "mob")
h <- function(y, x, start = NULL, weights = NULL, offset = NULL, cluster = NULL, ...)
  glm(y ~ 0 + x, family = gaussian(), start = start, weights = weights, ...)
(airct_partykit_4 <- partykit::ctree(Ozone ~ Temp | Solar.R + Wind, 
    data = airq, cluster = month, ytrafo = h, 
    control = partykit::ctree_control(maxsurrogate = 3)))
airq$node <- factor(predict(airct_partykit_4, type = "node"))
summary(m <- glm(Ozone ~ node + node:Temp - 1, data = airq))
mean((predict(m) - airq$Ozone)^2)

## ----airq-mob, echo = TRUE------------------------------------------
airq_lmtree <- partykit::lmtree(Ozone ~ Temp | Solar.R + Wind, 
                                data = airq, cluster = month)
info_node(node_party(airq_lmtree))
mean((predict(airq_lmtree, newdata = airq) - airq$Ozone)^2)

## ----closing, echo = FALSE, results = "hide"------------------------
detach(package:party)

## ----bib, echo = FALSE----------------------------------------------
thisdir <- getwd()
bibfile <- system.file("REFERENCES.bib", package = "partykit")
### bibfile may contain spaces LaTeX is unable to deal with on MacOS it seems
if (file.copy(bibfile, to = thisdir, overwrite = TRUE)) {
    bibfile <- "REFERENCES.bib"
} else {
    ### hope for the best
    bibfile <- file.path("..", "inst", "REFERENCES.bib")
}

