## ----setup, echo = FALSE, results = "hide", message = FALSE, warnings = FALSE----
suppressWarnings(RNGversion("3.5.2"))
options(width = 70)
library("partykit")
set.seed(290875)

## ----weather-data---------------------------------------------------
data("WeatherPlay", package = "partykit")
WeatherPlay

## ----weather-plot0, echo=FALSE, fig.height=5, fig.width=7.5---------
py <- party(
  partynode(1L,
    split = partysplit(1L, index = 1:3),
    kids = list(
      partynode(2L,
        split = partysplit(3L, breaks = 75),
        kids = list(
          partynode(3L, info = "yes"),
          partynode(4L, info = "no"))),
      partynode(5L, info = "yes"),
      partynode(6L,
        split = partysplit(4L, index = 1:2),
        kids = list(
          partynode(7L, info = "yes"),
          partynode(8L, info = "no"))))),
  WeatherPlay)
plot(py)

## ----weather-splits-------------------------------------------------
sp_o <- partysplit(1L, index = 1:3)
sp_h <- partysplit(3L, breaks = 75)
sp_w <- partysplit(4L, index = 1:2)

## ----weather-nodes--------------------------------------------------
pn <- partynode(1L, split = sp_o, kids = list(
  partynode(2L, split = sp_h, kids = list(
    partynode(3L, info = "yes"),
    partynode(4L, info = "no"))),
  partynode(5L, info = "yes"),
  partynode(6L, split = sp_w, kids = list(
    partynode(7L, info = "yes"),
    partynode(8L, info = "no")))))

## ----weather-nodes-print--------------------------------------------
pn

## ----weather-party--------------------------------------------------
py <- party(pn, WeatherPlay)
print(py)

## ----weather-plot, eval=FALSE---------------------------------------
# plot(py)

## ----weather-predict------------------------------------------------
predict(py, head(WeatherPlay))

## ----weather-methods-dim--------------------------------------------
length(py)
width(py)
depth(py)

## ----weather-methods-subset-----------------------------------------
py[6]

## ----weather-methods-names------------------------------------------
py2 <- py
names(py2)
names(py2) <- LETTERS[1:8]
py2

## ----weather-methods-nodeids----------------------------------------
nodeids(py)
nodeids(py, terminal = TRUE)

## ----weather-methods-nodeapply--------------------------------------
nodeapply(py, ids = c(1, 7), FUN = function(n) n$info)
nodeapply(py, ids = nodeids(py, terminal = TRUE),
  FUN = function(n) paste("Play decision:", n$info))

## ----weather-methods-predict----------------------------------------
predict(py, FUN = function(n) paste("Play decision:", n$info))

## ----weather-methods-print------------------------------------------
print(py, terminal_panel = function(n)
  c(", then the play decision is:", toupper(n$info)))

## ----weather-methods-plot, eval=FALSE-------------------------------
# plot(py, tp_args = list(FUN = function(i)
#   c("Play decision:", toupper(i))))

## ----weather-methods-plot1, echo=FALSE, fig.height=5, fig.width=7.5----
plot(py[6])

## ----weather-methods-plot2, echo=FALSE, fig.height=5, fig.width=7.5----
plot(py, tp_args = list(FUN = function(i) 
  c("Play decision:", toupper(i))))

## ----weather-prune--------------------------------------------------
nodeprune(py, 2)
nodeprune(py, c(2, 6))

## ----partysplit-1, echo = TRUE--------------------------------------
sp_h <- partysplit(3L, breaks = 75)
class(sp_h)

## ----partysplit-2, echo = TRUE--------------------------------------
unclass(sp_h)

## ----partysplit-3, echo = TRUE--------------------------------------
character_split(sp_h, data = WeatherPlay)

## ----partysplit-4, echo = TRUE--------------------------------------
kidids_split(sp_h, data = WeatherPlay)

## ----partysplit-5, echo = TRUE--------------------------------------
as.numeric(!(WeatherPlay$humidity <= 75)) + 1

## ----partysplit-6, echo = TRUE--------------------------------------
sp_o2 <- partysplit(1L, index = c(1L, 1L, 2L))
character_split(sp_o2, data = WeatherPlay)
table(kidids_split(sp_o2, data = WeatherPlay), WeatherPlay$outlook)

## ----partysplit-6a, echo = TRUE-------------------------------------
unclass(sp_o2)

## ----partysplit-7, echo = TRUE--------------------------------------
sp_o <- partysplit(1L, index = 1L:3L)
character_split(sp_o, data = WeatherPlay)

## ----partysplit-8, echo = TRUE--------------------------------------
sp_t <- partysplit(2L, breaks = c(69.5, 78.8), index = c(1L, 2L, 1L))
character_split(sp_t, data = WeatherPlay)
table(kidids_split(sp_t, data = WeatherPlay),
  cut(WeatherPlay$temperature, breaks = c(-Inf, 69.5, 78.8, Inf)))

## ----partynode-1, echo = TRUE---------------------------------------
n1 <- partynode(id = 1L)
is.terminal(n1)
print(n1)

## ----partynode-2, echo = TRUE---------------------------------------
n1 <- partynode(id = 1L, split = sp_o, kids = lapply(2L:4L, partynode))
print(n1, data = WeatherPlay)

## ----partynode-3, echo = TRUE---------------------------------------
fitted_node(n1, data = WeatherPlay)

## ----partynode-4, echo = TRUE---------------------------------------
kidids_node(n1, data = WeatherPlay)

## ----party-1a, echo = TRUE------------------------------------------
t1 <- party(n1, data = WeatherPlay)
t1

## ----party-1b, echo = TRUE------------------------------------------
party(n1, data = WeatherPlay[0, ])

## ----party-2, echo = TRUE-------------------------------------------
t2 <- party(n1, 
  data = WeatherPlay,
  fitted = data.frame(
    "(fitted)" = fitted_node(n1, data = WeatherPlay),
    "(response)" = WeatherPlay$play,
    check.names = FALSE),
  terms = terms(play ~ ., data = WeatherPlay),
)

## ----party-3, echo=TRUE---------------------------------------------
t2 <- as.constparty(t2)
t2

## ----constparty-plot, echo=FALSE, fig.height=5, fig.width=6---------
plot(t2, tnex = 1.5)

## ----party-4, echo=TRUE---------------------------------------------
nd <- data.frame(outlook = factor(c("overcast", "sunny"),
  levels = levels(WeatherPlay$outlook)))
predict(t2, newdata = nd, type = "response")
predict(t2, newdata = nd, type = "prob")
predict(t2, newdata = nd, type = "node")

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

