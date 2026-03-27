## System details:

#### Positron and OS details:

```
Positron Version: 2026.03.0 build 212
Code - OSS Version: 1.108.0
Commit: f3aae65e0a1a11d39226cd884520f49301daef82
Date: 2026-02-26T05:01:47.944Z
Electron: 39.2.7
Chromium: 142.0.7444.235
Node.js: 22.21.1
V8: 14.2.231.21-electron.0
OS: Linux x64 6.8.0-106-generic
```

#### Session details:

```
R 4.5.1
```

## Describe the issue:

Printing a `profvis` object (i.e., viewing profiling results) opens a new editor tab, but the tab frequently remains blank indefinitely. The behavior is unreliable — sometimes the profiling visualization renders correctly, but it is often slow and sometimes never loads at all. The tab stays open with no content and no error indication.

This issue does not occur in RStudio, where `profvis` output renders quickly and reliably every time.

Related: profvis tab support was added in #11067.

## Steps to reproduce the issue:

1. Open an R script in Positron.
2. Run the following code (adapted from the [profvis examples](https://profvis.r-lib.org/articles/examples.html)):

```r
library(profvis)

times <- 10
cols <- 150
data <- as.data.frame(x = matrix(rnorm(times * cols, mean = 5), ncol = cols))
data <- cbind(id = paste0("g", seq_len(times)), data)

test <- profvis({
  data1 <- data

  means <- apply(data1[, names(data1) != "id"], 2, mean)

  for (i in seq_along(means)) {
    data1[, names(data1) != "id"][, i] <- data1[, names(data1) != "id"][, i] -
      means[i]
  }
})

test
```

3. Observe that a new editor tab opens. The tab frequently remains blank with no profiling visualization rendered.
4. Repeat the above steps multiple times — the tab sometimes renders correctly but is often blank or very slow to load.

## Expected or desired behavior:

The `profvis` visualization should render reliably and promptly in the editor tab every time, consistent with the behavior in RStudio.

## Were there any error messages in the UI, Output panel, or Developer Tools console?

No error messages are displayed. The tab opens without errors — it simply remains blank with no content rendered.

---

_This issue was written with the assistance of Claude Opus 4.6._
