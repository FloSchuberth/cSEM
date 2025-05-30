% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/csem_resample.R
\name{resampleData}
\alias{resampleData}
\title{Resample data}
\usage{
resampleData(
 .object          = NULL,
 .data            = NULL,
 .resample_method = c("bootstrap", "jackknife", "permutation", 
                      "cross-validation"),
 .cv_folds        = 10,
 .id              = NULL,
 .R               = 499,
 .seed            = NULL
)
}
\arguments{
\item{.object}{An R object of class \link{cSEMResults} resulting from a call to \code{\link[=csem]{csem()}}.}

\item{.data}{A \code{data.frame}, a \code{matrix} or a \code{list} of data of either type.
Possible column types or classes of the data provided are:
"\code{logical}", "\code{numeric}" ("\code{double}" or "\code{integer}"), "\code{factor}" (ordered and unordered)
or a mix of several types. The data may also include
\strong{one} character column whose column name must be given to \code{.id}.
This column is assumed to contain group identifiers used to split the data into groups.
If \code{.data} is provided, \code{.object} is ignored. Defaults to \code{NULL}.}

\item{.resample_method}{Character string. The resampling method to use. One of:
"\emph{bootstrap}", "\emph{jackknife}", "\emph{permutation}", or "\emph{cross-validation}".
Defaults to "\emph{bootstrap}".}

\item{.cv_folds}{Integer. The number of cross-validation folds to use. Setting
\code{.cv_folds} to \code{N} (the number of observations) produces
leave-one-out cross-validation samples. Defaults to \code{10}.}

\item{.id}{Character string or integer. A character string giving the name or
an integer of the position of the column of \code{.data} whose levels are used
to split \code{.data} into groups. Defaults to \code{NULL}.}

\item{.R}{Integer. The number of bootstrap runs, permutation runs
or cross-validation repetitions to use. Defaults to \code{499}.}

\item{.seed}{Integer or \code{NULL}. The random seed to use. Defaults to \code{NULL} in which
case an arbitrary seed is chosen. Note that the scope of the seed is limited
to the body of the function it is used in. Hence, the global seed will
not be altered!}
}
\value{
The structure of the output depends on the type of input and the
resampling method:
\describe{
\item{Bootstrap}{If a \code{matrix} or \code{data.frame} without grouping variable
is provided (i.e., \code{.id = NULL}), the result is a list of length \code{.R}
(default \code{499}). Each element of that list is a bootstrap (re)sample.
If a grouping variable is specified or a list of data is provided
(where each list element is assumed to contain data for one group),
resampling is done by group. Hence,
the result is a list of length equal to the number of groups
with each list element containing \code{.R} bootstrap samples based on the
\code{N_g} observations of group \code{g}.}
\item{Jackknife}{If a \code{matrix} or \code{data.frame} without grouping variable
is provided (\code{.id = NULL}), the result is a list of length equal to the number
of observations/rows (\code{N}) of the data set provided.
Each element of that list is a jackknife (re)sample.
If a grouping variable is specified or a list of data is provided
(where each list element is assumed to contain data for one group),
resampling is done by group. Hence,
the result is a list of length equal to the number of group levels
with each list element containing \code{N} jackknife samples based on the
\code{N_g} observations of group \code{g}.}
\item{Permutation}{If a \code{matrix} or \code{data.frame} without grouping variable
is provided an error is returned as permutation will simply reorder the observations.
If a grouping variable is specified or a list of data is provided
(where each list element is assumed to contain data of one group),
group membership is permuted. Hence, the result is a list of length \code{.R}
where each element of that list is a permutation (re)sample.}
\item{Cross-validation}{If a \code{matrix} or \code{data.frame} without grouping variable
is provided a list of length \code{.R} is returned. Each list element
contains a list containing the \code{k} splits/folds subsequently
used as test and training data sets.
If a grouping variable is specified or a list of data is provided
(where each list element is assumed to contain data for one group),
cross-validation is repeated \code{.R} times for each group. Hence,
the result is a list of length equal to the number of groups,
each containing \code{.R} list elements (the repetitions) which in turn contain
the \code{k} splits/folds.
}
}
}
\description{
Resample data from a data set using common resampling methods.
For bootstrap or jackknife resampling, package users usually do not need to
call this function but directly use \code{\link[=resamplecSEMResults]{resamplecSEMResults()}} instead.
}
\details{
The function \code{resampleData()} is general purpose. It simply resamples data
from a data set according to the resampling method provided
via the \code{.resample_method} argument and returns a list of resamples.
Currently, \code{bootstrap}, \code{jackknife}, \code{permutation}, and  \code{cross-validation}
(both leave-one-out (LOOCV) and k-fold cross-validation) are implemented.

The user may provide the data set to resample either explicitly via the \code{.data}
argument or implicitly by providing a \link{cSEMResults} objects to \code{.object}
in which case the original data used in the call that created the
\link{cSEMResults} object is used for resampling.
If both, a \link{cSEMResults} object and a data set via \code{.data} are provided
the former is ignored.

As \code{\link[=csem]{csem()}} accepts a single data set, a list of data sets as well as data sets
that contain a column name used to split the data into groups,
the \link{cSEMResults} object may contain multiple data sets.
In this case, resampling is done by data set or group. Note that depending
on the number of data sets/groups provided this computation may be slower
as resampling will be repeated for each data set/group.

To split data provided via the \code{.data} argument into groups, the column name or
the column index of the column containing the group levels to split the data
must be given to \code{.id}. If data that contains grouping is taken from
a \link{cSEMResults} object, \code{.id} is taken from the object information. Hence,
providing  \code{.id} is redundant in this case and therefore ignored.

The number of bootstrap or permutation runs as well as the number of
cross-validation repetitions is given by \code{.R}. The default is
\code{499} but should be increased in real applications. See e.g.,
\insertCite{Hesterberg2015;textual}{cSEM}, p.380 for recommendations concerning
the bootstrap. For jackknife \code{.R} is ignored as it is based on the N leave-one-out data sets.

Choosing \code{resample_method = "permutation"} for ungrouped data causes an error
as permutation will simply reorder the observations which is usually not
meaningful. If a list of data is provided
each list element is assumed to represent the observations belonging to one
group. In this case, data is pooled and group adherence permuted.

For cross-validation the number of folds (\code{k}) defaults to \code{10}. It may be
changed via the \code{.cv_folds} argument. Setting \code{k = 2} (not 1!) splits
the data into a single training and test data set. Setting \code{k = N} (where \code{N} is the
number of observations) produces leave-one-out cross-validation samples.
Note: 1.) At least 2 folds required  (\code{k > 1}); 2.) \code{k} can not be larger than \code{N};
3.) If \code{N/k} is not not an integer the last fold will have less observations.

Random number generation (RNG) uses the L'Ecuyer-CRMR RGN stream as implemented
in the \href{https://github.com/futureverse/future.apply/}{future.apply package}
\insertCite{Bengtsson2018a}{cSEM}.
See \link[future.apply:future_lapply]{?future_lapply} for details. By default
a random seed is chosen.
}
\examples{
# ===========================================================================
# Using the raw data 
# ===========================================================================
### Bootstrap (default) -----------------------------------------------------

res_boot1 <- resampleData(.data = satisfaction)
str(res_boot1, max.level = 3, list.len = 3)

## To replicate a bootstrap draw use .seed:
res_boot1a <- resampleData(.data = satisfaction, .seed = 2364)
res_boot1b <- resampleData(.data = satisfaction, .seed = 2364)
                           
identical(res_boot1, res_boot1a) # TRUE

### Jackknife ---------------------------------------------------------------

res_jack <- resampleData(.data = satisfaction, .resample_method = "jackknife")
str(res_jack, max.level = 3, list.len = 3)

### Cross-validation --------------------------------------------------------
## Create dataset for illustration:
dat <- data.frame(
  "x1" = rnorm(100),
  "x2" = rnorm(100),
  "group" = sample(c("male", "female"), size = 100, replace = TRUE),
  stringsAsFactors = FALSE)

## 10-fold cross-validation (repeated 100 times)
cv_10a <- resampleData(.data = dat, .resample_method = "cross-validation", 
                      .R = 100)
str(cv_10a, max.level = 3, list.len = 3)

# Cross-validation can be done by group if a group identifyer is provided:
cv_10 <- resampleData(.data = dat, .resample_method = "cross-validation", 
                      .id = "group", .R = 100)

## Leave-one-out-cross-validation (repeated 50 times)
cv_loocv  <- resampleData(.data = dat[, -3], 
                          .resample_method = "cross-validation", 
                          .cv_folds = nrow(dat),
                          .R = 50)
str(cv_loocv, max.level = 2, list.len = 3)

### Permuation ---------------------------------------------------------------

res_perm <- resampleData(.data = dat, .resample_method = "permutation",
                         .id = "group")
str(res_perm, max.level = 2, list.len = 3)

# Forgetting to set .id causes an error
\dontrun{
res_perm <- resampleData(.data = dat, .resample_method = "permutation")
}

# ===========================================================================
# Using a cSEMResults object
# ===========================================================================

model <- "
# Structural model
QUAL ~ EXPE
EXPE ~ IMAG
SAT  ~ IMAG + EXPE + QUAL + VAL
LOY  ~ IMAG + SAT
VAL  ~ EXPE + QUAL

# Measurement model
EXPE =~ expe1 + expe2 + expe3 + expe4 + expe5
IMAG =~ imag1 + imag2 + imag3 + imag4 + imag5
LOY  =~ loy1  + loy2  + loy3  + loy4
QUAL =~ qual1 + qual2 + qual3 + qual4 + qual5
SAT  =~ sat1  + sat2  + sat3  + sat4
VAL  =~ val1  + val2  + val3  + val4
"
a <- csem(satisfaction, model)

# Create bootstrap and jackknife samples
res_boot <- resampleData(a, .resample_method = "bootstrap", .R = 499)
res_jack <- resampleData(a, .resample_method = "jackknife")

# Since `satisfaction` is the dataset used the following approaches yield
# identical results.
res_boot_data   <- resampleData(.data = satisfaction, .seed = 2364)
res_boot_object <- resampleData(a, .seed = 2364)

identical(res_boot_data, res_boot_object) # TRUE
}
\references{
\insertAllCited{}
}
\seealso{
\code{\link[=csem]{csem()}}, \link{cSEMResults}, \code{\link[=resamplecSEMResults]{resamplecSEMResults()}}
}
