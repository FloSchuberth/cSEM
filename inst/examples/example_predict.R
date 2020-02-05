\donttest{### Anime example taken from https://github.com/ISS-Analytics/pls-predict

# Load data
data(Anime) # data is similar to the Anime.csv found on 
            # https://github.com/ISS-Analytics/pls-predict but with irrelevant
            # columns removed

# Split into training and data the same way as it is done on 
# https://github.com/ISS-Analytics/pls-predict
set.seed(123)

index     <- sample.int(dim(Anime)[1], 83, replace = FALSE)
dat_train <- Anime[-index, ]
dat_test  <- Anime[index, ]

# Specify model
model <- "
# Structural model

ApproachAvoidance ~ PerceivedVisualComplexity + Arousal

# Measurement/composite model

ApproachAvoidance         =~ AA0 + AA1 + AA2 + AA3
PerceivedVisualComplexity <~ VX0 + VX1 + VX2 + VX3 + VX4
Arousal                   <~ Aro1 + Aro2 + Aro3 + Aro4
"

# Estimate (replicating the results of the `simplePLS()` function)
res <- csem(dat_train, 
            model, 
            .disattenuate = FALSE, # original PLS
            .iter_max = 300, 
            .tolerance = 1e-07, 
            .PLS_weight_scheme_inner = "factorial"
)

# Predict using a user-supplied training data set
pp <- predict(res, .test_data = dat_test)
pp$Predictions_target[1:6, ]
pp

### Compute prediction metrics  ------------------------------------------------
res2 <- csem(Anime, # whole data set
            model, 
            .disattenuate = FALSE, # original PLS
            .iter_max = 300, 
            .tolerance = 1e-07, 
            .PLS_weight_scheme_inner = "factorial"
)

# Predict using 10-fold cross-validation with 5 repetitions
pp2 <- predict(res, .benchmark = "lm")
pp2
## There is a plot method available
plot(pp2)}

