# GSCA -------------------------------------------------------------------
gsca_model <- "
# Composite Model
NetworkingBehavior <~ Behavior1 + Behavior2 + Behavior3
HonestyHumility <~ Honesty1 + Honesty2 + Honesty3

# Structural Model
NetworkingBehavior ~ HonestyHumility
"

gsca_CCA_model <- "
# Composite Model
NetworkingBehavior <~ Behavior1 + Behavior2 + Behavior3
HonestyHumility <~ Honesty1 + Honesty2 + Honesty3

# Structural Model
NetworkingBehavior ~~ HonestyHumility
"

data(LeDang2022)

dat <- rbind(
  subset(LeDang2022, Gender == "Male")[1:10,],
  subset(LeDang2022, Gender == "Female")[1:10,]
)

debugonce(calculateWeightsGSCA)
gsca <- csem(
  .data = dat,
  gsca_model,
  # gsca_CCA_model,
  .approach_weights = "GSCA",
  .dominant_indicators = NULL,
  .tolerance = 0.0001,
  .conv_criterion = "sum_diff_absolute",
  # .GSCA_modes = "CCMP"
  .GSCA_modes = "NCMP"
)

gsca_mg <- csem(
  .data = dat,
  gsca_model,
  .approach_weights = "GSCA",
  .dominant_indicators = NULL,
  .tolerance = 0.0001,
  .conv_criterion = "sum_diff_absolute",
  .id = "Gender"
)

# debugonce(calculateFIT)
calculateFIT(gsca)
# debugonce(calculateFIT)
calculateFIT(gsca_mg)

# debugonce(calculateAFIT)
calculateAFIT(gsca)
# debugonce(calculateAFIT)
calculateAFIT(gsca_mg)

# debugonce(calculateFIT_m)
calculateFIT_m(gsca)
# debugonce(calculateFIT_m)
calculateFIT_m(gsca_mg)

# debugonce(calculateFIT_s)
calculateFIT_s(gsca)
# debugonce(calculateFIT_s)
calculateFIT_s(gsca_mg)

# GSCAM ------------------------------------------------------------------
gscam_model <- "
# Latent Varible Model
NetworkingBehavior =~ Behavior1 + Behavior2 + Behavior3
HonestyHumility =~ Honesty1 + Honesty2 + Honesty3

# Structural Model
NetworkingBehavior ~ HonestyHumility
"

gscam_CFA_model <- "
# Latent Varible Model
NetworkingBehavior =~ Behavior1 + Behavior2 + Behavior3
HonestyHumility =~ Honesty1 + Honesty2 + Honesty3

# Correlated Factors
NetworkingBehavior ~~ HonestyHumility
"

data(LeDang2022)

dat <- rbind(
  subset(LeDang2022, Gender == "Male")[1:10,],
  subset(LeDang2022, Gender == "Female")[1:10,]
)

debugonce(calculateWeightsGSCAm)
gscam <- csem(
  .data = dat,
  gscam_model,
  # gscam_CFA_model,
  .approach_weights = "GSCA",
  .dominant_indicators = NULL,
  .tolerance = 0.0001,
  .conv_criterion = "sum_diff_absolute"
)

gscam_mg <- csem(
  .data = dat,
  gscam_model,
  .approach_weights = "GSCA",
  .dominant_indicators = NULL,
  .tolerance = 0.0001,
  .conv_criterion = "sum_diff_absolute",
  .id = "Gender"
)

# debugonce(calculateFIT)
calculateFIT(gscam)
# debugonce(calculateFIT)
# debugonce(bdiagGSCA)
calculateFIT(gscam_mg)

# debugonce(calculateAFIT)
calculateAFIT(gscam)
# debugonce(calculateAFIT)
calculateAFIT(gscam_mg)

# debugonce(calculateFIT_m)
calculateFIT_m(gscam)
# debugonce(calculateFIT_m)
calculateFIT_m(gscam_mg)

# debugonce(calculateFIT_s)
calculateFIT_s(gscam)
calculateFIT_s(gscam_mg)

# IGSCA ------------------------------------------------------------------
igsca_model <- "
# Composite Model
NetworkingBehavior <~ Behavior1 + Behavior2 + Behavior3

# Latent Variable Model
HonestyHumility =~ Honesty1 + Honesty2 + Honesty3

# Structural Model
NetworkingBehavior ~ HonestyHumility
"

igsca_model_nostruct <- "
# Composite Model
NetworkingBehavior <~ Behavior1 + Behavior2 + Behavior3

# Latent Variable Model
HonestyHumility =~ Honesty1 + Honesty2 + Honesty3

# CCA
NetworkingBehavior ~~ HonestyHumility
"

debugonce(calculateWeightsIGSCA)
igsca <- csem(
  .data = dat,
  igsca_model,
  # igsca_model_nostruct,
  .approach_weights = "GSCA",
  .dominant_indicators = NULL,
  .tolerance = 0.0001,
  .conv_criterion = "sum_diff_absolute",
  .GSCA_modes = "CCMP"
  # .GSCA_modes = "NCMP"
)

igsca_mg <- csem(
  .data = dat,
  igsca_model,
  .approach_weights = "GSCA",
  .dominant_indicators = NULL,
  .tolerance = 0.0001,
  .conv_criterion = "sum_diff_absolute",
  .id = "Gender"
)

debugonce(calculateFIT)
calculateFIT(igsca)
debugonce(calculateFIT)
calculateFIT(igsca_mg)

# debugonce(calculateAFIT)
calculateAFIT(igsca)
# debugonce(calculateAFIT)
calculateAFIT(igsca_mg)

debugonce(calculateFIT_m)
calculateFIT_m(igsca)
calculateFIT_m(igsca_mg)

debugonce(calculateFIT_s)
calculateFIT_s(igsca)
calculateFIT_s(igsca_mg)



igsca_model_nostruct <- "
# Composite Model
NetworkingBehavior <~ Behavior1 + Behavior2 + Behavior3

# Latent Variable Model
HonestyHumility =~ Honesty1 + Honesty2 + Honesty3
NetworkingBehavior ~~ HonestyHumility
"

data(LeDang2022)

dat <- rbind(
  subset(LeDang2022, Gender == "Male")[1:10,],
  subset(LeDang2022, Gender == "Female")[1:10,]
)

# debugonce(calculateWeightsIGSCA)
# debug(calculateWeightsGSCA)
igsca_nostruct <- csem(
  .data = dat,
  igsca_model_nostruct,
  .approach_weights = "GSCA",
  .dominant_indicators = NULL,
  .tolerance = 0.0001,
  .conv_criterion = "sum_diff_absolute",
  .GSCA_modes = "NCMP"
)

igsca_mg_nostruct <- csem(
  .data = dat,
  igsca_model_nostruct,
  .approach_weights = "GSCA",
  .dominant_indicators = NULL,
  .tolerance = 0.0001,
  .conv_criterion = "sum_diff_absolute",
  .id = "Gender",
  .GSCA_modes = "NCMP"
)

# debugonce(calculateFIT)
calculateFIT(igsca_nostruct)
# debugonce(calculateFIT)
calculateFIT(igsca_mg_nostruct)

# debugonce(calculateFIT_s)
calculateFIT_s(igsca_nostruct)
# debugonce(calculateFIT_s)
calculateFIT_s(igsca_mg_nostruct)
