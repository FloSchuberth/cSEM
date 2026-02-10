# TODO: See if I can replicate the analysis
# Some of this code is directly from https://osf.io/9tm2y/files/wk5vz?view_only=59d9792bab994892a735e4efc763b511
library(seminr)
library(naniar)

data(corp_rep_data, package = "seminr")

x_clean <-  naniar::replace_with_na_all(corp_rep_data, ~.x  == -99) |> 
  na.omit()

model <- '
Quality <~ qual_1 + qual_2 + qual_3 + qual_4 + qual_5 + qual_6 + qual_7 + qual_8
Performance <~ perf_1 + perf_2 + perf_3 + perf_4 + perf_5
CorpSocResp <~ csor_1 + csor_2 + csor_3 + csor_4 + csor_5
Attractiveness <~ attr_1 + attr_2 + attr_3

Competence =~ comp_1 + comp_2 + comp_3
CustomerLoyality =~ cusl_1 + cusl_2 + cusl_3
Likeability =~ like_1 + like_2 + like_3
CustomerSatisfaction =~ cusa

Competence ~ Quality + Performance + CorpSocResp + Attractiveness
Likeability ~ Quality + Performance + CorpSocResp + Attractiveness
CustomerSatisfaction ~ Competence + Likeability
CustomerLoyality ~ CustomerSatisfaction + Competence + Likeability'

mod <- csem(x_clean, model, .approach_weights = "GSCA", .disattenuate = TRUE,
.dominant_indicators = NULL,
.tolerance = 0.0001,
.conv_criterion = "sum_diff_absolute"
    )

tidy(mod)

# TODO: If I can or can't, then ask about how I can add it to the tests. This might be overkill/meaningless