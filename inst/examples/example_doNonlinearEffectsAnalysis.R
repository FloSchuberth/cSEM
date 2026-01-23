\dontrun{
model_Int <- "
# Measurement models
INV =~ INV1 + INV2 + INV3 +INV4
SAT =~ SAT1 + SAT2 + SAT3
INT =~ INT1 + INT2

# Structrual model containing an interaction term.
INT ~ INV + SAT + INV.SAT
"
  
# Estimate model
out <- csem(.data = Switching, .model = model_Int,
            # ADANCO settings
            .PLS_weight_scheme_inner = 'factorial',
            .tolerance = 1e-06,
            .resample_method = 'bootstrap'
)
  
# Do nonlinear effects analysis
neffects <- doNonlinearEffectsAnalysis(out, 
                                       .dependent = 'INT',
                                       .moderator = 'INV',
                                       .independent = 'SAT') 

# Get an overview
neffects

# Simple effects plot
plot(neffects, .plot_type = 'simpleeffects')

# Surface plot using plotly
plot(neffects, .plot_type = 'surface', .plot_package = 'plotly')

# Surface plot using persp
plot(neffects, .plot_type = 'surface', .plot_package = 'persp')

# Floodlight analysis
plot(neffects, .plot_type = 'floodlight')
}