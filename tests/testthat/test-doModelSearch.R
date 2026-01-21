# test_that("Not providing a cSEMResults object causes an error", {
#   expect_error(verify(list(1:3)))
# })

# .prob_mutation larger than 1
probs = list(1.5,-0.5,c(0.5,0.7),'a')

lapply(probs[4],function(x){
  doModelSearch(.object= out_str1, 
                .pop_size = 20,
                .n_generations =  20,
                .prob_mutation = x,
                .prob_crossover = 0.8,
                .fbar = -100000,
                .ms_criterion = 'bic')
})

lapply(probs[4],function(x){
  doModelSearch(.object= out_str1, 
                .pop_size = 20,
                .n_generations =  20,
                .prob_mutation = 0.5,
                .prob_crossover = x,
                .fbar = -100000,
                .ms_criterion = 'bic')
})

size=probs = list(1.5,-0.5,c(0.5,0.7),'a')


doModelSearch(.object= out_str1, 
              .pop_size = 20,
              .n_generations =  20,
              .prob_mutation = 1.5,
              .prob_crossover = 0.8,
              .fbar = -100000,
              .ms_criterion = 'bic')

# .prob_mutation smaller than 0
doModelSearch(.object= out_str1, 
              .pop_size = 20,
              .n_generations =  20,
              .prob_mutation = -.5,
              .prob_crossover = 0.8,
              .fbar = -100000,
              .ms_criterion = 'bic')

# .prob_mutation is a vector 
doModelSearch(.object= out_str1, 
              .pop_size = 20,
              .n_generations =  20,
              .prob_mutation = c(.5,0.7),
              .prob_crossover = 0.8,
              .fbar = -100000,
              .ms_criterion = 'bic')

doModelSearch(.object= out_str1, 
              .pop_size = 20,
              .n_generations =  20,
              .prob_mutation = c(.5,0.7),
              .prob_crossover = 0.8,
              .fbar = -100000,
              .ms_criterion = 'bic')
