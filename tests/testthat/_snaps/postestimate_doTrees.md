# IGSCA Trees Conditional Test on Matrix of Residuals Runs as expected

    Code
      trees_mx
    Output
      
      Model formula:
      ~x31 + x32 + x33 + x11 + x12 + x13 + x41 + x42 + x43 + x44 + 
          x21 + x22 + x23 + (z_true + noise_1 + noise_2)
      
      Fitted party:
      [1] root
      |   [2] z_true in 1
      |   |   [3] noise_1 <= -1.09783: *
      |   |   [4] noise_1 > -1.09783: *
      |   [5] z_true in 2: *
      
      Number of inner nodes:    2
      Number of terminal nodes: 3

# IGSCA Trees Conditional Test on Vector of Residuals Runs as expected

    Code
      trees_vec
    Output
      
      Model formula:
      ~x31 + x32 + x33 + x11 + x12 + x13 + x41 + x42 + x43 + x44 + 
          x21 + x22 + x23 + (z_true + noise_1 + noise_2)
      
      Fitted party:
      [1] root
      |   [2] z_true in 1: *
      |   [3] z_true in 2: *
      
      Number of inner nodes:    1
      Number of terminal nodes: 2

