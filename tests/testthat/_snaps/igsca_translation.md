# Matlab_R_Compared_Names

    Code
      igsca_sim(z0 = igsca_sim_in$z0, W0 = igsca_sim_in$W0, C0 = igsca_sim_in$C0, B0 = igsca_sim_in$
        B0, lv_type = igsca_sim_in$lv_type, ov_type = igsca_sim_in$ov_type, ind_domi = igsca_sim_in$
        ind_domi, nbt = 0, devmode = TRUE, swap_step = "noswap", devmode_checkobj = TRUE)
    Message
      As manually checked on December 17/2023, I'm fairly sure that the non-overlapping names are OK. Most of the Matlab 'unique' names actually show-up in the R sub-functions. A small minority simply don't appear for parsimony (C_t, crindex) and some probably aren't legal/good-idea R names (t). The R-unique object names are intermediaries from functions
    Output
      $RObj_NotIn_Matlab
      [1] "devdir"            "devmode"           "devmode_checkobj" 
      [4] "flipped_signs"     "prepared_for_ALS"  "swap_step"        
      [7] "tot"               "updated_C_B_D"     "updated_X_weights"
      
      $MatObj_NotIn_R
       [1] "C_t"       "Delta"     "H1"        "H2"        "M1"        "M2"       
       [7] "XI"        "beta"      "bz0"       "crindex"   "e"         "in_dir"   
      [13] "out_dir"   "t"         "t1"        "t2"        "vecZDelta"
      

