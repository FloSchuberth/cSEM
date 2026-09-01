# What I am focusing on
# grow_tree(influence = "mat", splitter = "native", control = igsca_tree_control())

airq <- subset(airquality, !is.na(Ozone))
airct <- ctree(Ozone ~ ., data = airq) # TODO: What does thes influence function and trafo look like for this?
airct_list <- unclass(airct)
airct_list$trafo
# debug(ctree)
debug(partykit::ctree_control()$splitfun)
# boomer::rig_in_place(partykit:::.split)
debug(partykit:::.split)
# When you're inside .split run
  debug(FUN)
    debug(.ctree_test)
      debug(.ctree_test_1d)
        debug(.ctree_test_internal)
        # Passes the influence to LinStatExpCov, which goes to doTest, which gives you an index.  
        debug(doTest)

# debugonce(grow_tree)
debug(doTrees)

# debug(partykit::ctree)

# debug(partykit::extree_fit)
# debug(partykit:::.extree_node)

# debug(partykit:::.ctree_select)
debug(partykit:::.select)
  # debug(partykit:::.ctree_test)
  #   debug(partykit:::.ctree_test_2d) # If iy is not null, but when is that the case?
  #   debug(partykit:::.ctree_test_1d)
  #     debug(partykit:::.ctree_test_internal)
  #       debug(LinStatExpCov) # From libcoin package
  #       debug(doTest) # From libcoin package
          debug(R_QuadraticTest) # C code from libcoin package

# debug(partykit:::.ctree_split)
debug(partykit:::.split)
  debug(partykit:::.ctree_test)
    debug(partykit:::.ctree_test_1d)
      debug(partykit:::.ctree_test_internal)
# debug(partykit::party)


# For trying to understand how the split-point selection works
debug(partykit:::.extree_node)
  debug(argmax_split)
    debug(partition_stat)
      debug(split_max_fitdiff)
      debug(split_max_dli)
      debug(split_max_dgi)

# For how the list of different possible splitting partitions is constructed
  debug(node_group_data)
  debug(candidate_partitions)
