# debugonce(grow_tree)
# debug(doTrees)

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

