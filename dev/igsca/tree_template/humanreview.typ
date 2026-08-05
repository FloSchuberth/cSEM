#set document(
  title: "doTrees() Port Human Review",
  author: "Michael",
)
#set page(numbering: "1", margin: 2.4cm)
#set text(font: "New Computer Modern", size: 10.5pt)
#set heading(numbering: "1.1")
#show heading.where(level: 1): set text(size: 14pt, weight: "bold")
#show heading.where(level: 2): set text(size: 12pt, weight: "bold")


= On-Going Questions

+ Does it matter whether or not `coin_distribution = "asymptotic"` or `"approximate"`?
  - What is the default in ctree?
  - Do all the different options work as intended?
  - What about when I use the partition family path instead of conditional-test?
+ Why does the documentation seem catered to someone who already understands the codebase, as opposed to someone reading it for the first time?
+ Does igscaTree handle convergence failures in a similar way as semtree?
  - Within a resample?
  - Within a particular tested multigroup model?
  - Separation in convergence failure for variable selection versus split point selection
+ What is a kernel?
+ What helper functions or ways of interacting with the fitted tree need to be added/documented?
  - Split-point + p-value?
  - 


= helper_doTrees.R

== igsca_tree_control

+ How much does this overlap or contradict with `partykit::ctree_control()`? Is the functionality consistent? The documentation?
+ Do the default arguments make sense?

= test-postestimate_doTrees.R

+ What are reasonable defaults for `ctl_mixed` and `ctl_part`?

= postestimate_doTrees.R

== doTrees

+ What do each of the attributes of igsca_info mean?
+ What does the 'native' splitter using coin really mean?