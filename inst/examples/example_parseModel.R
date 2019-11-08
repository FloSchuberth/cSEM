# ===========================================================================
# Providing a model in lavaan syntax 
# ===========================================================================
model <- "
# Structural model
y1 ~ y2 + y3

# Measurement model
y1 =~ x1 + x2 + x3
y2 =~ x4 + x5
y3 =~ x6 + x7

# Error correlation
x1 ~~ x2
"

m <- parseModel(model)
m

# ===========================================================================
# Providing a complete model in cSEM format (class cSEMModel)
# ===========================================================================
# If the model is already a cSEMModel object, the model is returned as is:

identical(m, parseModel(m)) # TRUE

# ===========================================================================
# Providing a list 
# ===========================================================================
# It is possible to provide a list that contains at least the
# elements "structural" and "measurement". This is generally discouraged
# as this may cause unexpected errors.

m_incomplete <- m[c("structural", "measurement", "construct_type")]
parseModel(m_incomplete)

# Providing a list containing list names that are not part of a `cSEMModel`
# causes an error:

\dontrun{
m_incomplete[c("name_a", "name_b")] <- c("hello world", "hello universe")
parseModel(m_incomplete)
}

# Failing to provide "structural" or "measurement" also causes an error:

\dontrun{
m_incomplete <- m[c("structural", "construct_type")]
parseModel(m_incomplete)
}