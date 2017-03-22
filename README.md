# DriftGenerator
Generates data with drift

# Covariate
* Covariate distribution is initialised with a bayesian network.
* Covariate is drifted by linear interpolation of the conditional probability tables.

# Posterior
* Each unique covariate value combination is assigned a class via a tree with the following structure:
>
                            [x_1,x_2]
                                +
                                |
                                |
                                |
                                |
                       +----------------+
                       |        |       |
                  x_1=0|   x_1=1|  x_1=2|
                       v        v       v
                      y_1       |       |
                                |       |
              +-----------------+       +-----------------+
              |        |        |       |        |        |
         x_2=0|   x_2=1|   x_2=2|       |x_2=0   |x_2=1   |x_2=2
              v        v        v       v        v        v
             y_1      y_1      y_1     y_1      y_1      y_1
* Note that the first value of each attribute layer is always assigned a class
* The Posterior is drifted by changing the class assigned at the leaf nodes
