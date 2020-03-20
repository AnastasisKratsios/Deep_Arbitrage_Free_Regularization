### This SGD Code is Due to:
# Dennis Prangle
# NewCastle University

# Anastasis Kratsios Would like to thank them for this implementation! :)


# Reference class objects for gradient descent algorithms.
# The "step" method take an estimated gradient as argument
# and returns an increment to apply to the current parameter values.
# See end of file for an example.

##Algorithm which does no adaptation
noAdapt = setRefClass("noAdapt",
                      fields = list(),
                      methods = list(
                        step = function(gradientEst) {
                          return(0*gradientEst)
                        }
                      )
)

##Plain SGD
SGD = setRefClass("SGD",
                  fields = list(initialLearningRate="numeric", stepCount="numeric"),
                  methods = list(
                    step = function(gradientEst) {
                      inc = - gradientEst * initialLearningRate / stepCount
                      stepCount <<- stepCount + 1
                      return(inc)
                    },
                    initialize = function(initialLearningRate, stepCount=1) {
                      .self$initialLearningRate = initialLearningRate
                      .self$stepCount = stepCount
                    }
                  )
)

##Momentum
momentum = setRefClass("momentum",
                       fields = list(initialLearningRate="numeric", alpha="numeric", velocity="numeric", stepCount="numeric"),
                       methods = list(
                         step = function(gradientEst) {
                           velocity <<- velocity * alpha - gradientEst * initialLearningRate / stepCount
                           inc = velocity
                           stepCount <<- stepCount + 1
                           return(inc)
                         },
                         initialize = function(initialLearningRate, alpha=0.5, initialVelocity=0, stepCount=1) {
                           .self$initialLearningRate = initialLearningRate
                           .self$alpha = alpha
                           .self$velocity = initialVelocity
                           .self$stepCount = stepCount
                         }
                       )
)

##AdaGrad
adaGrad = setRefClass("adaGrad",
                      fields = list(learningRate="numeric", delta="numeric", alpha="numeric", r="numeric", velocity="numeric"),
                      methods = list(
                        step = function(gradientEst) {
                          r <<- r + c(gradientEst)^2
                          velocity <<- velocity * alpha - c(gradientEst) * learningRate / (delta + sqrt(r))
                          return(velocity)
                        },
                        initialize = function(learningRate, alpha=0, delta=1E-7) {
                          .self$learningRate = learningRate
                          .self$delta = delta
                          .self$alpha = alpha
                          .self$r = 0 ##nb r and velocity start off as scalars but become vectors of the right length on their first update
                          .self$velocity = 0
                        }
                      )
)

##adam
adam = setRefClass("adam",
                   fields = list(epsilon="numeric", rho1="numeric", rho2="numeric", delta="numeric", r="numeric", s="numeric", t="numeric"),
                   methods = list(
                     step = function(gradientEst) {
                       t <<- t + 1
                       s <<- rho1*s + (1-rho1)*c(gradientEst)
                       r <<- rho2*r + (1-rho2)*(c(gradientEst)^2)
                       s_hat = s / (1-rho1^t)
                       r_hat = r / (1-rho2^t)
                       inc = - epsilon * s_hat / (sqrt(r_hat) + delta)
                       return(inc)
                     },
                     initialize = function(epsilon=0.001, rho1=0.9, rho2=0.999, delta=1E-8) {
                       .self$epsilon = epsilon
                       .self$rho1 = rho1
                       .self$rho2 = rho2
                       .self$delta = delta
                       .self$r = 0 ##nb r and s start off as scalars but become vectors of the right length on their first update
                       .self$s = 0
                       .self$t = 0
                     }
                   )
)

##AdaDelta
adaDelta = setRefClass("adaDelta",
                       fields = list(epsilon="numeric", rho="numeric", acc_grad="numeric", acc_up="numeric"),
                       methods = list(
                         step = function(gradientEst) {
                           acc_grad <<- rho*acc_grad + (1-rho)*c(gradientEst)^2
                           inc = -gradientEst*sqrt(acc_up+epsilon)/sqrt(acc_grad+epsilon)
                           acc_up <<- rho*acc_up + (1-rho)*inc^2
                           return(inc)
                         },
                         initialize = function(epsilon=1E-6, rho=0.95) { ##Default values taken from Zeiler (2012) section 4.1
                           .self$epsilon = epsilon
                           .self$rho = rho
                           .self$acc_grad = 0 ##nb acc_grad and acc_up start off as scalars but become vectors of the right length on their first update
                           .self$acc_up = 0
                         }
                       )
)

##RMSprop
rmsProp = setRefClass("rmsProp",
                      fields = list(learningRate="numeric", rho="numeric", epsilon="numeric", acc_grad="numeric"),
                      methods = list(
                        step = function(gradientEst) {
                          acc_grad <<- rho*acc_grad + (1-rho)*c(gradientEst)^2
                          inc = -gradientEst*learningRate / sqrt(acc_grad+epsilon)
                          return(inc)
                        },
                        initialize = function(learningRate=1E-3, rho=0.9, epsilon=1E-6) { ##Defaults taken from a mix of places
                          .self$learningRate = learningRate
                          .self$rho = rho
                          .self$epsilon = epsilon
                          .self$acc_grad = 0 ##nb acc_grad starts off as a scalar but becomes a vector of the right length on its first update
                        }
                      )
)


# Many thanks to :Dennis Prangle of NewCastle University for this wonderful implementation of the SGD algorithm in R!