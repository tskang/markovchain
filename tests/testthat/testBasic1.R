library(markovchain)

#create basic markov chains
markov1<-new("markovchain", states=c("a","b","c"), transitionMatrix=
               matrix(c(0.2,0.5,0.3,
                        0,1,0,
                        0.1,0.8,0.1),nrow=3, byrow=TRUE, dimnames=list(c("a","b","c"),
                                                                       c("a","b","c"))
               ))

require(matlab)
mathematicaMatr <- zeros(5)
mathematicaMatr[1,] <- c(0, 1/3, 0, 2/3, 0)
mathematicaMatr[2,] <- c(1/2, 0, 0, 0, 1/2)
mathematicaMatr[3,] <- c(0, 0, 1/2, 1/2, 0)
mathematicaMatr[4,] <- c(0, 0, 1/2, 1/2, 0)
mathematicaMatr[5,] <- c(0, 0, 0, 0, 1)
statesNames <- letters[1:5]
mathematicaMc <- new("markovchain", transitionMatrix = mathematicaMatr,
                     name = "Mathematica MC", states = statesNames)
####end creating DTMC
context("Basic DTMC proprieties")

test_that("States are those that should be", {
  expect_equal(absorbingStates(markov1), "b")
  expect_equal(transientStates(markov1), c("a","c"))
  expect_equal(is.irreducible(mathematicaMc),FALSE)
  expect_equal(transientStates(mathematicaMc), c("a","b"))
  expect_equal(is.accessible(mathematicaMc, "a", "c"),TRUE)
#   expect_equal(firstPassage(markov1, "b", 2), matrix(c(0,1,0,0,0,0), nrow=2, byrow=TRUE, 
#                                                      dimnames=list(c("1","2"), c("a","b","c"))))
  expect_equal(summary(mathematicaMc), list(closedClasses = list(c("c", "d"), c("e")), 
                                            transientClasses = list(c("a", "b"))))
})

#first passage time
#@TAE: check this: 
mathematicaBMatr= matrix(c(  0, 1/2, 1/2,  1/2, 0, 1/2,  1/2, 1/2, 0),byrow=TRUE, nrow=3) ;
mathematicabMc<-as(mathematicaBMatr, "markovchain")
test_that("First passgae time should be", {
  expect_equal(firstPassage(mathematicabMc, "s3",3), matrix(c(0.5, 0.5, 0,
                                                            0.25, 0.25, 0.5,
                                                            0.125, 0.125, 0.25), nrow=3, byrow=TRUE,
                                                            dimnames=list(c("1", "2", "3"), c("s1", "s2", "s3"))
                                                            ))
  #if you repeat this more thime results change. No good
  #should behave like
  #http://www.wolfram.com/mathematica/new-in-9/markov-chains-and-queues/distribution-of-times-to-reach-a-target-state.html
})

###testing proper conversion of objeects
context("Conversion of objects")
provaMatr2Mc<-as(mathematicaMatr,"markovchain")

test_that("Conversion of objects", 
          {
            expect_equal(class(provaMatr2Mc)=="markovchain",TRUE)
          })

