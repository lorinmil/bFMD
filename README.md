# bFMD
Balanced Functional Module Detection

## Preliminaries
This program takes in a high-dimensional dataset, X, and a continuous outcome, y, and applies the balanced Functional Module Detection (bFMD) method. Essentially, the purpose of bFMD is to identify a subset of variables within X that are jointly correlated with each other AND related to the outcome y. These subset of variables indicate nodes within a semidirected balanced graph. Balance suggests that the directionality between module variables and the outcome do not contradict each other in order to improve biological interpretability (e.g. if variables x1 and x2 are positively correlated with each other, then sign(cor(x1, y)) should equal sign(cor(x2, y)). Details of this method may be found in

    Tritchler, David, Lorin M. Towle-Miller, and Jeffrey C. Miecznikowski. "Balanced Functional Module Detection in Genomic Data." bioRxiv (2020).

## STEP 0: Simulations (optional)
If data is needed to explore the bFMD method, the simDatHub function found in sim1typeData.R may be executed. This will return a dataset X and a vector y.

    #Import all functions for simulating data with a functional module
    source("sim1typeData.R")
    #Perform the simulation (refer to sim1typeData.R for details on parameters)
    dat <- simDatHub()
    #Pull out the simulated dataset and outcome
    X <- dat$x
    y <- dat$y

## STEP 1: Perform bFMD
To execute the bFMD method, use the FMD function from FMD1.R as shown below. Either use the simulated data from STEP 0, or input and replace the inputs into the FMD function to match the objects desired.

    #Perform bFMD
    bFMDResults <- FMD(x=X, y=y)
    
## STEP 2: Analyze bFMD results
To analyze the objects returned from the FMD method and extract the identified module variables within X, execute the code below.

    #Pull out the variables identified as part of the functional module. 1=variable is in the model and 0=variable is NOT in the model.
    moduleVariables <- bFMDResults$select
    #How balanced is the detected module? 0=perfectly inbalanced (bad) and 1=perfectly balanced (ideal)
    balance <- bFMDResults$balance
    #What is the matrix that was constructed to detect module variables?
    w <- bFMDResults$w
    #What are the loadings from the sparse PCA? Note that sparse PCA was applied to matrix w (refer to manuscript for additional details)
    loadings <- bFMDResults$loadings
    
