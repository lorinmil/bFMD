# bFMD
Balanced Functional Module Detection

## Preliminaries
This program takes in a high-dimensional dataset, X, and a continuous outcome, y, and applies the balanced Functional Module Detection (bFMD) method. Essentially, the purpose of bFMD is to identify a subset of variables within X that are jointly correlated with each other AND related to the outcome y. These subset of variables indicate nodes within a semidirected balanced graph. Balance suggests that the directionality between module variables and the outcome do not contradict each other in order to improve biological interpretability (e.g. if variables x1 and x2 are positively correlated with each other, then sign(cor(x1, y)) should equal sign(cor(x2, y)). Details of this method may be found in (INSERT PREPRINT HERE)

## STEP 0: Simulations (optional)
