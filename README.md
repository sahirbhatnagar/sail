Article: Variable Selection in Regression-based Estimation of Dynamic Treatment Regimes
Author: Zeyu Bian, Erica EM Moodie, Susan M Shortreed and Sahir Bhatnagar
Email: zeyu.bian@mail.mcgill.ca
 
The enclosed R.file demonstrates the implementation of the pdWOLS method for selecting the important tailoring variables and estimating optimal dynamic treatment regimens.  We present the code of the simulation studies in Section 3.4 of the associated paper are based, where a multi-stage DTR’s estimation was implemented.

The function ‘gen_data’ is used to generate the simulated data, where all the notations of the outputs follow the ones in the main paper: Y is the outcome, X1 and X2 are covariates at first and second stage, A1 and A2 are the binary treatments at stages 1 and 2, w1 and w2 are the pdWOLS weights. The fitting procedure used the R package sail: https://sahirbhatnagar.com/sail/, specifically, the function “cv.sail” is used to conduct the four-fold cross validation. The given code in line 2 can help to install and library the corresponding package. 
