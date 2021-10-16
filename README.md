# MATLAB Code for Knies, Lorca, and Melo (2021)
This repository contains the MATLAB code used for the simulations, MLE routines, and resulting figures and tables generated for Knies, Lorca, and Melo (2021), "A Recursive Logit Model with Choice Aversion and Its Application to Transportation Networks," forthcoming at Transportation Research Part B: Methodological (Link to working paper: https://arxiv.org/abs/2010.02398)

The code in this repository simulates route choice probabilities and welfare calculations across the relevant transportation network examples discussed in our main text and appendices and generates the subsequent figures and tables included in our paper. The code contains simulations under various parameterizations for the choice aversion model as well as other Path Size Logit (PSL) and Recursive Logit (RL) models and extensions that we discuss. While the choice aversion model and many other PSL models have a closed-form expression, the Adaptive PSL introduced by Duncan et al. (2020) requires a fixed-point iteration algorithm. The algorithm used in our code for the APSL model is the simple Fixed-Point Iteration Method (FPIM). The Nested Recursive Logit introduced by Mai et al. (2015) also requires a fixed-point iteration algorithm for the continuation value function, which we have employed in our simulations.

The MLE routines are modifications of existing RL code developed by Tien Mai (https://github.com/maitien86/RL-Tutorial). Alongside incorporating choice aversion into the existing RL code and incorporating the Borl&auml;nge data used in Fosgerau et al. (2013) (https://github.com/maitien86/RecursiveLogit.Classical.V2/tree/master/Input), we also make the following changes:
   -   Link 7 only has two links attached in the network depicted in 
       tutorial slides: 10 and 17. We break from the diagram and allow 19
       to connect from 7 to 29, since it would otherwise be missing. Our recreation of this network (Figure 9) reflects this.
   -   In the incidence matrix, link 20 appears to connect to itself. 
       We have removed this, so that link 20 only "connects" to 29.
   -   In loadData.m, the following code erased the connection between
       link 20 and 29:
       
                   icd(:,nbNetworkStates:nbTotalStates) = 0;
                   
       I have changed it to
       
                   icd(:,nbTotalStates) = 0;
                   
                   
To conduct the cross-validation performed in Appendix D and reproduce the random samples, one should follow these steps:
   -   First, run generateHoldoutSample.m to obtain random samples for training and testing. 
   -   Next, run prediction_Borlange.m to obtain estimates on each training sample and then use those estimates to obtain the predicted LL on the corresponding testing sample. Since prediction_Borlange.m is a function with output that gives you estimated and predicted LL for each random sample coupling, you can input the following code:
   
                  [EstimatedLL, PredictedLL] = prediction_Borlange(isLS,true);
                  
       If you want to include the LS term, set isLS = true. Otherwise, set isLS = false.

License for the MATLAB code in this repository: [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
