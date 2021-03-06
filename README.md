# MATLAB Code for Knies, Lorca, and Melo (2021)
This repository contains the MATLAB code used for the simulations and resulting figures and tables generated for Knies, Lorca, and Melo (2021), "A Recursive Logit Model with Choice Aversion and Its Application to Transportation Networks." (Link: https://arxiv.org/abs/2010.02398)

The code in this repository simulates route choice probabilities and welfare calculations across the relevant transportation network examples discussed in our main text and appendices and generates the subsequent figures and tables included in our paper. The code contains simulations under various parameterizations for the choice aversion model as well as other Path Size Logit (PSL) models and extensions that we discuss. While the choice aversion model and many other PSL models have a closed-form expression, the Adaptive PSL introduced by Duncan et al. (2020) requires a fixed-point iteration algorithm. The algorithm used in our code for the APSL model is the simple Fixed-Point Iteration Method (FPIM). 

License for the MATLAB code in this repository: [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
