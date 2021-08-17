------------------------------------------------------------
Recursive logit route choice model estimation 
------------------------------------------------------------
%   Modified by Austin Knies - Indiana University Bloomington (August 2021)
%   Incorporates choice aversion term as in Knies, Lorca, and Melo (2021)

%   Changes to original code:
%   -   Link 7 only has two links attached in the network depicted in 
%       tutorial slides: 10 and 17. We break from the diagram and allow 19
%       to connect from 7 to 29, since it would otherwise be missing.
%   -   In the incidence matrix, link 20 appears to connect to itself. 
%       We have removed this, so that link 20 only "connects" to 29.
%   -   In loadData.m, the following code erased the connection between
%       link 20 and 29:
%                   icd(:,nbNetworkStates:nbTotalStates) = 0;
%       I have changed it to
%                   icd(:,nbTotalStates) = 0;
