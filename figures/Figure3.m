%% Figure 3 from  Knies and Melo (2020)
% Calculates route choice probabilities for Figure 2 for varying kappa = 0 to 2.5
% Generates Figure 3
clear all
clc 

% Instantaneous edge costs
a1 = 1.9;
a2 = 2;
a3 = 0.1;
a4 = 0.1;

% Parameter of interest: Kappa, the choice aversion parameter
kap = 0:0.1:2.5;
theta = 1; % unused in our model, but just an edge cost scale parameter
% used more in other PSL models

% Route utilities (minimizing costs = maximizing negative costs)
u1 = -theta*(a1+a3);
u2 = -theta*(a1+a4);
u3 = -theta*(a2);

% Choice set size for each route
% Note Aji =\= the notation in the paper. 
% For simplicity of code, here Aji = product of choice set cardinality at
% each node along route i in the network. 
Aj1=2*2*1;
Aj2=2*2*1;
Aj3=2*1*1;


% Route probabilities
P1 = exp(u1)./(((Aj1./Aj1).^kap).*exp(u1)+((Aj1./Aj2).^kap).*exp(u2)+((Aj1./Aj3).^kap).*exp(u3));
P2 = exp(u2)./(((Aj2./Aj1).^kap).*exp(u1)+((Aj2./Aj2).^kap).*exp(u2)+((Aj2./Aj3).^kap).*exp(u3));
P3 = exp(u3)./(((Aj3./Aj1).^kap).*exp(u1)+((Aj3./Aj2).^kap).*exp(u2)+((Aj3./Aj3).^kap).*exp(u3));
P = P1+P2+P3;

%Figure for varying kappa
figure
r = plot(kap,P1,kap,P2,kap,P3);
r(1).LineWidth = 2;
r(2).LineWidth = 2;
r(2).LineStyle = ":";
r(3).LineWidth = 2;
legend({'Route 1','Route 2','Route 3'},'Location','northwest')
xlabel('\kappa')
ylabel('Path Choice Probability')