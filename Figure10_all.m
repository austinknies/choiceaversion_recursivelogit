%% Figures 10(A) and 10(B) from Knies and Melo (2020)
% Calculates route choice probabilities and welfare for Figure 4 
% Network originally featured in Fosgerau et al. (2013) (Figure 2)
% Generates Figure 10(A) and 10(B)

% Varying kappa_2
% Heterogeneous/Node-Specific Choice Aversion Model

% Route 1 : 12,23,35
% Route 2 : 12,23,34,45
% Route 3 : 12,24,45
% Route 4 : 15
clear all
clc

% Instantaneous edge costs
t12 = 1;
t23 = 1;
t35 = 2;
t34 = 1;
t45 = 1;
t24 = 2;
t15 = 4;

% Route/path costs
c1 = t12 + t23 + t35;
c2 = t12 + t23 + t34 + t45;
c3 = t12 + t24 + t45;
c4 = t15;

% Kappa is the choice aversion parameter at each node other than the sink
% (end) node (directed acyclical graph assumption)
kap1 = 1;
kap2 = 0:0.1:5; % variable of interest; what happens as we let this vary?
kap3 = 1;  
kap4 = 1;

% Theta is the parameter on the cost function
theta = 1; % unused in our model, but just an edge cost scale parameter
% used more in other PSL models


% Route utilities (minimizing costs = maximizing negative costs)
u1 = -theta*c1;
u2 = -theta*c2;
u3 = -theta*c3;
u4 = -theta*c4;

% Choice set size for each route
% Note Aji =\= the notation in the paper. 
% For simplicity of code, here Aji = product of each choice set cardinality at
% each node along route i in the network, raised to the node-specific choice aversion parameter. 

Aj1=(2.^kap1).*(2.^kap2).*(2.^kap3);
Aj2=(2.^kap1).*(2.^kap2).*(2.^kap3).*(1.^kap4);
Aj3=(2.^kap1).*(2.^kap2).*(1.^kap4);
Aj4=(2.^kap1); 


% Route probabilities
P1 = exp(u1)./((Aj1./Aj1).*exp(u1)+(Aj1./Aj2).*exp(u2)+(Aj1./Aj3).*exp(u3)+(Aj1./Aj4).*exp(u4));
P2 = exp(u2)./((Aj2./Aj1).*exp(u1)+(Aj2./Aj2).*exp(u2)+(Aj2./Aj3).*exp(u3)+(Aj2./Aj4).*exp(u4));
P3 = exp(u3)./((Aj3./Aj1).*exp(u1)+(Aj3./Aj2).*exp(u2)+(Aj3./Aj3).*exp(u3)+(Aj3./Aj4).*exp(u4));
P4 = exp(u4)./((Aj4./Aj1).*exp(u1)+(Aj4./Aj2).*exp(u2)+(Aj4./Aj3).*exp(u3)+(Aj4./Aj4).*exp(u4));
P = P1+P2+P3+P4; % unused; just one easy check to help ensure you haven't messed up the choice probabilities

% Welfare calculation
Wca2 = log((exp(u1)./(Aj1)) + (exp(u2)./(Aj2)) + (exp(u3)./(Aj3)) + (exp(u4)./(Aj4)));

%Figure 10(A)
figure
p = plot(kap2,P1,kap2,P2,kap2,P3,kap2,P4);
p(1).LineWidth = 2;
p(2).LineWidth = 2;
p(3).LineWidth = 2;
p(2).LineStyle = ':';
p(4).LineWidth = 2;
legend({'Route 1','Route 2','Route 3','Route 4'},'Location','northwest')
xlabel('\kappa_2')
ylabel('Choice Probabilities P_i')

%Figure 10(B)
figure
q = plot(kap2,Wca2);
q(1).LineWidth = 2;
%legend({'Route 1','Route 2','Route 3','Route 4'},'Location','northeast')
xlabel('\kappa_2')
ylabel('Welfare')

%% Figures 10(C) and 10(D) from Knies and Melo (2020)
% Calculates route choice probabilities and welfare for Figure 4 
% Network originally featured in Fosgerau et al. (2013) (Figure 2)
% Generates Figure 10(C) and 10(D)

% Varying kappa_3
% Heterogeneous/Node-Specific Choice Aversion Model
% Route 1 : 12,23,35
% Route 2 : 12,23,34,45
% Route 3 : 12,24,45
% Route 4 : 15
clear all
clc

% Instantaneous edge costs
t12 = 1;
t23 = 1;
t35 = 2;
t34 = 1;
t45 = 1;
t24 = 2;
t15 = 4;

% Route/path costs
c1 = t12 + t23 + t35;
c2 = t12 + t23 + t34 + t45;
c3 = t12 + t24 + t45;
c4 = t15;

% Kappa is the choice aversion parameter at each node other than the sink
% (end) node (directed acyclical graph assumption)
kap1 = 1;
kap2 = 1; 
kap3 = 0:0.1:5; % variable of interest; what happens as we let this vary? 
kap4 = 1;

% Theta is the parameter on the cost function
theta = 1; % unused in our model, but just an edge cost scale parameter
% used more in other PSL models


% Route utilities (minimizing costs = maximizing negative costs)
u1 = -theta*c1;
u2 = -theta*c2;
u3 = -theta*c3;
u4 = -theta*c4;

% Choice set size for each route
% Note Aji =\= the notation in the paper. 
% For simplicity of code, here Aji = product of each choice set cardinality at
% each node along route i in the network, raised to the node-specific choice aversion parameter. 

Aj1=(2.^kap1).*(2.^kap2).*(2.^kap3);
Aj2=(2.^kap1).*(2.^kap2).*(2.^kap3).*(1.^kap4);
Aj3=(2.^kap1).*(2.^kap2).*(1.^kap4);
Aj4=(2.^kap1); 


% Route probabilities
P1 = exp(u1)./((Aj1./Aj1).*exp(u1)+(Aj1./Aj2).*exp(u2)+(Aj1./Aj3).*exp(u3)+(Aj1./Aj4).*exp(u4));
P2 = exp(u2)./((Aj2./Aj1).*exp(u1)+(Aj2./Aj2).*exp(u2)+(Aj2./Aj3).*exp(u3)+(Aj2./Aj4).*exp(u4));
P3 = exp(u3)./((Aj3./Aj1).*exp(u1)+(Aj3./Aj2).*exp(u2)+(Aj3./Aj3).*exp(u3)+(Aj3./Aj4).*exp(u4));
P4 = exp(u4)./((Aj4./Aj1).*exp(u1)+(Aj4./Aj2).*exp(u2)+(Aj4./Aj3).*exp(u3)+(Aj4./Aj4).*exp(u4));
P = P1+P2+P3+P4; % unused; just one easy check to help ensure you haven't messed up the choice probabilities

% Welfare calculation
Wca3 = log((exp(u1)./(Aj1)) + (exp(u2)./(Aj2)) + (exp(u3)./(Aj3)) + (exp(u4)./(Aj4)));

%Figure 10(C)
figure
p = plot(kap3,P1,kap3,P2,kap3,P3,kap3,P4);
p(1).LineWidth = 2;
p(2).LineWidth = 2;
p(3).LineWidth = 2;
p(2).LineStyle = ':';
p(4).LineWidth = 2;
legend({'Route 1','Route 2','Route 3','Route 4'},'Location','northwest')
xlabel('\kappa_3')
ylabel('Choice Probabilities P_i')

%Figure 10(D)
figure
q = plot(kap3,Wca3);
q(1).LineWidth = 2;
%legend({'Route 1','Route 2','Route 3','Route 4'},'Location','northeast')
xlabel('\kappa_3')
ylabel('Welfare')
