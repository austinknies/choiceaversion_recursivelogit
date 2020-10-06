%% Part 1 - Figure 8 from Knies and Melo (2020)
% MUST RUN BEFORE FIGURE8PART2.M (otherwise will not have proper .mat file)
% Calculates welfare across various PSL models and choice aversion model
% for Figure 7(A)


%% Choice Aversion Model
clear all
clc 

% Based on Braess_A.pdf
% Route 1  = a1 + a2
% Route 2 = a3 + a4
% a1 = a4 = x, a2 = a3 = 1
% Letting x vary from 0.5 to 3.
x = 0:0.1:3.0;
x(1) = 10^(-16);
% Kappa is the penalty term for the log of outgoing edge set cardinality
kap = 1;
% Theta is the parameter on the cost function
theta = 1;

% Route costs
c1 = x+1;
c2 = 1+x;

% Route utilities
u1 = -theta*c1;
u2 = -theta*c2;

% Choice set size for each route

Aj1=2*1*1;
Aj2=2*1*1;


% Route probabilities
P1 = exp(u1)./(((Aj1/Aj1)^kap).*exp(u1)+((Aj1/Aj2)^kap).*exp(u2));
P2 = exp(u2)./(((Aj2/Aj1)^kap).*exp(u1)+((Aj2/Aj2)^kap).*exp(u2));
P = P1+P2;

% Welfare calculation
Vca = log((exp(u1)./((Aj1)^kap)) + (exp(u2)./((Aj2)^kap)));
save Figure8part1.mat x Vca
% % plots route choice probabilities
% figure
% p = plot(x,P1,x,P2);
% p(1).LineWidth = 2;
% p(2).LineWidth = 2;
% p(2).LineStyle = ':';
% legend({'Route 1','Route 2'},'Location','southwest')
% xlabel('x')
% ylabel('Choice Probabilities P_i')

%% Choice Aversion (kappa = 2)

clear all
clc 

% Based on Braess_A.pdf
% Route 1  = a1 + a2
% Route 2 = a3 + a4
% a1 = a4 = x, a2 = a3 = 1
% Letting x vary from 0.5 to 3.
x = 0:0.1:3.0;
x(1) = 10^(-16);
% Kappa is the penalty term for the log of outgoing edge set cardinality
kap = 2;
% Theta is the parameter on the cost function
theta = 1;

% Route costs
c1 = x+1;
c2 = 1+x;

% Route utilities
u1 = -theta*c1;
u2 = -theta*c2;

% Choice set size for each route

Aj1=2*1*1;
Aj2=2*1*1;


% Route probabilities
P1 = exp(u1)./(((Aj1/Aj1)^kap).*exp(u1)+((Aj1/Aj2)^kap).*exp(u2));
P2 = exp(u2)./(((Aj2/Aj1)^kap).*exp(u1)+((Aj2/Aj2)^kap).*exp(u2));
P = P1+P2;

% Welfare calculation
Vca5 = log((exp(u1)./((Aj1)^kap)) + (exp(u2)./((Aj2)^kap)));
save('Figure8part1.mat','Vca5','-append')
% % plots route choice probabilities
% figure
% p = plot(x,P1,x,P2);
% p(1).LineWidth = 2;
% p(2).LineWidth = 2;
% p(2).LineStyle = ':';
% legend({'Route 1','Route 2'},'Location','southwest')
% xlabel('x')
% ylabel('Choice Probabilities P_i')

%% MNL model

clear all
clc 

% Based on Braess_A.pdf
% Route 1  = a1 + a2
% Route 2 = a3 + a4
% a1 = a4 = x, a2 = a3 = 1
% Letting x vary from 0.5 to 3.
x = 0:0.1:3.0;
x(1) = 10^(-16);
% Kappa is the penalty term for the log of outgoing edge set cardinality
kap = 0; % This makes the choice aversion model equivalent to MNL model
% Theta is the parameter on the cost function
theta = 1;
% Route costs
c1 = x+1;
c2 = 1+x;

% Route utilities
u1 = -theta*c1;
u2 = -theta*c2;

% Choice set size for each route

Aj1=2*1*1;
Aj2=2*1*1;


% Route probabilities
P1 = exp(u1)./(((Aj1/Aj1)^kap).*exp(u1)+((Aj1/Aj2)^kap).*exp(u2));
P2 = exp(u2)./(((Aj2/Aj1)^kap).*exp(u1)+((Aj2/Aj2)^kap).*exp(u2));
P = P1+P2;

% Welfare calculation
Vmnl = log((exp(u1)./((Aj1)^kap)) + (exp(u2)./((Aj2)^kap)));
save('Figure8part1.mat','Vmnl','-append')

% % plots route choice probabilities
% figure
% p = plot(x,P1,x,P2);
% p(1).LineWidth = 2;
% p(2).LineWidth = 2;
% p(2).LineStyle = ':';
% legend({'Route 1','Route 2'},'Location','southwest')
% xlabel('x')
% ylabel('Choice Probabilities P_i')

%% Path Size Logit Model

clear all
clc 

% Based on Braess_A.pdf
% Route 1  = a1 + a2
% Route 2 = a3 + a4
% a1 = a4 = x, a2 = a3 = 1
% Letting x vary from 0.5 to 3.
x = 0:0.1:3.0;
x(1) = 10^(-16);
% Kappa is the penalty term for the log of outgoing edge set cardinality
kap = 1;
% Theta is the parameter on the cost function
theta = 1;

% Route costs
c1 = x+1;
c2 = 1+x;

% Route utilities
u1 = -theta*c1;
u2 = -theta*c2;

% Gamma terms
% Travel time for each edge a
ta1 = x;
ta2 = 1;
ta3 = 1;
ta4 = x;
% Incidence for each edge in a route
delta1 = 1;
delta2 = 1;
delta3 = 1;
delta4 = 1;


Gamma1=(ta1./c1).*(1/delta1)+(ta2./c1).*(1/delta2);
Gamma2=(ta3./c2).*(1/delta3)+(ta4./c2).*(1/delta4);


% Route probabilities
P1 = exp(u1)./(((Gamma1./Gamma1).^kap).*exp(u1)+((Gamma2./Gamma1).^kap).*exp(u2));
P2 = exp(u2)./(((Gamma1./Gamma2).^kap).*exp(u1)+((Gamma2./Gamma2).^kap).*exp(u2));
P = P1+P2;

% Welfare calculations
Vpsl = log((exp(u1).*((Gamma1).^kap)) + (exp(u2).*((Gamma2).^kap)));
save('Figure8part1.mat','Vpsl','-append')

% % plots route choice probabilities
% figure
% p = plot(x,P1,x,P2);
% p(1).LineWidth = 2;
% p(2).LineWidth = 2;
% p(2).LineStyle = ':';
% legend({'Route 1','Route 2'},'Location','southwest')
% xlabel('x')
% ylabel('Choice Probabilities P_i')

%% Alternative PSL (PSL') Model

clear all
clc 

% Based on Braess_A.pdf
% Route 1  = a1 + a2
% Route 2 = a3 + a4
% a1 = a4 = x, a2 = a3 = 1
% Letting x vary from 0.5 to 3.
x = 0:0.1:3.0;
x(1) = 10^(-16);
% Kappa is the penalty term for the log of outgoing edge set cardinality
kap = 1;
% Theta is the parameter on the cost function
theta = 1;

% Route costs
c1 = x+1;
c2 = 1+x;

% Route utilities
u1 = -theta*c1;
u2 = -theta*c2;

% Gamma terms
% Travel time for each edge a
ta1 = x;
ta2 = 1;
ta3 = 1;
ta4 = x;
% Incidence for each edge in a route
% Includes a weighting term min(cl:l in R)/ck
for i=1:length(c1)
    minc(i) = min([c1(i),c2(i)]);
end % Finds the minimum cost between 3 routes for each value of x
delta1 = (minc./c1);
delta2 = (minc./c1);
delta3 = (minc./c2);
delta4 = (minc./c2);


Gamma1=(ta1./c1).*(1./delta1)+(ta2./c1).*(1./delta2);
Gamma2=(ta3./c2).*(1./delta3)+(ta4./c2).*(1./delta4);


% Route probabilities
P1 = exp(u1)./(((Gamma1./Gamma1).^kap).*exp(u1)+((Gamma2./Gamma1).^kap).*exp(u2));
P2 = exp(u2)./(((Gamma1./Gamma2).^kap).*exp(u1)+((Gamma2./Gamma2).^kap).*exp(u2));
P = P1+P2;

% Welfare calculations
Valtpsl = log((exp(u1).*((Gamma1).^kap)) + (exp(u2).*((Gamma2).^kap)));
save('Figure8part1.mat','Valtpsl','-append')

% % plots route choice probabilities
% figure
% p = plot(x,P1,x,P2);
% p(1).LineWidth = 2;
% p(2).LineWidth = 2;
% p(2).LineStyle = ':';
% legend({'Route 1','Route 2'},'Location','southwest')
% xlabel('x')
% ylabel('Choice Probabilities P_i')

%% GPSL Model

clear all
clc 

% Based on Braess_A.pdf
% Route 1  = a1 + a2
% Route 2 = a3 + a4
% a1 = a4 = x, a2 = a3 = 1
% Letting x vary from 0.5 to 3.
x = 0:0.1:3.0;
x(1) = 10^(-16);
% Kappa is the penalty term for the log of outgoing edge set cardinality
kap = 1;
% Theta is the parameter on the cost function
theta = 1;

% Route costs
c1 = x+1;
c2 = 1+x;

% Route utilities
u1 = -theta*c1;
u2 = -theta*c2;

% Gamma terms
% Travel time for each edge a
ta1 = x;
ta2 = 1;
ta3 = 1;
ta4 = x;
% Incidence for each edge in a route
lambda = 3; % this is similar to Duncan et al. (2020) parameterization
delta1 = ((1./c1).^lambda);
delta2 = ((1./c1).^lambda);
delta3 = ((1./c2).^lambda);
delta4 = ((1./c2).^lambda);


Gamma1=(ta1./c1).*(1./((c1.^lambda).*delta1))+(ta2./c1).*(1./((c1.^lambda).*delta2));
Gamma2=(ta3./c2).*(1./((c2.^lambda).*delta3))+(ta4./c2).*(1./((c2.^lambda).*delta4));


% Route probabilities
P1 = exp(u1)./(((Gamma1./Gamma1).^kap).*exp(u1)+((Gamma2./Gamma1).^kap).*exp(u2));
P2 = exp(u2)./(((Gamma1./Gamma2).^kap).*exp(u1)+((Gamma2./Gamma2).^kap).*exp(u2));
P = P1+P2;

% Welfare calculations
Vgpsl = log((exp(u1).*((Gamma1).^kap)) + (exp(u2).*((Gamma2).^kap)));
save('Figure8part1.mat','Vgpsl','-append')

% % plots route choice probabilities
% figure
% p = plot(x,P1,x,P2);
% p(1).LineWidth = 2;
% p(2).LineWidth = 2;
% p(2).LineStyle = ':';
% legend({'Route 1','Route 2'},'Location','southwest')
% xlabel('x')
% ylabel('Choice Probabilities P_i')

%% GPSL' Model

clear all
clc 

% Based on Braess_A.pdf
% Route 1  = a1 + a2
% Route 2 = a3 + a4
% a1 = a4 = x, a2 = a3 = 1
% Letting x vary from 0.5 to 3.
x = 0:0.1:3.0;
x(1) = 10^(-16);
% Kappa is the penalty term for the log of outgoing edge set cardinality
kap = 1;
% Theta is the parameter on the cost function
theta = 1;

% Route costs
c1 = x+1;
c2 = 1+x;

% Route utilities
u1 = -theta*c1;
u2 = -theta*c2;

% Gamma terms
% Travel time for each edge a
ta1 = x;
ta2 = 1;
ta3 = 1;
ta4 = x;
% Incidence for each edge in a route
lambda = theta; 
delta1 = exp(-lambda.*c1);
delta2 = exp(-lambda.*c1);
delta3 = exp(-lambda.*c2);
delta4 = exp(-lambda.*c2);


Gamma1=(ta1./c1).*(1./((exp(lambda.*c1)).*delta1))+(ta2./c1).*(1./((exp(lambda.*c1)).*delta2));
Gamma2=(ta3./c2).*(1./((exp(lambda.*c2)).*delta3))+(ta4./c2).*(1./((exp(lambda.*c2)).*delta4));


% Route probabilities
P1 = exp(u1)./(((Gamma1./Gamma1).^kap).*exp(u1)+((Gamma2./Gamma1).^kap).*exp(u2));
P2 = exp(u2)./(((Gamma1./Gamma2).^kap).*exp(u1)+((Gamma2./Gamma2).^kap).*exp(u2));
P = P1+P2;

% Welfare calculations
Valtgpsl = log((exp(u1).*((Gamma1).^kap)) + (exp(u2).*((Gamma2).^kap)));
save('Figure8part1.mat','Valtgpsl','-append')

% % plots route choice probabilities
% figure
% p = plot(x,P1,x,P2);
% p(1).LineWidth = 2;
% p(2).LineWidth = 2;
% p(2).LineStyle = ':';
% legend({'Route 1','Route 2'},'Location','southwest')
% xlabel('x')
% ylabel('Choice Probabilities P_i')

%% Adaptive PSL Model

clear all
clc 

% Based on Braess_A.pdf
% Route 1  = a1 + a2
% Route 2 = a3 + a4
% a1 = a4 = x, a2 = a3 = 1
% Letting x vary from 0 to 3.
% Can't use x here because I use it as a variable for the functions defined
% below
ex = 10^(-16); %as an approximation for zero
% Beta is the penalty term for gamma
beta = 1;
% Theta is the parameter on the cost function
theta = 1;

% Initialization
tau = 10^(-16);
N = 2;

i = 1;
error = 8; %% Duncan et al. (2020) states "predetermined convergence parameter"; cut error in half for easier convergence
I = 1000;
j = 1;

for j=1:length(0:0.1:3.0)
% Cost of routes
c(1) = ex+1;
c(2) = 1+ex;

% Route utilities
u(1) = -theta*c(1);
u(2) = -theta*c(2);

% Gamma terms
% Travel time for each edge a
ta(1) = ex;
ta(2) = 1;
ta(3) = 1;
ta(4) = ex;

% Initial guess
P0(1:2) = 1/N;


syms gamma1(x,y) gamma2(x,y) gamma3(x,y) g1(x,y) g2(x,y) g3(x,y) G1(x,y) G2(x,y) G3(x,y)

gamma1(x,y) = (ta(1)/c(1))*(x/(x))+(ta(2)/c(1))*(x/(x));
gamma2(x,y) = (ta(3)/c(2))*(y/(y))+(ta(4)/c(2))*(y/(y));

g1(x,y) = (((gamma1(x,y))^beta)*exp(-theta*c(1)))/(((gamma1(x,y))^beta)*exp(-theta*c(1))+((gamma2(x,y))^beta)*exp(-theta*c(2)));
g2(x,y) = (((gamma2(x,y))^beta)*exp(-theta*c(2)))/(((gamma1(x,y))^beta)*exp(-theta*c(1))+((gamma2(x,y))^beta)*exp(-theta*c(2)));

G1(x,y) = tau + (1-N*tau)*g1(x,y);
G2(x,y) = tau + (1-N*tau)*g2(x,y);





    
while i <= I
    
    P(1) = G1(P0(1),P0(2));
    P(2) = G2(P0(1),P0(2));
    P=double(P);
    pcheck = log(sum(abs(P-P0)));
    ptest = log(10^(-error));
    if log(sum(abs(P-P0)))<log(10^(-error))
        %fprintf('Solution has been found!')
        break
    end
i = i+1;
P0=P;
end
if i == I
    fprintf('Solution did not converge')
end
ChoiceProbabilities(j,:) = P(1:N);
Vadaptive(j) = log((exp(u(1)).*((double(gamma1(P(1),P(2)))).^beta)) + (exp(u(2)).*((double(gamma2(P(1),P(2)))).^beta)));
% j = j+1;
ex = ex + 0.1;

end

save('Figure8part1.mat','Vadaptive','-append')

% % plots route choice probabilities
% figure
% r = plot(0:length(ChoiceProbabilities)-1,ChoiceProbabilities);
% r(1).LineWidth = 2;
% r(2).LineWidth = 2;
% r(2).LineStyle = ':';
% legend({'Route 1','Route 2'},'Location','southwest')
% xlabel('x')
% xticklabels({'0','0.5','1','1.5','2','2.5','3'})
% ylabel('Choice Probability P_i')

%% Welfare (plots just the direct welfare calculation)

% clear
% load Figure8part1.mat
% 
% figure
% q = plot(x,Vca,x,Vca5,x,Vmnl,x,Vpsl,x,Valtpsl,x,Vgpsl,x,Valtgpsl,x,Vadaptive);
% q(1).LineWidth = 2;
% q(2).LineWidth = 2;
% q(3).LineWidth = 2;
% q(4).LineWidth = 2;
% q(5).LineWidth = 2;
% q(6).LineWidth = 2;
% q(7).LineWidth = 2;
% q(8).LineWidth = 2;
% q(8).LineStyle = ':';
% legend({'Choice Aversion, \kappa = 1','Choice Aversion, \kappa = 5','MNL','PSL','Alt. PSL','GPSL, \lambda = 3','Alt. GPSL, \lambda = \theta = 1','Adaptive PSL'},'Location','southwest')
% xlabel('x')
% ylabel('Welfare')


