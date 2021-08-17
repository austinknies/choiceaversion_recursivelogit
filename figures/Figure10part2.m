%% Part 2 - Figure 10 from Knies, Lorca, and Melo (2021)
% MUST RUN AFTER FIGURE10PART1.M (otherwise will not have proper .mat file)
% Calculates welfare across various PSL models and choice aversion model
% for Figure 8(B)
% Last Updated: August 2021

%% Choice Aversion Model
clear all
clc 

% Based on Braess_B.pdf
% Route 1  = a1 + a2
% Route 2 = a3 + a4
% Route 3 = a1 + a5 + a4
% a1 = a4 = x, a2 = a3 = 1, a5 = 0
% Letting x vary from 0 to 3.
% If x = 0.5, c3 = 1 < 1.5 = c1 = c2
% If x = 1, c3 = 2 = c1 = c2
% If x = 1.5, c3 = 3 > 2.5 = c1 = c2
% If x = 3, c3 = 6 > 4 = c1 = c2
x = 0:0.1:3.0;
x(1) = 10^(-16);
% Kappa is the penalty term for the log of outgoing edge set cardinality
kap = 1;
% Theta is the parameter on the cost function
theta = 1;

% Route costs
c1 = x+1;
c2 = 1+x;
c3 = x+0+x;

% Route utilities
u1 = -theta*c1;
u2 = -theta*c2;
u3 = -theta*c3;

% Choice set size for each route

Aj1=2*2*1;
Aj2=2*1*1;
Aj3=2*2*1;


% Route probabilities
P1 = exp(u1)./(((Aj1/Aj1)^kap).*exp(u1)+((Aj1/Aj2)^kap).*exp(u2)+((Aj1/Aj3)^kap).*exp(u3));
P2 = exp(u2)./(((Aj2/Aj1)^kap).*exp(u1)+((Aj2/Aj2)^kap).*exp(u2)+((Aj2/Aj3)^kap).*exp(u3));
P3 = exp(u3)./(((Aj3/Aj1)^kap).*exp(u1)+((Aj3/Aj2)^kap).*exp(u2)+((Aj3/Aj3)^kap).*exp(u3));
P = P1+P2+P3;

% Welfare calculation
Wca = log((exp(u1)./((Aj1)^kap)) + (exp(u2)./((Aj2)^kap)) + (exp(u3)./((Aj3)^kap)));
save Figure10part2.mat x Wca
% % Plots route choice probabilities
% figure
% p = plot(x,P1,x,P2,x,P3);
% p(1).LineWidth = 2;
% p(2).LineWidth = 2;
% p(3).LineWidth = 2;
% legend({'Route 1','Route 2','Route 3'},'Location','southwest')
% xlabel('x')
% ylabel('Choice Probabilities P_i')

%% Choice Aversion (kappa = 2)
clear all
clc 

% Based on Braess_B.pdf
% Route 1  = a1 + a2
% Route 2 = a3 + a4
% Route 3 = a1 + a5 + a4
% a1 = a4 = x, a2 = a3 = 1, a5 = 0
% Letting x vary from 0 to 3.
% If x = 0.5, c3 = 1 < 1.5 = c1 = c2
% If x = 1, c3 = 2 = c1 = c2
% If x = 1.5, c3 = 3 > 2.5 = c1 = c2
% If x = 3, c3 = 6 > 4 = c1 = c2
x = 0:0.1:3.0;
x(1) = 10^(-16);
% Kappa is the penalty term for the log of outgoing edge set cardinality
kap = 2;
% Theta is the parameter on the cost function
theta = 1;

% Route costs
c1 = x+1;
c2 = 1+x;
c3 = x+0+x;

% Route utilities
u1 = -theta*c1;
u2 = -theta*c2;
u3 = -theta*c3;

% Choice set size for each route

Aj1=2*2*1;
Aj2=2*1*1;
Aj3=2*2*1;


% Route probabilities
P1 = exp(u1)./(((Aj1/Aj1)^kap).*exp(u1)+((Aj1/Aj2)^kap).*exp(u2)+((Aj1/Aj3)^kap).*exp(u3));
P2 = exp(u2)./(((Aj2/Aj1)^kap).*exp(u1)+((Aj2/Aj2)^kap).*exp(u2)+((Aj2/Aj3)^kap).*exp(u3));
P3 = exp(u3)./(((Aj3/Aj1)^kap).*exp(u1)+((Aj3/Aj2)^kap).*exp(u2)+((Aj3/Aj3)^kap).*exp(u3));
P = P1+P2+P3;

% Welfare calculation
Wca5 = log((exp(u1)./((Aj1)^kap)) + (exp(u2)./((Aj2)^kap)) + (exp(u3)./((Aj3)^kap)));
save('Figure10part2.mat','Wca5','-append')
% % Plots route choice probabilities
% figure
% p = plot(x,P1,x,P2,x,P3);
% p(1).LineWidth = 2;
% p(2).LineWidth = 2;
% p(3).LineWidth = 2;
% legend({'Route 1','Route 2','Route 3'},'Location','southwest')
% xlabel('x')
% ylabel('Choice Probabilities P_i')

%% MNL model
clear all
clc 

% Based on Braess_B.pdf
% Route 1  = a1 + a2
% Route 2 = a3 + a4
% Route 3 = a1 + a5 + a4
% a1 = a4 = x, a2 = a3 = 1, a5 = 0
% Letting x vary from 0 to 3.
% If x = 0.5, c3 = 1 < 1.5 = c1 = c2
% If x = 1, c3 = 2 = c1 = c2
% If x = 1.5, c3 = 3 > 2.5 = c1 = c2
% If x = 3, c3 = 6 > 4 = c1 = c2
x = 0:0.1:3.0;
x(1) = 10^(-16);
% Kappa is the penalty term for the log of outgoing edge set cardinality
kap = 0; % This makes the choice aversion model equivalent to MNL model
% Theta is the parameter on the cost function
theta = 1;

% Route costs
c1 = x+1;
c2 = 1+x;
c3 = x+0+x;

% Route utilities
u1 = -theta*c1;
u2 = -theta*c2;
u3 = -theta*c3;

% Choice set size for each route

Aj1=2*2*1;
Aj2=2*1*1;
Aj3=2*2*1;


% Route probabilities
P1 = exp(u1)./(((Aj1/Aj1)^kap).*exp(u1)+((Aj1/Aj2)^kap).*exp(u2)+((Aj1/Aj3)^kap).*exp(u3));
P2 = exp(u2)./(((Aj2/Aj1)^kap).*exp(u1)+((Aj2/Aj2)^kap).*exp(u2)+((Aj2/Aj3)^kap).*exp(u3));
P3 = exp(u3)./(((Aj3/Aj1)^kap).*exp(u1)+((Aj3/Aj2)^kap).*exp(u2)+((Aj3/Aj3)^kap).*exp(u3));
P = P1+P2+P3;

% Welfare calculation
Wmnl = log((exp(u1)./((Aj1)^kap)) + (exp(u2)./((Aj2)^kap)) + (exp(u3)./((Aj3)^kap)));
save('Figure10part2.mat','Wmnl','-append')

% %Plots route choice probabilities
% figure
% p = plot(x,P1,x,P2,x,P3);
% p(1).LineWidth = 2;
% p(2).LineWidth = 2;
% p(2).LineStyle = ':';
% p(3).LineWidth = 2;
% legend({'Route 1','Route 2','Route 3'},'Location','southwest')
% xlabel('x')
% ylabel('Choice Probabilities P_i')

%% Path Size Logit Model

clear all
clc 

% Based on Braess_B.pdf
% Route 1  = a1 + a2
% Route 2 = a3 + a4
% Route 3 = a1 + a5 + a4
% a1 = a4 = x, a2 = a3 = 1, a5 = 0
% Letting x vary from 0 to 3.
% If x = 0.5, c3 = 1 < 1.5 = c1 = c2
% If x = 1, c3 = 2 = c1 = c2
% If x = 1.5, c3 = 3 > 2.5 = c1 = c2
% If x = 3, c3 = 6 > 4 = c1 = c2
x = 0:0.1:3.0;
x(1) = 10^(-16);
% Kappa is the penalty term for the log of outgoing edge set cardinality
kap = 1;
% Theta is the parameter on the cost function
theta = 1;

% Route costs
c1 = x+1;
c2 = 1+x;
c3 = x+0+x;

% Route utilities
u1 = -theta*c1;
u2 = -theta*c2;
u3 = -theta*c3;

% Gamma terms
% Travel time for each edge a
ta1 = x;
ta2 = 1;
ta3 = 1;
ta4 = x;
ta5 = 0;
% Incidence for each edge in a route
delta1 = 2;
delta2 = 1;
delta3 = 1;
delta4 = 2;
delta5 = 1;


Gamma1=(ta1./c1).*(1/delta1)+(ta2./c1).*(1/delta2);
Gamma2=(ta3./c2).*(1/delta3)+(ta4./c2).*(1/delta4);
Gamma3=(ta1./c3).*(1/delta1)+(ta5./c3).*(1/delta5)+(ta4./c3).*(1/delta4);


% Route probabilities
P1 = exp(u1)./(((Gamma1./Gamma1).^kap).*exp(u1)+((Gamma2./Gamma1).^kap).*exp(u2)+((Gamma3./Gamma1).^kap).*exp(u3));
P2 = exp(u2)./(((Gamma1./Gamma2).^kap).*exp(u1)+((Gamma2./Gamma2).^kap).*exp(u2)+((Gamma3./Gamma2).^kap).*exp(u3));
P3 = exp(u3)./(((Gamma1./Gamma3).^kap).*exp(u1)+((Gamma2./Gamma3).^kap).*exp(u2)+((Gamma3./Gamma3).^kap).*exp(u3));
P = P1+P2+P3;

% Welfare calculations
Wpsl = log((exp(u1).*((Gamma1).^kap)) + (exp(u2).*((Gamma2).^kap)) + (exp(u3).*((Gamma3).^kap)));
save('Figure10part2.mat','Wpsl','-append')

% %Plots route choice probabilities
% figure
% p = plot(x,P1,x,P2,x,P3);
% p(1).LineWidth = 2;
% p(2).LineWidth = 2;
% p(2).LineStyle = ':';
% p(3).LineWidth = 2;
% legend({'Route 1','Route 2','Route 3'},'Location','southwest')
% xlabel('x')
% ylabel('Choice Probabilities P_i')

%% Alternative PSL (PSL') Model

clear all
clc 

% Based on Braess_B.pdf
% Route 1  = a1 + a2
% Route 2 = a3 + a4
% Route 3 = a1 + a5 + a4
% a1 = a4 = x, a2 = a3 = 1, a5 = 0
% Letting x vary from 0 to 3.
% If x = 0.5, c3 = 1 < 1.5 = c1 = c2
% If x = 1, c3 = 2 = c1 = c2
% If x = 1.5, c3 = 3 > 2.5 = c1 = c2
% If x = 3, c3 = 6 > 4 = c1 = c2
x = 0:0.1:3.0;
x(1) = 10^(-16);
% Kappa is the penalty term for the log of outgoing edge set cardinality
kap = 1;
% Theta is the parameter on the cost function
theta = 1;

% Route costs
c1 = x+1;
c2 = 1+x;
c3 = x+0+x;

% Route utilities
u1 = -theta*c1;
u2 = -theta*c2;
u3 = -theta*c3;

% Gamma terms
% Travel time for each edge a
ta1 = x;
ta2 = 1;
ta3 = 1;
ta4 = x;
ta5 = 0;
% Incidence for each edge in a route
% Includes a weighting term min(cl:l in R)/ck
for i=1:length(c1)
    minc(i) = min([c1(i),c2(i),c3(i)]);
end % Finds the minimum cost between 3 routes for each value of x
delta1 = (minc./c1)+(minc./c3);
delta2 = (minc./c1);
delta3 = (minc./c2);
delta4 = (minc./c2)+(minc./c3);
delta5 = (minc./c3);


Gamma1=(ta1./c1).*(1./delta1)+(ta2./c1).*(1./delta2);
Gamma2=(ta3./c2).*(1./delta3)+(ta4./c2).*(1./delta4);
Gamma3=(ta1./c3).*(1./delta1)+(ta5./c3).*(1./delta5)+(ta4./c3).*(1./delta4);


% Route probabilities
P1 = exp(u1)./(((Gamma1./Gamma1).^kap).*exp(u1)+((Gamma2./Gamma1).^kap).*exp(u2)+((Gamma3./Gamma1).^kap).*exp(u3));
P2 = exp(u2)./(((Gamma1./Gamma2).^kap).*exp(u1)+((Gamma2./Gamma2).^kap).*exp(u2)+((Gamma3./Gamma2).^kap).*exp(u3));
P3 = exp(u3)./(((Gamma1./Gamma3).^kap).*exp(u1)+((Gamma2./Gamma3).^kap).*exp(u2)+((Gamma3./Gamma3).^kap).*exp(u3));
P = P1+P2+P3;

% Welfare calculations
Waltpsl = log((exp(u1).*((Gamma1).^kap)) + (exp(u2).*((Gamma2).^kap)) + (exp(u3).*((Gamma3).^kap)));
save('Figure10part2.mat','Waltpsl','-append')

% % Plots route choice probabilities
% figure
% p = plot(x,P1,x,P2,x,P3);
% p(1).LineWidth = 2;
% p(2).LineWidth = 2;
% p(2).LineStyle = ':';
% p(3).LineWidth = 2;
% legend({'Route 1','Route 2','Route 3'},'Location','southwest')
% xlabel('x')
% ylabel('Choice Probabilities P_i')

%% GPSL Model

clear all
clc 

% Based on Braess_B.pdf
% Route 1  = a1 + a2
% Route 2 = a3 + a4
% Route 3 = a1 + a5 + a4
% a1 = a4 = x, a2 = a3 = 1, a5 = 0
% Letting x vary from 0 to 3.
% If x = 0.5, c3 = 1 < 1.5 = c1 = c2
% If x = 1, c3 = 2 = c1 = c2
% If x = 1.5, c3 = 3 > 2.5 = c1 = c2
% If x = 3, c3 = 6 > 4 = c1 = c2
x = 0:0.1:3.0;
x(1) = 10^(-16);
% Kappa is the penalty term for the log of outgoing edge set cardinality
kap = 1;
% Theta is the parameter on the cost function
theta = 1;

% Route costs
c1 = x+1;
c2 = 1+x;
c3 = x+0+x;

% Route utilities
u1 = -theta*c1;
u2 = -theta*c2;
u3 = -theta*c3;

% Gamma terms
% Travel time for each edge a
ta1 = x;
ta2 = 1;
ta3 = 1;
ta4 = x;
ta5 = 0;
% Incidence for each edge in a route
lambda = 3; % this is similar to Duncan et al. (2020) parameterization
delta1 = ((1./c1).^lambda)+((1./c3).^lambda);
delta2 = ((1./c1).^lambda);
delta3 = ((1./c2).^lambda);
delta4 = ((1./c2).^lambda)+((1./c3).^lambda);
delta5 = ((1./c3).^lambda);


Gamma1=(ta1./c1).*(1./((c1.^lambda).*delta1))+(ta2./c1).*(1./((c1.^lambda).*delta2));
Gamma2=(ta3./c2).*(1./((c2.^lambda).*delta3))+(ta4./c2).*(1./((c2.^lambda).*delta4));
Gamma3=(ta1./c3).*(1./((c3.^lambda).*delta1))+(ta5./c3).*(1./((c3.^lambda).*delta5))+(ta4./c3).*(1./((c3.^lambda).*delta4));


% Route probabilities
P1 = exp(u1)./(((Gamma1./Gamma1).^kap).*exp(u1)+((Gamma2./Gamma1).^kap).*exp(u2)+((Gamma3./Gamma1).^kap).*exp(u3));
P2 = exp(u2)./(((Gamma1./Gamma2).^kap).*exp(u1)+((Gamma2./Gamma2).^kap).*exp(u2)+((Gamma3./Gamma2).^kap).*exp(u3));
P3 = exp(u3)./(((Gamma1./Gamma3).^kap).*exp(u1)+((Gamma2./Gamma3).^kap).*exp(u2)+((Gamma3./Gamma3).^kap).*exp(u3));
P = P1+P2+P3;

% Welfare calculations
Wgpsl = log((exp(u1).*((Gamma1).^kap)) + (exp(u2).*((Gamma2).^kap)) + (exp(u3).*((Gamma3).^kap)));
save('Figure10part2.mat','Wgpsl','-append')

% % Plots route choice probabilities
% figure
% p = plot(x,P1,x,P2,x,P3);
% p(1).LineWidth = 2;
% p(2).LineWidth = 2;
% p(2).LineStyle = ':';
% p(3).LineWidth = 2;
% legend({'Route 1','Route 2','Route 3'},'Location','southwest')
% xlabel('x')
% ylabel('Choice Probabilities P_i')

%% GPSL' Model

clear all
clc 

% Based on Braess_B.pdf
% Route 1  = a1 + a2
% Route 2 = a3 + a4
% Route 3 = a1 + a5 + a4
% a1 = a4 = x, a2 = a3 = 1, a5 = 0
% Letting x vary from 0 to 3.
% If x = 0.5, c3 = 1 < 1.5 = c1 = c2
% If x = 1, c3 = 2 = c1 = c2
% If x = 1.5, c3 = 3 > 2.5 = c1 = c2
% If x = 3, c3 = 6 > 4 = c1 = c2
x = 0:0.1:3.0;
x(1) = 10^(-16);
% Kappa is the penalty term for the log of outgoing edge set cardinality
kap = 1;
% Theta is the parameter on the cost function
theta = 1;

% Route costs
c1 = x+1;
c2 = 1+x;
c3 = x+0+x;

% Route utilities
u1 = -theta*c1;
u2 = -theta*c2;
u3 = -theta*c3;

% Gamma terms
% Travel time for each edge a
ta1 = x;
ta2 = 1;
ta3 = 1;
ta4 = x;
ta5 = 0;
% Incidence for each edge in a route
lambda = theta; 
delta1 = exp(-lambda.*c1)+exp(-lambda.*c3);
delta2 = exp(-lambda.*c1);
delta3 = exp(-lambda.*c2);
delta4 = exp(-lambda.*c2)+exp(-lambda.*c3);
delta5 = exp(-lambda.*c3);


Gamma1=(ta1./c1).*(1./((exp(lambda.*c1)).*delta1))+(ta2./c1).*(1./((exp(lambda.*c1)).*delta2));
Gamma2=(ta3./c2).*(1./((exp(lambda.*c2)).*delta3))+(ta4./c2).*(1./((exp(lambda.*c2)).*delta4));
Gamma3=(ta1./c3).*(1./((exp(lambda.*c3)).*delta1))+(ta5./c3).*(1./((exp(lambda.*c3)).*delta5))+(ta4./c3).*(1./((exp(lambda.*c3)).*delta4));


% Route probabilities
P1 = exp(u1)./(((Gamma1./Gamma1).^kap).*exp(u1)+((Gamma2./Gamma1).^kap).*exp(u2)+((Gamma3./Gamma1).^kap).*exp(u3));
P2 = exp(u2)./(((Gamma1./Gamma2).^kap).*exp(u1)+((Gamma2./Gamma2).^kap).*exp(u2)+((Gamma3./Gamma2).^kap).*exp(u3));
P3 = exp(u3)./(((Gamma1./Gamma3).^kap).*exp(u1)+((Gamma2./Gamma3).^kap).*exp(u2)+((Gamma3./Gamma3).^kap).*exp(u3));
P = P1+P2+P3;

% Welfare calculations
Waltgpsl = log((exp(u1).*((Gamma1).^kap)) + (exp(u2).*((Gamma2).^kap)) + (exp(u3).*((Gamma3).^kap)));
save('Figure10part2.mat','Waltgpsl','-append')

% % Plots route choice probabilities
% figure
% p = plot(x,P1,x,P2,x,P3);
% p(1).LineWidth = 2;
% p(2).LineWidth = 2;
% p(2).LineStyle = ':';
% p(3).LineWidth = 2;
% legend({'Route 1','Route 2','Route 3'},'Location','southwest')
% xlabel('x')
% ylabel('Choice Probabilities P_i')

%% Adaptive PSL Model

clear all
clc 

% Based on Braess_B.pdf
% Route 1  = a1 + a2
% Route 2 = a3 + a4
% Route 3 = a1 + a5 + a4
% a1 = a4 = x, a2 = a3 = 1, a5 = 0
% Letting x vary from 0.5 to 3.
% If x = 0.5, c3 = 1 < 1.5 = c1 = c2
% If x = 1, c3 = 2 = c1 = c2
% If x = 1.5, c3 = 3 > 2.5 = c1 = c2
% If x = 3, c3 = 6 > 4 = c1 = c2
% Can't use x here because I use it as a variable for the functions defined
% below
ex = 10^(-16); %as an approximation for zero
% Beta is the penalty term for gamma
beta = 1;
% Theta is the parameter on the cost function
theta = 1;

% Initialization
tau = 10^(-16);
N = 3;

i = 1;
error = 8; %% Duncan et al. (2020) states "predetermined convergence parameter"; cut error in half for easier convergence
I = 1000;
j = 1;

for j=1:length(0:0.1:3.0)
% Cost of routes
c(1) = ex+1;
c(2) = 1+ex;
c(3) = ex+ex;

% Route utilities
u(1) = -theta*c(1);
u(2) = -theta*c(2);
u(3) = -theta*c(3);

% Gamma terms
% Travel time for each edge a
ta(1) = ex;
ta(2) = 1;
ta(3) = 1;
ta(4) = ex;
ta(5) = 0;

% Initial guess
P0(1:3) = 1/N;


syms gamma1(x,y,z) gamma2(x,y,z) gamma3(x,y,z) g1(x,y,z) g2(x,y,z) g3(x,y,z) G1(x,y,z) G2(x,y,z) G3(x,y,z)

gamma1(x,y,z) = (ta(1)/c(1))*(x/(x+z))+(ta(2)/c(1))*(x/(x));
gamma2(x,y,z) = (ta(3)/c(2))*(y/(y))+(ta(4)/c(2))*(y/(y+z));
gamma3(x,y,z) = (ta(1)/c(3))*(z/(x+z))+(ta(5)/c(3))*(z/(z))+(ta(4)/c(3))*(z/(y+z));

g1(x,y,z) = (((gamma1(x,y,z))^beta)*exp(-theta*c(1)))/(((gamma1(x,y,z))^beta)*exp(-theta*c(1))+((gamma2(x,y,z))^beta)*exp(-theta*c(2))+((gamma3(x,y,z))^beta)*exp(-theta*c(3)));
g2(x,y,z) = (((gamma2(x,y,z))^beta)*exp(-theta*c(2)))/(((gamma1(x,y,z))^beta)*exp(-theta*c(1))+((gamma2(x,y,z))^beta)*exp(-theta*c(2))+((gamma3(x,y,z))^beta)*exp(-theta*c(3)));
g3(x,y,z) = (((gamma3(x,y,z))^beta)*exp(-theta*c(3)))/(((gamma1(x,y,z))^beta)*exp(-theta*c(1))+((gamma2(x,y,z))^beta)*exp(-theta*c(2))+((gamma3(x,y,z))^beta)*exp(-theta*c(3)));

G1(x,y,z) = tau + (1-N*tau)*g1(x,y,z);
G2(x,y,z) = tau + (1-N*tau)*g2(x,y,z);
G3(x,y,z) = tau + (1-N*tau)*g3(x,y,z);





    
while i <= I
    
    P(1) = G1(P0(1),P0(2),P0(3));
    P(2) = G2(P0(1),P0(2),P0(3));
    P(3) = G3(P0(1),P0(2),P0(3));
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
Wadaptive(j) = log((exp(u(1)).*((double(gamma1(P(1),P(2),P(3)))).^beta)) + (exp(u(2)).*((double(gamma2(P(1),P(2),P(3)))).^beta)) + (exp(u(3)).*((double(gamma3(P(1),P(2),P(3))).^beta))));
% j = j+1;
ex = ex + 0.1;

end

save('Figure10part2.mat','Wadaptive','-append')

% % plots route choice probabilities
% figure
% r = plot(0:length(ChoiceProbabilities)-1,ChoiceProbabilities);
% r(1).LineWidth = 2;
% r(2).LineWidth = 2;
% r(2).LineStyle = ':';
% r(3).LineWidth = 2;
% legend({'Route 1','Route 2','Route 3'},'Location','southwest')
% xlabel('x')
% xticklabels({'0','0.5','1','1.5','2','2.5','3'})
% ylabel('Choice Probability P_i')



%% Welfare (plots just the direct welfare calculation)
% clear
% load Figure10part2.mat 
% 
% figure
% q = plot(x,Wca,x,Wca5,x,Wmnl,x,Wpsl,x,Waltpsl,x,Wgpsl,x,Waltgpsl,x,Wadaptive);
% ylim([-10 5])
% q(1).LineWidth = 2;
% q(2).LineWidth = 2;
% q(3).LineWidth = 2;
% q(4).LineWidth = 2;
% q(5).LineWidth = 2;
% q(6).LineWidth = 2;
% q(7).LineWidth = 2;
% q(8).LineWidth = 2;
% q(8).LineStyle = ':';
% legend({'Choice Aversion, \kappa = 1','Choice Aversion, \kappa = 2','MNL','PSL','Alt. PSL','GPSL, \lambda = 3','Alt. GPSL, \lambda = \theta = 1','Adaptive PSL'},'Location','southwest')
% xlabel('x')
% ylabel('Welfare')

%% Difference in Welfare
clear
load Figure10part1.mat
load Figure10part2.mat

% Differences
Difca = Wca - Vca;
Difca5 = Wca5 - Vca5;
Difmnl = Wmnl - Vmnl;
Difpsl = Wpsl - Vpsl;
Difaltpsl = Waltpsl - Valtpsl;
Difgpsl = Wgpsl - Vgpsl;
Difaltgpsl = Waltgpsl - Valtgpsl;
Difadaptive = Wadaptive - Vadaptive;

% Plots Figure 10
figure
r = plot(x,Difca,x,Difca5,x,Difmnl,x,Difpsl,x,Difaltpsl,x,Difgpsl,x,Difaltgpsl,x,Difadaptive);
ylim([-1 1])
r(1).LineWidth = 2;
r(2).LineWidth = 2;
r(3).LineWidth = 2;
r(4).LineWidth = 2;
r(5).LineWidth = 2;
r(6).LineWidth = 2;
r(7).LineWidth = 2;
r(8).LineWidth = 2;
r(8).LineStyle = ':';
legend({'Choice Aversion, \kappa = 1','Choice Aversion, \kappa = 2','MNL','PSL','Alt. PSL','GPSL, \lambda = 3','Alt. GPSL, \lambda = \theta = 1','Adaptive PSL, \beta = \theta = 1'},'Location','southwest')
xlabel('x')
ylabel('Difference in Welfare')