%% Figure 9(B) from Knies and Melo (2020)
% Calculates welfare for Figures 7(A) and 7(B) using the 
% APSL model for varying kappa at particular values of x

%% Welfare for 7(A)
% Using APSL Model
clear all
clc 

% Based on Braess_A.pdf
% Route 1  = a1 + a2
% Route 2 = a3 + a4
% a1 = a4 = x, a2 = a3 = 1
% Letting x vary from 0 to 3.
% Can't use x here because I use it as a variable for the functions defined
% below

% Beta is the penalty term for gamma
beta = 0;
% Theta is the parameter on the cost function
theta = 1;
Vadaptive = zeros(length(0:0.5:3.0),length(0:0.1:10));

% Initialization
tau = 10^(-16);
N = 2;

i = 1;
error = 16; % "predetermined convergence rate"
I = 10000;
j = 1;

for k=1:length(0:0.1:10)
ex = 10^(-16); %as an approximation for zero
for j=1:length(0:0.5:3.0)
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
Vadaptive(j,k) = log((exp(u(1)).*((double(gamma1(P(1),P(2)))).^beta)) + (exp(u(2)).*((double(gamma2(P(1),P(2)))).^beta)));
% j = j+1;
ex = ex + 0.5;

fprintf('ex = %3.2f\n',ex);
end
beta=beta+0.1;
fprintf('beta = %3.2f \n',beta);
end


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

save Figure9B.mat Vadaptive

%% Welfare for Figure 7(B)
% Using APSL Model
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

% Beta is the penalty term for gamma
beta = 0;
% Theta is the parameter on the cost function
theta = 1;
Wadaptive=zeros(length(0:0.5:3.0),length(0:0.1:10));
% Initialization
tau = 10^(-16);
N = 3;

i = 1;
error = 16; % "predetermined convergence parameter"
I = 10000;
j = 1;
for k=1:length(0:0.1:10)
    ex = 10^(-16); %as an approximation for zero
for j=1:length(0:0.5:3.0)
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
Wadaptive(j,k) = log((exp(u(1)).*((double(gamma1(P(1),P(2),P(3)))).^beta)) + (exp(u(2)).*((double(gamma2(P(1),P(2),P(3)))).^beta)) + (exp(u(3)).*((double(gamma3(P(1),P(2),P(3))).^beta))));
% j = j+1;
ex = ex + 0.5;
fprintf('ex = %3.2f\n',ex);
end
beta=beta+0.1;
fprintf('beta = %3.2f \n',beta);
end


% % Plots route choice probabilities
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

save('Figure9B.mat','Wadaptive','-append')

%% Welfare Change
load Figure9B.mat
Difference = Wadaptive - Vadaptive;
PercentChange = Difference./abs(Vadaptive);
kap=0:0.1:10; %didn't have a kappa to save this time; this is just for the plot

% Figure 9(B)
figure
r = plot(kap,Difference(1,:),kap,Difference(2,:),kap,Difference(3,:),kap,Difference(4,:),kap,Difference(5,:),kap,Difference(6,:),kap,Difference(7,:));
r(1).LineWidth = 2;
r(2).LineWidth = 2;
r(3).LineWidth = 2;
r(4).LineWidth = 2;
r(5).LineWidth = 2;
r(6).LineWidth = 2;
r(7).LineWidth = 2;
legend({'x=0','x=0.5','x=1','x=1.5','x=2','x=2.5','x=3'},'Location','northeast')
xlabel('\beta')
ylabel('Difference in Welfare')


