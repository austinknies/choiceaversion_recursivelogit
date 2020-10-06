%% Figure 5(A) from Knies and Melo (2020)
% Calculates route choice probabilities for Figure 4 
% Network originally featured in Fosgerau et al. (2013) (Figure 2)
% Generates Figure 5(A)
% Choice Aversion Model
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

% Kappa is the penalty term for the log of outgoing edge set cardinality
kap = 0:0.1:10;
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
% For simplicity of code, here Aji = product of choice set cardinality at
% each node along route i in the network. 
Aj1=2*2*2*1;
Aj2=2*2*2*1;
Aj3=2*2*1*1;
Aj4=2*1; 


% Route probabilities
P1 = exp(u1)./(((Aj1/Aj1).^kap).*exp(u1)+((Aj1/Aj2).^kap).*exp(u2)+((Aj1/Aj3).^kap).*exp(u3)+((Aj1/Aj4).^kap).*exp(u4));
P2 = exp(u2)./(((Aj2/Aj1).^kap).*exp(u1)+((Aj2/Aj2).^kap).*exp(u2)+((Aj2/Aj3).^kap).*exp(u3)+((Aj2/Aj4).^kap).*exp(u4));
P3 = exp(u3)./(((Aj3/Aj1).^kap).*exp(u1)+((Aj3/Aj2).^kap).*exp(u2)+((Aj3/Aj3).^kap).*exp(u3)+((Aj3/Aj4).^kap).*exp(u4));
P4 = exp(u4)./(((Aj4/Aj1).^kap).*exp(u1)+((Aj4/Aj2).^kap).*exp(u2)+((Aj4/Aj3).^kap).*exp(u3)+((Aj4/Aj4).^kap).*exp(u4));
P = P1+P2+P3+P4; % unused; just one easy check to help ensure you haven't messed up the choice probabilities

% % Welfare calculation (unused here, but easy to calculate)
% Wca = log((exp(u1)./((Aj1).^kap)) + (exp(u2)./((Aj2).^kap)) + (exp(u3)./((Aj3).^kap))+ (exp(u4)./((Aj4).^kap)));


%Figure 5(A)
figure
p = plot(kap,P1,kap,P2,kap,P3,kap,P4);
p(1).LineWidth = 2;
p(2).LineWidth = 2;
p(3).LineWidth = 2;
p(2).LineStyle = ':';
p(4).LineWidth = 2;
legend({'Route 1','Route 2','Route 3','Route 4'},'Location','northeast')
xlabel('\kappa')
ylabel('Choice Probabilities P_i')

%% Figure 5(B) from Knies and Melo (2020)
% Calculates route choice probabilities for Figure 4 
% Network originally featured in Fosgerau et al. (2013) (Figure 2)
% Generates Figure 5(B)
% Adaptive PSL Model 

 alreadydone = false; % Make sure to reverse the commenting if you want
% this code to not run, isn't necessarily quick
%alreadydone = true;

% Route 1 : 12,23,35
% Route 2 : 12,23,34,45
% Route 3 : 12,24,45
% Route 4 : 15
if alreadydone~=true
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

% Theta is the parameter on the cost function
theta = 1;

% Route utilities
u1 = -theta*c1;
u2 = -theta*c2;
u3 = -theta*c3;
u4 = -theta*c4;

% Beta is the penalty term for link-path incidence (gamma), will go from 0 to 10
beta = 0;

% Initialization
tau = 10^(-16);
N = 4;

i = 1;
error = 16; % Duncan et al. (2020) states "predetermined convergence parameter"
I = 10000;
% Preallocation for speed
ChoiceProbabilities=zeros(length(0:0.1:10),N);
%Wadaptive=zeros(length(0:0.1:10)); %this is if you want welfare
%                                     calculations along the way
for j=1:length(0:0.1:10)

% Initial guess
P0(1:N) = 1/N;


syms gamma1(w,x,y,z) gamma2(w,x,y,z) gamma3(w,x,y,z) gamma4(w,x,y,z) g1(w,x,y,z) g2(w,x,y,z) g3(w,x,y,z) g4(w,x,y,z) G1(w,x,y,z) G2(w,x,y,z) G3(w,x,y,z) G4(w,x,y,z)

gamma1(w,x,y,z) = (t12/c1)*(w/(w+x+y))+(t23/c1)*(w/(w+x))+(t35/c1)*(w/(w));
gamma2(w,x,y,z) = (t12/c2)*(x/(w+x+y))+(t23/c2)*(x/(w+x))+(t34/c2)*(x/(x))+(t45/c2)*(x/(x+y));
gamma3(w,x,y,z) = (t12/c3)*(y/(w+x+y))+(t24/c3)*(y/(y))+(t45/c3)*(y/(x+y));
gamma4(w,x,y,z) = (t15/c4)*(z/z);

g1(w,x,y,z) = (((gamma1(w,x,y,z))^beta)*exp(u1))/(((gamma1(w,x,y,z))^beta)*exp(u1)+((gamma2(w,x,y,z))^beta)*exp(u2)+((gamma3(w,x,y,z))^beta)*exp(u3)+((gamma4(w,x,y,z))^beta)*exp(u4));
g2(w,x,y,z) = (((gamma2(w,x,y,z))^beta)*exp(u2))/(((gamma1(w,x,y,z))^beta)*exp(u1)+((gamma2(w,x,y,z))^beta)*exp(u2)+((gamma3(w,x,y,z))^beta)*exp(u3)+((gamma4(w,x,y,z))^beta)*exp(u4));
g3(w,x,y,z) = (((gamma3(w,x,y,z))^beta)*exp(u3))/(((gamma1(w,x,y,z))^beta)*exp(u1)+((gamma2(w,x,y,z))^beta)*exp(u2)+((gamma3(w,x,y,z))^beta)*exp(u3)+((gamma4(w,x,y,z))^beta)*exp(u4));
g4(w,x,y,z) = (((gamma4(w,x,y,z))^beta)*exp(u4))/(((gamma1(w,x,y,z))^beta)*exp(u1)+((gamma2(w,x,y,z))^beta)*exp(u2)+((gamma3(w,x,y,z))^beta)*exp(u3)+((gamma4(w,x,y,z))^beta)*exp(u4));

G1(w,x,y,z) = tau + (1-N*tau)*g1(w,x,y,z);
G2(w,x,y,z) = tau + (1-N*tau)*g2(w,x,y,z);
G3(w,x,y,z) = tau + (1-N*tau)*g3(w,x,y,z);
G4(w,x,y,z) = tau + (1-N*tau)*g4(w,x,y,z);




    
while i <= I
    
    P(1) = G1(P0(1),P0(2),P0(3),P0(4));
    P(2) = G2(P0(1),P0(2),P0(3),P0(4));
    P(3) = G3(P0(1),P0(2),P0(3),P0(4));
    P(4) = G4(P0(1),P0(2),P0(3),P0(4));
    P=double(P);
    pcheck = log(sum(abs(P-P0)));
    ptest = log(10^(-error));
    if pcheck<ptest
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
%Wadaptive(j) = log((exp(u1).*((double(gamma1(P(1),P(2),P(3),P(4)))).^beta)) + (exp(u2).*((double(gamma2(P(1),P(2),P(3),P(4)))).^beta)) + (exp(u3).*((double(gamma3(P(1),P(2),P(3),P(4))).^beta))) + (exp(u4).*((double(gamma1(P(1),P(2),P(3),P(4)))).^beta)));
% j = j+1;
beta=beta+0.1;
fprintf('beta = %3.2f \n',beta);
end


%Figure
figure
r = plot(0:length(ChoiceProbabilities)-1,ChoiceProbabilities);
r(1).LineWidth = 2;
r(2).LineWidth = 2;
r(3).LineStyle = ':';
r(3).LineWidth = 2;
r(4).LineWidth = 2;
legend({'Route 1','Route 2','Route 3','Route 4'},'Location','northwest')
xlabel('\beta')
xticklabels({'0','2','4','6','8','10'})
ylabel('Choice Probability P_i')

end


