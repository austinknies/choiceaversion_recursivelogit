%% Figure 11(A) from Knies, Lorca, and Melo (2021)
% Calculates welfare for Figures 8(A) and 8(B) using the choice aversion
% model for varying kappa at particular values of x
% Last Updated: August 2021

%% Welfare for Figure 8(A)
% Using Choice Aversion Model
clear all
clc 

% Based on Braess_A.pdf
% Route 1  = a1 + a2
% Route 2 = a3 + a4
% a1 = a4 = x, a2 = a3 = 1

% Kappa is the penalty term for the log of outgoing edge set cardinality
kap = 0:0.1:10;
sizekap = length(kap);
samplex = 0:0.5:3.0;
sizex = length(samplex);
Welfarebefore = zeros(sizex,sizekap);
% Theta is the parameter on the cost function
theta = 1;
% Letting x vary from 0 to 3.
x=0;
i=1;
while x<=3.0


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
P1 = exp(u1)./(((Aj1/Aj1).^kap).*exp(u1)+((Aj1/Aj2).^kap).*exp(u2));
P2 = exp(u2)./(((Aj2/Aj1).^kap).*exp(u1)+((Aj2/Aj2).^kap).*exp(u2));
P = P1+P2;

% Welfare calculation
Welfarebefore(i,:) = log((exp(u1)./((Aj1).^kap)) + (exp(u2)./((Aj2).^kap)));
x = x+0.5;
i = i+1;
end
% % plots route choice probabilities
% figure
% p = plot(x,P1,x,P2);
% p(1).LineWidth = 2;
% p(2).LineWidth = 2;
% p(2).LineStyle = ':';
% legend({'Route 1','Route 2'},'Location','southwest')
% xlabel('x')
% ylabel('Choice Probabilities P_i')

save Figure11A.mat Welfarebefore

%% Welfare for Figure 8(B)
% Using Choice Aversion Model
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
% Kappa is the penalty term for the log of outgoing edge set cardinality
kap = 0:0.1:10;
sizekap = length(kap);
samplex = 0:0.5:3.0;
sizex = length(samplex);
Welfareafter = zeros(sizex,sizekap);
% Theta is the parameter on the cost function
theta = 1;
% Letting x vary from 0 to 3.
x=0;
i=1;
while x<=3.0

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
P1 = exp(u1)./(((Aj1/Aj1).^kap).*exp(u1)+((Aj1/Aj2).^kap).*exp(u2)+((Aj1/Aj3).^kap).*exp(u3));
P2 = exp(u2)./(((Aj2/Aj1).^kap).*exp(u1)+((Aj2/Aj2).^kap).*exp(u2)+((Aj2/Aj3).^kap).*exp(u3));
P3 = exp(u3)./(((Aj3/Aj1).^kap).*exp(u1)+((Aj3/Aj2).^kap).*exp(u2)+((Aj3/Aj3).^kap).*exp(u3));
P = P1+P2+P3;

% Welfare calculation
Welfareafter(i,:) = log((exp(u1)./((Aj1).^kap)) + (exp(u2)./((Aj2).^kap)) + (exp(u3)./((Aj3).^kap)));
x = x+0.5;
i = i+1;
end
% % plots route choice probabilities
% figure
% p = plot(x,P1,x,P2,x,P3);
% p(1).LineWidth = 2;
% p(2).LineWidth = 2;
% p(3).LineWidth = 2;
% legend({'Route 1','Route 2','Route 3'},'Location','southwest')
% xlabel('x')
% ylabel('Choice Probabilities P_i')

save('Figure11A.mat','Welfareafter','-append')

%% Welfare Change
load Figure11A.mat
Difference = Welfareafter - Welfarebefore;
PercentChange = Difference./abs(Welfarebefore);

% Figure 11(A)
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
xlabel('\kappa')
ylabel('Difference in Welfare')


