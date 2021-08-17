%% Table 1 (Columns 2 and 3) from Knies, Lorca, and Melo (2021)
% Calculates CA model route choice probabilities for Figure 6
% Network originally featured in Mai et al. (2015), Figure 4
% Last Updated: August 2021

% To get column 2, keep simple==0, ensure lines 51-55 are commented
% To get column 3, keep simple==0, uncomment lines 51-55

clear all
clc 

simple = 0; % uncomment to run base version of CA (Table 1, Column 2 or 3)
% simple = 1; % uncomment for plot as mu_a and mu_b vary


if simple == 0


% Instantaneous edge costs
o = 1;
a = 1;
b = 2;
a1 = 1;
a2 = 2;
a3 = 2;
b1 = 1;
b2 = 1.5;
b3 = 1;
e = 1;
d = 0;
f = 2; % Not used as a path, but included in the overlapping calculation

% Parameters that are different: mu_a and mu_b 
mua = 0.5; % Mai et al. (2015) calibration is mua = 0.5, mub = 0.8
mub = 0.8; %
% all other scale parameters set to 1
muo = 1;
mua1 = 1;
mua2 = 1;
mua3 = 1;
mub1 = 1;
mub2 = 1;
mub3 = 1;
mue = 1;
muf = 1; % note: same as muo

% Choice aversion parameters
% Here, we set kappa = -omega_OL from Mai et al. (2015)'s estimation of mu
kapa = -(log(muo))./2;
kapb = -(log(mua))./3;
kapc = -(log(mub))./3;
kape = -(log(mua3))./2; % could also replace mua3 for mub1
kapd = -(log(mue))./2; % could replace with mua1/2 or mub2/3
kap = 1; % in heterogeneous case, we code kappa differently, so current kap in choice probabilities needs to be neutralized
theta = 1; % unused in our model, but just an edge cost scale parameter
% used more in other PSL models

% kapa = 1; % uncomment for kappa = 1 (baseline case)
% kapb = 1;
% kapc = 1;
% kape = 1; 
% kapd = 1;

% Route utilities (minimizing costs = maximizing negative costs)
u1 = -theta*(o+a+a1+d);
u2 = -theta*(o+a+a2+d);
u3 = -theta*(o+a+a3+e+d);
u4 = -theta*(o+b+b1+e+d);
u5 = -theta*(o+b+b2+d);
u6 = -theta*(o+b+b3+d);

% Choice set size for each route
% Note Aji =\= the notation in the paper. 
% For simplicity of code, here Aji = product of choice set cardinality at
% each node along route i in the network. 
Aj1=(2.^kapa).*((exp(kapb.*log(3)))).*(1.^kapd);
Aj2=(2.^kapa).*((exp(kapb.*log(3)))).*(1.^kapd);
Aj3=(2.^kapa).*((exp(kapb.*log(3)))).*(2.^kape)*(1.^kapd);
Aj4=(2.^kapa).*((exp(kapc.*log(3)))).*(2.^kape)*(1.^kapd);
Aj5=(2.^kapa).*((exp(kapc.*log(3)))).*(1.^kapd);
Aj6=(2.^kapa).*((exp(kapc.*log(3)))).*(1.^kapd);


% Route probabilities
P1 = exp(u1)./(((Aj1./Aj1).^kap).*exp(u1)+((Aj1./Aj2).^kap).*exp(u2)+((Aj1./Aj3).^kap).*exp(u3)+((Aj1./Aj4).^kap).*exp(u4)+((Aj1./Aj5).^kap).*exp(u5)+((Aj1./Aj6).^kap).*exp(u6));
P2 = exp(u2)./(((Aj2./Aj1).^kap).*exp(u1)+((Aj2./Aj2).^kap).*exp(u2)+((Aj2./Aj3).^kap).*exp(u3)+((Aj2./Aj4).^kap).*exp(u4)+((Aj2./Aj5).^kap).*exp(u5)+((Aj2./Aj6).^kap).*exp(u6));
P3 = exp(u3)./(((Aj3./Aj1).^kap).*exp(u1)+((Aj3./Aj2).^kap).*exp(u2)+((Aj3./Aj3).^kap).*exp(u3)+((Aj3./Aj4).^kap).*exp(u4)+((Aj3./Aj5).^kap).*exp(u5)+((Aj3./Aj6).^kap).*exp(u6));
P4 = exp(u4)./(((Aj4./Aj1).^kap).*exp(u1)+((Aj4./Aj2).^kap).*exp(u2)+((Aj4./Aj3).^kap).*exp(u3)+((Aj4./Aj4).^kap).*exp(u4)+((Aj4./Aj5).^kap).*exp(u5)+((Aj4./Aj6).^kap).*exp(u6));
P5 = exp(u5)./(((Aj5./Aj1).^kap).*exp(u1)+((Aj5./Aj2).^kap).*exp(u2)+((Aj5./Aj3).^kap).*exp(u3)+((Aj5./Aj4).^kap).*exp(u4)+((Aj5./Aj5).^kap).*exp(u5)+((Aj5./Aj6).^kap).*exp(u6));
P6 = exp(u6)./(((Aj6./Aj1).^kap).*exp(u1)+((Aj6./Aj2).^kap).*exp(u2)+((Aj6./Aj3).^kap).*exp(u3)+((Aj6./Aj4).^kap).*exp(u4)+((Aj6./Aj5).^kap).*exp(u5)+((Aj6./Aj6).^kap).*exp(u6));
P = P1+P2+P3+P4+P5+P6;
end
%% Varying mua and mub
% Calculates NRL model route choice probabilities for Figure 6
% Network originally featured in Mai et al. (2015), Figure 4
% Last Updated: August 2021
if simple==1

% Instantaneous edge costs
o = 1;
a = 1;
b = 2;
a1 = 1;
a2 = 2;
a3 = 2;
b1 = 1;
b2 = 1.5;
b3 = 1;
e = 1;
d = 0;
f = 2; % Not used as a path, but included in the overlapping calculation

% all main scale parameters set to 1
muo = 1;
mua1 = 1;
mua2 = 1;
mua3 = 1;
mub1 = 1;
mub2 = 1;
mub3 = 1;
mue = 1;
muf = 1; % note: same as muo

% Parameters that are different: mu_a and mu_b 
mua = 0.01; % Mai et al. (2015) calibration is mua = 0.5, mub = 0.8
    j = 1; % iteration for outer loop
while mua <=1
    mub = 0.01; %
        k = 1; % iteration for inner loop
    while mub<=1
    % Choice aversion parameters
    % Here, we set kappa = -omega_OL from Mai et al. (2015)'s estimation of mu
    kapa = -(log(muo))./2;
    kapb = -(log(mua))./3;
    kapc = -(log(mub))./3;
    kape = -(log(mua3))./2; % could also replace mua3 for mub1
    kapd = -(log(mue))./2; % could replace with mua1/2 or mub2/3
    kap = 1; % in heterogeneous case, we code kappa differently, so current kap in choice probabilities needs to be neutralized
    theta = 1; % unused in our model, but just an edge cost scale parameter
    % used more in other PSL models

    % Route utilities (minimizing costs = maximizing negative costs)
    u1 = -theta*(o+a+a1+d);
    u2 = -theta*(o+a+a2+d);
    u3 = -theta*(o+a+a3+e+d);
    u4 = -theta*(o+b+b1+e+d);
    u5 = -theta*(o+b+b2+d);
    u6 = -theta*(o+b+b3+d);

    % Choice set size for each route
    % Note Aji =\= the notation in the paper. 
    % For simplicity of code, here Aji = product of choice set cardinality at
    % each node along route i in the network. 
    Aj1=(2^kapa)*((exp(kapb*log(3))))*(1^kapd);
    Aj2=(2^kapa)*((exp(kapb*log(3))))*(1^kapd);
    Aj3=(2^kapa)*((exp(kapb*log(3))))*(2^kape)*(1^kapd);
    Aj4=(2^kapa)*((exp(kapc*log(3))))*(2^kape)*(1^kapd);
    Aj5=(2^kapa)*((exp(kapc*log(3))))*(1^kapd);
    Aj6=(2^kapa)*((exp(kapc*log(3))))*(1^kapd);


    % Route probabilities
    P1(j,k) = exp(u1)./((Aj1./Aj1).*exp(u1)+(Aj1./Aj2).*exp(u2)+(Aj1./Aj3).*exp(u3)+(Aj1./Aj4).*exp(u4)+(Aj1./Aj5).*exp(u5)+(Aj1./Aj6).*exp(u6));
    P2(j,k) = exp(u2)./(((Aj2./Aj1).^kap).*exp(u1)+((Aj2./Aj2).^kap).*exp(u2)+((Aj2./Aj3).^kap).*exp(u3)+((Aj2./Aj4).^kap).*exp(u4)+((Aj2./Aj5).^kap).*exp(u5)+((Aj2./Aj6).^kap).*exp(u6));
    P3(j,k) = exp(u3)./(((Aj3./Aj1).^kap).*exp(u1)+((Aj3./Aj2).^kap).*exp(u2)+((Aj3./Aj3).^kap).*exp(u3)+((Aj3./Aj4).^kap).*exp(u4)+((Aj3./Aj5).^kap).*exp(u5)+((Aj3./Aj6).^kap).*exp(u6));
    P4(j,k) = exp(u4)./(((Aj4./Aj1).^kap).*exp(u1)+((Aj4./Aj2).^kap).*exp(u2)+((Aj4./Aj3).^kap).*exp(u3)+((Aj4./Aj4).^kap).*exp(u4)+((Aj4./Aj5).^kap).*exp(u5)+((Aj4./Aj6).^kap).*exp(u6));
    P5(j,k) = exp(u5)./(((Aj5./Aj1).^kap).*exp(u1)+((Aj5./Aj2).^kap).*exp(u2)+((Aj5./Aj3).^kap).*exp(u3)+((Aj5./Aj4).^kap).*exp(u4)+((Aj5./Aj5).^kap).*exp(u5)+((Aj5./Aj6).^kap).*exp(u6));
    P6(j,k) = exp(u6)./(((Aj6./Aj1).^kap).*exp(u1)+((Aj6./Aj2).^kap).*exp(u2)+((Aj6./Aj3).^kap).*exp(u3)+((Aj6./Aj4).^kap).*exp(u4)+((Aj6./Aj5).^kap).*exp(u5)+((Aj6./Aj6).^kap).*exp(u6));
    P(j,k) = P1(j,k)+P2(j,k)+P3(j,k)+P4(j,k)+P5(j,k)+P6(j,k);
    
    mub = mub + 0.01
    k=k+1;
    end

mua = mua + 0.01
j=j+1;

end

% Figure for varying mu (high mu, low kappa)
figure
surf(P1,'FaceColor','r','EdgeAlpha', 0.2)
alpha 0.4
hold on
surf(P2,'FaceColor','g','EdgeAlpha', 0.2)
alpha 0.4
hold on
surf(P3,'FaceColor','b','EdgeAlpha', 0.2)
alpha 0.4
hold off
% colormap hsv
% shading flat
% alphamap('vup')
legend({'Route 1','Route 2','Route 3'},'Location','northwest')
ylabel('\mu_a')
yticks([1 50 100])
yticklabels({'0.01','0.5','1'})
xlabel('\mu_b')
xticks([1 50 100])
xticklabels({'0.01','0.5','1'})
zlabel('Path Choice Probability')
title('Choice Aversion Model - Paths Using Node B')
grid on


figure
surf(P4,'FaceColor','c','EdgeAlpha', 0.2)
alpha 0.4
hold on
surf(P5,'FaceColor','m','EdgeAlpha', 0.2)
alpha 0.4
hold on
surf(P6,'FaceColor','y','EdgeAlpha', 0.2)
alpha 0.4
hold off
% colormap hsv
% shading flat
legend({'Route 4','Route 5','Route 6'},'Location','northwest')
ylabel('\mu_a')
yticks([1 50 100])
yticklabels({'0.01','0.5','1'})
xlabel('\mu_b')
xticks([1 50 100])
xticklabels({'0.01','0.5','1'})
zlabel('Path Choice Probability')
title('Choice Aversion Model - Paths Using Node C')
grid on
end