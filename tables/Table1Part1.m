%% Table 1 (Columns 1) from Knies, Lorca, and Melo (2021)
% Calculates NRL model route choice probabilities for Figure 6
% Network originally featured in Mai et al. (2015), Figure 4
% Last Updated: August 2021

clear all
clc

simple = 0; % uncomment to run base version of NRL (Table 1, Column 1)
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
f = 2;
d = 0; % just being thorough


% Parameters that are different: mu_a and mu_b 
mua = 0.5; % Mai et al. (2015) calibration
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

theta = 1; % unused in model, but just an edge cost scale parameter

% Route utilities (minimizing costs = maximizing negative costs)
va = -theta*(a);
vb = -theta*(b);
va1 = -theta*(a1);
va2 = -theta*(a2);
va3 = -theta*(a3);
vb1 = -theta*(b1);
vb2 = -theta*(b2);
vb3 = -theta*(b3);
ve = -theta*(e);
vf = -theta*(f);
vd = -theta*(d); % just being thorough


% Value functions not necessary for the loop
Vdd = 0;
Vda1 = mua1.*log(exp((1./mua1).*(vd+Vdd))); % just being thorough
Vda2 = mua2.*log(exp((1./mua2).*(vd+Vdd)));
Vdb2 = mub2.*log(exp((1./mub2).*(vd+Vdd)));
Vdb3 = mub3.*log(exp((1./mub3).*(vd+Vdd)));
Vde = mue.*log(exp((1./mue).*(vd+Vdd)));

% Fixed point loop for Vdf
i = 1;
error = 16; % "predetermined convergence rate"
I = 10000;
% Dealing with edge f, which creates a cycle. F is a potential option on 
% paths 3 and 4, so f's continuation value includes the entire network. 

% When muf = muo, then Vdf = Vdo. 
Vdf = 0; % Initial guess for Vdf

while i<=I

    % All other value functions depending on Vdf
    Vda3 = mua3.*log(exp((1./mua3).*(ve+Vde))+exp((1./mua3).*(vf+Vdf)));
    Vdb1 = mub1.*log(exp((1./mub1).*(ve+Vde))+exp((1./mub1).*(vf+Vdf)));
    Vda = mua.*log(exp((1./mua).*(va1+Vda1))+exp((1./mua).*(va2+Vda2))+exp((1./mua).*(va3+Vda3)));
    Vdb = mub.*log(exp((1./mub).*(vb1+Vdb1))+exp((1./mub).*(vb2+Vdb2))+exp((1./mub).*(vb3+Vdb3)));
    Vdo = muo.*log(exp((1./muo).*(va+Vda))+exp((1./muo).*(vb+Vdb)));
    
    % Choice prob numerators and denom
    
    % Path 1: o, a, a1, d
    num1 = va./muo + va1./mua + vd./mua1 + (mua./muo - 1).*(Vda./mua) + (mua1./mua - 1).*(Vda1./mua1) + Vdd./mua1;
    % Path 2: o, a, a2, d
    num2 = va./muo + va2./mua + vd./mua2 + (mua./muo - 1).*(Vda./mua) + (mua2./mua - 1).*(Vda2./mua2) + Vdd./mua2;
    % Path 3: o, a, a3, e, d
    num3 = va./muo + va3./mua + ve./mua3 + vd./mue + (mua./muo - 1).*(Vda./mua) + (mua3./mua - 1).*(Vda3./mua3) + (mue./mua3 - 1).*(Vde./mue) + Vdd./mue;
    % Path 4: o, b, b1, e, d
    num4 = vb./muo + vb1./mub + ve./mub1 + vd./mue + (mub./muo - 1).*(Vdb./mub) + (mub1./mub - 1).*(Vdb1./mub1) + (mue./mub1 - 1).*(Vde./mue) + Vdd./mue;
    % Path 5: o, b, b2, d
    num5 = vb./muo + vb2./mub + vd./mub2 + (mub./muo - 1).*(Vdb./mub) + (mub2./mub - 1).*(Vdb2./mub2) + Vdd./mub2;
    % Path 6: o, b, b3, d
    num6 = vb./muo + vb3./mub + vd./mub3 + (mub./muo - 1).*(Vdb./mub) + (mub3./mub - 1).*(Vdb3./mub3) + Vdd./mub3;
    
    
    
    den = Vdo./muo;

    % Route probabilities
    P1 = exp(num1)./exp(den);
    P2 = exp(num2)./exp(den);
    P3 = exp(num3)./exp(den);
    P4 = exp(num4)./exp(den);
    P5 = exp(num5)./exp(den);
    P6 = exp(num6)./exp(den);
    P = P1+P2+P3+P4+P5+P6;

%    if (log(abs(P-1))<log(10^(-error)))&&(log(abs(Vdo-Vdf))<log(10^(-error)))
   if (log(abs(Vdo-Vdf))<log(10^(-error)))
    %fprintf('Solution has been found!')
        break
   end
    i = i+1;
    Vdf = Vdo;
end
if i > I
    fprintf('Solution did not converge')
else
    fprintf('Solution has been found! Iterations: %d',i)
    P
end


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

theta = 1; % unused in model, but just an edge cost scale parameter


% Route utilities (minimizing costs = maximizing negative costs)
va = -theta*(a);
vb = -theta*(b);
va1 = -theta*(a1);
va2 = -theta*(a2);
va3 = -theta*(a3);
vb1 = -theta*(b1);
vb2 = -theta*(b2);
vb3 = -theta*(b3);
ve = -theta*(e);
vf = -theta*(f);
vd = -theta*(d); % just being thorough

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

% Value functions not necessary for the loop
Vdd = 0;
Vda1 = mua1.*log(exp((1./mua1).*(vd+Vdd))); % just being thorough
Vda2 = mua2.*log(exp((1./mua2).*(vd+Vdd)));
Vdb2 = mub2.*log(exp((1./mub2).*(vd+Vdd)));
Vdb3 = mub3.*log(exp((1./mub3).*(vd+Vdd)));
Vde = mue.*log(exp((1./mue).*(vd+Vdd)));

% Parameters that are different: mu_a and mu_b 
mua = 0.01; % Mai et al. (2015) calibration is mua = 0.5, mub = 0.8
    j = 1; % iteration for outer loop
while mua <=1
    mub = 0.01; %
        k = 1; % iteration for inner loop
    while mub<=1
        
 	      % Fixed point loop for Vdf
            i = 1;
            error = 16; % "predetermined convergence rate"
            I = 10000;
            % Dealing with edge f, which creates a cycle. F is a potential option on 
            % paths 3 and 4, so f's continuation value includes the entire network. 

            % When muf = muo, then Vdf = Vdo. 
            Vdf = 0; % Initial guess for Vdf
        
            while i<=I

            % All other value functions depending on Vdf
            Vda3 = mua3.*log(exp((1./mua3).*(ve+Vde))+exp((1./mua3).*(vf+Vdf)));
            Vdb1 = mub1.*log(exp((1./mub1).*(ve+Vde))+exp((1./mub1).*(vf+Vdf)));
            Vda = mua.*log(exp((1./mua).*(va1+Vda1))+exp((1./mua).*(va2+Vda2))+exp((1./mua).*(va3+Vda3)));
            Vdb = mub.*log(exp((1./mub).*(vb1+Vdb1))+exp((1./mub).*(vb2+Vdb2))+exp((1./mub).*(vb3+Vdb3)));
            Vdo = muo.*log(exp((1./muo).*(va+Vda))+exp((1./muo).*(vb+Vdb)));
        
                % Choice prob numerators and denom
    
                % Path 1: o, a, a1, d
                num1 = va./muo + va1./mua + vd./mua1 + (mua./muo - 1).*(Vda./mua) + (mua1./mua - 1).*(Vda1./mua1) + Vdd./mua1;
                % Path 2: o, a, a2, d
                num2 = va./muo + va2./mua + vd./mua2 + (mua./muo - 1).*(Vda./mua) + (mua2./mua - 1).*(Vda2./mua2) + Vdd./mua2;
                % Path 3: o, a, a3, e, d
                num3 = va./muo + va3./mua + ve./mua3 + vd./mue + (mua./muo - 1).*(Vda./mua) + (mua3./mua - 1).*(Vda3./mua3) + (mue./mua3 - 1).*(Vde./mue) + Vdd./mue;
                % Path 4: o, b, b1, e, d
                num4 = vb./muo + vb1./mub + ve./mub1 + vd./mue + (mub./muo - 1).*(Vdb./mub) + (mub1./mub - 1).*(Vdb1./mub1) + (mue./mub1 - 1).*(Vde./mue) + Vdd./mue;
                % Path 5: o, b, b2, d
                num5 = vb./muo + vb2./mub + vd./mub2 + (mub./muo - 1).*(Vdb./mub) + (mub2./mub - 1).*(Vdb2./mub2) + Vdd./mub2;
                % Path 6: o, b, b3, d
                num6 = vb./muo + vb3./mub + vd./mub3 + (mub./muo - 1).*(Vdb./mub) + (mub3./mub - 1).*(Vdb3./mub3) + Vdd./mub3;
    
    
    
                den = Vdo./muo;
                   
                
                
               % Route probabilities
               P1(j,k) = exp(num1)./exp(den);
               P2(j,k) = exp(num2)./exp(den);
               P3(j,k) = exp(num3)./exp(den);
               P4(j,k) = exp(num4)./exp(den);
               P5(j,k) = exp(num5)./exp(den);
               P6(j,k) = exp(num6)./exp(den);
               P(j,k) = P1(j,k)+P2(j,k)+P3(j,k)+P4(j,k)+P5(j,k)+P6(j,k);

            %    if (log(abs(P-1))<log(10^(-error)))&&(log(abs(Vdo-Vdf))<log(10^(-error)))
                if (log(abs(Vdo-Vdf))<log(10^(-error)))
                %fprintf('Solution has been found!')
                    break
                end
                i = i+1;
                Vdf = Vdo;
                end
                if i > I
                    fprintf('Solution did not converge. Iterations: %d',i, j, k)
%                 else
%                     fprintf('Solution has been found! Iterations: %d',i);
%                     P(j,k);
                end
    mub = mub + 0.01
    k=k+1;
    end

mua = mua + 0.01
j=j+1;

end

% Figure for varying mu 
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
title('NRL Model - Paths Using Edge a')
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
title('NRL Model - Paths Using Edge b')
grid on
end