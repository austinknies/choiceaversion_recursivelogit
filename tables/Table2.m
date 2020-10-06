%% Table 2 from Knies and Melo (2020)
% Calculates route choice probabilities for Figure 6
% Produces Table 2, Column 2
% Network originally featured in Mai et al. (2015) (Figure 3)
% Choice Aversion Model (Heterogeneous/Node-Specific Choice Aversion Parameter Extension)


% Route 1 : [o,a,a1,d]
% Route 2 : [o,a,a2,d]
% Route 3 : [o,a,a3,d]
% Route 4 : [o,b,b1,d]
% Route 5 : [o,b,b2,d]
% Route 6 : [o,b,b3,d]

% kappa_c = 2 and kappa_b = 1
clear all
clc
% Link/edge costs
to = 1;
ta = 1;
tb = 2;
ta1 = 1;
ta2 = 2;
ta3 = 3;
tb1 = 2;
tb2 = 1.5;
tb3 = 1;
td = 0;

% Route/path costs
c1 = to + ta + ta1 + td;
c2 = to + ta + ta2 + td;
c3 = to + ta + ta3 + td;
c4 = to + tb + tb1 + td;
c5 = to + tb + tb2 + td;
c6 = to + tb + tb3 + td;

% Kappa is the choice aversion parameter at each node other than the sink
% (end) node (directed acyclical graph assumption)
kapa = 1;
kapb = 1; 
kapc = 2; 
kapd = 1;

% Theta is the parameter on the cost function
theta = 1;


% Route utilities
u1 = -theta*c1;
u2 = -theta*c2;
u3 = -theta*c3;
u4 = -theta*c4;
u5 = -theta*c5;
u6 = -theta*c6;

% Choice set size for each route

Aj1=(2.^kapa).*(3.^kapb).*(1.^kapd);
Aj2=(2.^kapa).*(3.^kapb).*(1.^kapd);
Aj3=(2.^kapa).*(3.^kapb).*(1.^kapd);
Aj4=(2.^kapa).*(3.^kapc).*(1.^kapd);
Aj5=(2.^kapa).*(3.^kapc).*(1.^kapd);
Aj6=(2.^kapa).*(3.^kapc).*(1.^kapd);


% Route probabilities
P = zeros(6,1);
P(1) = exp(u1)./((Aj1./Aj1).*exp(u1)+(Aj1./Aj2).*exp(u2)+(Aj1./Aj3).*exp(u3)+(Aj1./Aj4).*exp(u4)+(Aj1./Aj5).*exp(u5)+(Aj1./Aj6).*exp(u6));
P(2) = exp(u2)./((Aj2./Aj1).*exp(u1)+(Aj2./Aj2).*exp(u2)+(Aj2./Aj3).*exp(u3)+(Aj2./Aj4).*exp(u4)+(Aj2./Aj5).*exp(u5)+(Aj2./Aj6).*exp(u6));
P(3) = exp(u3)./((Aj3./Aj1).*exp(u1)+(Aj3./Aj2).*exp(u2)+(Aj3./Aj3).*exp(u3)+(Aj3./Aj4).*exp(u4)+(Aj3./Aj5).*exp(u5)+(Aj3./Aj6).*exp(u6));
P(4) = exp(u4)./((Aj4./Aj1).*exp(u1)+(Aj4./Aj2).*exp(u2)+(Aj4./Aj3).*exp(u3)+(Aj4./Aj4).*exp(u4)+(Aj4./Aj5).*exp(u5)+(Aj4./Aj6).*exp(u6));
P(5) = exp(u5)./((Aj5./Aj1).*exp(u1)+(Aj5./Aj2).*exp(u2)+(Aj5./Aj3).*exp(u3)+(Aj5./Aj4).*exp(u4)+(Aj5./Aj5).*exp(u5)+(Aj5./Aj6).*exp(u6));
P(6) = exp(u6)./((Aj6./Aj1).*exp(u1)+(Aj6./Aj2).*exp(u2)+(Aj6./Aj3).*exp(u3)+(Aj6./Aj4).*exp(u4)+(Aj6./Aj5).*exp(u5)+(Aj6./Aj6).*exp(u6));

% saving first round of route choice probs
Prob = zeros(6,5);
Prob(:,1) = P;
save Table2.mat Prob

%% Table 2 from Knies and Melo (2020)
% Calculates route choice probabilities for Figure 6
% Produces Table 2, Column 3
% Network originally featured in Mai et al. (2015) (Figure 3)
% Choice Aversion Model (Heterogeneous/Node-Specific Choice Aversion Parameter Extension)

% removes a1
% Route 1 : [o,a,a1,d]% gone
% Route 2 : [o,a,a2,d]
% Route 3 : [o,a,a3,d]
% Route 4 : [o,b,b1,d]
% Route 5 : [o,b,b2,d]
% Route 6 : [o,b,b3,d]


clear all
clc
% Link/edge costs
to = 1;
ta = 1;
tb = 2;
%ta1 = 1;
ta2 = 2;
ta3 = 3;
tb1 = 2;
tb2 = 1.5;
tb3 = 1;
td = 0;

% Route/path costs
%c1 = to + ta + ta1 + td;
c2 = to + ta + ta2 + td;
c3 = to + ta + ta3 + td;
c4 = to + tb + tb1 + td;
c5 = to + tb + tb2 + td;
c6 = to + tb + tb3 + td;

% Kappa is the choice aversion parameter at each node other than the sink
% (end) node (directed acyclical graph assumption)
kapa = 1;
kapb = 1; 
kapc = 2; 
kapd = 1;

% Theta is the parameter on the cost function
theta = 1;


% Route utilities
%u1 = -theta*c1;
u2 = -theta*c2;
u3 = -theta*c3;
u4 = -theta*c4;
u5 = -theta*c5;
u6 = -theta*c6;

% Choice set size for each route

%Aj1=(2.^kapa).*(3.^kapb).*(1.^kapd);
Aj2=(2.^kapa).*(2.^kapb).*(1.^kapd);
Aj3=(2.^kapa).*(2.^kapb).*(1.^kapd);
Aj4=(2.^kapa).*(3.^kapc).*(1.^kapd);
Aj5=(2.^kapa).*(3.^kapc).*(1.^kapd);
Aj6=(2.^kapa).*(3.^kapc).*(1.^kapd);


% Route probabilities
P = zeros(6,1);
P(1) = NaN; %exp(u1)./((Aj1./Aj1).*exp(u1)+(Aj1./Aj2).*exp(u2)+(Aj1./Aj3).*exp(u3)+(Aj1./Aj4).*exp(u4)+(Aj1./Aj5).*exp(u5)+(Aj1./Aj6).*exp(u6));
P(2) = exp(u2)./((Aj2./Aj2).*exp(u2)+(Aj2./Aj3).*exp(u3)+(Aj2./Aj4).*exp(u4)+(Aj2./Aj5).*exp(u5)+(Aj2./Aj6).*exp(u6));
P(3) = exp(u3)./((Aj3./Aj2).*exp(u2)+(Aj3./Aj3).*exp(u3)+(Aj3./Aj4).*exp(u4)+(Aj3./Aj5).*exp(u5)+(Aj3./Aj6).*exp(u6));
P(4) = exp(u4)./((Aj4./Aj2).*exp(u2)+(Aj4./Aj3).*exp(u3)+(Aj4./Aj4).*exp(u4)+(Aj4./Aj5).*exp(u5)+(Aj4./Aj6).*exp(u6));
P(5) = exp(u5)./((Aj5./Aj2).*exp(u2)+(Aj5./Aj3).*exp(u3)+(Aj5./Aj4).*exp(u4)+(Aj5./Aj5).*exp(u5)+(Aj5./Aj6).*exp(u6));
P(6) = exp(u6)./((Aj6./Aj2).*exp(u2)+(Aj6./Aj3).*exp(u3)+(Aj6./Aj4).*exp(u4)+(Aj6./Aj5).*exp(u5)+(Aj6./Aj6).*exp(u6));


load Table2.mat
% saving second round of route choice probs

Prob(:,2) = P;
save Table2.mat Prob

%% Table 2 from Knies and Melo (2020)
% Calculates route choice probabilities for Figure 6
% Produces Table 2, Column 4
% Network originally featured in Mai et al. (2015) (Figure 3)
% Choice Aversion Model (Heterogeneous/Node-Specific Choice Aversion Parameter Extension)

% removes a2
% Route 1 : [o,a,a1,d]
% Route 2 : [o,a,a2,d] % gone
% Route 3 : [o,a,a3,d]
% Route 4 : [o,b,b1,d]
% Route 5 : [o,b,b2,d]
% Route 6 : [o,b,b3,d]


clear all
clc
% Link/edge costs
to = 1;
ta = 1;
tb = 2;
ta1 = 1;
%ta2 = 2;
ta3 = 3;
tb1 = 2;
tb2 = 1.5;
tb3 = 1;
td = 0;

% Route/path costs
c1 = to + ta + ta1 + td;
%c2 = to + ta + ta2 + td;
c3 = to + ta + ta3 + td;
c4 = to + tb + tb1 + td;
c5 = to + tb + tb2 + td;
c6 = to + tb + tb3 + td;

% Kappa is the choice aversion parameter at each node other than the sink
% (end) node (directed acyclical graph assumption)
kapa = 1;
kapb = 1; 
kapc = 2; 
kapd = 1;

% Theta is the parameter on the cost function
theta = 1;


% Route utilities
u1 = -theta*c1;
%u2 = -theta*c2;
u3 = -theta*c3;
u4 = -theta*c4;
u5 = -theta*c5;
u6 = -theta*c6;

% Choice set size for each route

Aj1=(2.^kapa).*(2.^kapb).*(1.^kapd);
%Aj2=(2.^kapa).*(3.^kapb).*(1.^kapd);
Aj3=(2.^kapa).*(2.^kapb).*(1.^kapd);
Aj4=(2.^kapa).*(3.^kapc).*(1.^kapd);
Aj5=(2.^kapa).*(3.^kapc).*(1.^kapd);
Aj6=(2.^kapa).*(3.^kapc).*(1.^kapd);


% Route probabilities
P = zeros(6,1);
P(1) = exp(u1)./((Aj1./Aj1).*exp(u1)+(Aj1./Aj3).*exp(u3)+(Aj1./Aj4).*exp(u4)+(Aj1./Aj5).*exp(u5)+(Aj1./Aj6).*exp(u6));
P(2) = NaN; %exp(u2)./((Aj2./Aj1).*exp(u1)+(Aj2./Aj2).*exp(u2)+(Aj2./Aj3).*exp(u3)+(Aj2./Aj4).*exp(u4)+(Aj2./Aj5).*exp(u5)+(Aj2./Aj6).*exp(u6));
P(3) = exp(u3)./((Aj3./Aj1).*exp(u1)+(Aj3./Aj3).*exp(u3)+(Aj3./Aj4).*exp(u4)+(Aj3./Aj5).*exp(u5)+(Aj3./Aj6).*exp(u6));
P(4) = exp(u4)./((Aj4./Aj1).*exp(u1)+(Aj4./Aj3).*exp(u3)+(Aj4./Aj4).*exp(u4)+(Aj4./Aj5).*exp(u5)+(Aj4./Aj6).*exp(u6));
P(5) = exp(u5)./((Aj5./Aj1).*exp(u1)+(Aj5./Aj3).*exp(u3)+(Aj5./Aj4).*exp(u4)+(Aj5./Aj5).*exp(u5)+(Aj5./Aj6).*exp(u6));
P(6) = exp(u6)./((Aj6./Aj1).*exp(u1)+(Aj6./Aj3).*exp(u3)+(Aj6./Aj4).*exp(u4)+(Aj6./Aj5).*exp(u5)+(Aj6./Aj6).*exp(u6));

load Table2.mat
% saving third round of route choice probs

Prob(:,3) = P;
save Table2.mat Prob

%% Table 2 from Knies and Melo (2020)
% Calculates route choice probabilities for Figure 6
% Produces Table 2, Column 5
% Network originally featured in Mai et al. (2015) (Figure 3)
% Choice Aversion Model (Heterogeneous/Node-Specific Choice Aversion Parameter Extension)

% removes b1
% Route 1 : [o,a,a1,d]
% Route 2 : [o,a,a2,d]
% Route 3 : [o,a,a3,d]
% Route 4 : [o,b,b1,d]% gone
% Route 5 : [o,b,b2,d]
% Route 6 : [o,b,b3,d]


clear all
clc
% Link/edge costs
to = 1;
ta = 1;
tb = 2;
ta1 = 1;
ta2 = 2;
ta3 = 3;
%tb1 = 2;
tb2 = 1.5;
tb3 = 1;
td = 0;

% Route/path costs
c1 = to + ta + ta1 + td;
c2 = to + ta + ta2 + td;
c3 = to + ta + ta3 + td;
%c4 = to + tb + tb1 + td;
c5 = to + tb + tb2 + td;
c6 = to + tb + tb3 + td;

% Kappa is the choice aversion parameter at each node other than the sink
% (end) node (directed acyclical graph assumption)
kapa = 1;
kapb = 1; 
kapc = 2; 
kapd = 1;

% Theta is the parameter on the cost function
theta = 1;


% Route utilities
u1 = -theta*c1;
u2 = -theta*c2;
u3 = -theta*c3;
%u4 = -theta*c4;
u5 = -theta*c5;
u6 = -theta*c6;

% Choice set size for each route

Aj1=(2.^kapa).*(3.^kapb).*(1.^kapd);
Aj2=(2.^kapa).*(3.^kapb).*(1.^kapd);
Aj3=(2.^kapa).*(3.^kapb).*(1.^kapd);
%Aj4=(2.^kapa).*(3.^kapc).*(1.^kapd);
Aj5=(2.^kapa).*(2.^kapc).*(1.^kapd);
Aj6=(2.^kapa).*(2.^kapc).*(1.^kapd);


% Route probabilities
P = zeros(6,1);
P(1) = exp(u1)./((Aj1./Aj1).*exp(u1)+(Aj1./Aj2).*exp(u2)+(Aj1./Aj3).*exp(u3)+(Aj1./Aj5).*exp(u5)+(Aj1./Aj6).*exp(u6));
P(2) = exp(u2)./((Aj2./Aj1).*exp(u1)+(Aj2./Aj2).*exp(u2)+(Aj2./Aj3).*exp(u3)+(Aj2./Aj5).*exp(u5)+(Aj2./Aj6).*exp(u6));
P(3) = exp(u3)./((Aj3./Aj1).*exp(u1)+(Aj3./Aj2).*exp(u2)+(Aj3./Aj3).*exp(u3)+(Aj3./Aj5).*exp(u5)+(Aj3./Aj6).*exp(u6));
P(4) = NaN; %exp(u4)./((Aj4./Aj1).*exp(u1)+(Aj4./Aj2).*exp(u2)+(Aj4./Aj3).*exp(u3)+(Aj4./Aj4).*exp(u4)+(Aj4./Aj5).*exp(u5)+(Aj4./Aj6).*exp(u6));
P(5) = exp(u5)./((Aj5./Aj1).*exp(u1)+(Aj5./Aj2).*exp(u2)+(Aj5./Aj3).*exp(u3)+(Aj5./Aj5).*exp(u5)+(Aj5./Aj6).*exp(u6));
P(6) = exp(u6)./((Aj6./Aj1).*exp(u1)+(Aj6./Aj2).*exp(u2)+(Aj6./Aj3).*exp(u3)+(Aj6./Aj5).*exp(u5)+(Aj6./Aj6).*exp(u6));

load Table2.mat
% saving fourth round of route choice probs

Prob(:,4) = P;
save Table2.mat Prob

%% Table 2 from Knies and Melo (2020)
% Calculates route choice probabilities for Figure 6
% Produces Table 2, Column 6
% Network originally featured in Mai et al. (2015) (Figure 3)
% Choice Aversion Model (Heterogeneous/Node-Specific Choice Aversion Parameter Extension)

% removes b2
% Route 1 : [o,a,a1,d]
% Route 2 : [o,a,a2,d]
% Route 3 : [o,a,a3,d]
% Route 4 : [o,b,b1,d]
% Route 5 : [o,b,b2,d]% gone
% Route 6 : [o,b,b3,d]


clear all
clc
% Link/edge costs
to = 1;
ta = 1;
tb = 2;
ta1 = 1;
ta2 = 2;
ta3 = 3;
tb1 = 2;
%tb2 = 1.5;
tb3 = 1;
td = 0;

% Route/path costs
c1 = to + ta + ta1 + td;
c2 = to + ta + ta2 + td;
c3 = to + ta + ta3 + td;
c4 = to + tb + tb1 + td;
%c5 = to + tb + tb2 + td;
c6 = to + tb + tb3 + td;

% Kappa is the choice aversion parameter at each node other than the sink
% (end) node (directed acyclical graph assumption)
kapa = 1;
kapb = 1; 
kapc = 2; 
kapd = 1;

% Theta is the parameter on the cost function
theta = 1;


% Route utilities
u1 = -theta*c1;
u2 = -theta*c2;
u3 = -theta*c3;
u4 = -theta*c4;
%u5 = -theta*c5;
u6 = -theta*c6;

% Choice set size for each route

Aj1=(2.^kapa).*(3.^kapb).*(1.^kapd);
Aj2=(2.^kapa).*(3.^kapb).*(1.^kapd);
Aj3=(2.^kapa).*(3.^kapb).*(1.^kapd);
Aj4=(2.^kapa).*(2.^kapc).*(1.^kapd);
%Aj5=(2.^kapa).*(3.^kapc).*(1.^kapd);
Aj6=(2.^kapa).*(2.^kapc).*(1.^kapd);


% Route probabilities
P = zeros(6,1);
P(1) = exp(u1)./((Aj1./Aj1).*exp(u1)+(Aj1./Aj2).*exp(u2)+(Aj1./Aj3).*exp(u3)+(Aj1./Aj4).*exp(u4)+(Aj1./Aj6).*exp(u6));
P(2) = exp(u2)./((Aj2./Aj1).*exp(u1)+(Aj2./Aj2).*exp(u2)+(Aj2./Aj3).*exp(u3)+(Aj2./Aj4).*exp(u4)+(Aj2./Aj6).*exp(u6));
P(3) = exp(u3)./((Aj3./Aj1).*exp(u1)+(Aj3./Aj2).*exp(u2)+(Aj3./Aj3).*exp(u3)+(Aj3./Aj4).*exp(u4)+(Aj3./Aj6).*exp(u6));
P(4) = exp(u4)./((Aj4./Aj1).*exp(u1)+(Aj4./Aj2).*exp(u2)+(Aj4./Aj3).*exp(u3)+(Aj4./Aj4).*exp(u4)+(Aj4./Aj6).*exp(u6));
P(5) = NaN;%exp(u5)./((Aj5./Aj1).*exp(u1)+(Aj5./Aj2).*exp(u2)+(Aj5./Aj3).*exp(u3)+(Aj5./Aj4).*exp(u4)+(Aj5./Aj5).*exp(u5)+(Aj5./Aj6).*exp(u6));
P(6) = exp(u6)./((Aj6./Aj1).*exp(u1)+(Aj6./Aj2).*exp(u2)+(Aj6./Aj3).*exp(u3)+(Aj6./Aj4).*exp(u4)+(Aj6./Aj6).*exp(u6));

load Table2.mat
% saving fifth round of route choice probs

Prob(:,5) = P;
save Table2.mat Prob

% produces matrix of percent change from the baseline case
PercentChange = (Prob - Prob(:,1))./(Prob(:,1));
