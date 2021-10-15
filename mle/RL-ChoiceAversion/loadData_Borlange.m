%   Load route choice data

%% Read files
disp('Loading data ...')

global incidenceFull;
global Atts;
global Obs;
global nbobs;
global LSatt;
global isLinkSizeInclusive;

file_AttEstimatedtime='./Input/Borlange/ATTRIBUTEestimatedtime.txt';
file_incidence='./Input/Borlange/linkIncidence.txt';
file_observations='./Input/Borlange/observationsForEstimBAI.txt';

% get attributes and incidence data from files
incidenceFull = spconvert(load(file_incidence));

% number of network links = nbNetworkStates
% number of network + dummy links = nbTotalStates
[nbNetworkStates,nbTotalStates]=size(incidenceFull);

EstimatedTime = spconvert(load(file_AttEstimatedtime));
EstimatedTime(nbNetworkStates,nbTotalStates)=0;

Obs = spconvert(load(file_observations));
nbobs = length(Obs(:,1));

 

%% Define link attributes 

% choose appropriate columns from attributes data matrix

ChoiceSize=log(sum(incidenceFull,2)); % choice aversion term -> log(choice set cardinality)

%% Define link pairs attributes

% extend attribute value to absorbing links as 0
% attribute value must always be 0 for absorbing links!
ChoiceSize(nbTotalStates) = 0;
icd=incidenceFull;
icd(:,nbTotalStates) = 0;

% create link pair matrix x(a|k) for each attribute
% First index is k, second is a, attribute is for link a
icdChoiceSize = incidenceFull * spdiags(ChoiceSize,0,nbTotalStates,nbTotalStates);
icdEstimatedTime = (incidenceFull .* EstimatedTime);
% put everything in Atts variables
Atts  = objArray(1);
Atts(1).value = icd; %link dummy
Atts(2).value = icdEstimatedTime; % travel time
Atts(3).value = icdChoiceSize; % choice aversion term

%% Load link size attribute if needed

if isLinkSizeInclusive == true
    %either compute the link size attribute
    ExpV_is_ok = getLinkSizeAtt();
    if ExpV_is_ok ==0
        disp('Warning: failed to compute link size attribute')
    end
    %or load an already computed and saved link size attribute
%     load('./Input/Example1/LSatt.mat')
end