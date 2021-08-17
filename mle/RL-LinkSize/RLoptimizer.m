....
%   Link-based network route choice model with unrestricted choice set
%   Optimization algorithm
%   Chaire CN - DIRO - Université de Montréal

%   Modified by Austin Knies - Indiana University Bloomington (August 2021)

%   Changes to original code:
%   -   Link 7 only has two links attached in the network depicted in 
%       tutorial slides: 10 and 17. We break from the diagram and allow 19
%       to connect from 7 to 29, since it would otherwise be missing.
%       Our recreation of this network (Figure 9) reflects this.
%   -   In the incidence matrix, link 20 appears to connect to itself. 
%       We have removed this, so that link 20 only "connects" to 29.
%   -   In loadData.m, the following code erased the connection between
%       link 20 and 29:
%                   icd(:,nbNetworkStates:nbTotalStates) = 0;
%       I have changed it to
%                   icd(:,nbTotalStates) = 0;


%   MAIN PROGRAM
%   ---------------------------------------------------
%%
Credits;

% declare global variables
globalVar; 

% set the link size attribute
isLinkSizeInclusive = true;
saveResults = true;

% load network attributes and observations
loadData;

% set the optimization parameters
Op = Op_structure;
initialize_optimization_structure();

Gradient = zeros(nbobs,Op.n);

%---------------------------
%Starting optimization
tic ;
disp('Start Optimizing ....')
[Op.value, Op.grad ] = LL(Op.x);
PrintOut(Op);
% print result to string text
header = [sprintf('%s \n',file_observations) Op.Optim_Method];
header = [header sprintf('\nNumber of observations = %d \n', nbobs)];
header = [header sprintf('Hessian approx methods = %s \n', OptimizeConstant.getHessianApprox(Op.Hessian_approx))];
resultsTXT = header;
%------------------------------------------------
while (true)    
  Op.k = Op.k + 1;
  if strcmp(Op.Optim_Method,OptimizeConstant.LINE_SEARCH_METHOD);
    ok = line_search_iterate();
    if ok == true
        PrintOut(Op);
    else
        disp(' Unsuccessful process ...')
        break;
    end
  else
    ok = btr_interate();
    PrintOut(Op);
  end
  [isStop, Stoppingtype, isSuccess] = CheckStopping(Op);  
  %----------------------------------------
  % Check stopping criteria
  if isStop == true
      isSuccess
      fprintf('The algorithm stops, due to %s', Stoppingtype);
      resultsTXT = [resultsTXT sprintf('The algorithm stops, due to %s \n', Stoppingtype)];
      break;
  end
end

% Compute variance-covariance matrix
PrintOut(Op);
disp(' Calculating VAR-COV ...');
global Stdev;
Stdev = zeros(1,Op.n);
getCov;

%   Finishing ...
ElapsedTime = toc
resultsTXT = [resultsTXT sprintf('\n Number of function evaluation %d \n', Op.nFev)];
resultsTXT = [resultsTXT sprintf('\n Estimated time %d \n', ElapsedTime)];

if saveResults == true
    SaveResults(Stoppingtype, ElapsedTime);
end
