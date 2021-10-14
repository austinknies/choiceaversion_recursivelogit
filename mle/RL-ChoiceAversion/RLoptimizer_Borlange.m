....
%   Link-based network route choice model with unrestricted choice set
%   Optimization algorithm
%   Chaire CN - DIRO - Université de Montréal

%   Modified by Austin Knies - Indiana University Bloomington (August 2021)
%   Incorporates choice aversion term as in Knies, Lorca, and Melo (2021)
%   
%   Runs Borlange data used in Fosgerau et al. (2013) 


%   MAIN PROGRAM
%   ---------------------------------------------------
%%
Credits;

% declare global variables
globalVar; 

% set the link size attribute
isLinkSizeInclusive = false;
% When running Link Size, to determine whether you want to explore the
% interaction between the choice aversion penalization and link size, you
% also need to modify:
%           - Lines 14&15 and 21&22 in getLinkSizeAtt.m
%           - Lines 14&15 in initialize_optimization_structure_Borlange.m
%           - Lines 158&159 in 172-183 in getLL.m
saveResults = true;

% set the optimization parameters
Op = Op_structure;
initialize_optimization_structure_Borlange();

% load network attributes and observations
loadData_Borlange;



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
