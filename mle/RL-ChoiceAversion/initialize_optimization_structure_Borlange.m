%   Initialize optimization structure
%%
function [] = initialize_optimization_structure_Borlange()
    global Op;
    global isLinkSizeInclusive;
    Op.Optim_Method = OptimizeConstant.TRUST_REGION_METHOD;
    Op.Hessian_approx = OptimizeConstant.BFGS; 
    Op.ETA1 = 0.05;
    Op.ETA2 = 0.75;
    Op.maxIter = 150;
    Op.k = 0;
    Op.n = 3; % includes choice aversion
    if isLinkSizeInclusive == true
%         Op.n = Op.n + 1; % just + 1 if no interaction term
        Op.n = Op.n + 2; % for interaction term
    end
    Op.x = -ones(Op.n,1) * 1.5;
    Op.tol = 1e-6;
    Op.radius = 1.0;
    Op.Ak = eye(Op.n);
    Op.H = eye(Op.n);
end