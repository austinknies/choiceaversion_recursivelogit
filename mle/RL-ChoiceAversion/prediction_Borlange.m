%% Compute Prediction

function [LogL] = prediction_Borlange(isLS, isObs);

    globalVar;
    global TXT; 
    global SampleObs;
    
    saveResults = true;

    %% Data 

    file_AttEstimatedtime='./Input/Borlange/ATTRIBUTEestimatedtime.txt';
    file_incidence='./Input/Borlange/linkIncidence.txt';
    file_observations='./Input/Borlange/observationsForEstimBAI.txt';

    isLinkSizeInclusive = isLS;
    isFixedUturn = false;
    loadData_Borlange;
    
    if isObs ==  true
        PredSample = spconvert(load('./Input/Borlange/PredSampleObs.txt'));
        note = 'OBS';
    else
        PredSample = spconvert(load('./Input/Borlange/PredSampleODs.txt'));
        note = 'ODS';
    end
    TXT = ['RL prediction:',note,':'];
    TXT = [TXT sprintf('\n Link size = %d \n', isLinkSizeInclusive)];
    nTest = round(size(PredSample,1) / 2);
    %% Estimation for training samples
    for ii = 1: nTest
        
        Op = Op_structure;
        initialize_optimization_structure_Borlange();
        
        train = PredSample(ii*2-1,:);
        test  = PredSample(ii*2,:);
        SampleObs = train; 
        
        %% RLoptimizer part - Routine for training sample
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
        
        
        %% Now getting LL for holdout sample
        [LogL(ii),~] = getPLL(test);
        disp('Check')
        disp(ii)
        disp(LogL(ii))
        resultsTXT = [resultsTXT sprintf('\n Predicted LL %d \n', LogL(ii))];
    end    
end