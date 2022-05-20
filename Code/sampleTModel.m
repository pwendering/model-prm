function [fluxSamples,concSamples,norm1Obj] = sampleTModel(tModel,minFlux,maxFlux,n,threads)
%% [fluxSamples,concSamples] = sampleTModel(tModel,minFlux,maxFlux,n,threads)
% Projection sampling for net fluxes in a thermodynamically constraint
% metabolic model [matTFA, Salvy et al. (2019)].
% The sampling is done by creating a random vector of fluxes between within
% the feasible ranges as determined by variability analysis, followed by
% minimizing the first norm of differences to the predicted flux vector.
%
% INPUT:
%   struct tModel:              thermodynamically constraint metabolic
%                               model (matTFA toolbox)
%   double minFlux,maxFlux:     feasible ranges of reaction fluxes
%                               determined by TFA
%   double n:                   number of samples
%   double threads:             number of threads to use
%
% OUTPUT:
%   double fluxSamples:         matrix of sampled reaction fluxes (|R| x n)
%   double concSamples:         matrix of sampled concentrations, which are
%                               associated with each of the sampled points
%   double norm1Obj:            objective values from norm 1 minimzation to
%                               random flux vector

rng('default')

if any(minFlux > maxFlux)
    error('At least one minimum flux value is greater than the maximum')
end

if threads > 1
    [status,ME] = setupParpool(threads);
    if status
        warning('Error: %s\nParallel pool could not be initialized, proceeding without parallelization',ME.message)
        threads = 1;
    end
end
nfIdx = startsWith(tModel.varNames,'NF_');
lcIdx = startsWith(tModel.varNames,'LC_');

% initialize output
fluxSamples = sparse(sum(nfIdx),n);
concSamples = sparse(sum(lcIdx),n);
norm1Obj = sparse(n,1);

% add first norm minimization to model
tModel_norm1 = addMinNorm1Const(tModel);
clear tModel

if threads == 1
    
    printMsg = '-- Done with 0 samples\n';
    
    for i=1:n
        
        v = createRandomFluxVector(minFlux,maxFlux);
        projSol = performProjectionCPLEX(tModel_norm1,v);
        
        if ~isempty(projSol.x)
            [fluxCol,concCol,objVal] = deal(projSol.x(nfIdx),projSol.x(lcIdx),projSol.val);
            fluxSamples(:,i) = fluxCol;
            concSamples(:,i) = concCol;
            norm1Obj(i) = objVal;
        else
            fprintf(printMsg)
        end
        
        if mod(i,100)==0
            if i>100
                fprintf(repmat('\b',1,length(printMsg)))
            end
            printMsg = regexprep(printMsg,'\d+',num2str(i));
            fprintf(printMsg)
        end
    end
    
else
    
    environment = getEnvironment;
    solver = getCobraSolver('LP');
    
    parfor i=1:n
        
        restoreEnvironment(environment);
        changeCobraSolver(solver, 'LP', 0, -1);
        
        v = createRandomFluxVector(minFlux,maxFlux);
        projSol = performProjectionCPLEX(tModel_norm1,v);
        
        if ~isempty(projSol.x)
            [fluxCol,concCol,objVal] = deal(projSol.x(nfIdx),projSol.x(lcIdx),projSol.val);
            fluxSamples(:,i) = fluxCol;
            concSamples(:,i) = exp(concCol);
            norm1Obj(i) = objVal;
        end
    end
end

    function tModel_norm1 = addMinNorm1Const(tModel)
        nfIdx = startsWith(tModel.varNames,'NF_');
        [nRow,nVar] = size(tModel.A);
        nNF = sum(nfIdx);
        nDelta = 2*nNF;
        % create constraint matrix for norm 1
        norm1Mat = sparse(nNF,nVar+nDelta);
        norm1Mat(:,nfIdx) = eye(nNF); % v(net flux)
        norm1Mat(:,nVar+1:nVar+nNF) = -eye(nNF); %delta+
        norm1Mat(:,nVar+nNF+1:nVar+nDelta) = eye(nNF); %delta-
        % add constraints to model
        tModel_norm1 = struct;
        tModel_norm1.A = [tModel.A zeros(nRow,nDelta); norm1Mat];
        tModel_norm1.rhs = sparse([tModel.rhs; zeros(nNF,1)]);
        tModel_norm1.var_lb = sparse([tModel.var_lb; zeros(nDelta,1)]);
        tModel_norm1.var_ub = [tModel.var_ub; 1000*ones(nDelta,1)];
        tModel_norm1.vartypes = [tModel.vartypes; repmat({'C'},nDelta,1)];
        tModel_norm1.constraintType = [tModel.constraintType; repmat({'='},nNF,1)];
        tModel_norm1.constraintNames = [tModel.constraintNames; strcat('norm1_',tModel.varNames(nfIdx))];
        tModel_norm1.f = [zeros(nVar,1); ones(nDelta,1)];
        tModel_norm1.objtype = 1;
        tModel_norm1.description = tModel.description;
    end


end