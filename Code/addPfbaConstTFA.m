function pfba_model = addPfbaConstTFA(model, rhsVector)
%% pfba_model = addPfbaConstTFA(model)
% This function adds constraints to a TFA problem to enable minimization of
% absolute net fluxes. This is done by decomposing reaction net flux into a
% postitive (delta+) and negative(delta-) component. The objective is then
% set to minimization of the sum of delta+ and delta-.
% 
% INTPUT:
%   struct model:           TFA model (matTFA)
% OUTPUT:
%   struct pfba_model:      TFA model with pFBA constraints

% find net flux variables
nf_idx = startsWith(model.varNames, 'NF_');
% exlude biomass reaction
nf_idx(ismember(model.varNames, strcat('N',model.varNames(model.f==1)))) = 0;
n_nf = sum(nf_idx);

if nargin < 2 || isempty(rhsVector)
    rhsVector = zeros(n_nf, 1);
elseif numel(rhsVector) ~= n_nf
    error('Size of provided RHS vector does not match the number of NF variables')
end

% construct additional matrix with delta variables to be appended to the
% existing constraint matrix
delta_mat = [zeros(n_nf, size(model.A,2)) -eye(n_nf) eye(n_nf)];
% fill in ones for NF reaction indices
delta_mat(:, nf_idx) = eye(n_nf);

% append new matrix to existing matrix of TFA problem
pfba_model = model;
pfba_model.A = [[pfba_model.A zeros(size(pfba_model.A, 1), 2*n_nf)]; delta_mat];

% add additional information for constraints and variables
pfba_model.rhs = [pfba_model.rhs; rhsVector];
pfba_model.constraintNames = [pfba_model.constraintNames; strcat('pfba_', model.varNames(nf_idx))];
pfba_model.constraintType = [pfba_model.constraintType; repelem({'='}, n_nf, 1)];

pfba_model.varNames = [pfba_model.varNames; strcat('delta_plus_', pfba_model.varNames(nf_idx));...
    strcat('delta_minus_', pfba_model.varNames(nf_idx))];
pfba_model.var_lb = [pfba_model.var_lb; zeros(2*n_nf, 1)];
pfba_model.var_ub = [pfba_model.var_ub; repelem(1000, 2*n_nf, 1)];
pfba_model.vartypes = [pfba_model.vartypes; repmat({'C'}, 2*n_nf, 1)];
% set objective to minimization of deltas
pfba_model.f = [zeros(size(model.f)); ones(2*n_nf, 1)];

% set objective sense to minimization
pfba_model.objtype = 1;

end