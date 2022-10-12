function pfba_model = addPfbaConstTFA(model)
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
nf_idx = startswith(model.varNames, 'NF_');
n_nf = sum(nf_idx);

% construct additional matrix with delta variables to be appended to the
% existing constraint matrix
delta_mat = [zeros(n_nf, size(model.A,2)) -eye(n_nf) eye(n_nf)];
% fill in ones for NF reaction indices
delta_mat(1:n_nf, nf_idx) = 1;

% append new matrix to existing matrix of TFA problem
pfba_model = model;
pfba_model.A = [[pfba_model.A zeros(size(pfba_model.A, 1), 2*n_nf)]; delta_mat];

% add additional information for constraints and variables
pfba_model.rhs = [pfba_model.rhs; zeros(n_nf, 1)];
pfba_model.constraintNames = [pfba_model.constraintNames; strcat(model.varNames(nf_idx), '_pfba')];
pfba_model.constraintType = [pfba_model.constraintType; repelem('=', n_nf, 1)];

pfba_model.varNames = [pfba_model.varNames; strcat(pfba_model.varNames(nf_idx), '_delta_plus');...
    strcat(pfba_model.varNames(nf_idx), '_delta_minus')];
pfba_model.var_lb = [pfba_model.var_lb; zeros(2*n_nf, 1)];
pfba_model.var_ub = [pfba_model.var_ub; repelem(1000, 2*n_nf, 1)];

% set objective to minimization of deltas
pfba_model.f = [zeros(size(model.f)) ones(2*n_nf, 1)];

% set objective sense to minimization
pfba_model.objtype = 1;

end