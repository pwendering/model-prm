function projSol = performProjectionCPLEX(model, v)
% This function updates the right-hand side of constraints that specify the
% minimization of absolute distances between the predicted and a random
% flux vector. These constraints have the prefix 'norm_NF_' for each net flux
% reaction in the TFA model.
% INPUT
%   struct model:           TFA model created with the matTFA toolbox;
%                           contains constraints for the minimization of
%                           the first norm of fluxes through net flux
%                           reactions.
%   double v:               vector of random flux values within ranges
%                           determined by flux variability analysis with
%                           thermodynamic constraints.
% OUTPUT
%   struct projSol:         solution structure

% find minimization constraints
norm1_const_idx = startsWith(model.constraintNames,'norm1_NF_');
% add RHS values
model.rhs(norm1_const_idx) = v;
% solve problem
projSol = solveTFAmodelCplex(model,30);
end