function projSol = performProjectionCPLEX(tModel,v)
norm1_const_idx = startsWith(tModel.constraintNames,'norm1_NF_');
tModel.rhs(norm1_const_idx) = v;
projSol = solveTFAmodelCplex(tModel,30);
end