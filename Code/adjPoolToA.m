function error = adjPoolToA(model,PoolUB,A_ref)

model.ub(strcmp(model.rxns,'prot_pool_exchange')) = PoolUB;
params.feasTol = 1e-9;

s = optimizeCbModel(model,'max','one',1,params);

A = [1 -0.5 -1] * s.x(findRxnIDs(model,{'arm_RBC_h' 'arm_RBO_h' 'Tr_CO2m_REV'}));
error = abs(A-A_ref);

end