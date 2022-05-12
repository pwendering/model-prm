%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Modeling of photorespiratory mutants using metabolite abundance ratios
% Philipp Wendering, University of Potsdam (philipp.wendering@gmail.com)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize COBRA
initCobraToolbox(false)

% set up parallel pool
% delete(gcp('nocreate'));
% parpool(20);

% define output directory
res_dir = '../Results';

%% Read data

data_dir = '../Data';
data_dir_struct = dir(data_dir);
data_files = {data_dir_struct.name};

% metabolite abundance ranges from multiple experiments [g/gDW]
met_ranges = readtable(...
    fullfile(data_dir, 'met_ranges.csv'),...
    'ReadRowNames', true...
    );
met_range_av = mean(met_ranges{:,2:end},2,'omitnan');
met_range_sd = std(met_ranges{:,2:end},0,2,'omitnan');
met_range_sd(met_range_sd==0) = 0.3*met_range_av(met_range_sd==0);
range_met_names = met_ranges.Properties.RowNames;
range_ids = name2id(range_met_names);

min_conc_limit = min(min(met_ranges{:,2:end}));
max_conc_limit = max(max(met_ranges{:,2:end}));

% metabolic model (GECKO-formatted model of AraCore model)
% model reference: Arnold and Nikoloski (2014), doi: 10.1104/pp.114.235358
model = readCbModel(fullfile(data_dir,'raw_batch_ecModel.mat'));
RXN_IDX = find(~startsWith(model.rxns,{'draw_','prot_'}));
DRAW_RXN_IDX = startsWith(model.rxns,'draw_');
ENZ_MET_IDX = find(...
    startsWith(model.mets,'prot_') &...
    ~contains(model.mets, 'prot_pool')...
    );
fprintf('Setting sigma to 1 for GECKO model\n')
model.ub(findRxnIDs(model,'prot_pool_exchange')) = ...
    2*model.ub(findRxnIDs(model,'prot_pool_exchange'));

% set objective to light-limiting biomass reaction
model.c(:) = 0;
model.c(findRxnIDs(model,'Bio_opt')) = 1;
model.ub(findRxnIDs(model,{'Bio_CLim','Bio_NLim'})) = 0;

% dry weights
dw = readtable(...
    fullfile(data_dir, data_files{contains(data_files, 'dw_')}),...
    'ReadRowNames', true...
    );

% net CO2 assimilation rates
A_ml = readtable(...
    fullfile(data_dir, data_files{contains(data_files, 'A_ml')})...
    );

A_fl = readtable(...
    fullfile(data_dir, data_files{contains(data_files, 'A_fl')})...
    );

A_max = max([A_fl.mean_A; A_ml.mean_A]);

% photon uptake
LMA = 22.2; % [gDW/m2] (Hummel et al. 2010, doi: 10.1104/pp.110.157008)
I_ML = 3600 * 200 / LMA / 1000; % [mmol/gDW/h]
I_FL = 3600 * 700 / LMA / 1000; % [mmol/gDW/h] (use kinetic model to integrate fluctuations?)

% growth duration
t_growth = 34 * 24; % [h]

% calculate growth rate
dw_0 = 20 / 1e6;
mu = (log(dw.mean_DW) - log(dw_0)) ./ t_growth;
mu_wt = mu(strcmp(dw.Properties.RowNames,'Col0'));

mu(mu>mu_wt) = mu_wt;

% define mutant names and reaction knock-outs
mutants = {...
    'ggt1-1'...
    'ggt1-2'...
    'ggt2'...
    'hpr1-1'...
    'hpr1-2'...
    };

ko_rxns = {...
    findRxnsFromMets(model,'prot_Q9LR30')...
    findRxnsFromMets(model,'prot_Q9LR30')...
    findRxnsFromMets(model,'prot_Q9S7E9')...
    findRxnsFromMets(model,'prot_Q9C9W5')...
    findRxnsFromMets(model,'prot_Q9C9W5')...
    };

% light conditions
l_conds = {'ml','fl','ml_fl','fl_ml'};

% time points
tp = [-21 3 27 75 123];

%% Step 1 - minimize total flux with GECKO model
fprintf('Setting photon uptake to 200 umol/m2/s\n')
model.ub(findRxnIDs(model,'Im_hnu')) = I_ML;

% solve FBA problem and minimize sum of fluxes using first norm
params.feasTol = 1e-9;
step_1_solution = optimizeCbModel(model,'max','one',1,params);
% perform FVA for wild type at 90% of optimal growth
[minFluxWT,maxFluxWT] = fluxVariability(model,...
    'optPercentage',90,...
    'solverParams',params...
    );

% calculate C^wt = v / (kcat * E) for each reaction (per enzyme)
fprintf('Calculating C_wt...\n')
C_wt = calculateVBykE(model,step_1_solution.x);

l_idx = 1;

for t_idx = 1:numel(tp)
    cmdsz = matlab.desktop.commandwindow.size;
    fprintf([repmat('-',1,cmdsz(1)) '\n'])
    fprintf('Time after shift: %d\n', tp(t_idx))
    for m_idx = 1:numel(mutants)
        
        cmdsz = matlab.desktop.commandwindow.size;
        fprintf([repmat('-',1,floor(cmdsz(1)/2)) '\n'])
        fprintf('>>> Processing genotype %s ...\n', mutants{m_idx})
        
        % test if growth is possible with knock-outs
        ko_model = removeRxns(model,ko_rxns{m_idx});
        g_ko = optimizeCbModel(ko_model).f;
        g_wt = optimizeCbModel(model).f;
        gt_bool = g_ko > 0;
        if gt_bool
            fprintf('Growth test passed\n')
            fprintf('FBA results: WT: %.4e, %s: %.4e\n', g_ko, mutants{m_idx}, g_wt)
        else
            fprintf('KO disables growth prediction - skipping\n')
            continue
        end
        
        % determine ratio of growth rates (wt/mutant)
        biomass_ratio = mu_wt / ...
            mu(strcmp(dw.Properties.RowNames,mutants{m_idx}));
        
        % calculate B (derived from Michaelis-Menten rate law for multiple substrates)
        met_av = readtable(fullfile(data_dir, ['met_' l_conds{l_idx} '.csv']));
        met_sd = readtable(fullfile(data_dir, ['met_' l_conds{l_idx} '_sd.csv']));
        
        meas_names = met_av.Properties.VariableNames(3:end);
        meas_ids = name2id(meas_names);
        
        exp_data = struct;
        exp_data.meas_ids = meas_ids;
        exp_data.met_av = met_av;
        exp_data.met_sd = met_sd;
        
        range_data = struct;
        range_data.range_ids = range_ids;
        range_data.met_range_av = met_range_av;
        range_data.met_range_sd = met_range_sd;
        range_data.min_conc_limit = min_conc_limit;
        range_data.max_conc_limit = max_conc_limit;
        
        [B_wt_min,B_wt_mean,B_wt_max,B_mut_min,B_mut_mean,B_mut_max] = calculateB(...
            model, exp_data, range_data, mutants{m_idx}, tp(t_idx));
        
        ignore_rxns = {'protein','carbohydrate','lipid','Bio_AA'};
        B_wt_min(findRxnIDs(model,ignore_rxns)) = NaN;
        B_wt_mean(findRxnIDs(model,ignore_rxns)) = NaN;
        B_wt_max(findRxnIDs(model,ignore_rxns)) = NaN;
        B_mut_min(findRxnIDs(model,ignore_rxns)) = NaN;
        B_mut_mean(findRxnIDs(model,ignore_rxns)) = NaN;
        B_mut_max(findRxnIDs(model,ignore_rxns)) = NaN;
        
        %% Step 2 - calculate range for C^mut
        C_mut_min = nan(numel(RXN_IDX),1);
        C_mut_mean = nan(numel(RXN_IDX),1);
        C_mut_max = nan(numel(RXN_IDX),1);
        for i=1:numel(RXN_IDX)
            if C_wt(i) < 1
                C_mut_min(i) = C_wt(i) / ...
                    ( ( B_wt_min(i) / B_mut_min(i) ) * ( 1 - C_wt(i) ) + C_wt(i) );
                C_mut_mean(i) = C_wt(i) / ...
                    ( ( B_wt_mean(i) / B_mut_mean(i) ) * ( 1 - C_wt(i) ) + C_wt(i) );
                C_mut_max(i) = C_wt(i) / ...
                    ( ( B_wt_max(i) / B_mut_max(i) ) * ( 1 - C_wt(i) ) + C_wt(i) );
            end
        end
        
        % find reactions for which minimum is greater than maximum
        err_rxns = model.rxns(RXN_IDX(C_mut_min > C_mut_max));
        for i=1:numel(err_rxns)
            fprintf('%s has  C_mut_min > C_mut_max\n',err_rxns{i})
        end
        
        %% Step 3 - obtain flux distributions for mutant and wild type
        [RDIM,CDIM] = size(model.S);
        
        % constraints
        
        % 2-norm minimization (v^mut - z = v^wt)
        common_rxns = setdiff(model.rxns, ko_rxns{m_idx}, 'stable');
        COM_IDX = findRxnIDs(model,common_rxns);
        
        norm_2_lhs = zeros(numel(COM_IDX),CDIM+numel(COM_IDX));
        norm_2_lhs(:,COM_IDX) = eye(numel(COM_IDX));
        norm_2_lhs(:,CDIM+1:end) = -eye(numel(COM_IDX));
        
        norm_2_rhs = step_1_solution.x(COM_IDX);
        
        norm_2_sense = repelem('=', numel(COM_IDX),1);
        
        % steady state
        st_st_lhs = [model.S zeros(RDIM,numel(COM_IDX))];
        
        st_st_rhs = zeros(RDIM,1);
        
        st_st_sense = repelem('=',RDIM,1);
        
        % biomass ratio wt / mutant
        bio_ratio_lhs = zeros(1,CDIM+numel(COM_IDX));
        bio_ratio_lhs(model.c==1) = biomass_ratio;
        
        bio_ratio_rhs = step_1_solution.f;
        
        bio_ratio_sense = '=';
        
        % contrain v
        nan_c_max = isnan(C_mut_max);
        C_mut_max(nan_c_max&~isnan(C_mut_min)) = 1;
        C_mut_min(isnan(C_mut_min)&~nan_c_max) = 0;
        
        valid_idx = ...
            ~isnan(C_mut_max) & ...
            ~isnan(C_mut_min) & ...
            C_mut_max > 0 & ...
            C_mut_max > C_mut_min; % & abs(C_mut_max-C_mut_min) >= 0.1;
        
        fprintf('%d valid values for C_mut\n',sum(valid_idx))
        
        C_mut_lhs = zeros(2*sum(valid_idx),CDIM+numel(COM_IDX));
        row_counter = 0;
        for i=1:numel(RXN_IDX)
            if valid_idx(i)
                % find all associated kcat values and enzymes
                tmp_rxn_enz_idx = find(model.S(ENZ_MET_IDX,RXN_IDX(i))<0);
                tmp_enz_ids = erase(model.mets(ENZ_MET_IDX(tmp_rxn_enz_idx)),'prot_');
                tmp_draw_rxn_idx = find(endsWith(model.rxns,tmp_enz_ids));
                for j=1:numel(tmp_rxn_enz_idx)
                    
                    % lower bound
                    row_counter = row_counter + 1;
                    C_mut_lhs(row_counter, [RXN_IDX(i) tmp_draw_rxn_idx(j)]) = ...
                        [model.S(ENZ_MET_IDX(tmp_rxn_enz_idx(j)),RXN_IDX(i)) C_mut_min(i)];
                    
                    % upper bound
                    row_counter = row_counter + 1;
                    C_mut_lhs(row_counter, [RXN_IDX(i) tmp_draw_rxn_idx(j)]) = ...
                        -[model.S(ENZ_MET_IDX(tmp_rxn_enz_idx(j)),RXN_IDX(i)) C_mut_max(i)];
                end
            end
        end
        
        C_mut_rhs = zeros(row_counter,1);
        
        C_mut_sense = repelem('<',row_counter,1);
        
        % construct QCP
        
        % 1) without constraints from metabolite abundances
        problem = struct;
        
        problem.A = [st_st_lhs; bio_ratio_lhs; norm_2_lhs];
        problem.rhs = [st_st_rhs; bio_ratio_rhs; norm_2_rhs];
        problem.sense = [st_st_sense; bio_ratio_sense; norm_2_sense];
        
        problem.lb = [model.lb; -repelem(1000,numel(COM_IDX),1)];
        problem.lb(findRxnIDs(model,ko_rxns{m_idx})) = 0;
        
        problem.ub = [model.ub; repelem(1000,numel(COM_IDX),1)];
        problem.ub(findRxnIDs(model,ko_rxns{m_idx})) = 0;
        
        % objective
        problem.Q = sparse(CDIM+numel(COM_IDX),CDIM+numel(COM_IDX));
        problem.Q(CDIM+1:end, CDIM+1:end) = speye(numel(COM_IDX));
        
        % set solver parameters
        params = struct(...
            'Threads', 1,...
            'OutputFlag', 0,...
            'FeasibilityTol', 1e-9,...
            'NumericFocus', 3,...
            'BarCorrectors', 10000 ...
            );
        
        % solve LP using gurobi
        moma_solution = gurobi(problem,params);
        
        if contains(moma_solution.status,'OPTIMAL')
            % calculate ratio B_wt / B_mut
            % B = C / (1-C)
            B_wt_pfba = C_wt ./ (1-C_wt);
            
            C_mut_moma = calculateVBykE(model,moma_solution.x);
            B_mut_moma = C_mut_moma ./ (1-C_mut_moma);
            
            % compare C_mut (MOMA) with C_mut limits from metbolite data
            fig = figure('Visible','off');
            tiledlayout('flow')
            
            nexttile
            hold on
            scatter(C_mut_moma,C_mut_min,'^','k','filled')
            scatter(C_mut_moma,C_mut_max,'v','k','filled');
            scatter(C_mut_moma,C_mut_mean,'.','r')
            xlabel('C_{mut} (MOMA)', 'FontSize', 14)
            ylabel('C_{mut} (calc.)','FontSize', 14)
            legend({'C_{mut}^{min}', 'C_{mut}^{max}', 'C_{mut}^{mean}'},'Location','best',...
                'Box', 'off', 'FontSize', 12)
            
            nexttile
            scatter(C_mut_moma,C_wt,20,'k','filled')
            xlabel('C_{mut} (MOMA)', 'FontSize', 14)
            ylabel('C_{WT} (pFBA)','FontSize', 14)
            
            exportgraphics(fig,[res_dir filesep 'C_values_' mutants{m_idx} ...
                't_' num2str(tp(t_idx)) ...
                '_' l_conds{l_idx} '.png']);
        end
        % 2) including constraints from metabolite abundances
        problem = struct;
        
        problem.A = [st_st_lhs; bio_ratio_lhs; norm_2_lhs; C_mut_lhs];
        problem.rhs = [st_st_rhs; bio_ratio_rhs; norm_2_rhs; C_mut_rhs];
        problem.sense = [st_st_sense; bio_ratio_sense; norm_2_sense; C_mut_sense];
        
        problem.lb = [model.lb; -repelem(1000,numel(COM_IDX),1)];
        problem.lb(findRxnIDs(model,ko_rxns{m_idx})) = 0;
        
        problem.ub = [model.ub; repelem(1000,numel(COM_IDX),1)];
        problem.ub(findRxnIDs(model,ko_rxns{m_idx})) = 0;
        
        % objective
        problem.Q = sparse(CDIM+numel(COM_IDX),CDIM+numel(COM_IDX));
        problem.Q(CDIM+1:end, CDIM+1:end) = speye(numel(COM_IDX));
        
        % set solver parameters
        params = struct(...
            'Threads', 1,...
            'OutputFlag', 0,...
            'FeasibilityTol', 1e-9,...
            'NumericFocus', 3,...
            'BarCorrectors', 10000 ...
            );
        
        % solve LP using gurobi
        step_2_solution = gurobi(problem,params);
        
        if ~contains(step_2_solution.status,'OPTIMAL')
            fprintf('Solver status: %s\n', step_2_solution.status)
        else
            fprintf('Solution found (%s).\n',step_2_solution.status)
            phot_rxns = {'arm_RBC_h' 'arm_RBO_h' 'Tr_CO2m_REV'};
            A_pred_wt = [1 -0.5 -1] * step_1_solution.x(findRxnIDs(model,phot_rxns));
            A_pred_mt = [1 -0.5 -1] * step_2_solution.x(findRxnIDs(model,phot_rxns));
            fprintf(...
                'CO2 assimilation rates [umol/m2/s]: WT: %.4g %s: %.4g\n',...
                1000*LMA*A_pred_wt / 3600,... mmol/gDW/h => umol/m2/s
                mutants{m_idx},...
                1000*LMA*A_pred_mt / 3600 ...
                )
            fprintf('Oxygenation fluxes: [mmol/gDW/h]: WT: %.4g %s: %.4g\n',...
                step_1_solution.x(findRxnIDs(model,'arm_RBO_h')),...
                mutants{m_idx},...
                step_2_solution.x(findRxnIDs(model,'arm_RBO_h'))...
                )
            % rank reactions by distance
            rxn_dist = step_2_solution.x(end-numel(COM_IDX)+1:end);
            [rxn_dist_sorted,rxn_dist_order] = sort(rxn_dist,'descend');
            rxns_by_dist = model.rxns(COM_IDX(rxn_dist_order));
            subsyst_by_dist = model.subSystems(COM_IDX(rxn_dist_order));
            top_50_rxns = rxns_by_dist(1:50);
            top_50_subsyst = subsyst_by_dist(1:50);
            top_50_subsyst = [top_50_subsyst{:}]';
            
            % write results to file
            writetable(...
                [cell2table(...
                [model.rxns(COM_IDX(rxn_dist_order)),...
                model.rxnNames(COM_IDX(rxn_dist_order)),...
                model.subSystems(COM_IDX(rxn_dist_order))],...
                'VariableNames', {'RXN_ID', 'RXN_NAME', 'SUBSYSTEM'}...
                ),...
                array2table(...
                [rxn_dist_sorted,...
                step_1_solution.x(COM_IDX(rxn_dist_order)),...
                step_2_solution.x(COM_IDX(rxn_dist_order))],...
                'VariableNames', {'dist','wt_flux','mut_flux'}...
                )],...
                [res_dir filesep 'flux_dist_norm_2_' mutants{m_idx} ...
                't_' num2str(tp(t_idx)) ...
                '_' l_conds{l_idx} '.csv'],...
                'WriteRowNames', true...
                );
        end
        
        if contains(step_2_solution.status,'OPTIMAL')
            %% STEP 4 - Run flux variability analysis at optimal flux sum
            % remove quadratic objective
            problem_va = rmfield(problem,'Q');
            % initialize linear objective
            problem_va.obj = zeros(size(problem.A,2),1);
            % add quadratic constraint to limit the sum of squared differences
            % to the (sub-)optimal solution
            problem_va.quadcon.Qc = sparse(size(problem.A,2),size(problem.A,2));
            problem_va.quadcon.Qc(CDIM+1:end,CDIM+1:end) = speye(numel(COM_IDX));
            problem_va.quadcon.rhs = sparse(1.1*step_2_solution.objval);
            problem_va.quadcon.q = sparse(size(problem.A,2),1);
            
            % set LB for growth to 90% of optimum in previous step
            problem_va.lb(model.c==1) = 0.9*step_2_solution.x(model.c==1);
            
            minFluxMUT = zeros(numel(COM_IDX),1);
            maxFluxMUT = zeros(numel(COM_IDX),1);
            fprintf('Starting variability analysis...\n')
            parfor i=1:numel(COM_IDX)
                if i > 1 && mod(i,1000)==1
                    fprintf('Processed %d reactions ...\n',i-1)
                end
                tmp_problem = problem_va;
                % maximization
                tmp_problem.obj(COM_IDX(i)) = -1;
                sol = gurobi(tmp_problem,params);
                if contains(sol.status,'OPTIMAL')
                    maxFluxMUT(i) = sol.x(COM_IDX(i));
                end
                
                % minimization
                tmp_problem.obj(COM_IDX(i)) = 1;
                sol = gurobi(tmp_problem,params);
                if contains(sol.status,'OPTIMAL')
                    minFluxMUT(i) = sol.x(COM_IDX(i));
                end
            end
            
            % find non-overlapping reactions
            non_overlapping = false(numel(COM_IDX),1);
            for i=1:numel(COM_IDX)
                if minFluxMUT(i) > maxFluxWT(COM_IDX(i)) + 1e-8
                    non_overlapping(i) = 1;
                elseif minFluxWT(COM_IDX(i)) > maxFluxMUT(i) + 1e-8
                    non_overlapping(i) = 1;
                end
            end
            
            writetable(...
                [cell2table(...
                [model.rxns(COM_IDX(non_overlapping)),...
                model.rxnNames(COM_IDX(non_overlapping)),...
                model.subSystems(COM_IDX(non_overlapping))],...
                'VariableNames', {'ReactionID','ReactionName','ReactionSubSystem'}...
                ),...
                array2table(...
                [minFluxWT(COM_IDX(non_overlapping)) maxFluxWT(COM_IDX(non_overlapping)),...
                minFluxMUT(non_overlapping) maxFluxMUT(non_overlapping)],...
                'VariableNames', {...
                'minFluxWT','maxFluxWT','minFluxMUT','maxFluxMUT'...
                })],...
                [res_dir filesep 'diff_rxns_' mutants{m_idx} ...
                't_' num2str(tp(t_idx)) ...
                '_' l_conds{l_idx} '.csv']...
                );
        end
    end
end