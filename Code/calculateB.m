function [B_wt_min,B_wt_mean,B_wt_max,B_mut_min,B_mut_mean,B_mut_max] = calculateB(model, exp_data, range_data, mut_id, tp)
%% [B_wt_min,B_wt_mean,B_wt_max,B_mut_min,B_mut_mean,B_mut_max] = calculateB(model, exp_data, range_data, mut_id, tp)
% Calculate ratio based on Michaelis-Menten kinetics for each reaction i:
% 
% B = prod( (x_ij / KM_ij ) ^ alpha_ij), 
% 
% where x is the metabolite concentration, KM is the Michaelis constant and
% alpha is the stoichiometric coefficient for each metabolite j in each
% reaction i. 
% 
% INPUT:
%   struct model:
%       metabolic model
%   struct exp_data:
%       metabolite abundances of current experiment 
%           .met_ids:   metabolite IDs (as in model, without compartment info)
%           .met_av:    average abundance
%           .met_sd:    standaerd deviation
%   struct range_data:
%           .range_ids:         metabolite IDs (as in model, without compartment info)
%           .met_range_av:      average of metabolite abundances over multiple
%                               experiments
%           .met_range_sd:      standard deviation of metabolite abundances over multiple
%                               experiments
%           .min_conc_limit:    overall minimum abundance for metabolite
%                               (used if unmeasured and no range available)
%           .max_conc_limit:    overall maximum abundance for metabolite
%                               (used if unmeasured and no range available)
%   char mut_id:
%       mutant ID (lookup in metabolite data)
%   double tp:
%       time point after shift [-21 3 27 75 123]
% OUTPUT:
%   double B_wt_min, B_wt_mean, B_wt_max, B_mut_min, B_mut_mean, B_mut_max: 
%       minimum and maximum values B for wild-type and mutant

% Col-0 is wild type
wt_id = 'Col-0';

% dummy value for KM
KM = 1;

% create variables from metabolite data structs
[meas_ids, met_av, met_sd] = deal(exp_data.meas_ids, exp_data.met_av, exp_data.met_sd);
[range_ids, met_range_av, met_range_sd, min_conc_limit, max_conc_limit] = ...
    deal(range_data.range_ids, range_data.met_range_av, range_data.met_range_sd,...
    range_data.max_conc_limit, range_data.min_conc_limit);

% get indices of metabolic reactions
rxn_idx = getRxnEnzIdx(model);

% initialize output arrays
B_wt_min = zeros(numel(rxn_idx),1);
B_wt_mean = zeros(numel(rxn_idx),1);
B_wt_max = zeros(numel(rxn_idx),1);
B_mut_min = zeros(numel(rxn_idx),1);
B_mut_mean = zeros(numel(rxn_idx),1);
B_mut_max = zeros(numel(rxn_idx),1);
for i=1:numel(rxn_idx)
    
    % determine indices of metabolites taking part in current reaction
    rxn_met_idx = find(...
        model.S(:,rxn_idx(i)) ~= 0 & ...
        ~contains(model.mets, 'prot_')...
        );
    
    % IDs of metabolites in current reaction
    rxn_mets = model.mets(rxn_met_idx);
    
    % take care of arm reactions
    if any(contains(rxn_mets,'pmet_'))
        
        % get index and coefficient of pseudometabolite
        tmp_pmet_idx = rxn_met_idx(contains(rxn_mets,'pmet_'));
        pmet_coeff = model.S(tmp_pmet_idx,rxn_idx(i));
        
        if pmet_coeff < 0
            % pmet is consumed, find substrates of arm reaction
            arm_rxn_idx = find(model.S(tmp_pmet_idx,:) > 0);
            adj_mets = model.mets(any(model.S(:,arm_rxn_idx) < 0,2));
        else
            % pmet is produced, find products of consuming reaction
            arm_rxn_idx = find(model.S(tmp_pmet_idx,:) < 0);
            adj_mets = model.mets(any(model.S(:,arm_rxn_idx) > 0,2));
        end
    else
        arm_rxn_idx = 0;
        adj_mets = {};
    end
    
    % remove pmet and add ID / index of adjacent metabolites
    rxn_met_idx = rxn_met_idx(~contains(rxn_mets,'pmet_'));
    rxn_met_idx = [rxn_met_idx; findMetIDs(model,adj_mets)];
    
    rxn_mets = rxn_mets(~contains(rxn_mets,'pmet_'));
    rxn_mets = [rxn_mets; adj_mets];
    
    % stoichiometric coefficients
    rxn_met_coeff = full(model.S(rxn_met_idx,rxn_idx(i)));
    
    % find match in data from current experiment
    trunc_mets = regexprep(rxn_mets,'_\w$','');
    meas_match = ismember(trunc_mets, meas_ids);
    
    % find match in metabolite abundance range dataset
    range_match = ismember(trunc_mets, range_ids);
    range_match(meas_match) = false;
    
    % initialize concentrations
    concs = nan(numel(rxn_mets),1);
    sd = nan(numel(rxn_mets),1);
    
    % wild type (Col-0)
    concs(meas_match) = met_av{...
        met_av.time==tp & strcmp(met_av.genotype,wt_id),...
        2+cellfun(@(x)find(ismember(meas_ids, x)),trunc_mets(meas_match))...
        };
    sd(meas_match) = met_sd{...
        met_av.time==tp & strcmp(met_av.genotype,wt_id),...
        2+cellfun(@(x)find(ismember(meas_ids, x)),trunc_mets(meas_match))...
        };
    
    concs(range_match) = met_range_av(...
        cellfun(@(x)find(ismember(range_ids, x)),trunc_mets(range_match))...
        );
    sd(range_match) = met_range_sd(...
        cellfun(@(x)find(ismember(range_ids, x)),trunc_mets(range_match))...
        );
    
    conc_min = concs-sd;
    conc_min(isnan(conc_min)|conc_min<0) = min_conc_limit;
    
    conc_max = concs+sd;
    conc_max(isnan(conc_max)) = max_conc_limit;
    
    B_wt_min(i) = prod((conc_min ./ KM ) .^ rxn_met_coeff);
    B_wt_mean(i) = prod((concs ./ KM ) .^ rxn_met_coeff);
    B_wt_max(i) = prod((conc_max ./ KM ) .^ rxn_met_coeff);
    
    if arm_rxn_idx > 0
        B_wt_min(arm_rxn_idx) = B_wt_min(i);
        B_wt_mean(arm_rxn_idx) = B_wt_mean(i);
        B_wt_max(arm_rxn_idx) = B_wt_max(i);
    end
    
    % mutant
    
    % initialize concentrations
    concs = nan(numel(rxn_mets),1);
    sd = nan(numel(rxn_mets),1);
    
    concs(meas_match) = met_av{...
        met_av.time==tp & strcmp(met_av.genotype,mut_id),...
        2+cellfun(@(x)find(ismember(meas_ids, x)),trunc_mets(meas_match))...
        };
    sd(meas_match) = met_sd{...
        met_av.time==tp & strcmp(met_av.genotype,mut_id),...
        2+cellfun(@(x)find(ismember(meas_ids, x)),trunc_mets(meas_match))...
        };
    
    concs(range_match) = met_range_av(...
        cellfun(@(x)find(ismember(range_ids, x)),trunc_mets(range_match))...
        );
    sd(range_match) = met_range_sd(...
        cellfun(@(x)find(ismember(range_ids, x)),trunc_mets(range_match))...
        );
    
    conc_min = concs-sd;
    conc_min(isnan(conc_min)|conc_min<0) = min_conc_limit;
    
    conc_max = concs+sd;
    conc_max(isnan(conc_max)) = max_conc_limit;
    
    B_mut_min(i) = prod((conc_min ./ KM ) .^ rxn_met_coeff);
    B_mut_mean(i) = prod((concs ./ KM ) .^ rxn_met_coeff);
    B_mut_max(i) = prod((conc_max ./ KM ) .^ rxn_met_coeff);
    
    if arm_rxn_idx > 0
        B_mut_min(arm_rxn_idx) = B_mut_min(i);
        B_mut_mean(arm_rxn_idx) = B_mut_mean(i);
        B_mut_max(arm_rxn_idx) = B_mut_max(i);
    end
end
end