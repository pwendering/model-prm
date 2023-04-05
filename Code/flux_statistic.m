% Statistical comparison of average flux within and across genotypes and
% light conditions

clear;clc;close all

data_dir = '../Data';
res_dir = '../Results';

% read model and find reactions of interest
model_file = fullfile(data_dir,'AraCore-updated-rev.mat');
load(model_file);

% update GPR for GGT1 and HPR1; add HPR2 reaction (cytosol)
model = updateModelHprGgt(model);

% specify reactions for which heatmaps should be plotted
model.subSystems{findRxnIDs(model, 'GCEADH_c')} = {'photorespiration'};

pr_rxns = model.rxns(cellfun(@(x)contains(x, 'photorespiration'),...
    model.subSystems));

tca_cycle_rxns = model.rxns(cellfun(@(x)contains(x, 'tricarboxylic acid cycle'),...
    model.subSystems));

rxns = [pr_rxns; tca_cycle_rxns; {'GlyHMT_c'; 'GCEADH_h'; 'GGAT_h';...
    'GCAO_h'; 'MalDH3_h'; 'GluSFd_h'; 'GluSNAD_h'; 'GlnS_h'}];
rxn_idx = findRxnIDs(model, rxns);

% flux sampling data specifications
light_conditions = {'ml', 'fl'};
ph_ub = [286 191];
genotypes = {'Col-0', 'ggt1-1', 'ggt1-2', 'hpr1-1', 'hpr1-2'};
timepoints = {'-21'};

fluxes = cell(1, numel(rxns));

for i = 1:numel(rxns)
    fluxes{i} = nan(numel(light_conditions)*1000, numel(genotypes));
end

for i = 1:numel(genotypes)
    
    geno = genotypes{i};
    
    for j = 1:numel(light_conditions)
        lc = light_conditions{j};
        
        
        file_name = fullfile(res_dir, lc,...
            ['flux_samples_', geno, '_' lc '_t_', timepoints{1},...
            '_pHUB_' num2str(ph_ub(j)) '.csv']);
        
        
        if isfile(file_name)
            
            % read flux samples
            tmp_flux = readmatrix(file_name);
            tmp_flux(tmp_flux>0&tmp_flux<1e-9) = 0;
            fprintf('Number of failed samplings: %d\n',sum(~any(tmp_flux)))
            
            for k = 1:numel(rxns)
                fluxes{k}((j-1)*1000+1:j*1000,i) = tmp_flux(rxn_idx(k), :);
            end
        end
    end
end

%% Pairwise comparisons of genotypes within each light condition
ml_row_idx = 1:1000;
fl_row_idx = 1001:2000;
p_values = cell(numel(rxns),1);
for i=1:numel(rxns)
    disp(rxns{i})
    p_mat = nan(numel(genotypes));
    for j=1:numel(genotypes)-1
        for k=j+1:numel(genotypes)
            % ML (CL)
            f1 = fluxes{i}(ml_row_idx,j);
            f2 = fluxes{i}(ml_row_idx,k);
            if any(~isnan(f1)) && any(~isnan(f2))
                p_mat(j,k) = ranksum(f1, f2);
            end
            % FL
            f1 = fluxes{i}(fl_row_idx,j);
            f2 = fluxes{i}(fl_row_idx,k);
            if any(~isnan(f1)) && any(~isnan(f2))
                p_mat(k,j) = ranksum(f1, f2);
            end
        end
    end
    p_values{i} = p_mat;
end

% correct p-values
[~,~,~,p_values_adj] = fdr_bh(cell2mat(p_values));

% set all diagonal values to one
p_values_adj(repmat(logical(eye(numel(genotypes))), numel(rxns), 1)) = 1;

% round P-values to 6th digit
p_values_adj = round(p_values_adj, 6);

fid = fopen(fullfile('..', 'Results', 'stats_fluxes_pairwise.tsv'), 'w+');
fprintf(fid, ['Reaction ID/name/EC number\tFL\\CL\t' strjoin(genotypes, '\t') '\n']);
for i=1:numel(rxns)
    rxn_idx = findRxnIDs(model, rxns(i));
    rxn_name = model.rxnNames{rxn_idx};
    rxn_ec = model.rxnECNumbers{rxn_idx};
    for j=1:numel(genotypes)
        if j==1
            fprintf(fid, '%s\t', rxns{i});
        elseif j==2
            fprintf(fid, '%s\t', rxn_name);
        elseif j==3
            fprintf(fid, '%s\t', rxn_ec);
        else
            fprintf(fid, '\t');
        end
        fprintf(fid, ['%s\t' strjoin(repelem({'%.4g'}, numel(genotypes), 1), '\t') '\n'],...
            genotypes{j}, p_values_adj((i-1)*numel(genotypes)+j,:));
    end
end
fclose(fid);

%% Test for difference in flux in different light conditions for each genotype
ml_row_idx = 1:1000;
fl_row_idx = 1001:2000;
p_values = cell(numel(rxns),1);
for i=1:numel(rxns)
    disp(rxns{i})
    p_mat = nan(1, numel(genotypes));
    for j=1:numel(genotypes)
            % ML (CL)
            f1 = fluxes{i}(ml_row_idx,j);
            % FL
            f2 = fluxes{i}(fl_row_idx,j);
            if any(~isnan(f1)) && any(~isnan(f2))
                p_mat(j) = ranksum(f1, f2);
            end
    end
    p_values{i} = p_mat;
end

% correct p-values
[~,~,~,p_values_adj] = fdr_bh(cell2mat(p_values));

% round P-values to 6th digit
p_values_adj = round(p_values_adj, 6);

fid = fopen(fullfile('..', 'Results', 'stats_fluxes_between_conditions.tsv'), 'w+');
fprintf(fid, 'Reaction ID\tReaction Name\tEC Number\t\tGenotype\t\t\t\n');
fprintf(fid, ['\t\t\t' strjoin(genotypes, '\t') '\n']);
for i=1:numel(rxns)
    rxn_idx = findRxnIDs(model, rxns(i));
    rxn_name = model.rxnNames(rxn_idx);
    rxn_ec = model.rxnECNumbers(rxn_idx);
    fprintf(fid, strjoin([rxns(i) rxn_name rxn_ec], '\t'));
    fprintf(fid, ['\t' strjoin(repelem({'%.4g'}, numel(genotypes), 1), '\t') '\n'],...
        p_values_adj(i,:));
    fprintf(fid, '\n');
end
fclose(fid);


