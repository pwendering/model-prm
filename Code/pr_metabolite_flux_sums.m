% This script calculates flux sums for all metabolites that are associated
% with photorespiratory reactions in the AraCore model.
clear;clc;close all

data_dir = '../Data';
res_dir = '../Results';

% read model and find reactions of interest
model_file = fullfile(data_dir,'AraCore-updated-rev.mat');
load(model_file);

% update GPR for GGT1 and HPR1; add HPR2 reaction (cytosol)
model = updateModelHprGgt(model);

% flux sampling data specifications
light_conditions = {'ml', 'fl'};
ph_ub = [286 191];
genotypes = {'Col-0', 'ggt1-1', 'ggt1-2', 'hpr1-1', 'hpr1-2'};
timepoints = {'-21'};

% take subset of measured metabolites that are involved in photorespiration
pr_rxns = model.rxns(cellfun(@(x)contains(x, 'photorespiration'),...
    model.subSystems));
mets = unique(findMetsFromRxns(model, pr_rxns));
met_names = model.metNames(findMetIDs(model, mets));
mets = strtok(mets, '[');

% initialize arrays for average flux sums and standard deviations
met_fs_av = cell(1, numel(mets));
met_fs_sd = cell(1, numel(mets));
met_turnover_raw = cell(1, numel(mets));

for i=1:numel(mets)
    % initialize cell array for current metabolite
    met_fs_av{i} = cell(1, numel(light_conditions));
    met_fs_sd{i} = cell(1, numel(light_conditions));
    
    % initialize flux sums arrays for each light condition
    for j=1:numel(light_conditions)
        met_fs_av{i}{j} = nan(numel(genotypes), numel(timepoints));
        met_fs_sd{i}{j} = nan(numel(genotypes), numel(timepoints));
    end
end

for l = 1:numel(light_conditions)
    % metabolite abundances
    met_abundance_file = fullfile(data_dir,['met_' light_conditions{l} '.csv']);
    
    for i = 1:numel(genotypes)
        
        geno = genotypes{i};
        
        for j = 1:numel(timepoints)
            
            tp = timepoints{j};
            
            flux_sampling_file = fullfile(res_dir, light_conditions{l},...
                ['flux_samples_', geno, '_' light_conditions{l} '_t_', tp,...
                '_pHUB_' num2str(ph_ub(l)) '.csv']);
            
            if isfile(flux_sampling_file)
                
                % read flux sampling results
                tmp_flux = readmatrix(flux_sampling_file);
                tmp_flux(tmp_flux>0&tmp_flux<1e-9) = 0;
                fprintf('Number of failed samplings: %d\n',sum(~any(tmp_flux)))
                tmp_flux(:, ~any(tmp_flux)) = [];
                
                for m = 1:numel(mets)
                    tmp_model_mets = reduceCell(regexp(model.mets, ['^' mets{m} '\[.\]$'], 'match'));
                    tmp_met_idx = findMetIDs(model, tmp_model_mets);
                    tmp_rxns = findRxnsFromMets(model, tmp_model_mets);
                    tmp_rxn_idx = findRxnIDs(model, tmp_rxns);
                                        
                    tmp_flux_sum_mat = zeros(size(tmp_flux,2),numel(tmp_met_idx));
                    for I=1:numel(tmp_met_idx)
                        for J=1:numel(tmp_rxn_idx)
                            for K=1:size(tmp_flux,2)
                                tmp_flux_sum_mat(K,I) = tmp_flux_sum_mat(K,I) + abs(model.S(tmp_met_idx(I), tmp_rxn_idx(J)) * tmp_flux(tmp_rxn_idx(J),K));
                            end
                        end
                    end
                    tmp_flux_sums = 0.5 * sum(tmp_flux_sum_mat,2);
                    
                    met_fs_av{m}{l}(i,j) = mean(tmp_flux_sums, 'omitnan');
                    met_fs_sd{m}{l}(i,j) = std(tmp_flux_sums, 'omitnan');
                end 
            end
        end
    end
end

%% Plot results as barplots
for i=1:numel(mets)
    createBarplot(met_fs_av{i}, met_fs_sd{i}, genotypes,...
        timepoints, ['Flux sum of ' met_names{i}, ' [mmol/gDW/h'],...
        light_conditions, [res_dir filesep lower(met_names{i}) '_flux_sum.jpg'])    
end

function createBarplot(phenotype_av, phenotype_sd, genotypes, timepoints, y_label,...
    light_conditions, fig_filename)
fig = figure;
set(gcf, 'OuterPosition', [164.3333  175.6667  416.6667  444.6667])
n_lc = numel(light_conditions);

if isempty(phenotype_sd)
    eb = false;
else
    eb = true;
end

% get maximum value to set common y-axis limit for all barplots
y_max = 0;
y_log_flag = false;
for l = 1:n_lc
    
    av = 60.*phenotype_av{l};
    sd = 60.*phenotype_sd{l};

    if any(any(av<1e-4)) || any(any(av>1e4))
        y_log_flag = true;
    end
    
    if eb
        tmp_max = max(max(av+sd));
    else
        tmp_max = max(max(av));
    end
    
    y_max = max(y_max, tmp_max);
end
y_max = 1.2*y_max;

for l = 1:n_lc
    subplot(n_lc,1,l)
    
    av = 60.*phenotype_av{l};
    sd = 60.*phenotype_sd{l};

    bars = bar(av', 'facecolor', 'flat');
    bars.CData = lines(numel(genotypes));
    
    if eb
        hold on
        addErrorBars(bars,timepoints,av,sd)
    end
    
    xticklabels(timepoints)
    
    if l == n_lc
        xticklabels(genotypes)
    else
        xticklabels({''})
    end
    
    ylim([0 y_max])
    
    if y_log_flag
        set(gca, 'YScale', 'log')
        pause(3)
        n_y_ticks = numel(yticks);
        if n_y_ticks < 2
            yticks(logspace(round(log10(min(min(av)))),round(log10(y_max)),5))
            yticklabels(arrayfun(@(x)sprintf('%.1e',x),logspace(round(log10(min(min(av)))),round(log10(y_max)),5),'un',0))
        end
    end
    
    text(.007, 0.92, strrep(light_conditions{l}, '_', '->'), 'units', 'normalized', 'FontSize', 10)
    
end

han = axes(fig, 'visible', 'off');
han.Title.Visible = 'on';
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
ylabel(han, y_label, 'FontSize', 14);

if y_log_flag
    han2 = get(gca,'YLabel');
    set(han2,'Position', get(han2,'Position') - [0.01 0 0])
end

exportgraphics(gcf, fig_filename)

end


function addErrorBars(bars,timepoints,av,sd)
for i=1:numel(timepoints)
    for j=1:numel(bars)
        errorbar(i+bars(j).XOffset,...
            av(j,i),...
            sd(j,i),...
            'Color', 'k')
    end
end

ylim = get(gca, 'YLim');
if ylim(1) < 0
    ylim(1) = 0;
    set(gca, 'YLim', ylim)
end
end
