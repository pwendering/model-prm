% This script generates heatmaps for reactions associated with
% photorespiration, TCA cycle, and additional reactions in the AraCore model.
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

plot_rxns = [pr_rxns; tca_cycle_rxns; {'GlyHMT_c'; 'GCEADH_h'; 'GGAT_h'; 'GCAO_h'}];
plot_rxn_idx = findRxnIDs(model, plot_rxns);

% flux sampling data specifications
light_conditions = {'ml', 'fl'};
ph_ub = [286 191];
genotypes = {'Col-0', 'ggt1-1', 'ggt1-2', 'hpr1-1', 'hpr1-2'};
timepoints = {'-21'};

flux_av = cell(1, numel(plot_rxns));

for i = 1:numel(plot_rxns)
    flux_av{i} = cell(1, numel(light_conditions));
    for j = 1:numel(light_conditions)
        flux_av{i}{j} = nan(numel(genotypes), numel(timepoints));
    end
end

flux_sd = flux_av;

for l = 1:numel(light_conditions)
    fprintf('%s\t', light_conditions{l})
    for i = 1:numel(genotypes)
        
        geno = genotypes{i};
        fprintf('%s\t', geno)
        
        for j = 1:numel(timepoints)
            
            tp = timepoints{j};
            fprintf('%s\t', tp)
            
            file_name = fullfile(res_dir, light_conditions{l},...
                ['flux_samples_', geno, '_' light_conditions{l} '_t_', tp,...
                '_pHUB_' num2str(ph_ub(l)) '.csv']);
            
            if isfile(file_name)
                
                % read flux samples
                tmp_flux = readmatrix(file_name);
                tmp_flux(tmp_flux>0&tmp_flux<1e-9) = 0;
                fprintf('Number of failed samplings: %d\n',sum(~any(tmp_flux)))
                tmp_flux(:, ~any(tmp_flux)) = [];
                
                for k = 1:numel(plot_rxns)
                    rxn_fluxes = tmp_flux(plot_rxn_idx(k), :);
                    flux_av{k}{l}(i,j) = mean(rxn_fluxes, 'omitnan');
                    flux_sd{k}{l}(i,j) = std(rxn_fluxes, 'omitnan');
                end
            end
        end
    end
end

%% Create heatmaps
% get color range
tmp1 = [flux_av{:}];
tmp2 = round([tmp1{:}], 9);
tmp2(tmp2==-0) = 0;
crange = [floor(min(min(tmp2))) ceil(max(max(tmp2)))];

% load custom colormap
load('colormap_blue_red.mat')
for i=1:numel(plot_rxns)

    mat = round([flux_av{i}{1}'; flux_av{i}{2}'], 9);
    mat(mat==-0) = 0;
    mat(mat==0) = NaN;
    
    heatmap(mat, 'GridVisible', 'off', 'CellLabelColor','none',...
        'Colormap', colormap_blue_red,...
        'ColorLimits', crange,...
        'XDisplayLabels', repmat({''}, size(mat,2), 1),...
        'YDisplayLabels', repmat({''}, size(mat,1), 1),...
        'fontsize',14,...
        'colorbarvisible', 'off')
    set(gcf, 'OuterPosition', [353.6667  265.0000  658.6666  434.0000])
    exportgraphics(gcf, [res_dir filesep 'heatmaps' filesep plot_rxns{i} '.png'])
end

% only plot colorbar
plot(NaN,NaN)
colormap(colormap_blue_red)
cb = colorbar('north');
set(gca, 'Visible', 'off')
pause(.1)  % to prevent xticks to be placed on top of the cb for some reason
set(cb, 'Position', [0.1325    0.5    0.6605    0.1])
set(cb, 'Ticks', 0:.25:1)
set(cb, 'TickLabels', round(linspace(floor(crange(1)), ceil(crange(2)), 5),2))
set(cb, 'TickDirection', 'out')
set(cb, 'FontSize', 12)
cb.Title.String = 'log_2 flux [mmol/gDW/h]';
exportgraphics(gcf, [res_dir filesep 'heatmaps' filesep 'colorbar.jpg'])