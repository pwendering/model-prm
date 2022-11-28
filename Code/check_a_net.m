% check if net CO2 assimilation constraints are satisfied
clear;clc;close all

data_dir = '../Data';
res_dir = '../Results';

% read model and find reactions of interest
model_file = fullfile(data_dir,'AraCore-updated-rev.mat');
load(model_file);

add_hpr2_reaction

v_o_rxn = 'RBO_h';
v_o_idx = findRxnIDs(model, v_o_rxn);
v_c_rxn = 'RBC_h';
v_c_idx = findRxnIDs(model, v_c_rxn);

rd_rxn = 'Tr_CO2m';
rd_idx = findRxnIDs(model, rd_rxn);

% flux sampling data specifications
light_conditions = {'ml', 'fl'};
ph_ub = [286 191];
genotypes = {'Col-0', 'ggt1-1', 'ggt1-2', 'hpr1-1', 'hpr1-2'};
timepoints = {'-21'};


a_net_mean = cell(1, numel(light_conditions));
a_net_sd = cell(1, numel(light_conditions));

% load experimental data
a_net_ml = readtable(fullfile(data_dir, 'A_ml.csv'));
a_net_ml.Genotype(ismember(a_net_ml.Genotype, 'Col0')) = {'Col-0'};
a_net_ml = a_net_ml(cellfun(@(x)find(ismember(a_net_ml.Genotype,x)), genotypes),:);
a_net_fl = readtable(fullfile(data_dir, 'A_fl.csv'));
a_net_fl.Genotype(ismember(a_net_fl.Genotype, 'Col0')) = {'Col-0'};
a_net_fl = a_net_fl(cellfun(@(x)find(ismember(a_net_fl.Genotype,x)), genotypes),:);


for l = 1:numel(light_conditions)
    disp(light_conditions{l})
    % initialize matrices for metabolic phenotypes (genotype x time after shift)
    a_net_mean{l} = nan(numel(genotypes), numel(timepoints));
    a_net_sd{l} = nan(numel(genotypes), numel(timepoints));
    
    
    for i = 1:numel(genotypes)
        
        geno = genotypes{i};
        disp(geno)
        for j = 1:numel(timepoints)
            
            tp = timepoints{j};
            disp(tp)
            
            file_name = fullfile(res_dir, light_conditions{l},...
                ['flux_samples_', geno, '_' light_conditions{l} '_t_', tp,...
                '_pHUB_' num2str(ph_ub(l)) '.csv']);
            
            if isfile(file_name)
                
                % read flux samples
                tmp_flux = readmatrix(file_name);
                tmp_flux(tmp_flux>0&tmp_flux<1e-9) = 0;
                fprintf('Number of failed samplings: %d\n',sum(~any(tmp_flux)))
                tmp_flux(:, ~any(tmp_flux)) = [];
                
                % Net CO2 assimilation rate
                a_net_mean{l}(i,j) = mean(sum([1 -0.5 1]'.*tmp_flux([v_c_idx v_o_idx rd_idx],:)), 'omitnan');
                a_net_sd{l}(i,j) = std(sum([1 -0.5 1]'.*tmp_flux([v_c_idx v_o_idx rd_idx],:)), 'omitnan');
            end
        end
    end
end

%% scatter plot with error bars
colors = lines(numel(genotypes));

figure
tiledlayout(1,2)
nexttile
hold on

a_rel_exp = a_net_ml.mean_A / max(a_net_ml.mean_A);
a_rel_exp_sd = a_net_ml.sd_A / max(a_net_ml.mean_A);
a_rel_pred = a_net_mean{1} / max(a_net_mean{1});
a_rel_pred_sd = a_net_sd{1} / max(a_net_mean{1});

line([0 1.2],...
    [0 1.2], 'Color', [.5 .5 .5])
for i=1:numel(genotypes)
    h(i) = errorbar(a_rel_exp(i), a_rel_pred(i), a_rel_pred_sd(i),...
        a_rel_pred_sd(i), a_rel_exp_sd(i), a_rel_exp_sd(i),...
        '.', 'Color', colors(i,:), 'MarkerSize', 20);
end

xlabel('A_{net} experimental (scaled)')
ylabel('A_{net} predicted (scaled)')
disp('Spearman Correlation (ML)')
rho_s = corr(a_net_ml.mean_A,  a_net_mean{1}, 'type', 'Spearman');
disp(rho_s)
text(0.05, 0.93, 'C', 'Units', 'normalized', 'FontSize', 14)
set(gca, 'FontSize', 14, 'XLim', [.4 1.1], 'YLim', [.4 1.1], 'Box', 'on',...
    'LineWidth', 1.3)

nexttile
hold on

a_rel_exp = a_net_fl.mean_A / max(a_net_fl.mean_A);
a_rel_exp_sd = a_net_fl.sd_A / max(a_net_fl.mean_A);
a_rel_pred = a_net_mean{2} / max(a_net_mean{2});
a_rel_pred_sd = a_net_sd{2} / max(a_net_mean{2});

line([0 1.2],...
    [0 1.2], 'Color', [.5 .5 .5])
for i=1:numel(genotypes)
    h(i) = errorbar(a_rel_exp(i), a_rel_pred(i), a_rel_pred_sd(i), a_rel_pred_sd(i), a_rel_exp_sd(i),a_rel_exp_sd(i),...
        '.', 'Color', colors(i,:), 'MarkerSize', 20);
end

legend(h, genotypes, 'Box', 'off', 'Location', 'southeast')

xlabel('A_{net} experimental (scaled)')

disp('Spearman Correlation (FL)')
rho_s = corr(a_net_fl.mean_A,  a_net_mean{2}, 'type', 'Spearman');
disp(rho_s)
text(0.05, 0.93, 'D', 'Units', 'normalized', 'FontSize', 14)
set(gca, 'FontSize', 14, 'XLim', [.4 1.1], 'YLim', [.4 1.1], 'Box', 'on',...
    'LineWidth', 1.3)
set(gcf, 'OuterPosition', [134.3333  161.0000  935.3333  475.3333])
exportgraphics(gcf, fullfile(res_dir, 'a_net_check.png'))


%% Compare theoretical A_net for mutant with actual flux values
fprintf('\nChecking ratios ML\n')
for i=2:numel(genotypes)
    fprintf('Mutant: %s\n', genotypes{i})
    % get biomass ratio
    a_net_ratio = a_net_ml.mean_A(1) / ...
        a_net_ml.mean_A(strcmp(a_net_ml.Genotype, genotypes{i}));
    wt_a_netflux = a_net_mean{1}(1);
    
    theoretical_mutant_a_net = wt_a_netflux / a_net_ratio;
    actual_mutant_a_net = a_net_mean{1}(i);
    
    fprintf('Theoretical value: %.4g; Sampling mean: %.4g\n',...
        theoretical_mutant_a_net, actual_mutant_a_net)
    fprintf('Allowed range: %.4g - %.4g\n',....
        theoretical_mutant_a_net*(1-1e-3), theoretical_mutant_a_net*(1+1e-3))
    
end
