% check for violations of growth ratio constraints
clear;clc;close all

% define data and output directory
data_dir = '../Data';
res_dir = '../Results';

% read model and find reactions of interest
model_file = fullfile(data_dir,'AraCore-updated-rev.mat');
load(model_file);

% add HPR2 reaction (cytosol) and update gene association 
model = updateModelHprGgt(model);

bio_rxn_id = 'Bio_opt';
bio_rxn_idx = findRxnIDs(model, bio_rxn_id);

% dry weights
dw_ml = readtable(...
    fullfile(data_dir, 'dw_ml_34.csv'),...
    'ReadRowNames', true...
    );
das_ml = 34;
dw_fl = readtable(...
    fullfile(data_dir, 'dw_fl_43.csv'),...
    'ReadRowNames', true...
    );
das_fl = 43;


% growth duration
t_growth = das_ml * 24; % [h]

% calculate growth rate
dw_0 = 20 / 1e6;
mu_ml = (log(dw_ml.mean_DW) - log(dw_0)) ./ t_growth;
mu_wt_ml = mu_ml(~cellfun(@isempty, regexp(dw_ml.Properties.RowNames,'Col-?0')));

mu_ml(mu_ml>mu_wt_ml) = mu_wt_ml;

% min/max
mu_ml_min = (log(dw_ml.mean_DW-dw_ml.sd_DW) - log(dw_0)) ./ t_growth;
mu_wt_ml_min = mu_ml_min(~cellfun(@isempty, regexp(dw_ml.Properties.RowNames,'Col-?0')));

mu_ml_max = (log(dw_ml.mean_DW+dw_ml.sd_DW) - log(dw_0)) ./ t_growth;
mu_wt_ml_max = mu_ml_max(~cellfun(@isempty, regexp(dw_ml.Properties.RowNames,'Col-?0')));

% sd
mu_ml_sd = (log(dw_ml.sd_DW) - log(dw_0)) ./ t_growth;

% growth duration
t_growth = das_fl * 24; % [h]

% calculate growth rate
dw_0 = 20 / 1e6;
mu_fl = (log(dw_fl.mean_DW) - log(dw_0)) ./ t_growth;
mu_wt_fl = mu_fl(~cellfun(@isempty, regexp(dw_fl.Properties.RowNames,'Col-?0')));

mu_fl(mu_fl>mu_wt_fl) = mu_wt_fl;

% min/max
mu_fl_min = (log(dw_fl.mean_DW-dw_fl.sd_DW) - log(dw_0)) ./ t_growth;
mu_wt_fl_min = mu_fl_min(~cellfun(@isempty, regexp(dw_fl.Properties.RowNames,'Col-?0')));

mu_fl_max = (log(dw_fl.mean_DW+dw_fl.sd_DW) - log(dw_0)) ./ t_growth;
mu_wt_fl_max = mu_fl_max(~cellfun(@isempty, regexp(dw_fl.Properties.RowNames,'Col-?0')));

% sd
mu_fl_sd = (log(dw_fl.sd_DW) - log(dw_0)) ./ t_growth;

% define genotypes
wt = 'Col-0';
mutants = {'ggt1-1', 'ggt1-2', 'hpr1-1', 'hpr1-2'};

% flux sampling data specifications
light_conditions = {'ml', 'fl'};
ph_ub = [286 191];
timepoints = {'-21'};

% calculate experimental ratios between growth rates
rgr_ratio_exp_mean = {...
    mu_ml(cellfun(@(x)find(ismember(dw_ml.Properties.RowNames, x)), mutants)) / mu_wt_ml,...
   	mu_fl(cellfun(@(x)find(ismember(dw_fl.Properties.RowNames, x)), mutants)) / mu_wt_fl};

sd_fl = mu_fl/mu_wt_fl-mu_fl_min/mu_wt_fl_max;
sd_ml = mu_ml/mu_wt_ml-mu_ml_min/mu_wt_ml_max;
rgr_ratio_exp_sd = {...
    sd_ml(cellfun(@(x)find(ismember(dw_ml.Properties.RowNames, x)), mutants)),...
    sd_fl(cellfun(@(x)find(ismember(dw_ml.Properties.RowNames, x)), mutants))};

% initialize cell arrays for relative growth rates
rgr_ratio_est_mean = cell(1, numel(light_conditions));
rgr_ratio_est_sd = cell(1, numel(light_conditions));

% calculate min/max ratios in rgr between ML and FL
min_ratio = mu_wt_ml_min / mu_wt_fl_max;
max_ratio = mu_wt_ml_max / mu_wt_fl_min;

rgr_wt = zeros(numel(light_conditions), 1);

for l = 1:numel(light_conditions)
    
    % initialize matrix for metabolic phenotype for current light condition
    rgr_ratio_est_mean{l} = nan(numel(mutants), numel(timepoints));
    rgr_ratio_est_sd{l} = nan(numel(mutants), numel(timepoints));
    
    for i = 1:numel(timepoints)
        
        tp = timepoints{i};
        
        file_name = fullfile(res_dir, light_conditions{l},...
            ['flux_samples_', wt, '_' light_conditions{l} '_t_', tp,...
            '_pHUB_' num2str(ph_ub(l)) '.csv']);
        wt_flux = readmatrix(file_name);
        
        wt_rgr = wt_flux(bio_rxn_idx, :);
        wt_rgr(wt_rgr==0) = NaN;
        rgr_wt(l) = mean(wt_rgr, 'omitnan');
        
        for j = 1:numel(mutants)
            
            mut = mutants{j};
            
            file_name = fullfile(res_dir, light_conditions{l},...
                ['flux_samples_', mutants{j}, '_' light_conditions{l} '_t_', tp,...
                '_pHUB_' num2str(ph_ub(l)) '.csv']);
            
            
            if isfile(file_name)
                
                % read flux samples
                mut_flux = readmatrix(file_name);
                
                mut_rgr = mut_flux(bio_rxn_idx, :);
                mut_rgr(mut_rgr==0) = NaN;
                
                rgr_ratio_est_mean{l}(j,i) = mean(mut_rgr ./ wt_rgr, 'omitnan');
                rgr_ratio_est_sd{l}(j,i) = std(mut_rgr ./ wt_rgr, 'omitnan');
            end
        end
    end
end


%% Create scatterplots
colors = lines(numel(mutants)+1);
figure
tiledlayout(1,2)
nexttile
hold on
line([0 1.2], [0 1.2],'color', [.5 .5 .5])
for i=1:numel(mutants)
    h(i) = errorbar(rgr_ratio_exp_mean{1}(i), rgr_ratio_est_mean{1}(i), rgr_ratio_est_sd{1}(i),...
        rgr_ratio_est_sd{1}(i),rgr_ratio_exp_sd{1}(i), rgr_ratio_exp_sd{1}(i),...
        '.', 'Color', colors(i+1,:), 'MarkerSize', 20);
end

xlabel('RGR ratio to Col-0 (experimental)', 'FontSize', 14)
ylabel('RGR ratio to Col-0 (predicted)', 'FontSize', 14)
text(0.05, 0.93, 'a', 'Units', 'normalized', 'FontSize', 14)
set(gca, 'FontSize', 14, 'XLim', [.4 1.1], 'YLim', [.4 1.1], 'Box', 'on',...
    'LineWidth', 1.3)

nexttile
hold on
line([0 1.2], [0 1.2], 'color', [.5 .5 .5])
for i=1:numel(mutants)
    h(i) = errorbar(rgr_ratio_exp_mean{2}(i), rgr_ratio_est_mean{2}(i), rgr_ratio_est_sd{2}(i),...
        rgr_ratio_est_sd{2}(i), rgr_ratio_exp_sd{2}(i), rgr_ratio_exp_sd{2}(i),...
        '.', 'Color', colors(i+1,:), 'MarkerSize', 20);
end
xlabel('RGR ratio to Col-0 (experimental)', 'FontSize', 14)
% legend(h, mutants, 'Location', 'southeast', 'box', 'off')
text(0.05, 0.93, 'b', 'Units', 'normalized', 'FontSize', 14)
set(gca, 'FontSize', 14, 'XLim', [.4 1.1], 'YLim', [.4 1.1], 'Box', 'on',...
    'LineWidth', 1.3)
set(gcf, 'OuterPosition', [134.3333  161.0000  935.3333  475.3333])
exportgraphics(gcf, fullfile(res_dir, 'rgr_check.png'))

fprintf('Ratio between FL and ML conditions (Col-0)\n')
fprintf('Expected\tPredicted\n')
fprintf('%.2f\t%.2f\n', mu_wt_fl / mu_wt_ml, rgr_wt(2) / rgr_wt(1))
