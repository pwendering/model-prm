% calculate relative growth rates
clear;clc

data_dir = fullfile('..', 'Data');
res_dir = fullfile('..', 'Results');

% define input file and get sheet names of excel file (i.e. light
% conditions)
dw_file_name = fullfile(data_dir, '221116_AllWeight_PRmutants.xlsx');
conditions = sheetnames(dw_file_name);

% output file name
out_filename = fullfile(res_dir, 'rgr_all.xlsx');
LETTERS = upper('abcdefghijklmnopqrstuvwxyz');

% read seed weights from file
seed_weights = readtable(fullfile(data_dir, 'seed_weights.txt'));
seed_weights.weight_g = seed_weights.weight_ug / 1e6;

% plot colors
color_tab = array2table([[89 90 92]; [75 82 160]; [183 212 235]; [37 168 220];
    [241 123 46]; [192 35 39]; [166 208 139]; [13 148 69]; [113 58 151]]/255);
color_tab.genotype = {'Col-0', 'ggt2', 'ggt1-1', 'ggt1-2', 'hpr1-1', 'hpr1-2',...
    'gldt1-1', 'PsL', 'PGLP'}';

% initialize growth rates and genotype names
growth_rates_mean = zeros(10, numel(conditions));
growth_rates_median = zeros(10, numel(conditions));
growth_rates_p25 = zeros(10, numel(conditions));
growth_rates_p75 = zeros(10, numel(conditions));
growth_rates_wmin = zeros(10, numel(conditions));
growth_rates_wmax = zeros(10, numel(conditions));
growth_rates_sd = zeros(10, numel(conditions));
outliers = cell(1, numel(conditions));

for i = 1:numel(conditions)
    fprintf('Condition %s\n', conditions{i})
    
    % find columns with info, genotype, and dry weights
    dw_tab = readtable(dw_file_name,...
        'Sheet', conditions{i});
    dw_cols = 3:4:size(dw_tab, 2);
    geno_cols = 2:4:size(dw_tab, 2);
    info_cols = 1:4:size(dw_tab, 2);
    info = dw_tab{3, info_cols};
    exp_ids = dw_tab.Properties.VariableNames(1:4:size(dw_tab, 2));
    age_d = cellfun(@(x)str2double(regexp(x, '\d+', 'match')), info);
    
    % loop over experiments
    for j = 1:numel(exp_ids)
        
        age_h = age_d(j) * 24;
        dw_mg = dw_tab{:, dw_cols(j)};
        dw_g = dw_mg / 1000;
        
        tmp_genos = strrep(dw_tab{:, geno_cols(j)}, 'Col0', 'Col-0');
        tmp_genos = strrep(tmp_genos, 'PSL', 'PsL');
        tmp_genos = strtok(tmp_genos, ' ');
        tmp_genos = tmp_genos(~cellfun('isempty', tmp_genos));
        dw_g = dw_g(~cellfun('isempty', tmp_genos));
        
        tab = cell2table(...
            [[exp_ids(j) {'Genotype', 'RGR'}];...
            [repmat({''}, size(mu)) tmp_genos num2cell(mu)]]);
        
        col_idx = (j-1)*size(tab,2)+1:j*size(tab,2);
        writetable(tab, out_filename, 'Sheet', conditions{i},...
            'Range', [LETTERS(col_idx(1)) ':' LETTERS(col_idx(end))],...
            'WriteVariableNames', false)
    end
end