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
        
%         tmp_genos_uniq = unique(tmp_genos);
%         growth_rates_mean = nan(numel(tmp_genos_uniq), 1);
%         growth_rates_sd = nan(numel(tmp_genos_uniq), 1);
        
        sw = cellfun(@(x)seed_weights.weight_g(ismember(seed_weights.genotype,x)), tmp_genos);
        mu = (log(dw_g) - log(sw)) ./ age_h;
        
%         for k = 1:numel(tmp_genos_uniq)
%             dw = dw_g(ismember(tmp_genos, tmp_genos_uniq(k)));
%             sw = seed_weights.weight_g(ismember(seed_weights.genotype, tmp_genos_uniq(k)));
%             if isempty(sw)
%                 sw = mean(seed_weights.weight_g);
%             end
%             mu = (log(dw) - log(sw)) ./ age_h;
%             
%             growth_rates_mean(k) = mean(mu, 'omitnan');
%             growth_rates_sd(k) = std(mu, 'omitnan');
%         end
        
        %         tab = cell2table(...
        %             [[exp_ids(j) {'Genotype', 'RGR_mean', 'RGR_sd'}];...
        %             [repmat({''}, size(growth_rates_mean)) tmp_genos_uniq ...
        %             num2cell(growth_rates_mean) num2cell(growth_rates_sd)]]);
        
        tab = cell2table(...
            [[exp_ids(j) {'Genotype', 'RGR'}];...
            [repmat({''}, size(mu)) tmp_genos num2cell(mu)]]);
        
        col_idx = (j-1)*size(tab,2)+1:j*size(tab,2);
        writetable(tab, out_filename, 'Sheet', conditions{i},...
            'Range', [LETTERS(col_idx(1)) ':' LETTERS(col_idx(end))],...
            'WriteVariableNames', false)
    end
end

%{
rem_idx = ~any(growth_rates_mean, 2);
growth_rates_mean(rem_idx, :) = [];
growth_rates_sd(rem_idx, :) = [];
growth_rates_median(rem_idx, :) = [];
growth_rates_p25(rem_idx, :) = [];
growth_rates_p75(rem_idx, :) = [];
growth_rates_wmin(rem_idx, :) = [];
growth_rates_wmax(rem_idx, :) = [];

tab = array2table(growth_rates_mean,...
    'VariableNames', conditions,...
    'RowNames', genotypes);

%% Boxplots
% re-order genotypes according to order in color table
new_order = cellfun(@(x)find(ismember(genotypes, x)), color_tab.genotype);
growth_rates_mean = growth_rates_mean(new_order, :);
growth_rates_sd = growth_rates_sd(new_order, :);
growth_rates_median = growth_rates_median(new_order, :);
growth_rates_p25 = growth_rates_p25(new_order, :);
growth_rates_p75 = growth_rates_p75(new_order, :);
growth_rates_wmin = growth_rates_wmin(new_order, :);
growth_rates_wmax = growth_rates_wmax(new_order, :);
genotypes = genotypes(new_order);

for i = 1:numel(conditions)
    outliers{i} = outliers{i}(new_order);
end


bp = bar(growth_rates_mean');
xoffset = zeros(size(genotypes));
for i = 1:numel(genotypes)
    xoffset(i) = bp(i).XOffset;
end
close

figure
hold on
for i = 1:numel(conditions)
    h = zeros(numel(genotypes), 1);
    for j = 1:numel(genotypes)
        ol = outliers{i}{j};
        if ~isempty(ol)
            scatter(repelem(i+xoffset(j), 1, numel(ol)), ol,...
                5, 'filled', 'MarkerFaceColor', [.4 .4 .4],...
                'MarkerEdgeColor', [.4 .4 .4])
        end
        line([i+xoffset(j) i+xoffset(j)],...
            [growth_rates_wmin(j,i) growth_rates_wmax(j,i)],...
            'color', 'k',...
            'linewidth', .8)
        h(j) = line(...
            [i+xoffset(j) i+xoffset(j)],...
            [growth_rates_p25(j,i) growth_rates_p75(j,i)],...
            'LineWidth', 8,...
            'Color', color_tab{ismember(color_tab.genotype, genotypes(j)),1:3});
        scatter(i+xoffset(j), growth_rates_median(j,i),...
            'Marker', '_',...
            'LineWidth', 1,...
            'MarkerFaceColor', 'k',...
            'MarkerEdgeColor', 'k')
    end
end
hold off

xticks(1:numel(conditions))
xticklabels(conditions)
ylabel('Relative growth rate [h^{-1}]')
genotypes(~ismember(genotypes, 'Col-0')) = cellfun(@(x)['{\it' x '}'],...
    genotypes(~ismember(genotypes, 'Col-0')), 'un', 0);
legend(h, genotypes,...
    'Location', 'northwest',...
    'NumColumns', 5,...
    'Box', 'off')
xlabel('Light condition')
set(gca, 'FontSize', 14)
set(gcf, 'OuterPosition', 1000*[0.0190    0.1770    1.1973    0.5340])
print(fullfile('..', 'Results', 'rgr_boxplot.png'), '-dpng', '-r300')
%}