%% Metabolite data
filename = fullfile('..','Data','220404_AllData_Photorespiration_TB.xlsx');

% ML control
sheet = 'Metabolites ML control';
[met_table_ml_ctl,av_per_genotype_ml_ctl,genotype_ml_ctl,time_in_shift_ml_ctl] = parseMetaboliteData(filename,sheet);

% FL control
sheet = 'Metabolites FL control';
[met_table_fl_ctl,av_per_genotype_fl_ctl,genotype_fl_ctl,time_in_shift_fl_ctl] = parseMetaboliteData(filename,sheet);

% ML to FL
sheet = 'Metabolites ML to FL shift';
[met_table_ml_fl,av_per_genotype_ml_fl,genotype_ml_fl,time_in_shift_ml_fl] = parseMetaboliteData(filename,sheet);

% FL to ML
sheet = 'Metabolites FL to ML shift';
[met_table_fl_ml,av_per_genotype_fl_ml,genotype_fl_ml,time_in_shift_fl_ml] = parseMetaboliteData(filename,sheet);

clear sheet


function [met_table,av_per_genotype,genotype,shifttime] = parseMetaboliteData(filename,sheet)

% read table from file
met_table = readtable(filename,'Sheet',sheet);

% sort by genotype and time point
[~,idx_sort_geno] = sort(met_table.Genotype);
met_table = met_table(idx_sort_geno,:);

% create unique strings for genotype and timeshift
geno_time = strcat(...
    met_table.Genotype,...
    num2str(met_table.TimeInLightShift_h_)...
    );
[uniq_geno_time,ia] = unique(geno_time,'stable');
genotype = met_table.Genotype(ia);
shifttime = met_table.TimeInLightShift_h_(ia);

% take average over replicates (genotype and timeshift)
av_per_genotype = cell2mat(...
    cellfun(@(x)...
        mean(met_table{ismember(geno_time,x),4:end},1,'omitnan'),...
        uniq_geno_time,...
        'un',0 ...
        )...
    );
end
