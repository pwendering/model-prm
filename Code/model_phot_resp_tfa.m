%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Modeling of photorespiratory mutants using TFA
% using the matTFA toolbox
%   Reference: Salvy et al. (2019), doi: 10.1093/bioinformatics/bty499
%
% Github Fork: https://github.com/pwendering/matTFA
%
% Philipp Wendering, University of Potsdam (philipp.wendering@gmail.com)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc

% number of workers for parallel computing

N_CPU = 1;%20;
N_SAMPLES = 1000;

% define data and output directory
data_dir = '../Data';
res_dir = '../Results';

% set up parallel pool
PAR_FLAG = 0;
if N_CPU > 1
    [status,ME]=setupParpool(N_CPU);
    if status
        error('Error: %s\nParallel pool could not be initialized, proceeding without parallelization',ME.message)
    else
        PAR_FLAG = 1;
    end
end

% initialize COBRA
initCobraToolbox(false)

% fix setfield error
curr_path = pwd;
cd ~/cobratoolbox/external/analysis/PolytopeSamplerMatlab/code/utils
system('mv setfield.m Setfield.m');
cd(curr_path)

% CPLEX path and COBRA solver
cplexPath = '~/bin/ibm/ILOG/CPLEX_Studio129/cplex/matlab/x86-64_linux/';
% cplexPath = fullfile('C:/Program Files/IBM/ILOG/CPLEX_Studio129/cplex/matlab/x64_win64/');
addpath(genpath(cplexPath))
changeCobraSolver('ibm_cplex','all');

% path to matTFA toolbox
mattfa_path = fullfile('..','..','matTFA');
addpath(genpath(fullfile(mattfa_path,'matTFA')))
addpath(genpath(fullfile(mattfa_path,'thermoDatabases')))
addpath(genpath(fullfile(mattfa_path,'models')))
addpath(genpath(fullfile(mattfa_path,'plotting')))

%% Load the model
% metabolic model (updated metabolite names and added inchi keys)
% model reference: Arnold and Nikoloski (2014), doi: 10.1104/pp.114.235358
load(fullfile(data_dir,'AraCore-updated-rev.mat'));
model.metCompSymbol = regexprep(model.mets,'.*\[(?<comp>\w)\]$','$<comp>');
tmp = load(fullfile(mattfa_path,'models','smallEcoli.mat'));
eco_model = tmp.smallEcoli;

% change compartment symbol for Peroxisome from p to x
model.comps(strcmp(model.comps,'p')) = {'x'};
model.mets = strrep(model.mets,'[p]','[x]');
model.rxns = regexprep(model.rxns,'_p$','_x');

% update compartment data
model.CompartmentData = eco_model.CompartmentData;

% add chloroplast
model.CompartmentData.compSymbolList{end+1} = 'h';
model.CompartmentData.compMaxConc(end+1) = 0.05;
model.CompartmentData.compMinConc(end+1) = 1e-6;
model.CompartmentData.compNameList{end+1} = 'Chloroplast';
model.CompartmentData.ionicStr(end+1) = 0;
model.CompartmentData.membranePot(:,end+1) = zeros(size(model.CompartmentData.membranePot,1),1);
model.CompartmentData.membranePot(end+1,:) = zeros(1,size(model.CompartmentData.membranePot,2));
model.CompartmentData.pH(end+1) = 8; % https://www.pnas.org/doi/10.1073/pnas.0403709101

% add lumen
model.CompartmentData.compSymbolList{end+1} = 'l';
model.CompartmentData.compMaxConc(end+1) = 0.05;
model.CompartmentData.compMinConc(end+1) = 1e-6;
model.CompartmentData.compNameList{end+1} = 'Lumen';
model.CompartmentData.ionicStr(end+1) = 0;
model.CompartmentData.membranePot(:,end+1) = zeros(size(model.CompartmentData.membranePot,1),1);
model.CompartmentData.membranePot(end+1,:) = zeros(1,size(model.CompartmentData.membranePot,2));
model.CompartmentData.pH(end+1) = 5; % http://www.esalq.usp.br/lepse/imgs/conteudo_thumb/What-is-the-pH-of-the-lumen-of-each-intracellular-organelle.pdf

% add Intermembrane space
model.CompartmentData.compSymbolList{end+1} = 'i';
model.CompartmentData.compMaxConc(end+1) = 0.05;
model.CompartmentData.compMinConc(end+1) = 1e-6;
model.CompartmentData.compNameList{end+1} = 'Intermembrane Space';
model.CompartmentData.ionicStr(end+1) = 0;
model.CompartmentData.membranePot(:,end+1) = zeros(size(model.CompartmentData.membranePot,1),1);
model.CompartmentData.membranePot(end+1,:) = zeros(1,size(model.CompartmentData.membranePot,2));
model.CompartmentData.pH(end+1) = 7; % https://www.sciencedirect.com/science/article/pii/B9780123821638000219

% update compartment pH
model.CompartmentData.pH(2) = 8; % https://www.sciencedirect.com/science/article/pii/B9780123821638000219
model.CompartmentData.pH(strcmp(model.CompartmentData.compSymbolList,'c')) = 7.3; % https://www.sciencedirect.com/science/article/pii/S1674205214602173
model.CompartmentData.pH(strcmp(model.CompartmentData.compSymbolList,'n')) = 7.2; % https://www.sciencedirect.com/science/article/pii/S1674205214602173
model.CompartmentData.pH(strcmp(model.CompartmentData.compSymbolList,'x')) = 8.4; % https://www.sciencedirect.com/science/article/pii/S1674205214602173

% increase default metabolite concentration range
model.CompartmentData.compMinConc(:) = 1e-8;
model.CompartmentData.compMaxConc(:) = 0.05;

% convert to reverible and add rev field
model.rev = model.lb<0 & model.ub>0;

% Limit the bounds of the fluxes that are higher than 100 or lower than
% -100 mmol/(gDW * h)
if any(model.lb<-100) || any(model.ub>100)
    model.lb(model.lb<-100) = -100;
    model.ub(model.ub>+100) = +100;
end

%% Load the thermodynamics database
tmp = load(fullfile(data_dir, 'thermo_data.mat'));
ReactionDB = tmp.DB_AlbertyUpdate;
clear tmp

%% translate InChI to SEED ID
% translate InchI-Key to SEED ID
SEEDID = repmat({''},numel(model.mets),1);
for i=1:numel(model.mets)
    tmp_inchi_key = strsplit(model.metisinchikeyID{i},'|');
    match_idx = 0;
    j = 0;
    while j<numel(tmp_inchi_key) && sum(match_idx)==0
        j = j + 1;
        match_idx = startsWith(ReactionDB.compound.InChIKey,regexprep(tmp_inchi_key{j},'-\w$',''));
        if sum(match_idx)==1
            SEEDID{i} = ReactionDB.compound.ID{match_idx};
        elseif sum(match_idx) > 1
            fprintf('%s: multiple matches:\n',model.metNames{i})
            cellfun(@(x)fprintf('\t%s\n',x),ReactionDB.compound.ID(match_idx));
            SEEDID{i} = ReactionDB.compound.ID{find(match_idx,1)};
        end
    end
end
model.metSEEDID = SEEDID;

unmatched_met_translation = [
    {'1-Pyrroline-5-carboxylate'               } {'cpd02431'}
    {'10-Formyltetrahydrofolate'               } {'cpd00201'}
    {'2,3-Dihydroxy-3-methylbutanoate'         } {'cpd19030'}
    {'2,3-Dihydroxy-3-methylpentanoate'        } {'cpd19046'}
    {'2-Deoxycytidine 5-triphosphate'          } {'cpd00356'}
    {'5,10-Methylenetetrahydrofolate'          } {'cpd00125'}
    {'5,6,7,8-Tetrahydrofolate'                } {'cpd00087'}
    {'5-Methyltetrahydrofolate'                } {'cpd00345'}
    {'5-Phospho-D-ribosyl-N-formylglycineamide'} {'cpd02678'}
    {'5-Phospho-D-ribosylglycineamide'         } {'cpd02394'}
    {'Acetylglutamate'                         } {'cpd00477'}
    {'Acyl-carrier protein'                    } {'cpd11493'}
    {'Aminomethyldihydrolipoylprotein'         } {'cpd11830'}
    {'Argininosuccinate'                       } {'cpd02152'}
    {'Cellulose, 1 repeat units'               } {'cpd11746'}
    {'Cellulose, 2 repeat units'               } {'cpd11746'}
    {'Cellulose, 3 repeat units'               } {'cpd11746'}
    {'Coenzyme A'                              } {'cpd00010'}
    {'D-Fructose 1,6-bisphosphate'             } {'cpd00290'}
    {'D-Glucose 1-phosphate'                   } {'cpd00501'}
    {'D-Ribose'                                } {'cpd00105'}
    {'D-ribulose 1,5-bisphosphate'             } {'cpd00871'}
    {'Diaminopimelate'                         } {'cpd00504'}
    {'Dihydrolipolprotein'                     } {'cpd12225'}
    {'Dihydroorotate'                          } {'cpd00282'}
    {'Erythrose 4-phosphate'                   } {'cpd00236'}
    {'Glyceraldehyde 3-phosphate'              } {'cpd00102'}
    {'H+, proton'                              } {'cpd00067'}
    {'Lipoylprotein'                           } {'cpd12005'}
    {'Malonyl-Coenzyme A'                      } {'cpd00070'}
    {'Malonyl-acyl-carrier-protein'            } {'cpd11492'}
    {'N6-(1,2-Dicarboxyethyl)-AMP'             } {'cpd02375'}
    {'Oxidized cytochrome c'                   } {'cpd27754'}
    {'Oxidized ferredoxin'                     } {'cpd11621'}
    {'Oxidized plastocyanin'                   } {'cpd27746'}
    {'Oxidized plastoquinone'                  } {'cpd30517'}
    {'Oxidized thioredoxin'                    } {'cpd27735'}
    {'Photon'                                  } {'cpd30510'}
    {'Prephenate'                              } {'cpd00219'}
    {'Reduced cytochrome c'                    } {'cpd28079'}
    {'Reduced ferredoxin'                      } {'cpd28082'}
    {'Reduced plastocyanin'                    } {'cpd12239'}
    {'Reduced plastoquinone'                   } {'cpd30519'}
    {'Reduced thioredoxin'                     } {'cpd28060'}
    {'Sedoheptulose 1,7-bisphosphate'          } {'cpd00349'}
    {'Sedoheptulose 7-phosphate'               } {'cpd00238'}
    {'Starch, X+1 glucose units'               } {'cpd11657'}
    {'Starch, X+2 glucose units'               } {'cpd11657'}
    {'Starch, X+3 glucose units'               } {'cpd11657'}
    {'Starch, X+5 glucose units'               } {'cpd11657'}
    {'Succinyl-CoA'                            } {'cpd00078'}
    {'Succinyldihydrolipoamide'                } {'cpd00860'}
    {'UDP-Glucose'                             } {'cpd00026'}
    {'Ubiquinol'                               } {'cpd25914'}
    ];

for i=1:size(unmatched_met_translation,1)
    model.metSEEDID(ismember(model.metNames,unmatched_met_translation{i,1})) = ...
        unmatched_met_translation(i,2);
end

%% photon uptake
fprintf('Setting photon uptake to 200 umol/m2/s\n')
LMA = 22.2; % [gDW/m2] (Hummel et al. 2010, doi: 10.1104/pp.110.157008)
I_ML = 3600 * 200 / LMA / 1000; % [mmol/gDW/h]
I_FL = 3600 * 700 / LMA / 1000; % [mmol/gDW/h] (use kinetic model to integrate fluctuations?)
% I_ML = 100*200/700; % scale by assuming the 700 uE is saturating
% I_FL = 100;

model.ub(findRxnIDs(model,'Im_hnu')) = I_ML;

%% Set objective to light-limiting biomass reaction
model.c(:) = 0;
model.c(findRxnIDs(model,'Bio_opt')) = 1;
model.ub(findRxnIDs(model,{'Bio_CLim','Bio_NLim'})) = 0;

%% Test simple fba
solFBA = optimizeCbModel(model);
% We can set a lower bound for growth (e.g. 50% of maximal growth)
min_obj = roundsd(0.4*solFBA.f, 2, 'floor');
model.lb(model.c==1) = min_obj;

%% Perform FVA at 90% of the optimum
old_ub = model.lb(model.c==1);
model.lb(model.c==1) = 0.9*solFBA.f;
fva_wt = runMinMax(model,model.rxns,'runParallel',PAR_FLAG);
model.lb(model.c==1) = old_ub;
% - logical vector with indices for bidirectional reactions (flux ranges crossing zero)
n = @(x) x(:,1)<-1e-9 & x(:,2)>1e-9;
is_bd_fva_wt = (n(fva_wt));

% Are there any blocked reactions?
% solver tolerance is 1e-9
SolTol = 1e-9;
id_Blocked_in_FBA = find( (fva_wt(:,1)>-SolTol & fva_wt(:,1)<SolTol) & ...
    (fva_wt(:,2)>-SolTol & fva_wt(:,2)<SolTol) );
% If there exist block reactions
while ~isempty(id_Blocked_in_FBA)
    % remove them
    model = removeRxns(model, model.rxns(id_Blocked_in_FBA));
    fva_wt = runMinMax(model,'runParallel',PAR_FLAG);
    id_Blocked_in_FBA = find( (fva_wt(:,1)>-SolTol & fva_wt(:,1)<SolTol) & ...
        (fva_wt(:,2)>-SolTol & fva_wt(:,2)<SolTol) );
end

%% Prepare for TFA
%need field for description
prepped_m = prepModelforTFA(model, ReactionDB, model.CompartmentData);

%% Convert to TFA
tmp = convToTFA(prepped_m, ReactionDB, [], 'DGo', [], min_obj);

% Add net flux variables, which are equal to forwards flux - backwards flux
% NF_rxn = F_rxn - B_rxn
this_tmodel = addNetFluxVariables(tmp);
NF_idx = getAllVar(this_tmodel,{'NF'});

%% Perform FVA at 90% of the optimum
old_ub = model.lb(model.c==1);
model.lb(model.c==1) = 0.9*solFBA.f;
fva_wt = runMinMax(model,model.rxns,'runParallel',PAR_FLAG);
model.lb(model.c==1) = old_ub;

%% Read experimental data
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

% dry weights
dw_ml = readtable(...
    fullfile(data_dir, data_files{contains(data_files, 'dw_ml_34')}),...
    'ReadRowNames', true...
    );
das_ml = 34;
dw_fl = readtable(...
    fullfile(data_dir, data_files{contains(data_files, 'dw_fl_43')}),...
    'ReadRowNames', true...
    );
das_fl = 43;


% net CO2 assimilation rates
A_ml = readtable(...
    fullfile(data_dir, data_files{contains(data_files, 'A_ml')})...
    );

A_fl = readtable(...
    fullfile(data_dir, data_files{contains(data_files, 'A_fl')})...
    );

A_max = max([A_fl.mean_A; A_ml.mean_A]);


% growth duration
t_growth = das_ml * 24; % [h]

% calculate growth rate
dw_0 = 20 / 1e6;
mu_ml = (log(dw_ml.mean_DW) - log(dw_0)) ./ t_growth;
mu_wt_ml = mu_ml(~cellfun(@isempty, regexp(dw.Properties.RowNames,'Col?-0')));

mu_ml(mu_ml>mu_wt_ml) = mu_wt_ml;


% growth duration
t_growth = das_fl * 24; % [h]

% calculate growth rate
dw_0 = 20 / 1e6;
mu_fl = (log(dw_fl.mean_DW) - log(dw_0)) ./ t_growth;
mu_wt_fl = mu_fl(~cellfun(@isempty, regexp(dw.Properties.RowNames,'Col?-0')));

mu_fl(mu_fl>mu_wt_fl) = mu_wt_fl;

fprintf('Ratio between growth rates in ML and FL conditions: %.3f\n', mu_wt_ml / mu_wt_ml)

% define mutant names and reaction knock-outs
mutants = {...
    'ggt1-1'...
    'ggt1-2'...
    'ggt2'...
    'hpr1-1'...
    'hpr1-2'...
    };

ko_rxns = {...
    {'GGAT_p','GGAT_h','AlaTA_h','AlaTA_p'}...
    {'GGAT_p','GGAT_h','AlaTA_h','AlaTA_p'}...
    {'GGAT_p','GGAT_h','AlaTA_h','AlaTA_p'}...
    {'GCEADH_p','GCEADH_h'}...
    {'GCEADH_p','GCEADH_h'}...
    };

% light conditions and photon uptake upper bounds
% ML must be the first condition because ML reference growth must be
% predicted fist 
l_cond = {'ml','fl','ml_fl','fl_ml'};

ph_ub = I_ML;

% time points
tp = [-21 3 27 75 123];

%% dry weight to volume conversion
% weight increases linearly with volume
% Chen2016: https://doi.org/10.3389/fpls.2016.00392 (lettuce)
% Huxley1971: https://doi.org/10.2307/2402133 (several species)
% Aoyagi1992: https://doi.org/10.1016/0922-338X(92)90144-J

fw_per_dw = 10; % [gFW/gDW]
% fw_per_ml = 29.26 / 1000; % [gFW/ml]; Fig 8 (Chen2016)
fw_per_ml = 1 / (1.192-0.0380); % [gFW/ml]; Tab 1 (Huxley1971) mean a + mean b
% fw_per_ml = mean([1.010 1.028]); % [gFW/ml]; Abstract (Aoyagi1992)
fw_per_l = fw_per_ml * 1000; % [gFW/l]
dw_per_l = fw_per_l / fw_per_dw; % [gDW/l]

% Winter et al. 1993: total leaf volume: 902 uL / mg Chlorophyll (Barley)
% Anika: 0.8 mg Chl. per gFW
% fw_per_dw = 10; % [gFW/gDW]
% chl_per_fw = 0.8e-3; % [gChl / gFW]
% vol_per_chl = 0.902; % [mL / mgChl] = [L / gChl]
% fw_per_l = 1/(vol_per_chl * chl_per_fw);
% dw_per_l = fw_per_l / fw_per_dw;

%% loop over light conditions
for lc_idx = 1:numel(l_cond)
    cmdsz = matlab.desktop.commandwindow.size;
    fprintf([repmat('-',1,cmdsz(1)) '\n'])
    fprintf('Light Condition: %s\n', l_cond{lc_idx})

    %% loop over timepoints
    for t_idx = 1:numel(tp)
        cmdsz = matlab.desktop.commandwindow.size;
        fprintf([repmat('-',1,cmdsz(1)) '\n'])
        fprintf('Time after shift: %d\n', tp(t_idx))
        
        %% get metabolite concentration ranges for the wild type
        met_av = readtable(fullfile(data_dir, ['met_' l_cond{lc_idx} '.csv']));
        met_sd = readtable(fullfile(data_dir, ['met_' l_cond{lc_idx} '_sd.csv']));
        meas_names = met_av.Properties.VariableNames(3:end);
        
        % match names in file to model metabolite IDs
        meas_ids = name2id(meas_names);
        nm_idx = cellfun(@isempty,meas_ids);
        meas_ids = meas_ids(~nm_idx);
        
        met_av_wt = met_av{strcmp(met_av.genotype,'Col-0') & met_av.time==tp(t_idx),3:end};
        met_sd_wt = met_sd{strcmp(met_av.genotype,'Col-0') & met_av.time==tp(t_idx),3:end};
        
        % remove data for unmatched metabolites
        met_av_wt = met_av_wt(~nm_idx);
        met_sd_wt = met_sd_wt(~nm_idx);
        
        meas_mw = zeros(size(meas_ids));
        for i=1:numel(meas_ids)
            first_idx = find(contains(model.mets,[meas_ids{i} '[']),1);
            meas_mw(i) = getMolecularMass(model.metFormulas(first_idx)); % [mol/gDW]
        end
        
        meas_conc_min_wt = met_av_wt - met_sd_wt;
        meas_conc_max_wt = met_av_wt + met_sd_wt;
        
        meas_conc_min_wt = meas_conc_min_wt' * dw_per_l ./ meas_mw; % [M]
        meas_conc_max_wt = meas_conc_max_wt' * dw_per_l ./ meas_mw; % [M]
        
        meas_conc_min_wt(meas_conc_min_wt<0) = min(model.CompartmentData.compMinConc);
        
        % match measured concentrations to LC variables in TFA model
        varnames = strcat(...
            'LC_',...
            regexprep(...
            model.mets(startsWith(model.mets,strcat(meas_ids(~cellfun(@isempty,meas_ids)),'['))),...
            '(.*)\[(\w)\]$',...
            '$1_$2'...
            )...
            );
        
        LC_varNames = varnames(ismember(varnames,this_tmodel.varNames));
        C_lb_wt = zeros(size(LC_varNames));
        C_ub_wt = zeros(size(LC_varNames));
        
        for i=1:numel(meas_ids)
            meas_id_match = ismember(regexprep(LC_varNames,'LC_(.*)_\w','$1'),meas_ids(i));
            if sum(meas_id_match)>0
                
                C_lb_wt(meas_id_match) = meas_conc_min_wt(i);
                C_ub_wt(meas_id_match) = meas_conc_max_wt(i);
            else
                fprintf('%s not found in model variables\n',meas_ids{i})
            end
        end
        
        
        for m_idx = 1:numel(mutants)
            
            cmdsz = matlab.desktop.commandwindow.size;
            fprintf([repmat('-',1,floor(cmdsz(1)/2)) '\n'])
            fprintf('>>> Processing genotype %s ...\n', mutants{m_idx})
            
            % determine ratio of growth rates (wt/mutant)
            biomass_ratio = mu_wt / ...
                mu(strcmp(dw.Properties.RowNames,mutants{m_idx}));
            
            %% get metabolite concentration ranges for the mutant
            met_av_mut = met_av{strcmp(met_av.genotype,mutants{m_idx})&met_av.time==tp(t_idx),3:end};
            met_sd_mut = met_sd{strcmp(met_av.genotype,mutants{m_idx})&met_av.time==tp(t_idx),3:end};
            
            met_av_mut = met_av_mut(~nm_idx);
            met_sd_mut = met_sd_mut(~nm_idx);
            
            meas_conc_min_mut = met_av_mut - met_sd_mut;
            meas_conc_max_mut = met_av_mut + met_sd_mut;
            
            meas_conc_min_mut = meas_conc_min_mut' * dw_per_l ./ meas_mw; % [M]
            meas_conc_max_mut = meas_conc_max_mut' * dw_per_l ./ meas_mw; % [M]
            
            meas_conc_min_mut(meas_conc_min_mut<0) = min(model.CompartmentData.compMinConc);
            
            C_lb_mut = zeros(size(LC_varNames));
            C_ub_mut = zeros(size(LC_varNames));
            
            for i=1:numel(meas_ids)
                meas_id_match = ismember(regexprep(LC_varNames,'LC_(.*)_\w','$1'),meas_ids(i));
                if sum(meas_id_match)>0
                    
                    C_lb_mut(meas_id_match) = meas_conc_min_mut(i);
                    C_ub_mut(meas_id_match) = meas_conc_max_mut(i);
                    
                else
                    fprintf('%s not found in model variables\n',meas_ids{i})
                end
            end
            
            % find the indices of these variables in the variable names of the tfa
            id_LC_varNames = find_cell(LC_varNames, this_tmodel.varNames);
            
            for l_idx = 1:numel(ph_ub)
                
                fprintf('[[ Photon uptake = %.2g mmol/gDW/h ]]\n',ph_ub(l_idx))
                model.ub(findRxnIDs(model,'Im_hnu')) = ph_ub(l_idx);
                this_tmodel.var_ub(strcmp(this_tmodel.varNames,'F_Im_hnu')) = ph_ub(l_idx);
                
                % test if growth is possible with knock-outs
                ko_model = removeRxns(model,ko_rxns{m_idx});
                g_ko = optimizeCbModel(ko_model).f;
                g_wt = optimizeCbModel(model).f;
                gt_bool = g_ko > 1e-6;
                if gt_bool
                    fprintf('Growth test passed\n')
                    fprintf('FBA results: WT: %.4e, %s: %.4e\n', g_ko, mutants{m_idx}, g_wt)
                else
                    fprintf('KO disables growth prediction - skipping\n')
                    continue
                end
                
                %% Solve TFA for wild type model
                wt_tmodel = this_tmodel;
                
                % if calculated growth rate for FL applies, constrain ratio
                % of growth rate between ML and FL
                if contains(l_cond(lc_idx), 'fl') && ~contains(l_cond(lc_idx), 'fl_')
                    wt_tmodel.var_ub(wt_tmodel.f==1) = ref_ml_growth * mu_wt_fl / mu_wt_ml + 1e-10;
                    wt_tmodel.var_lb(wt_tmodel.f==1) = ref_ml_growth * mu_wt_fl / mu_wt_ml - 1e-10;
                end
                
                % add relaxed metabolite concentration ranges to measured
                % metabolites
                wt_tmodel = addRelaxedMetConcRanges(wt_tmodel,LC_varNames,log(C_lb_wt),log(C_ub_wt));
                
                % add minimization objective for relaxation variables
                bio_obj = wt_tmodel.f;
                wt_tmodel.f(startsWith(wt_tmodel.varNames,'EPS_')) = -1;
                
                % solve TFA
                tfa_wt = solveTFAmodelCplex(wt_tmodel);
                
                % reset model objective
                wt_tmodel.f = bio_obj;
                
                if isempty(tfa_wt.x)
                    fprintf('WT thermo model cannot be solved with relaxed metabolite concentration ranges\n')
                    continue
                end
                
                % positive and negative concentration relaxations
                eps_minus_wt = tfa_wt.x(startsWith(wt_tmodel.varNames,'EPS_MINUS'));
                eps_plus_wt = tfa_wt.x(startsWith(wt_tmodel.varNames,'EPS_PLUS'));
                ch_conc_idx = eps_minus_wt>1e-4 | eps_plus_wt>1e-4;
                
                % (relaxed) TFA metabolite concentrations
                met_conc_wt = tfa_wt.x(cellfun(@(x)find(ismember(wt_tmodel.varNames,x)),...
                    LC_varNames));
                
                % create figure with updated metabolite concentration ranges
                tmp_fig = figure(...
                    'Visible','off',...
                    'units','normalized',...
                    'outerposition',[0 0 1 1]);
                tiledlayout(2,1)
                nexttile
                labels = categorical(strrep(erase(LC_varNames,'LC_'),'_','\_'));
                arrayfun(@(i)line([labels(i) labels(i)],[C_lb_wt(i) C_ub_wt(i)],'linewidth',4,'color','k'),1:numel(LC_varNames))
                hold on
                h = zeros(2,1);
                h(1)=scatter(labels(~ch_conc_idx),exp(met_conc_wt(~ch_conc_idx)),15,'filled');
                h(2)=scatter(labels(ch_conc_idx),exp(met_conc_wt(ch_conc_idx)),15,'filled');
                hold off
                set(gca,'YScale','log','FontSize',10)
                legend(h,{'unchanged','relaxed'},'FontSize',14,'box','off',...
                    'location','southeast')
                ylabel('TFA metabolite concentration [M]','FontSize',14)
                text(0.01,0.98,'Col-0','units','normalized','fontweight','bold')
                
                % fix metabolite measured metabolite concentrations to newly
                % obtained values with 10% tolerance
                wt_tmodel.var_lb(cellfun(@(x)find(ismember(wt_tmodel.varNames,x)),...
                    LC_varNames)) = met_conc_wt-abs(.1*(met_conc_wt));
                wt_tmodel.var_ub(cellfun(@(x)find(ismember(wt_tmodel.varNames,x)),...
                    LC_varNames)) = met_conc_wt+abs(.1*(met_conc_wt));
                
                % WT optimal growth rate
                tfa_wt_2 = solveTFAmodelCplex(wt_tmodel);
                if ~isempty(tfa_wt_2.x)
                    wt_opt = tfa_wt_2.x(wt_tmodel.f==1);
                else
                    fprintf('Growth optimization with fixed metabolite concentrations not successful\n')
                    wt_opt = tfa_wt.x(wt_tmodel.f==1);
                end
                
                % save optimal growth rate for ML condition
                if lc_idx==1
                    ref_ml_growth = wt_opt;
                end
                
                % set lower bound for biomass to 90% of optimal value
                wt_tmodel.var_lb(wt_tmodel.f==1) = 0.9*wt_opt;
                
                % Run tva with the data
                fprintf('Running variability analysis\n')
                tva_wt = runTMinMax(wt_tmodel, wt_tmodel.varNames(NF_idx),...
                    'runParallel',PAR_FLAG);
                
                % Run flux sampling
                if m_idx == 1
                    fprintf('Running sampling\n')
                    tic
                    [fluxSamples_wt,concSamples_wt] = sampleTModel(wt_tmodel,tva_wt(:,1),tva_wt(:,2),N_SAMPLES,N_CPU);
                    t_end = toc;
                    fprintf('Sampling time: %.2f s\n',t_end)
                end
                
                %% solve TFA for mutant model
                mut_tmodel = this_tmodel;
                mut_tmodel.var_lb(ismember(mut_tmodel.varNames,strcat('NF_', ko_rxns{m_idx}))) = 0;
                mut_tmodel.var_ub(ismember(mut_tmodel.varNames,strcat('NF_', ko_rxns{m_idx}))) = 0;
                
                % run FVA for mutant model
                mut_model = model;
                mut_model.lb(ismember(mut_model.rxns,ko_rxns{m_idx})) = 0;
                mut_model.ub(ismember(mut_model.rxns,ko_rxns{m_idx})) = 0;
                mutFBASol = optimizeCbModel(mut_model);
                mut_model.lb(mut_model.c==1) = 0.9*mutFBASol.f;
                fva_mut = runMinMax(mut_model,mut_model.rxns,'runParallel',PAR_FLAG);
                is_bd_fva_mut = (n(fva_mut));
                clear mut_model mutFBASol
                
                % add biomass ratio constraint
                mut_tmodel.var_ub(mut_tmodel.f==1) = wt_opt / biomass_ratio + 1e-10;
                mut_tmodel.var_lb(mut_tmodel.f==1) = wt_opt / biomass_ratio - 1e-10;
                
                % add relaxed metabolite concentration ranges to measured
                % metabolites
                mut_tmodel = addRelaxedMetConcRanges(mut_tmodel,LC_varNames,log(C_lb_mut),log(C_ub_mut));
                
                % add minimization objective for relaxation variables
                bio_obj = mut_tmodel.f;
                mut_tmodel.f(startsWith(mut_tmodel.varNames,'EPS_')) = -1e-3;
                
                % solve TFA
                tfa_mut = solveTFAmodelCplex(mut_tmodel);
                
                % reset model objective
                mut_tmodel.f = bio_obj;
                
                if isempty(tfa_mut.x)
                    fprintf('Mutant thermo model cannot be solved with relaxed metabolite concentration ranges\n')
                    continue
                end
                
                % positive and negative concentration relaxations
                eps_minus_mut = tfa_mut.x(startsWith(mut_tmodel.varNames,'EPS_MINUS'));
                eps_plus_mut = tfa_mut.x(startsWith(mut_tmodel.varNames,'EPS_PLUS'));
                ch_conc_idx = eps_minus_mut>1e-4 | eps_plus_mut>1e-4;
                
                % (relaxed) TFA metabolite concentrations
                met_conc_mut = tfa_mut.x(cellfun(@(x)find(ismember(mut_tmodel.varNames,x)),...
                    LC_varNames));
                
                % add second panel to figure
                nexttile
                labels = categorical(strrep(erase(LC_varNames,'LC_'),'_','\_'));
                arrayfun(@(i)line([labels(i) labels(i)],[C_lb_mut(i) C_ub_mut(i)],'linewidth',4,'color','k'),1:numel(LC_varNames))
                hold on
                h = zeros(2,1);
                h(1)=scatter(labels(~ch_conc_idx),exp(met_conc_mut(~ch_conc_idx)),15,'filled');
                h(2)=scatter(labels(ch_conc_idx),exp(met_conc_mut(ch_conc_idx)),15,'filled');
                hold off
                set(gca,'YScale','log','FontSize',10)
                legend(h,{'unchanged','relaxed'},'FontSize',14,'box','off',...
                    'location','southeast')
                ylabel('TFA metabolite concentration [M]','FontSize',14)
                text(0.01,0.98,['{\it ' mutants{m_idx} '}'],'units','normalized','fontweight','bold')
                
                exportgraphics(tmp_fig,[res_dir filesep 'metabolite_concentration_' mutants{m_idx} ...
                    '_' l_cond '_timepoint_' num2str(tp(t_idx)) '_phUB_' num2str(ph_ub(l_idx),3) '.png'])
                delete(tmp_fig)
                
                % fix metabolite measured metabolite concentrations to newly
                % obtained values with 10% tolerance
                mut_tmodel.var_lb(cellfun(@(x)find(ismember(mut_tmodel.varNames,x)),...
                    LC_varNames)) = met_conc_mut-abs(.1*(met_conc_mut));
                mut_tmodel.var_ub(cellfun(@(x)find(ismember(mut_tmodel.varNames,x)),...
                    LC_varNames)) = met_conc_mut+abs(.1*(met_conc_mut));
                
                % fix biomass flux between 90% and 100% of mutant growth
                mut_tmodel.var_lb(mut_tmodel.f==1) = 0.9*(wt_opt / biomass_ratio) - 1e-10;
                
                % Run tva with the data
                fprintf('Running variability analysis\n')
                tva_mut = runTMinMax(mut_tmodel, wt_tmodel.varNames(NF_idx),...
                    'runParallel',PAR_FLAG);
                
                % Run flux sampling
                fprintf('Running sampling\n')
                tic
                [fluxSamples_mut,concSamples_mut] = sampleTModel(mut_tmodel,tva_mut(:,1),tva_mut(:,2),N_SAMPLES,N_CPU);
                t_end = toc;
                fprintf('Sampling time: %.2f s\n',t_end)
                
                %% find non-overlapping reactions between wild type and mutant
                nf_idx = find(startsWith(wt_tmodel.varNames,'NF_'));
                non_overlapping = false(numel(nf_idx),1);
                for i=1:numel(nf_idx)
                    if tva_mut(i,1) > tva_wt(i,2) + 1e-9
                        non_overlapping(i) = 1;
                    elseif tva_wt(i,1) > tva_mut(i,2) + 1e-9
                        non_overlapping(i) = 1;
                    end
                end
                
                rxn_idx = cellfun(@(x)findRxnIDs(model,x),erase(wt_tmodel.varNames(nf_idx),'NF_'));
                
                writetable(...
                    [cell2table(...
                    [model.rxns(rxn_idx(non_overlapping)),...
                    model.rxnNames(rxn_idx(non_overlapping)),...
                    model.subSystems(rxn_idx(non_overlapping))],...
                    'VariableNames', {'ReactionID','ReactionName','ReactionSubSystem'}...
                    ),...
                    array2table(...
                    [tva_wt(non_overlapping,1) tva_wt(non_overlapping,2),...
                    tva_mut(non_overlapping,1) tva_mut(non_overlapping,2)],...
                    'VariableNames', {...
                    'minFluxWT','maxFluxWT','minFluxMUT','maxFluxMUT'...
                    })],...
                    [res_dir filesep 'diff_rxns_tfa_' mutants{m_idx} ...
                    '_' l_cond{lc_idx} '_t_' num2str(tp(t_idx)) ...
                    '_phUB_' num2str(ph_ub(l_idx),3) '.csv']...
                    );
                
                
                %% find reactions that show differences in flexibility between wild type and mutant
                is_bd_tva_wt = (n(tva_wt));
                is_bd_tva_mut = (n(tva_mut));
                diff_dir_idx = xor(is_bd_tva_wt, is_bd_tva_mut);
                
                dir_tab = [...
                    cell2table([...
                    model.rxns(diff_dir_idx),...
                    model.rxnNames(diff_dir_idx),...
                    model.subSystems(diff_dir_idx)],...
                    'VariableNames', {'ReactionID','ReactionName','ReactionSubSystem'}...
                    ),...
                    array2table([...
                    is_bd_fva_wt(diff_dir_idx),...
                    is_bd_fva_mut(diff_dir_idx),...
                    is_bd_tva_wt(diff_dir_idx),...
                    is_bd_tva_mut(diff_dir_idx),...
                    fva_wt(diff_dir_idx,1),fva_wt(diff_dir_idx,2),...
                    fva_mut(diff_dir_idx,1),fva_mut(diff_dir_idx,2),...
                    tva_wt(diff_dir_idx,1),tva_wt(diff_dir_idx,2),...
                    tva_mut(diff_dir_idx,1),tva_mut(diff_dir_idx,2)],...
                    'VariableNames',{...
                    'IS_BD_FVA_WT','IS_BD_FVA_MUT','IS_BD_TVA_WT','IS_BD_TVA_Mut',...
                    'minFlux_FVA_WT','maxFlux_FVA_WT','minFlux_FVA_Mut','maxFlux_FVA_Mut',...
                    'minFlux_TVA_WT','maxFlux_TVA_WT','minFlux_TVA_Mut','maxFlux_TVA_Mut'}...
                    )];
                
                writetable(dir_tab,[res_dir filesep 'rxn_flexibility_' ...
                    mutants{m_idx} ...
                    '_' l_cond{lc_idx} ...
                    '_t_' num2str(tp(t_idx)) ...
                    '_pHUB_' num2str(ph_ub(l_idx),3) '.csv'])
                
                %% write FVA and TVA results to file
                writetable(...
                    array2table([...
                    fva_wt fva_mut tva_wt tva_mut ...
                    ],...
                    'VariableNames',{...
                    'minFlux_FVA_WT','maxFlux_FVA_WT','minFlux_FVA_Mut','maxFlux_FVA_Mut',...
                    'minFlux_TVA_WT','maxFlux_TVA_WT','minFlux_TVA_Mut','maxFlux_TVA_Mut'},...
                    'RowNames',model.rxns),...
                    [res_dir filesep 'va_results_' mutants{m_idx} ...
                    '_' l_cond{lc_idx} ...
                    '_t_' num2str(tp(t_idx)) ...
                    '_pHUB_' num2str(ph_ub(l_idx),3) '.csv'],...
                    'WriteRowNames', true)
                
                %% write relaxed metabolite concentrations to file
                writetable(...
                    array2table([...
                    C_lb_wt C_ub_wt exp(met_conc_wt) eps_minus_wt eps_plus_wt ...
                    C_lb_mut C_ub_mut exp(met_conc_mut) eps_minus_mut eps_plus_mut],...
                    'VariableNames',...
                    {'CONC_LB_WT','CONC_UB_WT','TFA_CONC_WT','E_MINUS_LOG_WT','E_PLUS_LOG_WT',...
                    'CONC_LB_MUT','CONC_UB_MUT','TFA_CONC_MUT','E_MINUS_LOG_MUT','E_PLUS_LOG_MUT'},...
                    'RowNames', erase(LC_varNames,'LC_')),...
                    [res_dir filesep 'met_conc_' mutants{m_idx} ...
                    '_' l_cond{lc_idx} ...
                    '_t_' num2str(tp(t_idx)) ...
                    '_pHUB_' num2str(ph_ub(l_idx),3) '.csv'],...
                    'WriteRowNames',true)
                
                %% write sampling results to file
                if m_idx == 1
                    % wild type (only once per time point and light intensity
                    
                    % flux
                    writematrix(...
                        fluxSamples_wt,...
                        [res_dir filesep 'flux_samples_Col-0_' ...
                        l_cond{lc_idx} ...
                        '_t_' num2str(tp(t_idx)) ...
                        '_pHUB_' num2str(ph_ub(l_idx),3) '.csv'])
                    
                    % concentration
                    writematrix(...
                        concSamples_wt,...
                        [res_dir filesep 'conc_samples_Col-0_' ...
                        l_cond{lc_idx} ...
                        't_' num2str(tp(t_idx)) ...
                        '_pHUB_' num2str(ph_ub(l_idx),3) '.csv'])
                end
                
                % mutant
                
                % flux
                writematrix(...
                    fluxSamples_mut,...
                    [res_dir filesep 'flux_samples_' mutants{m_idx} ...
                    '_' l_cond{lc_idx} ...
                    '_t_' num2str(tp(t_idx)) ...
                    '_pHUB_' num2str(ph_ub(l_idx),3) '.csv'])
                
                % concentration
                writematrix(...
                    concSamples_mut,...
                    [res_dir filesep 'conc_samples_' mutants{m_idx} ...
                    '_' l_cond{lc_idx} ...
                    '_t_' num2str(tp(t_idx)) ...
                    '_pHUB_' num2str(ph_ub(l_idx),3) '.csv'])
                
            end
        end
    end
end