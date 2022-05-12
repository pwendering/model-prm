%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Modeling of photorespiratory mutants using TFA
% Philipp Wendering, University of Potsdam (philipp.wendering@gmail.com)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc

% basic variables
NCPU=4;

% set up parallel pool
delete(gcp('nocreate'));
parpool(NCPU);

% initialize COBRA
initCobraToolbox(false)

% CPLEX path and COBRA solver
cplexPath = 'C:\Program Files\IBM\ILOG\CPLEX_Studio129\cplex\matlab\x64_win64';
addpath(genpath(cplexPath))
changeCobraSolver('ibm_cplex')

% path to matTFA toolbox
mattfa_path = fullfile('..','..','mattfa');

addpath(genpath(fullfile(mattfa_path,'matTFA')))
addpath(genpath(fullfile(mattfa_path,'thermoDatabases')))
addpath(genpath(fullfile(mattfa_path,'models')))
addpath(genpath(fullfile(mattfa_path,'plotting')))

% define output directory
res_dir = '../Results';

%% Load the model
% metabolic model (updated metabolite names and added inchi keys)
% model reference: Arnold and Nikoloski (2014), doi: 10.1104/pp.114.235358
load(fullfile('..','Data','AraCore-updated-rev.mat'));
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
tmp = load('thermo_data.mat');
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
model.ub(findRxnIDs(model,'Im_hnu')) = I_ML;

%% Test simple fba
solFBA = optimizeCbModel(model);
% We can set a lower bound for growth (e.g. 50% of maximal growth)
min_obj = roundsd(0.5*solFBA.f, 2, 'floor');
model.lb(model.c==1) = min_obj;

%% Perform FVA
fva_wt = runMinMax(model,model.rxns,false);
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
    fva_wt = runMinMax(model);
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

%% Read experimental data
data_dir = '../Data';
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

% set objective to light-limiting biomass reaction
model.c(:) = 0;
model.c(findRxnIDs(model,'Bio_opt')) = 1;
model.ub(findRxnIDs(model,{'Bio_CLim','Bio_NLim'})) = 0;

% dry weights
dw = readtable(...
    fullfile(data_dir, data_files{contains(data_files, 'dw_')}),...
    'ReadRowNames', true...
    );

% net CO2 assimilation rates
A_ml = readtable(...
    fullfile(data_dir, data_files{contains(data_files, 'A_ml')})...
    );

A_fl = readtable(...
    fullfile(data_dir, data_files{contains(data_files, 'A_fl')})...
    );

A_max = max([A_fl.mean_A; A_ml.mean_A]);


% growth duration
t_growth = 34 * 24; % [h]

% calculate growth rate
dw_0 = 20 / 1e6;
mu = (log(dw.mean_DW) - log(dw_0)) ./ t_growth;
mu_wt = mu(strcmp(dw.Properties.RowNames,'Col0'));

mu(mu>mu_wt) = mu_wt;

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

% light conditions
l_conds = {'ml','fl','ml_fl','fl_ml'};

% time points
tp = [-21 3 27 75 123];

l_idx = 1;

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

%% loop over timepoints
for t_idx = 1:numel(tp)
    cmdsz = matlab.desktop.commandwindow.size;
    fprintf([repmat('-',1,cmdsz(1)) '\n'])
    fprintf('Time after shift: %d\n', tp(t_idx))
    for m_idx = 1:numel(mutants)
        
        cmdsz = matlab.desktop.commandwindow.size;
        fprintf([repmat('-',1,floor(cmdsz(1)/2)) '\n'])
        fprintf('>>> Processing genotype %s ...\n', mutants{m_idx})
        
        % test if growth is possible with knock-outs
        ko_model = removeRxns(model,ko_rxns{m_idx});
        g_ko = optimizeCbModel(ko_model).f;
        g_wt = optimizeCbModel(model).f;
        gt_bool = g_ko > 0;
        if gt_bool
            fprintf('Growth test passed\n')
            fprintf('FBA results: WT: %.4e, %s: %.4e\n', g_ko, mutants{m_idx}, g_wt)
        else
            fprintf('KO disables growth prediction - skipping\n')
            continue
        end
        
        % determine ratio of growth rates (wt/mutant)
        biomass_ratio = mu_wt / ...
            mu(strcmp(dw.Properties.RowNames,mutants{m_idx}));
        
        %% get metabolite concentration ranges
        met_av = readtable(fullfile(data_dir, ['met_' l_conds{l_idx} '.csv']));
        met_av_mut = met_av{strcmp(met_av.genotype,mutants{m_idx})&met_av.time==tp(t_idx),3:end};
        met_av_wt = met_av{strcmp(met_av.genotype,'Col-0')&met_av.time==tp(t_idx),3:end};
        met_sd = readtable(fullfile(data_dir, ['met_' l_conds{l_idx} '_sd.csv']));
        met_sd_mut = met_sd{strcmp(met_av.genotype,mutants{m_idx})&met_av.time==tp(t_idx),3:end};
        met_sd_wt= met_sd{strcmp(met_av.genotype,'Col-0')&met_av.time==tp(t_idx),3:end};
        
        meas_names = met_av.Properties.VariableNames(3:end);
        meas_ids = name2id(meas_names);
        
        % remove data for unmatched metabolites
        nm_idx = cellfun(@isempty,meas_ids);
        met_av_mut = met_av_mut(~nm_idx);
        met_av_wt = met_av_wt(~nm_idx);
        met_sd_mut = met_sd_mut(~nm_idx);
        met_sd_wt = met_sd_wt(~nm_idx);
        meas_ids = meas_ids(~nm_idx);
        
        meas_mw = zeros(size(meas_ids));
        for i=1:numel(meas_ids)
            first_idx = find(contains(model.mets,[meas_ids{i} '[']),1);
            meas_mw(i) = getMolecularMass(model.metFormulas(first_idx)); % [mol/gDW]
        end
        
        meas_conc_min_mut = met_av_mut - met_sd_mut;
        meas_conc_min_wt = met_av_wt - met_sd_wt;
        meas_conc_max_mut = met_av_mut + met_sd_mut;
        meas_conc_max_wt = met_av_wt + met_sd_wt;
        
        meas_conc_min_mut = meas_conc_min_mut' * dw_per_l ./ meas_mw; % [M]
        meas_conc_min_wt = meas_conc_min_wt' * dw_per_l ./ meas_mw; % [M]
        meas_conc_max_mut = meas_conc_max_mut' * dw_per_l ./ meas_mw; % [M]
        meas_conc_max_wt = meas_conc_max_wt' * dw_per_l ./ meas_mw; % [M]
        
        
        meas_conc_min_mut(meas_conc_min_mut<0) = min(model.CompartmentData.compMinConc);
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
        C_lb_mut = zeros(size(LC_varNames));
        C_lb_wt = zeros(size(LC_varNames));
        C_ub_mut = zeros(size(LC_varNames));
        C_ub_wt = zeros(size(LC_varNames));
        
        for i=1:numel(meas_ids)
            meas_id_match = ismember(regexprep(LC_varNames,'LC_(.*)_\w','$1'),meas_ids(i));
            if sum(meas_id_match)>0
                
                C_lb_mut(meas_id_match) = meas_conc_min_mut(i);
                C_lb_wt(meas_id_match) = meas_conc_min_wt(i);
                
                C_ub_mut(meas_id_match) = meas_conc_max_mut(i);
                C_ub_wt(meas_id_match) = meas_conc_max_wt(i);
            else
                fprintf('%s not found in model variables\n',meas_ids{i})
            end
        end
        
        % find the indices of these variables in the variable names of the tfa
        id_LC_varNames = find_cell(LC_varNames, this_tmodel.varNames);
        
        %% Solve TFA for wild type model
        wt_tmodel = this_tmodel;
        
        % add relaxed metabolite concentration ranges to measured
        % metabolites
        wt_tmodel = addRelaxedMetConcRanges(wt_tmodel,LC_varNames,log(C_lb_wt),log(C_ub_wt));
        
        % add minimization objective for relaxation variables
        wt_tmodel.f(startsWith(wt_tmodel.varNames,'EPS_')) = -1;
        
        % solve TFA
        tfa_wt = solveTFAmodelCplex(wt_tmodel);
        
        % WT optimal growth rate
        wt_opt = tfa_wt.x(wt_tmodel.f==1);
        
        % positive and negative concentration relaxations
        eps_minus = tfa_wt.x(startsWith(wt_tmodel.varNames,'EPS_MINUS'));
        eps_plus = tfa_wt.x(startsWith(wt_tmodel.varNames,'EPS_PLUS'));
        ch_conc_idx = eps_minus>1e-4 | eps_plus>1e-4;
        
        % (relaxed) TFA metabolite concentrations
        met_conc = tfa_wt.x(cellfun(@(x)find(ismember(wt_tmodel.varNames,x)),...
            LC_varNames));
        
        % create figure with updated metabolite concentration ranges
        tmp_fig = figure('Visible','off');
        tiledlayout(2,1)
        nexttile
        labels = categorical(erase(strrep(LC_varNames,'_','\_'),'LC_'));
        arrayfun(@(i)line([labels(i) labels(i)],[C_lb_wt(i) C_ub_wt(i)],'linewidth',4,'color','k'),1:numel(LC_varNames))
        hold on
        h = zeros(2,1);
        h(1)=scatter(labels(~ch_conc_idx),exp(met_conc(~ch_conc_idx)),15,'filled');
        h(2)=scatter(labels(ch_conc_idx),exp(met_conc(ch_conc_idx)),15,'filled');
        hold off
        set(gca,'YScale','log','FontSize',10)
        legend(h,{'unchanged','relaxed'},'FontSize',14,'box','off',...
            'location','southeast')
        ylabel('TFA metabolite concentration [M]','FontSize',14)
        text(0.01,0.95,'Col-0','units','normalized','fontweight','bold')
        
        % fix metabolite measured metabolite concentrations to newly
        % obtained values with 10% tolerance
        wt_tmodel.var_lb(cellfun(@(x)find(ismember(wt_tmodel.varNames,x)),...
            LC_varNames)) = met_conc-abs(.1*(met_conc));
        wt_tmodel.var_ub(cellfun(@(x)find(ismember(wt_tmodel.varNames,x)),...
            LC_varNames)) = met_conc+abs(.1*(met_conc));
        
        % set lower bound for biomass to 90% of optimal value
        wt_tmodel.var_lb(wt_tmodel.f==1) = 0.9*wt_opt;
        
        % Run tva with the data
        if ~isempty(tfa_wt.x)
            tva_wt = runTMinMax(wt_tmodel, wt_tmodel.varNames(NF_idx));
        end
        
        %% solve TFA for mutant model
        mut_tmodel = this_tmodel;
        mut_tmodel.var_lb(ismember(mut_tmodel.varNames,strcat('NF_', ko_rxns{m_idx}))) = 0;
        mut_tmodel.var_ub(ismember(mut_tmodel.varNames,strcat('NF_', ko_rxns{m_idx}))) = 0;
        
        % run FVA for mutant model
        mut_model = model;
        mut_model.lb(ismember(mut_model.rxns,ko_rxns{m_idx})) = 0;
        mut_model.ub(ismember(mut_model.rxns,ko_rxns{m_idx})) = 0;
        fva_mut = runMinMax(mut_model);
        is_bd_fva_mut = (n(fva_mut));
        
        % add biomass ratio constraint
        mut_tmodel.var_ub(mut_tmodel.f==1) = wt_opt / biomass_ratio + 1e-10;
        mut_tmodel.var_lb(mut_tmodel.f==1) = wt_opt / biomass_ratio - 1e-10;
        
        % add relaxed metabolite concentration ranges to measured
        % metabolites
        mut_tmodel = addRelaxedMetConcRanges(mut_tmodel,LC_varNames,log(C_lb_wt),log(C_ub_wt));
        
        % add minimization objective for relaxation variables
        mut_tmodel.f(startsWith(mut_tmodel.varNames,'EPS_')) = -1e-3;
        
        % solve TFA
        tfa_mut = solveTFAmodelCplex(mut_tmodel);
        
        % positive and negative concentration relaxations
        eps_minus = tfa_mut.x(startsWith(mut_tmodel.varNames,'EPS_MINUS'));
        eps_plus = tfa_mut.x(startsWith(mut_tmodel.varNames,'EPS_PLUS'));
        ch_conc_idx = eps_minus>1e-4 | eps_plus>1e-4;
        
        % (relaxed) TFA metabolite concentrations
        met_conc = tfa_mut.x(cellfun(@(x)find(ismember(mut_tmodel.varNames,x)),...
            LC_varNames));
        
        % add second panel to figure
        nexttile
        labels = categorical(erase(strrep(LC_varNames,'_','\_'),'LC_'));
        arrayfun(@(i)line([labels(i) labels(i)],[C_lb_mut(i) C_ub_mut(i)],'linewidth',4,'color','k'),1:numel(LC_varNames))
        hold on
        h = zeros(2,1);
        h(1)=scatter(labels(~ch_conc_idx),exp(met_conc(~ch_conc_idx)),15,'filled');
        h(2)=scatter(labels(ch_conc_idx),exp(met_conc(ch_conc_idx)),15,'filled');
        hold off
        set(gca,'YScale','log','FontSize',10)
        legend(h,{'unchanged','relaxed'},'FontSize',14,'box','off',...
            'location','southeast')
        ylabel('TFA metabolite concentration [M]','FontSize',14)
        text(0.01,0.95,['{\it ' mutants{m_idx} '}'],'units','normalized','fontweight','bold')
        set(tmp_fig,'Outerposition',1000*[1.5297   -0.2063    1.5507    0.9347])
        
        exportgraphics(tmp_fig,[res_dir filesep 'metabolite_concentration_' mutants{m_idx} ...
            '_timepoint_' num2str(tp(t_idx)) '.png'])
        
        % fix metabolite measured metabolite concentrations to newly
        % obtained values with 10% tolerance
        mut_tmodel.var_lb(cellfun(@(x)find(ismember(mut_tmodel.varNames,x)),...
            LC_varNames)) = met_conc-abs(.1*(met_conc));
        mut_tmodel.var_ub(cellfun(@(x)find(ismember(mut_tmodel.varNames,x)),...
            LC_varNames)) = met_conc+abs(.1*(met_conc));
        
        % fix biomass flux between 90% and 100% of mutant growth
        mut_tmodel.var_lb(mut_tmodel.f==1) = 0.9*(wt_opt / biomass_ratio) - 1e-10;
        
        % Run tva with the data
        if ~isempty(tfa_mut.x)
            tva_mut = runTMinMax(mut_tmodel, mut_tmodel.varNames(NF_idx));
        end
        
        %% find non-overlapping reactions between wild type and mutant
        nf_idx = find(startsWith(wt_tmodel.varNames,'NF_'));
        non_overlapping = false(numel(nf_idx),1);
        for i=1:numel(nf_idx)
            if tva_mut(i,1) > tva_wt(i,2) + 1e-8
                non_overlapping(i) = 1;
            elseif tva_wt(i,1) > tva_mut(i,2) + 1e-8
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
            't_' num2str(tp(t_idx)) ...
            '_' l_conds{l_idx} '.csv']...
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
            fva_mut(diff_dir_idx,1),fva_mut(diff_dir_idx,1),...
            tva_wt(diff_dir_idx,1),tva_wt(diff_dir_idx,2),...
            tva_mut(diff_dir_idx,1),tva_mut(diff_dir_idx,2)],...
            'VariableNames',{...
            'IS_BD_FVA_WT','IS_BD_FVA_MUT','IS_BD_TVA_WT','IS_BD_TVA_Mut',...
            'minFlux_FVA_WT','maxFlux_FVA_WT','minFlux_FVA_MUT','maxFlux_FVA_MUT',...
            'minFlux_WT','maxFlux_WT','minFlux_Mut','maxFlux_Mut'}...
            )];
        
        writetable(dir_tab,[res_dir filesep 'rxn_flx_' mutants{m_idx} '_t_' ...
            num2str(tp(t_idx)) '.csv'])
        
    end
end