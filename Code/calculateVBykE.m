function C = calculateVBykE(model,v)
%% C = calculateVBykE(model,v,rxn_idx,enz_idx)
% Calculate ratio of reaction flux divided by the draw reaction fluxes of
% the associated enzymes and kcat values:
%                    v
%           C = ------------ .
%                k_cat * e
% 
% In case of multiple unique values obtained from k_cat*E for one reaction,
% the average value is assigned to the reaction.
% 
% INPUT:
%     struct model:         enzyme-constrained metabolic model (GECKO format)
%     double v:             flux vector
% OUTPUT:
%     double C:             calculated ratio between v and k_cat * E

% get indices of metabolic reactions and enzyme draw reactions
[rxn_idx, enz_idx] = getRxnEnzIdx(model);

% initialize C as cell array since we can get multiple values for the same
% reaction
C = cell(numel(rxn_idx),1);

for i=1:numel(C)
    
    % get current reaction index
    tmp_rxn_idx = rxn_idx(i);
    
    if v(tmp_rxn_idx) ~= 0
        
        % determine consumed enzyme metabolites for current reaction
        tmp_enz_idx = find(model.S(enz_idx,tmp_rxn_idx) < 0);
        
        % array to store all C values for current reaction
        tmp_c = zeros(numel(tmp_enz_idx),1);
        
        for j=1:numel(tmp_enz_idx)
            
            % get Swissprit IDs of current enzyme metabolite
            tmp_enz_id = erase(...
                model.mets(enz_idx(tmp_enz_idx(j))),...
                'prot_'...
                );
            
            % find draw reaction associated with current enzyme
            tmp_draw_rxn_idx = endsWith(model.rxns,tmp_enz_id);
            
            % calculate C for current reaction-enzyme pair
            tmp_c(j) = ...
                v(tmp_rxn_idx) / ... % v [mmol/gDW/h]
                (-1/model.S(enz_idx(tmp_enz_idx(j)),tmp_rxn_idx)) / ... % kcat [/h]
                v(tmp_draw_rxn_idx); % E [mmol/gDW]
            
        end
        
        % determine unique C
        uniq_tmp_c = unique(tmp_c);
        
        if numel(uniq_tmp_c)>1
            fprintf('%s: %d unique fractions\n Assigning the average of %.2g\n',...
                char(tmp_enz_id), numel(uniq_tmp_c), mean(uniq_tmp_c))
            C{i} = mean(uniq_tmp_c);
        else
            C{i} = uniq_tmp_c;
        end
        
    end
end

fprintf('-\tC could not be determined for %d reactions (%d%%)\n',...
    sum(cellfun(@isempty,C)),...
    round(100*sum(cellfun(@isempty,C))/numel(C),0))
fprintf('-\t%d reactions have C == 0\n',...
    sum(cell2mat(cellfun(@(x)x==0,C,'un',0))))
fprintf('-\t%d reactions have C == 1\n',...
    sum(cell2mat(cellfun(@(x)x==1,C,'un',0))))
fprintf('-\t%d reactions have 0 < C < 1\n',...
    sum(cell2mat(cellfun(@(x)x>0&x<1,C,'un',0))))

% assign NaN to reactions that did not carry flux in given flux distribution
C(cellfun(@isempty,C)) = {NaN};

% convert to double matrix 
C = cell2mat(C);

end