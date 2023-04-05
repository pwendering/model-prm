function model = updateModelHprGgt(model)
% This function adds the cytosolic hydroxypyruvate dehydrogenase reaction 
% catalysed by the HPR2 protein. Therefore, a transport reaction for 
% hydroxypyruvate from the peroxisome to the cytosol compartment is added 
% as well.
% Further, HPR1 and GGT1 are removed from the gene association of chloroplastic
% reactions.

% correct associations for HPR genes, add HPR2 reaction in cytosol
model = addMetabolite(model, 'HPR[c]', 'metName', 'Hydroxypyruvate',...
    'metFormula', 'C3H3O4');
model.metisinchikeyID{end} = 'HHDDCCUIIUWNGJ-UHFFFAOYSA-N';
model = addReaction(model, 'Tr_HPR1',...
    'reactionName', 'Hydroxypyruvate transporter',...
    'reactionFormula', 'HPR[p] <=> HPR[c]');
model = addReaction(model, 'GCEADH_c',...
    'reactionName', 'GCEA dehydrogenase',...
    'reactionFormula', 'HPR[c] + NADPH[c] + H[c] -> GCEA[c] + NADP[c]',...
    'geneRule', 'AT1G79870');

% remove HPR1 from gene associations of chloroplastic reactions
model = changeGeneAssociation(model, 'GCEADH_h', 'AT1G12550');

% remove GGT1 from gene associations of chloroplastic reactions
model = changeGeneAssociation(model, 'AlaTA_h', 'AT1G70580 or AT1G17290');
model = changeGeneAssociation(model, 'GGAT_h', 'AT1G70580');

end