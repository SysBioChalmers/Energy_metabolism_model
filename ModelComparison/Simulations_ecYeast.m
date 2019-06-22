% ecYeast
% The model was downloaded from GECKO v1.3.0 (ecYeastGEM_batch_7.6_v2.2)
% model = readCbModel('ecYeastGEM_batch.xml');

%% Simulations
% Ref
load('ecYeastGEM.mat');

model = changeRxnBounds(model,'r_1714_REV',0,'b'); 
model = changeRxnBounds(model,'r_1714',-1000,'l'); 
model = changeObjective(model,'r_2111');
sol_ref = optimizeCbModel(model,'max','one');
mu_ref = sol_ref.f;

% Change protein cost one by one

c = 0.5; % change the protein cost to 50% of the original value

[~, txt, ~] = xlsread('RxnList.xlsx','Yeast');
rxn_list = txt(:,2);
clear txt;

mu_list = zeros(length(rxn_list),1);
for i = 1:length(rxn_list)
    model_tmp = model;
    rxnid = rxn_list(i);
    rxn_idx = find(ismember(model_tmp.rxns,rxnid));
    coef_idx = model_tmp.S(:,rxn_idx);
	met_idx = find(coef_idx);
    metlist = model_tmp.mets(coef_idx ~= 0);
    coeflist = model_tmp.S(met_idx,rxn_idx);
    
    if ~any(contains(metlist,'prot_'))
        mu_list(i) = NaN;
    else
        newcoeflist = coeflist;
        newcoeflist(contains(metlist,'prot_')) = newcoeflist(contains(metlist,'prot_'))*c;
        model_tmp.S(met_idx,rxn_idx) = newcoeflist;
        sol_tmp = optimizeCbModel(model_tmp,'max','one');
        mu_list(i) = sol_tmp.f;
    end
end

rel_values = mu_list./mu_ref;