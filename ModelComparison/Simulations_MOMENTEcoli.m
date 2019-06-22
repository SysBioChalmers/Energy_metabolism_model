% MOMENT model of E.coli 
% load('Ecoli_model.mat');
% momentEcoli = reformulateEcoliModel(model);

prot_cost_ecoli = struct();
[num, txt, ~] = xlsread('Model_ecoli.xlsx','Protein_cost_info');
prot_cost_ecoli.id = txt(2:end,1);
prot_cost_ecoli.value = num;
clear num txt;
prot_cost_ecoli.id(ismember(prot_cost_ecoli.id,'CYTBO34pp')) = {'CYTBO3_4pp'};

%% Simulations
load('Ecoli_model.mat');
% Change gene ID in the model.rules into that in the model.genes.
for i = 1:length(model.genes)
    old = strcat('x(',num2str(i),')');
    new = model.genes{i};
    model.rules = cellfun(@(x) strrep(x,old,new),...
                            model.rules,'UniformOutput',false);
end

% calculate protein costs
prot_cost_info = momentProteinCost(model);
% % change the protein costs to manually determined values for energy part
% for i = 1:length(prot_cost_ecoli.id)
%     id_tmp = prot_cost_ecoli.id(i);
%     if ismember(id_tmp,prot_cost_info.id)
%         prot_cost_info.value(ismember(prot_cost_info.id,id_tmp)) = prot_cost_ecoli.value(i);
%     end
% end

model.protCost = prot_cost_info.value;

% split reverse reactions
model = reformulateEcoliModel(model);
model = changeRxnBounds(model, 'EX_glc__D_e', -1000, 'l');

tot_prot_weight = 0.15;
prot_cost_info = model.protCost;

% solve LP
sol = solveMOMENT(model,'BIOMASS_Ec_iML1515_core_75p37M','max',prot_cost_info,tot_prot_weight);
glc = -sol.fluxes(strcmp(model.rxns,'EX_glc__D_e'));
ac = sol.fluxes(strcmp(model.rxns,'EX_ac_e'));
fluxes_ref_ecoli = [sol.fluxes(strcmp(model.rxns,'BIOMASS_Ec_iML1515_core_75p37M')),glc,ac,(ac*2)/(glc*6)];

c = 0.5; % change the protein cost to 50% of the original value

fluxes_sim_ecoli = zeros(length(prot_cost_ecoli.id),4);
for i = 1:length(prot_cost_ecoli.id)
    model_tmp = model;
    rxnid = prot_cost_ecoli.id(i);
    model_tmp.protCost(ismember(model_tmp.rxns,rxnid)) = model_tmp.protCost(ismember(model_tmp.rxns,rxnid))*c;
    sol = solveMOMENT(model_tmp,'BIOMASS_Ec_iML1515_core_75p37M','max',model_tmp.protCost,tot_prot_weight);
    glc = -sol.fluxes(strcmp(model_tmp.rxns,'EX_glc__D_e'));
    ac = sol.fluxes(strcmp(model_tmp.rxns,'EX_ac_e'));
	fluxes_sim_ecoli(i,:) = [sol.fluxes(strcmp(model_tmp.rxns,'BIOMASS_Ec_iML1515_core_75p37M')),glc,ac,(ac*2)/(glc*6)];
end

% % solve LP (considering energy metabolism constraint)
% 
% prot_cost_EM = prot_cost_info;
% for i = 1:length(model.rxns)
%     rxn_tmp = model.rxns{i};
%     if contains(rxn_tmp,'_fwd') || contains(rxn_tmp,'_rvs')
%         rxn_tmp = rxn_tmp(1:end-4);
%     end
%     if ~ismember(rxn_tmp,prot_cost_ecoli.id)
%         prot_cost_EM(i,1) = 0;
%     end
% end
% 
% EM_weight = 0.04;
% sol = solveMOMENTmod(model,'BIOMASS_Ec_iML1515_core_75p37M','max',prot_cost_info,tot_prot_weight,prot_cost_EM,EM_weight);
% 
% % sol = solveMOMENT(model,'BIOMASS_Ec_iML1515_core_75p37M','max',prot_cost_info,tot_prot_weight);
% glc = -sol.fluxes(strcmp(model.rxns,'EX_glc__D_e'));
% ac = sol.fluxes(strcmp(model.rxns,'EX_ac_e'));
% fluxes_ref_ecoli = [sol.fluxes(strcmp(model.rxns,'BIOMASS_Ec_iML1515_core_75p37M')),glc,ac,(ac*2)/(glc*6)];
% 
% c = 0.5; % change the protein cost to 50% of the original value
% 
% fluxes_sim_ecoli = zeros(length(prot_cost_ecoli.id),4);
% for i = 1:length(prot_cost_ecoli.id)
%     model_tmp = model;
%     rxnid = prot_cost_ecoli.id(i);
%     model_tmp.protCost(ismember(model_tmp.rxns,rxnid)) = model_tmp.protCost(ismember(model_tmp.rxns,rxnid))*c;
%     sol = solveMOMENTmod(model,'BIOMASS_Ec_iML1515_core_75p37M','max',prot_cost_info,tot_prot_weight,prot_cost_EM,EM_weight);
%     glc = -sol.fluxes(strcmp(model_tmp.rxns,'EX_glc__D_e'));
%     ac = sol.fluxes(strcmp(model_tmp.rxns,'EX_ac_e'));
% 	fluxes_sim_ecoli(i,:) = [sol.fluxes(strcmp(model_tmp.rxns,'BIOMASS_Ec_iML1515_core_75p37M')),glc,ac,(ac*2)/(glc*6)];
% end

value = fluxes_sim_ecoli(:,1)/fluxes_ref_ecoli(1);
value = round(value,2);