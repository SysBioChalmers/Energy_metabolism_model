% Calculate protein cost for each pathway
% Related to Fig 2A
% Related to Fig S1
% The unit of protein cost is g/gCDW/flux.

% cd ../cobratoolbox;
% initCobraToolbox;

addpath('Functions/');

%% Modify model
% E.coli
model = xls2model('Model_ecoli.xlsx');
model.S(contains(model.mets,'h2o[c]'),contains(model.rxns,'Biomass')) = 0;
model.S(contains(model.mets,'atp[c]'),contains(model.rxns,'Biomass')) = 0;
model.S(contains(model.mets,'pi[c]'),contains(model.rxns,'Biomass')) = 0;
model.S(contains(model.mets,'adp[c]'),contains(model.rxns,'Biomass')) = 0;
model.S(contains(model.mets,'h[c]'),contains(model.rxns,'Biomass')) = 0;
model = changeRxnBounds(model, 'ATPM', 0, 'l');
model = changeRxnBounds(model, 'ATPM', 1000, 'u');
model_ecoli = model;

% Yeast
model = xls2model('Model_yeast.xlsx');
model.S(contains(model.mets,'h2o[c]'),contains(model.rxns,'Biomass')) = 0;
model.S(contains(model.mets,'atp[c]'),contains(model.rxns,'Biomass')) = 0;
model.S(contains(model.mets,'pi[c]'),contains(model.rxns,'Biomass')) = 0;
model.S(contains(model.mets,'adp[c]'),contains(model.rxns,'Biomass')) = 0;
model.S(contains(model.mets,'h[c]'),contains(model.rxns,'Biomass')) = 0;
model = changeRxnBounds(model, 'ATPM', 0, 'l');
model = changeRxnBounds(model, 'ATPM', 1000, 'u');
model_yeast = model;

clear model;

%% Import protein cost data
prot_cost_info = struct();
[num, txt, ~] = xlsread('Model_ecoli.xlsx','Protein_cost_info');
prot_cost_info.id = txt(2:end,1);
prot_cost_info.value = num;
clear num txt;
prot_cost_ecoli = prot_cost_info;

prot_cost_info = struct();
[num, txt, ~] = xlsread('Model_yeast.xlsx','Protein_cost_info');
prot_cost_info.id = txt(2:end,1);
prot_cost_info.value = num;
clear num txt;
prot_cost_yeast = prot_cost_info;

clear prot_cost_info;

%% E.coli
% HY
model = model_ecoli;
model = changeRxnBounds(model, 'GLCtex_HY', 1, 'b');
model = changeRxnBounds(model, 'GLCtex_LY', 0, 'b');
sol = solveModel(model,'ATPM','max',prot_cost_ecoli,1000);
[HYp_ecoli, ~, ~] = calculateProtein(model, prot_cost_ecoli, sol.fluxes);
% LY
model = model_ecoli;
model = changeRxnBounds(model, 'GLCtex_HY', 0, 'b');
model = changeRxnBounds(model, 'GLCtex_LY', 1, 'b');
sol = solveModel(model,'ATPM','max',prot_cost_ecoli,1000);
[~, LYp_ecoli, ~] = calculateProtein(model, prot_cost_ecoli, sol.fluxes);
% Bio
model = model_ecoli;
model = changeRxnBounds(model, 'EXglc', -1, 'b');
model = changeRxnBounds(model, 'GLCtex_HY', 0, 'b');
model = changeRxnBounds(model, 'GLCtex_LY', 0, 'b');
sol = solveModel(model,'Biomass','max',prot_cost_ecoli,1000);
[~, ~, Biop_ecoli] = calculateProtein(model, prot_cost_ecoli, sol.fluxes);

clear model sol;

%% Yeast
% HY
model = model_yeast;
model = changeRxnBounds(model, 'GLCt1_HY', 1, 'b');
model = changeRxnBounds(model, 'GLCt1_LY', 0, 'b');
sol = solveModel(model,'ATPM','max',prot_cost_yeast,1000);
[HYp_yeast, ~, ~] = calculateProtein(model, prot_cost_yeast, sol.fluxes);
% LY
model = model_yeast;
model = changeRxnBounds(model, 'GLCt1_HY', 0, 'b');
model = changeRxnBounds(model, 'GLCt1_LY', 1, 'b');
sol = solveModel(model,'ATPM','max',prot_cost_yeast,1000);
[~, LYp_yeast, ~] = calculateProtein(model, prot_cost_yeast, sol.fluxes);
% Bio
model = model_yeast;
model = changeRxnBounds(model, 'EXglc', -1, 'b');
model = changeRxnBounds(model, 'GLCt1_HY', 0, 'b');
model = changeRxnBounds(model, 'GLCt1_LY', 0, 'b');
sol = solveModel(model,'Biomass','max',prot_cost_yeast,1000);
[~, ~, Biop_yeast] = calculateProtein(model, prot_cost_yeast, sol.fluxes);

clear model sol;

clear model_ecoli model_yeast prot_cost_ecoli prot_cost_yeast;