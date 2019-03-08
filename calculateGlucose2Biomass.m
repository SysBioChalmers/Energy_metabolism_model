% Calculate flux distribution with fixing glucose uptake rate at 1 and
% maximizing growth rate.

% This is for determining biomass formation pathway. 

% cd ../cobratoolbox;
% initCobraToolbox;

%% E.coli
load('Ecoli_model.mat');
% The model is iML1515.

mu = 1;

model_test = changeObjective(model, 'EX_glc__D_e');
model_test = changeRxnBounds(model_test, 'EX_glc__D_e', -1000, 'l'); %glucose exchange
model_test = changeRxnBounds(model_test, 'EX_glc__D_e', 0, 'u'); %glucose exchange
model_test = changeRxnBounds(model_test, 'BIOMASS_Ec_iML1515_core_75p37M', mu, 'b'); %growth rate

model_test = changeRxnBounds(model_test,'DHAPT',0,'b');
model_test = changeRxnBounds(model_test,'GLCt2pp',0,'b');
model_test = changeRxnBounds(model_test,'HEX1',0,'b');
model_test = changeRxnBounds(model_test,'GLCabcpp',0,'b');

model_test = changeRxnBounds(model_test,'PTAr',0,'l');

model_test = changeRxnBounds(model_test,'ATPM',-75.55*mu,'b');

sol_ecoli = optimizeCbModel(model_test,'max','one');

clear model model_test mu;

%% Yeast
load('Yeast_model'); %The NGAM reaction has been added.

mu = 1;

model_test = changeObjective(model, 'r_1714');
model_test = changeRxnBounds(model_test, 'r_1714', -1000, 'l'); %glucose exchange
model_test = changeRxnBounds(model_test, 'r_1714', 0, 'u'); %glucose exchange
model_test = changeRxnBounds(model_test, 'r_2111', mu, 'b'); %growth rate

model_test = changeRxnBounds(model_test,'NGAM',-59.276*mu,'b');

sol_yeast = optimizeCbModel(model_test,'max','one');

clear model model_test mu;