% Calculate flux distribution with fixing glucose uptake rate at 1 and
% maximizing glycerol production rate.

% This is for determining fluxes towards glycerol in yeast model. 

load('Yeast_model'); %The NGAM reaction has been added.

model_test = changeObjective(model, 'r_1808');
model_test = changeRxnBounds(model_test, 'r_1714', -1, 'b'); %glucose exchange

sol_yeast = optimizeCbModel(model_test,'max','one');

clear model model_test mu;