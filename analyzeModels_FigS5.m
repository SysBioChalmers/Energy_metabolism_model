% Analyse models
% Change ratio of LY/HY
% Related to Fig S5

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
model = changeRxnBounds(model, 'EXbiomass', 0, 'b');
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
model = changeRxnBounds(model, 'EXbiomass', 0, 'b');
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

%% Simulations 
% Use adjusted ratio, i.e., ecoli is 2.31 yeast is 1.66.

% E.coli
tot_prot_list_ecoli = 0.1;
glc_list_ecoli = 0.1:0.1:20;
fluxes_ecoli = zeros(5,length(glc_list_ecoli),length(tot_prot_list_ecoli));
for i = 1:length(tot_prot_list_ecoli)
    tot_prot = tot_prot_list_ecoli(i);
    
    for j = 1:length(glc_list_ecoli)
        glc_in = glc_list_ecoli(j);
        model = changeRxnBounds(model_ecoli, 'EXglc', -glc_in, 'b');
        sol = solveModel(model,'ATPM','max',prot_cost_ecoli,tot_prot,0.5883);
        if sol.exitflag == 1
            HY = sol.fluxes(strcmp(model.rxns,'GLCtex_HY'))*23.5;
            LY = sol.fluxes(strcmp(model.rxns,'GLCtex_LY'))*11;
            ATP = sol.fluxes(strcmp(model.rxns,'ATPM'));
            [HYp, LYp, Biop] = calculateProtein(model, prot_cost_ecoli, sol.fluxes,0.5883);
            prot = HYp + LYp + Biop;
            fluxes_ecoli(:,j,i) = [HY; LY; ATP; glc_in; prot];
        end
    end
    
end
fluxes_ecoli(fluxes_ecoli<0.0001) = 0;
clear ATP glc_in HY i j LY model sol tot_prot HYp LYp Biop prot;

% Yeast
tot_prot_list_yeast = 0.1;
glc_list_yeast = 0.1:0.1:40;
fluxes_yeast = zeros(5,length(glc_list_yeast),length(tot_prot_list_yeast));
for i = 1:length(tot_prot_list_yeast)
    tot_prot = tot_prot_list_yeast(i);
    
    for j = 1:length(glc_list_yeast)
        glc_in = glc_list_yeast(j);
        model = changeRxnBounds(model_yeast, 'EXglc', -glc_in, 'b');
        sol = solveModel(model,'ATPM','max',prot_cost_yeast,tot_prot,0.6914);
        if sol.exitflag == 1
            HY = sol.fluxes(strcmp(model.rxns,'GLCt1_HY'))*22;
            LY = sol.fluxes(strcmp(model.rxns,'GLCt1_LY'))*2;
            ATP = sol.fluxes(strcmp(model.rxns,'ATPM'));
            [HYp, LYp, Biop] = calculateProtein(model, prot_cost_yeast, sol.fluxes,0.6914);
            prot = HYp + LYp + Biop;
            fluxes_yeast(:,j,i) = [HY; LY; ATP; glc_in; prot];
        end
    end
    
end
fluxes_yeast(fluxes_yeast<0.0001) = 0;
clear ATP glc_in HY i j LY model sol tot_prot HYp LYp Biop prot;

figure('Name','ecoli_2.31');

data = fluxes_ecoli(:,:,1);
num = length(find(data(3,:) > 0));
data = data(:,1:num);
x = data(4,:);
y1_tot = data(3,:);
y1_ly = data(2,:);

hold on;
plot(x,y1_tot,'k','LineWidth',1.5);
h = fill([x,fliplr(x)],[y1_ly,zeros(1,length(x))],[215,25,28]/256);
set(h,'edgealpha',0,'facealpha',0.7);
h = fill([x,fliplr(x)],[y1_tot,fliplr(y1_ly)],[43,131,186]/256);
set(h,'edgealpha',0,'facealpha',0.7);
xlim([0 13]);
ylim([0 150]);
xlabel(['Glucose uptake rate',char(13,10)','(mmol/gCDW/h)'],'FontSize',12,'FontName','Helvetica');
ylabel(['ATP production rate',char(13,10)','(mmol/gCDW/h)'],'FontSize',12,'FontName','Helvetica');
legend({'Total ATP';'ATP from Low-yield pathway';'ATP from High-yield pathway'},'FontSize',8,'FontName','Helvetica','location','ne');

set(gcf,'position',[500 200 180 130]);
set(gca,'position',[0.27 0.33 0.7 0.63]);

figure('Name','yeast_1.66');

data = fluxes_yeast(:,:,1);
num = length(find(data(3,:) > 0));
data = data(:,1:num);
x = data(4,:);
y1_tot = data(3,:);
y1_ly = data(2,:);

hold on;
plot(x,y1_tot,'k','LineWidth',1.5);
h = fill([x,fliplr(x)],[y1_ly,zeros(1,length(x))],[215,25,28]/256);
set(h,'edgealpha',0,'facealpha',0.7);
h = fill([x,fliplr(x)],[y1_tot,fliplr(y1_ly)],[43,131,186]/256);
set(h,'edgealpha',0,'facealpha',0.7);
xlim([0 25]);
ylim([0 50]);
xlabel(['Glucose uptake rate',char(13,10)','(mmol/gCDW/h)'],'FontSize',12,'FontName','Helvetica');
ylabel(['ATP production rate',char(13,10)','(mmol/gCDW/h)'],'FontSize',12,'FontName','Helvetica');
legend({'Total ATP';'ATP from Low-yield pathway';'ATP from High-yield pathway'},'FontSize',8,'FontName','Helvetica','location','ne');

set(gcf,'position',[700 200 180 130]);
set(gca,'position',[0.27 0.33 0.7 0.63]);

%% Simulations 
% Change raio, i.e., for ecoli from 2.31 to 1.66, for yeast from 1.66 to 2.31

% E.coli
tot_prot_list_ecoli = 0.1;
glc_list_ecoli = 0.1:0.1:20;
fluxes_ecoli = zeros(5,length(glc_list_ecoli),length(tot_prot_list_ecoli));
for i = 1:length(tot_prot_list_ecoli)
    tot_prot = tot_prot_list_ecoli(i);
    
    for j = 1:length(glc_list_ecoli)
        glc_in = glc_list_ecoli(j);
        model = changeRxnBounds(model_ecoli, 'EXglc', -glc_in, 'b');
        sol = solveModel(model,'ATPM','max',prot_cost_ecoli,tot_prot,0.5883*2.31/1.66);
        if sol.exitflag == 1
            HY = sol.fluxes(strcmp(model.rxns,'GLCtex_HY'))*23.5;
            LY = sol.fluxes(strcmp(model.rxns,'GLCtex_LY'))*11;
            ATP = sol.fluxes(strcmp(model.rxns,'ATPM'));
            [HYp, LYp, Biop] = calculateProtein(model, prot_cost_ecoli, sol.fluxes,0.5883*2.31/1.66);
            prot = HYp + LYp + Biop;
            fluxes_ecoli(:,j,i) = [HY; LY; ATP; glc_in; prot];
        end
    end
    
end
fluxes_ecoli(fluxes_ecoli<0.0001) = 0;
clear ATP glc_in HY i j LY model sol tot_prot HYp LYp Biop prot;

% Yeast
tot_prot_list_yeast = 0.1;
glc_list_yeast = 0.1:0.1:40;
fluxes_yeast = zeros(5,length(glc_list_yeast),length(tot_prot_list_yeast));
for i = 1:length(tot_prot_list_yeast)
    tot_prot = tot_prot_list_yeast(i);
    
    for j = 1:length(glc_list_yeast)
        glc_in = glc_list_yeast(j);
        model = changeRxnBounds(model_yeast, 'EXglc', -glc_in, 'b');
        sol = solveModel(model,'ATPM','max',prot_cost_yeast,tot_prot,0.6914*1.66/2.31);
        if sol.exitflag == 1
            HY = sol.fluxes(strcmp(model.rxns,'GLCt1_HY'))*22;
            LY = sol.fluxes(strcmp(model.rxns,'GLCt1_LY'))*2;
            ATP = sol.fluxes(strcmp(model.rxns,'ATPM'));
            [HYp, LYp, Biop] = calculateProtein(model, prot_cost_yeast, sol.fluxes,0.6914*1.66/2.31);
            prot = HYp + LYp + Biop;
            fluxes_yeast(:,j,i) = [HY; LY; ATP; glc_in; prot];
        end
    end
    
end
fluxes_yeast(fluxes_yeast<0.0001) = 0;
clear ATP glc_in HY i j LY model sol tot_prot HYp LYp Biop prot;

figure('Name','ecoli_1.66');

data = fluxes_ecoli(:,:,1);
num = length(find(data(3,:) > 0));
data = data(:,1:num);
x = data(4,:);
y1_tot = data(3,:);
y1_ly = data(2,:);

hold on;
plot(x,y1_tot,'k','LineWidth',1.5);
h = fill([x,fliplr(x)],[y1_ly,zeros(1,length(x))],[215,25,28]/256);
set(h,'edgealpha',0,'facealpha',0.7);
h = fill([x,fliplr(x)],[y1_tot,fliplr(y1_ly)],[43,131,186]/256);
set(h,'edgealpha',0,'facealpha',0.7);
xlim([0 13]);
ylim([0 150]);
xlabel(['Glucose uptake rate',char(13,10)','(mmol/gCDW/h)'],'FontSize',12,'FontName','Helvetica');
ylabel(['ATP production rate',char(13,10)','(mmol/gCDW/h)'],'FontSize',12,'FontName','Helvetica');
% legend({'Total ATP';'ATP from Low-yield pathway';'ATP from High-yield pathway'},'FontSize',8,'FontName','Helvetica','location','ne');

set(gcf,'position',[100 200 180 130]);
set(gca,'position',[0.27 0.33 0.7 0.63]);

figure('Name','yeast_2.31');

data = fluxes_yeast(:,:,1);
num = length(find(data(3,:) > 0));
data = data(:,1:num);
x = data(4,:);
y1_tot = data(3,:);
y1_ly = data(2,:);

hold on;
plot(x,y1_tot,'k','LineWidth',1.5);
h = fill([x,fliplr(x)],[y1_ly,zeros(1,length(x))],[215,25,28]/256);
set(h,'edgealpha',0,'facealpha',0.7);
h = fill([x,fliplr(x)],[y1_tot,fliplr(y1_ly)],[43,131,186]/256);
set(h,'edgealpha',0,'facealpha',0.7);
xlim([0 25]);
ylim([0 50]);
xlabel(['Glucose uptake rate',char(13,10)','(mmol/gCDW/h)'],'FontSize',12,'FontName','Helvetica');
ylabel(['ATP production rate',char(13,10)','(mmol/gCDW/h)'],'FontSize',12,'FontName','Helvetica');
% legend({'Total ATP';'ATP from Low-yield pathway';'ATP from High-yield pathway'},'FontSize',8,'FontName','Helvetica','location','ne');

set(gcf,'position',[200 200 180 130]);
set(gca,'position',[0.27 0.33 0.7 0.63]);

clear data h num x y1_ly y1_tot;