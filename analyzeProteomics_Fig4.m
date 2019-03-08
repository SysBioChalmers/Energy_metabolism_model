% Calculate protein mass based on proteomics data.
% Use chemostat data to calculate minimal protein allocation and apparent
% saturation.

% Related to Fig 4

% cd ../cobratoolbox;
% initCobraToolbox;

addpath('Functions/');

%% Exp data

[num_yeast, ~, ~] = xlsread('Protein_abundance.xlsx','Yeast');
[num_ecoli, ~, ~] = xlsread('Protein_abundance.xlsx','Ecoli');

% Yeast dataset 1 (PMID: 28365149)
tot_prot = sum(num_yeast(3:120,3:12));
glc_rate = num_yeast(122,3:12);

figure('Name','yeast_exp_1');
hold on;
box on;
scatter(glc_rate,tot_prot,40,'ko','filled');
set(gca,'FontSize',10,'FontName','Helvetica');
xlabel('Glucose uptake rate (mmol/gCDW/h)','FontSize',12,'FontName','Helvetica');
ylabel(['Measured protein allocation',char(13,10)','(g/gCDW)'],'FontSize',12,'FontName','Helvetica');
ylim([0 0.15]);
set(gcf,'position',[500 300 300 150]);
set(gca,'position',[0.2 0.22 0.76 0.71]);

tot_prot = sum(num_yeast(3:120,3:21));
glc_rate = num_yeast(122,3:21);
figure('Name','yeast_exp_2');
hold on;
box on;
scatter(glc_rate,tot_prot,40,'ko','filled');
set(gca,'FontSize',10,'FontName','Helvetica');
xlabel('Glucose uptake rate (mmol/gCDW/h)','FontSize',12,'FontName','Helvetica');
ylabel(['Measured protein allocation',char(13,10)','(g/gCDW)'],'FontSize',12,'FontName','Helvetica');
ylim([0 0.15]);
set(gcf,'position',[500 400 300 150]);
set(gca,'position',[0.2 0.22 0.76 0.71]);

% Ecoli dataset 1 (PMID: 25712329)
tot_prot = sum(num_ecoli(3:75,3:24));
glc_rate = num_ecoli(77,3:24);


figure('Name','ecoli_exp');
hold on;
box on;
scatter(glc_rate,tot_prot,40,'ko','filled');
set(gca,'FontSize',10,'FontName','Helvetica');
xlabel('Glucose uptake rate (mmol/gCDW/h)','FontSize',12,'FontName','Helvetica');
ylabel(['Measured protein allocation',char(13,10)','(g/gCDW)'],'FontSize',12,'FontName','Helvetica');
ylim([0 0.15]);
set(gcf,'position',[500 500 300 150]);
set(gca,'position',[0.2 0.22 0.76 0.71]);



%% Yeast CEN.PK 113-7D

% Chemostat data (DOI: 10.1038/s42255-018-0006-7)
% Data were obtained from Fig 3 of (DOI: 10.1038/s42255-018-0006-7), which 
% were originally reported in (PMID: 21354323) and (PMID: 9603825).
chm_data = [0.02	0.05	0.05	0.10	0.15	0.20	0.20	0.25	0.28	0.30	0.30	0.32	0.35	0.35	0.36	0.38	0.38    % Growth rate
            0.25	0.57	0.60	1.10	1.67	2.15	2.22	2.83	3.24	3.70	4.67	5.44	7.62	8.09	8.33	10.23	13.18   % Glucose rate
            0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.51	3.24	4.42	8.72	6.91	6.71	14.91	12.50]; % Ethanol rate

% Batch data
bch_data = [0.40        0.37        0.46  % Growth rate
            19.9        19.6        22.4  % Glucose rate
            29.6        30.5        36.1];% Ethanol rate
% PMID:     19684065    11157958    20199578

exp_data = [chm_data bch_data];

model_yeast = xls2model('Model_yeast.xlsx');

prot_cost_yeast = struct();
[num, txt, ~] = xlsread('Model_yeast.xlsx','Protein_cost_info');
prot_cost_yeast.id = txt(2:end,1);
prot_cost_yeast.value = num;
clear num txt;

model = model_yeast;
% model = changeRxnBounds(model, 'ATPM', 0, 'l');
% model = changeRxnBounds(model, 'ATPM', 1000, 'u');

res = zeros(8,length(exp_data));

for i = 1:length(exp_data)
    mu = exp_data(1,i);
    glc = exp_data(2,i);
    eth = exp_data(3,i);
    model = changeRxnBounds(model, 'EXbiomass', mu, 'b');
    model = changeRxnBounds(model, 'EXglc', -glc, 'b');
    model = changeRxnBounds(model, 'EXetoh', eth, 'b');
    [minProt, fluxes] = minimizeProtein(model, prot_cost_yeast, 'ATPM', 'max',0.6914);
    if ~isempty(fluxes)
        glc_new = -fluxes(strcmp(model.rxns,'EXglc'));
        HY = fluxes(strcmp(model.rxns,'GLCt1_HY'))/glc_new;
        LY = fluxes(strcmp(model.rxns,'GLCt1_LY'))/glc_new;
        Bio = fluxes(strcmp(model.rxns,'GLCt1_Bio'))/glc_new;
        mu_new = fluxes(strcmp(model.rxns,'EXbiomass'));
        [HYp, LYp, Biop] = calculateProtein(model, prot_cost_yeast, fluxes,0.6914);
        res(:,i) = [HY; LY; Bio; HYp; LYp; Biop; glc_new; minProt];
    else
        [minProt, fluxes] = obtainFeasibleSol(model,prot_cost_yeast,'Yeast',0.01,'ATPM','max',0.6914);
        glc_new = -fluxes(strcmp(model.rxns,'EXglc'));
        HY = fluxes(strcmp(model.rxns,'GLCt1_HY'))/glc_new;
        LY = fluxes(strcmp(model.rxns,'GLCt1_LY'))/glc_new;
        Bio = fluxes(strcmp(model.rxns,'GLCt1_Bio'))/glc_new;
        mu_new = fluxes(strcmp(model.rxns,'EXbiomass'));
        [HYp, LYp, Biop] = calculateProtein(model, prot_cost_yeast, fluxes,0.6914);
        res(:,i) = [HY; LY; Bio; HYp; LYp; Biop; glc_new; minProt];
    end
end

figure('Name','ecoli_estimated_protein');
hold on;
box on;
plot(res(7,1:17),res(8,1:17),'o','LineWidth',1.5,'Color','k','MarkerSize',6);
plot(res(7,18:20),res(8,18:20),'^','LineWidth',1.5,'Color','k','MarkerSize',7);
set(gca,'FontSize',10,'FontName','Helvetica');
xlabel('Glucose uptake rate (mmol/gCDW/h)','FontSize',12,'FontName','Helvetica');
ylabel(['Minimal protein allocation',char(13,10)','(g/gCDW)'],'FontSize',12,'FontName','Helvetica');
legend({'Chemostat data',...
        'Batch data'},'FontSize',9,'FontName','Helvetica','location','se');

set(gcf,'position',[0 0 300 150]);
set(gca,'position',[0.2 0.22 0.76 0.71]);

% Plot saturation
max_prot = median(res(8,10:20));
sat = res(8,1:20) ./ max_prot * 100;

figure('Name','yeast_saturation');
hold on;
box on;
plot(res(7,1:17),sat(1:17),'o','LineWidth',1.5,'Color','k','MarkerSize',6);
plot(res(7,18:20),sat(18:20),'^','LineWidth',1.5,'Color','k','MarkerSize',7);
set(gca,'FontSize',10,'FontName','Helvetica');
xlabel('Glucose uptake rate (mmol/gCDW/h)','FontSize',12,'FontName','Helvetica');
ylabel('Apparent saturation (%)','FontSize',12,'FontName','Helvetica');
legend({'Chemostat data',...
        'Batch data'},'FontSize',9,'FontName','Helvetica','location','se');
ylim([0 125]);
set(gcf,'position',[100 400 300 150]);
set(gca,'position',[0.2 0.22 0.76 0.71]);


%% Ecoli MG1655

% Chemostat data (DOI: 10.1038/s42255-018-0006-7)
% Data were obtained from Fig 5 of (DOI: 10.1038/s42255-018-0006-7), which 
% were originally reported in (PMID: 16672514) (PMID: 15838044) (PMID: 21122111)
% (PMID: 16461663) and (PMID: 25712329).
chm_data = [0.04	0.05	0.05	0.06	0.08	0.10	0.11	0.11	0.12	0.19	0.21	0.21	0.21	0.24	0.28	0.29	0.30	0.30	0.32	0.34	0.40	0.40	0.41	0.46	0.51	0.51	0.52	0.54	0.55	0.59	0.65  % Growth rate
            0.89	0.78	1.14	1.14	1.64	1.67	1.70	1.96	1.90	2.81	3.05	2.66	3.08	3.03	3.48	4.22	3.62	4.36	4.59	4.14	4.80	5.05	5.60	5.52	6.58	6.79	6.40	6.76	7.26	7.85	7.60  % Glucose rate
            0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.16	0.00	0.00	0.00	0.16	0.68	1.33	0.47	1.66	0.58	1.27	1.99	2.51	6.11	4.80];% Acetate rate
% Batch data
bch_data = [0.65        0.69        0.61        0.66  % Growth rate
            9.65        8.59        7.54        8.6  % Glucose rate
            6.76        3.91        2.72        4.3];% Acetate rate
% PMID:     27136056    25304508    24249002    21481254

exp_data = [chm_data bch_data];

model_ecoli = xls2model('Model_ecoli.xlsx');

prot_cost_ecoli = struct();
[num, txt, ~] = xlsread('Model_ecoli.xlsx','Protein_cost_info');
prot_cost_ecoli.id = txt(2:end,1);
prot_cost_ecoli.value = num;
clear num txt;

model = model_ecoli;
% model = changeRxnBounds(model, 'ATPM', 0, 'l');
% model = changeRxnBounds(model, 'ATPM', 1000, 'u');

res = zeros(8,length(exp_data));

for i = 1:length(exp_data)
    mu = exp_data(1,i);
    glc = exp_data(2,i);
    ac = exp_data(3,i);
    model = changeRxnBounds(model, 'EXbiomass', mu, 'b');
    model = changeRxnBounds(model, 'EXglc', -glc, 'b');
    model = changeRxnBounds(model, 'EXac', ac, 'b');
    [minProt, fluxes] = minimizeProtein(model, prot_cost_ecoli, 'ATPM', 'max',0.5883);
    if ~isempty(fluxes)
        glc_new = -fluxes(strcmp(model.rxns,'EXglc'));
        HY = fluxes(strcmp(model.rxns,'GLCtex_HY'))/glc_new;
        LY = fluxes(strcmp(model.rxns,'GLCtex_LY'))/glc_new;
        Bio = fluxes(strcmp(model.rxns,'GLCptspp_Bio'))/glc_new;
        mu_new = fluxes(strcmp(model.rxns,'EXbiomass'));
        [HYp, LYp, Biop] = calculateProtein(model, prot_cost_ecoli, fluxes,0.5883);
        res(:,i) = [HY; LY; Bio; HYp; LYp; Biop; glc_new; minProt];
    else
        [minProt, fluxes] = obtainFeasibleSol(model,prot_cost_ecoli,'Ecoli',0.01,'ATPM','max',0.5883);
        glc_new = -fluxes(strcmp(model.rxns,'EXglc'));
        HY = fluxes(strcmp(model.rxns,'GLCtex_HY'))/glc_new;
        LY = fluxes(strcmp(model.rxns,'GLCtex_LY'))/glc_new;
        Bio = fluxes(strcmp(model.rxns,'GLCptspp_Bio'))/glc_new;
        mu_new = fluxes(strcmp(model.rxns,'EXbiomass'));
        [HYp, LYp, Biop] = calculateProtein(model, prot_cost_ecoli, fluxes,0.5883);
        res(:,i) = [HY; LY; Bio; HYp; LYp; Biop; glc_new; minProt];
    end
end

figure('Name','ecoli_estimated_protein');
hold on;
box on;
plot(res(7,1:31),res(8,1:31),'o','LineWidth',1.5,'Color','k','MarkerSize',6);
plot(res(7,32:35),res(8,32:35),'^','LineWidth',1.5,'Color','k','MarkerSize',7);
set(gca,'FontSize',10,'FontName','Helvetica');
xlabel('Glucose uptake rate (mmol/gCDW/h)','FontSize',12,'FontName','Helvetica');
ylabel(['Minimal protein allocation',char(13,10)','(g/gCDW)'],'FontSize',12,'FontName','Helvetica');
legend({'Chemostat data',...
        'Batch data'},'FontSize',9,'FontName','Helvetica','location','se');
ylim([0 0.12]);
set(gcf,'position',[0 400 300 150]);
set(gca,'position',[0.2 0.22 0.76 0.71]);


% Plot saturation
max_prot = median(res(8,21:35));
sat = res(8,1:35) ./ max_prot * 100;

figure('Name','ecoli_saturation');
hold on;
box on;
plot(res(7,1:31),sat(1:31),'o','LineWidth',1.5,'Color','k','MarkerSize',6);
plot(res(7,32:35),sat(32:35),'^','LineWidth',1.5,'Color','k','MarkerSize',7);
set(gca,'FontSize',10,'FontName','Helvetica');
xlabel('Glucose uptake rate (mmol/gCDW/h)','FontSize',12,'FontName','Helvetica');
ylabel('Apparent saturation (%)','FontSize',12,'FontName','Helvetica');
legend({'Chemostat data',...
        'Batch data'},'FontSize',9,'FontName','Helvetica','location','se');
ylim([0 125]);
set(gcf,'position',[200 400 300 150]);
set(gca,'position',[0.2 0.22 0.76 0.71]);

