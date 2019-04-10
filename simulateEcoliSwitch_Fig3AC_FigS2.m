% Simulate E. coli metabolic switch
% Search for the ratio of LY/HY for best-fit simulations

% Related to Fig 3AC
% Related to Fig S2A

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

mu_list = 0.01:0.01:1;

f_hy_ecoli = 1;
f_ly_ecoli = 1;

% Estimate protein allocation
prot_list = zeros(1,length(bch_data));
for i = 1:length(bch_data)
    model = model_ecoli;
    model = changeRxnBounds(model, 'EXbiomass', bch_data(1,i), 'b');
    model = changeRxnBounds(model, 'EXglc', -bch_data(2,i), 'b');
    model = changeRxnBounds(model, 'EXac', bch_data(3,i), 'b');
    [minProt, fluxes] = minimizeProtein(model, prot_cost_ecoli, 'ATPM', 'max', f_hy_ecoli,f_ly_ecoli);
    if ~isempty(fluxes)
        prot_list(1,i) = minProt;
    else
        [minProt, ~] = obtainFeasibleSol(model,prot_cost_ecoli,'Ecoli',0.01,'ATPM','max',f_hy_ecoli,f_ly_ecoli);
        prot_list(1,i) = minProt;
    end
end
min_prot_ecoli = median(prot_list);
clear fluxes minProt model prot_list i;

% Generate figure

fluxes_sim_ecoli = zeros(3,length(mu_list));
for i = 1:length(mu_list)
    mu = mu_list(i);
    model = changeRxnBounds(model_ecoli, 'EXbiomass', mu, 'b');
    sol = solveModel(model,'EXglc','max',prot_cost_ecoli,min_prot_ecoli,f_hy_ecoli,f_ly_ecoli);
    if sol.exitflag == 1
        glc = -sol.fluxes(strcmp(model.rxns,'EXglc'));
        ac = sol.fluxes(strcmp(model.rxns,'EXac'));
        fluxes_sim_ecoli(:,i) = [mu; glc; ac];
    end
end
fluxes_sim_ecoli = fluxes_sim_ecoli(:,1:length(find(fluxes_sim_ecoli(1,:))));
clear ac glc i min_prot_list mu sol;

figure('Name','ecoli_1');
hold on;
box on;
plot(chm_data(2,:),chm_data(3,:),'o','LineWidth',1.5,'Color','k','MarkerSize',6);
plot(bch_data(2,:),bch_data(3,:),'^','LineWidth',1.5,'Color','k','MarkerSize',7);
plot(fluxes_sim_ecoli(2,:),fluxes_sim_ecoli(3,:),'-','LineWidth',1.5,'Color','k');

set(gca,'FontSize',10,'FontName','Helvetica');
xlabel('Glucose uptake (mmol/gCDW/h)','FontSize',12,'FontName','Helvetica');
ylabel(['Acetate production',char(13,10)','(mmol/gCDW/h)'],'FontSize',12,'FontName','Helvetica');

xlim([0 12]);

legend({'Chemostat data',...
        'Batch data',...
        'Predicted data'},'FontSize',9,'FontName','Helvetica','location','nw');

set(gcf,'position',[300 400 240 185]);
set(gca,'position',[0.2 0.18 0.76 0.8]);

clear bch_data chm_data exp_data fluxes_sim_ecoli model mu_list;

clear f_hy_ecoli f_ly_ecoli min_prot_ecoli;

%% Figures for chemostats

% E.coli (PMID: 16672514/21694717)
fluxes_exp_ecoli = [0.04  0.1   0.19  0.29  0.41  0.49  0.6;   % mu
                    0.89  1.22  2.20  3.34  4.71  5.93  7.86;  % glucose
                    0     0     0     0     0     1.81  6.02;  % acetate
                    2.45  2.94  4.87  10.12 10.57 15.25 12.94];% o2
% Batch data (PMID: 27136056)
bch_data = [0.65  % Growth rate
            9.65  % Glucose rate
            6.76];% Acetate rate

exp_data = [fluxes_exp_ecoli(1:3,:) bch_data];

x = exp_data(1,6:end);
y = exp_data(3,6:end);
p_exp_yx = polyfit(y,x,1);
crtc_exp = polyval(p_exp_yx,0);
p_exp_xy = polyfit(x,y,1);
slope_exp = p_exp_xy(1);
clear x y p_exp_xy p_exp_yx;

% Original simulations
mu_list = 0.01:0.01:1;

f_hy_org = 1;
f_ly_org = 1;

model = model_ecoli;
model = changeRxnBounds(model, 'EXbiomass', bch_data(1), 'b');
model = changeRxnBounds(model, 'EXglc', -bch_data(2), 'b');
model = changeRxnBounds(model, 'EXac', bch_data(3), 'b');
[minProt, fluxes] = minimizeProtein(model, prot_cost_ecoli, 'ATPM', 'max', f_hy_org, f_ly_org);
if ~isempty(fluxes)
	min_prot_org = minProt;
else
	[minProt, ~] = obtainFeasibleSol(model,prot_cost_ecoli,'Ecoli',0.01,'ATPM','max',f_hy_org,f_ly_org);
	min_prot_org = minProt;
end
clear fluxes i minProt model;

fluxes_sim_ecoli = zeros(3,length(mu_list));
for i = 1:length(mu_list)
    mu = mu_list(i);
    model = changeRxnBounds(model_ecoli, 'EXbiomass', mu, 'b');
    sol = solveModel(model,'EXglc','max',prot_cost_ecoli,min_prot_org,f_hy_org,f_ly_org);
    if sol.exitflag == 1
        glc = -sol.fluxes(strcmp(model.rxns,'EXglc'));
        ac = sol.fluxes(strcmp(model.rxns,'EXac'));
        fluxes_sim_ecoli(:,i) = [mu; glc; ac];
    end
end
fluxes_sim_ecoli = fluxes_sim_ecoli(:,1:length(find(fluxes_sim_ecoli(1,:))));
clear ac glc i mu sol;

start_tmp = find(fluxes_sim_ecoli(3,:), 1);
x = fluxes_sim_ecoli(1,start_tmp:end);
y = fluxes_sim_ecoli(3,start_tmp:end);
p_sim_yx = polyfit(y,x,1);
crtc_sim = polyval(p_sim_yx,0);
p_sim_xy = polyfit(x,y,1);
slope_sim = p_sim_xy(1);
clear start_tmp x y p_sim_xy p_sim_yx;

% Test protein efficiency
f_list = 0.2:0.2:1;
color_list = [200,200,200;150,150,150;100,100,100;50,50,50;0,0,0]/255;
figure('Name','ecoli_test');
% Decrease the protein efficiency of HY pathway
for j = 1:length(f_list)
    f = f_list(j);
    model = model_ecoli;
    model = changeRxnBounds(model, 'EXbiomass', bch_data(1), 'b');
    model = changeRxnBounds(model, 'EXglc', -bch_data(2), 'b');
    model = changeRxnBounds(model, 'EXac', bch_data(3), 'b');
    [minProt, fluxes] = minimizeProtein(model, prot_cost_ecoli, 'ATPM', 'max', f, 1);
    if ~isempty(fluxes)
        min_prot_org = minProt;
    else
        [minProt, ~] = obtainFeasibleSol(model,prot_cost_ecoli,'Ecoli',0.01,'ATPM','max',f,1);
        min_prot_org = minProt;
    end
    

    fluxes_sim_ecoli = zeros(3,length(mu_list));
    for i = 1:length(mu_list)
        mu = mu_list(i);
        model = changeRxnBounds(model_ecoli, 'EXbiomass', mu, 'b');
        sol = solveModel(model,'EXglc','max',prot_cost_ecoli,min_prot_org,f,1);
        if sol.exitflag == 1
            glc = -sol.fluxes(strcmp(model.rxns,'EXglc'));
            ac = sol.fluxes(strcmp(model.rxns,'EXac'));
            fluxes_sim_ecoli(:,i) = [mu; glc; ac];
        end
    end
    fluxes_sim_ecoli = fluxes_sim_ecoli(:,1:length(find(fluxes_sim_ecoli(1,:))));
    hold on;
    plot(fluxes_sim_ecoli(1,:),fluxes_sim_ecoli(3,:),'Color',color_list(j,:),'LineWidth',2);

end
set(gca,'FontSize',10,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',12,'FontName','Helvetica');
ylabel('Acetate flux (mmol/gCDW/h)','FontSize',12,'FontName','Helvetica');
set(gcf,'position',[200 0 220 185]);
set(gca,'position',[0.15 0.18 0.8 0.8]);

clear f fluxes_sim_ecoli j min_prot_org ac glc i mu sol fluxes  minProt model f_list;


% Change protein efficiency of HY to fit critical point
while abs(crtc_sim-crtc_exp) > 0.001
    if crtc_sim > crtc_exp
        f_hy_org = f_hy_org * 0.95;
    else
        f_hy_org = f_hy_org * 1.05;
    end
    

        model = model_ecoli;
        model = changeRxnBounds(model, 'EXbiomass', bch_data(1), 'b');
        model = changeRxnBounds(model, 'EXglc', -bch_data(2), 'b');
        model = changeRxnBounds(model, 'EXac', bch_data(3), 'b');
        [minProt, fluxes] = minimizeProtein(model, prot_cost_ecoli, 'ATPM', 'max', f_hy_org, f_ly_org);
        if ~isempty(fluxes)
            min_prot_org = minProt;
        else
            [minProt, ~] = obtainFeasibleSol(model,prot_cost_ecoli,'Ecoli',0.01,'ATPM','max',f_hy_org,f_ly_org);
            min_prot_org = minProt;
        end

    
    fluxes_sim_ecoli = zeros(3,length(mu_list));
    for i = 1:length(mu_list)
        mu = mu_list(i);
        model = changeRxnBounds(model_ecoli, 'EXbiomass', mu, 'b');
        sol = solveModel(model,'EXglc','max',prot_cost_ecoli,min_prot_org,f_hy_org,f_ly_org);
        if sol.exitflag == 1
            glc = -sol.fluxes(strcmp(model.rxns,'EXglc'));
            ac = sol.fluxes(strcmp(model.rxns,'EXac'));
            fluxes_sim_ecoli(:,i) = [mu; glc; ac];
        end
    end
    fluxes_sim_ecoli = fluxes_sim_ecoli(:,1:length(find(fluxes_sim_ecoli(1,:))));

    start_tmp = find(fluxes_sim_ecoli(3,:), 1);
    x = fluxes_sim_ecoli(1,start_tmp:end);
    y = fluxes_sim_ecoli(3,start_tmp:end);
    p_sim_yx = polyfit(y,x,1);
    crtc_sim = polyval(p_sim_yx,0);
    p_sim_xy = polyfit(x,y,1);
    slope_sim = p_sim_xy(1);
end


clear ac fluxes fluxes_sim_ecoli glc i minProt model mu min_prot_org;
clear p_sim_xy p_sim_yx sol start_tmp x y;

f_hy_ecoli = f_hy_org;
f_ly_ecoli = f_ly_org;

% Estimate protein allocation

    model = model_ecoli;
    model = changeRxnBounds(model, 'EXbiomass', bch_data(1), 'b');
    model = changeRxnBounds(model, 'EXglc', -bch_data(2), 'b');
    model = changeRxnBounds(model, 'EXac', bch_data(3), 'b');
    [minProt, fluxes] = minimizeProtein(model, prot_cost_ecoli, 'ATPM', 'max', f_hy_ecoli,f_ly_ecoli);
    if ~isempty(fluxes)
        min_prot_ecoli = minProt;
    else
        [minProt, ~] = obtainFeasibleSol(model,prot_cost_ecoli,'Ecoli',0.01,'ATPM','max',f_hy_ecoli,f_ly_ecoli);
        min_prot_ecoli = minProt;
    end

clear fluxes minProt model i;

mu_list = 0.01:0.01:1;
fluxes_sim_ecoli = zeros(8,length(mu_list));
for i = 1:length(mu_list)
    mu = mu_list(i);
    model = changeRxnBounds(model_ecoli, 'EXbiomass', mu, 'b');
    sol = solveModel(model,'EXglc','max',prot_cost_ecoli,min_prot_ecoli,f_hy_ecoli,f_ly_ecoli);
    if sol.exitflag == 1
        glc = abs(sol.fluxes(strcmp(model.rxns,'EXglc')));
        ac = sol.fluxes(strcmp(model.rxns,'EXac'));
        o2 = abs(sol.fluxes(strcmp(model.rxns,'EXo2')));
        [HYp, LYp, Biop] = calculateProtein(model, prot_cost_ecoli, sol.fluxes, f_hy_ecoli,f_ly_ecoli);
        fluxes_sim_ecoli(:,i) = [mu; glc; ac; o2; HYp; LYp; Biop; HYp+LYp+Biop];
    end
end
fluxes_sim_ecoli = fluxes_sim_ecoli(:,1:length(find(fluxes_sim_ecoli(1,:))));
clear ac glc i model mu mu_list o2 sol HYp LYp Biop;
clear f_hy_org f_ly_org;

figure('Name','ecoli_chemostat');
hold on;
box on;
plot(fluxes_exp_ecoli(1,:),fluxes_exp_ecoli(2,:),'o','LineWidth',1.5,'Color',[55,126,184]/255,'MarkerSize',6);
plot(fluxes_exp_ecoli(1,:),fluxes_exp_ecoli(3,:),'o','LineWidth',1.5,'Color',[228,26,28]/255,'MarkerSize',6);
plot(fluxes_exp_ecoli(1,:),fluxes_exp_ecoli(4,:),'o','LineWidth',1.5,'Color',[77,175,74]/255,'MarkerSize',6);
plot(fluxes_sim_ecoli(1,:),fluxes_sim_ecoli(2,:),'-','LineWidth',1.5,'Color',[55,126,184]/255);
plot(fluxes_sim_ecoli(1,:),fluxes_sim_ecoli(3,:),'-','LineWidth',1.5,'Color',[228,26,28]/255);
plot(fluxes_sim_ecoli(1,:),fluxes_sim_ecoli(4,:),'-','LineWidth',1.5,'Color',[77,175,74]/255);

set(gca,'FontSize',10,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',12,'FontName','Helvetica');
ylabel('Flux (mmol/gCDW/h)','FontSize',12,'FontName','Helvetica');
legend({'Glucose uptake',...
        'Acetate production'...
        'O2 uptake',},'FontSize',9,'FontName','Helvetica','location','nw');
xlim([0 0.8]);
set(gcf,'position',[0 400 240 185]);
set(gca,'position',[0.2 0.18 0.76 0.8]);