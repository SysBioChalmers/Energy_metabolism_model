% Simulate yeast metabolic switch
% Search for the ratio of LY/HY for best-fit simulations

% Related to Fig 3BD
% Related to Fig S2B

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

mu_list = 0.01:0.005:0.5;

f_hy_yeast = 1;
f_ly_yeast = 1;

% Estimate protein allocation
prot_list = zeros(1,length(bch_data));
for i = 1:length(bch_data)
    model = model_yeast;
    model = changeRxnBounds(model, 'EXbiomass', bch_data(1,i), 'b');
    model = changeRxnBounds(model, 'EXglc', -bch_data(2,i), 'b');
    model = changeRxnBounds(model, 'EXetoh', bch_data(3,i), 'b');
    [minProt, fluxes] = minimizeProtein(model, prot_cost_yeast, 'ATPM', 'max', f_hy_yeast,f_ly_yeast);
    if ~isempty(fluxes)
        prot_list(1,i) = minProt;
    else
        [minProt, ~] = obtainFeasibleSol(model,prot_cost_yeast,'Yeast',0.01,'ATPM','max',f_hy_yeast,f_ly_yeast);
        prot_list(1,i) = minProt;
    end
end
min_prot_yeast = median(prot_list);
clear fluxes minProt model prot_list i;

% Generate figure

fluxes_sim_yeast = zeros(3,length(mu_list));
for i = 1:length(mu_list)
    mu = mu_list(i);
    model = changeRxnBounds(model_yeast, 'EXbiomass', mu, 'b');
    sol = solveModel(model,'EXglc','max',prot_cost_yeast,min_prot_yeast,f_hy_yeast,f_ly_yeast);
    if sol.exitflag == 1
        glc = -sol.fluxes(strcmp(model.rxns,'EXglc'));
        eth = sol.fluxes(strcmp(model.rxns,'EXetoh'));
        fluxes_sim_yeast(:,i) = [mu; glc; eth];
    end
end
fluxes_sim_yeast = fluxes_sim_yeast(:,1:length(find(fluxes_sim_yeast(1,:))));
clear eth glc i min_prot_list mu sol;

figure('Name','yeast_1');
hold on;
box on;
plot(chm_data(2,:),chm_data(3,:),'o','LineWidth',0.75,'Color','k','MarkerSize',8);
plot(bch_data(2,:),bch_data(3,:),'^','LineWidth',0.75,'Color','k','MarkerSize',10);
plot(fluxes_sim_yeast(2,:),fluxes_sim_yeast(3,:),'-','LineWidth',0.75,'Color','k');

xlim([0 25]);

set(gca,'FontSize',12,'FontName','Helvetica');
xlabel('Glucose uptake (mmol/gCDW/h)','FontSize',14,'FontName','Helvetica');
ylabel(['Ethanol production',char(13,10)','(mmol/gCDW/h)'],'FontSize',14,'FontName','Helvetica');

legend({'Chemostat data',...
        'Batch data',...
        'Predicted data'},'FontSize',12,'FontName','Helvetica','location','nw');


set(gcf,'position',[300 400 240 185]);
set(gca,'position',[0.2 0.18 0.76 0.8]);

clear bch_data chm_data exp_data fluxes_sim_yeast model mu_list;
clear f_hy_yeast f_ly_yeast min_prot_yeast;

%% Figures for chemostats

% Yeast (PMID: 9603825)
fluxes_exp_yeast = [0.1  0.15  0.2  0.25  0.28  0.3   0.32  0.35  0.36  0.38 ; % mu
                    1.1  1.67  2.15 2.83  3.24  3.7   5.44  8.09  8.33  10.23; % glucose
                    0    0     0    0     0     0.51  4.42  6.91  6.71  14.91; % ethanol
                    2.73 2.5   5.07 6.8   8.3   8.8   6.83  6.6   7.1   4.19];% o2
% Batch data (PMID: 20199578)
bch_data = [0.46  % Growth rate
            22.4  % Glucose rate
            36.1];% Ethanol rate

exp_data = [fluxes_exp_yeast(1:3,:) bch_data];

x = exp_data(1,6:end);
y = exp_data(3,6:end);
p_exp_yx = polyfit(y,x,1);
crtc_exp = polyval(p_exp_yx,0);
p_exp_xy = polyfit(x,y,1);
slope_exp = p_exp_xy(1);
clear x y p_exp_xy p_exp_yx;

% Original simulations
mu_list = 0.01:0.005:0.5;

f_hy_org = 1;
f_ly_org = 1;


model = model_yeast;
model = changeRxnBounds(model, 'EXbiomass', bch_data(1), 'b');
model = changeRxnBounds(model, 'EXglc', -bch_data(2), 'b');
model = changeRxnBounds(model, 'EXetoh', bch_data(3), 'b');
[minProt, fluxes] = minimizeProtein(model, prot_cost_yeast, 'ATPM', 'max', f_hy_org, f_ly_org);
if ~isempty(fluxes)
	min_prot_org = minProt;
else
	[minProt, ~] = obtainFeasibleSol(model,prot_cost_yeast,'Yeast',0.01,'ATPM','max',f_hy_org,f_ly_org);
	min_prot_org = minProt;
end

clear fluxes i minProt model;

fluxes_sim_yeast = zeros(3,length(mu_list));
for i = 1:length(mu_list)
    mu = mu_list(i);
    model = changeRxnBounds(model_yeast, 'EXbiomass', mu, 'b');
    sol = solveModel(model,'EXglc','max',prot_cost_yeast,min_prot_org,f_hy_org,f_ly_org);
    if sol.exitflag == 1
        glc = -sol.fluxes(strcmp(model.rxns,'EXglc'));
        eth = sol.fluxes(strcmp(model.rxns,'EXetoh'));
        fluxes_sim_yeast(:,i) = [mu; glc; eth];
    end
end
fluxes_sim_yeast = fluxes_sim_yeast(:,1:length(find(fluxes_sim_yeast(1,:))));
clear ac glc i mu sol;

start_tmp = find(fluxes_sim_yeast(3,:), 1);
x = fluxes_sim_yeast(1,start_tmp:end);
y = fluxes_sim_yeast(3,start_tmp:end);
p_sim_yx = polyfit(y,x,1);
crtc_sim = polyval(p_sim_yx,0);
p_sim_xy = polyfit(x,y,1);
slope_sim = p_sim_xy(1);
clear start_tmp x y p_sim_xy p_sim_yx;

% Test protein efficiency
f_list = 0.2:0.2:1;
color_list = [200,200,200;150,150,150;100,100,100;50,50,50;0,0,0]/255;
figure('Name','yeast_test');
% Decrease the protein efficiency of HY pathway
for j = 1:length(f_list)
    f = f_list(j);
    model = model_yeast;
    model = changeRxnBounds(model, 'EXbiomass', bch_data(1), 'b');
    model = changeRxnBounds(model, 'EXglc', -bch_data(2), 'b');
    model = changeRxnBounds(model, 'EXetoh', bch_data(3), 'b');
    [minProt, fluxes] = minimizeProtein(model, prot_cost_yeast, 'ATPM', 'max', f, 1);
    if ~isempty(fluxes)
        min_prot_org = minProt;
    else
        [minProt, ~] = obtainFeasibleSol(model,prot_cost_yeast,'Yeast',0.01,'ATPM','max',f,1);
        min_prot_org = minProt;
    end
    

    fluxes_sim_yeast = zeros(3,length(mu_list));
    for i = 1:length(mu_list)
        mu = mu_list(i);
        model = changeRxnBounds(model_yeast, 'EXbiomass', mu, 'b');
        sol = solveModel(model,'EXglc','max',prot_cost_yeast,min_prot_org,f,1);
        if sol.exitflag == 1
            glc = -sol.fluxes(strcmp(model.rxns,'EXglc'));
            eth = sol.fluxes(strcmp(model.rxns,'EXetoh'));
            fluxes_sim_yeast(:,i) = [mu; glc; eth];
        end
    end
    fluxes_sim_yeast = fluxes_sim_yeast(:,1:length(find(fluxes_sim_yeast(1,:))));
    hold on;
    plot(fluxes_sim_yeast(1,:),fluxes_sim_yeast(3,:),'Color',color_list(j,:),'LineWidth',2);

end
set(gca,'FontSize',10,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',12,'FontName','Helvetica');
ylabel('Ethanol flux (mmol/gCDW/h)','FontSize',12,'FontName','Helvetica');

set(gcf,'position',[200 0 220 185]);
set(gca,'position',[0.15 0.18 0.8 0.8]);

clear f fluxes_sim_yeast j min_prot_org eth glc i mu sol fluxes  minProt model f_list;


% Change protein efficiency of HY to fit critical point
while abs(crtc_sim-crtc_exp) > 0.001
    if crtc_sim > crtc_exp
        f_hy_org = f_hy_org * 0.95;
    else
        f_hy_org = f_hy_org * 1.05;
    end
    

        model = model_yeast;
        model = changeRxnBounds(model, 'EXbiomass', bch_data(1), 'b');
        model = changeRxnBounds(model, 'EXglc', -bch_data(2), 'b');
        model = changeRxnBounds(model, 'EXetoh', bch_data(3), 'b');
        [minProt, fluxes] = minimizeProtein(model, prot_cost_yeast, 'ATPM', 'max', f_hy_org, f_ly_org);
        if ~isempty(fluxes)
            min_prot_org = minProt;
        else
            [minProt, ~] = obtainFeasibleSol(model,prot_cost_yeast,'Yeast',0.01,'ATPM','max',f_hy_org,f_ly_org);
            min_prot_org = minProt;
        end

    
    fluxes_sim_yeast = zeros(3,length(mu_list));
    for i = 1:length(mu_list)
        mu = mu_list(i);
        model = changeRxnBounds(model_yeast, 'EXbiomass', mu, 'b');
        sol = solveModel(model,'EXglc','max',prot_cost_yeast,min_prot_org,f_hy_org,f_ly_org);
        if sol.exitflag == 1
            glc = -sol.fluxes(strcmp(model.rxns,'EXglc'));
            eth = sol.fluxes(strcmp(model.rxns,'EXetoh'));
            fluxes_sim_yeast(:,i) = [mu; glc; eth];
        end
    end
    fluxes_sim_yeast = fluxes_sim_yeast(:,1:length(find(fluxes_sim_yeast(1,:))));

    start_tmp = find(fluxes_sim_yeast(3,:), 1);
    x = fluxes_sim_yeast(1,start_tmp:end);
    y = fluxes_sim_yeast(3,start_tmp:end);
    p_sim_yx = polyfit(y,x,1);
    crtc_sim = polyval(p_sim_yx,0);
    p_sim_xy = polyfit(x,y,1);
    slope_sim = p_sim_xy(1);
end

clear eth fluxes fluxes_sim_yeast glc i minProt model mu min_prot_org;
clear p_sim_xy p_sim_yx sol start_tmp x y;

f_hy_yeast = f_hy_org;
f_ly_yeast = f_ly_org;

% Estimate protein allocation
    model = model_yeast;
    model = changeRxnBounds(model, 'EXbiomass', bch_data(1), 'b');
    model = changeRxnBounds(model, 'EXglc', -bch_data(2), 'b');
    model = changeRxnBounds(model, 'EXetoh', bch_data(3), 'b');
    [minProt, fluxes] = minimizeProtein(model, prot_cost_yeast, 'ATPM', 'max', f_hy_yeast,f_ly_yeast);
    if ~isempty(fluxes)
        min_prot_yeast = minProt;
    else
        [minProt, ~] = obtainFeasibleSol(model,prot_cost_yeast,'Yeast',0.01,'ATPM','max',f_hy_yeast,f_ly_yeast);
        min_prot_yeast = minProt;
    end
clear fluxes minProt model i;

mu_list = 0.01:0.005:0.5;
fluxes_sim_yeast = zeros(8,length(mu_list));
for i = 1:length(mu_list)
    mu = mu_list(i);
    model = changeRxnBounds(model_yeast, 'EXbiomass', mu, 'b');
    sol = solveModel(model,'EXglc','max',prot_cost_yeast,min_prot_yeast,f_hy_yeast,f_ly_yeast);
    if sol.exitflag == 1
        glc = abs(sol.fluxes(strcmp(model.rxns,'EXglc')));
        eth = sol.fluxes(strcmp(model.rxns,'EXetoh'));
        o2 = abs(sol.fluxes(strcmp(model.rxns,'EXo2')));
        [HYp, LYp, Biop] = calculateProtein(model, prot_cost_yeast, sol.fluxes, f_hy_yeast,f_ly_yeast);
        fluxes_sim_yeast(:,i) = [mu; glc; eth; o2; HYp; LYp; Biop; HYp+LYp+Biop];
    end
end
fluxes_sim_yeast = fluxes_sim_yeast(:,1:length(find(fluxes_sim_yeast(1,:))));
clear eth glc i model mu mu_list o2 sol HYp LYp Biop;
clear f_hy_org f_ly_org;

figure('Name','yeast_chemostat');
hold on;
box on;
plot(fluxes_exp_yeast(1,:),fluxes_exp_yeast(2,:),'o','LineWidth',0.75,'Color',[55,126,184]/255,'MarkerSize',8);
plot(fluxes_exp_yeast(1,:),fluxes_exp_yeast(3,:),'o','LineWidth',0.75,'Color',[255,127,0]/255,'MarkerSize',8);
plot(fluxes_exp_yeast(1,:),fluxes_exp_yeast(4,:),'o','LineWidth',0.75,'Color',[77,175,74]/255,'MarkerSize',8);
plot(fluxes_sim_yeast(1,:),fluxes_sim_yeast(2,:),'-','LineWidth',0.75,'Color',[55,126,184]/255);
plot(fluxes_sim_yeast(1,:),fluxes_sim_yeast(3,:),'-','LineWidth',0.75,'Color',[255,127,0]/255);
plot(fluxes_sim_yeast(1,:),fluxes_sim_yeast(4,:),'-','LineWidth',0.75,'Color',[77,175,74]/255);

xlim([0 0.5]);

set(gca,'FontSize',12,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
ylabel('Flux (mmol/gCDW/h)','FontSize',14,'FontName','Helvetica');
legend({'Glucose uptake',...
        'Ethanol production'...
        'O2 uptake',},'FontSize',12,'FontName','Helvetica','location','nw');

set(gcf,'position',[0 400 240 185]);
set(gca,'position',[0.2 0.18 0.76 0.8]);