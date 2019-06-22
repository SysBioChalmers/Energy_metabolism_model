% Change protein cost (increase and decrease) for each protein to see how
% phenotype changes

% Related to Fig 6

c = 2; % change the protein cost to 50% of the original value

%% Ecoli

model_ecoli = xls2model('Model_ecoli.xlsx');

prot_cost_ecoli = struct();
[num, txt, ~] = xlsread('Model_ecoli.xlsx','Protein_cost_info');
prot_cost_ecoli.id = txt(2:end,1);
prot_cost_ecoli.value = num;
clear num txt;

min_prot_ecoli = 0.07;
f_hy_ecoli = 0.5883;
f_ly_ecoli = 1;

model = changeRxnBounds(model_ecoli, 'EXglc', -1000, 'l');
sol = solveModel(model,'EXbiomass','max',prot_cost_ecoli,min_prot_ecoli,f_hy_ecoli,f_ly_ecoli);
glc = -sol.fluxes(strcmp(model.rxns,'EXglc'));
ac = sol.fluxes(strcmp(model.rxns,'EXac'));
fluxes_ref_ecoli = [sol.fluxes(strcmp(model.rxns,'EXbiomass')),glc,ac,(ac*2)/(glc*6)];

fluxes_sim_ecoli = zeros(length(prot_cost_ecoli.id),4);
for i = 1:length(prot_cost_ecoli.id)
    prot_cost_ecoli_tmp = prot_cost_ecoli;
    prot_cost_ecoli_tmp.value(i) = prot_cost_ecoli_tmp.value(i) * c;
    sol = solveModel(model,'EXbiomass','max',prot_cost_ecoli_tmp,min_prot_ecoli,f_hy_ecoli,f_ly_ecoli);
	glc = -sol.fluxes(strcmp(model.rxns,'EXglc'));
	ac = sol.fluxes(strcmp(model.rxns,'EXac'));
	fluxes_sim_ecoli(i,:) = [sol.fluxes(strcmp(model.rxns,'EXbiomass')),glc,ac,(ac*2)/(glc*6)];
end

%% Yeast

model_yeast = xls2model('Model_yeast.xlsx');

prot_cost_yeast = struct();
[num, txt, ~] = xlsread('Model_yeast.xlsx','Protein_cost_info');
prot_cost_yeast.id = txt(2:end,1);
prot_cost_yeast.value = num;
clear num txt;

min_prot_yeast = 0.07;
f_hy_yeast = 0.6914;
f_ly_yeast = 1;

model = changeRxnBounds(model_yeast, 'EXglc', -1000, 'l');
sol = solveModel(model,'EXbiomass','max',prot_cost_yeast,min_prot_yeast,f_hy_yeast,f_ly_yeast);
glc = -sol.fluxes(strcmp(model.rxns,'EXglc'));
eth = sol.fluxes(strcmp(model.rxns,'EXetoh'));
fluxes_ref_yeast = [sol.fluxes(strcmp(model.rxns,'EXbiomass')),glc,eth,(eth*2)/(glc*6)];

fluxes_sim_yeast = zeros(length(prot_cost_yeast.id),4);
for i = 1:length(prot_cost_yeast.id)
    prot_cost_yeast_tmp = prot_cost_yeast;
    prot_cost_yeast_tmp.value(i) = prot_cost_yeast_tmp.value(i) * c;
    sol = solveModel(model,'EXbiomass','max',prot_cost_yeast_tmp,min_prot_yeast,f_hy_yeast,f_ly_yeast);
	glc = -sol.fluxes(strcmp(model.rxns,'EXglc'));
	eth = sol.fluxes(strcmp(model.rxns,'EXetoh'));
	fluxes_sim_yeast(i,:) = [sol.fluxes(strcmp(model.rxns,'EXbiomass')),glc,eth,(eth*2)/(glc*6)];
end

%% Figures
figure();
set(gcf,'position',[450 50 450 300]);

color_ecoli = [142,1,82]/255;
color_yeast = [39,100,25]/255;

% The rxns with assumed kcats will be not displayed in the figures
remove_id_ecoli = {'ACt2rpp'};
remove_id_yeast = {'PYRt2m';'AKGDam';'AKGDbm';'GCC2cm';'SUCOASm';'ATPtmH';'PIt2m';'H2Ot'};

% Ecoli
id = prot_cost_ecoli.id;
value = fluxes_sim_ecoli(:,1)/fluxes_ref_ecoli(1);
idx1 = ~ismember(id,remove_id_ecoli);
id = id(idx1);
value = value(idx1);
value = round(value,2);
[value, idx2] = sort(value,'ascend');
id = id(idx2);

subplot(2,1,1);
h1 = bar(1:length(id),value,'FaceColor',color_ecoli,'FaceAlpha',0.3,'EdgeColor',color_ecoli,'LineWidth',0.5);
set(gca,'XTick',1:1:length(id));
set(gca,'XTickLabel',id);
set(gca,'FontSize',10,'FontName','Helvetica');
ylim([0.6 1.05]);
ylabel('Change in growth','FontSize',12,'FontName','Helvetica','Color','k');
xtickangle(90);
box off;

% Yeast
id = prot_cost_yeast.id;
value = fluxes_sim_yeast(:,1)/fluxes_ref_yeast(1);
idx1 = ~ismember(id,remove_id_yeast);
id = id(idx1);
value = value(idx1);
value = round(value,2);
[value, idx2] = sort(value,'ascend');
id = id(idx2);

subplot(2,1,2);
h2 = bar(1:length(id),value,'FaceColor',color_yeast,'FaceAlpha',0.3,'EdgeColor',color_yeast,'LineWidth',0.5);
set(gca,'XTick',1:1:length(id));
set(gca,'XTickLabel',id);
set(gca,'FontSize',10,'FontName','Helvetica');
ylim([0.6 1.05]);
ylabel('Change in growth','FontSize',12,'FontName','Helvetica','Color','k');
xtickangle(90);
box off;








