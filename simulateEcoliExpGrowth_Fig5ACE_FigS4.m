% Various E. coli strains under exponential growth 
% Related to Fig 5ACE
% Related to Fig S4A

% Data from (PMID: 25304508 27135538)

% cd ../../../GitHub/cobratoolbox;
% initCobraToolbox;

addpath('Functions/');

model_ecoli = xls2model('Model_ecoli.xlsx');

prot_cost_info = struct();
[num, txt, ~] = xlsread('Model_ecoli.xlsx','Protein_cost_info');
prot_cost_info.id = txt(2:end,1);
prot_cost_info.value = num;
clear num txt;

%% Modify model

% Remove GAM
model = model_ecoli;
model.S(contains(model.mets,'h2o[c]'),contains(model.rxns,'Biomass')) = 0;
model.S(contains(model.mets,'atp[c]'),contains(model.rxns,'Biomass')) = 0;
model.S(contains(model.mets,'pi[c]'),contains(model.rxns,'Biomass')) = 0;
model.S(contains(model.mets,'adp[c]'),contains(model.rxns,'Biomass')) = 0;
model.S(contains(model.mets,'h[c]'),contains(model.rxns,'Biomass')) = 0;

% Unlimit NGAM
model = changeRxnBounds(model, 'ATPM', 0, 'l');
model = changeRxnBounds(model, 'ATPM', 1000, 'u');


%% Wildtype
model_wt = changeRxnBounds(model, 'EXbiomass', 0.69, 'b');
model_wt = changeRxnBounds(model_wt, 'EXglc', -8.59, 'b');
model_wt = changeRxnBounds(model_wt, 'EXac', 3.91, 'b');

[minProt_wt, fluxes] = minimizeProtein(model_wt, prot_cost_info, 'EXbiomass', 'max',0.5883);
mu = fluxes(strcmp(model_wt.rxns,'EXbiomass'));
totGlc = abs(fluxes(strcmp(model_wt.rxns,'EXglc')));
HY = fluxes(strcmp(model_wt.rxns,'GLCtex_HY'))/totGlc;
LY = fluxes(strcmp(model_wt.rxns,'GLCtex_LY'))/totGlc;
Bio = 1 - HY - LY;
ATP = fluxes(strcmp(model_wt.rxns,'ATPM'));
[HYp, LYp, Biop] = calculateProtein(model_wt, prot_cost_info, fluxes, 0.5883);
[Glcp, TCAp, OXPp, Prdp] = calculateProteome(model_wt, prot_cost_info, fluxes, 'Ecoli', 0.5883);
wt = [minProt_wt; mu; HY; LY; Bio; ATP; HYp; LYp; Biop; Glcp; TCAp; OXPp; Prdp];

clear minProt_wt mu HY LY Bio fluxes model_wt totGlc ATP Glcp TCAp OXPp Prdp;

%% Evoloved strains
data = [0.98 -13.51 8.43;
        0.96 -12.19 7.89;
        0.93 -12.77 7.11;
        1.01 -13.13 5.12;
        0.97 -11.01 3.97;
        0.92 -10.43 2.36;
        0.89 -12.59 5.05;
        0.92 -13.13 6.99;
        0.95 -13.98 9.27;
        0.85 -10.23 4.75; % rpoBE546V (PMID: 27135538)
        0.88 -10.44 4.87];% rpoBE762K (PMID: 27135538)
    %   mu   q_glc  q_ac
strainid = {'3';'4';'6';'7';'7A';'7B';'8';'9';'10';'MT1';'MT2'};

evolved = zeros(13,11);
for i = 1:length(strainid)
    model_evo = changeRxnBounds(model, 'EXbiomass', data(i,1), 'b');
    model_evo = changeRxnBounds(model_evo, 'EXglc', data(i,2), 'b');
    model_evo = changeRxnBounds(model_evo, 'EXac', data(i,3), 'b');
    
    [minProt, fluxes] = minimizeProtein(model_evo, prot_cost_info, 'EXbiomass', 'max',0.5883);
    
    mu = fluxes(strcmp(model_evo.rxns,'EXbiomass'));
    totGlc = abs(fluxes(strcmp(model_evo.rxns,'EXglc')));
    HY = fluxes(strcmp(model_evo.rxns,'GLCtex_HY'))/totGlc;
    LY = fluxes(strcmp(model_evo.rxns,'GLCtex_LY'))/totGlc;
    Bio = 1 - HY - LY;
    ATP = fluxes(strcmp(model_evo.rxns,'ATPM'));
    [HYp, LYp, Biop] = calculateProtein(model_evo, prot_cost_info, fluxes,0.5883);
    [Glcp, TCAp, OXPp, Prdp] = calculateProteome(model_evo, prot_cost_info, fluxes, 'Ecoli', 0.5883);
    evolved(:,i) = [minProt; mu; HY; LY; Bio; ATP; HYp; LYp; Biop; Glcp; TCAp; OXPp; Prdp];
end

clear i minProt data mu HY LY Bio fluxes model_evo totGlc ATP Glcp TCAp OXPp Prdp;

%% Figure

res = [wt evolved];

% Flux distribution among HY, LY and Biomass pathways
figure('Name','1');

hold on;
yyaxis left;
c = categorical(transpose(['WT';strainid]));
c = reordercats(c,transpose(['WT';strainid]));
b = bar(c,transpose(res(3:5,:)),'stacked');
b(1).FaceColor = [43,131,186]/256;
b(1).FaceAlpha = 0.6;
b(2).FaceColor = [215,25,28]/256;
b(2).FaceAlpha = 0.6;
b(3).FaceColor = [146,39,143]/256;
b(3).FaceAlpha = 0.6;
ylim([0 1]);
set(gca,'FontSize',10,'FontName','Helvetica');
set(gca,'ycolor','k');
ylabel('Fraction in total glucose fluxes','FontSize',12,'FontName','Helvetica','Color','k');

yyaxis right;
plot(c,res(2,:),'ko','LineWidth',1,'MarkerSize',6);
set(gca,'ycolor','k');
ylim([0.66 1.04]);
ylabel('Growth rate (/h)','FontSize',12,'FontName','Helvetica','Color','k');

legend({'High-yield flux',...
        'Low-yield flux',...
        'Biomass flux',...
        'Growth rate'},'FontSize',9,'FontName','Helvetica','location','nw');

set(gcf,'position',[0 0 250 220]);
set(gca,'position',[0.15 0.13 0.65 0.85]);

clear b Biop c HYp LYp;

figure('Name','1_s');

hold on;
yyaxis left;
c = categorical(transpose(['WT';strainid]));
c = reordercats(c,transpose(['WT';strainid]));
b = bar(c,transpose(res(7:9,:)./res(1,:)),'stacked');
b(1).FaceColor = [43,131,186]/256;
b(1).FaceAlpha = 0.6;
b(2).FaceColor = [215,25,28]/256;
b(2).FaceAlpha = 0.6;
b(3).FaceColor = [146,39,143]/256;
b(3).FaceAlpha = 0.6;
ylim([0 1]);
set(gca,'FontSize',10,'FontName','Helvetica');
set(gca,'ycolor','k');
ylabel('Fraction in proteome','FontSize',12,'FontName','Helvetica','Color','k');

yyaxis right;
plot(c,res(2,:),'ko','LineWidth',1,'MarkerSize',6);
set(gca,'ycolor','k');
ylim([0.66 1.04]);
ylabel('Growth rate (/h)','FontSize',12,'FontName','Helvetica','Color','k');

legend({'High-yield',...
        'Low-yield',...
        'Biomass',...
        'Growth rate'},'FontSize',9,'FontName','Helvetica','location','nw');

set(gcf,'position',[10 0 250 220]);
set(gca,'position',[0.15 0.13 0.65 0.85]);

clear b Biop c HYp LYp;

figure('Name','1_ss');

hold on;
yyaxis left;
c = categorical(transpose(['WT';strainid]));
c = reordercats(c,transpose(['WT';strainid]));
b = bar(c,transpose(res(10:13,:)./res(1,:)),'stacked');
b(1).FaceColor = [27,158,119]/256;
b(1).FaceAlpha = 0.6;
b(2).FaceColor = [217,95,2]/256;
b(2).FaceAlpha = 0.6;
b(3).FaceColor = [117,112,179]/256;
b(3).FaceAlpha = 0.6;
b(4).FaceColor = [231,41,138]/256;
b(4).FaceAlpha = 0.6;
ylim([0 1]);
set(gca,'FontSize',10,'FontName','Helvetica');
set(gca,'ycolor','k');
ylabel('Fraction in proteome','FontSize',12,'FontName','Helvetica','Color','k');

yyaxis right;
plot(c,res(2,:),'ko','LineWidth',1,'MarkerSize',6);
set(gca,'ycolor','k');
ylim([0.66 1.04]);
ylabel('Growth rate (/h)','FontSize',12,'FontName','Helvetica','Color','k');

legend({'Glycolysis',...
        'TCA cycle',...
        'OXPHOS',...
        'Acetate',...
        'Growth rate'},'FontSize',9,'FontName','Helvetica','location','nw');

set(gcf,'position',[20 0 250 220]);
set(gca,'position',[0.15 0.13 0.65 0.85]);

clear b Biop c HYp LYp;


% Growth rate and protein requirement
figure('Name','2');
hold on;
box on;
x = res(2,:);
y = res(1,:);
scatter(x,y,50,'o','k');
ylim([0.03 0.1]);
set(gca,'FontSize',10,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',11,'FontName','Helvetica');
ylabel(['Protein allocation',char(13,10)','(g/gCDW)'],'FontSize',11,'FontName','Helvetica');

p = polyfit(x,y,1);
p_y = polyval(p,x);
plot(x,p_y,'k-','LineWidth',1.5);

R = corrcoef(x,y);
corr_coef = ['R^2 = ',num2str(round(R(1,2)^2,3))];
text(0.65,0.13,corr_coef,'FontSize',13,'FontName','Helvetica','Color','k');

xlim([0.6 1.1]);
ylim([0.07 0.2]);

set(gcf,'position',[200 0 230 130]);
set(gca,'position',[0.23 0.25 0.72 0.7]);

clear x y corr_coef p p_y R;

% Growth rate and ATP production rate
figure('Name','3');
hold on;
box on;
x = res(2,:);
y = res(6,:);
scatter(x,y,50,'o','k');
ylim([40 120]);
set(gca,'FontSize',10,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',11,'FontName','Helvetica');
ylabel(['ATP production rate',char(13,10)','(mmol/gCDW/h)'],'FontSize',11,'FontName','Helvetica');

p = polyfit(x,y,1);
p_y = polyval(p,x);
plot(x,p_y,'k-','LineWidth',1.5);

R = corrcoef(x,y);
corr_coef = ['R^2 = ',num2str(round(R(1,2)^2,3))];
text(0.65,100,corr_coef,'FontSize',13,'FontName','Helvetica','Color','k');

xlim([0.6 1.1]);

set(gcf,'position',[400 0 230 130]);
set(gca,'position',[0.23 0.25 0.72 0.7]);


clear x y corr_coef p p_y R;

% ATP production rate and protein requirements
figure('Name','4');
hold on;
box on;

x = res(1,:);
y = res(6,:);
s = scatter(x,y,80,res(2,:),'filled');
s.MarkerEdgeColor = 'k';
h = colorbar;
colormap Gray;
set(get(h,'label'),'string','Growth rate (/h)');

set(gca,'FontSize',10,'FontName','Helvetica');
xlabel('Protein allocation (g/gCDW)','FontSize',12,'FontName','Helvetica');
ylabel('ATP production rate (mmol/gCDW/h)','FontSize',12,'FontName','Helvetica');

xlim([0.07 0.19]);
ylim([42 102]);

p = polyfit(res(1,:),res(6,:),1);
p_y = polyval(p,res(1,:));
plot(res(1,:),p_y,'k-','LineWidth',1.5);

R = corrcoef(res(1,:),res(6,:));
corr_coef = ['R^2 = ',num2str(round(R(1,2)^2,3))];
text(0.14,75,corr_coef,'FontSize',13,'FontName','Helvetica','Color','k');

set(gcf,'position',[600 400 310 250]);
set(gca,'position',[0.13 0.15 0.63 0.8]);

% Change heatmap color: Edit -> Figure properties -> Colormap -> Custom
% Lowest value: [255 255 255]   Highest value: [142,1,82]

clear x y s h;

clear corr_coef p p_y R;


% Protein allocation and ATP requirement for biomass formation
figure('Name','5');
hold on;
box on;
x = res(1,:);
y = res(6,:) ./ res(2,:);
s = scatter(x,y,80,res(2,:),'filled');
s.MarkerEdgeColor = 'k';
h = colorbar;
colormap Gray;
set(get(h,'label'),'string','Growth rate (/h)');

set(gca,'FontSize',10,'FontName','Helvetica');
xlabel('Protein allocation (g/gCDW)','FontSize',12,'FontName','Helvetica');
ylabel('ATP demand (mmol/gCDW)','FontSize',12,'FontName','Helvetica');

xlim([0.07 0.19]);
ylim([47 113]);

set(gcf,'position',[800 400 310 250]);
set(gca,'position',[0.13 0.15 0.63 0.8]);

% Change heatmap color: Edit -> Figure properties -> Colormap -> Custom
% Lowest value: [255 255 255]   Highest value: [142,1,82]

clear x y s h;

%% Supplementary figures
% % ATP production rate and protein requirements
% figure('Name','4s');
% hold on;
% box on;
% 
% model_hy = changeRxnBounds(model, 'GLCtex_LY', 0, 'b');
% glc_data = -2:-2:-10;
% hy = zeros(2,length(glc_data));
% for i = 1:length(glc_data)
%     glc_in = glc_data(i);
%     model_hy = changeRxnBounds(model_hy, 'EXglc', glc_in, 'b');
%     model_hy = changeObjective(model_hy, 'ATPM');
%     sol = optimizeCbModel(model_hy, 'max', 'one');
%     [HYp, ~, ~] = calculateProtein(model, prot_cost_info, sol.x,0.5883);
%     hy(1,i) = HYp;
%     hy(2,i) = sol.f;
% end
% 
% model_ly = changeRxnBounds(model, 'GLCtex_HY', 0, 'b');
% glc_data = -4:-4:-20;
% ly = zeros(2,length(glc_data));
% for i = 1:length(glc_data)
%     glc_in = glc_data(i);
%     model_ly = changeRxnBounds(model_ly, 'EXglc', glc_in, 'b');
%     model_ly = changeObjective(model_ly, 'ATPM');
%     sol = optimizeCbModel(model_ly, 'max', 'one');
%     [~, LYp, ~] = calculateProtein(model, prot_cost_info, sol.x,0.5883);
%     ly(1,i) = LYp;
%     ly(2,i) = sol.f;
% end
% 
% clear i Biop HYp LYp glc_data glc_in model_hy model_ly sol;
% 
% plot([0, hy(1,:)],[0, hy(2,:)],'-','LineWidth',2,'Color',[43,131,186]/255);
% plot([0, ly(1,:)],[0, ly(2,:)],'-','LineWidth',2,'Color',[215,25,28]/255);
% 
% x = res(7,:) + res(8,:);
% y = res(6,:);
% s = scatter(x,y,80,res(2,:),'filled');
% s.MarkerEdgeColor = 'k';
% h = colorbar;
% colormap Gray;
% set(get(h,'label'),'string','Growth rate (/h)');
% 
% set(gca,'FontSize',10,'FontName','Helvetica');
% xlabel('Protein allocation (g/gCDW)','FontSize',12,'FontName','Helvetica');
% ylabel('ATP production rate (mmol/gCDW/h)','FontSize',12,'FontName','Helvetica');
% 
% % xlim([0 0.095]);
% % ylim([0 105]);
% 
% p = polyfit(res(7,:)+res(8,:),res(6,:),1);
% p_y = polyval(p,res(7,:)+res(8,:));
% plot(res(7,:)+res(8,:),p_y,'k-','LineWidth',1.5);
% 
% R = corrcoef(res(7,:)+res(8,:),res(6,:));
% corr_coef = ['R^2 = ',num2str(round(R(1,2)^2,3))];
% text(0.065,60,corr_coef,'FontSize',13,'FontName','Helvetica','Color','k');
% 
% set(gcf,'position',[600 0 310 250]);
% set(gca,'position',[0.13 0.15 0.63 0.8]);
% 
% % Change heatmap color: Edit -> Figure properties -> Colormap -> Custom
% % Lowest value: [255 255 255]   Highest value: [146 39 143]
% 
% clear x y s h;
% 
% clear corr_coef hy ly p p_y R wt evolved;
% 
% % % heatmap
% % imagesc(res(2,:));
% 
% % Protein allocation and ATP requirement for biomass formation
% figure('Name','5s');
% hold on;
% box on;
% x = res(7,:) + res(8,:);
% y = res(6,:) ./ res(2,:);
% s = scatter(x,y,80,res(2,:),'filled');
% s.MarkerEdgeColor = 'k';
% h = colorbar;
% colormap Gray;
% set(get(h,'label'),'string','Growth rate (/h)');
% 
% set(gca,'FontSize',10,'FontName','Helvetica');
% xlabel('Protein allocation (g/gCDW)','FontSize',12,'FontName','Helvetica');
% ylabel('ATP demand (mmol/gCDW)','FontSize',12,'FontName','Helvetica');
% 
% % xlim([0.02 0.1]);
% % ylim([45 115]);
% 
% set(gcf,'position',[800 0 310 250]);
% set(gca,'position',[0.13 0.15 0.63 0.8]);
% 
% % Change heatmap color: Edit -> Figure properties -> Colormap -> Custom
% % Lowest value: [255 255 255]   Highest value: [146 39 143]
% 
% clear x y s h;






