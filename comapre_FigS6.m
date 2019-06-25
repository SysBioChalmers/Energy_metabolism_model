% Compare E. coli and yeast
% Related to Fig S6

%% Import data
[num, ~, ~] = xlsread('Comparison.xlsx','ecoli_raw');
mw_ecoli = num(:,1);
kcat_ecoli = num(:,2);
pc_ecoli = num(:,3);
clear num;

[num, ~, ~] = xlsread('Comparison.xlsx','yeast_raw');
mw_yeast = num(:,1);
kcat_yeast = num(:,2);
pc_yeast = num(:,3);
clear num;

[num, txt, ~] = xlsread('Comparison.xlsx','comparison_rxn');
abs_ecoli = num(:,1);
abs_yeast = num(:,2);
fc = num(:,3);
rxnID = txt(2:end,1);
clear num txt;

[num, txt, ~] = xlsread('Comparison.xlsx','comparison_pathway');
pw_ecoli = num(:,1);
pw_yeast = num(:,2);
fcpw = num(:,3);
pwID = txt(2:end,1);
clear num txt;

%% Figures

figure('Name','mw');
hold on;
box on;
[f,x] = ecdf(mw_ecoli);
plot(x,f,'LineWidth',1.5,'Color',[197,27,138]/256);
clear x f;
[f,x] = ecdf(mw_yeast);
plot(x,f,'LineWidth',1.5,'Color',[49,163,84]/256);
clear x f;
xlim([0 1000]);
legend({'E. coli','S. cerevisiae',},'FontSize',12,'FontName','Helvetica','location','se');
legend('boxoff');
xlabel('Molecular weight (kDa)','FontSize',14,'FontName','Helvetica');
ylabel('Cumulative frequency','FontSize',14,'FontName','Helvetica');
set(gcf,'position',[100 300 270 150]);
set(gca,'position',[0.2 0.22 0.76 0.71]);

figure('Name','kcat');
[f,x] = ecdf(kcat_ecoli);
semilogx(x,f,'LineWidth',1.5,'Color',[197,27,138]/256);
clear x f;
hold on;
[f,x] = ecdf(kcat_yeast);
semilogx(x,f,'LineWidth',1.5,'Color',[49,163,84]/256);
clear x f;
xlim([10 10000]);
legend({'E. coli','S. cerevisiae',},'FontSize',12,'FontName','Helvetica','location','se');
legend('boxoff');
xlabel('kcat (/s)','FontSize',14,'FontName','Helvetica');
ylabel('Cumulative frequency','FontSize',14,'FontName','Helvetica');
set(gcf,'position',[300 300 270 150]);
set(gca,'position',[0.2 0.22 0.76 0.71]);

figure('Name','pc');
hold on;
box on;
[f,x] = ecdf(pc_ecoli);
plot(x,f,'LineWidth',1.5,'Color',[197,27,138]/256);
clear x f;
[f,x] = ecdf(pc_yeast);
plot(x,f,'LineWidth',1.5,'Color',[49,163,84]/256);
clear x f;
legend({'E. coli','S. cerevisiae',},'FontSize',12,'FontName','Helvetica','location','se');
legend('boxoff');
xlabel('MW/kcat','FontSize',14,'FontName','Helvetica');
ylabel('Cumulative frequency','FontSize',14,'FontName','Helvetica');
set(gcf,'position',[500 300 270 150]);
set(gca,'position',[0.2 0.22 0.76 0.71]);

figure('Name','hm_abs');
cdata = [abs_ecoli abs_yeast];
xvalues = {'Eco','Sce'};
yvalues = rxnID;
map = colormap(heatmap(xvalues,yvalues,cdata,'Colormap',gray));
map = sort(map,'descend');
h = heatmap(xvalues,yvalues,cdata,'Colormap',map,'CellLabelColor','none');
h.Title = 'Protein cost';
set(gca,'FontSize',12,'FontName','Helvetica');
set(gcf,'position',[600 300 100 600]);
set(gca,'position',[0.2 0.1 0.4 0.75]);

figure('Name','hm_fc');
cdata = log2(fc);
maxc = [208,28,139];
minc = [49,163,84];
step_max = round(abs(max(cdata)),2)*100;
itv_max = ([255,255,255] - maxc)/step_max;
step_min = round(abs(min(cdata)),2)*100;
itv_min = (minc - [255,255,255])/step_min;
color_map_max = transpose([maxc(1):itv_max(1):255;maxc(2):itv_max(2):255;maxc(3):itv_max(3):255;]);
color_map_min = transpose([255:itv_min(1):minc(1);255:itv_min(2):minc(2);255:itv_min(3):minc(3);]);
color_map = [sort(color_map_min,'ascend');sort(color_map_max,'descend')]/255;
h = heatmap(log2(fc),'Colormap',color_map,'ColorMethod','count','CellLabelColor','none');
h.Title = 'log2(FC)';
set(gca,'FontSize',12,'FontName','Helvetica');
set(gcf,'position',[200 300 100 600]);
set(gca,'position',[0.2 0.1 0.4 0.75]);
