% Proteomics data analysis.
% Related to Fig S3

[num_yeast, ~, ~] = xlsread('Protein_abundance.xlsx','Yeast');
[num_ecoli, ~, ~] = xlsread('Protein_abundance.xlsx','Ecoli');


%% Different conditions

% Add Violinplot path
% addpath(genpath('../Violinplot-Matlab'));

% Yeast dataset 1 (PMID: 28365149)

ref = num_yeast(3:120,3);
others = num_yeast(3:120,4:12);
result = log2(others./ref);

tot_ref = sum(ref);
tot_others = sum(others);
tot_result = log2(tot_others./tot_ref);

figure('Name','1');
hold on;

plot([1;length(tot_result)],[0;0],'k-.','LineWidth',1.5);
% violins = violinplot(result,'','ViolinAlpha',0.2,'ViolinColor',[240,2,127]/255);
h = boxplot(result,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[240,2,127]/255);
set(h,{'linew'},{1});
scatter(1:length(tot_result),tot_result,40,'k','o','filled','MarkerFaceAlpha',0.7);


set(gca,'FontSize',10,'FontName','Helvetica');
ylabel('log2 Foldchange','FontSize',12,'FontName','Helvetica');

clear h others ref result tot_others tot_ref tot_result;

set(gcf,'position',[0 0 170 190]);
set(gca,'position',[0.19 0.1 0.78 0.85]);

% E.coli dataset 1 (PMID: 25712329)

ref = num_ecoli(3:75,3);
others = num_ecoli(3:75,4:24);
result = log2(others./ref);

tot_ref = sum(ref);
tot_others = sum(others);
tot_result = log2(tot_others./tot_ref);

figure('Name','2');
hold on;

plot([1;length(tot_result)],[0;0],'k-.','LineWidth',1.5);
% violins = violinplot(result,'','ViolinAlpha',0.2,'ViolinColor',[56,108,176]/255);
h = boxplot(result,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[56,108,176]/255);
set(h,{'linew'},{1});
scatter(1:length(tot_result),tot_result,40,'k','o','filled','MarkerFaceAlpha',0.7);


set(gca,'FontSize',10,'FontName','Helvetica');
ylabel('log2 Foldchange','FontSize',12,'FontName','Helvetica');

clear h others ref result tot_others tot_ref tot_result;

set(gcf,'position',[200 0 350 190]);
set(gca,'position',[0.11 0.1 0.87 0.85]);

% E.coli dataset 2 (PMID: 26641532)

ref = num_ecoli(3:76,25);
others = num_ecoli(3:76,26:46);
result = log2(others./ref);

tot_ref = sum(ref);
tot_others = sum(others);
tot_result = log2(tot_others./tot_ref);

figure('Name','3');
hold on;

plot([1;length(tot_result)],[0;0],'k-.','LineWidth',1.5);
% violins = violinplot(result,'','ViolinAlpha',0.2,'ViolinColor',[102,166,30]/255);
h = boxplot(result,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[102,166,30]/255);
set(h,{'linew'},{1});
scatter(1:length(tot_result),tot_result,40,'k','o','filled','MarkerFaceAlpha',0.7);


set(gca,'FontSize',10,'FontName','Helvetica');
ylabel('log2 Foldchange','FontSize',12,'FontName','Helvetica');

clear h others ref result tot_others tot_ref tot_result;

set(gcf,'position',[500 0 350 190]);
set(gca,'position',[0.11 0.1 0.87 0.85]);

clear num_ecoli num_yeast;
