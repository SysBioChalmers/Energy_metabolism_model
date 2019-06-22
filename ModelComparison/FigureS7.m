figure();
set(gcf,'position',[450 50 450 650]);

color_ecoli = [142,1,82]/255;
color_yeast = [39,100,25]/255;

% The rxns with assumed kcats will be not displayed in the figures
remove_id_ecoli = {'ACt2rpp'};
remove_id_yeast = {'PYRt2m';'AKGDam';'AKGDbm';'GCC2cm';'SUCOASm';'ATPtmH';'PIt2m';'H2Ot'};

% Ecoli protein cost = 0.5
[value, id, ~] = xlsread('Benchmarking_result.xlsx','Ecoli0.5');
idx1 = ~ismember(id,remove_id_ecoli);
id = id(idx1);
value = value(idx1);
value = round(value,2);
[value, idx2] = sort(value,'descend');
id = id(idx2);

subplot(4,1,1);
h1 = bar(1:length(id),value,'FaceColor',color_ecoli,'FaceAlpha',0.3,'EdgeColor',color_ecoli,'LineWidth',0.5);
set(gca,'XTick',1:1:length(id));
set(gca,'XTickLabel',id);
set(gca,'FontSize',10,'FontName','Helvetica');
title('2-fold increase in activity (MOMENT-Ecoli)','FontSize',12,'FontName','Helvetica','Color','k');
ylim([0.9 1.1]);
ylabel('Change in growth','FontSize',12,'FontName','Helvetica','Color','k');
xtickangle(90);
box off;

% Yeast protein cost = 0.5
[value, id, ~] = xlsread('Benchmarking_result.xlsx','Yeast0.5');
idx1 = ~ismember(id,remove_id_yeast);
id = id(idx1);
value = value(idx1);
value = round(value,2);
[value, idx2] = sort(value,'descend');
id = id(idx2);

subplot(4,1,2);
h2 = bar(1:length(id),value,'FaceColor',color_yeast,'FaceAlpha',0.3,'EdgeColor',color_yeast,'LineWidth',0.5);
set(gca,'XTick',1:1:length(id));
set(gca,'XTickLabel',id);
set(gca,'FontSize',10,'FontName','Helvetica');
title('2-fold increase in activity (ec-Yeast)','FontSize',12,'FontName','Helvetica','Color','k');
ylim([0.9 1.1]);
ylabel('Change in growth','FontSize',12,'FontName','Helvetica','Color','k');
xtickangle(90);
box off;

% Ecoli protein cost = 0.1
[value, id, ~] = xlsread('Benchmarking_result.xlsx','Ecoli0.1');
idx1 = ~ismember(id,remove_id_ecoli);
id = id(idx1);
value = value(idx1);
value = round(value,2);
[value, idx2] = sort(value,'descend');
id = id(idx2);

subplot(4,1,3);
h3 = bar(1:length(id),value,'FaceColor',color_ecoli,'FaceAlpha',0.3,'EdgeColor',color_ecoli,'LineWidth',0.5);
set(gca,'XTick',1:1:length(id));
set(gca,'XTickLabel',id);
set(gca,'FontSize',10,'FontName','Helvetica');
title('10-fold increase in activity (MOMENT-Ecoli)','FontSize',12,'FontName','Helvetica','Color','k');
ylim([0.9 1.2]);
ylabel('Change in growth','FontSize',12,'FontName','Helvetica','Color','k');
xtickangle(90);
box off;

% Yeast protein cost = 0.1
[value, id, ~] = xlsread('Benchmarking_result.xlsx','Yeast0.1');
idx1 = ~ismember(id,remove_id_yeast);
id = id(idx1);
value = value(idx1);
value = round(value,2);
[value, idx2] = sort(value,'descend');
id = id(idx2);

subplot(4,1,4);
h4 = bar(1:length(id),value,'FaceColor',color_yeast,'FaceAlpha',0.3,'EdgeColor',color_yeast,'LineWidth',0.5);
set(gca,'XTick',1:1:length(id));
set(gca,'XTickLabel',id);
set(gca,'FontSize',10,'FontName','Helvetica');
title('10-fold increase in activity (ec-Yeast)','FontSize',12,'FontName','Helvetica','Color','k');
ylim([0.9 1.2]);
ylabel('Change in growth','FontSize',12,'FontName','Helvetica','Color','k');
xtickangle(90);
box off;