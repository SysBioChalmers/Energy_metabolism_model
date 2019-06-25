% Various yeast strains under exponential growth 
% Related to Fig 5 and Fig S5

%% Experimental data (PMID: 27229477/25278608)

fluxes_exp = [0.29   0.12   0.22   0.22   0.22   0.22   0.23 % growth rate
              11.7   11.9   16.1   15.9   16.3   16.2   15.8 % glucose
              16.7   17.1   21.8   21.5   22.1   21.9   21.9 % ethanol
              0.8    2.8    3.4    2.9    3.6    3.8    3.5];% glycerol

strainid = {'WT(30)' 'WT' 'TT11' 'TT21' 'TT22' 'TT31' 'TT33'};


%% Simulations
model_yeast = xls2model('Model_yeast.xlsx');

prot_cost_yeast = struct();
[num, txt, ~] = xlsread('Model_yeast.xlsx','Protein_cost_info');
prot_cost_yeast.id = txt(2:end,1);
prot_cost_yeast.value = num;
clear num txt;

model = model_yeast;
model.S(contains(model.mets,'h2o[c]'),contains(model.rxns,'Biomass')) = 0;
model.S(contains(model.mets,'atp[c]'),contains(model.rxns,'Biomass')) = 0;
model.S(contains(model.mets,'pi[c]'),contains(model.rxns,'Biomass')) = 0;
model.S(contains(model.mets,'adp[c]'),contains(model.rxns,'Biomass')) = 0;
model.S(contains(model.mets,'h[c]'),contains(model.rxns,'Biomass')) = 0;
model = changeRxnBounds(model, 'ATPM', 0, 'l');
model = changeRxnBounds(model, 'ATPM', 1000, 'u');

res = zeros(12,length(fluxes_exp));

for i = 1:length(fluxes_exp)
    mu = fluxes_exp(1,i);
    glc = fluxes_exp(2,i);
    eth = fluxes_exp(3,i);
    glyc = fluxes_exp(4,i);
    model = changeRxnBounds(model, 'EXbiomass', mu, 'b');
    model = changeRxnBounds(model, 'EXglc', -glc+glyc/1.64, 'b');
    model = changeRxnBounds(model, 'EXetoh', eth, 'b');
    [minProt, fluxes] = minimizeProtein(model, prot_cost_yeast, 'ATPM', 'max',0.6914);
    minProt_adj = minProt + 0.005 * glyc/1.64;
    if ~isempty(fluxes)
        glc_new = -fluxes(strcmp(model.rxns,'EXglc'));
        HY = fluxes(strcmp(model.rxns,'GLCt1_HY'))/(glc_new+glyc/1.64);
        LY = fluxes(strcmp(model.rxns,'GLCt1_LY'))/(glc_new+glyc/1.64);
        Bio = fluxes(strcmp(model.rxns,'GLCt1_Bio'))/(glc_new+glyc/1.64);
        ATPM = fluxes(strcmp(model.rxns,'ATPM'));
        [HYp, LYp, Biop] = calculateProtein(model, prot_cost_yeast, fluxes,0.6914);
        [Glcp, TCAp, OXPp, Prdp] = calculateProteome(model, prot_cost_yeast, fluxes, 'Yeast', 0.6914);
        res(:,i) = [HY; LY; Bio; HYp; LYp; Biop; ATPM; minProt_adj; Glcp; TCAp; OXPp; Prdp];
    else
        [minProt, fluxes] = obtainFeasibleSol(model,prot_cost_yeast,'Yeast',0.01,'ATPM','max',0.6914);
        minProt_adj = minProt + 0.005 * glyc/1.64;
        glc_new = -fluxes(strcmp(model.rxns,'EXglc'));
        HY = fluxes(strcmp(model.rxns,'GLCt1_HY'))/(glc_new+glyc/1.64);
        LY = fluxes(strcmp(model.rxns,'GLCt1_LY'))/(glc_new+glyc/1.64);
        Bio = fluxes(strcmp(model.rxns,'GLCt1_Bio'))/(glc_new+glyc/1.64);
        ATPM = fluxes(strcmp(model.rxns,'ATPM'));
        [HYp, LYp, Biop] = calculateProtein(model, prot_cost_yeast, fluxes,0.6914);
        [Glcp, TCAp, OXPp, Prdp] = calculateProteome(model, prot_cost_yeast, fluxes, 'Yeast', 0.6914);
        res(:,i) = [HY; LY; Bio; HYp; LYp; Biop; ATPM; minProt_adj; Glcp; TCAp; OXPp; Prdp];
    end
    
end

clear mu glc eth glyc model fluxes HY LY Bio ATPM glc_new;
clear HYp LYp Biop i minProt minProt_adj Glcp TCAp OXPp Prdp;


%% Figure

figure('Name','1');

hold on;
yyaxis left;
c = categorical(strainid);
c = reordercats(c,strainid);
b = bar(c,transpose(res(1:3,:)),'stacked');
b(1).FaceColor = [43,131,186]/256;
b(1).FaceAlpha = 0.6;
b(2).FaceColor = [215,25,28]/256;
b(2).FaceAlpha = 0.6;
b(3).FaceColor = [146,39,143]/256;
b(3).FaceAlpha = 0.6;
ylim([0 1]);
set(gca,'FontSize',12,'FontName','Helvetica');
set(gca,'ycolor','k');
ylabel('Fraction in total glucose fluxes','FontSize',14,'FontName','Helvetica','Color','k');

yyaxis right;
plot(c,fluxes_exp(1,:),'ko','LineWidth',1,'MarkerSize',5);
set(gca,'ycolor','k');
ylim([0.1 0.3]);

ylabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica','Color','k');

% legend({'High-yield flux',...
%         'Low-yield flux',...
%         'Biomass flux',...
%         'Growth rate'},'FontSize',9,'FontName','Helvetica','location','nw');

set(gcf,'position',[0 0 250 220]);
set(gca,'position',[0.15 0.13 0.65 0.85]);

clear b c;


figure('Name','1_s');

hold on;
yyaxis left;
c = categorical(strainid);
c = reordercats(c,strainid);
b = bar(c,transpose(res(4:6,:)./res(8,:)),'stacked');
b(1).FaceColor = [43,131,186]/256;
b(1).FaceAlpha = 0.6;
b(2).FaceColor = [215,25,28]/256;
b(2).FaceAlpha = 0.6;
b(3).FaceColor = [146,39,143]/256;
b(3).FaceAlpha = 0.6;
ylim([0 1]);
set(gca,'FontSize',12,'FontName','Helvetica');
set(gca,'ycolor','k');
ylabel('Fraction in proteome','FontSize',14,'FontName','Helvetica','Color','k');

yyaxis right;
plot(c,fluxes_exp(1,:),'ko','LineWidth',1,'MarkerSize',5);
set(gca,'ycolor','k');
ylim([0.1 0.3]);

ylabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica','Color','k');

legend({'High-yield',...
        'Low-yield',...
        'Biomass',...
        'Growth rate'},'FontSize',9,'FontName','Helvetica','location','nw');

set(gcf,'position',[10 0 250 220]);
set(gca,'position',[0.15 0.13 0.65 0.85]);

clear b c;

figure('Name','1_ss');

hold on;
yyaxis left;
c = categorical(strainid);
c = reordercats(c,strainid);
b = bar(c,transpose(res(9:12,:)./res(8,:)),'stacked');
b(1).FaceColor = [27,158,119]/256;
b(1).FaceAlpha = 0.6;
b(2).FaceColor = [217,95,2]/256;
b(2).FaceAlpha = 0.6;
b(3).FaceColor = [117,112,179]/256;
b(3).FaceAlpha = 0.6;
b(4).FaceColor = [231,41,138]/256;
b(4).FaceAlpha = 0.6;
ylim([0 1]);
set(gca,'FontSize',12,'FontName','Helvetica');
set(gca,'ycolor','k');
ylabel('Fraction in proteome','FontSize',14,'FontName','Helvetica','Color','k');

yyaxis right;
plot(c,fluxes_exp(1,:),'ko','LineWidth',1,'MarkerSize',5);
set(gca,'ycolor','k');
ylim([0.1 0.3]);

ylabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica','Color','k');

legend({'Glycolysis',...
        'TCA cycle',...
        'OXPHOS',...
        'Ethanol',...
        'Growth rate'},'FontSize',9,'FontName','Helvetica','location','nw');

set(gcf,'position',[20 0 250 220]);
set(gca,'position',[0.15 0.13 0.65 0.85]);

clear b c;


% Protein allocation and ATP production rate
figure('Name','2');
hold on;
box on;

x = res(8,:);
y = res(7,:);
s = scatter(x,y,80,fluxes_exp(1,:),'filled');
s.MarkerEdgeColor = 'k';

h = colorbar;
colormap Gray;
set(get(h,'label'),'string','Growth rate (/h)');

set(gca,'FontSize',12,'FontName','Helvetica');
xlabel('Protein allocation (g/gCDW)','FontSize',14,'FontName','Helvetica');
ylabel('ATP production rate (mmol/gCDW/h)','FontSize',14,'FontName','Helvetica');

xlim([0.08 0.22]);
ylim([26 61]);

p = polyfit(res(8,:),res(7,:),1);
p_y = polyval(p,res(8,:));
plot(res(8,:),p_y,'k-','LineWidth',1.5);

R = corrcoef(res(8,:),res(7,:));
corr_coef = ['R^2 = ',num2str(round(R(1,2)^2,3))];
text(0.15,43,corr_coef,'FontSize',14,'FontName','Helvetica','Color','k');

set(gcf,'position',[600 400 310 250]);
set(gca,'position',[0.13 0.15 0.63 0.8]);

% Protein allocation and ATP requirement for biomass formation
figure('Name','3');
hold on;
box on;
x = res(8,:);
y = res(7,:) ./ fluxes_exp(1,:);
s = scatter(x,y,80,fluxes_exp(1,:),'filled');
s.MarkerEdgeColor = 'k';
h = colorbar;
colormap Gray;
set(get(h,'label'),'string','Growth rate (/h)');

set(gca,'FontSize',12,'FontName','Helvetica');
xlabel('Protein allocation (g/gCDW)','FontSize',14,'FontName','Helvetica');
ylabel('ATP demand (mmol/gCDW)','FontSize',14,'FontName','Helvetica');

xlim([0.08 0.22]);
% ylim([47 113]);

set(gcf,'position',[800 400 310 250]);
set(gca,'position',[0.13 0.15 0.63 0.8]);

% Change heatmap color: Edit -> Figure properties -> Colormap -> Custom
% Lowest value: [255 255 255]   Highest value: [39 100 25]

clear x y s h;

% Growth rate and protein requirement
figure('Name','4');
hold on;
box on;
x = fluxes_exp(1,2:end);
y = res(8,2:end);
scatter(x,y,50,'o','k');
set(gca,'FontSize',12,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
ylabel(['Protein allocation',char(13,10)','(g/gCDW)'],'FontSize',14,'FontName','Helvetica');

p = polyfit(x,y,1);
p_y = polyval(p,x);
plot(x,p_y,'k-','LineWidth',1.5);

R = corrcoef(x,y);
corr_coef = ['R^2 = ',num2str(round(R(1,2)^2,3))];
text(0.13,0.16,corr_coef,'FontSize',14,'FontName','Helvetica','Color','k');

xlim([0.1 0.25]);
ylim([0.08 0.22]);

set(gcf,'position',[200 0 230 130]);
set(gca,'position',[0.23 0.25 0.72 0.7]);

clear x y corr_coef p p_y R;

% Growth rate and ATP production rate
figure('Name','5');
hold on;
box on;
x = fluxes_exp(1,2:end);
y = res(7,2:end);
scatter(x,y,50,'o','k');
ylim([40 120]);
set(gca,'FontSize',12,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
ylabel(['ATP production rate',char(13,10)','(mmol/gCDW/h)'],'FontSize',14,'FontName','Helvetica');

p = polyfit(x,y,1);
p_y = polyval(p,x);
plot(x,p_y,'k-','LineWidth',1.5);

R = corrcoef(x,y);
corr_coef = ['R^2 = ',num2str(round(R(1,2)^2,3))];
text(0.12,50,corr_coef,'FontSize',14,'FontName','Helvetica','Color','k');

xlim([0.1 0.25]);
ylim([28 62]);

set(gcf,'position',[400 0 230 130]);
set(gca,'position',[0.23 0.25 0.72 0.7]);


clear x y corr_coef p p_y R;
