% Analyse models
% Related to Fig 2BC

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

% E.coli
tot_prot_list_ecoli = 0.05:0.05:0.15;
glc_list_ecoli = 0.1:0.1:20;
fluxes_ecoli = zeros(5,length(glc_list_ecoli),length(tot_prot_list_ecoli));
for i = 1:length(tot_prot_list_ecoli)
    tot_prot = tot_prot_list_ecoli(i);
    
    for j = 1:length(glc_list_ecoli)
        glc_in = glc_list_ecoli(j);
        model = changeRxnBounds(model_ecoli, 'EXglc', -glc_in, 'b');
        sol = solveModel(model,'ATPM','max',prot_cost_ecoli,tot_prot);
        if sol.exitflag == 1
            HY = sol.fluxes(strcmp(model.rxns,'GLCtex_HY'))*23.5;
            LY = sol.fluxes(strcmp(model.rxns,'GLCtex_LY'))*11;
            ATP = sol.fluxes(strcmp(model.rxns,'ATPM'));
            [HYp, LYp, Biop] = calculateProtein(model, prot_cost_ecoli, sol.fluxes);
            prot = HYp + LYp + Biop;
            fluxes_ecoli(:,j,i) = [HY; LY; ATP; glc_in; prot];
        end
    end
    
end
fluxes_ecoli(fluxes_ecoli<0.0001) = 0;

clear ATP glc_in HY i j LY model sol tot_prot HYp LYp Biop prot;

% Yeast
tot_prot_list_yeast = 0.05:0.05:0.15;
glc_list_yeast = 0.1:0.1:40;
fluxes_yeast = zeros(5,length(glc_list_yeast),length(tot_prot_list_yeast));
for i = 1:length(tot_prot_list_yeast)
    tot_prot = tot_prot_list_yeast(i);
    
    for j = 1:length(glc_list_yeast)
        glc_in = glc_list_yeast(j);
        model = changeRxnBounds(model_yeast, 'EXglc', -glc_in, 'b');
        sol = solveModel(model,'ATPM','max',prot_cost_yeast,tot_prot);
        if sol.exitflag == 1
            HY = sol.fluxes(strcmp(model.rxns,'GLCt1_HY'))*22;
            LY = sol.fluxes(strcmp(model.rxns,'GLCt1_LY'))*2;
            ATP = sol.fluxes(strcmp(model.rxns,'ATPM'));
            [HYp, LYp, Biop] = calculateProtein(model, prot_cost_yeast, sol.fluxes);
            prot = HYp + LYp + Biop;
            fluxes_yeast(:,j,i) = [HY; LY; ATP; glc_in; prot];
        end
    end
    
end
fluxes_yeast(fluxes_yeast<0.0001) = 0;

clear ATP glc_in HY i j LY model sol tot_prot HYp LYp Biop prot;

%% Figure

figure('Name','ecoli_3d');
% for i = 1:length(tot_prot_list_ecoli)
%     tot_prot = tot_prot_list_ecoli(i);
%     ATP_tmp = fluxes_ecoli(:,:,i);
%     num = length(find(ATP_tmp(3,:) > 0));
%     prot_data = tot_prot * ones(1,num);
%     ATP_data = ATP_tmp(:,1:num);
%     
%     plot3(ATP_data(4,:),prot_data,ATP_data(3,:),'k','LineWidth',1.5);
%     hold on;
%     h = fill3([ATP_data(4,:),fliplr(ATP_data(4,:))],[prot_data,prot_data],[ATP_data(2,:),zeros(1,length(ATP_data(4,:)))],[215,25,28]/256);
% 	
%     set(h,'edgealpha',0,'facealpha',0.7);
%     h = fill3([ATP_data(4,:),fliplr(ATP_data(4,:))],[prot_data,prot_data],[ATP_data(3,:),fliplr(ATP_data(2,:))],[43,131,186]/256);
% 	set(h,'edgealpha',0,'facealpha',0.7);
% 
% end
% 
% grid on;
% set(gca,'FontSize',10,'FontName','Helvetica');
% zlim([0 220]);
% ylabel(['Protein allocation',char(13,10)','(g/gCDW)'],'FontSize',12,'FontName','Helvetica');
% xlabel(['Glucose uptake rate',char(13,10)','(mmol/gCDW/h)'],'FontSize',12,'FontName','Helvetica');
% zlabel(['ATP production rate',char(13,10)','(mmol/gCDW/h)'],'FontSize',12,'FontName','Helvetica');
% legend({'Total ATP';'ATP from Low-yield pathway';'ATP from High-yield pathway'},'FontSize',8,'FontName','Helvetica','location','ne');
% 
% az = -17;
% el = 53;
% view(az, el);
% 
% set(gcf,'position',[0 0 350 300]);
% set(gca,'position',[0.25 0.25 0.65 0.65]);

x1 = fluxes_ecoli(4,1,1);
y1 = 0.05;
z1 = fluxes_ecoli(3,1,1);
idx_tmp2 = find(fluxes_ecoli(2,:,1),1)-1;
x2 = fluxes_ecoli(4,idx_tmp2,1);
y2 = 0.05;
z2 = fluxes_ecoli(3,idx_tmp2,1);
x3 = fluxes_ecoli(4,1,3);
y3 = 0.15;
z3 = fluxes_ecoli(3,1,3);
idx_tmp4 = find(fluxes_ecoli(2,:,3),1);
x4 = fluxes_ecoli(4,idx_tmp4,3);
y4 = 0.15;
z4 = fluxes_ecoli(3,idx_tmp4,3);
C(:,:,1) = [0 0;0 0];
C(:,:,2) = [0 0;0 0];
C(:,:,3) = [0 0;0 0];
surf([x1 x2; x3 x4],[y1 y2;y3 y4],[z1 z2;z3 z4],C,'FaceAlpha',0.5,'EdgeColor','none');
hold on;
x5 = max(fluxes_ecoli(4,:,1));
y5 = 0.05;
z5 = max(fluxes_ecoli(3,:,1));
x6 = max(fluxes_ecoli(4,:,3));
y6 = 0.15;
z6 = max(fluxes_ecoli(3,:,3));
surf([x2 x5;x4 x6],[y2 y5;y4 y6],[z2 z5;z4 z6],C,'FaceAlpha',0.2,'EdgeColor','none');
clear x1 y1 z1 x2 y2 z2 x3 y3 z3 x4 y4 z4 x5 y5 z5 x6 y6 z6 C idx_tmp2 idx_tmp4;

for i = 1:length(tot_prot_list_ecoli)
    tot_prot = tot_prot_list_ecoli(i);
    ATP_tmp = fluxes_ecoli(:,:,i);
    num = length(find(ATP_tmp(3,:) > 0));
    prot_data = tot_prot * ones(1,num);
    ATP_data = ATP_tmp(:,1:num);
    
    if i == 2
        plot3(ATP_data(4,:),prot_data,ATP_data(3,:),'k','LineWidth',1);
        hold on;
        plot3([0,0],[tot_prot,tot_prot],[0,ATP_data(1,1)],'k-.','LineWidth',1);
        plot3([ATP_data(4,end),ATP_data(4,end)],[tot_prot,tot_prot],[0,ATP_data(3,end)],'k-.','LineWidth',1);
        plot3([0,ATP_data(4,end)],[tot_prot,tot_prot],[0,0],'k-.','LineWidth',1);
    else
        plot3(ATP_data(4,:),prot_data,ATP_data(3,:),'k','LineWidth',0.5);
        hold on;
%         plot3([0,0],[tot_prot,tot_prot],[0,ATP_data(1,1)],'k','LineWidth',0.5);
%         plot3([ATP_data(4,end),ATP_data(4,end)],[tot_prot,tot_prot],[0,ATP_data(3,end)],'k','LineWidth',0.5);
%         plot3([0,ATP_data(4,end)],[tot_prot,tot_prot],[0,0],'k','LineWidth',0.5);
    end

end

grid on;
set(gca,'FontSize',12,'FontName','Helvetica');
zlim([0 220]);
ylabel(['Protein allocation',char(13,10)','(g/gCDW)'],'FontSize',14,'FontName','Helvetica');
xlabel(['Glucose uptake rate',char(13,10)','(mmol/gCDW/h)'],'FontSize',14,'FontName','Helvetica');
zlabel(['ATP production rate',char(13,10)','(mmol/gCDW/h)'],'FontSize',14,'FontName','Helvetica');

az = -45;
el = 30;
view(az, el);

set(gcf,'position',[0 500 350 300]);
set(gca,'position',[0.25 0.25 0.65 0.65]);

clear i tot_prot ATP_tmp num prot_data ATP_data h az el;

figure('Name','yeast_3d');
% for i = 1:length(tot_prot_list_yeast)
%     tot_prot = tot_prot_list_yeast(i);
%     ATP_tmp = fluxes_yeast(:,:,i);
%     num = length(find(ATP_tmp(3,:) > 0));
%     prot_data = tot_prot * ones(1,num);
%     ATP_data = ATP_tmp(:,1:num);
%     
%     plot3(ATP_data(4,:),prot_data,ATP_data(3,:),'k','LineWidth',1.5);
%     hold on;
%     h = fill3([ATP_data(4,:),fliplr(ATP_data(4,:))],[prot_data,prot_data],[ATP_data(2,:),zeros(1,length(ATP_data(4,:)))],[215,25,28]/256);
% 	
%     set(h,'edgealpha',0,'facealpha',0.7);
%     h = fill3([ATP_data(4,:),fliplr(ATP_data(4,:))],[prot_data,prot_data],[ATP_data(3,:),fliplr(ATP_data(2,:))],[43,131,186]/256);
% 	set(h,'edgealpha',0,'facealpha',0.7);
% 
% end
% grid on;
% set(gca,'FontSize',10,'FontName','Helvetica');
% xlim([0 37]);
% zlim([0 80]);
% ylabel(['Protein allocation',char(13,10)','(g/gCDW)'],'FontSize',12,'FontName','Helvetica');
% xlabel(['Glucose uptake rate',char(13,10)','(mmol/gCDW/h)'],'FontSize',12,'FontName','Helvetica');
% zlabel(['ATP production rate',char(13,10)','(mmol/gCDW/h)'],'FontSize',12,'FontName','Helvetica');
% legend({'Total ATP';'ATP from Low-yield pathway';'ATP from High-yield pathway'},'FontSize',8,'FontName','Helvetica','location','ne');
% 
% az = -17;
% el = 53;
% view(az, el);
% 
% set(gcf,'position',[500 0 350 300]);
% set(gca,'position',[0.25 0.25 0.65 0.65]);

x1 = fluxes_yeast(4,1,1);
y1 = 0.05;
z1 = fluxes_yeast(3,1,1);
idx_tmp2 = find(fluxes_yeast(2,:,1),1)-1;
x2 = fluxes_yeast(4,idx_tmp2,1);
y2 = 0.05;
z2 = fluxes_yeast(3,idx_tmp2,1);
x3 = fluxes_yeast(4,1,3);
y3 = 0.15;
z3 = fluxes_yeast(3,1,3);
idx_tmp4 = find(fluxes_yeast(2,:,3),1);
x4 = fluxes_yeast(4,idx_tmp4,3);
y4 = 0.15;
z4 = fluxes_yeast(3,idx_tmp4,3);
C(:,:,1) = [0 0;0 0];
C(:,:,2) = [0 0;0 0];
C(:,:,3) = [0 0;0 0];
surf([x1 x2; x3 x4],[y1 y2;y3 y4],[z1 z2;z3 z4],C,'FaceAlpha',0.5,'EdgeColor','none');
hold on;
x5 = max(fluxes_yeast(4,:,1));
y5 = 0.05;
z5 = max(fluxes_yeast(3,:,1));
x6 = max(fluxes_yeast(4,:,3));
y6 = 0.15;
z6 = max(fluxes_yeast(3,:,3));
surf([x2 x5;x4 x6],[y2 y5;y4 y6],[z2 z5;z4 z6],C,'FaceAlpha',0.2,'EdgeColor','none');
clear x1 y1 z1 x2 y2 z2 x3 y3 z3 x4 y4 z4 x5 y5 z5 x6 y6 z6 C idx_tmp2 idx_tmp4;

for i = 1:length(tot_prot_list_yeast)
    tot_prot = tot_prot_list_yeast(i);
    ATP_tmp = fluxes_yeast(:,:,i);
    num = length(find(ATP_tmp(3,:) > 0));
    prot_data = tot_prot * ones(1,num);
    ATP_data = ATP_tmp(:,1:num);
    
    if i == 2
        plot3(ATP_data(4,:),prot_data,ATP_data(3,:),'k','LineWidth',1);
        hold on;
        plot3([0,0],[tot_prot,tot_prot],[0,ATP_data(1,1)],'k-.','LineWidth',1);
        plot3([ATP_data(4,end),ATP_data(4,end)],[tot_prot,tot_prot],[0,ATP_data(3,end)],'k-.','LineWidth',1);
        plot3([0,ATP_data(4,end)],[tot_prot,tot_prot],[0,0],'k-.','LineWidth',1);
    else
        plot3(ATP_data(4,:),prot_data,ATP_data(3,:),'k','LineWidth',0.5);
        hold on;
%         plot3([0,0],[tot_prot,tot_prot],[0,ATP_data(1,1)],'k','LineWidth',0.5);
%         plot3([ATP_data(4,end),ATP_data(4,end)],[tot_prot,tot_prot],[0,ATP_data(3,end)],'k','LineWidth',0.5);
%         plot3([0,ATP_data(4,end)],[tot_prot,tot_prot],[0,0],'k','LineWidth',0.5);
    end

end

grid on;
set(gca,'FontSize',12,'FontName','Helvetica');
xlim([0 37]);
zlim([0 80]);
ylabel(['Protein allocation',char(13,10)','(g/gCDW)'],'FontSize',14,'FontName','Helvetica');
xlabel(['Glucose uptake rate',char(13,10)','(mmol/gCDW/h)'],'FontSize',14,'FontName','Helvetica');
zlabel(['ATP production rate',char(13,10)','(mmol/gCDW/h)'],'FontSize',14,'FontName','Helvetica');

az = -45;
el = 30;
view(az, el);

set(gcf,'position',[0 0 350 300]);
set(gca,'position',[0.25 0.25 0.65 0.65]);


clear i tot_prot ATP_tmp num prot_data ATP_data h az el;




figure('Name','ecoli_2d');

data = fluxes_ecoli(:,:,2);
num = length(find(data(3,:) > 0));
data = data(:,1:num);
x = data(4,:);
y1_tot = data(3,:);
y1_ly = data(2,:);
y2 = data(5,:);

hold on;
plot(x,y1_tot,'k','LineWidth',1);
h = fill([x,fliplr(x)],[y1_ly,zeros(1,length(x))],[215,25,28]/256);
set(h,'edgealpha',0,'facealpha',0.7);
h = fill([x,fliplr(x)],[y1_tot,fliplr(y1_ly)],[43,131,186]/256);
set(h,'edgealpha',0,'facealpha',0.7);
xlim([0 13]);
ylim([0 150]);
xlabel(['Glucose uptake rate',char(13,10)','(mmol/gCDW/h)'],'FontSize',14,'FontName','Helvetica');
ylabel(['ATP production rate',char(13,10)','(mmol/gCDW/h)'],'FontSize',14,'FontName','Helvetica');
legend({'Total ATP';'ATP from Low-yield pathway';'ATP from High-yield pathway'},'FontSize',8,'FontName','Helvetica','location','ne');

set(gcf,'position',[500 200 180 130]);
set(gca,'position',[0.27 0.33 0.7 0.63]);

figure('Name','ecoli_2d_prot');
hold on;
plot(x,y2,'LineWidth',1,'Color',[197,27,138]/255);
plot(x,0.1*ones(1,length(x)),':','LineWidth',1,'Color',[197,27,138]/255);
ylim([0 0.1]);
set(gca,'ytick',0:0.1:0.1);
xlim([0 13]);
set(gcf,'position',[500 0 180 50]);
set(gca,'position',[0.27 0.35 0.7 0.5]);


clear data h num x y1_ly y1_tot y2;


figure('Name','yeast_2d');

data = fluxes_yeast(:,:,2);
num = length(find(data(3,:) > 0));
data = data(:,1:num);
x = data(4,:);
y1_tot = data(3,:);
y1_ly = data(2,:);
y2 = data(5,:);

hold on;
plot(x,y1_tot,'k','LineWidth',1);
h = fill([x,fliplr(x)],[y1_ly,zeros(1,length(x))],[215,25,28]/256);
set(h,'edgealpha',0,'facealpha',0.7);
h = fill([x,fliplr(x)],[y1_tot,fliplr(y1_ly)],[43,131,186]/256);
set(h,'edgealpha',0,'facealpha',0.7);
xlim([0 25]);
ylim([0 50]);
xlabel(['Glucose uptake rate',char(13,10)','(mmol/gCDW/h)'],'FontSize',14,'FontName','Helvetica');
ylabel(['ATP production rate',char(13,10)','(mmol/gCDW/h)'],'FontSize',14,'FontName','Helvetica');
legend({'Total ATP';'ATP from Low-yield pathway';'ATP from High-yield pathway'},'FontSize',8,'FontName','Helvetica','location','ne');

set(gcf,'position',[700 200 180 130]);
set(gca,'position',[0.27 0.33 0.7 0.63]);

figure('Name','yeast_2d_prot');
hold on;
plot(x,y2,'LineWidth',1,'Color',[197,27,138]/255);
plot(x,0.1*ones(1,length(x)),':','LineWidth',1,'Color',[197,27,138]/255);
ylim([0 0.1]);
set(gca,'ytick',0:0.1:0.1);
xlim([0 25]);
set(gcf,'position',[700 0 180 50]);
set(gca,'position',[0.27 0.35 0.7 0.5]);

clear data h num x y1_ly y1_tot y2;