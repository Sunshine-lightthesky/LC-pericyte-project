close all;clc;clear;
%% Load workspace;
% Load NVC and non-NVC vessel diameter dynamics, ID of vessels, for CBF simulation of LPC intensities;
load('2024_workspace_CBFmodeling_LPCintensity.mat');
% load('../data/2024_workspace_CBFmodeling_LPCintensity.mat');
%% Plot of LPC vessel dynamics, under different LC-NE activity intensities;
colorlist=[[0.0863,0.5216,0.4392];[0.10,0.84,0.73];[1.00,0.30,0.00];[1.00,0.67,0.39];[0.6471,0.1490,0.3922];[0.9216,0.4980,0.6863]];
t=-49:200;

figure,for i=1:6
   p(i)=plot(t,modi_LPC_intensity{1}(:,i),'linewidth',1,'color',colorlist(i,:));hold on
end
ax=gca,xlabel('Time (s)','FontSize',16);ylabel('Vessel diam. change (%)','FontSize',16);ylim([-12,4]);
hold on,plot([0,0],[ax.YLim(1),ax.YLim(2)],'linewidth',0.5,'color',[0 0 0],'linestyle','--','Marker','none');
xlim([-50 200]);ax.YTick=-12:4:4;ylim([-12 4]);
legend([p(1) p(2) p(3) p(4) p(5) p(6)],{'Fore. Cap.','Fore. Art.','Mid. Cap.','Mid. Art.','Hind. Cap.','Hind. Art.'},'Location','SouthEast','FontSize',16);box off;
set(gca,'FontName','Arial','FontSize',20);set(gcf,'color', [1 1 1]);box off;set(gcf,'color', [1 1 1]);set(gca,'Position', [.15 .18 .75 .75]);
% saveas(gcf,'../results/mildLPC_vessel_dynamics.png');

figure,for i=1:6
   p(i)=plot(t,modi_LPC_intensity{2}(:,i),'linewidth',1,'color',colorlist(i,:));hold on
end
ax=gca,xlabel('Time (s)','FontSize',16);ylabel('Vessel diam. change (%)','FontSize',16);ylim([-12,4]);
hold on,plot([0,0],[ax.YLim(1),ax.YLim(2)],'linewidth',0.5,'color',[0 0 0],'linestyle','--','Marker','none');
xlim([-50 200]);ax.YTick=-12:4:4;ylim([-12 4]);
set(gca,'FontName','Arial','FontSize',20);set(gcf,'color', [1 1 1]);box off;set(gcf,'color', [1 1 1]);set(gca,'Position', [.15 .18 .75 .75]);
% saveas(gcf,'../results/mediumLPC_vessel_dynamics.png');

figure,for i=1:6
   p(i)=plot(t,modi_LPC_intensity{3}(:,i),'linewidth',1,'color',colorlist(i,:));hold on
end
ax=gca,xlabel('Time (s)','FontSize',16);ylabel('Vessel diam. change (%)','FontSize',16);ylim([-12,4]);
hold on,plot([0,0],[ax.YLim(1),ax.YLim(2)],'linewidth',0.5,'color',[0 0 0],'linestyle','--','Marker','none');
xlim([-50 200]);ax.YTick=-12:4:4;ylim([-12 4]);
set(gca,'FontName','Arial','FontSize',20);set(gcf,'color', [1 1 1]);box off;set(gcf,'color', [1 1 1]);set(gca,'Position', [.15 .18 .75 .75]);
% saveas(gcf,'../results/strongLPC_vessel_dynamics.png');
%% Plot of NVC vessel dynamics, under different LC-NE activity intensities;
for i=1:3
figure,plot(t,NVC_kernel_art{i}*100,'linewidth',1,'color',[0.8196,0.2078,0.3490]);ylabel('Diameter change (%)');xlabel('Time (s)');
hold on,ax=gca,ylim([-1 4]);ax=gca,ax.YTick=-1:1:4;xlim([-50 200]);
plot([0,0],[ax.YLim(1),ax.YLim(2)],'color',[0.1882,0.1882,0.1882],'linestyle','--','linewidth',0.5,'Marker','none');
set(gca,'FontName','Arial','FontSize',20);box off;set(gcf,'color', [1 1 1]);set(gca,'Position', [.15 .18 .75 .75]); 
end
for i=1:3
figure,plot(t,NVC_kernel_cap{i}*100,'linewidth',1,'color',[0.9490,0.7529,0.3020]);ylabel('Diameter change (%)');xlabel('Time (s)');
hold on,ax=gca,ylim([-1 4]);ax=gca,ax.YTick=-1:1:4;xlim([-50 200]);
plot([0,0],[ax.YLim(1),ax.YLim(2)],'color',[0.1882,0.1882,0.1882],'linestyle','--','linewidth',0.5,'Marker','none');
set(gca,'FontName','Arial','FontSize',20);box off;set(gcf,'color', [1 1 1]);set(gca,'Position', [.15 .18 .75 .75]); 
end
%% Global LPC dynamics, under different LC-NE activity intensities;
snr=50;px_dBW=0;
t=1:1:250;raw_radii=repmat(vessel_radii,[length(t),1]);
sim_radii=cell(1,3);
for k=1:3
    sim_radii{k}=raw_radii;
    for j=1:length(Global_LPcoupling)
    for i=1:length(Global_LPcoupling{j})
    sim_radii{k}(:,Global_LPcoupling{j}(i))=vessel_radii(Global_LPcoupling{j}(i))*(1+0.01*modi_LPC_intensity{k}(:,j));          
    end
    end
    for m=1:length(vessel_radii)
    A=awgn(sim_radii{k}(:,m),snr,px_dBW);% add noise;
    sim_radii{k}(:,m)=smoothdata(A,'gaussian',10);    
    end
end
Global_sim_radii=sim_radii;
%% CBF simulation, LPC intensities;
vesselnum=length(vessel_radii);

para=[];tic
for arbi=1:length(NVC_art_group)
% NVC vessel ID;
NVC_arteriole_id=NVC_art_group{arbi};
NVC_capillary_id=NVC_cap_group{arbi};
%
NVC_arteriole_idnew=[];
for i=1:length(NVC_arteriole_id)
    NVC_arteriole_idnew(i)=find(pathid==NVC_arteriole_id(i));
end
NVC_capillary_idnew=[];
for i=1:length(NVC_capillary_id)
    NVC_capillary_idnew(i)=find(pathid==NVC_capillary_id(i));
end
NVC_idnew=[NVC_arteriole_idnew,NVC_capillary_idnew];
%% NVC blood flow dynamics in the absence of LPC;
snr=50;px_dBW=0;% SNR and power of Gaussian noise;
raw_radii=repmat(vessel_radii,[length(t),1]);
sim_radii=raw_radii;smooth_radii=[];
for i=1:vesselnum
   if ismember(i,NVC_capillary_idnew) 
   sim_radii(:,i)=awgn(vessel_radii(i)*(1+pure_NVC_cap),snr,px_dBW);
   end
   if ismember(i,NVC_arteriole_idnew)
   sim_radii(:,i)=awgn(vessel_radii(i)*(1+pure_NVC_art),snr,px_dBW);
   end
   if ~ismember(i,NVC_capillary_idnew) & ~ismember(i,NVC_arteriole_idnew)
   sim_radii(:,i)=awgn(raw_radii(:,i),snr,px_dBW);    
   end
   smooth_radii(:,i)=smoothdata(sim_radii(:,i),'gaussian',20);
end
control_Local_sim_radii=smooth_radii;
% Blood flow simulation;
[control_Local_flow_eachvessel,control_Local_flow,control_Local_pressure]=bfmodel(control_Local_sim_radii);
control_Local_dfdf_flow=dfvsfo(abs(control_Local_flow_eachvessel));

% response traces of blood flow;
control_Local_dfdf_flow_NVC=control_Local_dfdf_flow(:,NVC_idnew);

control_Local_amp=[];control_Local_FWHM=[];control_Local_AUC=[];
[control_Local_amp,control_Local_FWHM,control_Local_AUC] = fwhm(control_Local_dfdf_flow);

control_Local_amp_NVC = [];control_Local_FWHM_NVC = [];control_Local_AUC_NVC = [];
control_Local_amp_NVC = control_Local_amp(:,NVC_idnew)';
control_Local_FWHM_NVC = control_Local_FWHM(:,NVC_idnew)';
control_Local_AUC_NVC = control_Local_AUC(:,NVC_idnew)';
%% NVC blood flow dynamics with LPC;
% LPC vessel dynamics;
LPCintensitynum=length(Global_sim_radii);
snr=50;px_dBW=0;
raw_radii=repmat(vessel_radii,[length(t),1]);

Local_sim_radii=cell(1,LPCintensitynum);
for k=1:LPCintensitynum
    for i=1:vesselnum
    if ismember(i,NVC_capillary_idnew) 
    Local_sim_radii{k}(:,i)=awgn(vessel_radii(i)*(1+NVC_kernel_cap{k}),snr,px_dBW);
    end
    if ismember(i,NVC_arteriole_idnew) 
    Local_sim_radii{k}(:,i)=awgn(vessel_radii(i)*(1+NVC_kernel_art{k}),snr,px_dBW);
    end    
    if ~ismember(i,NVC_capillary_idnew) & ~ismember(i,NVC_arteriole_idnew)
    Local_sim_radii{k}(:,i)=awgn(raw_radii(:,i),snr,px_dBW);    
    end    
    Local_sim_radii{k}(:,i)=smoothdata(Local_sim_radii{k}(:,i),'gaussian',20);
    end
end

% LPCintensitynum=length(Global_sim_radii);
Localplusglobal=[];
for ii=1:LPCintensitynum
    Localplusglobal{ii}=Global_sim_radii{ii};
    for i=1:length(NVC_arteriole_idnew)
    Localplusglobal{ii}(:,NVC_arteriole_idnew(i))=Local_sim_radii{ii}(:,NVC_arteriole_idnew(i));
    end
    for i=1:length(NVC_capillary_idnew)
    Localplusglobal{ii}(:,NVC_capillary_idnew(i))=Local_sim_radii{ii}(:,NVC_capillary_idnew(i));
    end    
end

datanum=length(Localplusglobal);% 3 intensity levels;
% Blood flow simulation;
Localplusglobal_flow_eachvessel=cell(1,datanum);Localplusglobal_flow=cell(1,datanum);Localplusglobal_pressure=cell(1,datanum);
Localplusglobal_dfdf_flow=cell(1,datanum);tic
for k=1:datanum
    [Localplusglobal_flow_eachvessel{k},Localplusglobal_flow{k},Localplusglobal_pressure{k}]=bfmodel(Localplusglobal{k});
    Localplusglobal_dfdf_flow{k}=dfvsfo(abs(Localplusglobal_flow_eachvessel{k}));
end
toc

% Calculation of blood flow parameters;
localplusglobal_amp=[];localplusglobal_FWHM=[];localplusglobal_AUC=[];
for i=1:datanum
   [localplusglobal_amp(:,i),localplusglobal_FWHM(:,i),localplusglobal_AUC(:,i)] = fwhm(Localplusglobal_dfdf_flow{i});
end

amp_NVC = [];FWHM_NVC = [];AUC_NVC = [];
amp_NVC = localplusglobal_amp(NVC_idnew,:);
FWHM_NVC = localplusglobal_FWHM(NVC_idnew,:);
AUC_NVC = localplusglobal_AUC(NVC_idnew,:);

% response traces of blood flow;
Localplusglobal_dfdf_flow_NVC=cell(1,datanum);
for i=1:datanum
    Localplusglobal_dfdf_flow_NVC{i}=Localplusglobal_dfdf_flow{i}(:,NVC_idnew);
end
% save blood flow parameters;
i=arbi;
para(i).NVC_arteriole_id=NVC_arteriole_id;
para(i).NVC_capillary_id=NVC_capillary_id;
para(i).NVC_idnew=NVC_idnew;
para(i).control_Local_sim_radii=control_Local_sim_radii;
para(i).control_Local_flow_eachvessel=control_Local_flow_eachvessel;
para(i).control_Local_flow=control_Local_flow;
para(i).control_Local_dfdf_flow=control_Local_dfdf_flow;
para(i).control_Local_dfdf_flow_NVC=control_Local_dfdf_flow_NVC;
para(i).control_Local_amp_NVC=control_Local_amp_NVC;
para(i).control_Local_FWHM_NVC=control_Local_FWHM_NVC;
para(i).control_Local_AUC_NVC=control_Local_AUC_NVC;
para(i).Localplusglobal=Localplusglobal;
para(i).Localplusglobal_flow_eachvessel=Localplusglobal_flow_eachvessel;
para(i).Localplusglobal_flow=Localplusglobal_flow;
para(i).Localplusglobal_pressure=Localplusglobal_pressure;
para(i).Localplusglobal_dfdf_flow=Localplusglobal_dfdf_flow;
para(i).Localplusglobal_dfdf_flow_NVC=Localplusglobal_dfdf_flow_NVC;
para(i).localplusglobal_amp=localplusglobal_amp;
para(i).localplusglobal_FWHM=localplusglobal_FWHM;
para(i).localplusglobal_AUC=localplusglobal_AUC;
para(i).amp_NVC=amp_NVC;
para(i).FWHM_NVC=FWHM_NVC;
para(i).AUC_NVC=AUC_NVC;
end
toc
%% data from multiple simulation trials;
% blood flow traces;
simtrialnum=length(para);% number of simulation trials;
All.Localplusglobal_dfdf_flow_NVC=cell(1,datanum);
for i=1:datanum
for j=1:simtrialnum
   All.Localplusglobal_dfdf_flow_NVC{i}=[All.Localplusglobal_dfdf_flow_NVC{i},para(j).Localplusglobal_dfdf_flow_NVC{i}];   
end
end
All.control_Local_dfdf_flow_NVC=[];
for j=1:simtrialnum
   All.control_Local_dfdf_flow_NVC=[All.control_Local_dfdf_flow_NVC,para(j).control_Local_dfdf_flow_NVC];
end

All.amp_NVC=[];All.FWHM_NVC=[];All.AUC_NVC=[];
for i=1:simtrialnum
    All.amp_NVC=cat(1,All.amp_NVC,para(i).amp_NVC);
    All.FWHM_NVC=cat(1,All.FWHM_NVC,para(i).FWHM_NVC);
    All.AUC_NVC=cat(1,All.AUC_NVC,para(i).AUC_NVC);
end
All.control_amp_NVC=[];All.control_FWHM_NVC=[];All.control_AUC_NVC=[];
for i=1:simtrialnum
    All.control_amp_NVC=cat(1,All.control_amp_NVC,para(i).control_Local_amp_NVC);
    All.control_FWHM_NVC=cat(1,All.control_FWHM_NVC,para(i).control_Local_FWHM_NVC);
    All.control_AUC_NVC=cat(1,All.control_AUC_NVC,para(i).control_Local_AUC_NVC);
end
%% plot of blood flow dynamics;
trace=All.Localplusglobal_dfdf_flow_NVC;
trace{4}=All.control_Local_dfdf_flow_NVC*100;% change to percent;
for i=1:3
    trace{i}=trace{i}*100;% change to percent;
end
mean_data=[];std_data=[];sem_data=[];
for i=1:4%datanum+1
   mean_data(:,i)=mean(trace{i}');
   std_data(:,i)=std(trace{i}');
   sem_data(:,i)=std_data(:,i)/sqrt(size(trace{i},2));
end

t = -49:200;
colorlist=[[0.2201,0.3433,0.5494];[0.1195,0.6075,0.5402];[0.5258,0.8335,0.2881];[0.4431,0.4471,0.4784]];
figure,for i=1:4%datanum
   shadedErrorBar(t,mean_data(:,i)',sem_data(:,i),'lineprops',{'color',colorlist(i,:)},'transparent',1,'patchSaturation', 0.05);hold on;    
end
hold on,for i=1:4%datanum
   p(i)=plot(t,mean_data(:,i),'linewidth',1,'color',colorlist(i,:),'linestyle','-','Marker','none');hold on
end
ylim([-5 40]);ax=gca,ax.YTick=-20:10:40;
hold on;plot([0,0],[ax.YLim(1),ax.YLim(2)],'linewidth',0.5,'color',[0.1882,0.1882,0.1882],'linestyle','--','Marker','none');
xlim([-50 200]);xlabel( 'Time (s) ','FontSize',18,'FontName','Arial');ylabel('Local CBF change (%)','FontSize',18,'FontName','Arial');
legend([p(1) p(2) p(3) p(4)],{'NVC + miLPC','NVC + meLPC','NVC + stLPC','NVC'},'FontSize',16);set(gca,'Position', [.15 .18 .75 .75]);
set(gca,'FontName','Arial','FontSize',20);set(gcf,'color', [1 1 1]);box off;
% saveas(gcf,'../results/Local_CBF_change_LPCintensity_trace.png');
%% Summary data, Local NVC dynamics under different LPC intensities;
colorface=[[0.5686,0.5961,0.6824];[0.2012,0.3837,0.5543];[0.1263,0.6441,0.5253];[0.6266,0.8546,0.2234]];
coloredge=[[0.4431,0.4471,0.4784];[0.2718,0.2093,0.5044];[0.1534,0.4970,0.5577];[0.2669,0.7488,0.4406]];

A=[All.control_amp_NVC,All.amp_NVC]*100;% Peak amplitude, change to percent;
barData1=A;p1=[];
[p1(1),h,stats] = signrank(barData1(:,1),barData1(:,2), 'alpha', 0.05, 'tail', 'both');
[p1(2),h,stats] = signrank(barData1(:,2),barData1(:,3), 'alpha', 0.05, 'tail', 'both');
[p1(3),h,stats] = signrank(barData1(:,3),barData1(:,4), 'alpha', 0.05, 'tail', 'both');
FDR=mafdr(p1,'BHFDR', true);% P correction;
Y=barData1;X=ones(size(Y)).*(1:size(Y,2));group_num=size(barData1,2);
figure,for i=1:group_num
    barHdl(i)=bar(i,nanmean(barData1(:,i)),'FaceColor',colorface(i,:),'LineWidth',2,'EdgeColor',coloredge(i,:));hold on
end
hold on,plot(X',Y','Color',[0,0,0,0.1],'Marker','o','LineWidth',1,'MarkerEdgeColor',[1,1,1].*.6,'MarkerSize',5,'LineWidth',.1);hold on
for i=1:group_num % plot error bars;
   errorbar(i,nanmean(barData1(:,i)),nanstd(barData1(:,i),1)./sqrt(length(find(~isnan(barData1(:,i))))),...
    'vertical','LineStyle','none','LineWidth',2,'Color',coloredge(i,:),'CapSize',18);hold on
end
barBaseLineHdl=barHdl(1).BaseLine;
barBaseLineHdl.LineWidth=1.2;hold on;ax=gca,
ax.Color='none';ax.LineWidth=1.5;ax.TickDir='in';ax.FontSize=12;
ax.XTick=1:1:group_num;ax.Title.FontWeight='bold';
ax.Title.FontSize=18;ax.YDir='normal';box off;
ax.YTick=0:20:120;ax.YLim=[0,120];
ylabel('Peak (%)','FontName','Arial','Fontsize',16);box off %血流增加的最大幅值;
plot([1.1,1.9],[70,70],'Color',[0,0,0],'LineWidth',1.5);% 添加显著性标志
plot([2.1,2.9],[70,70],'Color',[0,0,0],'LineWidth',1.5);
plot([3.1,3.9],[70,70],'Color',[0,0,0],'LineWidth',1.5);
text(1.5,75,num2str(p1(1)),'FontSize',12,'HorizontalAlignment','center','FontWeight','normal');
text(2.5,75,num2str(p1(2)),'FontSize',12,'HorizontalAlignment','center','FontWeight','normal');
text(3.5,75,num2str(p1(3)),'FontSize',12,'HorizontalAlignment','center','FontWeight','normal');
set(gca,'FontName','Arial','FontSize',20);set(gcf,'color', [1 1 1]);box off;
dr_name = {'','',''};set(ax(1),'XTickLabel',dr_name,'xtick',1:1:3);box off;
legend('NVC','NVC + miLPC','NVC + meLPC','NVC+stLPC','FontSize',14,'Location','NorthWest');
% saveas(gcf,'../results/Local_CBF_change_LPCintensity_peakamplitude.png');
%
A=[All.control_FWHM_NVC,All.FWHM_NVC];% FWHM;
barData1=A;p1=[];
[p1(1),h,stats] = signrank(barData1(:,1),barData1(:,2), 'alpha', 0.05, 'tail', 'both');
[p1(2),h,stats] = signrank(barData1(:,2),barData1(:,3), 'alpha', 0.05, 'tail', 'both');
[p1(3),h,stats] = signrank(barData1(:,3),barData1(:,4), 'alpha', 0.05, 'tail', 'both');
FDR=mafdr(p1,'BHFDR', true);% P correction;
Y=barData1;X=ones(size(Y)).*(1:size(Y,2));
figure,for i=1:group_num
    barHdl(i)=bar(i,nanmean(barData1(:,i)),'FaceColor',colorface(i,:),'LineWidth',2,'EdgeColor',coloredge(i,:));hold on
end
hold on,plot(X',Y','Color',[0,0,0,0.1],'Marker','o','LineWidth',1,'MarkerEdgeColor',[1,1,1].*.6,'MarkerSize',5,'LineWidth',.1);hold on
for i=1:group_num % plot error bars;
   errorbar(i,nanmean(barData1(:,i)),nanstd(barData1(:,i),1)./sqrt(length(find(~isnan(barData1(:,i))))),...
    'vertical','LineStyle','none','LineWidth',2,'Color',coloredge(i,:),'CapSize',18);hold on
end
barBaseLineHdl=barHdl(1).BaseLine;
barBaseLineHdl.LineWidth=1.2;hold on;ax=gca,
ax.Color='none';ax.LineWidth=1.5;ax.TickDir='in';ax.FontSize=12;
ax.XTick=1:1:group_num;ax.Title.FontWeight='bold';
ax.Title.FontSize=18;ax.YDir='normal';box off;
ax.YTick=0:50:300;ax.YLim=[0,200];
ylabel('FWHM (s)','FontName','Arial','Fontsize',16);box off %在历次event中，有多大概率有响应;
plot([1.1,1.9],[160,160],'Color',[0,0,0],'LineWidth',1.5);% 添加显著性标志
plot([2.1,2.9],[160,160],'Color',[0,0,0],'LineWidth',1.5);
plot([3.1,3.9],[160,160],'Color',[0,0,0],'LineWidth',1.5);
text(1.5,167,num2str(p1(1)),'FontSize',12,'HorizontalAlignment','center','FontWeight','normal');
text(2.5,167,num2str(p1(2)),'FontSize',12,'HorizontalAlignment','center','FontWeight','normal');
text(3.5,167,num2str(p1(3)),'FontSize',12,'HorizontalAlignment','center','FontWeight','normal');
set(gca,'FontName','Arial','FontSize',20);set(gcf,'color', [1 1 1]);box off;
dr_name = {'','',''};set(ax(1),'XTickLabel',dr_name,'xtick',1:1:3);box off;
% saveas(gcf,'../results/Local_CBF_change_LPCintensity_FWHM.png');
%
A=[All.control_AUC_NVC,All.AUC_NVC];% AUC;
barData1=A;p1=[];
[p1(1),h,stats] = signrank(barData1(:,1),barData1(:,2), 'alpha', 0.05, 'tail', 'both');
[p1(2),h,stats] = signrank(barData1(:,2),barData1(:,3), 'alpha', 0.05, 'tail', 'both');
[p1(3),h,stats] = signrank(barData1(:,3),barData1(:,4), 'alpha', 0.05, 'tail', 'both');
FDR=mafdr(p1,'BHFDR', true);% P correction;
Y=barData1;X=ones(size(Y)).*(1:size(Y,2));
figure,for i=1:group_num
    barHdl(i)=bar(i,nanmean(barData1(:,i)),'FaceColor',colorface(i,:),'LineWidth',2,'EdgeColor',coloredge(i,:));hold on
end
hold on,plot(X',Y','Color',[0,0,0,0.1],'Marker','o','LineWidth',1,'MarkerEdgeColor',[1,1,1].*.6,'MarkerSize',5,'LineWidth',.1);hold on
for i=1:group_num % plot error bars;
   errorbar(i,nanmean(barData1(:,i)),nanstd(barData1(:,i),1)./sqrt(length(find(~isnan(barData1(:,i))))),...
    'vertical','LineStyle','none','LineWidth',2,'Color',coloredge(i,:),'CapSize',18);hold on
end
barBaseLineHdl=barHdl(1).BaseLine;
barBaseLineHdl.LineWidth=1.2;hold on;ax=gca,
ax.Color='none';ax.LineWidth=1.5;ax.TickDir='in';ax.FontSize=12;
ax.XTick=1:1:group_num;ax.Title.FontWeight='bold';
ax.Title.FontSize=18;ax.YDir='normal';box off;
ax.YTick=-20:20:100;ax.YLim=[0,100];
ylabel('Total flow (AUC)','FontName','Arial','Fontsize',16);box off %在历次event中，有多大概率有响应;
plot([1.1,1.9],[80,80],'Color',[0,0,0],'LineWidth',1.5);% 添加显著性标志
plot([2.1,2.9],[80,80],'Color',[0,0,0],'LineWidth',1.5);
plot([3.1,3.9],[80,80],'Color',[0,0,0],'LineWidth',1.5);
text(1.5,83,num2str(p1(1)),'FontSize',12,'HorizontalAlignment','center','FontWeight','normal');
text(2.5,83,num2str(p1(2)),'FontSize',12,'HorizontalAlignment','center','FontWeight','normal');
text(3.5,83,num2str(p1(3)),'FontSize',12,'HorizontalAlignment','center','FontWeight','normal');
set(gca,'FontName','Arial','FontSize',20);set(gcf,'color', [1 1 1]);box off;
dr_name = {'','',''};set(ax(1),'XTickLabel',dr_name,'xtick',1:1:3);box off;
% saveas(gcf,'../results/Local_CBF_change_LPCintensity_AUC.png');
%%
% save('../results/CBFsimulation_LPCintensity');
save('results_CBFsimulation_LPCintensity');
