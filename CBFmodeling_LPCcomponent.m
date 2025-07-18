clear all;close all;clc
%% Load workspace;
% Load NVC and non-NVC vessel diameter dynamics, ID of vessels, for CBF simulation of LPC components;
load('2025_workspace_CBFmodeling.mat');
%% CBF simulation for the contribution of direct of indirect LPC to blood flow regulation;
% NVC_art_group;NVC_cap_group;% ID of NVC vessels in each simulation;
vesselnum=length(vessel_radii);

tic
para=[];
for arbi=1:length(NVC_art_group)
NVC_arteriole_id=NVC_art_group{arbi};
NVC_capillary_id=NVC_cap_group{arbi};

NVC_arteriole_idnew=[];
for i=1:length(NVC_arteriole_id)
    NVC_arteriole_idnew(i)=find(pathid==NVC_arteriole_id(i));
end
NVC_capillary_idnew=[];
for i=1:length(NVC_capillary_id)
    NVC_capillary_idnew(i)=find(pathid==NVC_capillary_id(i));
end
NVC_idnew=[NVC_arteriole_idnew,NVC_capillary_idnew];
%% NVC blood flow dynamics without LPC;
% NVC dynamics with LPC;
snr=50;px_dBW=0;% SNR and power of Gaussian noise;
t=-49:200;
raw_radii=repmat(vessel_radii,[length(t),1]);
sim_radii=raw_radii;smooth_radii=[];
for i=1:vesselnum
   if ismember(i,NVC_capillary_idnew) 
   sim_radii(:,i)=awgn(vessel_radii(i)*(1+control_NVC_cap),snr,px_dBW);
   end
   if ismember(i,NVC_arteriole_idnew)
   sim_radii(:,i)=awgn(vessel_radii(i)*(1+control_NVC_art),snr,px_dBW);
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

% NVC dynamics without LPC (NER blockade);
sim_radii=raw_radii;smooth_radii=[];
for i=1:vesselnum
   if ismember(i,NVC_capillary_idnew) 
   sim_radii(:,i)=awgn(vessel_radii(i)*(1+NERblockade_NVC_cap),snr,px_dBW);
   end
   if ismember(i,NVC_arteriole_idnew)
   sim_radii(:,i)=awgn(vessel_radii(i)*(1+NERblockade_NVC_art),snr,px_dBW);
   end
   if ~ismember(i,NVC_capillary_idnew) & ~ismember(i,NVC_arteriole_idnew)
   sim_radii(:,i)=awgn(raw_radii(:,i),snr,px_dBW);    
   end
   smooth_radii(:,i)=smoothdata(sim_radii(:,i),'gaussian',20);
end
NERblockade_Local_sim_radii=smooth_radii;
% Blood flow simulation;
[NERblockade_Local_flow_eachvessel,NERblockade_Local_flow,NERblockade_Local_pressure]=bfmodel(NERblockade_Local_sim_radii);
NERblockade_Local_dfdf_flow=dfvsfo(abs(NERblockade_Local_flow_eachvessel));
%% NVC blood flow dynamics with LPC;
Localplusglobal=Global_sim_radii;
control_groupid=[1,3];LCablation_groupid=[2];mix_groupid=[4];
for k=1:length(control_groupid)
    for i=1:length(NVC_arteriole_idnew)
    Localplusglobal{control_groupid(k)}(:,NVC_arteriole_idnew(i))=control_Local_sim_radii(:,NVC_arteriole_idnew(i));        
    end
    for i=1:length(NVC_capillary_idnew)
    Localplusglobal{control_groupid(k)}(:,NVC_capillary_idnew(i))=control_Local_sim_radii(:,NVC_capillary_idnew(i));        
    end    
end

for k=1:length(LCablation_groupid)
    for i=1:length(NVC_arteriole_idnew)
    Localplusglobal{LCablation_groupid(k)}(:,NVC_arteriole_idnew(i))=NERblockade_Local_sim_radii(:,NVC_arteriole_idnew(i));  
    end
    for i=1:length(NVC_capillary_idnew)
    Localplusglobal{LCablation_groupid(k)}(:,NVC_capillary_idnew(i))=NERblockade_Local_sim_radii(:,NVC_capillary_idnew(i));  
    end       
end
for k=1:length(mix_groupid)% glia ablation, vasoconstrictive drive on capillaries but not arterioles;
    for i=1:length(NVC_arteriole_idnew)
    Localplusglobal{mix_groupid(k)}(:,NVC_arteriole_idnew(i))=NERblockade_Local_sim_radii(:,NVC_arteriole_idnew(i));  
    end    
    for i=1:length(NVC_capillary_idnew)
    Localplusglobal{mix_groupid(k)}(:,NVC_capillary_idnew(i))=control_Local_sim_radii(:,NVC_capillary_idnew(i));        
    end        
end

datanum=length(Localplusglobal);% 6 experimental contexts;
% Blood flow simulation;
Localplusglobal_flow_eachvessel=cell(1,datanum);Localplusglobal_flow=cell(1,datanum);Localplusglobal_pressure=cell(1,datanum);
Localplusglobal_dfdf_flow=cell(1,datanum);
for k=1:datanum
    [Localplusglobal_flow_eachvessel{k},Localplusglobal_flow{k},Localplusglobal_pressure{k}]=bfmodel(Localplusglobal{k});
    Localplusglobal_dfdf_flow{k}=dfvsfo(abs(Localplusglobal_flow_eachvessel{k}));
end
%% Calculation of blood flow parameters;
localplusglobal_amp=[];localplusglobal_FWHM=[];localplusglobal_AUC=[];
for i=1:datanum
   [localplusglobal_amp(:,i),localplusglobal_FWHM(:,i),localplusglobal_AUC(:,i)] = fwhm(Localplusglobal_dfdf_flow{i});
end

amp_NVC = [];FMHM_NVC = [];AUC_NVC = [];
amp_NVC = localplusglobal_amp(NVC_idnew,:);
FMHM_NVC = localplusglobal_FWHM(NVC_idnew,:);
AUC_NVC = localplusglobal_AUC(NVC_idnew,:);

% response traces of blood flow;
Localplusglobal_dfdf_flow_NVC=cell(1,datanum);
for i=1:datanum
    Localplusglobal_dfdf_flow_NVC{i}=Localplusglobal_dfdf_flow{i}(:,NVC_idnew);
end
%% save blood flow parameters;
i=arbi;
para(i).NVC_arteriole_id=NVC_arteriole_id;
para(i).NVC_capillary_id=NVC_capillary_id;
para(i).NVC_idnew=NVC_idnew;
para(i).control_Local_sim_radii=control_Local_sim_radii;
para(i).control_Local_flow_eachvessel=control_Local_flow_eachvessel;
para(i).control_Local_flow=control_Local_flow;
para(i).control_Local_dfdf_flow=control_Local_dfdf_flow;
para(i).NERblockade_Local_sim_radii=NERblockade_Local_sim_radii;
para(i).NERblockade_Local_flow_eachvessel=NERblockade_Local_flow_eachvessel;
para(i).NERblockade_Local_flow=NERblockade_Local_flow;
para(i).NERblockade_Local_dfdf_flow=NERblockade_Local_dfdf_flow;
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
para(i).FMHM_NVC=FMHM_NVC;
para(i).AUC_NVC=AUC_NVC;
end
toc
%% data from multiple simulation trials;
% Traces of blood flow dynamics;
simtrialnum=length(para);% number of simulations;
All.Localplusglobal_dfdf_flow_NVC=cell(1,datanum);
for i=1:datanum
for j=1:simtrialnum
   All.Localplusglobal_dfdf_flow_NVC{i}=[All.Localplusglobal_dfdf_flow_NVC{i},para(j).Localplusglobal_dfdf_flow_NVC{i}];   
end
end
% figure,plot(mean(All.Localplusglobal_dfdf_flow_NVC{2}'));hold on,
% plot(mean(All.Localplusglobal_dfdf_flow_NVC{4}'));hold on,
% plot(mean(All.Localplusglobal_dfdf_flow_NVC{3}'));

All.amp_NVC=[];All.FMHM_NVC=[];All.AUC_NVC=[];
for i=1:simtrialnum
    All.amp_NVC=cat(1,All.amp_NVC,para(i).amp_NVC);
    All.FMHM_NVC=cat(1,All.FMHM_NVC,para(i).FMHM_NVC);
    All.AUC_NVC=cat(1,All.AUC_NVC,para(i).AUC_NVC);
end
%% traces of CBF responses;
t=-49:200;
% LC ablation;
A=All.Localplusglobal_dfdf_flow_NVC{2}*100;B=All.Localplusglobal_dfdf_flow_NVC{1}*100;
coloredge=[[0.4431,0.4471,0.4784];[0.6510,0.3333,0.6157]];
mean_A=mean(A');std_A=std(A');sem_A=std_A/sqrt(size(A,2));mean_B=mean(B');std_B=std(B');sem_B=std_B/sqrt(size(B,2));
figure,shadedErrorBar(t,mean_A',sem_A,'lineprops',{'color',coloredge(1,:)},'transparent',1,'patchSaturation', 0.1);hold on;    
hold on;shadedErrorBar(t,mean_B',sem_B,'lineprops',{'color',coloredge(2,:)},'transparent',1,'patchSaturation', 0.1);hold on;    
hold on;p(1)=plot(t,mean_A,'linewidth',1,'color',coloredge(1,:),'linestyle','-','Marker','none');
hold on;p(2)=plot(t,mean_B,'linewidth',1,'color',coloredge(2,:),'linestyle','-','Marker','none');
xlim([-50 200]);ax=gca,ylim([-5 30]);ax.YTick=-10:10:30;
hold on;plot([0,0],[ax.YLim(1),ax.YLim(2)],'linewidth',1,'color',[0.1882,0.1882,0.1882],'linestyle','--','Marker','none');
xlabel( 'Time (s)','FontSize',18,'FontName','Arial');ylabel('Blood flow (%)','FontSize',18,'FontName','Arial');
legend([p(1) p(2)],{'NVC','NVC + d-LPC + i-LPC'},'FontSize',18);title('LC ablation','FontSize',18);
set(gca,'FontName','Arial','FontSize',20);set(gcf,'color', [1 1 1]);box off;set(gcf,'color', [1 1 1]);%ylim([-0.3 0.4]);
% saveas(gcf,'../results/Local_CBF_change_LCablation.png');

% RG ablation;
A=All.Localplusglobal_dfdf_flow_NVC{4}*100;B=All.Localplusglobal_dfdf_flow_NVC{3}*100;
coloredge=[[0.3647,0.6392,0.6157];[0.6510,0.3333,0.6157]];
mean_A=mean(A');std_A=std(A');sem_A=std_A/sqrt(size(A,2));mean_B=mean(B');std_B=std(B');sem_B=std_B/sqrt(size(B,2));
figure,shadedErrorBar(t,mean_A',sem_A,'lineprops',{'color',coloredge(1,:)},'transparent',1,'patchSaturation', 0.1);hold on;    
hold on;shadedErrorBar(t,mean_B',sem_B,'lineprops',{'color',coloredge(2,:)},'transparent',1,'patchSaturation', 0.1);hold on;    
hold on;p(1)=plot(t,mean_A,'linewidth',1,'color',coloredge(1,:),'linestyle','-','Marker','none');
hold on;p(2)=plot(t,mean_B,'linewidth',1,'color',coloredge(2,:),'linestyle','-','Marker','none');
xlim([-50 200]);ax=gca,ylim([-5 30]);ax.YTick=-10:10:30;
hold on;plot([0,0],[ax.YLim(1),ax.YLim(2)],'linewidth',1,'color',[0.1882,0.1882,0.1882],'linestyle','--','Marker','none');
xlabel( 'Time (s)','FontSize',18,'FontName','Arial');ylabel('Blood flow (%)','FontSize',18,'FontName','Arial');
legend([p(1) p(2)],{'NVC + d-LPC','NVC + d-LPC + i-LPC'},'FontSize',18);title('RG ablation','FontSize',18);
set(gca,'FontName','Arial','FontSize',20);set(gcf,'color', [1 1 1]);box off;set(gcf,'color', [1 1 1]);%ylim([-0.3 0.4]);
% saveas(gcf,'../results/Local_CBF_change_RGablation.png');

% NVC, NVC + d-LPC, NVC + i-LPC put together, Fig. 7c;
A=All.Localplusglobal_dfdf_flow_NVC{2}*100;B=All.Localplusglobal_dfdf_flow_NVC{4}*100;C=All.Localplusglobal_dfdf_flow_NVC{3}*100;
coloredge=[[0.4431,0.4471,0.4784];[0.3647,0.6392,0.6157];[0.6510,0.3333,0.6157]];
mean_A=mean(A');std_A=std(A');sem_A=std_A/sqrt(size(A,2));
mean_B=mean(B');std_B=std(B');sem_B=std_B/sqrt(size(B,2));
mean_C=mean(C');std_C=std(C');sem_C=std_C/sqrt(size(C,2));

figure,shadedErrorBar(t,mean_A',sem_A,'lineprops',{'color',coloredge(1,:)},'transparent',1,'patchSaturation', 0.1);hold on;    
hold on;shadedErrorBar(t,mean_B',sem_B,'lineprops',{'color',coloredge(2,:)},'transparent',1,'patchSaturation', 0.1);hold on;    
hold on;shadedErrorBar(t,mean_C',sem_C,'lineprops',{'color',coloredge(3,:)},'transparent',1,'patchSaturation', 0.1);hold on;    
hold on;p(1)=plot(t,mean_A,'linewidth',1,'color',coloredge(1,:),'linestyle','-','Marker','none');
hold on;p(2)=plot(t,mean_B,'linewidth',1,'color',coloredge(2,:),'linestyle','-','Marker','none');
hold on;p(3)=plot(t,mean_C,'linewidth',1,'color',coloredge(3,:),'linestyle','-','Marker','none');
xlim([-50 200]);ax=gca,ylim([-5 30]);ax.YTick=-10:10:30;
hold on;plot([0,0],[ax.YLim(1),ax.YLim(2)],'linewidth',1,'color',[0.1882,0.1882,0.1882],'linestyle','--','Marker','none');
xlim([-50 200]);xlabel( 'Time (s)','FontSize',18,'FontName','Arial');ylabel('Local CBF change (%)','FontSize',18,'FontName','Arial');
legend([p(1) p(2) p(3)],{'NVC','NVC + d-LPC','NVC + d-LPC + i-LPC'},'FontSize',18);
set(gca,'FontName','Arial','FontSize',20);set(gcf,'color', [1 1 1]);box off;set(gcf,'color', [1 1 1]);%ylim([-0.3 0.4]);
% saveas(gcf,'../results/Local_CBF_change_NVC_d-LPC_i-LPC.png');
%% Summary data of CBF parameters, Extended Data Fig. 11;
%% LC ablation;
coloredge=[[0.4431,0.4471,0.4784];[0.6510,0.3333,0.6157]];colorface=[[0.5686,0.5961,0.6824];[0.8078,0.5765,0.7490]];colorerrorbar=coloredge;

A=All.amp_NVC*100;barData1=A(:,[2,1]);% peak;
[p1,h,stats] = signrank(barData1(:,1),barData1(:,2), 'alpha', 0.05, 'tail', 'both');group_num=size(barData1,2);
Y=barData1;X=ones(size(Y)).*(1:size(Y,2));figure,for i=1:group_num
    barHdl(i)=bar(i,nanmean(barData1(:,i)),0.8,'FaceColor',colorface(i,:),'FaceAlpha',0.5,'LineWidth',2,'EdgeColor',coloredge(i,:));hold on
end
barBaseLineHdl=barHdl(1).BaseLine;barBaseLineHdl.LineWidth=1.2;
hold on,plot(X',Y','Color',[0,0,0,0.1],'Marker','o','LineWidth',1,'MarkerEdgeColor',[1,1,1].*.6,'MarkerSize',5,'LineWidth',.1);hold on
for i=1:group_num % plot error bars;
   errorbar(i,nanmean(barData1(:,i)),nanstd(barData1(:,i),1)./sqrt(length(find(~isnan(barData1(:,i))))),...
    'vertical','LineStyle','none','LineWidth',2,'Color',colorerrorbar(i,:),'CapSize',18);hold on
end
hold on;ax=gca,ax.Color='none';ax.LineWidth=1.5;ax.TickDir='in';ax.FontSize=14;ax.XTick=1:1:2;box off
ylabel('Peak (%)','FontName','Arial','Fontsize',16);box off;hold on,
ax.Title.FontWeight='bold';ax.Title.FontSize=18;ax.YDir='normal';
ax.YTick=0:20:80;ylim([0,80]);
plot([1.1,1.9],[60,60],'Color',[0,0,0],'LineWidth',1.5);hold on,
set(gca,'FontName','Arial','FontSize',18);set(gcf,'color', [1 1 1]);
dr_name = {'w/o LPC','w LPC'};set(ax(1),'XTickLabel',dr_name,'xtick',1:1:2);box off;%,'XTickLabelRotation',30
text(1.5,64,num2str(p1),'FontSize',20,'HorizontalAlignment','center','FontWeight','normal');
% text(1.5,62,'****','FontSize',24,'HorizontalAlignment','center','FontWeight','normal');
% saveas(gcf,'../results/Peak_LCablation.png');

A=All.FMHM_NVC;barData1=A(:,[2,1]);% Full width at half maximum;
[p1,h,stats] = signrank(barData1(:,1),barData1(:,2), 'alpha', 0.05, 'tail', 'both');group_num=size(barData1,2);
Y=barData1;X=ones(size(Y)).*(1:size(Y,2));figure,for i=1:group_num
    barHdl(i)=bar(i,nanmean(barData1(:,i)),0.8,'FaceColor',colorface(i,:),'FaceAlpha',0.5,'LineWidth',2,'EdgeColor',coloredge(i,:));hold on
end
barBaseLineHdl=barHdl(1).BaseLine;barBaseLineHdl.LineWidth=1.2;
hold on,plot(X',Y','Color',[0,0,0,0.1],'Marker','o','LineWidth',1,'MarkerEdgeColor',[1,1,1].*.6,'MarkerSize',5,'LineWidth',.1);hold on
for i=1:group_num % plot error bars;
   errorbar(i,nanmean(barData1(:,i)),nanstd(barData1(:,i),1)./sqrt(length(find(~isnan(barData1(:,i))))),...
    'vertical','LineStyle','none','LineWidth',2,'Color',colorerrorbar(i,:),'CapSize',18);hold on
end
hold on;ax=gca,ax.Color='none';ax.LineWidth=1.5;ax.TickDir='in';ax.FontSize=14;ax.XTick=1:1:2;box off
ylabel('FWHM (s)','FontName','Arial','Fontsize',16);box off;hold on,
ax.Title.FontWeight='bold';ax.Title.FontSize=18;ax.YDir='normal';
ax.YTick=0:50:200;ylim([0,200]);
plot([1.1,1.9],[150,150],'Color',[0,0,0],'LineWidth',1.5);hold on,
set(gca,'FontName','Arial','FontSize',18);set(gcf,'color', [1 1 1]);
dr_name = {'w/o LPC','w LPC'};set(ax(1),'XTickLabel',dr_name,'xtick',1:1:2);box off;
text(1.5,160,num2str(p1),'FontSize',20,'HorizontalAlignment','center','FontWeight','normal');
% text(1.5,156,'****','FontSize',24,'HorizontalAlignment','center','FontWeight','normal');
% saveas(gcf,'../results/FWHM_LCablation.png');

A=All.AUC_NVC;barData1=A(:,[2,1]);% AUC;
[p1,h,stats] = signrank(barData1(:,1),barData1(:,2), 'alpha', 0.05, 'tail', 'both');group_num=size(barData1,2);
Y=barData1;X=ones(size(Y)).*(1:size(Y,2));figure,for i=1:group_num
    barHdl(i)=bar(i,nanmean(barData1(:,i)),0.8,'FaceColor',colorface(i,:),'FaceAlpha',0.5,'LineWidth',2,'EdgeColor',coloredge(i,:));hold on
end
barBaseLineHdl=barHdl(1).BaseLine;barBaseLineHdl.LineWidth=1.2;
hold on,plot(X',Y','Color',[0,0,0,0.1],'Marker','o','LineWidth',1,'MarkerEdgeColor',[1,1,1].*.6,'MarkerSize',5,'LineWidth',.1);hold on
for i=1:group_num % plot error bars;
   errorbar(i,nanmean(barData1(:,i)),nanstd(barData1(:,i),1)./sqrt(length(find(~isnan(barData1(:,i))))),...
    'vertical','LineStyle','none','LineWidth',2,'Color',colorerrorbar(i,:),'CapSize',18);hold on
end
hold on;ax=gca,ax.Color='none';ax.LineWidth=1.5;ax.TickDir='in';ax.FontSize=14;ax.XTick=1:1:2;box off
ylabel('Total change (AUC)','FontName','Arial','Fontsize',16);box off;hold on,
ax.Title.FontWeight='bold';ax.Title.FontSize=18;ax.YDir='normal';
ax.YTick=0:20:80;ylim([0,80]);
plot([1.1,1.9],[60,60],'Color',[0,0,0],'LineWidth',1.5);hold on,
set(gca,'FontName','Arial','FontSize',18);set(gcf,'color', [1 1 1]);
dr_name = {'w/o LPC','w LPC'};set(ax(1),'XTickLabel',dr_name,'xtick',1:1:2);box off;
text(1.5,64,num2str(p1),'FontSize',20,'HorizontalAlignment','center','FontWeight','normal');
% text(1.5,62,'*','FontSize',24,'HorizontalAlignment','center','FontWeight','normal');
% saveas(gcf,'../results/AUC_LCablation.png');
%% RG ablation;
coloredge=[[0.3647,0.6392,0.6157];[0.6510,0.3333,0.6157]];colorface=[[0.5333,0.7490,0.7216];[0.8078,0.5765,0.7490]];colorerrorbar=coloredge;

A=All.amp_NVC*100;barData1=A(:,[4,3]);% peak;
[p1,h,stats] = signrank(barData1(:,1),barData1(:,2), 'alpha', 0.05, 'tail', 'both');group_num=size(barData1,2);
Y=barData1;X=ones(size(Y)).*(1:size(Y,2));
figure,for i=1:group_num
    barHdl(i)=bar(i,nanmean(barData1(:,i)),0.8,'FaceColor',colorface(i,:),'FaceAlpha',0.5,'LineWidth',2,'EdgeColor',coloredge(i,:));hold on
end
barBaseLineHdl=barHdl(1).BaseLine;barBaseLineHdl.LineWidth=1.2;
hold on,plot(X',Y','Color',[0,0,0,0.1],'Marker','o','LineWidth',1,'MarkerEdgeColor',[1,1,1].*.6,'MarkerSize',5,'LineWidth',.1);hold on
for i=1:group_num
   errorbar(i,nanmean(barData1(:,i)),nanstd(barData1(:,i),1)./sqrt(length(find(~isnan(barData1(:,i))))),...
    'vertical','LineStyle','none','LineWidth',2,'Color',colorerrorbar(i,:),'CapSize',18);hold on
end
hold on;ax=gca,ax.Color='none';ax.LineWidth=1.5;ax.TickDir='in';ax.FontSize=14;ax.XTick=1:1:2;box off
ylabel('Peak (%)','FontName','Arial','Fontsize',16);box off;hold on,
ax.Title.FontWeight='bold';ax.Title.FontSize=18;ax.YDir='normal';
ax.YTick=0:20:80;ylim([0,80]);
plot([1.1,1.9],[60,60],'Color',[0,0,0],'LineWidth',1.5);hold on,
set(gca,'FontName','Arial','FontSize',18);set(gcf,'color', [1 1 1]);
dr_name = {'w/o RG','w RG'};
set(ax(1),'XTickLabel',dr_name,'xtick',1:1:2);box off;
text(1.5,64,num2str(p1),'FontSize',20,'HorizontalAlignment','center','FontWeight','normal');
% text(1.5,62,'***','FontSize',24,'HorizontalAlignment','center','FontWeight','normal');
% saveas(gcf,'../results/Peak_RGablation.png');

A=All.FMHM_NVC;barData1=A(:,[4,3]);% FWHM;
[p1,h,stats] = signrank(barData1(:,1),barData1(:,2), 'alpha', 0.05, 'tail', 'both');group_num=size(barData1,2);
Y=barData1;X=ones(size(Y)).*(1:size(Y,2));
figure,for i=1:group_num
    barHdl(i)=bar(i,nanmean(barData1(:,i)),0.8,'FaceColor',colorface(i,:),'FaceAlpha',0.5,'LineWidth',2,'EdgeColor',coloredge(i,:));hold on
end
barBaseLineHdl=barHdl(1).BaseLine;barBaseLineHdl.LineWidth=1.2;
hold on,plot(X',Y','Color',[0,0,0,0.1],'Marker','o','LineWidth',1,'MarkerEdgeColor',[1,1,1].*.6,'MarkerSize',5,'LineWidth',.1);hold on
for i=1:group_num
   errorbar(i,nanmean(barData1(:,i)),nanstd(barData1(:,i),1)./sqrt(length(find(~isnan(barData1(:,i))))),...
    'vertical','LineStyle','none','LineWidth',2,'Color',colorerrorbar(i,:),'CapSize',18);hold on
end
hold on;ax=gca,ax.Color='none';ax.LineWidth=1.5;ax.TickDir='in';ax.FontSize=14;ax.XTick=1:1:2;box off
ylabel('FWHM (s)','FontName','Arial','Fontsize',16);box off;hold on,
ax.Title.FontWeight='bold';ax.Title.FontSize=18;ax.YDir='normal';
ax.YTick=0:50:200;ylim([0,200]);
plot([1.1,1.9],[150,150],'Color',[0,0,0],'LineWidth',1.5);hold on,
set(gca,'FontName','Arial','FontSize',18);set(gcf,'color', [1 1 1]);
dr_name = {'w/o RG','w RG'};set(ax(1),'XTickLabel',dr_name,'xtick',1:1:2);box off;
text(1.5,160,num2str(p1),'FontSize',20,'HorizontalAlignment','center','FontWeight','normal');
% text(1.5,156,'****','FontSize',24,'HorizontalAlignment','center','FontWeight','normal');
% saveas(gcf,'../results/FWHM_RGablation.png');

A=All.AUC_NVC;barData1=A(:,[4,3]);%  AUC;
[p1,h,stats] = signrank(barData1(:,1),barData1(:,2), 'alpha', 0.05, 'tail', 'both');group_num=size(barData1,2);
Y=barData1;X=ones(size(Y)).*(1:size(Y,2));
figure,for i=1:group_num
    barHdl(i)=bar(i,nanmean(barData1(:,i)),0.8,'FaceColor',colorface(i,:),'FaceAlpha',0.5,'LineWidth',2,'EdgeColor',coloredge(i,:));hold on
end
barBaseLineHdl=barHdl(1).BaseLine;barBaseLineHdl.LineWidth=1.2;
hold on,plot(X',Y','Color',[0,0,0,0.1],'Marker','o','LineWidth',1,'MarkerEdgeColor',[1,1,1].*.6,'MarkerSize',5,'LineWidth',.1);hold on
for i=1:group_num
   errorbar(i,nanmean(barData1(:,i)),nanstd(barData1(:,i),1)./sqrt(length(find(~isnan(barData1(:,i))))),...
    'vertical','LineStyle','none','LineWidth',2,'Color',colorerrorbar(i,:),'CapSize',18);hold on
end
hold on;ax=gca,ax.Color='none';ax.LineWidth=1.5;ax.TickDir='in';ax.FontSize=14;ax.XTick=1:1:2;box off
ylabel('Total flow (AUC)','FontName','Arial','Fontsize',16);box off;hold on,
ax.Title.FontWeight='bold';ax.Title.FontSize=18;ax.YDir='normal';
ax.YTick=0:20:80;ylim([0,80]);
plot([1.1,1.9],[60,60],'Color',[0,0,0],'LineWidth',1.5);hold on,
set(gca,'FontName','Arial','FontSize',18);set(gcf,'color', [1 1 1]);
dr_name = {'w/o RG','w RG'};ylim([-10,70]);set(ax(1),'XTickLabel',dr_name,'xtick',1:1:2);box off;
text(1.5,64,num2str(p1),'FontSize',20,'HorizontalAlignment','center','FontWeight','normal');
% text(1.5,62,'****','FontSize',24,'HorizontalAlignment','center','FontWeight','normal');
% saveas(gcf,'../results/AUC_RGablation.png');
%% Plot of CBF parameters, NVC, NVC + d-LPC, NVC + i-LPC put together, Fig. 7c;
coloredge=[[0.4431,0.4471,0.4784];[0.3647,0.6392,0.6157];[0.6510,0.3333,0.6157]];
colorface=[[0.5686,0.5961,0.6824];[0.5333,0.7490,0.7216];[0.8078,0.5765,0.7490]];
colorerrorbar=coloredge;

A=All.amp_NVC*100;barData1=A(:,[2,4,3]);% peak;
[p1(1),h,stats] = signrank(barData1(:,1),barData1(:,2), 'alpha', 0.05, 'tail', 'both');
[p1(2),h,stats] = signrank(barData1(:,2),barData1(:,3), 'alpha', 0.05, 'tail', 'both');
[p1(3),h,stats] = signrank(barData1(:,1),barData1(:,3), 'alpha', 0.05, 'tail', 'both');
FDR=mafdr(p1,'BHFDR', true);% P correction;
group_num=size(barData1,2);
Y=barData1;X=ones(size(Y)).*(1:size(Y,2));
figure,for i=1:group_num
    barHdl(i)=bar(i,nanmean(barData1(:,i)),0.8,'FaceColor',colorface(i,:),'FaceAlpha',0.5,'LineWidth',2,'EdgeColor',coloredge(i,:));hold on
end
barBaseLineHdl=barHdl(1).BaseLine;barBaseLineHdl.LineWidth=1.2;
hold on,plot(X',Y','Color',[0,0,0,0.1],'Marker','o','LineWidth',1,'MarkerEdgeColor',[1,1,1].*.6,'MarkerSize',5,'LineWidth',.1);hold on
for i=1:group_num% plot error bars;
   errorbar(i,nanmean(barData1(:,i)),nanstd(barData1(:,i),1)./sqrt(length(find(~isnan(barData1(:,i))))),...
    'vertical','LineStyle','none','LineWidth',2,'Color',colorerrorbar(i,:),'CapSize',18);hold on
end
hold on;ax=gca,ax.Color='none';ax.LineWidth=1.5;ax.TickDir='in';ax.FontSize=14;ax.XTick=1:1:2;box off
ylabel('Peak (%)','FontName','Arial','Fontsize',16);box off;hold on,
ax.Title.FontWeight='bold';ax.Title.FontSize=18;ax.YDir='normal';
ax.YTick=0:20:80;ylim([0,80]);
set(gca,'FontName','Arial','FontSize',18);set(gcf,'color', [1 1 1]);
plot([1.1,1.9],[60,60],'Color',[0,0,0],'LineWidth',1.5);hold on,
plot([2.1,2.9],[60,60],'Color',[0,0,0],'LineWidth',1.5);hold on,
text(1.5,64,num2str(p1(1)),'FontSize',16,'HorizontalAlignment','center','FontWeight','normal');
% text(1.5,62,'****','FontSize',24,'HorizontalAlignment','center','FontWeight','normal');
text(2.5,64,num2str(p1(2)),'FontSize',16,'HorizontalAlignment','center','FontWeight','normal');
% text(2.5,62,'***','FontSize',24,'HorizontalAlignment','center','FontWeight','normal');
plot([1.3,2.7],[70,70],'Color',[0,0,0],'LineWidth',1.5);hold on,
text(2,74,num2str(p1(3)),'FontSize',16,'HorizontalAlignment','center','FontWeight','normal');
% text(2,72,'****','FontSize',24,'HorizontalAlignment','center','FontWeight','normal');
dr_name = {'','',''};set(ax(1),'XTickLabel',dr_name,'xtick',1:1:3);box off;
legend('NVC','NVC + d-LPC','NVC + d-LPC + i-LPC','FontSize',14,'Location','NorthWest');
% saveas(gcf,'../results/Peak_NVC_d-LPC_i-LPC.png');

A=All.FMHM_NVC;barData1=A(:,[2,4,3]);% FWHM;
[p1(1),h,stats] = signrank(barData1(:,1),barData1(:,2), 'alpha', 0.05, 'tail', 'both');
[p1(2),h,stats] = signrank(barData1(:,2),barData1(:,3), 'alpha', 0.05, 'tail', 'both');
[p1(3),h,stats] = signrank(barData1(:,1),barData1(:,3), 'alpha', 0.05, 'tail', 'both');
FDR=mafdr(p1,'BHFDR', true);% P correction;
group_num=size(barData1,2);
Y=barData1;X=ones(size(Y)).*(1:size(Y,2));
figure,for i=1:group_num
    barHdl(i)=bar(i,nanmean(barData1(:,i)),0.8,'FaceColor',colorface(i,:),'FaceAlpha',0.5,'LineWidth',2,'EdgeColor',coloredge(i,:));hold on
end
barBaseLineHdl=barHdl(1).BaseLine;barBaseLineHdl.LineWidth=1.2;
hold on,plot(X',Y','Color',[0,0,0,0.1],'Marker','o','LineWidth',1,'MarkerEdgeColor',[1,1,1].*.6,'MarkerSize',5,'LineWidth',.1);hold on
for i=1:group_num% plot error bars;
   errorbar(i,nanmean(barData1(:,i)),nanstd(barData1(:,i),1)./sqrt(length(find(~isnan(barData1(:,i))))),...
    'vertical','LineStyle','none','LineWidth',2,'Color',colorerrorbar(i,:),'CapSize',18);hold on
end
hold on;ax=gca,ax.Color='none';ax.LineWidth=1.5;ax.TickDir='in';ax.FontSize=14;ax.XTick=1:1:2;box off
ylabel('FWHM (s)','FontName','Arial','Fontsize',16);box off;hold on,
ax.Title.FontWeight='bold';ax.Title.FontSize=18;ax.YDir='normal';
ax.YTick=0:50:200;ylim([0,200]);
set(gca,'FontName','Arial','FontSize',18);set(gcf,'color', [1 1 1]);
plot([1.1,1.9],[160,160],'Color',[0,0,0],'LineWidth',1.5);hold on,
plot([2.1,2.9],[160,160],'Color',[0,0,0],'LineWidth',1.5);hold on,
text(1.5,168,num2str(p1(1)),'FontSize',16,'HorizontalAlignment','center','FontWeight','normal');
text(2.5,168,num2str(p1(2)),'FontSize',16,'HorizontalAlignment','center','FontWeight','normal');
plot([1.3,2.7],[180,180],'Color',[0,0,0],'LineWidth',1.5);hold on,
text(2,188,num2str(p1(3)),'FontSize',16,'HorizontalAlignment','center','FontWeight','normal');
dr_name = {'','',''};set(ax(1),'XTickLabel',dr_name,'xtick',1:1:3);box off;
% saveas(gcf,'../results/FWHM_NVC_d-LPC_i-LPC.png');

A=All.AUC_NVC;barData1=A(:,[2,4,3]);% AUC;
[p1(1),h,stats] = signrank(barData1(:,1),barData1(:,2), 'alpha', 0.05, 'tail', 'both');
[p1(2),h,stats] = signrank(barData1(:,2),barData1(:,3), 'alpha', 0.05, 'tail', 'both');
[p1(3),h,stats] = signrank(barData1(:,1),barData1(:,3), 'alpha', 0.05, 'tail', 'both');
FDR=mafdr(p1,'BHFDR', true);% P correction;
group_num=size(barData1,2);
Y=barData1;X=ones(size(Y)).*(1:size(Y,2));
figure,for i=1:group_num
    barHdl(i)=bar(i,nanmean(barData1(:,i)),0.8,'FaceColor',colorface(i,:),'FaceAlpha',0.5,'LineWidth',2,'EdgeColor',coloredge(i,:));hold on
end
barBaseLineHdl=barHdl(1).BaseLine;barBaseLineHdl.LineWidth=1.2;
hold on,plot(X',Y','Color',[0,0,0,0.1],'Marker','o','LineWidth',1,'MarkerEdgeColor',[1,1,1].*.6,'MarkerSize',5,'LineWidth',.1);hold on
for i=1:group_num% plot error bars;
   errorbar(i,nanmean(barData1(:,i)),nanstd(barData1(:,i),1)./sqrt(length(find(~isnan(barData1(:,i))))),...
    'vertical','LineStyle','none','LineWidth',2,'Color',colorerrorbar(i,:),'CapSize',18);hold on
end
hold on;ax=gca,ax.Color='none';ax.LineWidth=1.5;ax.TickDir='in';ax.FontSize=14;ax.XTick=1:1:2;box off
ylabel('Total flow (AUC)','FontName','Arial','Fontsize',16);box off;hold on,
ax.Title.FontWeight='bold';ax.Title.FontSize=18;ax.YDir='normal';
ax.YTick=0:20:80;ylim([-5,80]);
set(gca,'FontName','Arial','FontSize',18);set(gcf,'color', [1 1 1]);
plot([1.1,1.9],[45,45],'Color',[0,0,0],'LineWidth',1.5);hold on,
plot([2.1,2.9],[45,45],'Color',[0,0,0],'LineWidth',1.5);hold on,
text(1.5,49,num2str(p1(1)),'FontSize',16,'HorizontalAlignment','center','FontWeight','normal');
text(2.5,49,num2str(p1(2)),'FontSize',16,'HorizontalAlignment','center','FontWeight','normal');
plot([1.3,2.7],[53,53],'Color',[0,0,0],'LineWidth',1.5);hold on,
text(2,57,num2str(p1(3)),'FontSize',16,'HorizontalAlignment','center','FontWeight','normal');
dr_name = {'','',''};set(ax(1),'XTickLabel',dr_name,'xtick',1:1:3);box off;
% saveas(gcf,'../results/AUC_NVC_d-LPC_i-LPC.png');
%%
% save('../results/CBF_simulation');
save('results_CBF_simulation');
