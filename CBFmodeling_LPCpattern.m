close all;clc;clear;
%% Load workspace;
% Load NVC and non-NVC vessel diameter dynamics, ID of vessels, for CBF simulation of LPC pattern;
load('2024_workspace_CBFmodeling_LPCpattern.mat');
% load('../data/20240319_workspace_CBFmodeling_LPCpattern.mat');
%% LPC kernel;shift of cap to art latency;
t=1:1:250;
% align all (arteriole) vessel diameter to the same constriction onset;
timealign_kernel=modigliaintact_kernel;
artregion_id=[2,4,6];
art_delay=zeros(1,6)*NaN;art_delay(2)=28;art_delay(4)=6;art_delay(6)=38;
snr=50;px_dBW=0;
for k=1:length(artregion_id)
    timealign_kernel([51:length(t)-art_delay(artregion_id(k))],artregion_id(k))=timealign_kernel([51+art_delay(artregion_id(k)):length(t)],artregion_id(k));
    for k=3
    timealign_kernel([length(t)-art_delay(artregion_id(k))+1:length(t)],artregion_id(k))=awgn(zeros(art_delay(artregion_id(k)),1),snr,px_dBW);
    timealign_kernel(:,artregion_id(k))=smoothdata(timealign_kernel(:,artregion_id(k)),'gaussian',15);
    end
end
figure,plot(timealign_kernel);
% shift of Cap. to Art. latency;
snr=30;px_dBW=0;winshift=5;shiftnum=13;
disrupted_kernel=cell(1,shiftnum);disrupted_kernel{1}=timealign_kernel;
for ii=2:shiftnum
    disrupted_kernel{ii}=timealign_kernel;
    shiftonset=50+winshift*(ii-1);
    for k=1:length(artregion_id)
    disrupted_kernel{ii}([51:shiftonset],artregion_id(k))=awgn(ones(winshift*(ii-1),1)*timealign_kernel(30,artregion_id(k)),snr,px_dBW);
    disrupted_kernel{ii}([shiftonset+1:length(t)],artregion_id(k))=timealign_kernel([51:length(t)-winshift*(ii-1)],artregion_id(k));
    end
    disrupted_kernel{ii}=smoothdata(disrupted_kernel{ii},'gaussian',10);
end
% Plot of disrupted LPC vessel diameter dynamics;
for i=1:13
    figure,plot(disrupted_kernel{i});
end
%% LPC vessel dynamics, disrupted latency;
snr=50;px_dBW=0;
t=1:1:250;raw_radii=repmat(vessel_radii,[length(t),1]);
sim_radii=cell(1,shiftnum);
for k=1:shiftnum
    sim_radii{k}=raw_radii;
    for j=1:length(Global_LPcoupling)
    for i=1:length(Global_LPcoupling{j})
    sim_radii{k}(:,Global_LPcoupling{j}(i))=vessel_radii(Global_LPcoupling{j}(i))*(1+0.01*disrupted_kernel{k}(:,j));          
    end
    end
    for m=1:length(vessel_radii)
    A=awgn(sim_radii{k}(:,m),snr,px_dBW);% add noise;
    sim_radii{k}(:,m)=smoothdata(A,'gaussian',10);    
    end
end
Global_sim_radii=sim_radii;
%% CBF simulation, in the presence of LPC of different Cap. to Art. latency;
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

% NVC blood flow dynamics in the presence of LPC;
snr=50;px_dBW=0;% SNR and power of Gaussian noise;
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
% CBF simulation;
[control_Local_flow_eachvessel,control_Local_flow,control_Local_pressure]=bfmodel(control_Local_sim_radii);
control_Local_dfdf_flow=dfvsfo(abs(control_Local_flow_eachvessel));
% response traces of blood flow;
control_Local_dfdf_flow_NVC=control_Local_dfdf_flow(:,NVC_idnew);

% Local plus global,vessel dynamics;
Localplusglobal=[];
for ii=1:shiftnum
    Localplusglobal{ii}=Global_sim_radii{ii};
    for i=1:length(NVC_arteriole_idnew)
    Localplusglobal{ii}(:,NVC_arteriole_idnew(i))=control_Local_sim_radii(:,NVC_arteriole_idnew(i));
    end
    for i=1:length(NVC_capillary_idnew)
    Localplusglobal{ii}(:,NVC_capillary_idnew(i))=control_Local_sim_radii(:,NVC_capillary_idnew(i));
    end    
end

datanum=length(Localplusglobal);
% CBF simulation;
Localplusglobal_flow_eachvessel=cell(1,datanum);Localplusglobal_flow=cell(1,datanum);Localplusglobal_pressure=cell(1,datanum);
Localplusglobal_dfdf_flow=cell(1,datanum);tic
for k=1:datanum
    [Localplusglobal_flow_eachvessel{k},Localplusglobal_flow{k},Localplusglobal_pressure{k}]=bfmodel(Localplusglobal{k});
    Localplusglobal_dfdf_flow{k}=dfvsfo(abs(Localplusglobal_flow_eachvessel{k}));k
end
toc

% Calculation of blood flow parameters;
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
simtrialnum=length(para);% number of simulation trials;
All.control_Local_dfdf_flow_NVC=[];
for j=1:simtrialnum
    All.control_Local_dfdf_flow_NVC=[All.control_Local_dfdf_flow_NVC,para(j).control_Local_dfdf_flow_NVC];
end
All.Localplusglobal_dfdf_flow_NVC=cell(1,datanum);
for i=1:datanum
for j=1:simtrialnum
%    All.control_Local_dfdf_flow{i}=[All.control_Local_dfdf_flow{i},para(j).control_Local_dfdf_flow_NVC];
   All.Localplusglobal_dfdf_flow_NVC{i}=[All.Localplusglobal_dfdf_flow_NVC{i},para(j).Localplusglobal_dfdf_flow_NVC{i}];   
end
end
All.amp_NVC=[];All.FMHM_NVC=[];All.AUC_NVC=[];
for i=1:simtrialnum
    All.amp_NVC=cat(1,All.amp_NVC,para(i).amp_NVC);
    All.FMHM_NVC=cat(1,All.FMHM_NVC,para(i).FMHM_NVC);
    All.AUC_NVC=cat(1,All.AUC_NVC,para(i).AUC_NVC);
end
% plot of blood flow dynamics;
trace=All.Localplusglobal_dfdf_flow_NVC;
for i=1:length(trace)
   trace{i}=trace{i}*100; 
end
mean_data=[];std_data=[];sem_data=[];
for i=1:datanum
   mean_data(:,i)=mean(trace{i}');
   std_data(:,i)=std(trace{i}');
   sem_data(:,i)=std_data(:,i)/sqrt(size(trace{i},2));
end
% gradient=slanCM('magma');
colorlist=gradient;colorratio=16;
t=-49:200;figure,
for i=1:datanum
   shadedErrorBar(t,mean_data(:,i)',sem_data(:,i),'lineprops',{'color',colorlist((i+1)*colorratio,:)},'transparent',1,'patchSaturation', 0.05);hold on;    
end
hold on,for i=1:datanum
   p(i)=plot(t,mean_data(:,i),'linewidth',1,'color',colorlist(i*colorratio,:),'linestyle','-','Marker','none');hold on
end
ylim([-5 30]);ax=gca,ax.YTick=-20:10:30;%xlim([0 250]);
hold on;plot([0,0],[ax.YLim(1),ax.YLim(2)],'linewidth',0.5,'color',[0.1882,0.1882,0.1882],'linestyle','--','Marker','none');
xlim([-50 200]);xlabel( 'Time (s) ','FontSize',18,'FontName','Arial');ylabel('Blood flow (%)','FontSize',18,'FontName','Arial');
set(gca,'FontName','Arial','FontSize',20);set(gcf,'color', [1 1 1]);box off;set(gca,'Position', [.15 .18 .8 .75]); 
% saveas(gcf,'../results/Local_CBF_change_LPCpattern_trace.png');

figure,clim=[0 60];colormap(colorlist([1:colorratio:colorratio*shiftnum],:));% colorbar;
h=colorbar();caxis(clim);set(h,'FontSize',18);
%% Summary data, NVC flow in the presence of LPC with different Cap. to Art. latency;
A=All.amp_NVC*100;barData1=A;% peak; change to percent;
p1=[];
[p1(1),h,stats] = signrank(barData1(:,1),barData1(:,7), 'alpha', 0.05, 'tail', 'both');
[p1(2),h,stats] = signrank(barData1(:,7),barData1(:,13), 'alpha', 0.05, 'tail', 'both');
FDR=mafdr(p1,'BHFDR', true);% P correction;
group_num=size(barData1,2);
figure,for i=1:group_num
    barHdl(i)=bar(i,nanmean(barData1(:,i)),'FaceColor',gradient(i*colorratio+40,:),'FaceAlpha',0.8,'LineWidth',2,'EdgeColor',gradient(i*colorratio,:));hold on
end
barBaseLineHdl=barHdl(1).BaseLine;barBaseLineHdl.LineWidth=1.2;
for i=1:group_num
    scatter(i+rand(1,length(find(~isnan(barData1(:,i))))).*.6-.3,barData1(find(~isnan(barData1(:,i))),i),15,'CData',gradient(i*colorratio,:),'linewidth',0.5,'MarkerEdgeAlpha',.8);hold on
end
for i=1:group_num
   errorbar(i,nanmean(barData1(:,i)),nanstd(barData1(:,i),1)./sqrt(length(find(~isnan(barData1(:,i))))),...
    'vertical','LineStyle','none','LineWidth',2,'Color',gradient(i*colorratio,:),'CapSize',12);hold on
end
hold on;ax=gca,ax.Color='none';
ax.LineWidth=1.5;ax.TickDir='in';
ax.FontSize=12;ax.XTick=1:1:group_num;
ax.Title.FontWeight='bold';ax.Title.FontSize=18;
ax.YDir='normal';shiftlatency=[0:12]*5;
dr_name = [];for i=1:shiftnum
   dr_name{i} = num2str(shiftlatency(i));
end
set(ax(1),'XTickLabel',dr_name,'xtick',1:1:group_num);box off;xtickangle(60);
set(gca,'FontName','Arial','FontSize',20);box off;set(gcf,'color', [1 1 1]);set(gca,'Position', [.15 .18 .8 .75]); 
ax=gca,ax.YTick=0:20:100;% ax.Title.String='Pericyte response ratio';
ax.YLim=[0,80];%truncAxis('y',[33;52]);
ylabel('Peak CBF increase (%)','FontName','Arial','Fontsize',16);box off %在历次event中，有多大概率有响应;
plot([1.1,6.9],[60,60],'Color',[0,0,0],'LineWidth',1);hold on,
text(4,64,num2str(p1(1)),'FontSize',16,'HorizontalAlignment','center','FontWeight','normal');
plot([7.1,12.9],[60,60],'Color',[0,0,0],'LineWidth',1);hold on,
text(10,64,num2str(p1(2)),'FontSize',16,'HorizontalAlignment','center','FontWeight','normal');
set(gca,'FontName','Arial','FontSize',20);box off;set(gcf,'color', [1 1 1]);set(gca,'Position', [.15 .18 .8 .75]); 
% saveas(gcf,'../results/Local_CBF_change_LPCpattern_peak.png');

A=All.FMHM_NVC;barData1=A;% FWHM;
p1=[];
[p1(1),h,stats] = signrank(barData1(:,1),barData1(:,7), 'alpha', 0.05, 'tail', 'both');
[p1(2),h,stats] = signrank(barData1(:,7),barData1(:,13), 'alpha', 0.05, 'tail', 'both');
FDR=mafdr(p1,'BHFDR', true);% P correction;
group_num=size(barData1,2);
figure,for i=1:group_num
    barHdl(i)=bar(i,nanmean(barData1(:,i)),'FaceColor',gradient(i*colorratio+40,:),'FaceAlpha',0.8,'LineWidth',2,'EdgeColor',gradient(i*colorratio,:));hold on
end
barBaseLineHdl=barHdl(1).BaseLine;barBaseLineHdl.LineWidth=1.2;
for i=1:group_num
    scatter(i+rand(1,length(find(~isnan(barData1(:,i))))).*.6-.3,barData1(find(~isnan(barData1(:,i))),i),15,'CData',gradient(i*colorratio,:),'linewidth',0.5,'MarkerEdgeAlpha',.8);hold on
end
for i=1:group_num
   errorbar(i,nanmean(barData1(:,i)),nanstd(barData1(:,i),1)./sqrt(length(find(~isnan(barData1(:,i))))),...
    'vertical','LineStyle','none','LineWidth',2,'Color',gradient(i*colorratio,:),'CapSize',12);hold on
end
hold on;ax=gca,ax.Color='none';
ax.LineWidth=1.5;ax.TickDir='in';
ax.FontSize=12;ax.XTick=1:1:group_num;
ax.Title.FontWeight='bold';ax.Title.FontSize=18;
ax.YDir='normal';shiftlatency=[0:12]*5;
dr_name = [];for i=1:shiftnum
   dr_name{i} = num2str(shiftlatency(i));
end
set(ax(1),'XTickLabel',dr_name,'xtick',1:1:group_num);box off;xtickangle(60);
set(gca,'FontName','Arial','FontSize',20);box off;set(gcf,'color', [1 1 1]);set(gca,'Position', [.15 .18 .8 .75]); 
ax.YTick=0:50:200;% ax.Title.String='Pericyte response ratio';
ax.YLim=[0,200];
ylabel('FWHM (s)','FontName','Arial','Fontsize',16);box off %在历次event中，有多大概率有响应;
plot([1.1,6.9],[160,160],'Color',[0,0,0],'LineWidth',1);hold on,
text(4,168,num2str(p1(1)),'FontSize',16,'HorizontalAlignment','center','FontWeight','normal');
plot([7.1,12.9],[160,160],'Color',[0,0,0],'LineWidth',1);hold on,
text(10,168,num2str(p1(2)),'FontSize',16,'HorizontalAlignment','center','FontWeight','normal');
set(gca,'FontName','Arial','FontSize',20);box off;set(gcf,'color', [1 1 1]);set(gca,'Position', [.15 .18 .8 .75]); 
% saveas(gcf,'../results/Local_CBF_change_LPCpattern_FWHM.png');

A=All.AUC_NVC;barData1=A;% AUC;
p1=[];
[p1(1),h,stats] = signrank(barData1(:,1),barData1(:,7), 'alpha', 0.05, 'tail', 'both');
[p1(2),h,stats] = signrank(barData1(:,7),barData1(:,13), 'alpha', 0.05, 'tail', 'both');
FDR=mafdr(p1,'BHFDR', true);% P correction;
group_num=size(barData1,2);
figure,for i=1:group_num
    barHdl(i)=bar(i,nanmean(barData1(:,i)),'FaceColor',gradient(i*colorratio+40,:),'FaceAlpha',0.8,'LineWidth',2,'EdgeColor',gradient(i*colorratio,:));hold on
end
barBaseLineHdl=barHdl(1).BaseLine;barBaseLineHdl.LineWidth=1.2;
for i=1:group_num
    scatter(i+rand(1,length(find(~isnan(barData1(:,i))))).*.6-.3,barData1(find(~isnan(barData1(:,i))),i),15,'CData',gradient(i*colorratio,:),'linewidth',0.5,'MarkerEdgeAlpha',.8);hold on
end
for i=1:group_num
   errorbar(i,nanmean(barData1(:,i)),nanstd(barData1(:,i),1)./sqrt(length(find(~isnan(barData1(:,i))))),...
    'vertical','LineStyle','none','LineWidth',2,'Color',gradient(i*colorratio,:),'CapSize',12);hold on
end
hold on;ax=gca,ax.Color='none';
ax.LineWidth=1.5;ax.TickDir='in';
ax.FontSize=12;ax.XTick=1:1:group_num;
ax.Title.FontWeight='bold';ax.Title.FontSize=18;
ax.YDir='normal';shiftlatency=[0:12]*5;
dr_name = [];for i=1:shiftnum
   dr_name{i} = num2str(shiftlatency(i));
end
set(ax(1),'XTickLabel',dr_name,'xtick',1:1:group_num);box off;xtickangle(60);
set(gca,'FontName','Arial','FontSize',20);box off;set(gcf,'color', [1 1 1]);set(gca,'Position', [.15 .18 .8 .75]); 
set(gca,'FontName','Arial','FontSize',20);box off;set(gcf,'color', [1 1 1]);set(gca,'Position', [.15 .18 .8 .75]); 
ax.YTick=0:20:80;ax.YLim=[-5,70];% ax.Title.String='Pericyte response ratio';
ylabel('Total flow (AUC)','FontName','Arial','Fontsize',16);box off %在历次event中，有多大概率有响应;
plot([1.1,6.9],[35,35],'Color',[0,0,0],'LineWidth',1);hold on,
text(4,38,num2str(p1(1)),'FontSize',16,'HorizontalAlignment','center','FontWeight','normal');
plot([7.1,12.9],[35,35],'Color',[0,0,0],'LineWidth',1);hold on,
text(10,38,num2str(p1(2)),'FontSize',16,'HorizontalAlignment','center','FontWeight','normal');
set(gca,'FontName','Arial','FontSize',20);box off;set(gcf,'color', [1 1 1]);set(gca,'Position', [.15 .18 .8 .75]); 
% saveas(gcf,'../results/Local_CBF_change_LPCpattern_AUC.png');
%% 
% save('../results/CBF_simulation_LPCpattern');
save('results_CBF_simulation_LPCpattern');




