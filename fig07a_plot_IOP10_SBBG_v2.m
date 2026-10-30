% script to read and plot SBBG data collected during IOP10 for SWEX 2022
% GJ Duine, ERI UC Santa Barbara, June 2022

clear;clc;tic

% read data file
% load('Z:\SBARC\duine\sundowners\swex2022\data\sbbg\lidar\SBBG_VAD_116_IOP10.mat')
load('/home/voyager-sbarc/duine/sundowners/swex2022/data/sbbg/lidar/SBBG_VAD_116_IOP10.mat')
timePDTMat=repmat(tnumPDTSBBG,[size(wspdSBBG,1) 1]);

IOPstart=datenum(2022,5,12,10,0,0);
IOPend=datenum(2022,5,13,10,0,0);

xTicks=[IOPstart:datenum(0,0,0,4,0,0):IOPend];

% wspdSBBGSNR=wspdSBBG;
% wspdSBBGSNR(SNRSBBG(:,:)<1.01)=NaN;
wspdSBBG(wspdSBBG>12 & heightSBBG>1000) =NaN; % to get rid of bade values abover 1000 m



%% start plotting

close
fig=figGD('off');
sp1=subplot_tight(2,1,2,[0.05 0.1]);
cmap=jet(12);
colormap(cmap)

p1=pcolor(timePDTMat(2:end,:),heightSBBG(2:end,:),wspdSBBG(2:end,:));
set(gca,'box','on'); hold on;
ylabel('Elevation [m]');
cb1=colorbar; shading interp;
% title(['SBBG IOP 10 stationary Montecito fire station 2 :',num2str(latSBBG(1,1)),'N ',num2str(lonSBBG(1,1)),'E'])
ylim([0 1500])
xlim([IOPstart IOPend])
cb1.Ticks=[0:5:20];
caxis([0 20])
ylabel(cb1,'Wind speed [m s^{-1}]')
set(sp1,'TickDir','out')
% datetick('x','HH:MM','keeplimits'); grid on
set(gca,'TickDir','out','XTick',xTicks,'XTickLabels',datestr(xTicks,'HH:MM'))        
text(IOPstart+0.04,1400,['(a) SBBG - OBS Lidar'],'FontWeight','bold','BackgroundColor','w')

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20)

xlabel(['Time ',datestr(tnumPDTSBBG(1,1),'mmm dd'),' to ',datestr(tnumPDTSBBG(1,end),'mmm dd yyyy'),' [hours PDT]'])
set(sp1,'Position',[0.2000 0.1000 0.6500 0.40])


export_fig('fig07a_plot_IOP10_SBBG_shading_interp_v2','-png')


% first period strong winds
datestr(timePDTMat(1,117:123))
wspdSBBG(1:10,117:123)
heightSBBG(1:10,117:123)

% second period of the strong winds (00:37 to 03:07)
datestr(timePDTMat(1,128:138))


toc