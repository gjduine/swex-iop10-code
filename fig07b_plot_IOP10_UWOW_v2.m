% script to read and plot UWOW data collected during IOP10 for SWEX 2022
% GJ Duine, ERI UC Santa Barbara, June 2022

% IOP10 had a different approach with UWOW at Montecito fire station 2

clear;clc;tic

% read data file
load('/home/sbarc/duine/sundowners/swex2022/data/uwow/prelim/IOP10_stationary/UWOW_processed/UWOW_IOP10.mat')

timePDTMat=repmat(tnumPDTUwow,[size(wspdUwow,1) 1]);

IOPstart=datenum(2022,5,12,10,0,0);
IOPend=datenum(2022,5,13,10,0,0);

xTicks=[IOPstart:datenum(0,0,0,4,0,0):IOPend];

wspdUwowSNR=wspdUwow;
wspdUwowSNR(SNRUwow(:,:)<1.01)=NaN;

%% start plotting

close
fig=figGD('off');
sp1=subplot_tight(2,1,2,[0.05 0.1]);
cmap=jet(12);
colormap(cmap)

p1=pcolor(timePDTMat(3:end,:),heightUwow(3:end,:),wspdUwowSNR(3:end,:));
set(gca,'box','on'); hold on;
ylabel('Elevation [m]');
cb1=colorbar; shading interp;
% title(['UWOW IOP 10 stationary Montecito fire station 2 :',num2str(latUwow(1,1)),'N ',num2str(lonUwow(1,1)),'E'])
ylim([0 1500])
xlim([IOPstart IOPend])
cb1.Ticks=[0:5:20];
caxis([0 20])
ylabel(cb1,'Wind speed [m s^{-1}]')
set(sp1,'TickDir','out')
% datetick('x','HH:MM','keeplimits'); grid on
set(gca,'TickDir','out','XTick',xTicks,'XTickLabels',datestr(xTicks,'HH:MM'))        
text(IOPstart+0.04,1400,['(b) MFS2 - OBS Lidar'],'FontWeight','bold','BackgroundColor','w')

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20)

xlabel(['Time ',datestr(tnumPDTUwow(1,1),'mmm dd'),' to ',datestr(tnumPDTUwow(1,end),'mmm dd yyyy'),' [hours PDT]'])
set(sp1,'Position',[0.2000 0.1000 0.6500 0.40])


export_fig('plot_IOP10_UWOW_shading_interp_v2','-png')


toc