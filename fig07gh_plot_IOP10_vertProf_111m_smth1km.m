% script to read and plot UWOW data collected during IOP10 for SWEX 2022
% and combine these with WRF model output such that the same locations are 
% shown in the time height plots. 
% In this case, the WRF output was a forecast run during SWEX. The output
% is every 10 minutes and 1 km grid spacing.

% The wrf model output is structured such that the array mataches the array
% from UWOW. Therefore, only need to find the right times ate ach time and
% then reconstruct the model output matrix

% GJ Duine, ERI UC Santa Barbara, June 2022


clear;clc;tic

shdng={'interp','flat'};

% some "parameters" that determine the type of run
runX='run_1km_545_24_z107_2way'; % coarse grid input
IOPnr='IOP10';
tryX='2Dom_smth1km'; % domain settings, specific to case study, likely
doX='d04'; % d04 or d05
res='111m'; % 111m, or XXm
LESnest='ndown'; % type ndown or online
eVertStr='107';
km_opt='km_opt2';
runXstr=strrep(runX,'_','-');

% read data files
% UWOW data
load(['/home/sbarc/duine/sundowners/swex2022/data/uwow/prelim/',IOPnr,'_stationary/UWOW_processed/UWOW_',IOPnr,'.mat'])
% WRF data that has only vertical profiles in the UWOW location 
load(['/home/sbarc/duine/sundowners/swex2022/analysis/analysis_',IOPnr,'/model/era5_10min_v451/timeHeight/UWOW_WRF_ERA5_',IOPnr,'_10min.mat']);
%      /home/sbarc/duine/sundowners/swex2022/analysis/analysis_IOP10/model/era5_10min/timeHeight/wrf_timeHeight_UWOW.m

% WRF locations
folder='/home/sbarc/duine/sundowners/swex2022/data/model/IOP10/';
f1=['swex2022_',IOPnr,'_LES_',runX,'_',doX,'_',LESnest,'_',res,'_',eVertStr,'_',tryX,'_',km_opt,'_dt30.mat'];

% f1=['swex2022_',IOPnr,'_vertProfs_moreLocs_LES_',doX,'_',res,'_',tryX,'_dt30.mat'];
load(strcat(folder,f1));

timePDTMat=repmat(tnumPDTUwow,[size(wspdUwow,1) 1]);

%% start reconstructing a 2-D array from WRF model that matches the UWOW drive locations and timing
% 
% for t=1:length(tnumPDTUwow)
%     [minVal(t),id(t)]=min(abs(tnumPDTUwow(t)-LES.tnumPDT));
% end

% runs={'run_545_ERA5_1km';... % important to keep this as first, as it refers to different folder than others (v451 instead of v460)
%     'run_1km_141_24_z107';'run_1km_242_24_z107';'run_1km_242_24_z107Auto';...
%     'run_1km_540_25_z107';'run_1km_545_24_z107'};
% runsTi={'MYNN v4.5 z-55';...
%     'YSU';'MYJ';'MYJ-zAuto';...
%     'SMS3DTKE';'MYNN v4.6'};
IOPstart=datenum(2022,5,12,10,0,0);
IOPend=datenum(2022,5,13,10,0,0);

xTicks=[IOPstart:datenum(0,0,0,4,0,0):IOPend];

clear('wWRF','wspdWRF','wdirWRF','heightWRF','tnumPDTWRF')
nDims=size(LES.z);

% locations of interest are st=44, st=45, st=48, st=49

stName{19}='MFS-2';
stName{21}='SBBG';
panel{19}='(h)';
panel{21}='(g)';

for st=[19,21]%[2:4,8,10:21]

    for t=1:length(LES.tnumPDT)
        wspdWRF(:,t)=squeeze(LES.wspd(st,:,(t)));
        wdirWRF(:,t)=squeeze(LES.wdir(st,:,(t)));        
        wWRF(:,t)=squeeze(LES.W(st,:,(t)));
        uWRF(:,t)=squeeze(LES.U(st,:,(t)));
        vWRF(:,t)=squeeze(LES.V(st,:,(t)));
        qWRF(:,t)=squeeze(LES.Q(st,:,(t)))*1000;
        bvWRF(:,t)=squeeze(LES.BV(st,:,(t)));
        thetaWRF(:,t)=squeeze(LES.thetaV(st,:,(t)));
        heightWRF(:,t)=squeeze(LES.z(st,:,(t)));
        tnumPDTWRF(:,t)=repmat(LES.tnumPDT(t),[nDims(2) 1]);
    end

    wspdUwowSNR=wspdUwow;
    wspdUwowSNR(SNRUwow(:,:)<1.01)=NaN;

    wdirUwowSNR=wdirUwow;
    wdirUwowSNR(SNRUwow(:,:)<1.01)=NaN;

    thetaMin=floor(min(min(thetaWRF)));
    thetaMax=ceil(max(max(thetaWRF)));

    %% start plotting. First wind speed relative to ground level


    titleStr=[IOPnr,' ',runXstr,' DX:',res,'; km opt',num2str(LES.km_opt),'; ',LES.obs_stations{st},' ',num2str(LES.hgt_st(st),4),' m'];
    textStr=[panel{st}, ' WRF - ',res, ' smth1km ',stName{st}];

    % select shading
    sh=1;
    heightUwowMSL=heightUwow+terElevUwow;
    
    %%%%% brunt vaisala plus pbl depth %%%%%%
        close
    fig=figGD('off');

    % heightWRFUwowMSL=heightWRFUwow+LES.hgt_st(st)';
    sp2=subplot_tight(2,1,2,[0.05 0.1]);
    cmap=flipud(gray(30));
    colormap(cmap)
    clims2=[-0.002 0.03];
    c1=contourf(tnumPDTWRF,heightWRF,real(bvWRF),[clims2(1):0.002:clims2(2)],'edgecolor','none'); % first fill color
    set(gca,'box','on'); hold on;

    hold on
    %         cmap2=flipud(hot(30));
    % cmap2=flipud(gray(30));
    colormap(sp2,cmap);   caxis(clims2);   set(sp2,'CLim',clims2)

    % pblh
    line(tnumPDTWRF(1,:),LES.PBLH_st(st,:),'color','r','LineWidth',3)
    line([tnumPDTWRF(1,1) tnumPDTWRF(1,end)],...
        [LES.zFroude_st(st) LES.zFroude_st(st)],...
        'color','b','LineWidth',3,'LineStyle','--')
    % line([WRF.(runX).tnumPST(indSt) WRF.(runX).tnumPST(indEnd)],...
    %     [WRF.(runX).hgt_st(st) WRF.(runX).hgt_st(st)],...
    %     'color','k','LineWidth',1,'LineStyle',':') % terrain elevation

    % p2a=line(tnumPDTUwow,repmat(LES.hgt_st(st),[length(tnumPDTUwow) 1]),'LineWidth',2,'color','k');
    ylabel('Elevation [m]');
    cb2=colorbar; shading(shdng{sh});
    % title(titleStr)
    ylim([0 1500])
    xlim([IOPstart IOPend])
    % caxis([0 20])
    ylabel(cb2,'N [s^{-1}]')
    % datetick('x','HH:MM','keeplimits'); grid on
    set(gca,'TickDir','out','XTick',xTicks,'XTickLabels',datestr(xTicks,'HH:MM'))        
    text(IOPstart+0.04,1400,textStr,'FontWeight','bold','BackgroundColor','w')
    xlabel(['Time ',datestr(tnumPDTUwow(1,1),'mmm dd'),' to ',datestr(tnumPDTUwow(1,end),'mmm dd yyyy'),' [hours PDT]'])
    % set(sp1,'Position',[0.2000 0.5500 0.6500 0.40])
    set(sp2,'Position',[0.2000 0.1000 0.6500 0.40])

    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20)

    export_fig(['plot_',IOPnr,'_vertProf_',res,'_smth1km_N_',LES.obs_stations{st},'_shading_',shdng{sh}],'-png')





    %%%%% wind speed and theta %%%%%

    close
    fig=figGD('off');


    % heightWRFUwowMSL=heightWRFUwow+LES.hgt_st(st)';
    sp2=subplot_tight(2,1,2,[0.05 0.1]);
    cmap=jet(12);
    colormap(cmap)

    p2=pcolor(tnumPDTWRF,heightWRF,wspdWRF);
    set(gca,'box','on'); hold on;
    [c,h]=contour(tnumPDTWRF,heightWRF,thetaWRF,...
        [thetaMin:3:thetaMax],'k','Linewidth',2,'ShowText','on');

    contour(tnumPDTWRF,heightWRF,thetaWRF,...
        [thetaMin:1:thetaMax],'k','Linewidth',1,'ShowText','off')

    % windbarbm_NH(tnumPDTWRF,heightWRF,vWRF,vWRF,...
    %     2.5,'color','k','LineWidth',2);

    % p2a=line(tnumPDTUwow,repmat(LES.hgt_st(st),[length(tnumPDTUwow) 1]),'LineWidth',2,'color','k');
    ylabel('Elevation [m]');
    cb2=colorbar; shading(shdng{sh});
    % title(titleStr)
    ylim([0 1500])
    xlim([IOPstart IOPend])
    caxis([0 20])
    ylabel(cb2,'Wind speed [m s^{-1}]')%,'Ticks',[0:90:360],'TickLabels',{'N','E','S','W','N'})
    % set(gca,'TickDir','out','XTickLabel','')
    % datetick('x','HH:MM','keeplimits'); grid on
    set(gca,'TickDir','out','XTick',xTicks,'XTickLabels',datestr(xTicks,'HH:MM'))        
    text(IOPstart+0.04,1400,textStr,'FontWeight','bold','BackgroundColor','w')

    xlabel(['Time ',datestr(tnumPDTUwow(1,1),'mmm dd'),' to ',datestr(tnumPDTUwow(1,end),'mmm dd yyyy'),' [hours PDT]'])
    % set(sp1,'Position',[0.2000 0.5500 0.6500 0.40])
    set(sp2,'Position',[0.2000 0.1000 0.6500 0.40])

    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20)
    clabel(c,h,'FontSize',16)

    export_fig(['plot_',IOPnr,'_vertProf_',res,'_smth1km_wspd_theta_',LES.obs_stations{st},'_shading_',shdng{sh}],'-png')



    %%%%% wind speed only %%%%%

    close
    fig=figGD('off');


    % heightWRFUwowMSL=heightWRFUwow+LES.hgt_st(st)';
    sp2=subplot_tight(2,1,2,[0.05 0.1]);
    cmap=jet(12);
    colormap(cmap)

    p2=pcolor(tnumPDTWRF,heightWRF,wspdWRF);
    set(gca,'box','on'); hold on;
    % p2a=line(tnumPDTUwow,repmat(LES.hgt_st(st),[length(tnumPDTUwow) 1]),'LineWidth',2,'color','k');
    ylabel('Elevation [m]');
    cb2=colorbar; shading(shdng{sh});
    % title(titleStr)
    ylim([0 1500])
    xlim([IOPstart IOPend])
    caxis([0 20])
    ylabel(cb2,'Wind speed [m s^{-1}]')%,'Ticks',[0:90:360],'TickLabels',{'N','E','S','W','N'})
    % set(gca,'TickDir','out','XTickLabel','')
    % datetick('x','HH:MM','keeplimits'); grid on
    set(gca,'TickDir','out','XTick',xTicks,'XTickLabels',datestr(xTicks,'HH:MM'))        
    text(IOPstart+0.04,1400,textStr,'FontWeight','bold','BackgroundColor','w')

    xlabel(['Time ',datestr(tnumPDTUwow(1,1),'mmm dd'),' to ',datestr(tnumPDTUwow(1,end),'mmm dd yyyy'),' [hours PDT]'])
    % set(sp1,'Position',[0.2000 0.5500 0.6500 0.40])
    set(sp2,'Position',[0.2000 0.1000 0.6500 0.40])

    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20)


    export_fig(['plot_',IOPnr,'_vertProf_',res,'_smth1km_wspd_',LES.obs_stations{st},'_shading_',shdng{sh}],'-png')


    %% %%%% Q-vapor %%%%%%

    close
    fig=figGD('off');


    % heightWRFUwowMSL=heightWRFUwow+WRF.hgt_st(st)';
    sp3=subplot_tight(2,1,2,[0.05 0.1]);
    cmap=flipud(jet(30));
    colormap(cmap);
    p2=pcolor(tnumPDTWRF,heightWRF,qWRF);
    set(gca,'box','on'); hold on;
    % p2a=line(tnumPDTUwow,repmat(WRF.hgt_st(st),[length(tnumPDTUwow) 1]),'LineWidth',2,'color','k');
    ylabel('Elevation [m]');
    cb2=colorbar; shading(shdng{sh});
    % title(titleStr)

    ylim([0 1500])
    xlim([IOPstart IOPend])
    caxis([1 5])
    ylabel(cb2,'Mixing ratio [g kg^{-1}]')%,'Ticks',[0:90:360],'TickLabels',{'N','E','S','W','N'})
    set(gca,'TickDir','out','XTickLabel','')
    datetick('x','HH:MM','keeplimits'); grid on
    text(IOPstart+0.04,1400,textStr,'FontWeight','bold','BackgroundColor','w')

    xlabel(['Time ',datestr(tnumPDTUwow(1,1),'mmm dd'),' to ',datestr(tnumPDTUwow(1,end),'mmm dd yyyy'),' [hours PDT]'])
    % set(sp1,'Position',[0.2000 0.5500 0.6500 0.40])
    set(sp3,'Position',[0.2000 0.1000 0.6500 0.40])

    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20)


    export_fig(['plot_',IOPnr,'_vertProf_',res,'_smth1km_Q_',LES.obs_stations{st},'_shading_',shdng{sh}],'-png')


    %% wind direction with terrain elevation

    close
    fig=figGD('off');


    sp2=subplot_tight(2,1,2,[0.05 0.1]);
    colormap('jet');
    p2=pcolor(tnumPDTWRF,heightWRF,wdirWRF);
    set(gca,'box','on'); hold on;
    % p2a=line(tnumPDTUwow,LES.hgt_st(st),'LineWidth',2,'color','k');
    ylabel('Elevation [m]');
    cb2=colorbar; shading(shdng{sh});
    % title(titleStr)
    ylim([0 1500])
    xlim([IOPstart IOPend])
    caxis([0 360])
    ylabel(cb2,'Wind direction [^\circN]')%,'Ticks',[0:90:360],'TickLabels',{'N','E','S','W','N'})
    cb2.Ticks=[0:90:360];
    cb2.TickLabels={'N','E','S','W','N'};
    % set(gca,'TickDir','out','XTickLabel','')
    % datetick('x','HH:MM','keeplimits'); grid on
    set(gca,'TickDir','out','XTick',xTicks,'XTickLabels',datestr(xTicks,'HH:MM'))        
    text(IOPstart+0.04,1400,textStr,'FontWeight','bold','BackgroundColor','w')

    xlabel(['Time ',datestr(tnumPDTUwow(1,1),'mmm dd'),' to ',datestr(tnumPDTUwow(1,end),'mmm dd yyyy'),' [hours PDT]'])
    % set(sp1,'Position',[0.2000 0.5500 0.6500 0.40])
    set(sp2,'Position',[0.2000 0.1000 0.6500 0.40])

    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20)

    export_fig(['plot_',IOPnr,'_vertProf_',res,'_smth1km_wdir_',LES.obs_stations{st},'_shading_',shdng{sh}],'-png')




    %% W-component

    close
    fig=figGD('off');


    sp2=subplot_tight(2,1,2,[0.05 0.1]);
    p2=pcolor(tnumPDTWRF,heightWRF,wWRF);
    set(gca,'box','on'); hold on;
    bwr=flipud(bluewhitered(30));
    colormap(bwr);    % p2a=line(tnumPDTUwow,LES.hgt_st(st),'LineWidth',2,'color','k');
    ylabel('Elevation [m]');
    cb2=colorbar; shading(shdng{sh});
    % title(titleStr)
    ylim([0 1500])
    xlim([IOPstart IOPend])
    set(gca,'CLim',[-5 5])
    caxis([-5 5])

    hold on
    ylabel(cb2,'W-component [m s^{-1}]')%,'Ticks',[0:90:360],'TickLabels',{'N','E','S','W','N'})
    cb2.Ticks=[-5:1:5];
    % cb2.TickLabels={'N','E','S','W','N'};
    % set(gca,'TickDir','out','XTickLabel','')
    % datetick('x','HH:MM','keeplimits'); grid on
    set(gca,'TickDir','out','XTick',xTicks,'XTickLabels',datestr(xTicks,'HH:MM'))        
    text(IOPstart+0.04,1400,textStr,'FontWeight','bold','BackgroundColor','w')

    % need to set it a second time
    bwr=flipud(bluewhitered(30));
    colormap(bwr);
    caxis([-5 5])
    box on
    xlabel(['Time ',datestr(tnumPDTUwow(1,1),'mmm dd'),' to ',datestr(tnumPDTUwow(1,end),'mmm dd yyyy'),' [hours PDT]'])
    % set(sp1,'Position',[0.2000 0.5500 0.6500 0.40])
    set(sp2,'Position',[0.2000 0.1000 0.6500 0.40])

    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20)

    export_fig(['plot_',IOPnr,'_vertProf_',res,'_smth1km_Wcomp_',LES.obs_stations{st},'_shading_',shdng{sh}],'-png')


end % end of station loop




toc