% script to SWEX observations and compare with model output. 
% This script will find for each IOP and output in the respective IOP and
% EOP folders

% GJ. Duine, January 2023, ERI UCSB, 

clear;clc;tic

% load ISFS data
folder='/home/sbarc/duine/sundowners/swex2022/data/ISFS/';
f1='SWEX_ISFS_5min_v10_20220401_20220515'; 
% f1='SWEX_ISFS_v01_20220331_20220519'; % this was the preliminary data
load(strcat(folder,f1));
clear('folder','f1')

% load ISS data
folder='/home/sbarc/duine/sundowners/swex2022/data/ISS/sfc_stations/';
f1='SWEX_ISS_v01_20220316_20220515';
load(strcat(folder,f1));
clear('folder','f1')

% load RAWS and NWS data
folder='/home/sbarc/duine/sundowners/swex2022/data_full/obs/';
f1='Stations_SB_20220401_20220515_PDT_lvl1.mat';
load(strcat(folder,f1));
clear('folder','f1')

% load radisounding data
folder='/home/duine/reports/articles/IOP10_LES/figs_progress/radiosonde_MFS2/';
f1='RS_allProfiles_lvl0_v2.mat';
load(strcat(folder,f1));
clear('folder','f1')

% load forecast. This is to select for stations first which is why there is
% IOP01 in the filename. The model data for each IOP data will be loaded in
% a loop
% folder='/home/sbarc/duine/sundowners/swex2022/data/model/';
% f1='WRF_fcst_SWEX2022_vertProfs_moreLocs_IOP01.mat';
% load(strcat(folder,f1));
% WRFfcst=WRF;
% clear('WRF')

% see comments forecast (above)

% 
folder='/home/sbarc/duine/sundowners/swex2022/data/model/IOP10/';
f1='WRF_SWEX2022_IOP10_sensitivity_vertProfs_moreLocs.mat';
load(strcat(folder,f1));
WRF1km=WRF;
clear('WRF')


% load file with sensitivity
folder='/home/sbarc/duine/sundowners/swex2022/data/model/IOP10/';
f1='WRF_SWEX2022_IOP10_sensitivity_vertProfs_moreLocs.mat';
load(strcat(folder,f1));
WRFsens=WRF;
clear('WRF')

folder='/home/sbarc/duine/sundowners/swex2022/data/model/IOP10/';
f1='swex2022_IOP10_LES_run_1km_540_25_z107_1way_d04_ndown_333m_107_try4_km_opt5.mat';
load(strcat(folder,f1));
LES333m=LES;
clear('LES')



folder='/home/sbarc/duine/sundowners/swex2022/data/model/IOP10/';
f1='swex2022_IOP10_LES_run_1km_545_24_z81_2way_d05_ndown_111m_81_try4_km_opt4_dt30.mat';
load(strcat(folder,f1));
LES111m=LES;
clear('LES')



[obs_U_lvl1,obs_V_lvl1]=wspdWdirToUV(obs_wspd_lvl1,obs_wdir_lvl1);
[isfs_U_10m,isfs_V_10m]=wspdWdirToUV(isfs_wspd_10m,isfs_wdir_10m);
[iss_U_10m,iss_V_10m]=wspdWdirToUV(iss_wspd_10m,iss_wdir_10m);

[isfs_Tdew_2m]=fdewpoint(isfs_temp_2m,isfs_RH_2m);


%% select the stations we want to plot. Make LES reference as these domains
% have the least spatial domains

stSelect=[1:21]; % the plots for these "stations" are shown


% find the stations in model data arrays
clear('ho','ha','hu','he','hi','idx')
iRaws=zeros(length(stSelect),1);
iNws=zeros(length(stSelect),1);
iIsfs=zeros(length(stSelect),1);
iIss=zeros(length(stSelect),1);
iMeso=zeros(length(stSelect),1);
iLong=zeros(length(stSelect),1);
for st=stSelect
    ho=strcmp(stations,LES111m.obs_stations{st}); %  check for RAWS or NWS stations existence
    ha=strcmp(isfs_st,LES111m.obs_stations{st}); % check for ISFS station existence
    hu=strcmp(iss_st,LES111m.obs_stations{st}); % check for ISFS station existence
    he=strcmp(WRFsens.obs_stations,LES111m.obs_stations{st}); % check WRF for "station" existence
    hi=strcmp(WRFsens.obs_stations,LES111m.obs_stations{st}); % check WRF for "station" existence
    if sum(ho)~=0 % if station is found in stations
        if strfind(stations{ho},'C1') % if C1 is in string, it's RAWS
            iRaws(st,:)=find(ho==1);
            iNws(st,:)=0; % zeros are dummies
            iIsfs(st,:)=0;
            iIss(st,:)=0;
        else % if no C1, it's likely ASOS
            iRaws(st,:)=0; % zeros are dummies
            iNws(st,:)=find(ho==1);
            iIsfs(st,:)=0;
            iIss(st,:)=0;
        end
    elseif sum(ha)~=0 % this is ISFS
        iRaws(st,:)=0; % zeros are dummies
        iNws(st,:)=0;
        iIsfs(st,:)=find(ha==1);
        iIss(st,:)=0;
    elseif sum(hu)~=0 % this is ISS
        iRaws(st,:)=0; % zeros are dummies
        iNws(st,:)=0;
        iIsfs(st,:)=0;
        iIss(st,:)=find(hu==1);
    end
    
    clear('ho','ha','hu')

    if sum(he)~=0 % this is mesoscale WRF sensitivity
        iMeso(st,:)=find(he==1);
    end
    clear('he')
    
    if sum(hi)~=0 % this is WRF long run
        iLong(st,:)=find(hi==1);
    else
        iLong(st,:)=0;
    end
    clear('hi')
end



%% plot the stufffffff

% some definitions
maxAlt=6000;

IOPs=['IOP01';'IOP02';'IOP03';'IOP04';'IOP05';'IOP06';'IOP07';'IOP08';...
      'IOP09';'IOP10';...
      'EOP01';'EOP02';'EOP03'];
IOPdates=['04-04';'04-05';'04-13';'04-18';'04-23';'04-28';'05-07';'05-08';... % reflects the start day of the IOP
    '05-10';'05-12';...
    '04-17';'04-25';'05-04'];
xlimsIOPs=[datenum(2022,4,4,10,0,0);datenum(2022,4,5,10,0,0);
           datenum(2022,4,13,10,0,0);datenum(2022,4,18,10,0,0);
           datenum(2022,4,23,10,0,0);datenum(2022,4,28,10,0,0);
           datenum(2022,5,7,10,0,0);datenum(2022,5,8,10,0,0);
           datenum(2022,5,10,10,0,0);datenum(2022,5,12,10,0,0);
           datenum(2022,4,17,10,0,0);datenum(2022,4,25,10,0,0);
           datenum(2022,5,4,10,0,0);];

clear('LES111m','LES333m')
% load the LES runs in a structure
folder='/home/sbarc/duine/sundowners/swex2022/data/model/IOP10/';
% runs={'swex2022_IOP10_LES_run_1km_540_25_z107_1way_d05_ndown_111m_107_try4_km_opt5_dt30';...
%       'swex2022_IOP10_LES_run_1km_545_24_z107_2way_d05_ndown_111m_107_try4_km_opt2_dt30';...
%       'swex2022_IOP10_LES_run_1km_545_24_z107_2way_d05_ndown_111m_107_try4_km_opt4_dt30';...
%       'swex2022_IOP10_LES_run_1km_545_24_z81_2way_d05_ndown_111m_81_try4_km_opt4_dt30';...
%       'swex2022_IOP10_LES_run_1km_545_24_z55_2way_d05_ndown_111m_55_try4_km_opt4_dt25'};

% runsTi={'SMS3DTKE z107';'MYNN+3D TKE 1.5';'MYNN+1D Sm. z107';'MYNN+1D Sm. z81';'MYNN+1D Sm. z55'};

runs={'swex2022_IOP10_LES_run_1km_545_24_z107_2way_d04_ndown_111m_107_2Dom_km_opt2_dt30';...
      'swex2022_IOP10_LES_run_1km_545_24_z107_2way_d04_ndown_111m_107_2Dom_smth1km_km_opt2_dt30'};%...
    % 'swex2022_IOP10_LES_run_1km_545_24_z107_2way_d05_ndown_111m_107_try4_km_opt4_dt30'};

runsTi={'111m';'111m-smth1km'};%;'MYNN+1D Sm. z107'};

% load the runs, but also shorten the filename to enable structure naming
for run=1:size(runs,1)
    runX=runs{run}(20:end);
    load(strcat(folder,runs{run}));
    LES111m.(runX)=LES;
    LES111m.(runX).Tdew2_st=fdewpoint(LES111m.(runX).T2_st,LES111m.(runX).RH2);

    clear('LES')
end
LES111m.(runX).obs_stations{19}='MFS-2';
LES111m.(runX).obs_stations{21}='SBBG';
% cmap=[1 0 0;1 0 0;1 0 0;1 0 0;1 0 0;1 0 0;1 0 0;1 0 0;1 0 0;1 0 0;1 0 0;1 0 0;1 0 0;];
cmap=distinguishable_colors(size(runs,1),'w');
fg=[0.13 0.54 0.13];
% c333m=[0.54 0.13 0.13];
% c111m=[0.54 0.54 0.13];
LineStyles={'-.';'--';':';'-.';'--';':';'-.'};
Markers={'x','o','^','+','>','<','v'};


for i=10%:length(IOPs)

    IOPx=IOPs(i,:);
    disp(['working on ',IOPx])
    figFolder=[cd,'/'];
    xTicks=[xlimsIOPs(i):datenum(0,0,0,4,0,0):xlimsIOPs(i)+1];
    % % load forecast. There are three forecasts cycles that cover a single
    % % IOP
    % folder='/home/sbarc/duine/sundowners/swex2022/data/model/';
    % f1=['WRF_fcst_SWEX2022_vertProfs_moreLocs_',IOPx,'.mat'];
    % load(strcat(folder,f1));
    % WRFfcst=WRF;
    % clear('WRF')

    % load simulation. Based on ERA5.
    % folder='/home/sbarc/duine/sundowners/swex2022/data/model/';
    % f1=['WRF_ERA5_SWEX2022_vertProfs_moreLocs_',IOPx,'.mat'];
    % load(strcat(folder,f1));
    % WRFlong=WRF;
    % clear('WRF')

    % load 1km simulation data
    folder='/home/sbarc/duine/sundowners/swex2022/data/model/IOP10/';
    f1=['WRF_SWEX2022_IOP10_sensitivity_vertProfs_moreLocs.mat'];
    load(strcat(folder,f1));
    WRF1km=WRF.run_1km_545_24_z107_2way;
    WRF1km.hgt_st=WRF.hgt_st;
    WRF1km.obs_lats=WRF.obs_lats;
    WRF1km.obs_lons=WRF.obs_lons;
    WRF1km.obs_stations=WRF.obs_stations;
    WRF1km.zFroude_st=WRF.zFroude_st;
    clear('WRF')

    WRF1km.Tdew2_st=fdewpoint(WRF1km.T2_st,WRF1km.RH2);

    % load file with sensitivity
    folder='/home/sbarc/duine/sundowners/swex2022/data/model/IOP10/';
    f1='WRF_SWEX2022_IOP10_sensitivity_vertProfs_moreLocs.mat';
    load(strcat(folder,f1));
    WRFsens=WRF;
    clear('WRF')

    %%%% actual plotting from here %%%%%
    for st=[19,21]%2,3,10,11,15,19,21]%stSelect



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % froude nr and V-component average only %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % froude number
        close
        fig2=figGD('off');
        sp1=subplot_tight(2,2,1,[0.1 0.15]);
        hold on; grid on; box on;
        ylimFr=6;
        % Define shading limits
        shade_y = [1, ylimFr]; % X-range for shading
        shade_x = [xlimsIOPs(i) xlimsIOPs(i)+1];    % Use current Y-axis limits for shading

        % Plot shading before other elements
        fill([shade_x(1) shade_x(2) shade_x(2) shade_x(1)], ...
            [shade_y(1) shade_y(1) shade_y(2) shade_y(2)], ...
            [0.8 0.8 0.8], ... % Gray color
            'EdgeColor', 'none', 'FaceAlpha', 0.5); % Transparent shading
        hold on;
        p1=plot(WRF1km.tnumPDT,WRF1km.fri(iLong(st),:),...
            'color',fg,'Linewidth',3,'LineStyle','-');
        % part for WRF sensitivity
        c=1;
        fn=fieldnames(LES111m);
        for ic=1:length(fn)
            p2(c)=plot(LES111m.(fn{ic}).tnumPDT,LES111m.(fn{ic}).fri(st,:),...
                'color',cmap(c,:),'Linewidth',2,'LineStyle',LineStyles(c));
            legStr{c}=[runsTi{c}];
            c=c+1;
        end
        clear('c')
        if st==19
            p3=line(tnumPDT_launch,fri,'LineStyle','none','Marker','x','color','k','LineWidth',10,'MarkerSize',15);
        end
        
        if st==19 % legend only for MFS-2
            leg=legend([p1;p2(:);p3],{'WRF 1km',legStr{:},'Radiosounding'},'Location','NorthWest');
        % else
        %     leg=legend([p1;p2(:)],{'WRF 1km',legStr{:}},'Location','NorthWest');
        end
        ylabel('Fr [-]')
        ylim([0 ylimFr])
        if st==19
            % text(datenum(xlimsIOPs(i)+0.5),4.5,'(a)','FontWeight','bold')
            title(['(a) ' LES111m.(fn{ic}).obs_stations{st}])
        elseif st ==21
            title(['(c) ' LES111m.(fn{ic}).obs_stations{st}])% text(datenum(xlimsIOPs(i)+0.5),4.5,'(b)','FontWeight','bold')
        end
        % xlabel(['Days/hours in ',datestr(WRF1km.tnumPDT(25),'mmm'),' 2022 [PDT]'])
        xlim([xlimsIOPs(i) xlimsIOPs(i)+1])
        datetick('x','dd/HH','keeplimits')
        XlabelsTickzzz=char('','14:00','18:00','22:00','02:00','06:00','10:00');
        set(gca,'XTick',xTicks,'XTickLabels','')        

        set(findall(gcf, '-property', 'FontSize'), 'FontSize',16)
        clear('legObs','p0','p1','p2')



        % V10
        sp3=subplot_tight(2,2,3,[0.1 0.15]);
        hold on; grid on; box on;
        p1=plot(WRF1km.tnumPDT,WRF1km.meanVzi(iLong(st),:),...
            'color',fg,'Linewidth',3,'LineStyle','-');

         % part for WRF sensitivity
        c=1;
        fn=fieldnames(LES111m);
        substring='run';
        fnCheck=contains(fn,substring);
        for ic=1:length(fn)
            if fnCheck(ic)==1
                p2(c)=plot(LES111m.(fn{ic}).tnumPDT,LES111m.(fn{ic}).meanVzi(st,:),...
                    'color',cmap(c,:),'Linewidth',2,'LineStyle',LineStyles(c));
                legStr{c}=[runsTi{c}];
                c=c+1;
            end
        end
        clear('c')
        if st==19
            p3=line(tnumPDT_launch,meanVzi,'LineStyle','none','Marker','x','color','k','LineWidth',10,'MarkerSize',15);
        end
        % title([LES111m.(fn{ic}).obs_stations{st}])
        % leg=legend([p1;p2(:)],{'1km',legStr{:}},'Location','');        
        ylabel('Mean V (zi) [m s^{-1}]')
        xlabel(['Days/hours in ',datestr(WRF1km.tnumPDT(25),'mmm'),' 2022 [PDT]'])
        % if st==19 % only a legend for MFS-2
        %     leg=legend([p1;p2(:);p3],{'WRF 1km',legStr{:},'Radiosounding'},'Location','SouthWest');
        % % else
        % %     leg=legend([p1;p2(:)],{'WRF 1km',legStr{:}},'Location','SouthWest');
        % end
        if st==19
            % text(datenum(xlimsIOPs(i)+0.5),4.5,'(a)','FontWeight','bold')
            title(['(b) ' LES111m.(fn{ic}).obs_stations{st}])
        elseif st ==21
            title(['(d) ' LES111m.(fn{ic}).obs_stations{st}])% text(datenum(xlimsIOPs(i)+0.5),4.5,'(b)','FontWeight','bold')
        end

        % if st==19
        %     text(datenum(xlimsIOPs(i)+0.5),2.5,'(c)','FontWeight','bold')
        % elseif st ==21
        %     text(datenum(xlimsIOPs(i)+0.5),2.5,'(d)','FontWeight','bold')
        % end

        ylim([-20 5])
        xlim([xlimsIOPs(i) xlimsIOPs(i)+1])
        datetick('x','dd/HH','keeplimits')
        set(gca,'XTick',xTicks,'XTickLabels',XlabelsTickzzz)        

        set(findall(gcf, '-property', 'FontSize'), 'FontSize',16)
        % set(sp1,'')
        set(sp3,'Position',[0.1500 0.1700 0.2750 0.3500])
        figname=strcat([figFolder,'fig10_meteogram_gprime_vmean_',IOPx,'_',LES111m.(fn{ic}).obs_stations{st}]);
        export_fig(figname,'-png')
        print(gcf, '-dpdf',figname)
        % export_fig(figname,'-pdf')
        clear('figname','legObs','p0','p1','p2','leg')







    end
end




stop








        % 10-m winds
        close
        fig2=figGD('off');

        subplot_tight(2,2,1,[0.1 0.05])
        obsCheck=0; % make it zero for starters, to test for legend purposes
        hold on; grid on; box on;

        % need to address RAWS, NWS, ISFS. All have different measurement
        % strategies or come from different data source, which are reflected below
        if iRaws(st)>0 % RAWS: 6.1 m
            obsCheck=1;
            legObs=stations{iRaws(st)};
            p0=plot(obs_tnum_PDT_lvl1(iRaws(st),:),obs_wspd_lvl1(iRaws(st),:),...
                'color','k','Linewidth',3,'LineStyle','-');
            p1=plot(WRF1km.tnumPDT,WRF1km.wspd6_st(iLong(st),:),...
                'color',fg,'Linewidth',3,'LineStyle','-');

        elseif iIsfs(st)>0 % ISFS: 10 m
            obsCheck=1;
            legObs=isfs_st{iIsfs(st)};
            p0=plot(isfs_tnumPDT(:),isfs_wspd_10m(iIsfs(st),:),...
                'color','k','Linewidth',3,'LineStyle','-');
            p1=plot(WRF1km.tnumPDT,WRF1km.wspd10_st(iLong(st),:),...
                'color',fg,'Linewidth',3,'LineStyle','-');

        elseif iIss(st)>0 % ISS: 10 m
            obsCheck=1;
            legObs=iss_st{iIss(st)};
            p0=plot(iss_tnumPDT(iIss(st),:),iss_wspd_10m(iIss(st),:),...
                'color','k','Linewidth',3,'LineStyle','-');
            p1=plot(WRF1km.tnumPDT,WRF1km.wspd10_st(iLong(st),:),...
                'color',fg,'Linewidth',3,'LineStyle','-');

        elseif iNws(st)>0 % NWS shold be 10 m.
            obsCheck=1;
            legObs=stations{iNws(st)};
            p0=plot(obs_tnum_PDT_lvl1(iNws(st),:),obs_wspd_lvl1(iNws(st),:),...
                'color','k','Linewidth',3,'LineStyle','-');
            p1=plot(WRF1km.tnumPDT,WRF1km.wspd10_st(iLong(st),:),...
                'color',fg,'Linewidth',3,'LineStyle','-');

        else % if there is no RAWS, no ISFS, no NWS, just model output
            p1=plot(WRF1km.tnumPDT,WRF1km.wspd10_st(iLong(st),:),...
                'color',fg,'Linewidth',3,'LineStyle','-');
        end

        % part for WRF sensitivity
        c=1;
        fn=fieldnames(LES111m);
        for ic=1:length(fn)
            if iRaws(st)~=0 % RAWS is 6.1 m
                p2(c)=plot(LES111m.(fn{ic}).tnumPDT,LES111m.(fn{ic}).wspd6_st(st,:),...
                    'color',cmap(c,:),'Linewidth',2,'LineStyle',LineStyles(c));
            else % ISFS, NWS and others
                p2(c)=plot(LES111m.(fn{ic}).tnumPDT,LES111m.(fn{ic}).wspd10_st(st,:),...
                    'color',cmap(c,:),'Linewidth',2,'LineStyle',LineStyles(c));
            end

            legStr{c}=[runsTi{c}];
            c=c+1;
        end
        clear('c')

        title([LES111m.(fn{ic}).obs_stations{st}])
        if obsCheck==1;   leg=legend([p0;p1;p2(:)],{['Obs ',legObs],'1km',legStr{:}},'Location','WestOutside');
        else;       leg=legend([p1;p2(:)],{'1km',legStr{:}},'Location','WestOutside');
        end
        ylabel('Near-surface wind speed [m s^{-1}]')
        xlabel(['Days/hours in ',datestr(WRF1km.tnumPDT(25),'mmm'),' 2022 [PDT]'])

        xlim([xlimsIOPs(i) xlimsIOPs(i)+1])
        datetick('x','dd/HH','keeplimits')
        clear('legObs','p0','p1','p2')


        subplot_tight(2,2,2,[0.1 0.05])
        obsCheck=0; % make it zero for starters, to test for legend purposes
        hold on; grid on; box on;

        % need to address RAWS, NWS, ISFS. All have different measurement
        % strategies or come from different data source, which are reflected below
        if iRaws(st)>0 % RAWS: 6.1 m
            obsCheck=1;
            legObs=stations{iRaws(st)};
            p0=plot(obs_tnum_PDT_lvl1(iRaws(st),:),obs_wdir_lvl1(iRaws(st),:),...
                'color','k','Linewidth',3,'LineStyle','none','Marker','+');
            p1=plot(WRF1km.tnumPDT,WRF1km.wdir6_st(iLong(st),:),...
                'color',fg,'Linewidth',3,'LineStyle','none','Marker','*');

        elseif iIsfs(st)>0 % ISFS: 10 m
            obsCheck=1;
            legObs=isfs_st{iIsfs(st)};
            p0=plot(isfs_tnumPDT(:),isfs_wdir_10m(iIsfs(st),:),...
                'color','k','Linewidth',3,'LineStyle','none','Marker','+');
            p1=plot(WRF1km.tnumPDT,WRF1km.wdir10_st(iLong(st),:),...
                'color',fg,'Linewidth',3,'LineStyle','none','Marker','*');
        elseif iIss(st)>0 % ISS: 10 m
            obsCheck=1;
            legObs=iss_st{iIss(st)};
            p0=plot(iss_tnumPDT(iIss(st),:),iss_wdir_10m(iIss(st),:),...
                'color','k','Linewidth',3,'LineStyle','none','Marker','+');
            p1=plot(WRF1km.tnumPDT,WRF1km.wdir10_st(iLong(st),:),...
                'color',fg,'Linewidth',3,'LineStyle','none','Marker','*');
        elseif iNws(st)>0 % NWS shold be 10 m.
            obsCheck=1;
            legObs=stations{iNws(st)};
            p0=plot(obs_tnum_PDT_lvl1(iNws(st),:),obs_wdir_lvl1(iNws(st),:),...
                'color','k','Linewidth',3,'LineStyle','none','Marker','+');
            p1=plot(WRF1km.tnumPDT,WRF1km.wdir10_st(iLong(st),:),...
                'color',fg,'Linewidth',3,'LineStyle','none','Marker','*');
        else % if there is no RAWS, no ISFS, no NWS, just model output
            p1=plot(WRF1km.tnumPDT,WRF1km.wdir10_st(iLong(st),:),...
                'color',fg,'Linewidth',3,'LineStyle','none','Marker','*');
        end


        % part for WRF sensitivity
        c=1;
        fn=fieldnames(LES111m);
        for ic=1:length(fn)
            if iRaws(st)~=0 % RAWS is 6.1 m
                p2(c)=plot(LES111m.(fn{ic}).tnumPDT,LES111m.(fn{ic}).wdir6_st(st,:),...
                    'color',cmap(c,:),'Linewidth',2,'LineStyle','none');
            else % ISFS, NWS and others
                p2(c)=plot(LES111m.(fn{ic}).tnumPDT,LES111m.(fn{ic}).wdir10_st(st,:),...
                    'color',cmap(c,:),'Linewidth',2,'LineStyle','none');
            end
            legStr{c}=[runsTi{c}];
            c=c+1;
        end
        clear('c')
        title([LES111m.(fn{ic}).obs_stations{st}])
        if obsCheck==1;   leg=legend([p0;p1;p2(:)],{['Obs ',legObs],'1km',legStr{:}},'Location','EastOutside');
        else;       leg=legend([p1;p2(:)],{'1km',legStr{:}},'Location','EastOutside');
        end
        ylabel('Wind direction [^{\circ}N]')
        set(gca,'YLim',[0 360],'YTick',[0:90:360],'YTickLabel',{'N';'E';'S';'W';'N'})
        xlabel(['Days/hours in ',datestr(WRF1km.tnumPDT(25),'mmm'),' 2022 [PDT]'])
        xlim([xlimsIOPs(i) xlimsIOPs(i)+1])
        datetick('x','dd/HH','keeplimits')
        clear('legObs','p0','p1','p2')


        % temperature here
        subplot_tight(2,2,3,[0.1 0.05])
        obsCheck=0; % make it zero for starters, to test for legend purposes
        hold on; grid on; box on;

        % need to address RAWS, NWS, ISFS. All have different measurement
        % strategies or come from different data source, which are reflected below
        % these were designed for wind that are observed at different
        % levels, however, for temperature this is probably different
        if iRaws(st)>0 % RAWS: 6.1 m
            obsCheck=1;
            legObs=stations{iRaws(st)};
            p0=plot(obs_tnum_PDT_lvl1(iRaws(st),:),obs_temp_lvl1(iRaws(st),:),...
                'color','k','Linewidth',3,'LineStyle','-');
            p1=plot(WRF1km.tnumPDT,WRF1km.T2_st(iLong(st),:)-273.15,...
                'color',fg,'Linewidth',3,'LineStyle','-');
        elseif iIsfs(st)>0 % ISFS: 2 m
            obsCheck=1;
            legObs=isfs_st{iIsfs(st)};
            p0=plot(isfs_tnumPDT(:),isfs_temp_2m(iIsfs(st),:),...
                'color','k','Linewidth',3,'LineStyle','-');
            p1=plot(WRF1km.tnumPDT,WRF1km.T2_st(iLong(st),:)-273.15,...
                'color',fg,'Linewidth',3,'LineStyle','-');
        elseif iIss(st)>0 % ISS: 10 m
            obsCheck=1;
            legObs=iss_st{iIss(st)};
            p0=plot(iss_tnumPDT(iIss(st),:),iss_Temp_2m(iIss(st),:),... % note that ISS is v01!!!
                'color','k','Linewidth',3,'LineStyle','-');
            p1=plot(WRF1km.tnumPDT,WRF1km.T2_st(iLong(st),:)-273.15,...
                'color',fg,'Linewidth',3,'LineStyle','-');

        elseif iNws(st)>0 % NWS 
            obsCheck=1;
            legObs=stations{iNws(st)};
            p0=plot(obs_tnum_PDT_lvl1(iNws(st),:),obs_temp_lvl1(iNws(st),:),...
                'color','k','Linewidth',3,'LineStyle','-');
            p1=plot(WRF1km.tnumPDT,WRF1km.T2_st(iLong(st),:)-273.15,...
                'color',fg,'Linewidth',3,'LineStyle','-');
        else % if there is no RAWS, no ISFS, no NWS, just model output
            p1=plot(WRF1km.tnumPDT,WRF1km.T2_st(iLong(st),:)-273.15,...
                'color',fg,'Linewidth',3,'LineStyle','-');
        end



        % part for WRF sensitivity
        c=1;
        fn=fieldnames(LES111m);
        for ic=1:length(fn)
            % if iRaws(st)~=0 % RAWS is 6.1 m
                p2(c)=plot(LES111m.(fn{ic}).tnumPDT,LES111m.(fn{ic}).T2_st(st,:)-273.15,...
                    'color',cmap(c,:),'Linewidth',2,'LineStyle',LineStyles(c));
            % else % ISFS, NWS and others
            %     p2(c)=plot(LES111m.(fn{ic}).tnumPDT,LES111m.(fn{ic}).T2_st(st,:)-273.15,...
            %         'color',cmap(c,:),'Linewidth',2,'LineStyle',LineStyles(c));
            % end
            legStr{c}=[runsTi{c}];
            c=c+1;
        end
        clear('c')
        title([LES111m.(fn{ic}).obs_stations{st}])
        if obsCheck==1;   leg=legend([p0;p1;p2(:)],{['Obs ',legObs],'1km',legStr{:}},'Location','WestOutside');
        else;       leg=legend([p1;p2(:)],{'1km',legStr{:}},'Location','WestOutside');
        end
        ylabel('Temperature 2-m [^{\circ}C]')
        xlabel(['Days/hours in ',datestr(WRF1km.tnumPDT(25),'mmm'),' 2022 [PDT]'])

        xlim([xlimsIOPs(i) xlimsIOPs(i)+1])
        datetick('x','dd/HH','keeplimits')
        clear('legObs','p0','p1','p2')




        % dew point temperature here
        subplot_tight(2,2,4,[0.1 0.05])
        obsCheck=0; % make it zero for starters, to test for legend purposes
        hold on; grid on; box on;

        % need to address RAWS, NWS, ISFS. All have different measurement
        % strategies or come from different data source, which are reflected below
        % these were designed for wind that are observed at different
        % levels, however, for temperature this is probably different
        if iRaws(st)>0 % RAWS: 6.1 m
            obsCheck=1;
            legObs=stations{iRaws(st)};
            p0=plot(obs_tnum_PDT_lvl1(iRaws(st),:),obs_temp_dew_lvl1(iRaws(st),:),...
                'color','k','Linewidth',3,'LineStyle','-');
            p1=plot(WRF1km.tnumPDT,WRF1km.Tdew2_st(iLong(st),:),...
                'color',fg,'Linewidth',3,'LineStyle','-');
        elseif iIsfs(st)>0 % ISFS: 2 m
            obsCheck=1;
            legObs=isfs_st{iIsfs(st)};
            p0=plot(isfs_tnumPDT(:),isfs_Tdew_2m(iIsfs(st),:),...
                'color','k','Linewidth',3,'LineStyle','-');
            p1=plot(WRF1km.tnumPDT,WRF1km.Tdew2_st(iLong(st),:),...
                'color',fg,'Linewidth',3,'LineStyle','-');
        elseif iIss(st)>0 % ISS: 10 m
            obsCheck=1;
            legObs=iss_st{iIss(st)};
            p0=plot(iss_tnumPDT(iIss(st),:),iss_Tdew_2m(iIss(st),:),... % note that ISS is v01!!!
                'color','k','Linewidth',3,'LineStyle','-');
            p1=plot(WRF1km.tnumPDT,WRF1km.Tdew2_st(iLong(st),:),...
                'color',fg,'Linewidth',3,'LineStyle','-');
          

        elseif iNws(st)>0 % NWS 
            obsCheck=1;
            legObs=stations{iNws(st)};
            p0=plot(obs_tnum_PDT_lvl1(iNws(st),:),obs_temp_dew_lvl1(iNws(st),:),...
                'color','k','Linewidth',3,'LineStyle','-');
            p1=plot(WRF1km.tnumPDT,WRF1km.Tdew2_st(iLong(st),:),...
                'color',fg,'Linewidth',3,'LineStyle','-');
        else % if there is no RAWS, no ISFS, no NWS, just model output
            p1=plot(WRF1km.tnumPDT,WRF1km.Tdew2_st(iLong(st),:),...
                'color',fg,'Linewidth',3,'LineStyle','-');
        end

          % part for WRF sensitivity
        c=1;
        fn=fieldnames(LES111m);
        for ic=1:length(fn)
            % if iRaws(st)~=0 % RAWS is 6.1 m
                p2(c)=plot(LES111m.(fn{ic}).tnumPDT,LES111m.(fn{ic}).Tdew2_st(st,:),...
                    'color',cmap(c,:),'Linewidth',2,'LineStyle',LineStyles(c));
            % else % ISFS, NWS and others
            %     p2(c)=plot(WRFsens.(fn{ic}).tnumPDT,WRFsens.(fn{ic}).Tdew2_st(iMeso(st),:),...
            %         'color',cmap(c,:),'Linewidth',2,'LineStyle',LineStyles(c));
            % end
            legStr{c}=[runsTi{c}];
            c=c+1;
        end
        clear('c')
        title([LES111m.(fn{ic}).obs_stations{st}])
        if obsCheck==1;   leg=legend([p0;p1;p2(:)],{['Obs ',legObs],'1km',legStr{:}},'Location','EastOutside');
        else;       leg=legend([p1;p2(:)],{'1km',legStr{:}},'Location','EastOutside');
        end
        ylabel('Dew point 2-m [^{\circ}C]')
        xlabel(['Days/hours in ',datestr(WRF1km.tnumPDT(25),'mmm'),' 2022 [PDT]'])

        xlim([xlimsIOPs(i) xlimsIOPs(i)+1])
        datetick('x','dd/HH','keeplimits')
        clear('legObs','p0','p1','p2')



        set(findall(gcf, '-property', 'FontSize'), 'FontSize',16)

        figname=strcat([figFolder,'/meteogram_v10_wspdSfc_wdir_temps_',IOPx,'_',LES111m.(fn{ic}).obs_stations{st}]);
        export_fig(figname,'-png')
        %     export_fig(figname,'-pdf')
        clear('p0','figname')
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 10-m V-component
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        close
        fig2=figGD('off');
        subplot_tight(2,2,1,[0.1 0.05])
        hold on; grid on; box on;

        obsCheck=0; % make it zero for starters, to test for legend purposes

        % need to address RAWS, NWS, ISFS. All have different measurement
        % strategies or come from different data source, which are reflected below
        if iRaws(st)>0 % RAWS: 6.1 m
            obsCheck=1;
            legObs=stations{iRaws(st)};
            p0=plot(obs_tnum_PDT_lvl1(iRaws(st),:),obs_V_lvl1(iRaws(st),:),...
                'color','k','Linewidth',3,'LineStyle','-');
            p1=plot(WRF1km.tnumPDT,WRF1km.V6_st(iLong(st),:),...
                'color',fg,'Linewidth',3,'LineStyle','-');

        elseif iIsfs(st)>0 % ISFS: 10 m
            obsCheck=1;
            legObs=isfs_st{iIsfs(st)};
            p0=plot(isfs_tnumPDT(:),isfs_V_10m(iIsfs(st),:),...
                'color','k','Linewidth',3,'LineStyle','-');
            p1=plot(WRF1km.tnumPDT,WRF1km.V10_st(iLong(st),:),...
                'color',fg,'Linewidth',3,'LineStyle','-');

        elseif iIss(st)>0 % ISS: 10 m
            obsCheck=1;
            legObs=iss_st{iIss(st)};
            p0=plot(iss_tnumPDT(iIss(st),:),iss_V_10m(iIss(st),:),...
                'color','k','Linewidth',3,'LineStyle','-');
            p1=plot(WRF1km.tnumPDT,WRF1km.V10_st(iLong(st),:),...
                'color',fg,'Linewidth',3,'LineStyle','-');

        elseif iNws(st)>0 % NWS shold be 10 m.
            obsCheck=1;
            legObs=stations{iNws(st)};
            p0=plot(obs_tnum_PDT_lvl1(iNws(st),:),obs_V_lvl1(iNws(st),:),...
                'color','k','Linewidth',3,'LineStyle','-');
            p1=plot(WRF1km.tnumPDT,WRF1km.V10_st(iLong(st),:),...
                'color',fg,'Linewidth',3,'LineStyle','-');

        else % if there is no RAWS, no ISFS, no NWS, just model output
            p1=plot(WRF1km.tnumPDT,WRF1km.V10_st(iLong(st),:),...
                'color',fg,'Linewidth',3,'LineStyle','-');

        end

        % part for WRF sensitivity
        c=1;
        fn=fieldnames(LES111m);
        substring='run';
        fnCheck=contains(fn,substring);
        for ic=1:length(fn)
            if iRaws(st)~=0 % RAWS is 6.1 m
                p2(c)=plot(LES111m.(fn{ic}).tnumPDT,LES111m.(fn{ic}).V6_st(st,:),...
                    'color',cmap(c,:),'Linewidth',2,'LineStyle',LineStyles(c));
            else % ISFS, NWS and others
                p2(c)=plot(LES111m.(fn{ic}).tnumPDT,LES111m.(fn{ic}).V10_st(st,:),...
                    'color',cmap(c,:),'Linewidth',2,'LineStyle',LineStyles(c));
            end
            legStr{c}=[runsTi{c}];
            c=c+1;

        end
        clear('c')

        title([LES111m.(fn{ic}).obs_stations{st}])
        if obsCheck==1;   leg=legend([p0;p1;p2(:)],{['Obs ',legObs],'1km',legStr{:}},'Location','EastOutside');
        else;       leg=legend([p1;p2(:)],{'1km',legStr{:}},'Location','WestOutside');
        end
        ylabel('Near-surface V-component [m s^{-1}]')
        xlabel(['Days/hours in ',datestr(WRF1km.tnumPDT(25),'mmm'),' 2022 [PDT]'])

        xlim([xlimsIOPs(i) xlimsIOPs(i)+1])
        datetick('x','dd/HH','keeplimits')

        set(findall(gcf, '-property', 'FontSize'), 'FontSize',16)

        figname=strcat([figFolder,'meteogram_v10_Vsfc_',IOPx,'_',LES111m.(fn{ic}).obs_stations{st}]);
        export_fig(figname,'-png')
        %     export_fig(figname,'-pdf')
        clear('p0','figname') % emove the one legend
        clear('legObs','p0','p1','p2')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % V-component 1000 m msl
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        close
        fig2=figGD('off');
        subplot_tight(1,3,1,[0.25 0.05])
        hold on; grid on; box on;

        p1=plot(WRF1km.tnumPDT,WRF1km.V1000(iLong(st),:),...
            'color',fg,'Linewidth',3,'LineStyle','-');

        % part for WRF sensitivity
        c=1;
        fn=fieldnames(LES111m);
        for ic=1:length(fn)
            p2(c)=plot(LES111m.(fn{ic}).tnumPDT,LES111m.(fn{ic}).V1000(st,:),...
                'color',cmap(c,:),'Linewidth',2,'LineStyle',LineStyles(c));
            legStr{c}=[runsTi{c}];
            c=c+1;
        end
        clear('c')

        title([LES111m.(fn{ic}).obs_stations{st}])
        leg=legend([p1;p2(:)],{'1km',legStr{:}},'Location','SouthOutside');        
        ylabel('V-component 1000 m asl [m s^{-1}]')
        xlabel(['Days/hours in ',datestr(WRF1km.tnumPDT(25),'mmm'),' 2022 [PDT]'])
        xlim([xlimsIOPs(i) xlimsIOPs(i)+1])
        datetick('x','dd/HH','keeplimits')
        
        set(findall(gcf, '-property', 'FontSize'), 'FontSize',16)
        clear('legObs','p0','p1','p2')

        


        % V-component 1500 m msl
        subplot_tight(1,3,2,[0.25 0.05])
        hold on; grid on; box on;

        p1=plot(WRF1km.tnumPDT,WRF1km.V1500(iLong(st),:),...
            'color',fg,'Linewidth',3,'LineStyle','-');
        % part for WRF sensitivity
        c=1;
        fn=fieldnames(LES111m);
        for ic=1:length(fn)
            p2(c)=plot(LES111m.(fn{ic}).tnumPDT,LES111m.(fn{ic}).V1500(st,:),...
                'color',cmap(c,:),'Linewidth',2,'LineStyle',LineStyles(c));
            legStr{c}=[runsTi{c}];
            c=c+1;
        end
        clear('c')

        title([LES111m.(fn{ic}).obs_stations{st}])
        leg=legend([p1;p2(:)],{'1km',legStr{:}},'Location','SouthOutside');        
        ylabel('V-component 1500 m asl [m s^{-1}]')
        xlabel(['Days/hours in ',datestr(WRF1km.tnumPDT(25),'mmm'),' 2022 [PDT]'])
        xlim([xlimsIOPs(i) xlimsIOPs(i)+1])
        datetick('x','dd/HH','keeplimits')

        set(findall(gcf, '-property', 'FontSize'), 'FontSize',16)
        clear('legObs','p0','p1','p2')


        % V-component ridge level
        subplot_tight(1,3,3,[0.25 0.05])
        hold on; grid on; box on;
        p1=plot(WRF1km.tnumPDT,WRF1km.VridgeLvl(iLong(st),:),...
            'color',fg,'Linewidth',3,'LineStyle','-');

         % part for WRF sensitivity
        c=1;
        fn=fieldnames(LES111m);
        for ic=1:length(fn)
            if fnCheck(ic)==1
                p2(c)=plot(LES111m.(fn{ic}).tnumPDT,LES111m.(fn{ic}).VridgeLvl(st,:),...
                    'color',cmap(c,:),'Linewidth',2,'LineStyle',LineStyles(c));
                legStr{c}=[runsTi{c}];
                c=c+1;
            end
        end
        clear('c')

        title([LES111m.(fn{ic}).obs_stations{st}])
        leg=legend([p1;p2(:)],{'1km',legStr{:}},'Location','SouthOutside');        
        ylabel('V-component ridge level [m s^{-1}]')
        xlabel(['Days/hours in ',datestr(WRF1km.tnumPDT(25),'mmm'),' 2022 [PDT]'])

        xlim([xlimsIOPs(i) xlimsIOPs(i)+1])
        datetick('x','dd/HH','keeplimits')

        set(findall(gcf, '-property', 'FontSize'), 'FontSize',16)

        figname=strcat([figFolder,'meteogram_v10_V1000_V1500_VridgeLvl_',IOPx,'_',LES111m.(fn{ic}).obs_stations{st}]);
        export_fig(figname,'-png')
        %     export_fig(figname,'-pdf')
        clear('figname','legObs','p0','p1','p2','leg')











        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PBLH, fr, V10
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % froude number
        close
        fig2=figGD('off');
        subplot_tight(2,2,1,[0.1 0.15])
        hold on; grid on; box on;
        ylimFr=5;
        % Define shading limits
        shade_y = [1, ylimFr]; % X-range for shading
        shade_x = [xlimsIOPs(i) xlimsIOPs(i)+1];    % Use current Y-axis limits for shading

        % Plot shading before other elements
        fill([shade_x(1) shade_x(2) shade_x(2) shade_x(1)], ...
            [shade_y(1) shade_y(1) shade_y(2) shade_y(2)], ...
            [0.8 0.8 0.8], ... % Gray color
            'EdgeColor', 'none', 'FaceAlpha', 0.5); % Transparent shading
        hold on;
        p1=plot(WRF1km.tnumPDT,WRF1km.fri(iLong(st),:),...
            'color',fg,'Linewidth',3,'LineStyle','-');
        % part for WRF sensitivity
        c=1;
        fn=fieldnames(LES111m);
        for ic=1:length(fn)
            p2(c)=plot(LES111m.(fn{ic}).tnumPDT,LES111m.(fn{ic}).fri(st,:),...
                'color',cmap(c,:),'Linewidth',2,'LineStyle',LineStyles(c));
            legStr{c}=[runsTi{c}];
            c=c+1;
        end
        clear('c')

        title([LES111m.(fn{ic}).obs_stations{st}])
        leg=legend([p1;p2(:)],{'1km',legStr{:}},'Location','NorthWest');
        
        ylabel('Fr [-]')
        ylim([0 ylimFr])
        xlabel(['Days/hours in ',datestr(WRF1km.tnumPDT(25),'mmm'),' 2022 [PDT]'])
        xlim([xlimsIOPs(i) xlimsIOPs(i)+1])
        datetick('x','dd/HH','keeplimits')
        
        set(findall(gcf, '-property', 'FontSize'), 'FontSize',16)
        clear('legObs','p0','p1','p2')

        


        % PBLHeight
        subplot_tight(2,2,2,[0.1 0.05])
        hold on; grid on; box on;

        p1=plot(WRF1km.tnumPDT,WRF1km.PBLH_st(iLong(st),:),...
            'color',fg,'Linewidth',3,'LineStyle','-');
        % part for WRF sensitivity
        c=1;
        fn=fieldnames(LES111m);
        for ic=1:length(fn)
            p2(c)=plot(LES111m.(fn{ic}).tnumPDT,LES111m.(fn{ic}).PBLH_st(st,:),...
                'color',cmap(c,:),'Linewidth',2,'LineStyle',LineStyles(c));
            legStr{c}=[runsTi{c}];
            c=c+1;
        end
        clear('c')

        title([LES111m.(fn{ic}).obs_stations{st}])
        leg=legend([p1;p2(:)],{'1km',legStr{:}},'Location','EastOutside');        
        ylabel('PBL height [m agl]')
        xlabel(['Days/hours in ',datestr(WRF1km.tnumPDT(25),'mmm'),' 2022 [PDT]'])
        xlim([xlimsIOPs(i) xlimsIOPs(i)+1])
        datetick('x','dd/HH','keeplimits')

        set(findall(gcf, '-property', 'FontSize'), 'FontSize',16)
        clear('legObs','p0','p1','p2')


        % V10
        subplot_tight(2,2,3,[0.1 0.15])
        hold on; grid on; box on;
        p1=plot(WRF1km.tnumPDT,WRF1km.meanVzi(iLong(st),:),...
            'color',fg,'Linewidth',3,'LineStyle','-');

         % part for WRF sensitivity
        c=1;
        fn=fieldnames(LES111m);
        for ic=1:length(fn)
            if fnCheck(ic)==1
                p2(c)=plot(LES111m.(fn{ic}).tnumPDT,LES111m.(fn{ic}).meanVzi(st,:),...
                    'color',cmap(c,:),'Linewidth',2,'LineStyle',LineStyles(c));
                legStr{c}=[runsTi{c}];
                c=c+1;
            end
        end
        clear('c')

        title([LES111m.(fn{ic}).obs_stations{st}])
        % leg=legend([p1;p2(:)],{'1km',legStr{:}},'Location','');        
        ylabel('Mean V (zi) [m s^{-1}]')
        xlabel(['Days/hours in ',datestr(WRF1km.tnumPDT(25),'mmm'),' 2022 [PDT]'])

        xlim([xlimsIOPs(i) xlimsIOPs(i)+1])
        datetick('x','dd/HH','keeplimits')


        % PBLHeight
        subplot_tight(2,2,4,[0.1 0.05])
        hold on; grid on; box on;

        p1=plot(WRF1km.tnumPDT,WRF1km.PBLH_st(iLong(st),:),...
            'color',fg,'Linewidth',3,'LineStyle','-');
        % part for WRF sensitivity
        c=1;
        fn=fieldnames(LES111m);
        for ic=1:length(fn)
            p2(c)=plot(LES111m.(fn{ic}).tnumPDT,LES111m.(fn{ic}).g_prime(st,:),...
                'color',cmap(c,:),'Linewidth',2,'LineStyle',LineStyles(c));
            legStr{c}=[runsTi{c}];
            c=c+1;
        end
        clear('c')

        title([LES111m.(fn{ic}).obs_stations{st}])
        leg=legend([p1;p2(:)],{'1km',legStr{:}},'Location','EastOutside');        
        ylabel('g`[m s^{-2}]')
        ylim([0 0.5])

        xlabel(['Days/hours in ',datestr(WRF1km.tnumPDT(25),'mmm'),' 2022 [PDT]'])
        xlim([xlimsIOPs(i) xlimsIOPs(i)+1])
        datetick('x','dd/HH','keeplimits')

        set(findall(gcf, '-property', 'FontSize'), 'FontSize',16)
        clear('legObs','p0','p1','p2')



        set(findall(gcf, '-property', 'FontSize'), 'FontSize',16)

        figname=strcat([figFolder,'meteogram_PBLH_Fr_gprime_v10_',IOPx,'_',LES111m.(fn{ic}).obs_stations{st}]);
        export_fig(figname,'-png')
        %     export_fig(figname,'-pdf')
        clear('figname','legObs','p0','p1','p2','leg')







        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % jet strength (V-component) and height
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        close
        fig2=figGD('off');

        % jet strength
        subplot_tight(2,2,1,[0.1 0.05])
        hold on; grid on; box on;

        % ERA5 longterm
        [minV,idZ]=min(WRF1km.V(iLong(st),:,:),[],2); % compute minimum in the jet and its height
        minV=squeeze(minV);
        zMinV=diag(squeeze(WRF1km.z(iLong(st),squeeze(idZ),:)));
        p1=plot(WRF1km.tnumPDT,minV,...
            'color',fg,'Linewidth',3,'LineStyle','-');




        % part for WRF sensitivity
        c=1;
        fn=fieldnames(LES111m);
        for ic=1:length(fn)
            % WRF sensitivity
            [minV,idZ]=min(LES111m.(fn{ic}).V(st,:,:),[],2); % compute minimum in the jet and its height
            minV=squeeze(minV);
            zMinV=diag(squeeze(LES111m.(fn{ic}).z(st,squeeze(idZ),:)));
            p2(c)=plot(LES111m.(fn{ic}).tnumPDT,minV,...
                'color',cmap(c,:),'Linewidth',2,'LineStyle',LineStyles(c));
            legStr{c}=[runsTi{c}];
            c=c+1;
        end
        clear('c')

        title([LES111m.(fn{ic}).obs_stations{st}])
        leg=legend([p1;p2(:)],{'1km',legStr{:}},'Location','WestOutside');
        ylabel('Min V-component of jet [m s^{-1}]')
        xlabel(['Days/hours in ',datestr(WRF1km.tnumPDT(25),'mmm'),' 2022 [PDT]'])

        xlim([xlimsIOPs(i) xlimsIOPs(i)+1])
        datetick('x','dd/HH','keeplimits')
        clear('leg','legObs','p0','p1','p2')

        % jet height
        subplot_tight(2,2,2,[0.1 0.05])
        hold on; grid on; box on;

        % first ERA5
        [minV,idZ]=min(WRF1km.V(iLong(st),:,:),[],2); % compute minimum in the jet and its height
        minV=squeeze(minV);
        zMinV=diag(squeeze(WRF1km.z(iLong(st),squeeze(idZ),:)));
        p1=plot(WRF1km.tnumPDT,zMinV,...
            'color',fg,'Linewidth',3,'LineStyle','-');



        % part for WRF sensitivity
        c=1;
        fn=fieldnames(LES111m);
        for ic=1:length(fn)

            % WRF sensitivity
            [minV,idZ]=min(LES111m.(fn{ic}).V(st,:,:),[],2); % compute minimum in the jet and its height
            minV=squeeze(minV);
            zMinV=diag(squeeze(LES111m.(fn{ic}).z(st,squeeze(idZ),:)));
            p2(c)=plot(LES111m.(fn{ic}).tnumPDT,zMinV,...
                'color',cmap(c,:),'Linewidth',2,'LineStyle',LineStyles(c));
            legStr{c}=[runsTi{c}];
            c=c+1;

        end
        clear('c')

        title([LES111m.(fn{ic}).obs_stations{st}])
        xlabel(['Days/hours in ',datestr(WRF1km.tnumPDT(25),'mmm'),' 2022 [PDT]'])
        ylabel('Elevation of peak in jet [m agl]')

        line([WRF1km.tnumPDT(1) WRF1km.tnumPDT(end)],[LES111m.(fn{ic}).zFroude_st(st) LES111m.(fn{ic}).zFroude_st(st)],...
            'LineStyle','--','color',fg,'LineWidth',1.5)
        xlim([xlimsIOPs(i) xlimsIOPs(i)+1])
        ylim([0 4000])
        datetick('x','dd/HH','keeplimits')
        leg=legend([p1;p2(:)],{'1km',legStr{:}},'Location','EastOutside');




        set(findall(gcf, '-property', 'FontSize'), 'FontSize',16)

        figname=strcat([figFolder,'meteogram_v10_minV_elevationJet_',IOPx,'_',LES111m.(fn{ic}).obs_stations{st}]);
        export_fig(figname,'-png')
        %     export_fig(figname,'-pdf')
        clear('figname','leg','legObs','p0','p1','p2')
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % jet strength (total wind speed) and height
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        close
        fig2=figGD('off');

        % jet strength
        subplot_tight(2,2,1,[0.1 0.05])
        obsCheck=0; % make it zero for starters, to test for legend purposes

        hold on; grid on; box on;

        if iRaws(st)>0 % plot RAWS gusts
            obsCheck=1;
            legObs=['Obs gust ',stations{iRaws(st)},' 6.1m'];
            p0=plot(obs_tnum_PDT_lvl1(iRaws(st),:),obs_wspd_gust_lvl1(iRaws(st),:),...
                'color','k','Linewidth',3,'LineStyle','-');
        elseif iNws(st)>0 % for NWS it might exist
            obsCheck=1;
            legObs=['Obs gust ',stations{iNws(st)},' 10m'];
            p0=plot(obs_tnum_PDT_lvl1(iNws(st),:),obs_wspd_gust_lvl1(iNws(st),:),...
                'color','k','Linewidth',3,'LineStyle','-');
        elseif iIss(st)>0 % for ISS there is max wind available, but it is at 3 m instead of 10m
            obsCheck=1;
            legObs=['Obs gust ',iss_st{iIss(st)},' 3m'];
            p0=plot(iss_tnumPDT(iIss(st),:),iss_wspdMax_3m(iIss(st),:),...
                'color','k','Linewidth',3,'LineStyle','-');
        end

        % ERA5
        [maxWspd,idZ]=max(WRF1km.wspd(iLong(st),1:10,:),[],2); % compute minimum in the jet and its height. Only consider lowest 4 km (z=1:16)
        maxWspd=squeeze(maxWspd);
        zMaxWspd=diag(squeeze(WRF1km.z(iLong(st),squeeze(idZ),:)));
        p1=plot(WRF1km.tnumPDT,maxWspd,...
            'color',fg,'Linewidth',3,'LineStyle','-');

 
         % part for WRF sensitivity
        c=1;
        fn=fieldnames(LES111m);
        substring='run';
        fnCheck=contains(fn,substring);
        for ic=1:length(fn)
            % WRF sensitivity
            [maxWspd,idZ]=max(LES111m.(fn{ic}).wspd(st,1:10,:),[],2); % compute minimum in the jet and its height. Only consider lowest 4 km (z=1:16)
            maxWspd=squeeze(maxWspd);
            zMaxWspd=diag(squeeze(LES111m.(fn{ic}).z(st,squeeze(idZ),:)));
            p2(c)=plot(LES111m.(fn{ic}).tnumPDT,maxWspd,...
                'color',cmap(c,:),'Linewidth',2,'LineStyle',LineStyles(c));
            legStr{c}=[runsTi{c}];
            c=c+1;
        end
        clear('c')

        title([LES111m.(fn{ic}).obs_stations{st}])
        if obsCheck==1;   leg=legend([p0;p1;p2(:)],{['Obs ',legObs],'1km',legStr{:}},'Location','WestOutside');
        else;       leg=legend([p1;p2(:)],{'1km',legStr{:}},'Location','WestOutside');
        end
        ylabel('Max wind of jet [m s^{-1}]')
        xlabel(['Days/hours in ',datestr(WRF1km.tnumPDT(25),'mmm'),' 2022 [PDT]'])
        xlim([xlimsIOPs(i) xlimsIOPs(i)+1])
        datetick('x','dd/HH','keeplimits')
        clear('leg','legObs','p0','p1','p2')
        

        % jet height
        subplot_tight(2,2,2,[0.1 0.05])
        hold on; grid on; box on;

        % first ERA5
        [maxWspd,idZ]=max(WRF1km.wspd(iLong(st),1:10,:),[],2); % compute maximum in the jet and its height in lowest 4 km
        maxWspd=squeeze(maxWspd);
        zMaxWspd=diag(squeeze(WRF1km.z(iLong(st),squeeze(idZ),:)));
        p1=plot(WRF1km.tnumPDT,zMaxWspd,...
            'color',fg,'Linewidth',3,'LineStyle','-');


        % part for WRF sensitivity
        c=1;
        fn=fieldnames(LES111m);
        for ic=1:length(fn)
            % WRF sensitivity
            [maxWspd,idZ]=max(LES111m.(fn{ic}).wspd(st,1:10,:),[],2); % compute minimum in the jet and its height in lowest 4 km
            maxWspd=squeeze(maxWspd);
            zMaxWspd=diag(squeeze(LES111m.(fn{ic}).z(st,squeeze(idZ),:)));
            p2(c)=plot(LES111m.(fn{ic}).tnumPDT,zMaxWspd,...
                'color',cmap(c,:),'Linewidth',2,'LineStyle',LineStyles(c));
            legStr{c}=[runsTi{c}];
            c=c+1;
        end
        clear('c')

        title([LES111m.(fn{ic}).obs_stations{st}])
        xlabel(['Days/hours in ',datestr(WRF1km.tnumPDT(25),'mmm'),' 2022 [PDT]'])
        ylabel('Elevation of peak in jet [m agl]')

        line([WRF1km.tnumPDT(1) WRF1km.tnumPDT(end)],[WRF1km.zFroude_st(st) WRF1km.zFroude_st(st)],...
            'LineStyle','--','color',fg,'LineWidth',1.5)
        xlim([xlimsIOPs(i) xlimsIOPs(i)+1])
        ylim([0 4000])
        datetick('x','dd/HH','keeplimits')
        leg=legend([p1;p2(:)],{'1km',legStr{:}},'Location','EastOutside');



        set(findall(gcf, '-property', 'FontSize'), 'FontSize',16)

        figname=strcat([figFolder,'/meteogram_v10_maxWspd_elevationJet_',IOPx,'_',LES111m.(fn{ic}).obs_stations{st}]);
        export_fig(figname,'-png')
        
        %     export_fig(figname,'-pdf')
        clear('figname','leg','legObs','p0','p1','p2')
%     end % end of station loop
%     clear('WRF','WRF1km','WRFsens','IOPx','figFolder')
% end % end of IOP loop
toc