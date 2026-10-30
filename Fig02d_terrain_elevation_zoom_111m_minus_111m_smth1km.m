% script to plot terrain of several domains

clear
clc
tic

% location 
pathfile='/home/voyager-sbarc/wrf/wrf422/WPS-4.2.2-duine/geo_em/SWEX2022/IOP10/111m/2Doms_smth1km/';
runpath1(1,:) = strcat(pathfile,'/geo_em.d02.nc');

pathfile='/home/voyager-sbarc/wrf/wrf422/WPS-4.2.2-duine/geo_em/SWEX2022/IOP10/111m/2Doms/';
runpath2(1,:) = strcat(pathfile,'/geo_em.d02.nc');

%% Import surface and other data data

for i = 1
    doX=strcat('d0',num2str(i));
    WRF1.(doX).Lons = double(ncread(runpath1(i,:),'XLONG_M'));
    WRF1.(doX).Lats = double(ncread(runpath1(i,:),'XLAT_M'));
    WRF1.(doX).HGT = double(ncread(runpath1(i,:),'HGT_M'));
    WRF1.(doX).LU = double(ncread(runpath1(i,:),'LU_INDEX'));
    
    % some constants of the domains
    WRF1.(doX).DX=ncreadatt(runpath1(i,:),'/','DX')/1000;
    WRF1.(doX).DY=ncreadatt(runpath1(i,:),'/','DY')/1000;
    WRF1.(doX).i_parentStart=ncreadatt(runpath1(i,:),'/','i_parent_start');
    WRF1.(doX).j_parentStart=ncreadatt(runpath1(i,:),'/','j_parent_start');
    WRF1.(doX).i_parentEnd=ncreadatt(runpath1(i,:),'/','i_parent_end');
    WRF1.(doX).j_parentEnd=ncreadatt(runpath1(i,:),'/','j_parent_end');
    WRF1.(doX).e_we=ncreadatt(runpath1(i,:),'/','WEST-EAST_GRID_DIMENSION');
    WRF1.(doX).e_sn=ncreadatt(runpath1(i,:),'/','SOUTH-NORTH_GRID_DIMENSION');
    
    
end


for i = 1
    doX=strcat('d0',num2str(i));
    WRF2.(doX).Lons = double(ncread(runpath2(i,:),'XLONG_M'));
    WRF2.(doX).Lats = double(ncread(runpath2(i,:),'XLAT_M'));
    WRF2.(doX).HGT = double(ncread(runpath2(i,:),'HGT_M'));
    WRF2.(doX).LU = double(ncread(runpath2(i,:),'LU_INDEX'));
    
    % some constants of the domains
    WRF2.(doX).DX=ncreadatt(runpath2(i,:),'/','DX')/1000;
    WRF2.(doX).DY=ncreadatt(runpath2(i,:),'/','DY')/1000;
    WRF2.(doX).i_parentStart=ncreadatt(runpath2(i,:),'/','i_parent_start');
    WRF2.(doX).j_parentStart=ncreadatt(runpath2(i,:),'/','j_parent_start');
    WRF2.(doX).i_parentEnd=ncreadatt(runpath2(i,:),'/','i_parent_end');
    WRF2.(doX).j_parentEnd=ncreadatt(runpath2(i,:),'/','j_parent_end');
    WRF2.(doX).e_we=ncreadatt(runpath2(i,:),'/','WEST-EAST_GRID_DIMENSION');
    WRF2.(doX).e_sn=ncreadatt(runpath2(i,:),'/','SOUTH-NORTH_GRID_DIMENSION');
    
    
end

% stLES  ={     'HP', 'Refugio',  'OXLFC',     'KIZA',    'GVTC1',       's9',         's3', 'Sedgwick'};
% latsLES=[  34.5187;  34.517;  34.4898;   34.6061;     34.4884;    34.4758;        34.4657;    34.6955];
% lonsLES=[-120.0789;-120.075;-120.0459; -120.0756;     -120.2357;  -120.2146;    -120.0275;  -120.0470];
% 
% stLES  ={'Refugio',  'OXLFC',     'KIZA',    'GVTC1',       's9',         's3', 'Sedgwick'};
% latsLES=[   34.517;  34.4898;   34.6061;     34.4884;    34.4758;        34.4657;    34.6955];
% lonsLES=[ -120.075;-120.0459; -120.0756;     -120.2357;  -120.2146;    -120.0275;  -120.0470];

load('/home/voyager-sbarc/duine/sundowners/swex2022/instruments/SWEX2022_instrument_locations_BAMS_v2.mat')


%% plotting

res={strcat(num2str(WRF1.d01.DX),' km');'';''};

% cmapsea=[1 1 1];
% HGT_d01(HGT_d01<1)=20;
% HGT_d01(LU_INDEX_d01==17)=-100;
xLims=[-119.75 -119.6];
yLims=[34.4 34.55];
cLims=[-250 250];
close
fig=figGD('on');
for d=1
    cmapsea=[1 1 1];

%     if (d==1 || d==2 || d==3 )
%         WRF.(doX).HGT(WRF.(doX).LU==17)=-20;
%     elseif d==4
%         WRF.(doX).HGT(WRF.(doX).HGT<1)=20;
%         WRF.(doX).HGT(WRF.(doX).LU==17)=-20;
%     end
    
    
    sp(d)=subplot_tight(1,1,1,[0.2 0.34]);
    
    doX=strcat('d0',num2str(d));
    bwr=bluewhitered(40);
    colormap(bwr);
    % demcmap(WRF.(doX).HGT,100,cmapsea);
    p1(d)=pcolor(squeeze(WRF1.(doX).Lons),squeeze(WRF1.(doX).Lats),squeeze(WRF2.(doX).HGT)-squeeze(WRF1.(doX).HGT));
    hold on;
    
    
%     p2=contour(WRF.(doX).Lons,WRF.(doX).Lats,WRF.(doX).HGT,[0:200:1200],'color','k','LineWidth',2);
    
    shading flat
    cb = colorbar;
    ylabel(cb,'Elevation above sea level (m)')
    
    if (d==1 || d==2 || d==3 )
        bordersGD('continental us','k','LineWidth',3)
%         bordersGD('countries','k','LineWidth',2)
    elseif d==4
        bordersGD('continental us','k','LineWidth',3)
    end
    box on
    


    for st=1:length(stLIDAR)
        m1=line(lonsLIDAR(st),latsLIDAR(st),...
            'LineStyle','None','Marker','^','color','k','LineWidth',3,'MarkerSize',12);
    end

    m1a=line(-119.6593,34.4402,'LineStyle','None','Marker','^','color','k','LineWidth',3,'MarkerSize',12);

    % for st=1:length(stISFS)
    %     m2=line(lonsISFS(st),latsISFS(st),...
    %         'LineStyle','None','Marker','o','color','k','LineWidth',3,'MarkerSize',12);
    % end
    % orange=[1 0.5 0];
    % for st=1:length(stRAWS)
    %     m3=line(lonsRAWS(st),latsRAWS(st),...
    %         'LineStyle','None','Marker','x','color','k','LineWidth',3,'MarkerSize',20);
    % end
    % ASOS, same as RAWS
    % for st=1:length(stASOS)
    %     m4=line(lonsASOS(st),latsASOS(st),...
    %         'LineStyle','None','Marker','x','color','k','LineWidth',3,'MarkerSize',20);
    % end
    % 
    % lightblue=[0 1 1];%[0.67 0.84 0.90];
    % for st=1:length(stCeilo)
    %     m5=line(lonsCeilo(st),latsCeilo(st),...
    %         'LineStyle','None','Marker','+','color',lightblue,'LineWidth',3,'MarkerSize',12);
    % end

    % for st=1:length(stWindProf915)
    %     m6=line(lonsWindProf915(st),latsWindProf915(st),...
    %         'LineStyle','None','Marker','x','color','k','LineWidth',3,'MarkerSize',12);
    % end
    % 
    % for st=1:length(stWindProf449)
    %     m6=line(lonsWindProf449(st),latsWindProf449(st),...
    %         'LineStyle','None','Marker','x','color','k','LineWidth',3,'MarkerSize',12);
    % end

    

    
    
    % leg=legend([m1 m2 m3 m5 m6],...
    %     {'lidar','ISFS','RAWS/ASOS','ceilometer','449/915 MHz','Radiosonde'},'Location','East');
    
    dx=WRF1.(doX).DX*1000;
    doX=strcat('d0',num2str(d));
    % text(-119.8598,34.5050,strcat(res(d,:)),'FontSize',16,...
    %     'BackGroundColor','w','FontWeight','bold');
    title(['(d) 111m minus 111m smth1km']);

    set(gca,'CLim',cLims)
    caxis(cLims)
     % these come from LES:
    set(gca,'XLim',xLims,'YLim',yLims,'CLim',cLims,...
        'XTick',[-119.8 -119.7 -119.6],'YTick',[34.4 34.5 34.6],'tickdir','out')
    set(gca,'Position',[0.3200 0.0900 0.2590 0.6000])
%     set(cb,'Position',[0.6550 0.2501 0.0111 0.4997])
    leg=legend([m1],{'lidar'},'Location','SouthWest');


    % need to set it a second time
    bwr=bluewhitered(40);
    colormap(bwr);
    caxis([cLims])

    xlabel('Longitude [deg W]')
    ylabel('Latitude [deg N]')
    % set(gca,'XTick',[],'YTick',[])
    
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 18)
    
    
end
export_fig(['terrain_elevation_zoom_111m_minus_111m_smth1km'],'-png')
stop


%% land use map
for i=1:17
    [LU_name{i},LU_plot(i)]=lu_z0(i,1);
end

close
figGD('on')


for d=1:4
    sp(d)=subplot_tight(2,2,d,[0.05 0.1]);
    doX=strcat('d0',num2str(d));
    
    % use standard matlab colormap
%     colormap(hsv(17));
%     colormap(colorcube(17));
    colormap(cmapModisWRF);
%     colormap(vga);
    p1=pcolor(squeeze(WRF1.(doX).Lons),squeeze(WRF1.(doX).Lats),squeeze(WRF1.(doX).LU));
    shading flat; hold on
    if d==2
        cbar=colorbar;
        set(cbar,'Limits',[1 18],'Ticks',[1.5:17.5],'TickLabels',LU_name)
        set(gca,'Clim',[1 18]);
        ylabel(cbar,'Landuse index')
    end
    
    set(gca,'Clim',[1 18]);

    xlabel('Longitude [^\circ]'); ylabel('Latitude [^\circ]')
    bordersGD('continental us','k')
    hold on

    
%     [c1,h1]=contour(squeeze(WRF.(doX).Lons),squeeze(WRF.(doX).Lats),squeeze(WRF.(doX).HGT),[0:200:1200],'color','k','LineWidth',2,'ShowText','on');
%     clabel(c1,h1,'color','k')
    if d==1
        for j=1:4
            doX=strcat('d0',num2str(d+1));
            % box ljnes
            line([WRF1.(doX).Lons(1,1) WRF1.(doX).Lons(1,end)],[WRF1.(doX).Lats(1,1) WRF1.(doX).Lats(1,end)],'color','k','LineWidth',2)
            line([WRF1.(doX).Lons(1,end) WRF1.(doX).Lons(end,end)],[WRF1.(doX).Lats(1,end) WRF1.(doX).Lats(end,end)],'color','k','LineWidth',2)
            line([WRF1.(doX).Lons(end,end) WRF1.(doX).Lons(end,1)],[WRF1.(doX).Lats(end,end) WRF1.(doX).Lats(end,1)],'color','k','LineWidth',2)
            line([WRF1.(doX).Lons(end,1) WRF1.(doX).Lons(1,1)],[WRF1.(doX).Lats(end,1) WRF1.(doX).Lats(1,1)],'color','k','LineWidth',2)
        end
    elseif d==2
        for j=2:3
            doX=strcat('d0',num2str(d+1));
            % box ljnes
            line([WRF1.(doX).Lons(1,1) WRF1.(doX).Lons(1,end)],[WRF1.(doX).Lats(1,1) WRF1.(doX).Lats(1,end)],'color','k','LineWidth',2)
            line([WRF1.(doX).Lons(1,end) WRF1.(doX).Lons(end,end)],[WRF1.(doX).Lats(1,end) WRF1.(doX).Lats(end,end)],'color','k','LineWidth',2)
            line([WRF1.(doX).Lons(end,end) WRF1.(doX).Lons(end,1)],[WRF1.(doX).Lats(end,end) WRF1.(doX).Lats(end,1)],'color','k','LineWidth',2)
            line([WRF1.(doX).Lons(end,1) WRF1.(doX).Lons(1,1)],[WRF1.(doX).Lats(end,1) WRF1.(doX).Lats(1,1)],'color','k','LineWidth',2)
        end
    elseif d==3
        j=3;
        doX=strcat('d0',num2str(d+1));
        % box ljnes
        line([WRF1.(doX).Lons(1,1) WRF1.(doX).Lons(1,end)],[WRF1.(doX).Lats(1,1) WRF1.(doX).Lats(1,end)],'color','k','LineWidth',2)
        line([WRF1.(doX).Lons(1,end) WRF1.(doX).Lons(end,end)],[WRF1.(doX).Lats(1,end) WRF1.(doX).Lats(end,end)],'color','k','LineWidth',2)
        line([WRF1.(doX).Lons(end,end) WRF1.(doX).Lons(end,1)],[WRF1.(doX).Lats(end,end) WRF1.(doX).Lats(end,1)],'color','k','LineWidth',2)
        line([WRF1.(doX).Lons(end,1) WRF1.(doX).Lons(1,1)],[WRF1.(doX).Lats(end,1) WRF1.(doX).Lats(1,1)],'color','k','LineWidth',2)
    end
    % set(gca,'YLim',[34.3 34.8])
    % set(gca,'XLim',[-120.6 -119.2])
    if d==3 || d==4
        for i=1:length(latsLES)
            line(lonsLES(i),latsLES(i),'color','k','LineStyle','none','Marker','x','MarkerSize',10,'LineWidth',3)
        end
    end
    
    box on
    
    
end
export_fig('sbareg_d01-d04_LU_nlcd2011_to_modis','-png')
% export_fig('WRF_stations_LU_lvl1_zoom','-pdf')



toc