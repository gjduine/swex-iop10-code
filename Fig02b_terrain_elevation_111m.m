% script to plot terrain of several domains

clear
clc
tic

% location 
pathfile='/home/voyager-sbarc/wrf/wrf422/WPS-4.2.2-duine/geo_em/SWEX2022/IOP10/111m/2Doms';

runpath(1,:) = strcat(pathfile,'/geo_em.d02.nc');

%% Import surface and other data data

for i = 1
    doX=strcat('d0',num2str(i));
    WRF.(doX).Lons = double(ncread(runpath(i,:),'XLONG_M'));
    WRF.(doX).Lats = double(ncread(runpath(i,:),'XLAT_M'));
    WRF.(doX).HGT = double(ncread(runpath(i,:),'HGT_M'));
    WRF.(doX).LU = double(ncread(runpath(i,:),'LU_INDEX'));
    
    % some constants of the domains
    WRF.(doX).DX=ncreadatt(runpath(i,:),'/','DX')/1000;
    WRF.(doX).DY=ncreadatt(runpath(i,:),'/','DY')/1000;
    WRF.(doX).i_parentStart=ncreadatt(runpath(i,:),'/','i_parent_start');
    WRF.(doX).j_parentStart=ncreadatt(runpath(i,:),'/','j_parent_start');
    WRF.(doX).i_parentEnd=ncreadatt(runpath(i,:),'/','i_parent_end');
    WRF.(doX).j_parentEnd=ncreadatt(runpath(i,:),'/','j_parent_end');
    WRF.(doX).e_we=ncreadatt(runpath(i,:),'/','WEST-EAST_GRID_DIMENSION');
    WRF.(doX).e_sn=ncreadatt(runpath(i,:),'/','SOUTH-NORTH_GRID_DIMENSION');
    
    
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

res={strcat(num2str(WRF.d01.DX),' km');'';''};
% cmapsea=[1 1 1];
% HGT_d01(HGT_d01<1)=20;
% HGT_d01(LU_INDEX_d01==17)=-100;

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
    % demcmap(WRF.(doX).HGT,100,cmapsea);
    cmap=cmapTerrElev(35:128);
    colormap(cmap);
    p1(d)=pcolor(squeeze(WRF.(doX).Lons),squeeze(WRF.(doX).Lats),squeeze(WRF.(doX).HGT));
    hold on;
    
    
%     p2=contour(WRF.(doX).Lons,WRF.(doX).Lats,WRF.(doX).HGT,[0:200:1200],'color','k','LineWidth',2);
    
    shading flat
    % cb = colorbar;
    % ylabel(cb,'Elevation above sea level (m)')
    
    if (d==1 || d==2 || d==3 )
        bordersGD('continental us','k','LineWidth',3)
%         bordersGD('countries','k','LineWidth',2)
    elseif d==4
        bordersGD('continental us','k','LineWidth',3)
    end
    box on
    


    for st=1:length(stLIDAR)
        m1=line(lonsLIDAR(st),latsLIDAR(st),...
            'LineStyle','None','Marker','^','color','r','LineWidth',3,'MarkerSize',12);
    end

    m1a=line(-119.6593,34.4402,'LineStyle','None','Marker','^','color','r','LineWidth',3,'MarkerSize',12);

    % for st=1:length(stISFS)
    %     m2=line(lonsISFS(st),latsISFS(st),...
    %         'LineStyle','None','Marker','o','color','k','LineWidth',3,'MarkerSize',10);
    % end
    % orange=[1 0.5 0];
    % for st=1:length(stRAWS)
    %     m3=line(lonsRAWS(st),latsRAWS(st),...
    %         'LineStyle','None','Marker','.','color',orange,'LineWidth',2,'MarkerSize',15);
    % end
    % ASOS, same as RAWS
    % for st=1:length(stASOS)
    %     m4=line(lonsASOS(st),latsASOS(st),...
    %         'LineStyle','None','Marker','.','color',orange,'LineWidth',2,'MarkerSize',15);
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
    
    
    doX=strcat('d0',num2str(d));
    % text(-119.8598,34.5050,strcat(res(d,:)),'FontSize',16,...
    %     'BackGroundColor','w','FontWeight','bold');
    title(strcat('(b) 111 m'))

    
     % these come from LES:
    set(gca,'XLim',[-119.8563 -119.5493],'YLim',[34.3495 34.62231],'CLim',[0 1400],...
        'XTick',[-119.8 -119.6],'YTick',[34.4 34.5 34.6],'YTickLabel','','tickdir','out')
    set(gca,'Position',[0.3200 0.0900 0.3090 0.6650])
%     set(cb,'Position',[0.6550 0.2501 0.0111 0.4997])


    xlabel('Longitude [deg W]')
    % ylabel('Latitude [deg N]')
    % set(gca,'XTick',[],'YTick',[])
    
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 18)
    
    
end
export_fig('terrain_elevation_111m','-png')
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
    p1=pcolor(squeeze(WRF.(doX).Lons),squeeze(WRF.(doX).Lats),squeeze(WRF.(doX).LU));
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
            line([WRF.(doX).Lons(1,1) WRF.(doX).Lons(1,end)],[WRF.(doX).Lats(1,1) WRF.(doX).Lats(1,end)],'color','k','LineWidth',2)
            line([WRF.(doX).Lons(1,end) WRF.(doX).Lons(end,end)],[WRF.(doX).Lats(1,end) WRF.(doX).Lats(end,end)],'color','k','LineWidth',2)
            line([WRF.(doX).Lons(end,end) WRF.(doX).Lons(end,1)],[WRF.(doX).Lats(end,end) WRF.(doX).Lats(end,1)],'color','k','LineWidth',2)
            line([WRF.(doX).Lons(end,1) WRF.(doX).Lons(1,1)],[WRF.(doX).Lats(end,1) WRF.(doX).Lats(1,1)],'color','k','LineWidth',2)
        end
    elseif d==2
        for j=2:3
            doX=strcat('d0',num2str(d+1));
            % box ljnes
            line([WRF.(doX).Lons(1,1) WRF.(doX).Lons(1,end)],[WRF.(doX).Lats(1,1) WRF.(doX).Lats(1,end)],'color','k','LineWidth',2)
            line([WRF.(doX).Lons(1,end) WRF.(doX).Lons(end,end)],[WRF.(doX).Lats(1,end) WRF.(doX).Lats(end,end)],'color','k','LineWidth',2)
            line([WRF.(doX).Lons(end,end) WRF.(doX).Lons(end,1)],[WRF.(doX).Lats(end,end) WRF.(doX).Lats(end,1)],'color','k','LineWidth',2)
            line([WRF.(doX).Lons(end,1) WRF.(doX).Lons(1,1)],[WRF.(doX).Lats(end,1) WRF.(doX).Lats(1,1)],'color','k','LineWidth',2)
        end
    elseif d==3
        j=3;
        doX=strcat('d0',num2str(d+1));
        % box ljnes
        line([WRF.(doX).Lons(1,1) WRF.(doX).Lons(1,end)],[WRF.(doX).Lats(1,1) WRF.(doX).Lats(1,end)],'color','k','LineWidth',2)
        line([WRF.(doX).Lons(1,end) WRF.(doX).Lons(end,end)],[WRF.(doX).Lats(1,end) WRF.(doX).Lats(end,end)],'color','k','LineWidth',2)
        line([WRF.(doX).Lons(end,end) WRF.(doX).Lons(end,1)],[WRF.(doX).Lats(end,end) WRF.(doX).Lats(end,1)],'color','k','LineWidth',2)
        line([WRF.(doX).Lons(end,1) WRF.(doX).Lons(1,1)],[WRF.(doX).Lats(end,1) WRF.(doX).Lats(1,1)],'color','k','LineWidth',2)
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