% Composite figure v1: VAD (left) + vertical stare (right) for HQ, SBBG, MFS-2
% Layout: subplot_tight(3,2,n) with explicit Position overrides
% Panels: (a) HQ VAD, (b) SBBG VAD, (c) MFS-2 VAD
%         (d) HQ stare, (e) SBBG stare, (f) MFS-2 stare
% GJ Duine, ERI UC Santa Barbara

clear; clc; tic; close all

%% ---- Settings ----
IOPstart  = datenum(2022,5,12,10,0,0);   % PDT
IOPend    = datenum(2022,5,13,10,0,0);
UTCoffset = 7/24;
xTicks    = IOPstart : datenum(0,0,0,4,0,0) : IOPend;
ylimVal   = [0 1500];
climVAD   = [0 20];
climW     = [-5 5];
tGrid     = IOPstart : 1/1440 : IOPend;

%% ---- 1. HQ VAD ----
dataVAD   = '/home/voyager-sbarc/SWEX2022_datasets/Lidars_ground/HQ/v10/VAD/';
hq_height = ncread([dataVAD 'VAD_42_20220512.nc'], 'height');

t12   = double(ncread([dataVAD 'VAD_42_20220512.nc'], 'time'));
ws12  = ncread([dataVAD 'VAD_42_20220512.nc'], 'wind_speed');
snr12 = ncread([dataVAD 'VAD_42_20220512.nc'], 'mean_snr');
ws12(snr12 < -28.5) = NaN;
tHQ12 = datenum(2022,5,12,0,0,0) + t12/86400 - UTCoffset;

t13   = double(ncread([dataVAD 'VAD_42_20220513.nc'], 'time'));
ws13  = ncread([dataVAD 'VAD_42_20220513.nc'], 'wind_speed');
snr13 = ncread([dataVAD 'VAD_42_20220513.nc'], 'mean_snr');
ws13(snr13 < -28.5) = NaN;
tHQ13 = datenum(2022,5,13,0,0,0) + t13/86400 - UTCoffset;

k12 = tHQ12 >= IOPstart & tHQ12 <= IOPend;
k13 = tHQ13 >= IOPstart & tHQ13 <= IOPend;
tHQ_vad     = [tHQ12(k12); tHQ13(k13)]';
wsHQ_vad    = [ws12(:,k12), ws13(:,k13)];
iAlt_HQ_vad = 1:51;

%% ---- 2. SBBG VAD ----
load('/home/voyager-sbarc/duine/sundowners/swex2022/data/sbbg/lidar/SBBG_VAD_116_IOP10.mat')
kSBBG      = tnumPDTSBBG >= IOPstart & tnumPDTSBBG <= IOPend;
tSBBG_vad  = tnumPDTSBBG(kSBBG);
wsSBBG_vad = wspdSBBG(:, kSBBG);
hSBBG_vad  = heightSBBG(:, kSBBG);
wsSBBG_vad(wsSBBG_vad>12 & hSBBG_vad>1000) =NaN; % to get rid of bade values abover 1000 m


%% ---- 3. MFS-2 VAD ----
load('/home/voyager-sbarc/duine/sundowners/swex2022/data/uwow/prelim/IOP10_stationary/UWOW_processed/UWOW_IOP10.mat')
wspdUwow(SNRUwow < 1.01) = NaN;
kUWOW      = tnumPDTUwow >= IOPstart & tnumPDTUwow <= IOPend;
tUWOW_vad  = tnumPDTUwow(kUWOW);
wsUWOW_vad = wspdUwow(:, kUWOW);
hUWOW_vad  = heightUwow(:, kUWOW);

%% ---- 4. HQ stare ----
dataVER12 = '/home/voyager-sbarc/SWEX2022_datasets/Lidars_ground/HQ/v10/VER/20220512/';
dataVER13 = '/home/voyager-sbarc/SWEX2022_datasets/Lidars_ground/HQ/v10/VER/20220513/';
files12   = dir([dataVER12 'cfrad.*.nc']);
files13   = dir([dataVER13 'cfrad.*.nc']);
allFiles  = [files12; files13];

fileTimesUTC = NaN(numel(allFiles), 1);
for k = 1:numel(allFiles)
    tok = regexp(allFiles(k).name, 'cfrad\.(\d{8})_(\d{6})_', 'tokens');
    if ~isempty(tok)
        d = tok{1}{1}; s = tok{1}{2};
        fileTimesUTC(k) = datenum(str2double(d(1:4)), str2double(d(5:6)), ...
            str2double(d(7:8)), str2double(s(1:2)), str2double(s(3:4)), str2double(s(5:6)));
    end
end
keep = fileTimesUTC >= (IOPstart+UTCoffset) & fileTimesUTC < (IOPend+UTCoffset);
allFiles = allFiles(keep);  fileTimesUTC = fileTimesUTC(keep);
[fileTimesUTC, ord] = sort(fileTimesUTC);
allFiles = allFiles(ord);

HQrange = ncread(fullfile(allFiles(1).folder, allFiles(1).name), 'range');
nGates  = numel(HQrange);
HQvertW = [];  HQtPDT = [];
for k = 1:numel(allFiles)
    fpath = fullfile(allFiles(k).folder, allFiles(k).name);
    tSec  = double(ncread(fpath, 'time'));
    W     = ncread(fpath, 'radial_wind_speed');
    cnr   = ncread(fpath, 'cnr');
    W(cnr < -22) = NaN;
    HQvertW = [HQvertW, W];
    HQtPDT  = [HQtPDT,  fileTimesUTC(k) - UTCoffset + tSec(:)'/86400];
end
iAlt_HQ_st = 1:29;

dt_hq     = diff(HQtPDT)*86400;
gapIdx_hq = find(dt_hq > 30);
stS_hq    = [1; gapIdx_hq(:)+1];
stE_hq    = [gapIdx_hq(:); numel(HQtPDT)];
WgridHQ   = NaN(nGates, numel(tGrid));
for k = 1:numel(stS_hq)
    cols = stS_hq(k):stE_hq(k);
    w_k  = mean(HQvertW(:,cols), 2, 'omitnan');
    inS  = tGrid >= HQtPDT(cols(1)) & tGrid <= HQtPDT(cols(end));
    WgridHQ(:, inS) = repmat(w_k, 1, sum(inS));
end

%% ---- 5. SBBG stare ----
load('/home/voyager-sbarc/duine/sundowners/swex2022/data/sbbg/lidar/SBBG_stare_IOP10.mat')
sbbgStareVertW(sbbgStareSNR < 1.005) = NaN;
iAlt_SBBG_st = 3:70;

dt_sb     = diff(sbbgStareTnumPDT)*86400;
gapIdx_sb = find(dt_sb > 30);
stS_sb    = [1; gapIdx_sb(:)+1];
stE_sb    = [gapIdx_sb(:); numel(sbbgStareTnumPDT)];
WgridSBBG = NaN(size(sbbgStareVertW,1), numel(tGrid));
for k = 1:numel(stS_sb)
    cols = stS_sb(k):stE_sb(k);
    w_k  = mean(sbbgStareVertW(:,cols), 2, 'omitnan');
    inS  = tGrid >= sbbgStareTnumPDT(cols(1)) & tGrid <= sbbgStareTnumPDT(cols(end));
    WgridSBBG(:, inS) = repmat(w_k, 1, sum(inS));
end

%% ---- 6. MFS-2 stare ----
load('/home/voyager-sbarc/duine/sundowners/swex2022/data/uwow/prelim/IOP10_stationary/stare/UWOW_stare_IOP10.mat')
uwowStareVertW(uwowStareSNR < 1.005) = NaN;
iAlt_UWOW_st = 3:70;

dt_uw     = diff(uwowStareTnumPDT)*86400;
gapIdx_uw = find(dt_uw > 30);
stS_uw    = [1; gapIdx_uw(:)+1];
stE_uw    = [gapIdx_uw(:); numel(uwowStareTnumPDT)];
WgridUWOW = NaN(size(uwowStareVertW,1), numel(tGrid));
for k = 1:numel(stS_uw)
    cols = stS_uw(k):stE_uw(k);
    w_k  = mean(uwowStareVertW(:,cols), 2, 'omitnan');
    inS  = tGrid >= uwowStareTnumPDT(cols(1)) & tGrid <= uwowStareTnumPDT(cols(end));
    WgridUWOW(:, inS) = repmat(w_k, 1, sum(inS));
end

%% ====================================================================
%% ---- Figure layout ----
% Margins and panel geometry (normalized figure units)
% Left col needs left margin for y-labels; centre gap holds VAD colorbar;
% right margin holds stare colorbar.  Panels are wide-short like individual plots.
lm = 0.08;   % left margin
cm = 0.06;   % centre gap  (VAD colorbar lives here)
rm = 0.08;   % right margin (stare colorbar lives here)
bm = 0.10;   % bottom margin
tm = 0.02;   % top margin
rg = 0.025;  % row gap

pw = (1 - lm - cm - rm) / 2;         % panel width  ≈ 0.345
ph = (1 - bm - tm - 2*rg) / 3;       % panel height ≈ 0.27
x1 = lm;                              % left column x
x2 = lm + pw + cm;                    % right column x
y3 = bm;                              % bottom row y
y2 = y3 + ph + rg;                    % middle row y
y1 = y2 + ph + rg;                    % top row y

fig = figGD('on');

% Create panels with subplot_tight, then override positions
gap = [0.05 0.09];
sp1 = subplot_tight(3,2,1,gap);  sp4 = subplot_tight(3,2,2,gap);
sp2 = subplot_tight(3,2,3,gap);  sp5 = subplot_tight(3,2,4,gap);
sp3 = subplot_tight(3,2,5,gap);  sp6 = subplot_tight(3,2,6,gap);

set(sp1,'Position',[x1 y1 pw ph]);  set(sp4,'Position',[x2 y1 pw ph]);
set(sp2,'Position',[x1 y2 pw ph]);  set(sp5,'Position',[x2 y2 pw ph]);
set(sp3,'Position',[x1 y3 pw ph]);  set(sp6,'Position',[x2 y3 pw ph]);

%% ---- Panel (a): HQ VAD ----
axes(sp1);
tM = repmat(tHQ_vad,                [numel(iAlt_HQ_vad) 1]);
hM = repmat(hq_height(iAlt_HQ_vad), [1 numel(tHQ_vad)]);
pcolor(tM, hM, wsHQ_vad(iAlt_HQ_vad,:));
shading interp;  caxis(climVAD);  colormap(sp1, jet(12));
set(sp1,'TickDir','out','Box','on','XTick',xTicks,'YLim',ylimVal,...
    'XLim',[IOPstart IOPend],'XTickLabels',{});
grid on;  ylabel('Elevation [m]');
text(IOPstart+0.015, 1380, '(a) SBCFD-HQ', 'FontWeight','bold','BackgroundColor','w');

%% ---- Panel (b): SBBG VAD ----
axes(sp2);
tM = repmat(tSBBG_vad, [size(wsSBBG_vad,1) 1]);
pcolor(tM, hSBBG_vad, wsSBBG_vad);
shading interp;  caxis(climVAD);  colormap(sp2, jet(12));
set(sp2,'TickDir','out','Box','on','XTick',xTicks,'YLim',ylimVal,...
    'XLim',[IOPstart IOPend],'XTickLabels',{});
grid on;  ylabel('Elevation [m]');
text(IOPstart+0.015, 1380, '(b) SBBG', 'FontWeight','bold','BackgroundColor','w');

%% ---- Panel (c): MFS-2 VAD ----
axes(sp3);
tM = repmat(tUWOW_vad, [size(wsUWOW_vad,1) 1]);
pcolor(tM, hUWOW_vad, wsUWOW_vad);
shading interp;  caxis(climVAD);  colormap(sp3, jet(12));
set(sp3,'TickDir','out','Box','on','XTick',xTicks,'YLim',ylimVal,...
    'XLim',[IOPstart IOPend],'XTickLabels',datestr(xTicks,'HH:MM'));
grid on;  ylabel('Elevation [m]');
xlabel(['Time ',datestr(IOPstart,'mmm dd'),' to ',datestr(IOPend,'mmm dd yyyy'),' [hours PDT]']);
text(IOPstart+0.015, 1380, '(c) MFS-2', 'FontWeight','bold','BackgroundColor','w');

%% ---- Panel (d): HQ stare ----
axes(sp4);
tM2 = repmat(tGrid,               [numel(iAlt_HQ_st) 1]);
hM2 = repmat(HQrange(iAlt_HQ_st), [1 numel(tGrid)]);
pcolor(tM2, hM2, WgridHQ(iAlt_HQ_st,:));
shading flat;  caxis(climW);
axes(sp4);  bwr = flipud(bluewhitered(30));  colormap(sp4, bwr);
set(sp4,'TickDir','out','Box','on','XTick',xTicks,'YLim',ylimVal,...
    'XLim',[IOPstart IOPend],'XTickLabels',{},'YTickLabels',{});
grid on;
text(IOPstart+0.015, 1380, '(d) SBCFD-HQ', 'FontWeight','bold','BackgroundColor','w');

%% ---- Panel (e): SBBG stare ----
axes(sp5);
tM2 = repmat(tGrid,                    [numel(iAlt_SBBG_st) 1]);
hM2 = repmat(sbbgElevation(iAlt_SBBG_st), [1 numel(tGrid)]);
pcolor(tM2, hM2, WgridSBBG(iAlt_SBBG_st,:));
shading flat;  caxis(climW);
axes(sp5);  bwr = flipud(bluewhitered(30));  colormap(sp5, bwr);
set(sp5,'TickDir','out','Box','on','XTick',xTicks,'YLim',ylimVal,...
    'XLim',[IOPstart IOPend],'XTickLabels',{},'YTickLabels',{});
grid on;
text(IOPstart+0.015, 1380, '(e) SBBG', 'FontWeight','bold','BackgroundColor','w');

%% ---- Panel (f): MFS-2 stare ----
axes(sp6);
tM2 = repmat(tGrid,                      [numel(iAlt_UWOW_st) 1]);
hM2 = repmat(uwowElevation(iAlt_UWOW_st), [1 numel(tGrid)]);
pcolor(tM2, hM2, WgridUWOW(iAlt_UWOW_st,:));
shading flat;  caxis(climW);
axes(sp6);  bwr = flipud(bluewhitered(30));  colormap(sp6, bwr);
set(sp6,'TickDir','out','Box','on','XTick',xTicks,'YLim',ylimVal,...
    'XLim',[IOPstart IOPend],'XTickLabels',datestr(xTicks,'HH:MM'),'YTickLabels',{});
grid on;
xlabel(['Time ',datestr(IOPstart,'mmm dd'),' to ',datestr(IOPend,'mmm dd yyyy'),' [hours PDT]']);
text(IOPstart+0.015, 1380, '(f) MFS-2', 'FontWeight','bold','BackgroundColor','w');

%% ---- Restore panel positions (pcolor/colormap can shift them) ----
set(sp1,'Position',[x1 y1 pw ph]);  set(sp4,'Position',[x2 y1 pw ph]);
set(sp2,'Position',[x1 y2 pw ph]);  set(sp5,'Position',[x2 y2 pw ph]);
set(sp3,'Position',[x1 y3 pw ph]);  set(sp6,'Position',[x2 y3 pw ph]);

%% ---- Colorbars ----
col_h  = y1+ph-y3;          % full column height
cb_h   = col_h / 2;         % colorbar = half column height, centred
cb_y   = y3 + col_h/4;      % bottom of colorbar
cb_w   = 0.009;              % narrow bar

% VAD: vertical bar placed in the centre gap
cb1 = colorbar(sp2, 'EastOutside');
set(sp2,'Position',[x1 y2 pw ph]);   % restore after colorbar shrinks it
cb1.Position = [x1+pw+0.008,  cb_y,  cb_w,  cb_h];
cb1.Ticks    = 0:5:20;
ylabel(cb1, 'Wind speed [m s^{-1}]');
colormap(sp2, jet(12));  caxis(sp2, climVAD);

% Stare: vertical bar placed in the right margin
axes(sp5);  bwr = flipud(bluewhitered(30));
cb2 = colorbar(sp5, 'EastOutside');
set(sp5,'Position',[x2 y2 pw ph]);
cb2.Position = [x2+pw+0.008,  cb_y,  cb_w,  cb_h];
cb2.Ticks    = climW(1):1:climW(2);
ylabel(cb2, 'W-component [m s^{-1}]');
colormap(sp5, bwr);  caxis(sp5, climW);

%% ---- Font and export ----
set(findall(fig,'-property','FontSize'),'FontSize',16);
export_fig('observations_lidar_VAD_stare_v1', '-png', '-r150');
export_fig('observations_lidar_VAD_stare_v1', '-pdf', '-r150');

toc
