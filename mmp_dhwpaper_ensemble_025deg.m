%	             mmp_dhwpaper_ensemble_025deg
%
%  This script creates ensemble figures from CMIP5 and CMIP6
%  sea surface temperature data for the WWF Mesoamerican
%  Reef project and corresponding paper with GISS Model-E.
%  This version uses 0.25-degree resolution SST data.
%
%  -- Creates a map of the mean change in sea surface temperatures
%  for the month of September from the GISS-E2.1-G model for the
%  2050-2059 period under SSP585.
%
%  -- Creates a map of the unadjusted mean annual maximum DHW value
%  from the GISS-E2.1-G model for the 2050-2059 period under SSP585.
%
%  -- Creates a map of the bias-adjusted mean annual maximum DHW
%  value from the GISS-E2.1-G model for the 2050-2059 period under
%  SSP585.
%
%  -- Creates ensemble mean maps of the change in mean annual maximum
%  DHW for the 2050-2059 period under each SSP245 and SSP585, plus a
%  version for each the CMIP5 and CMIP6 ensembles only.
%
%  -- Creates a map of the ensemble mean annual maximum DHW for the
%  baseline (1980-2009) period across all unadjusted GCMs, plus a version
%  for each the CMIP5 and CMIP6 ensembles only.
%
%  -- Creates a map of the ensemble mean change in the annual frequency
%  of HotSpots greater than 1C for the 2050-2059 period under each SSP245
%  and SSP585, plus a version for each the CMIP5 and CMIP6 ensembles only.
%
%  -- Creates a map of the ensemble mean change in the annual maximum
%  consecutive duration of HotSpots greater than 1C for the 2050-2059
%  period under each SSP245 and SSP585, plus a version for each the CMIP5
%  and CMIP6 ensembles only.
%
%  INPUTS: none
%  OUTPUTS: none
%
%
%
%		                             author: Meridel Phillips
%                                  mmp2192@columbia.edu
%				                       date: 4/12/2021
%
function mmp_dhwpaper_ensemble_025deg();
%--------------------------------------------------
%--------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Begin Debug
%% End Debug

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /Users/mmphill2/Documents/MATLAB/
cd /Users/mmphill2/datasets/MarineHeatWaves/

slices = {'2010s','2020s','2030s','2040s','2050s','2060s','2070s','2080s','2090s'};
slicestart = [2010 2020 2030 2040 2050 2060 2070 2080 2090];
sliceend = [2019 2029 2039 2049 2059 2069 2079 2089 2099];

allgcm = {'ACCESS1-0','GISS-E2-R','NorESM1-M','ACCESS-CM2','CanESM5','GFDL-CM4','GISS_E213','IPSL-CM6A-LR','MRI-ESM2-0','NorESM2-LM'};
allensemble = {'r1i1p1','r6i1p1','r1i1p1','r1i1p1f1','r1i1p1f1','r1i1p1f1','none','r1i1p1f1','r1i1p1f1','r1i1p1f1'};

cmip5gcm = {'ACCESS1-0','GISS-E2-R','NorESM1-M'};
cmip6gcm = {'ACCESS-CM2','CanESM5','GFDL-CM4','GISS_E213','IPSL-CM6A-LR','MRI-ESM2-0','NorESM2-LM'};

allscen = {'ssp245','ssp585'};
allscennice = {'SSP245','SSP585'};

cmipmask = importdata(['OceanMask_025degree.mat']);

targetlats = [-89.75:0.25:90];
targetlons = [0.25:0.25:360];
targetlat = ones(length(targetlats),length(targetlons))*NaN;
targetlon = ones(length(targetlats),length(targetlons))*NaN;
for ii=1:length(targetlons),
   targetlat(:,ii) = targetlats;
end;
for jj=1:length(targetlats),
   targetlon(jj,:) = targetlons';
end;

slat = 8.5;
nlat = 31.5;
wlon = 261.5;
elon = 293.5;

sind = dsearchn(targetlat(:,1),slat);
nind = dsearchn(targetlat(:,1),nlat);
wind = dsearchn((targetlon(1,:)'),wlon);
eind = dsearchn((targetlon(1,:)'),elon);

latidx = [min(sind,nind):max(sind,nind)];
lonidx = [min(wind,eind):max(wind,eind)];

sublat = targetlat(latidx,lonidx);
sublon = targetlon(latidx,lonidx);

clear targetlat
clear targetlon
clear targetlats
clear targetlons

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      GISS E213
%%%    delta September SSTs for 2050s SSP585

histyears = 1980:2014;
futyears = 2050:2059;

allyears = 1980:2100;

thisgcm = 'GISS_E213';
gcmtitle = 'GISS E213';

targetlats = [-89.75:0.25:90];
targetlons = [-179.75:0.25:180];
targetlat = ones(length(targetlats),length(targetlons))*NaN;
targetlon = ones(length(targetlats),length(targetlons))*NaN;
for ii=1:length(targetlons),
   targetlat(:,ii) = targetlats;
end;
for jj=1:length(targetlats),
   targetlon(jj,:) = targetlons';
end;

slat = 8.5;
nlat = 31.5;
wlon = 261.5-360;
elon = 293.5-360;

sind = dsearchn(targetlat(:,1),slat);
nind = dsearchn(targetlat(:,1),nlat);
wind = dsearchn((targetlon(1,:)'),wlon);
eind = dsearchn((targetlon(1,:)'),elon);

latidx = [min(sind,nind):max(sind,nind)];
lonidx = [min(wind,eind):max(wind,eind)];

sublat = targetlat(latidx,lonidx);
sublon = targetlon(latidx,lonidx);

clear targetlat
clear targetlon
clear targetlats
clear targetlons

regionslat = 14.21;
regionnlat = 23.55;
regionwlon = -91.20;
regionelon = -79.54;

histgissdata = ones(length(latidx),length(lonidx),length(histyears),365)*NaN;
for xyear=1:length(histyears),
    thisyear = histyears(xyear);
    if (thisyear<2015),
        filename = ['tos_MAR_HR_day_CMIP6_' thisgcm '_historical_' num2str(thisyear) '.mat'];
    else
        filename = ['tos_MAR_HR_day_CMIP6_' thisgcm '_ssp585_' num2str(thisyear) '.mat'];
    end;
    subdata = importdata(filename);
    histgissdata(:,:,xyear,:) = subdata;
end;

september = [244:273];
histseptdata = squeeze(nanmean(squeeze(histgissdata(:,:,:,september)),4));
clear histgissdata
histmap = squeeze(nanmean(squeeze(histseptdata(:,:,1:30)),3));
clear histseptdata

futgissdata = ones(length(latidx),length(lonidx),length(futyears),365)*NaN;
for xyear=1:length(futyears),
    thisyear = futyears(xyear);
    if (thisyear<2015),
        filename = ['tos_MAR_HR_day_CMIP6_' thisgcm '_historical_' num2str(thisyear) '.mat'];
    else
        filename = ['tos_MAR_HR_day_CMIP6_' thisgcm '_ssp585_' num2str(thisyear) '.mat'];
    end;
    subdata = importdata(filename);
    futgissdata(:,:,xyear,:) = subdata;
end;

september = [244:273];
futseptdata = squeeze(nanmean(squeeze(futgissdata(:,:,:,september)),4));
clear futgissdata
futmap = squeeze(nanmean(squeeze(futseptdata(:,:,1:10)),3));
clear futseptdata

deltamap = futmap-histmap;

p = importdata('/Users/mmphill2/Documents/MATLAB/cmap_yellowred.mat');
p2 = squeeze(p(1:3:end,:));
p2(1:2,:) = [];

f = figure; colormap(p);
acr_pcolormapr4_nocoast(deltamap.*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[0.5 3],[gcmtitle ' CMIP6 Change in Mean September SSTs 2050-2059 SSP585'],['degrees C']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',[thisgcm '_CMIP6_DeltaMeanSST_September_2050-2059_SSP585.png']);
print(f,'-depsc',[thisgcm '_CMIP6_DeltaMeanSST_September_2050-2059_SSP585.eps']);

close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      GISS E213
%%%    unadjusted DHW for 2050s SSP585

thisgcm = 'GISS_E213';
gcmtitle = 'GISS E213';

maxmonth = importdata([thisgcm '_MesoAmericanReef_MaxMeanMonth_1985-1993.mat']);

allyears = [2050:2059];
nhalfweeks = round(((length(allyears)*365)/7)*2);

for xyear=1:length(allyears),
   thisyear = allyears(xyear);

   mmstart = [1 32 60 91 121 152 182 213 244 274 305 335];
   mmlength = [31 28 31 30 31 30 31 31 30 31 30 31];
   yearlength = 365;
   filename = ['tos_MAR_HR_day_CMIP6_' thisgcm '_ssp585_' num2str(thisyear) '.mat'];
   subdata = importdata([filename]);
   halfweekidx = [1:3.5:size(subdata,3)];
   yyhalfweeklydata = ones(length(latidx),length(lonidx),length(halfweekidx))*NaN;
   counter = 1;

   for xweek=1:length(halfweekidx),
      halfweekstart = floor(halfweekidx(xweek));
      halfweekend = floor(halfweekidx(xweek))+3;
      if (halfweekend>size(subdata,3)),
          halfweekend = size(subdata,3);
      end;
      if ((halfweekend-halfweekstart)<1),
          ;
      else
          yyhalfweeklydata(:,:,counter) = (squeeze(nanmean(squeeze(subdata(:,:,halfweekstart:halfweekend)),3)));
          counter=counter+1;
      end;
   end;
   clear subdata
   save([thisgcm '_MesoAmericanReef_HalfWeeklySSTs_' num2str(thisyear) '_NonAdjusted.mat'],'yyhalfweeklydata');

   %%% Calculate hotspots by subtracting max monthly SST
   %%% climatology from each twice-weekly SST

   allhotspots = yyhalfweeklydata*NaN;
   for xweek=1:size(yyhalfweeklydata,3),
      allhotspots(:,:,xweek) = yyhalfweeklydata(:,:,xweek)-maxmonth;
   end;
   allhotspots(allhotspots < 0) = NaN;
   clear yyhalfweeklydata
   save([thisgcm '_MesoAmericanReef_HalfWeeklyHotSpots_' num2str(thisyear) '_NonAdjusted.mat'],'allhotspots');

   %%% Calculate Degree Heating Weeks using the following equation
   %%% DHW = 0.5 * summation of previous 24 twice-weekly hotspots > 1 degree

   alldhw = allhotspots*NaN;

   if (xyear>1),

      lastyear = thisyear-1;
      lastyyhotspots = importdata([thisgcm '_MesoAmericanReef_HalfWeeklyHotSpots_' num2str(lastyear) '_NonAdjusted.mat']);

      for xweek=1:size(allhotspots,3),
          startdate = (xweek-23);
          halfweeklyidx = [startdate:xweek];
          if (xweek<24),
            weeklymeananom = ones(length(latidx),length(lonidx),length(halfweeklyidx))*NaN;
            for iweek=1:length(halfweeklyidx),
                 if (halfweeklyidx(iweek)<1),
                    lastyyidx = size(lastyyhotspots,3)+halfweeklyidx(iweek);
                    weeklymeananom(:,:,iweek) = squeeze(lastyyhotspots(:,:,lastyyidx));
                 else
                    weeklymeananom(:,:,iweek) = squeeze(allhotspots(:,:,halfweeklyidx(iweek)));
                 end;
            end;
          else
              weeklymeananom = squeeze(allhotspots(:,:,halfweeklyidx));
          end;
          for ii=1:length(lonidx),
              for jj=1:length(latidx),
                 if (length(find(squeeze(weeklymeananom(jj,ii,:))>1))>0),
                    alldhw(jj,ii,xweek) = (nansum(squeeze(weeklymeananom(jj,ii,find(squeeze(weeklymeananom(jj,ii,:))>1)))))/2;
                 end;
              end;
          end;
      end;

    else      %%% if it is first year, start at week 12

         for xweek=24:size(allhotspots,3),
            startdate = (xweek-23);
            halfweeklyidx = [startdate:xweek];
            weeklymeananom = squeeze(allhotspots(:,:,halfweeklyidx));
            for ii=1:length(lonidx),
               for jj=1:length(latidx),
                  if (length(find(squeeze(weeklymeananom(jj,ii,:))>1))>0),
                     alldhw(jj,ii,xweek) = (nansum(squeeze(weeklymeananom(jj,ii,find(squeeze(weeklymeananom(jj,ii,:))>1)))))/2;
      	           end;
               end;
            end;
         end;

    end;
    save([thisgcm '_MesoAmericanReef_HalfWeeklyDHW_' num2str(thisyear) '_NonAdjusted.mat'],'alldhw');

    clear allhotspots
    clear alldhw
    clear lastyyhotspots

end;

annualmaxdhw = ones(length(latidx),length(lonidx),length(allyears))*NaN;
for xyear=1:length(allyears),
  thisyear = allyears(xyear);
  thisyydhw = importdata([thisgcm '_MesoAmericanReef_HalfWeeklyDHW_' num2str(thisyear) '_NonAdjusted.mat']);
  annualmaxdhw(:,:,xyear) = nanmax(permute(thisyydhw,[3 1 2]));
  clear thisyydhw
end;
allmeanmaxdhw = squeeze(nanmean(annualmaxdhw,3));
clear annualmaxdhw
allmeanmaxdhw(isnan(allmeanmaxdhw)) = 0;

hexcolors = importdata('dhwpaper_colorbar21.mat');
p = ones(length(hexcolors),3)*NaN;
for rgb=1:size(p,1),
   p(rgb,:) = hex2rgb(hexcolors(rgb));
end;

f = figure; colormap(p);
acr_pcolormapr4_nocoast(allmeanmaxdhw.*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[10.5 31.5],['GISS E213 Non-Adjusted Mean Annual Max Degree Heating Weeks SSP585 2050-2059'],['degree heating weeks (degrees C)']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['GISS_E213_ssp585_2050-2059_AvgMaxDHW_NonAdjusted.png']);
print(f,'-depsc',['GISS_E213_ssp585_2050-2059_AvgMaxDHW_NonAdjusted.eps']);
close all
save(['GISS_E213_ssp585_2050-2059_AvgMaxDHW_NonAdjusted.mat'],'allmeanmaxdhw');

clear allmeanmaxdhw

%% make MERRA-2 adjusted version with this colorbar bounds

allmeanmaxdhw = importdata(['GISS_E213_MERRA2_ssp585_2050-2059_AvgMaxDHW.mat']);
f = figure; colormap(p);
acr_pcolormapr4_nocoast(allmeanmaxdhw.*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[10.5 31.5],['GISS E213 (MERRA-2 Adjusted) Mean Annual Max Degree Heating Weeks SSP585 2050-2059'],['degree heating weeks (degrees C)']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['GISS_E213_ssp585_2050-2059_AvgMaxDHW_Adjusted.png']);
print(f,'-depsc',['GISS_E213_ssp585_2050-2059_AvgMaxDHW_Adjusted.eps']);
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Ensemble Mean
%%%    2050s mean annual max DHW change for 45 and 85

targetlats = [-89.75:0.25:90];
targetlons = [0.25:0.25:360];
targetlat = ones(length(targetlats),length(targetlons))*NaN;
targetlon = ones(length(targetlats),length(targetlons))*NaN;
for ii=1:length(targetlons),
   targetlat(:,ii) = targetlats;
end;
for jj=1:length(targetlats),
   targetlon(jj,:) = targetlons';
end;

slat = 8.5;
nlat = 31.5;
wlon = 261.5;
elon = 293.5;

sind = dsearchn(targetlat(:,1),slat);
nind = dsearchn(targetlat(:,1),nlat);
wind = dsearchn((targetlon(1,:)'),wlon);
eind = dsearchn((targetlon(1,:)'),elon);

latidx = [min(sind,nind):max(sind,nind)];
lonidx = [min(wind,eind):max(wind,eind)];

sublat = targetlat(latidx,lonidx);
sublon = targetlon(latidx,lonidx);

clear targetlat
clear targetlon
clear targetlats
clear targetlons

regionslat = 14.21;
regionnlat = 23.55;
regionwlon = (-91.20)+360;
regionelon = (-79.54)+360;

allgcmdhw = ones(length(latidx),length(lonidx),length(allgcm),2)*NaN;
cmip5gcmdhw = ones(length(latidx),length(lonidx),length(cmip5gcm),2)*NaN;
cmip6gcmdhw = ones(length(latidx),length(lonidx),length(cmip6gcm),2)*NaN;

for xgcm=1:length(cmip5gcm),
    thisgcm = allgcm{xgcm};
    data = importdata([thisgcm '_MERRA2_rcp45_2050-2059_DeltaAvgMaxDHW.mat']);
    cmip5gcmdhw(:,:,xgcm,1) = data;
    allgcmdhw(:,:,xgcm,1) = data;
    clear data
    data = importdata([thisgcm '_MERRA2_rcp85_2050-2059_DeltaAvgMaxDHW.mat']);
    cmip5gcmdhw(:,:,xgcm,2) = data;
    allgcmdhw(:,:,xgcm,2) = data;
    clear data
end;
for xgcm=(length(cmip5gcm)+1):length(allgcm),
    thisgcm = allgcm{xgcm};
    data = importdata([thisgcm '_MERRA2_ssp245_2050-2059_DeltaAvgMaxDHW.mat']);
    cmip6gcmdhw(:,:,(xgcm-(length(cmip5gcm))),1) = data;
    allgcmdhw(:,:,xgcm,1) = data;
    clear data
    data = importdata([thisgcm '_MERRA2_ssp585_2050-2059_DeltaAvgMaxDHW.mat']);
    cmip6gcmdhw(:,:,(xgcm-(length(cmip5gcm))),2) = data;
    allgcmdhw(:,:,xgcm,2) = data;
    clear data
end;

%%%%%% maps of ensemble DHW change

cmip5map = squeeze(nanmean(cmip5gcmdhw,3));
cmip6map = squeeze(nanmean(cmip6gcmdhw,3));
map = squeeze(nanmean(allgcmdhw,3));

hexcolors = importdata('dhwpaper_colorbar21.mat');
p = ones(length(hexcolors),3)*NaN;
for rgb=1:size(p,1),
   p(rgb,:) = hex2rgb(hexcolors(rgb));
end;

f = figure; colormap(p);
acr_pcolormapr4_nocoast(squeeze(cmip5map(:,:,1)).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[7.5 28.5],['CMIP5 Ensemble Mean Change in Mean Annual Max DHW RCP45 2050-2059'],['change in mean annual max degree heating weeks (degrees C) from baseline']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['CMIP5_EnsembleMean_RCP45_2050-2059_DeltaAvgMaxDHW_subset.png']);
print(f,'-depsc',['CMIP5_EnsembleMean_RCP45_2050-2059_DeltaAvgMaxDHW_subset.eps']);
close all

f = figure; colormap(p);
acr_pcolormapr4_nocoast(squeeze(cmip5map(:,:,2)).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[7.5 28.5],['CMIP5 Ensemble Mean Change in Mean Annual Max DHW RCP85 2050-2059'],['change in mean annual max degree heating weeks (degrees C) from baseline']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['CMIP5_EnsembleMean_RCP85_2050-2059_DeltaAvgMaxDHW_subset.png']);
print(f,'-depsc',['CMIP5_EnsembleMean_RCP85_2050-2059_DeltaAvgMaxDHW_subset.eps']);
close all

f = figure; colormap(p);
acr_pcolormapr4_nocoast(squeeze(cmip6map(:,:,1)).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[7.5 28.5],['CMIP6 Ensemble Mean Change in Mean Annual Max DHW SSP245 2050-2059'],['change in mean annual max degree heating weeks (degrees C) from baseline']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['CMIP6_EnsembleMean_SSP245_2050-2059_DeltaAvgMaxDHW_subset.png']);
print(f,'-depsc',['CMIP6_EnsembleMean_SSP245_2050-2059_DeltaAvgMaxDHW_subset.eps']);
close all

f = figure; colormap(p);
acr_pcolormapr4_nocoast(squeeze(cmip6map(:,:,2)).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[7.5 28.5],['CMIP6 Ensemble Mean Change in Mean Annual Max DHW SSP585 2050-2059'],['change in mean annual max degree heating weeks (degrees C) from baseline']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['CMIP6_EnsembleMean_SSP585_2050-2059_DeltaAvgMaxDHW_subset.png']);
print(f,'-depsc',['CMIP6_EnsembleMean_SSP585_2050-2059_DeltaAvgMaxDHW_subset.eps']);
close all

f = figure; colormap(p);
acr_pcolormapr4_nocoast(squeeze(map(:,:,1)).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[7.5 28.5],['Ensemble Mean Change in Mean Annual Max DHW RCP45/SSP245 2050-2059'],['change in mean annual max degree heating weeks (degrees C) from baseline']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['EnsembleMean_SSP245_2050-2059_DeltaAvgMaxDHW_subset.png']);
print(f,'-depsc',['EnsembleMean_SSP245_2050-2059_DeltaAvgMaxDHW_subset.eps']);
close all

f = figure; colormap(p);
acr_pcolormapr4_nocoast(squeeze(map(:,:,2)).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[7.5 28.5],['Ensemble Mean Change in Mean Annual Max DHW RCP85/SSP585 2050-2059'],['change in mean annual max degree heating weeks (degrees C) from baseline']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['EnsembleMean_SSP585_2050-2059_DeltaAvgMaxDHW_subset.png']);
print(f,'-depsc',['EnsembleMean_SSP585_2050-2059_DeltaAvgMaxDHW_subset.eps']);
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% individual GCMs

hexcolors = importdata('dhwpaper_colorbar21.mat');
p = ones(length(hexcolors),3)*NaN;
for rgb=1:size(p,1),
   p(rgb,:) = hex2rgb(hexcolors(rgb));
end;

for xgcm=1:length(allgcm),
   thisgcm = allgcm{xgcm};
   if (strcmp(thisgcm,'GISS_E213')),
      gcmtitle = 'GISS E213';
   else
      gcmtitle = thisgcm;
   end;

   if (xgcm<4),
      thisscen = 'rcp45';
      thisscennice = 'RCP45';
   else
      thisscen = 'ssp245';
      thisscennice = 'SSP245';
   end;

   thisgcmdhw = squeeze(allgcmdhw(:,:,xgcm,1));
   f = figure; colormap(p);
   acr_pcolormapr4_nocoast(thisgcmdhw.*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[7.5 28.5],['MERRA2 (' gcmtitle ' Adjusted) Change in Mean Annual Max DHW ' thisscennice ' 2050-2059'],['change in mean max degree heating weeks (degrees C) from baseline']);
   hold on;
   h1 = borders('countries','color','k','linewidth',1);
   hold on;
   plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
   print(f,'-dpng',[thisgcm '_MERRA2_' thisscen '_2050-2059_DeltaAvgMaxDHW.png']);
   print(f,'-depsc',[thisgcm '_MERRA2_' thisscen '_2050-2059_DeltaAvgMaxDHW.eps']);
   close all

   if (xgcm<4),
      thisscen = 'rcp85';
      thisscennice = 'RCP85';
   else
      thisscen = 'ssp585';
      thisscennice = 'SSP585';
   end;

   thisgcmdhw = squeeze(allgcmdhw(:,:,xgcm,2));
   f = figure; colormap(p);
   acr_pcolormapr4_nocoast(thisgcmdhw.*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[7.5 28.5],['MERRA2 (' gcmtitle ' Adjusted) Change in Mean Annual Max DHW ' thisscennice ' 2050-2059'],['change in mean max degree heating weeks (degrees C) from baseline']);
   hold on;
   h1 = borders('countries','color','k','linewidth',1);
   hold on;
   plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
   print(f,'-dpng',[thisgcm '_MERRA2_' thisscen '_2050-2059_DeltaAvgMaxDHW.png']);
   print(f,'-depsc',[thisgcm '_MERRA2_' thisscen '_2050-2059_DeltaAvgMaxDHW.eps']);
   close all

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Figure 1b
%%%    ensemble mean baseline DHW

allgcmbasedhw = ones(length(latidx),length(lonidx),length(allgcm))*NaN;
cmip5gcmbasedhw = ones(length(latidx),length(lonidx),length(cmip5gcm))*NaN;
cmip6gcmbasedhw = ones(length(latidx),length(lonidx),length(cmip6gcm))*NaN;

for xgcm=1:length(cmip5gcm),
    thisgcm = allgcm{xgcm};
    data = importdata([thisgcm '_1980-2009_AvgMaxDHW.mat']);
    cmip5gcmbasedhw(:,:,xgcm) = data;
    clear data
end;
for xgcm=(length(cmip5gcm)+1):length(allgcm),
    thisgcm = allgcm{xgcm};
    data = importdata([thisgcm '_1980-2009_AvgMaxDHW.mat']);
    cmip6gcmbasedhw(:,:,(xgcm-(length(cmip5gcm)))) = data;
    clear data
end;
for xgcm=1:length(allgcm),
    thisgcm = allgcm{xgcm};
    data = importdata([thisgcm '_1980-2009_AvgMaxDHW.mat']);
    allgcmbasedhw(:,:,xgcm) = data;
    clear data
end;

cmip5map = squeeze(nanmean(cmip5gcmbasedhw,3));
cmip6map = squeeze(nanmean(cmip6gcmbasedhw,3));
map = squeeze(nanmean(allgcmbasedhw,3));

p = ones(7,3)*NaN;
p(1,:) = [200 200 200];
p(2,:) = [0 0 128];
p(3,:) = [0 130 200];
p(4,:) = [60 180 75];
p(5,:) = [255 255 25];
p(6,:) = [245 130 48];
p(7,:) = [128 0 0];
p = p/255;

f = figure; colormap(p);
acr_pcolormapr4_nocoast(cmip5map.*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[-0.5 6.5],['CMIP5 Ensemble Mean Annual Max Degree Heating Weeks 1980-2009'],['degree heating weeks (degrees C)']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['CMIP5_EnsembleMean_1980-2009_AvgMaxDHW_subset.png']);
print(f,'-depsc',['CMIP5_EnsembleMean_1980-2009_AvgMaxDHW_subset.eps']);
close all

f = figure; colormap(p);
acr_pcolormapr4_nocoast(cmip6map.*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[-0.5 6.5],['CMIP6 Ensemble Mean Annual Max Degree Heating Weeks 1980-2009'],['degree heating weeks (degrees C)']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['CMIP6_EnsembleMean_1980-2009_AvgMaxDHW_subset.png']);
print(f,'-depsc',['CMIP6_EnsembleMean_1980-2009_AvgMaxDHW_subset.eps']);
close all

f = figure; colormap(p);
acr_pcolormapr4_nocoast(map.*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[-0.5 6.5],['Ensemble Mean Annual Max Degree Heating Weeks 1980-2009'],['degree heating weeks (degrees C)']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['EnsembleMean_1980-2009_AvgMaxDHW_subset.png']);
print(f,'-depsc',['EnsembleMean_1980-2009_AvgMaxDHW_subset.eps']);
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% individual GCMs

p = ones(7,3)*NaN;
p(1,:) = [200 200 200];
p(2,:) = [0 0 128];
p(3,:) = [0 130 200];
p(4,:) = [60 180 75];
p(5,:) = [255 255 25];
p(6,:) = [245 130 48];
p(7,:) = [128 0 0];
p = p/255;

for xgcm=1:length(allgcm),
   thisgcm = allgcm{xgcm};
   if (strcmp(thisgcm,'GISS_E213')),
      gcmtitle = 'GISS E213';
   else
      gcmtitle = thisgcm;
   end;

   thisgcmdhw = squeeze(allgcmbasedhw(:,:,xgcm));

   f = figure; colormap(p);
   acr_pcolormapr4_nocoast(thisgcmdhw.*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[-0.5 6.5],[gcmtitle ' Annual Max Degree Heating Weeks 1980-2009'],['degree heating weeks (degrees C)']);
   hold on;
   h1 = borders('countries','color','k','linewidth',1);
   hold on;
   plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
   print(f,'-dpng',[thisgcm '_1980-2009_AvgMaxDHW.png']);
   print(f,'-depsc',[thisgcm '_1980-2009_AvgMaxDHW.eps']);
   close all

end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Figure 4a, 4b, 4c, 4d --- RCP45
%%%    ensemble mean change in frequency/duration of hotspots

allgcmfreq = ones(length(latidx),length(lonidx),length(allgcm))*NaN;
cmip5gcmfreq = ones(length(latidx),length(lonidx),3)*NaN;
cmip6gcmfreq = ones(length(latidx),length(lonidx),7)*NaN;

for xgcm=1:3,
    thisgcm = allgcm{xgcm};
    data = importdata([thisgcm '_MERRA2_HotSpots_DeltaFrequency_rcp45_2050-2059.mat']);
    cmip5gcmfreq(:,:,xgcm) = data;
    allgcmfreq(:,:,xgcm) = data;
    clear data
end;
for xgcm=4:length(allgcm),
    thisgcm = allgcm{xgcm};
    data = importdata([thisgcm '_MERRA2_HotSpots_DeltaFrequency_ssp245_2050-2059.mat']);
    cmip6gcmfreq(:,:,(xgcm-3)) = data;
    allgcmfreq(:,:,xgcm) = data;
    clear data
end;

cmip5map = squeeze(nanmean(cmip5gcmfreq,3));
cmip6map = squeeze(nanmean(cmip6gcmfreq,3));
map = squeeze(nanmean(allgcmfreq,3));

hexcolors = importdata('dhwpaper_tealmagenta32.mat');
p = ones(length(hexcolors),3)*NaN;
for rgb=1:size(p,1),
   p(rgb,:) = hex2rgb(hexcolors(rgb));
end;

hexcolors = importdata('dhwpaper_colorbar32.mat');
p2 = ones(length(hexcolors),3)*NaN;
for rgb=1:size(p2,1),
   p2(rgb,:) = hex2rgb(hexcolors(rgb));
end;

f = figure; colormap(p);
acr_pcolormapr4_nocoast((cmip5map/2).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[3.5 35.5],['CMIP5 Ensemble Mean Change in Annual Frequency of Hotspots>1 RCP45 2050-2059'],['change in number of weeks with >1C anomaly from baseline']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['CMIP5_EnsembleMean_RCP45_2050-2059_HotSpots_DeltaFrequency_subset.png']);
print(f,'-depsc',['CMIP5_EnsembleMean_RCP45_2050-2059_HotSpots_DeltaFrequency_subset.eps']);
close all

f = figure; colormap(p2);
acr_pcolormapr4_nocoast((cmip5map/2).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[3.5 35.5],['CMIP5 Ensemble Mean Change in Annual Frequency of Hotspots>1 RCP45 2050-2059'],['change in number of weeks with >1C anomaly from baseline']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['CMIP5_EnsembleMean_RCP45_2050-2059_HotSpots_DeltaFrequency_subset_old.png']);
print(f,'-depsc',['CMIP5_EnsembleMean_RCP45_2050-2059_HotSpots_DeltaFrequency_subset_old.eps']);
close all

f = figure; colormap(p);
acr_pcolormapr4_nocoast((cmip6map/2).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[3.5 35.5],['CMIP6 Ensemble Mean Change in Annual Frequency of Hotspots>1 RCP45 2050-2059'],['change in number of weeks with >1C anomaly from baseline']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['CMIP6_EnsembleMean_RCP45_2050-2059_HotSpots_DeltaFrequency_subset.png']);
print(f,'-depsc',['CMIP6_EnsembleMean_RCP45_2050-2059_HotSpots_DeltaFrequency_subset.eps']);
close all

f = figure; colormap(p2);
acr_pcolormapr4_nocoast((cmip6map/2).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[3.5 35.5],['CMIP6 Ensemble Mean Change in Annual Frequency of Hotspots>1 RCP45 2050-2059'],['change in number of weeks with >1C anomaly from baseline']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['CMIP6_EnsembleMean_RCP45_2050-2059_HotSpots_DeltaFrequency_subset_old.png']);
print(f,'-depsc',['CMIP6_EnsembleMean_RCP45_2050-2059_HotSpots_DeltaFrequency_subset_old.eps']);
close all

f = figure; colormap(p);
acr_pcolormapr4_nocoast((map/2).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[3.5 35.5],['Ensemble Mean Change in Annual Frequency of Hotspots>1 RCP45/SSP245 2050-2059'],['change in number of weeks with >1C anomaly from baseline']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['EnsembleMean_SSP245_2050-2059_HotSpots_DeltaFrequency_subset.png']);
print(f,'-depsc',['EnsembleMean_SSP245_2050-2059_HotSpots_DeltaFrequency_subset.eps']);
close all

f = figure; colormap(p2);
acr_pcolormapr4_nocoast((map/2).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[3.5 35.5],['Ensemble Mean Change in Annual Frequency of Hotspots>1 RCP45/SSP245 2050-2059'],['change in number of weeks with >1C anomaly from baseline']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['EnsembleMean_SSP245_2050-2059_HotSpots_DeltaFrequency_subset_old.png']);
print(f,'-depsc',['EnsembleMean_SSP245_2050-2059_HotSpots_DeltaFrequency_subset_old.eps']);
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% individual GCMs

hexcolors = importdata('dhwpaper_tealmagenta32.mat');
p = ones(length(hexcolors),3)*NaN;
for rgb=1:size(p,1),
   p(rgb,:) = hex2rgb(hexcolors(rgb));
end;

for xgcm=1:length(allgcm),
   thisgcm = allgcm{xgcm};
   if (strcmp(thisgcm,'GISS_E213')),
      gcmtitle = 'GISS E213';
   else
      gcmtitle = thisgcm;
   end;

   if (xgcm<4),
      thisscen = 'rcp45';
      thisscennice = 'RCP45';
   else
      thisscen = 'ssp245';
      thisscennice = 'SSP245';
   end;

   thisgcmdhw = squeeze(allgcmfreq(:,:,xgcm));
   f = figure; colormap(p);
   acr_pcolormapr4_nocoast((thisgcmdhw/2).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[3.5 35.5],['MERRA2 (' gcmtitle ' Adjusted) Change in Mean Annual Frequency of Hotspots>1C ' thisscennice ' 2050-2059'],['change in number of weeks with >1C anomaly from baseline']);
   hold on;
   h1 = borders('countries','color','k','linewidth',1);
   hold on;
   plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
   print(f,'-dpng',[thisgcm '_MERRA2_HotSpots_DeltaFrequency_' thisscen '_2050-2059.png']);
   print(f,'-depsc',[thisgcm '_MERRA2_HotSpots_DeltaFrequency_' thisscen '_2050-2059.eps']);
   close all

end;

allgcmduration = ones(length(latidx),length(lonidx),length(allgcm))*NaN;
cmip5gcmduration = ones(length(latidx),length(lonidx),3)*NaN;
cmip6gcmduration = ones(length(latidx),length(lonidx),7)*NaN;

for xgcm=1:3,
    thisgcm = allgcm{xgcm};
    data = importdata([thisgcm '_MERRA2_HotSpots_DeltaMaxDuration_rcp45_2050-2059.mat']);
    cmip5gcmduration(:,:,xgcm) = data;
    allgcmduration(:,:,xgcm) = data;
    clear data
end;
for xgcm=4:length(allgcm),
    thisgcm = allgcm{xgcm};
    data = importdata([thisgcm '_MERRA2_HotSpots_DeltaMaxDuration_ssp245_2050-2059.mat']);
    cmip6gcmduration(:,:,(xgcm-3)) = data;
    allgcmduration(:,:,xgcm) = data;
    clear data
end;

cmip5map = squeeze(nanmean(cmip5gcmduration,3));
cmip6map = squeeze(nanmean(cmip6gcmduration,3));
map = squeeze(nanmean(allgcmduration,3));

hexcolors = importdata('dhwpaper_tealmagenta24.mat');
p = ones(length(hexcolors),3)*NaN;
for rgb=1:size(p,1),
   p(rgb,:) = hex2rgb(hexcolors(rgb));
end;

hexcolors = importdata('dhwpaper_colorbar24.mat');
p2 = ones(length(hexcolors),3)*NaN;
for rgb=1:size(p2,1),
   p2(rgb,:) = hex2rgb(hexcolors(rgb));
end;

f = figure; colormap(p);
acr_pcolormapr4_nocoast((cmip5map/2).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[4.5 28.5],['CMIP5 Ensemble Mean Change in Annual MaxDuration of Hotspots>1 RCP45 2050-2059'],['change from baseline in maximum consecutive weeks with >1C anomaly']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['CMIP5_EnsembleMean_RCP45_2050-2059_HotSpots_DeltaMaxDuration_subset.png']);
print(f,'-depsc',['CMIP5_EnsembleMean_RCP45_2050-2059_HotSpots_DeltaMaxDuration_subset.eps']);
close all

f = figure; colormap(p2);
acr_pcolormapr4_nocoast((cmip5map/2).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[4.5 28.5],['CMIP5 Ensemble Mean Change in Annual MaxDuration of Hotspots>1 RCP45 2050-2059'],['change from baseline in maximum consecutive weeks with >1C anomaly']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['CMIP5_EnsembleMean_RCP45_2050-2059_HotSpots_DeltaMaxDuration_subset_old.png']);
print(f,'-depsc',['CMIP5_EnsembleMean_RCP45_2050-2059_HotSpots_DeltaMaxDuration_subset_old.eps']);
close all

f = figure; colormap(p);
acr_pcolormapr4_nocoast((cmip6map/2).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[4.5 28.5],['CMIP6 Ensemble Mean Change in Annual MaxDuration of Hotspots>1 SSP245 2050-2059'],['change from baseline in maximum consecutive weeks with >1C anomaly']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['CMIP6_EnsembleMean_SSP245_2050-2059_HotSpots_DeltaMaxDuration_subset.png']);
print(f,'-depsc',['CMIP6_EnsembleMean_SSP245_2050-2059_HotSpots_DeltaMaxDuration_subset.eps']);
close all

f = figure; colormap(p2);
acr_pcolormapr4_nocoast((cmip6map/2).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[4.5 28.5],['CMIP6 Ensemble Mean Change in Annual MaxDuration of Hotspots>1 SSP245 2050-2059'],['change from baseline in maximum consecutive weeks with >1C anomaly']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['CMIP6_EnsembleMean_SSP245_2050-2059_HotSpots_DeltaMaxDuration_subset_old.png']);
print(f,'-depsc',['CMIP6_EnsembleMean_SSP245_2050-2059_HotSpots_DeltaMaxDuration_subset_old.eps']);
close all

f = figure; colormap(p);
acr_pcolormapr4_nocoast((map/2).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[4.5 28.5],['Ensemble Mean Change in Annual MaxDuration of Hotspots>1 RCP45/SSP245 2050-2059'],['change from baseline in maximum consecutive weeks with >1C anomaly']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['EnsembleMean_SSP245_2050-2059_HotSpots_DeltaMaxDuration_subset.png']);
print(f,'-depsc',['EnsembleMean_SSP245_2050-2059_HotSpots_DeltaMaxDuration_subset.eps']);
close all

f = figure; colormap(p2);
acr_pcolormapr4_nocoast((map/2).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[4.5 28.5],['Ensemble Mean Change in Annual MaxDuration of Hotspots>1 RCP45/SSP245 2050-2059'],['change from baseline in maximum consecutive weeks with >1C anomaly']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['EnsembleMean_SSP245_2050-2059_HotSpots_DeltaMaxDuration_subset_old.png']);
print(f,'-depsc',['EnsembleMean_SSP245_2050-2059_HotSpots_DeltaMaxDuration_subset_old.eps']);
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% individual GCMs

hexcolors = importdata('dhwpaper_tealmagenta24.mat');
p = ones(length(hexcolors),3)*NaN;
for rgb=1:size(p,1),
   p(rgb,:) = hex2rgb(hexcolors(rgb));
end;

for xgcm=1:length(allgcm),
   thisgcm = allgcm{xgcm};
   if (strcmp(thisgcm,'GISS_E213')),
      gcmtitle = 'GISS E213';
   else
      gcmtitle = thisgcm;
   end;

   if (xgcm<4),
      thisscen = 'rcp45';
      thisscennice = 'RCP45';
   else
      thisscen = 'ssp245';
      thisscennice = 'SSP245';
   end;

   thisgcmdhw = squeeze(allgcmduration(:,:,xgcm));
   f = figure; colormap(p);
   acr_pcolormapr4_nocoast((thisgcmdhw/2).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[4.5 28.5],['MERRA2 (' gcmtitle ' Adjusted) Change in Mean Annual Max Duration of Hotspots>1C ' thisscennice ' 2050-2059'],['change from baseline in maximum consecutive weeks with >1C anomaly']);
   hold on;
   h1 = borders('countries','color','k','linewidth',1);
   hold on;
   plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
   print(f,'-dpng',[thisgcm '_MERRA2_HotSpots_DeltaMaxDuration_' thisscen '_2050-2059.png']);
   print(f,'-depsc',[thisgcm '_MERRA2_HotSpots_DeltaMaxDuration_' thisscen '_2050-2059.eps']);
   close all

end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Figure 4a, 4b, 4c, 4d -- RCP85
%%%    ensemble mean change in frequency/duration of hotspots

allgcmfreq = ones(length(latidx),length(lonidx),length(allgcm))*NaN;
cmip5gcmfreq = ones(length(latidx),length(lonidx),3)*NaN;
cmip6gcmfreq = ones(length(latidx),length(lonidx),7)*NaN;

for xgcm=1:3,
    thisgcm = allgcm{xgcm};
    data = importdata([thisgcm '_MERRA2_HotSpots_DeltaFrequency_rcp85_2050-2059.mat']);
    cmip5gcmfreq(:,:,xgcm) = data;
    allgcmfreq(:,:,xgcm) = data;
    clear data
end;
for xgcm=4:length(allgcm),
    thisgcm = allgcm{xgcm};
    data = importdata([thisgcm '_MERRA2_HotSpots_DeltaFrequency_ssp585_2050-2059.mat']);
    cmip6gcmfreq(:,:,(xgcm-3)) = data;
    allgcmfreq(:,:,xgcm) = data;
    clear data
end;

cmip5map = squeeze(nanmean(cmip5gcmfreq,3));
cmip6map = squeeze(nanmean(cmip6gcmfreq,3));
map = squeeze(nanmean(allgcmfreq,3));

%p2 = importdata('/home/mmphill2/matlab/matfiles/cmap_rainbow_lowgrey.mat');
%p3 = importdata('/home/mmphill2/matlab/matfiles/cmap_rainbow_long.mat');
%hexcolors = importdata('/home/mmphill2/matlab/matfiles/dhw_frequency_colors.mat');
%p4 = p3*NaN;
%for rgb=1:size(p4,1),
%   p4(rgb,:) = hex2rgb(hexcolors(rgb));
%end;

hexcolors = importdata('dhwpaper_tealmagenta32.mat');
p = ones(length(hexcolors),3)*NaN;
for rgb=1:size(p,1),
   p(rgb,:) = hex2rgb(hexcolors(rgb));
end;

hexcolors = importdata('dhwpaper_colorbar32.mat');
p2 = ones(length(hexcolors),3)*NaN;
for rgb=1:size(p2,1),
   p2(rgb,:) = hex2rgb(hexcolors(rgb));
end;

f = figure; colormap(p);
acr_pcolormapr4_nocoast((cmip5map/2).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[3.5 35.5],['CMIP5 Ensemble Mean Change in Annual Frequency of Hotspots>1 RCP85 2050-2059'],['change in number of weeks with >1C anomaly from baseline']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['CMIP5_EnsembleMean_RCP85_2050-2059_HotSpots_DeltaFrequency_subset.png']);
print(f,'-depsc',['CMIP5_EnsembleMean_RCP85_2050-2059_HotSpots_DeltaFrequency_subset.eps']);
close all

f = figure; colormap(p2);
acr_pcolormapr4_nocoast((cmip5map/2).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[3.5 35.5],['CMIP5 Ensemble Mean Change in Annual Frequency of Hotspots>1 RCP85 2050-2059'],['change in number of weeks with >1C anomaly from baseline']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['CMIP5_EnsembleMean_RCP85_2050-2059_HotSpots_DeltaFrequency_subset_old.png']);
print(f,'-depsc',['CMIP5_EnsembleMean_RCP85_2050-2059_HotSpots_DeltaFrequency_subset_old.eps']);
close all

f = figure; colormap(p);
acr_pcolormapr4_nocoast((cmip6map/2).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[3.5 35.5],['CMIP6 Ensemble Mean Change in Annual Frequency of Hotspots>1 SSP585 2050-2059'],['change in number of weeks with >1C anomaly from baseline']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['CMIP6_EnsembleMean_SSP585_2050-2059_HotSpots_DeltaFrequency_subset.png']);
print(f,'-depsc',['CMIP6_EnsembleMean_SSP585_2050-2059_HotSpots_DeltaFrequency_subset.eps']);
close all

f = figure; colormap(p2);
acr_pcolormapr4_nocoast((cmip6map/2).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[3.5 35.5],['CMIP6 Ensemble Mean Change in Annual Frequency of Hotspots>1 SSP585 2050-2059'],['change in number of weeks with >1C anomaly from baseline']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['CMIP6_EnsembleMean_SSP585_2050-2059_HotSpots_DeltaFrequency_subset_old.png']);
print(f,'-depsc',['CMIP6_EnsembleMean_SSP585_2050-2059_HotSpots_DeltaFrequency_subset_old.eps']);
close all

f = figure; colormap(p);
acr_pcolormapr4_nocoast((map/2).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[3.5 35.5],['Ensemble Mean Change in Annual Frequency of Hotspots>1 RCP85/SSP585 2050-2059'],['change in number of weeks with >1C anomaly from baseline']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['EnsembleMean_SSP585_2050-2059_HotSpots_DeltaFrequency_subset.png']);
print(f,'-depsc',['EnsembleMean_SSP585_2050-2059_HotSpots_DeltaFrequency_subset.eps']);
close all

f = figure; colormap(p2);
acr_pcolormapr4_nocoast((map/2).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[3.5 35.5],['Ensemble Mean Change in Annual Frequency of Hotspots>1 RCP85/SSP585 2050-2059'],['change in number of weeks with >1C anomaly from baseline']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['EnsembleMean_SSP585_2050-2059_HotSpots_DeltaFrequency_subset_old.png']);
print(f,'-depsc',['EnsembleMean_SSP585_2050-2059_HotSpots_DeltaFrequency_subset_old.eps']);
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% individual GCMs

hexcolors = importdata('dhwpaper_tealmagenta32.mat');
p = ones(length(hexcolors),3)*NaN;
for rgb=1:size(p,1),
   p(rgb,:) = hex2rgb(hexcolors(rgb));
end;

for xgcm=1:length(allgcm),
   thisgcm = allgcm{xgcm};
   if (strcmp(thisgcm,'GISS_E213')),
      gcmtitle = 'GISS E213';
   else
      gcmtitle = thisgcm;
   end;

   if (xgcm<4),
      thisscen = 'rcp85';
      thisscennice = 'RCP85';
   else
      thisscen = 'ssp585';
      thisscennice = 'SSP585';
   end;

   thisgcmdhw = squeeze(allgcmfreq(:,:,xgcm));
   f = figure; colormap(p);
   acr_pcolormapr4_nocoast((thisgcmdhw/2).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[3.5 35.5],['MERRA2 (' gcmtitle ' Adjusted) Change in Mean Annual Frequency of Hotspots>1C ' thisscennice ' 2050-2059'],['change in number of weeks with >1C anomaly from baseline']);
   hold on;
   h1 = borders('countries','color','k','linewidth',1);
   hold on;
   plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
   print(f,'-dpng',[thisgcm '_MERRA2_HotSpots_DeltaFrequency_' thisscen '_2050-2059.png']);
   print(f,'-depsc',[thisgcm '_MERRA2_HotSpots_DeltaFrequency_' thisscen '_2050-2059.eps']);
   close all

end;

allgcmduration = ones(length(latidx),length(lonidx),length(allgcm))*NaN;
cmip5gcmduration = ones(length(latidx),length(lonidx),3)*NaN;
cmip6gcmduration = ones(length(latidx),length(lonidx),7)*NaN;

for xgcm=1:3,
    thisgcm = allgcm{xgcm};
    data = importdata([thisgcm '_MERRA2_HotSpots_DeltaMaxDuration_rcp85_2050-2059.mat']);
    cmip5gcmduration(:,:,xgcm) = data;
    allgcmduration(:,:,xgcm) = data;
    clear data
end;
for xgcm=4:length(allgcm),
    thisgcm = allgcm{xgcm};
    data = importdata([thisgcm '_MERRA2_HotSpots_DeltaMaxDuration_ssp585_2050-2059.mat']);
    cmip6gcmduration(:,:,(xgcm-3)) = data;
    allgcmduration(:,:,xgcm) = data;
    clear data
end;

cmip5map = squeeze(nanmean(cmip5gcmduration,3));
cmip6map = squeeze(nanmean(cmip6gcmduration,3));
map = squeeze(nanmean(allgcmduration,3));

hexcolors = importdata('dhwpaper_tealmagenta24.mat');
p = ones(length(hexcolors),3)*NaN;
for rgb=1:size(p,1),
   p(rgb,:) = hex2rgb(hexcolors(rgb));
end;

hexcolors = importdata('dhwpaper_colorbar24.mat');
p2 = ones(length(hexcolors),3)*NaN;
for rgb=1:size(p2,1),
   p2(rgb,:) = hex2rgb(hexcolors(rgb));
end;

f = figure; colormap(p);
acr_pcolormapr4_nocoast((cmip5map/2).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[4.5 28.5],['CMIP5 Ensemble Mean Change in Annual MaxDuration of Hotspots>1 RCP85 2050-2059'],['change from baseline in maximum consecutive weeks with >1C anomaly']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['CMIP5_EnsembleMean_RCP85_2050-2059_HotSpots_DeltaMaxDuration_subset.png']);
print(f,'-depsc',['CMIP5_EnsembleMean_RCP85_2050-2059_HotSpots_DeltaMaxDuration_subset.eps']);
close all

f = figure; colormap(p2);
acr_pcolormapr4_nocoast((cmip5map/2).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[4.5 28.5],['CMIP5 Ensemble Mean Change in Annual MaxDuration of Hotspots>1 RCP85 2050-2059'],['change from baseline in maximum consecutive weeks with >1C anomaly']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['CMIP5_EnsembleMean_RCP85_2050-2059_HotSpots_DeltaMaxDuration_subset_old.png']);
print(f,'-depsc',['CMIP5_EnsembleMean_RCP85_2050-2059_HotSpots_DeltaMaxDuration_subset_old.eps']);
close all

f = figure; colormap(p);
acr_pcolormapr4_nocoast((cmip6map/2).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[4.5 28.5],['CMIP6 Ensemble Mean Change in Annual MaxDuration of Hotspots>1 SSP585 2050-2059'],['change from baseline in maximum consecutive weeks with >1C anomaly']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['CMIP6_EnsembleMean_SSP585_2050-2059_HotSpots_DeltaMaxDuration_subset.png']);
print(f,'-depsc',['CMIP6_EnsembleMean_SSP585_2050-2059_HotSpots_DeltaMaxDuration_subset.eps']);
close all

f = figure; colormap(p2);
acr_pcolormapr4_nocoast((cmip6map/2).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[4.5 28.5],['CMIP6 Ensemble Mean Change in Annual MaxDuration of Hotspots>1 SSP585 2050-2059'],['change from baseline in maximum consecutive weeks with >1C anomaly']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['CMIP6_EnsembleMean_SSP585_2050-2059_HotSpots_DeltaMaxDuration_subset_old.png']);
print(f,'-depsc',['CMIP6_EnsembleMean_SSP585_2050-2059_HotSpots_DeltaMaxDuration_subset_old.eps']);
close all

f = figure; colormap(p);
acr_pcolormapr4_nocoast((map/2).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[4.5 28.5],['Ensemble Mean Change in Annual MaxDuration of Hotspots>1 RCP85/SSP585 2050-2059'],['change from baseline in maximum consecutive weeks with >1C anomaly']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['EnsembleMean_SSP585_2050-2059_HotSpots_DeltaMaxDuration_subset.png']);
print(f,'-depsc',['EnsembleMean_SSP585_2050-2059_HotSpots_DeltaMaxDuration_subset.eps']);
close all

f = figure; colormap(p2);
acr_pcolormapr4_nocoast((map/2).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[4.5 28.5],['Ensemble Mean Change in Annual MaxDuration of Hotspots>1 RCP85/SSP585 2050-2059'],['change from baseline in maximum consecutive weeks with >1C anomaly']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['EnsembleMean_SSP585_2050-2059_HotSpots_DeltaMaxDuration_subset_old.png']);
print(f,'-depsc',['EnsembleMean_SSP585_2050-2059_HotSpots_DeltaMaxDuration_subset_old.eps']);
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% individual GCMs

hexcolors = importdata('dhwpaper_tealmagenta24.mat');
p = ones(length(hexcolors),3)*NaN;
for rgb=1:size(p,1),
   p(rgb,:) = hex2rgb(hexcolors(rgb));
end;

for xgcm=1:length(allgcm),
   thisgcm = allgcm{xgcm};
   if (strcmp(thisgcm,'GISS_E213')),
      gcmtitle = 'GISS E213';
   else
      gcmtitle = thisgcm;
   end;

   if (xgcm<4),
      thisscen = 'rcp85';
      thisscennice = 'RCP85';
   else
      thisscen = 'ssp585';
      thisscennice = 'SSP585';
   end;

   thisgcmdhw = squeeze(allgcmduration(:,:,xgcm));
   f = figure; colormap(p);
   acr_pcolormapr4_nocoast((thisgcmdhw/2).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[4.5 28.5],['MERRA2 (' gcmtitle ' Adjusted) Change in Mean Annual Max Duration of Hotspots>1C ' thisscennice ' 2050-2059'],['change from baseline in maximum consecutive weeks with >1C anomaly']);
   hold on;
   h1 = borders('countries','color','k','linewidth',1);
   hold on;
   plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
   print(f,'-dpng',[thisgcm '_MERRA2_HotSpots_DeltaMaxDuration_' thisscen '_2050-2059.png']);
   print(f,'-depsc',[thisgcm '_MERRA2_HotSpots_DeltaMaxDuration_' thisscen '_2050-2059.eps']);
   close all

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% maximum monthly delta value for each decade, model and scenario

maxdeltas = ones(length(allgcm),length(slices),2)*NaN;
deltaid = ones(length(allgcm),length(slices),2)*NaN;

for xgcm = 1:length(allgcm),
   if (xgcm<4),
      allscen = {'rcp45','rcp85'};
   else
      allscen = {'ssp245','ssp585'};
   end;
   thisgcm = allgcm{xgcm};
   regionaldeltas = ones(length(slices),length(allscen),12)*NaN;

   for xscen=1:length(allscen),
      thisscen = allscen{xscen};
      for tt=1:length(slices),
         thisslice = slices{tt};
         data = importdata([thisgcm '_' thisscen '_' thisslice '_RegionalDeltas.mat']);
         regionaldeltas(tt,xscen,:) = data;
         clear data
      end;
   end;

   save([thisgcm '_RegionalDeltas.mat'],'regionaldeltas');

   for xscen=1:length(allscen),
      for xslice=1:length(slices),
         md = nanmax(squeeze(regionaldeltas(xslice,xscen,:)));
         mdix = find(squeeze(regionaldeltas(xslice,xscen,:)) == nanmax(squeeze(regionaldeltas(xslice,xscen,:))));

         maxdeltas(xgcm,xslice,xscen) = md;
         deltaid(xgcm,xslice,xscen) = mdix;
      end;
   end;

end;
