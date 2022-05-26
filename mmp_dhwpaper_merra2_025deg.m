%	            mmp_dhwpaper_merra2_025deg
%
%  This script calculates degree heating weeks and other
%  corresponding metrics from MERRA-2 sea surface
%  temperature outputs for the WWF Mesoamerican Reef
%  project and corresponding paper with Model-E. This
%  version uses 0.25-degree resolution SST data.
%
%  -- Creates a map of MERRA-2 mean annual maximum DHW
%  value for the 1980-2009 period.
%
%  -- Creates a map of the maximum DHW value for the entire
%  year of 1998 as well as the maximum over just the time
%  period associated with bleaching (September 18 - October 1).
%
%  -- Creates a map of the maximum DHW value for the entire
%  year of 2005 as well as the maximum over just the time
%  period associated with bleaching (July 15 - November 15).
%
%  -- Creates a map of MERRA-2 mean annual maximum DHW value
%  for the 2010-2019 period.
%
%  -- Creates a map of the MERRA-2 climatological maximum mean
%  month in the 1985-1993 period (excluding 1991 and 1992), the
%  value from which all HotSpots are calculated.
%
%  -- Creates a map of the MERRA-2 mean annual frequency of HotSpots
%  greater than 1C for the 1980-2009 period.
%
%  -- Creates a map of the MERRA-2 mean annual maximum consecutive
%  duration of HotSpots greater than 1C for the 1980-2009 period.
%
%  -- Creates a map of the MERRA-2 mean annual frequency of HotSpots
%  greater than 1C for the 2010-2019 period.
%
%  -- Creates a map of the MERRA-2 mean annual maximum consecutive
%  duration of HotSpots greater than 1C for the 2010-2019 period.
%
%  -- Creates a map of MERRA-2 climatological mean SST for each
%  month in the 1980-2009 period.
%
%  -- Creates a map of MERRA-2 SSTs for each day of the September
%  18 - October 1 1998 bleaching event.
%
%  INPUTS: none
%  OUTPUTS: none
%
%
%
%		                             author: Meridel Phillips
%                                            mmp2192@columbia.edu
%				             date: 9/21/2020
%
function mmp_dhwpaper_merra2_025deg();
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

months = {'Jan','Feb','March','April','May','June','July','Aug','Sept','Oct','Nov','Dec'};
allmm = {'01','02','03','04','05','06','07','08','09','10','11','12'};

%%% regrid onto common grid for all SSTs

lats = [-89.5:0.5:90];
lons = [-179.5:0.5:180];
lat = ones(length(lats),length(lons))*NaN;
lon = ones(length(lats),length(lons))*NaN;
for ii=1:length(lons),
   lat(:,ii) = lats;
end;
for jj=1:length(lats),
   lon(jj,:) = lons';
end;

allyears = [1980:2019];

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

%for yy=1:length(allyears),
%   thisyear = allyears(yy);
%   data = importdata(['merra2_05deg_sst_global_daily_' num2str(thisyear) '.mat']);
%   regriddata = ones(size(sublat,1),size(sublon,2),size(data,3))*NaN;
%   for dd=1:size(data,3),
%      sample = griddata(lat,lon,(squeeze(data(:,:,dd))),sublat,sublon);
%      regriddata(:,:,dd) = inpaintn(sample);
%      clear sample
%      disp(['MERRA-2 ' num2str(thisyear) ' day' num2str(dd)])
%   end;
%   clear data
%   save(['merra2_sst_MAR_HR_daily_' num2str(thisyear) '.mat'],'regriddata');
%end;

oceanmask = importdata('OceanMask_025degree.mat');

regionslat = 14.21;
regionnlat = 23.55;
regionwlon = -91.20;
regionelon = -79.54;
londist = cosd(sublat(:,1))*(111.3215);
latdist = 111.3215;
area = sublon*NaN;
for ii=1:length(sublon(1,:)),
   area(:,ii) = londist*latdist;
end;

regionsind = dsearchn(sublat(:,1),regionslat);
regionnind = dsearchn(sublat(:,1),regionnlat);
regionwind = dsearchn(sublon(1,:)',regionwlon);
regioneind = dsearchn(sublon(1,:)',regionelon);

thislatidx = [min(regionsind,regionnind):max(regionsind,regionnind)];
thislonidx = [min(regionwind,regioneind):max(regionwind,regioneind)];
smalloceanmask = squeeze(oceanmask(thislatidx,thislonidx));

clear lat
clear lon

%%% for DHW calculation:
%%% find monthly mean SST climatologies from 1985 to 1993,
%%% excluding 1991 and 1992 (lat x lon x 12 months x 7 years)

allyears = [1985 1986 1987 1988 1989 1990 1993];
mm = {'01','02','03','04','05','06','07','08','09','10','11','12'};
dd = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31'};

mmsstclim = ones(length(latidx),length(lonidx),12,7)*NaN;

for xyear=1:length(allyears),
   thisyear = allyears(xyear);
   if (thisyear == 1988),
      mmstart = [1 32 61 92 122 153 183 214 245 275 306 336];
      mmlength = [31 29 31 30 31 30 31 31 30 31 30 31];
   else
      mmstart = [1 32 60 91 121 152 182 213 244 274 305 335];
      mmlength = [31 28 31 30 31 30 31 31 30 31 30 31];
   end;
   filename = ['merra2_sst_MAR_HR_daily_' num2str(thisyear) '.mat'];
   data = importdata([filename]);
   for xmonth=1:12,
      subdata = squeeze(data(:,:,(mmstart(xmonth):(mmstart(xmonth)+mmlength(xmonth)-1))));
      mmsstclim(:,:,xmonth,xyear) = squeeze(nanmean(subdata,3));
   end;
end;

clear data
clear subdata

%%% Next find montly mean SST climatology over
%%% entire 7-year period (lat x lon x 12 months)

meanmmsstclim = squeeze(nanmean(mmsstclim,4));
clear mmsstclim

%%% Calculate the maximum monthly mean SST climatology
%%% over entire 7-year period (lat x lon)

maxmonth = squeeze(nanmax(squeeze(permute(meanmmsstclim,[3 1 2]))));

save(['MERRA2_MesoAmericanReef_MaxMeanMonth_1985-1993.mat'],'maxmonth');
clear meanmmsstclim

maxmonth = importdata(['MERRA2_MesoAmericanReef_MaxMeanMonth_1985-1993.mat']);

%%% Calculate twice-weekly SST from 1980 to 2019
%%% not using MERRA2 pre-1980 so first 12 weeks of 1980 will not calculate DHW

allyears = [1980:2019];
leapyears = [1980:4:2019];
nhalfweeks = round((((length(allyears)*365)+length(leapyears))/7)*2);

for xyear=1:length(allyears),
   thisyear = allyears(xyear);
   subdata = importdata(['merra2_sst_MAR_HR_daily_' num2str(thisyear) '.mat']);

   if (size(subdata,3)>365),
      subdata(:,:,60) = [];
   end;

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
   save(['MERRA2_MesoAmericanReef_HalfWeeklySSTs_' num2str(thisyear) '.mat'],'yyhalfweeklydata');

   %%% Calculate hotspots by subtracting max monthly SST
   %%% climatology from each twice-weekly SST

   allhotspots = yyhalfweeklydata*NaN;
   for xweek=1:size(yyhalfweeklydata,3),
      allhotspots(:,:,xweek) = yyhalfweeklydata(:,:,xweek)-maxmonth;
   end;
   allhotspots(allhotspots < 0) = NaN;
   clear yyhalfweeklydata
   save(['MERRA2_MesoAmericanReef_HalfWeeklyHotSpots_' num2str(thisyear) '.mat'],'allhotspots');

   %%% Calculate Degree Heating Weeks using the following equation
   %%% DHW = 0.5 * summation of previous 24 twice-weekly hotspots > 1 degree

   alldhw = allhotspots*NaN;

   if (xyear>1),

     lastyear = thisyear-1;
     lastyyhotspots = importdata(['MERRA2_MesoAmericanReef_HalfWeeklyHotSpots_' num2str(lastyear) '.mat']);

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

   else      %%% if it is 1980, start at week 12

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
   save(['MERRA2_MesoAmericanReef_HalfWeeklyDHW_' num2str(thisyear) '.mat'],'alldhw');

   %%%% create regional time series of MERRA-2 half weekly DHW

   smallmerradhw = squeeze(alldhw(thislatidx,thislonidx,:));
   for xweek=1:size(smallmerradhw,3),
      smallmerradhw(:,:,xweek) = squeeze(smallmerradhw(:,:,xweek)).*smalloceanmask;
   end;
   areaweights = smallmerradhw*NaN;
   for xweek=1:size(smallmerradhw,3),
      subarea = squeeze(area(thislatidx,thislonidx));
      subarea(isnan(squeeze(smallmerradhw(:,:,xweek)))) = NaN;
      areaweights(:,:,xweek) = (subarea./nansum(nansum(subarea)));
   end;

   smallmerradhwweighted = smallmerradhw*NaN;
   for xweek=1:size(smallmerradhw,3),
      smallmerradhwweighted(:,:,xweek) = (smallmerradhw(:,:,xweek)).*(areaweights(:,:,xweek));
   end;

   %% full half-weekly time series
   merraregionaldhw = squeeze(nansum(squeeze(nansum(smallmerradhwweighted))));

   smallmerrahotspots = squeeze(allhotspots(thislatidx,thislonidx,:));
   for xweek=1:size(smallmerrahotspots,3),
      smallmerrahotspots(:,:,xweek) = squeeze(smallmerrahotspots(:,:,xweek)).*smalloceanmask;
   end;
   areaweights = smallmerrahotspots*NaN;
   for xweek=1:size(smallmerrahotspots,3),
      subarea = squeeze(area(thislatidx,thislonidx));
      subarea(isnan(squeeze(smallmerrahotspots(:,:,xweek)))) = NaN;
      areaweights(:,:,xweek) = (subarea./nansum(nansum(subarea)));
   end;

   smallmerrahotspotsweighted = smallmerrahotspots*NaN;
   for xweek=1:size(smallmerrahotspots,3),
      smallmerrahotspotsweighted(:,:,xweek) = (smallmerrahotspots(:,:,xweek)).*(areaweights(:,:,xweek));
   end;

   %% full half-weekly time series of hotspots
   merraregionalhotspots = squeeze(nansum(squeeze(nansum(smallmerrahotspotsweighted))));

   save(['MERRA2_RegionalDHW_HalfWeekly_' num2str(thisyear) '.mat'],'merraregionaldhw');
   save(['MERRA2_RegionalHotSpots_HalfWeekly_' num2str(thisyear) '.mat'],'merraregionalhotspots');

   clear allhotspots
   clear alldhw
   clear lastyyhotspots

end;

%% annual max time series of DHW
merraregionaldhwyy = ones(1,length(allyears))*NaN;

for xyear=1:length(allyears),
   thisyear = allyears(xyear);
   thisannualdhw = importdata(['MERRA2_RegionalDHW_HalfWeekly_' num2str(thisyear) '.mat']);
   merraregionaldhwyy(xyear) = nanmax(thisannualdhw);
end;

save(['MERRA2_RegionalDHW_AnnualMax.mat'],'merraregionaldhwyy');

%%% only 1980-2009 for baseline

clear smallmerradhw
clear smallmerradhwweighted
clear smallmerrahotspots
clear smallmerrahotspotsweighted
clear weeklymeananom

subyears = [1980:2009];
annualmaxdhw = ones(length(latidx),length(lonidx),length(subyears))*NaN;
for xyear=1:length(subyears),
   thisyear = subyears(xyear);
   thisyydhw = importdata(['MERRA2_MesoAmericanReef_HalfWeeklyDHW_' num2str(thisyear) '.mat']);
   annualmaxdhw(:,:,xyear) = nanmax(permute(thisyydhw,[3 1 2]));
   clear thisyydhw
end;
allmeanmaxdhw = squeeze(nanmean(annualmaxdhw,3));
clear annualmaxdhw

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
acr_pcolormapr4_nocoast(allmeanmaxdhw.*oceanmask,sublat,sublon,[slat nlat],[wlon elon],[-0.5 6.5],['MERRA-2 Annual Max Degree Heating Weeks 1980-2009'],['degree heating weeks (degrees C)']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',1.4);
print(f,'-dpng',['MERRA2_1980-2009_AvgMaxDHW.png']);
print(f,'-depsc',['MERRA2_1980-2009_AvgMaxDHW.eps']);
close all
save(['MERRA2_1980-2009_AvgMaxDHW.mat'],'allmeanmaxdhw');

clear allmeanmaxdhw

%%%%% maximum DHW during 1998 bleaching
%%%%% event (September 18 - October 1)

thisyydhw = importdata(['MERRA2_MesoAmericanReef_HalfWeeklyDHW_1998.mat']);

% maximum across 1998
allmaxdhw = squeeze(nanmax(squeeze(permute(thisyydhw,[3 1 2]))));
allmaxdhw(isnan(allmaxdhw)) = 0;

rainbow = importdata('/Users/mmphill2/Documents/MATLAB/cmap_rainbow_lowgrey.mat');
% to create these in R use
% colorbar = colorRampPalette(c('#000075','#4363d8','#3cb44b','#ffe119','#f58231','#800000'))(n = 17)
% n can be the length of the colorbar
hexcolors = importdata('dhwpaper_colorbar17.mat');
p1 = ones(length(hexcolors),3)*NaN;
for rgb=1:size(p1,1),
   p1(rgb,:) = hex2rgb(hexcolors(rgb));
end;
p = ones(length(hexcolors)+4,3)*NaN;
p(1:4,:) = rainbow(1:4,:);
p(5:end,:) = p1;

f = figure; colormap(p);
acr_pcolormapr4_nocoast(allmaxdhw.*oceanmask,sublat,sublon,[slat nlat],[wlon elon],[-0.5 20.5],['MERRA-2 Maximum Degree Heating Week During 1998'],['degree heating weeks (degrees C)']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['MERRA2_1998_MaxDHW.png']);
print(f,'-depsc',['MERRA2_1998_MaxDHW.eps']);
close all
save(['MERRA2_1998_MaxDHW.mat'],'allmaxdhw');

%% julian day 261-274
%% which is 74:78 in half weeks

bleaching1998 = squeeze(thisyydhw(:,:,74:78));
allmaxdhw = squeeze(nanmax(squeeze(permute(bleaching1998,[3 1 2]))));
allmaxdhw(isnan(allmaxdhw)) = 0;

f = figure; colormap(p);
acr_pcolormapr4_nocoast(allmaxdhw.*oceanmask,sublat,sublon,[slat nlat],[wlon elon],[-0.5 20.5],['MERRA-2 Maximum Degree Heating Week During 1998 Bleaching Event'],['degree heating weeks (degrees C)']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['MERRA2_1998BleachingEvent_MaxDHW.png']);
print(f,'-depsc',['MERRA2_1998BleachingEvent_MaxDHW.eps']);
close all
save(['MERRA2_1998BleachingEvent_MaxDHW.mat'],'allmaxdhw');

%%%%% maximum DHW during 2005 bleaching
%%%%% event (July 15 - November 15)

thisyydhw = importdata(['MERRA2_MesoAmericanReef_HalfWeeklyDHW_2005.mat']);

% maximum across 2005
allmaxdhw = squeeze(nanmax(squeeze(permute(thisyydhw,[3 1 2]))));
allmaxdhw(isnan(allmaxdhw)) = 0;

f = figure; colormap(p);
acr_pcolormapr4_nocoast(allmaxdhw.*oceanmask,sublat,sublon,[slat nlat],[wlon elon],[-0.5 20.5],['MERRA-2 Maximum Degree Heating Week During 2005'],['degree heating weeks (degrees C)']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['MERRA2_2005_MaxDHW.png']);
print(f,'-depsc',['MERRA2_2005_MaxDHW.eps']);
close all
save(['MERRA2_2005_MaxDHW.mat'],'allmaxdhw');

%% julian day 196-319
%% which is 56:91 in half weeks

bleaching2005 = squeeze(thisyydhw(:,:,56:91));
allmaxdhw = squeeze(nanmax(squeeze(permute(bleaching2005,[3 1 2]))));
allmaxdhw(isnan(allmaxdhw)) = 0;

f = figure; colormap(p);
acr_pcolormapr4_nocoast(allmaxdhw.*oceanmask,sublat,sublon,[slat nlat],[wlon elon],[-0.5 20.5],['MERRA-2 Maximum Degree Heating Week During 2005 Bleaching Event'],['degree heating weeks (degrees C)']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['MERRA2_2005BleachingEvent_MaxDHW.png']);
print(f,'-depsc',['MERRA2_2005BleachingEvent_MaxDHW.eps']);
close all
save(['MERRA2_2005BleachingEvent_MaxDHW.mat'],'allmaxdhw');

clear bleaching1998
clear bleaching2005
clear allmaxdhw
clear thisyydhw

%%%%% 2010-2019 versions
%%%%% to compare avg annual max and frequency maps

decadeyears = [2010:2019];
annualmaxdhw = ones(length(latidx),length(lonidx),length(decadeyears))*NaN;
for xyear=1:length(decadeyears),
   thisyear = decadeyears(xyear);
   thisyydhw = importdata(['MERRA2_MesoAmericanReef_HalfWeeklyDHW_' num2str(thisyear) '.mat']);
   annualmaxdhw(:,:,xyear) = nanmax(permute(thisyydhw,[3 1 2]));
   clear thisyydhw
end;
allmeanmaxdhw = squeeze(nanmean(annualmaxdhw,3));
clear annualmaxdhw

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
acr_pcolormapr4_nocoast(allmeanmaxdhw.*oceanmask,sublat,sublon,[slat nlat],[wlon elon],[-0.5 6.5],['MERRA-2 Annual Max Degree Heating Weeks 2010-2019'],['degree heating weeks (degrees C)']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',1.4);
print(f,'-dpng',['MERRA2_2010-2019_AvgMaxDHW.png']);
print(f,'-depsc',['MERRA2_2010-2019_AvgMaxDHW.eps']);
close all
save(['MERRA2_2010-2019_AvgMaxDHW.mat'],'allmeanmaxdhw');

hexcolors = importdata('dhwpaper_colorbar20.mat');
p = ones(length(hexcolors),3)*NaN;
for rgb=1:size(p,1),
   p(rgb,:) = hex2rgb(hexcolors(rgb));
end;

f = figure; colormap(p);
acr_pcolormapr4_nocoast((maxmonth.*oceanmask)-273.15,sublat,sublon,[slat nlat],[wlon elon],[27.5 30],['MERRA-2 1985-1993 Climatological Max Mean Month SST'],['degrees C']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',2);
print(f,'-dpng',['MERRA2_1985-1993_MaxMonth.png']);
print(f,'-depsc',['MERRA2_1985-1993_MaxMonth.eps']);
close all

clear allmeanmaxdhw

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% FREQUENCY of hotspots above 1C

merradhwfreq = ones(length(latidx),length(lonidx),length(subyears))*NaN;
merradhwduration = ones(length(latidx),length(lonidx),length(subyears))*NaN;

for xyear=1:length(subyears),
   thisyear = subyears(xyear);
   thisyyhotspots = importdata(['MERRA2_MesoAmericanReef_HalfWeeklyHotSpots_' num2str(thisyear) '.mat']);
   for ii=1:length(lonidx),
      for jj=1:length(latidx),
         merrasstpoint = squeeze(thisyyhotspots(jj,ii,:));
         merradhwfreq(jj,ii,xyear) = length(find(merrasstpoint>1));

         clear streak
	       if (merrasstpoint(1) > 1),
	          streak(1) = 1;
	       else
	          streak(1) = 0;
	       end;
	       for dd=2:length(merrasstpoint),
	          if (merrasstpoint(dd) > 1),
	             streak(dd) = streak(dd-1) + 1;
	             streak(dd-1) = 0;
	          else
	             streak(dd) = 0;
	          end;
	        end;
	        merradhwduration(jj,ii,xyear) = nanmax(streak);
      end;
   end;
end;

merradhwduration = squeeze(nanmean(merradhwduration,3));
merradhwfreq = squeeze(nanmean(merradhwfreq,3));

% for this one in R use
% colorbar = brewer.pal(n = 8, name = "PuRd")
% writeMat("dhwpaper_purplered8.mat",colorbar=colorbar)
%hexcolors = importdata('dhwpaper_purplered8.mat');
%p = ones(length(hexcolors),3)*NaN;
%for rgb=1:size(p,1),
%   p(rgb,:) = hex2rgb(hexcolors(rgb));
%end;

% for this one in R use
% colorbar = colorRampPalette(c('#aaffc3','#42d4f4','#f032e6','#911eb4','#e6194B'))(n = 7)
% writeMat("dhwpaper_tealmagenta7.mat",colorbar=colorbar)
hexcolors = importdata('dhwpaper_tealmagenta7.mat');
p = ones(length(hexcolors),3)*NaN;
for rgb=1:size(p,1),
   p(rgb,:) = hex2rgb(hexcolors(rgb));
end;

f = figure; colormap(p);
% divide by 2 for weeks instead of half-weeks
acr_pcolormapr4_nocoast((merradhwfreq/2).*oceanmask,sublat,sublon,[slat nlat],[wlon elon],[-0.5 3],['MERRA-2 Mean Annual Frequency of Hotspots>1C 1980-2009'],['number of weeks with >1C anomaly from baseline']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',1.4);
print(f,'-dpng',['MERRA2_HotSpots_Frequency_1980-2009.png']);
print(f,'-depsc',['MERRA2_HotSpots_Frequency_1980-2009.eps']);
close all
save(['MERRA2_HotSpots_Frequency_1980-2009.mat'],'merradhwfreq');

%hexcolors = importdata('dhwpaper_purplered5.mat');
%p = ones(length(hexcolors),3)*NaN;
%for rgb=1:size(p,1),
%   p(rgb,:) = hex2rgb(hexcolors(rgb));
%end;

hexcolors = importdata('dhwpaper_tealmagenta5.mat');
p = ones(length(hexcolors),3)*NaN;
for rgb=1:size(p,1),
   p(rgb,:) = hex2rgb(hexcolors(rgb));
end;

f = figure; colormap(p);
acr_pcolormapr4_nocoast((merradhwduration/2).*oceanmask,sublat,sublon,[slat nlat],[wlon elon],[-0.5 2],['MERRA-2 Mean Annual Max Duration of Hotspots>1C 1980-2009'],['maximum consecutive weeks with >1C anomaly from baseline']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',1.4);
print(f,'-dpng',['MERRA2_HotSpots_MaxDuration_1980-2009.png']);
print(f,'-depsc',['MERRA2_HotSpots_MaxDuration_1980-2009.eps']);
close all
save(['MERRA2_HotSpots_MaxDuration_1980-2009.mat'],'merradhwduration');

%%%% 2010 to 2019 version

merradhwfreq = ones(length(latidx),length(lonidx),length(decadeyears))*NaN;
merradhwduration = ones(length(latidx),length(lonidx),length(decadeyears))*NaN;

for xyear=1:length(decadeyears),
   thisyear = decadeyears(xyear);
   thisyyhotspots = importdata(['MERRA2_MesoAmericanReef_HalfWeeklyHotSpots_' num2str(thisyear) '.mat']);
   for ii=1:length(lonidx),
      for jj=1:length(latidx),
         merrasstpoint = squeeze(thisyyhotspots(jj,ii,:));
         merradhwfreq(jj,ii,xyear) = length(find(merrasstpoint>1));

         clear streak
	       if (merrasstpoint(1) > 1),
	          streak(1) = 1;
	       else
	          streak(1) = 0;
	       end;
	       for dd=2:length(merrasstpoint),
	          if (merrasstpoint(dd) > 1),
	             streak(dd) = streak(dd-1) + 1;
	             streak(dd-1) = 0;
	          else
	             streak(dd) = 0;
	          end;
	        end;
	        merradhwduration(jj,ii,xyear) = nanmax(streak);
      end;
   end;
end;

merradhwduration = squeeze(nanmean(merradhwduration,3));
merradhwfreq = squeeze(nanmean(merradhwfreq,3));

%hexcolors = importdata('dhwpaper_purplered8.mat');
%p = ones(length(hexcolors),3)*NaN;
%for rgb=1:size(p,1),
%   p(rgb,:) = hex2rgb(hexcolors(rgb));
%end;

hexcolors = importdata('dhwpaper_tealmagenta7.mat');
p = ones(length(hexcolors),3)*NaN;
for rgb=1:size(p,1),
   p(rgb,:) = hex2rgb(hexcolors(rgb));
end;

f = figure; colormap(p);
acr_pcolormapr4_nocoast((merradhwfreq/2).*oceanmask,sublat,sublon,[slat nlat],[wlon elon],[-0.5 3],['MERRA-2 Mean Annual Frequency of Hotspots>1C 2010-2019'],['number of weeks with >1C anomaly from baseline']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',1.4);
print(f,'-dpng',['MERRA2_HotSpots_Frequency_2010-2019.png']);
print(f,'-depsc',['MERRA2_HotSpots_Frequency_2010-2019.eps']);
close all
save(['MERRA2_HotSpots_Frequency_2010-2019.mat'],'merradhwfreq');

%hexcolors = importdata('dhwpaper_purplered5.mat');
%p = ones(length(hexcolors),3)*NaN;
%for rgb=1:size(p,1),
%   p(rgb,:) = hex2rgb(hexcolors(rgb));
%end;

hexcolors = importdata('dhwpaper_tealmagenta5.mat');
p = ones(length(hexcolors),3)*NaN;
for rgb=1:size(p,1),
   p(rgb,:) = hex2rgb(hexcolors(rgb));
end;

f = figure; colormap(p);
acr_pcolormapr4_nocoast((merradhwduration/2).*oceanmask,sublat,sublon,[slat nlat],[wlon elon],[-0.5 2],['MERRA-2 Mean Annual Max Duration of Hotspots>1C 2010-2019'],['maximum consecutive weeks with >1C anomaly from baseline']);
hold on;
h1 = borders('countries','color','k','linewidth',1);
hold on;
plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',1.4);
print(f,'-dpng',['MERRA2_HotSpots_MaxDuration_2010-2019.png']);
print(f,'-depsc',['MERRA2_HotSpots_MaxDuration_2010-2019.eps']);
close all
save(['MERRA2_HotSpots_MaxDuration_2010-2019.mat'],'merradhwduration');

clear thisyyhotspots
clear merradhwfreq
clear merradhwduration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% climatological monthly mean
%%%%%%%%%%%%%%%%%%% 1980-2009

mmsstclim = ones(length(latidx),length(lonidx),12,length(subyears))*NaN;

for xyear=1:length(subyears),
   thisyear = subyears(xyear);
   if ((thisyear == 1980)||(thisyear == 1984)||(thisyear == 1988)||(thisyear == 1992)||(thisyear == 1996)||(thisyear == 2000)||(thisyear == 2004)||(thisyear == 2008)||(thisyear == 2012)||(thisyear == 2016)),
      mmstart = [1 32 61 92 122 153 183 214 245 275 306 336];
      mmlength = [31 29 31 30 31 30 31 31 30 31 30 31];
      yearlength = 366;
   else
      mmstart = [1 32 60 91 121 152 182 213 244 274 305 335];
      mmlength = [31 28 31 30 31 30 31 31 30 31 30 31];
      yearlength = 365;
   end;

   filename = ['merra2_sst_MAR_HR_daily_' num2str(thisyear) '.mat'];
   subdata = importdata([filename]);
   for xmonth=1:12,
      mmsstclim(:,:,xmonth,xyear) = squeeze(nanmean(squeeze(subdata(:,:,(mmstart(xmonth):((mmstart(xmonth)+mmlength(xmonth))-1)))),3));
   end;
   clear subdata
end;

months = {'Jan','Feb','March','April','May','June','July','Aug','Sept','Oct','Nov','Dec'};

z = ones(10,3)*NaN;
z(1,:) = [0 0 0.6250];
z(2,:) = [0 0.1250 1];
z(3,:) = [0 0.5625 1];
z(4,:) = [0 0.9375 1];
z(5,:) = [0.375 1 0.625];
z(6,:) = [0.8125 1 0.1875];
z(7,:) = [1 0.75 0];
z(8,:) = [1 0.3125 0];
z(9,:) = [0.875 0 0];
z(10,:) = [0.5 0 0];

for mm=1:12,
   thismonth = months{mm};
   map = squeeze(nanmean(squeeze(mmsstclim(:,:,mm,:)),3));
   f = figure; colormap(z);
   acr_pcolormapr4_nocoast((map-273.15).*oceanmask,sublat,sublon,[slat nlat],[wlon elon],[25 30],['MERRA-2 Mean ' thismonth ' SSTs 1980-2009'],['degrees C']);
   hold on;
   h1 = borders('countries','color','k','linewidth',1);
   hold on;
   plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',1.4);
   print(f,'-dpng',['MERRA2_MeanSST_' num2str(mm) thismonth '_1980-2009.png']);
   print(f,'-depsc',['MERRA2_MeanSST_' num2str(mm) thismonth '_1980-2009.eps']);
   close all
   save(['MERRA2_MeanSST_' num2str(mm) thismonth '_1980-2009.mat'],'map');
end;

%%%% create regional monthly climatological SST variable
smallmmsstclim = squeeze(nanmean(squeeze(mmsstclim(thislatidx,thislonidx,:,:)),4));
for xweek=1:size(smallmmsstclim,3),
   smallmmsstclim(:,:,xweek) = squeeze(smallmmsstclim(:,:,xweek)).*smalloceanmask;
end;

areaweights = smallmmsstclim*NaN;
for xweek=1:size(smallmmsstclim,3),
   subarea = squeeze(area(thislatidx,thislonidx));
   subarea(isnan(squeeze(smallmmsstclim(:,:,xweek)))) = NaN;
   areaweights(:,:,xweek) = (subarea./nansum(nansum(subarea)));
end;

smallmmsstclimweighted = smallmmsstclim*NaN;
for xweek=1:size(smallmmsstclim,3),
   smallmmsstclimweighted(:,:,xweek) = (smallmmsstclim(:,:,xweek)).*(areaweights(:,:,xweek));
end;

%% full monthly
merraregionalsstclim = squeeze(nansum(squeeze(nansum(smallmmsstclimweighted))));

save(['MERRA2_RegionalClimatologicalSST_Monthly.mat'],'merraregionalsstclim');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% SST maps of bleaching event in 1998

mmstart = [1 32 60 91 121 152 182 213 244 274 305 335];
mmlength = [31 28 31 30 31 30 31 31 30 31 30 31];

thisyydata = importdata(['merra2_sst_MAR_HR_daily_1998.mat']);

z = ones(10,3)*NaN;
z(1,:) = [0 0 0.6250];
z(2,:) = [0 0.1250 1];
z(3,:) = [0 0.5625 1];
z(4,:) = [0 0.9375 1];
z(5,:) = [0.375 1 0.625];
z(6,:) = [0.8125 1 0.1875];
z(7,:) = [1 0.75 0];
z(8,:) = [1 0.3125 0];
z(9,:) = [0.875 0 0];
z(10,:) = [0.5 0 0];

map = squeeze(thisyydata(:,:,261:274));
clear thisyydata
for dd=1:size(map,3),
   f = figure; colormap(z);
   thismap = squeeze(map(:,:,dd));
   acr_pcolormapr4_nocoast((thismap-273.15).*oceanmask,sublat,sublon,[slat nlat],[wlon elon],[27.5 32.5],['MERRA-2 1998 Bleaching Event Day ' num2str(dd)],['degrees C']);
   hold on;
   h1 = borders('countries','color','k','linewidth',1);
   hold on;
   plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',1.4);
   print(f,'-dpng',['MERRA2_SSTs_1998bleaching_day' num2str(dd) '.png']);
   print(f,'-depsc',['MERRA2_SSTs_1998bleaching_day' num2str(dd) '.eps']);
   close all
   save(['MERRA2_SSTs_1998bleaching_day' num2str(dd) '.mat'],'thismap');
end;
