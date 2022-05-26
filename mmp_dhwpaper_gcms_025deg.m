%	               mmp_dhwpaper_gcms_025deg
%
%  This script calculates degree heating weeks and other
%  corresponding metrics from various CMIP5 and CMIP6 model
%  sea surface temperature outputs for the WWF Mesoamerican
%  Reef project and corresponding paper with Model-E.
%  This version uses 0.25-degree resolution SST data.
%
%  -- Creates a map of each GCM’s climatological mean SST
%  for each month in the 1980-2009 period.
%
%  -- Creates a map of each GCM’s mean annual maximum DHW
%  value for the baseline 1980-2009 period.
%
%  -- Creates a map of each GCM’s mean annual maximum DHW
%  value for the 2050-2059 period under each rcp45/ssp245 and
%  rcp85/ssp585.
%
%  -- Creates a map of each GCM’s change between baseline and
%  future mean annual maximum DHW value for the 2050-2059 period
%  under each rcp45/ssp245 and rcp85/ssp585.
%
%  -- Creates a map of each GCM’s mean annual frequency of
%  HotSpots greater than 1C for the 2050-2059 period under each
%  rcp45/ssp245 and rcp85/ssp585.
%
%  -- Creates a map of each GCM’s change between baseline and
%  future mean annual frequency of HotSpots greater than 1C for
%  the 2050-2059 period under each rcp45/ssp245 and rcp85/ssp585.
%
%  -- Creates a map of each GCM’s mean annual maximum consecutive
%  duration of HotSpots greater than 1C for the 2050-2059 period
%  under each rcp45/ssp245 and rcp85/ssp585.
%
%  -- Creates a map of each GCM’s change between baseline and
%  future mean annual maximum consecutive duration of HotSpots
%  greater than 1C for the 2050-2059 period under each rcp45/ssp245
%  and rcp85/ssp585.
%
%  INPUTS: none
%  OUTPUTS: none
%
%
%
%		                             author: Meridel Phillips
%                                  mmp2192@columbia.edu
%				                       date: 9/22/2020
%
function mmp_dhwpaper_gcms_025deg();
%--------------------------------------------------
%--------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Begin Debug
%% End Debug

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /Users/mmphill2/Documents/MATLAB/
cd /Users/mmphill2/datasets/Andromeda/temp/

institution = {'ACCESS1-0','GISS-E2-R','NorESM1-M','ACCESS-CM2','CanESM5','GFDL-CM4','GISS_E213','IPSL-CM6A-LR','MRI-ESM2-0','NorESM2-LM'};
ensemble = {'r1i1p1','r6i1p1','r1i1p1','r1i1p1f1','r1i1p1f1','r1i1p1f1','none','r1i1p1f1','r1i1p1f1','r1i1p1f1'};

for xgcm = 1:length(institution),

   thisgcm = institution{xgcm};
   if (xgcm<4),
      histyears = [1980:2005];
      allscen = {'rcp45','rcp85'};
      longscen = {'RCP45','RCP85'};
   else
      histyears = [1980:2014];
      allscen = {'ssp245','ssp585'};
      longscen = {'SSP245','SSP585'};
   end;

   if (strcmp(thisgcm,'GISS_E213')),
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
       regionslat = 14.21;
       regionnlat = 23.55;
       regionwlon = -91.20;
       regionelon = -79.54;

   else

       gcmtitle = thisgcm;

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
       regionslat = 14.21;
       regionnlat = 23.55;
       regionwlon = (-91.20)+360;
       regionelon = (-79.54)+360;

   end;

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

   months = {'Jan','Feb','March','April','May','June','July','Aug','Sept','Oct','Nov','Dec'};
   allmm = {'01','02','03','04','05','06','07','08','09','10','11','12'};

   cmipmask = importdata('OceanMask_025degree.mat');
   smallcmipmask = squeeze(cmipmask(thislatidx,thislonidx));

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%% climatological monthly mean
   %%%%%%%%%%%%%%%%%%% 1980-2009

   allyears = [1980:2009];
   leapyears = [1976:4:2099];
   mmsstclim = ones(length(latidx),length(lonidx),12,length(allyears))*NaN;

   for xyear=1:length(allyears),
      thisyear = allyears(xyear);
      if ((strcmp(thisgcm,'ACCESS1-0'))||(strcmp(thisgcm,'CMCC-CMS'))||(strcmp(thisgcm,'MIROC-ESM'))||(strcmp(thisgcm,'MPI-ESM-LR'))||(strcmp(thisgcm,'ACCESS-CM2'))||(strcmp(thisgcm,'IPSL-CM6A-LR'))||(strcmp(thisgcm,'MIROC6'))||(strcmp(thisgcm,'MPI-ESM1-2-LR'))||(strcmp(thisgcm,'MRI-ESM2-0'))),
         if (ismember(thisyear,leapyears)),
            mmstart = [1 32 61 92 122 153 183 214 245 275 306 336];
            mmlength = [31 29 31 30 31 30 31 31 30 31 30 31];
            yearlength = 366;
         else
            mmstart = [1 32 60 91 121 152 182 213 244 274 305 335];
            mmlength = [31 28 31 30 31 30 31 31 30 31 30 31];
            yearlength = 365;
         end;
      else
         mmstart = [1 32 60 91 121 152 182 213 244 274 305 335];
         mmlength = [31 28 31 30 31 30 31 31 30 31 30 31];
         yearlength = 365;
      end;
      if (xgcm<4),
         if (xyear<27),
            filename = ['tos_MAR_HR_day_CMIP5_' thisgcm '_historical_' num2str(thisyear) '.mat'];
         else
            filename = ['tos_MAR_HR_day_CMIP5_' thisgcm '_rcp85_' num2str(thisyear) '.mat'];
         end;
      else
         filename = ['tos_MAR_HR_day_CMIP6_' thisgcm '_historical_' num2str(thisyear) '.mat'];
      end;
      subdata = importdata([filename]);
      for xmonth=1:12,
         mmsstclim(:,:,xmonth,xyear) = squeeze(nanmean(squeeze(subdata(:,:,(mmstart(xmonth):((mmstart(xmonth)+mmlength(xmonth))-1)))),3));
      end;
      clear subdata
   end;

   months = {'Jan','Feb','March','April','May','June','July','Aug','Sept','Oct','Nov','Dec'};

   save([thisgcm '_MesoamericanReef_AllMonthlySSTs_1980-2009.mat'],'mmsstclim');

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
      if (xgcm<4),
         acr_pcolormapr4_nocoast((map-273.15).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[25 30],[gcmtitle ' CMIP5 Mean ' thismonth ' SSTs 1980-2009'],['degrees C']);
      else
         acr_pcolormapr4_nocoast((map).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[25 30],[gcmtitle ' CMIP6 Mean ' thismonth ' SSTs 1980-2009'],['degrees C']);
      end;
      hold on;
      h1 = borders('countries','color','k','linewidth',1);
      hold on;
      plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',1.4);
      print(f,'-dpng',[thisgcm '_MeanSST_' num2str(mm) thismonth '_1980-2009.png']);
      print(f,'-depsc',[thisgcm '_MeanSST_' num2str(mm) thismonth '_1980-2009.eps']);
      close all
      save([thisgcm '_MeanSST_' num2str(mm) thismonth '_1980-2009.mat'],'map');
   end;

   %%%% create regional monthly climatological SST variable
   smallmmsstclim = squeeze(nanmean(squeeze(mmsstclim(thislatidx,thislonidx,:,:)),4));
   for xweek=1:size(smallmmsstclim,3),
      smallmmsstclim(:,:,xweek) = squeeze(smallmmsstclim(:,:,xweek)).*smallcmipmask;
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
   regionalsstclim = squeeze(nansum(squeeze(nansum(smallmmsstclimweighted))));

   save([thisgcm '_RegionalClimatologicalSST_Monthly.mat'],'regionalsstclim');

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%% DHW maps for baseline

   %%% for DHW calculation:
   %%% find monthly mean SST climatologies from 1985 to 1993,
   %%% excluding 1991 and 1992 (lat x lon x 12 months x 7 years)

   %% index of years 1985-1990 and 1993
   climyears = [6 7 8 9 10 11 14];
   submmsstclim = squeeze(mmsstclim(:,:,:,climyears));
   clear mmsstclim

   %%% Next find montly mean SST climatology over
   %%% entire 7-year period (lat x lon x 12 months)

   meanmmsstclim = squeeze(nanmean(submmsstclim,4));
   clear submmsstclim

   %%% Calculate the maximum monthly mean SST climatology
   %%% over entire 7-year period (lat x lon)

   maxmonth = squeeze(nanmax(squeeze(permute(meanmmsstclim,[3 1 2]))));
   clear meanmmsstclim
   save([thisgcm '_MesoAmericanReef_MaxMeanMonth_1985-1993.mat'],'maxmonth');

   maxmonth = importdata([thisgcm '_MesoAmericanReef_MaxMeanMonth_1985-1993.mat']);

   clear map
   clear smallmmsstclim
   clear smallmmsstclimweighted

   %%% Calculate twice-weekly SST from 1980 to 2009
   %%% not using pre-1980 so first 12 weeks of 1980 will not calculate DHW

   allyears = [1980:2009];
   leapyears = [1980:4:2009];
   nhalfweeks = round((((length(allyears)*365)+length(leapyears))/7)*2);

   for xyear=1:length(allyears),
      thisyear = allyears(xyear);

      if ((strcmp(thisgcm,'ACCESS1-0'))||(strcmp(thisgcm,'CMCC-CMS'))||(strcmp(thisgcm,'MIROC-ESM'))||(strcmp(thisgcm,'MPI-ESM-LR'))||(strcmp(thisgcm,'ACCESS-CM2'))||(strcmp(thisgcm,'IPSL-CM6A-LR'))||(strcmp(thisgcm,'MIROC6'))||(strcmp(thisgcm,'MPI-ESM1-2-LR'))||(strcmp(thisgcm,'MRI-ESM2-0'))),
         if (ismember(thisyear,leapyears)),
            mmstart = [1 32 61 92 122 153 183 214 245 275 306 336];
            mmlength = [31 29 31 30 31 30 31 31 30 31 30 31];
            yearlength = 366;
         else
            mmstart = [1 32 60 91 121 152 182 213 244 274 305 335];
            mmlength = [31 28 31 30 31 30 31 31 30 31 30 31];
            yearlength = 365;
         end;
      else
         mmstart = [1 32 60 91 121 152 182 213 244 274 305 335];
         mmlength = [31 28 31 30 31 30 31 31 30 31 30 31];
         yearlength = 365;
      end;

      if (xgcm<4),
         if (xyear<27),
            filename = ['tos_MAR_HR_day_CMIP5_' thisgcm '_historical_' num2str(thisyear) '.mat'];
         else
            filename = ['tos_MAR_HR_day_CMIP5_' thisgcm '_rcp85_' num2str(thisyear) '.mat'];
         end;
      else
         filename = ['tos_MAR_HR_day_CMIP6_' thisgcm '_historical_' num2str(thisyear) '.mat'];
      end;
      subdata = importdata([filename]);
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
      save([thisgcm '_MesoAmericanReef_HalfWeeklySSTs_' num2str(thisyear) '.mat'],'yyhalfweeklydata');

      %%% Calculate hotspots by subtracting max monthly SST
      %%% climatology from each twice-weekly SST

      allhotspots = yyhalfweeklydata*NaN;
      for xweek=1:size(yyhalfweeklydata,3),
         allhotspots(:,:,xweek) = yyhalfweeklydata(:,:,xweek)-maxmonth;
      end;
      allhotspots(allhotspots < 0) = NaN;
      clear yyhalfweeklydata
      save([thisgcm '_MesoAmericanReef_HalfWeeklyHotSpots_' num2str(thisyear) '.mat'],'allhotspots');

      %%% Calculate Degree Heating Weeks using the following equation
      %%% DHW = 0.5 * summation of previous 24 twice-weekly hotspots > 1 degree

      alldhw = allhotspots*NaN;

      if (xyear>1),

        lastyear = thisyear-1;
        lastyyhotspots = importdata([thisgcm '_MesoAmericanReef_HalfWeeklyHotSpots_' num2str(lastyear) '.mat']);

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
      save([thisgcm '_MesoAmericanReef_HalfWeeklyDHW_' num2str(thisyear) '.mat'],'alldhw');

      %%%% create regional time series of half weekly DHW

      smalldhw = squeeze(alldhw(thislatidx,thislonidx,:));
      for xweek=1:size(smalldhw,3),
         smalldhw(:,:,xweek) = squeeze(smalldhw(:,:,xweek)).*smallcmipmask;
      end;
      areaweights = smalldhw*NaN;
      for xweek=1:size(smalldhw,3),
         subarea = squeeze(area(thislatidx,thislonidx));
         subarea(isnan(squeeze(smalldhw(:,:,xweek)))) = NaN;
         areaweights(:,:,xweek) = (subarea./nansum(nansum(subarea)));
      end;

      smalldhwweighted = smalldhw*NaN;
      for xweek=1:size(smalldhw,3),
         smalldhwweighted(:,:,xweek) = (smalldhw(:,:,xweek)).*(areaweights(:,:,xweek));
      end;

      %% full half-weekly time series
      regionaldhw = squeeze(nansum(squeeze(nansum(smalldhwweighted))));

      smallhotspots = squeeze(allhotspots(thislatidx,thislonidx,:));
      for xweek=1:size(smallhotspots,3),
         smallhotspots(:,:,xweek) = squeeze(smallhotspots(:,:,xweek)).*smallcmipmask;
      end;
      areaweights = smallhotspots*NaN;
      for xweek=1:size(smallhotspots,3),
         subarea = squeeze(area(thislatidx,thislonidx));
         subarea(isnan(squeeze(smallhotspots(:,:,xweek)))) = NaN;
         areaweights(:,:,xweek) = (subarea./nansum(nansum(subarea)));
      end;

      smallhotspotsweighted = smallhotspots*NaN;
      for xweek=1:size(smallhotspots,3),
         smallhotspotsweighted(:,:,xweek) = (smallhotspots(:,:,xweek)).*(areaweights(:,:,xweek));
      end;

      %% full half-weekly time series
      regionalhotspots = squeeze(nansum(squeeze(nansum(smallhotspotsweighted))));

      save([thisgcm '_RegionalDHW_HalfWeekly_' num2str(thisyear) '.mat'],'regionaldhw');
      save([thisgcm '_RegionalHotSpots_HalfWeekly_' num2str(thisyear) '.mat'],'regionalhotspots');

      clear allhotspots
      clear alldhw
      clear lastyyhotspots

   end;

   %% annual max time series of DHW
   regionaldhwyy = ones(1,length(allyears))*NaN;

   for xyear=1:length(allyears),
      thisyear = allyears(xyear);
      thisannualdhw = importdata([thisgcm '_RegionalDHW_HalfWeekly_' num2str(thisyear) '.mat']);
      regionaldhwyy(xyear) = nanmax(thisannualdhw);
   end;

   save([thisgcm '_RegionalDHW_AnnualMax.mat'],'regionaldhwyy');

   clear smalldhw
   clear smalldhwweighted
   clear smallhotspots
   clear smallhotspotsweighted
   clear weeklymeananom

   annualmaxdhw = ones(length(latidx),length(lonidx),length(allyears))*NaN;
   for xyear=1:length(allyears),
      thisyear = allyears(xyear);
      thisyydhw = importdata([thisgcm '_MesoAmericanReef_HalfWeeklyDHW_' num2str(thisyear) '.mat']);
      annualmaxdhw(:,:,xyear) = nanmax(permute(thisyydhw,[3 1 2]));
      clear thisyydhw
   end;
   allmeanmaxdhw = squeeze(nanmean(annualmaxdhw,3));
   clear annualmaxdhw
   allmeanmaxdhw(isnan(allmeanmaxdhw)) = 0;

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
   acr_pcolormapr4_nocoast(allmeanmaxdhw.*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[-0.5 6.5],[gcmtitle ' Annual Max Degree Heating Weeks 1980-2009'],['degree heating weeks (degrees C)']);
   hold on;
   h1 = borders('countries','color','k','linewidth',1);
   hold on;
   plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',1.4);
   print(f,'-dpng',[thisgcm '_1980-2009_AvgMaxDHW.png']);
   print(f,'-depsc',[thisgcm '_1980-2009_AvgMaxDHW.eps']);
   close all
   save([thisgcm '_1980-2009_AvgMaxDHW.mat'],'allmeanmaxdhw');

   clear allmeanmaxdhw

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%    BIAS ADJUSTMENT

   histyears = [1980:2009];
   allyears = [2010:2100];
   slices = {'2010s','2020s','2030s','2040s','2050s','2060s','2070s','2080s','2090s'};
   slicestart = [2010 2020 2030 2040 2050 2060 2070 2080 2090];
   sliceend = [2019 2029 2039 2049 2059 2069 2079 2089 2099];
   idxstart = [1 11 21 31 41 51 61 71 81];
   idxend = [10 20 30 40 50 60 70 80 90];

   gcmhistorical = load([thisgcm '_MesoamericanReef_AllMonthlySSTs_1980-2009.mat']);
   gcmhistorical = gcmhistorical.mmsstclim;
   subgcmhistorical = squeeze(gcmhistorical(thislatidx,thislonidx,:,:));

   meanhist = squeeze(nanmean(gcmhistorical,4));
   smallmeanhist = squeeze(nanmean(subgcmhistorical,4));
   clear gcmhistorical
   clear subgcmhistorical

   for xscen=1:length(allscen),
      thisscen = allscen{xscen};
      thisscennice = longscen{xscen};

      for tt=1:length(slices),
         thisslice = slices{tt};
         sliceyrs = [slicestart(tt):sliceend(tt)];
         idxyears = [idxstart(tt):idxend(tt)];

         gcmfuture = ones(length(latidx),length(lonidx),12,length(sliceyrs))*NaN;

         for xyear=1:length(sliceyrs),
            thisyear = sliceyrs(xyear);
            if (xgcm<4),
               filename = ['tos_MAR_HR_day_CMIP5_' thisgcm '_' thisscen '_' num2str(thisyear) '.mat'];
            else
               if (thisyear<2015),
                  filename = ['tos_MAR_HR_day_CMIP6_' thisgcm '_historical_' num2str(thisyear) '.mat'];
               else
                  filename = ['tos_MAR_HR_day_CMIP6_' thisgcm '_' thisscen '_' num2str(thisyear) '.mat'];
               end;
            end;
            subdata = importdata([filename]);
   	        if (size(subdata,3)>365),
               mmstart = [1 32 61 92 122 153 183 214 245 275 306 336];
               mmlength = [31 29 31 30 31 30 31 31 30 31 30 31];
            else
               mmstart = [1 32 60 91 121 152 182 213 244 274 305 335];
               mmlength = [31 28 31 30 31 30 31 31 30 31 30 31];
            end;
            for xmonth=1:12,
               gcmfuture(:,:,xmonth,xyear) = squeeze(nanmean(squeeze(subdata(:,:,(mmstart(xmonth):((mmstart(xmonth)+mmlength(xmonth))-1)))),3));
            end;
            clear subdata
         end;
         subgcmfuture = squeeze(gcmfuture(thislatidx,thislonidx,:,:));
         meanfut = squeeze(nanmean(gcmfuture,4));
         clear gcmfuture
         smallmeanfut = squeeze(nanmean(subgcmfuture,4));
         clear subgcmfuture

         gcmdelta = meanfut*NaN;
         smallgcmdelta = smallmeanfut*NaN;
         for xmonth=1:12,
            gcmdelta(:,:,xmonth) = (squeeze(meanfut(:,:,xmonth)))-(squeeze(meanhist(:,:,xmonth)));
   	        smallgcmdelta(:,:,xmonth) = (squeeze(smallmeanfut(:,:,xmonth)))-(squeeze(smallmeanhist(:,:,xmonth)));
         end;

         %%%% area average gcm delta values for line plot
         for xweek=1:size(smallgcmdelta,3),
            smallgcmdelta(:,:,xweek) = squeeze(smallgcmdelta(:,:,xweek)).*smallcmipmask;
         end;
         areaweights = smallgcmdelta*NaN;
         for xmonth=1:size(smallgcmdelta,3),
            subarea = squeeze(area(thislatidx,thislonidx));
            subarea(isnan(squeeze(smallgcmdelta(:,:,xmonth)))) = NaN;
            areaweights(:,:,xmonth) = (subarea./nansum(nansum(subarea)));
         end;

         regionaldeltas = ones(1,length(months))*NaN;
         gcmdeltaweighted = smallgcmdelta*NaN;
         for xmonth=1:size(smallgcmdelta,3),
            gcmdeltaweighted(:,:,xmonth) = (smallgcmdelta(:,:,xmonth)).*(areaweights(:,:,xmonth));
         end;
         regionaldeltas(1,:) = squeeze(nansum(squeeze(nansum(gcmdeltaweighted))));
         save([thisgcm '_' thisscen '_' thisslice '_RegionalDeltas.mat'],'regionaldeltas');

         for xbaseyear=1:length(histyears),
            thisyear = histyears(xbaseyear);
            merradaily = importdata(['merra2_sst_MAR_HR_daily_' num2str(thisyear) '.mat']);
            if (size(merradaily,3)>365),
               merradaily(:,:,60) = []; %% sorry no leap years in future
            end;
            mmstart = [1 32 60 91 121 152 182 213 244 274 305 335];
            mmlength = [31 28 31 30 31 30 31 31 30 31 30 31];
            merrafuture = merradaily*NaN;
            for xmonth=1:12,
               thismonthidx = [(mmstart(xmonth)):((mmstart(xmonth)+mmlength(xmonth))-1)];
               for dd=1:length(thismonthidx),
                  dayidx = thismonthidx(dd);
                  merrafuture(:,:,dayidx) = (squeeze(gcmdelta(:,:,xmonth)))+(squeeze(merradaily(:,:,dayidx)));
               end; % day
            end;   % month

            save(['MERRA2_DailySSTs_' thisgcm '_' thisscen '_' num2str(slicestart(tt)) '-' num2str(sliceend(tt)) '_yy' num2str(xbaseyear) '.mat'],'merrafuture');

         end; % MERRA-2 year
      end;    % slice
   end;       % scen

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%% DHW maps: max, average annual max
   %%%%%%%%%%%%%%%%%%%  and change between baseline and future

   merraregionaldhwyy = importdata(['MERRA2_RegionalDHW_AnnualMax.mat']);
   merrameanmaxdhw = importdata(['MERRA2_1980-2009_AvgMaxDHW.mat']);
   maxmonth = importdata(['MERRA2_MesoAmericanReef_MaxMeanMonth_1985-1993.mat']);

   mmstart = [1 32 60 91 121 152 182 213 244 274 305 335];
   mmlength = [31 28 31 30 31 30 31 31 30 31 30 31];
   yearlength = 365;

   for xscen=1:length(allscen),
      thisscen = allscen{xscen};
      thisscennice = longscen{xscen};

      for tt=1:length(slices),
         thisslice = slices{tt};
         sliceyrs = [slicestart(tt):sliceend(tt)];

         for xyear=1:30,

            filename = ['MERRA2_DailySSTs_' thisgcm '_' thisscen '_' num2str(sliceyrs(1)) '-' num2str(sliceyrs(end)) '_yy' num2str(xyear) '.mat'];
            subdata = load([filename]);
            subdata = subdata.merrafuture;

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

            %%% Calculate hotspots by subtracting max monthly SST
            %%% climatology from each twice-weekly SST

            allhotspots = yyhalfweeklydata*NaN;
            for xweek=1:size(yyhalfweeklydata,3),
               allhotspots(:,:,xweek) = yyhalfweeklydata(:,:,xweek)-maxmonth;
            end;
            allhotspots(allhotspots < 0) = NaN;
            clear yyhalfweeklydata
            save([thisgcm '_MesoAmericanReef_HalfWeeklyHotSpots_' thisscen '_' num2str(sliceyrs(1)) '-' num2str(sliceyrs(end)) '_yy' num2str(xyear) '.mat'],'allhotspots');

            %%% Calculate Degree Heating Weeks using the following equation
            %%% DHW = 0.5 * summation of previous 24 twice-weekly hotspots > 1 degree

            alldhw = allhotspots*NaN;

            if (xyear>1),

              lastyear = xyear-1;
              lastyyhotspots = importdata([thisgcm '_MesoAmericanReef_HalfWeeklyHotSpots_' thisscen '_' num2str(sliceyrs(1)) '-' num2str(sliceyrs(end)) '_yy' num2str(lastyear) '.mat']);

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

            end; % alldhw
            save([thisgcm '_MesoAmericanReef_HalfWeeklyDHW_' thisscen '_' num2str(sliceyrs(1)) '-' num2str(sliceyrs(end)) '_yy' num2str(xyear) '.mat'],'alldhw');
            clear allhotspots

            %%%% create regional time series of half weekly DHW

            smalldhw = squeeze(alldhw(thislatidx,thislonidx,:));
            for xweek=1:size(smalldhw,3),
               smalldhw(:,:,xweek) = squeeze(smalldhw(:,:,xweek)).*smallcmipmask;
            end;
            areaweights = smalldhw*NaN;
            for xweek=1:size(smalldhw,3),
               subarea = squeeze(area(thislatidx,thislonidx));
               subarea(isnan(squeeze(smalldhw(:,:,xweek)))) = NaN;
               areaweights(:,:,xweek) = (subarea./nansum(nansum(subarea)));
            end;

            smalldhwweighted = smalldhw*NaN;
            for xweek=1:size(smalldhw,3),
               smalldhwweighted(:,:,xweek) = (smalldhw(:,:,xweek)).*(areaweights(:,:,xweek));
            end;
            %% full half-weekly time series
            regionaldhw = squeeze(nansum(squeeze(nansum(smalldhwweighted))));
            save([thisgcm '_RegionalDHW_HalfWeekly_' thisscen '_' num2str(sliceyrs(1)) '-' num2str(sliceyrs(end)) '_yy' num2str(xyear) '.mat'],'regionaldhw');

            clear alldhw
            clear lastyyhotspots

         end; % 30 slice years

         annualmaxdhw = ones(length(latidx),length(lonidx),30)*NaN;
         for xyear=1:30,
            thisyydhw = importdata([thisgcm '_MesoAmericanReef_HalfWeeklyDHW_' thisscen '_' num2str(sliceyrs(1)) '-' num2str(sliceyrs(end)) '_yy' num2str(xyear) '.mat']);
            annualmaxdhw(:,:,xyear) = nanmax(permute(thisyydhw,[3 1 2]));
            clear thisyydhw
         end;
         allmeanmaxdhw = squeeze(nanmean(annualmaxdhw,3));
         clear annualmaxdhw
         diffavgmaxdhw = allmeanmaxdhw-merrameanmaxdhw;

         if (tt == 5),

            hexcolors = importdata('dhwpaper_colorbar21.mat');
            p = ones(length(hexcolors),3)*NaN;
            for rgb=1:size(p,1),
               p(rgb,:) = hex2rgb(hexcolors(rgb));
            end;

            f = figure; colormap(p);
            acr_pcolormapr4_nocoast(allmeanmaxdhw.*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[7.5 28.5],['MERRA2 (' gcmtitle ' Adjusted) Mean Annual Max Degree Heating Weeks ' thisscennice ' ' num2str(sliceyrs(1)) '-' num2str(sliceyrs(end))],['degree heating weeks (degrees C)']);
            hold on;
            h1 = borders('countries','color','k','linewidth',1);
            hold on;
            plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',1.4);
            print(f,'-dpng',[thisgcm '_MERRA2_' thisscen '_' num2str(sliceyrs(1)) '-' num2str(sliceyrs(end)) '_AvgMaxDHW.png']);
            print(f,'-depsc',[thisgcm '_MERRA2_' thisscen '_' num2str(sliceyrs(1)) '-' num2str(sliceyrs(end)) '_AvgMaxDHW.eps']);
            close all

            f = figure; colormap(p);
            acr_pcolormapr4_nocoast(diffavgmaxdhw.*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[7.5 28.5],['MERRA2 (' gcmtitle ' Adjusted) Change in Mean Annual Max DHW ' thisscennice ' ' num2str(sliceyrs(1)) '-' num2str(sliceyrs(end))],['change in mean max degree heating weeks (degrees C) from baseline']);
            hold on;
            h1 = borders('countries','color','k','linewidth',1);
            hold on;
            plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',1.4);
            print(f,'-dpng',[thisgcm '_MERRA2_' thisscen '_' num2str(sliceyrs(1)) '-' num2str(sliceyrs(end)) '_DeltaAvgMaxDHW.png']);
            print(f,'-depsc',[thisgcm '_MERRA2_' thisscen '_' num2str(sliceyrs(1)) '-' num2str(sliceyrs(end)) '_DeltaAvgMaxDHW.eps']);
            close all

         end;

         save([thisgcm '_MERRA2_' thisscen '_' num2str(sliceyrs(1)) '-' num2str(sliceyrs(end)) '_AvgMaxDHW.mat'],'allmeanmaxdhw');
         save([thisgcm '_MERRA2_' thisscen '_' num2str(sliceyrs(1)) '-' num2str(sliceyrs(end)) '_DeltaAvgMaxDHW.mat'],'diffavgmaxdhw');

         %% annual max time series of DHW
         regionaldhwyy = ones(1,30)*NaN;

         for xyear=1:30,
            thisannualdhw = importdata([thisgcm '_RegionalDHW_HalfWeekly_' thisscen '_' num2str(sliceyrs(1)) '-' num2str(sliceyrs(end)) '_yy' num2str(xyear) '.mat']);
            regionaldhwyy(xyear) = nanmax(thisannualdhw);
         end;

         save([thisgcm '_RegionalDHW_AnnualMax_' thisscen '_' num2str(sliceyrs(1)) '-' num2str(sliceyrs(end)) '.mat'],'regionaldhwyy');

      end;  %% tt
   end;     %% xscen

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%% FREQUENCY of hotspots above 1C

   merradhwfreq = ones(length(latidx),length(lonidx),length(histyears))*NaN;
   merradhwduration = ones(length(latidx),length(lonidx),length(histyears))*NaN;

   for xyear=1:length(histyears),
      thisyear = histyears(xyear);
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

   for xscen=1:length(allscen),
      thisscen = allscen{xscen};
      thisscennice = longscen{xscen};

      for tt=1:length(slices),
         thisslice = slices{tt};
         sliceyrs = [slicestart(tt):sliceend(tt)];
         cmipdhwfreq = ones(length(latidx),length(lonidx),30)*NaN;
         cmipdhwduration = ones(length(latidx),length(lonidx),30)*NaN;

         for xyear=1:30,

            allhotspots = importdata([thisgcm '_MesoAmericanReef_HalfWeeklyHotSpots_' thisscen '_' num2str(sliceyrs(1)) '-' num2str(sliceyrs(end)) '_yy' num2str(xyear) '.mat']);

            for ii=1:length(lonidx),
               for jj=1:length(latidx),
                  cmipsstpoint = squeeze(allhotspots(jj,ii,:));
                  cmipdhwfreq(jj,ii,xyear) = length(find(cmipsstpoint>1));
                  clear streak
                  if (cmipsstpoint(1) > 1),
                     streak(1) = 1;
                  else
                     streak(1) = 0;
                  end;
                  for dd=2:length(cmipsstpoint),
                     if (cmipsstpoint(dd) > 1),
                        streak(dd) = streak(dd-1) + 1;
                        streak(dd-1) = 0;
                     else
                        streak(dd) = 0;
                     end;
                  end;
                  cmipdhwduration(jj,ii,xyear) = nanmax(streak);
               end;
            end;
         end;  %% xyear

         cmipdhwfreq = squeeze(nanmean(cmipdhwfreq,3));
         cmipdhwduration = squeeze(nanmean(cmipdhwduration,3));

         freqdiffmap = cmipdhwfreq-merradhwfreq;
         durationdiffmap = cmipdhwduration-merradhwduration;

         if (tt == 5),

            % for this one in R use
            % colorbar = colorRampPalette(c("#F7F4F9","#E7E1EF","#D4B9DA","#C994C7","#DF65B0","#E7298A","#CE1256","#91003F"))(n = 32)
            %hexcolors = importdata('dhwpaper_purplered32.mat');
            %p = ones(length(hexcolors),3)*NaN;
            %for rgb=1:size(p,1),
            %   p(rgb,:) = hex2rgb(hexcolors(rgb));
            %end;

            % colorbar = colorRampPalette(c('#aaffc3','#42d4f4','#f032e6','#911eb4','#e6194B'))(n = 32)
            hexcolors = importdata('dhwpaper_tealmagenta32.mat');
            p = ones(length(hexcolors),3)*NaN;
            for rgb=1:size(p,1),
               p(rgb,:) = hex2rgb(hexcolors(rgb));
            end;

            f = figure; colormap(p);
            acr_pcolormapr4_nocoast((cmipdhwfreq/2).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[3.5 35.5],['MERRA2 (' gcmtitle ' Adjusted) Mean Annual Frequency of Hotspots>1C ' thisscennice ' ' num2str(sliceyrs(1)) '-' num2str(sliceyrs(end))],['number of weeks with >1C anomaly from baseline']);
            hold on;
            h1 = borders('countries','color','k','linewidth',1);
            hold on;
            plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',1.4);
            print(f,'-dpng',[thisgcm '_MERRA2_HotSpots_Frequency_' thisscen '_' num2str(sliceyrs(1)) '-' num2str(sliceyrs(end)) '.png']);
            print(f,'-depsc',[thisgcm '_MERRA2_HotSpots_Frequency_' thisscen '_' num2str(sliceyrs(1)) '-' num2str(sliceyrs(end)) '.eps']);
            close all

            f = figure; colormap(p);
            acr_pcolormapr4_nocoast((freqdiffmap/2).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[3.5 35.5],['MERRA2 (' gcmtitle ' Adjusted) Change in Mean Annual Frequency of Hotspots>1C ' thisscennice ' ' num2str(sliceyrs(1)) '-' num2str(sliceyrs(end))],['change in number of weeks with >1C anomaly from baseline']);
            hold on;
            h1 = borders('countries','color','k','linewidth',1);
            hold on;
            plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',1.4);
            print(f,'-dpng',[thisgcm '_MERRA2_HotSpots_DeltaFrequency_' thisscen '_' num2str(sliceyrs(1)) '-' num2str(sliceyrs(end)) '.png']);
            print(f,'-depsc',[thisgcm '_MERRA2_HotSpots_DeltaFrequency_' thisscen '_' num2str(sliceyrs(1)) '-' num2str(sliceyrs(end)) '.eps']);
            close all

         end;

         save([thisgcm '_MERRA2_HotSpots_Frequency_' thisscen '_' num2str(sliceyrs(1)) '-' num2str(sliceyrs(end)) '.mat'],'cmipdhwfreq');
         save([thisgcm '_MERRA2_HotSpots_DeltaFrequency_' thisscen '_' num2str(sliceyrs(1)) '-' num2str(sliceyrs(end)) '.mat'],'freqdiffmap');

         if (tt == 5),

            hexcolors = importdata('dhwpaper_tealmagenta24.mat');
            p = ones(length(hexcolors),3)*NaN;
            for rgb=1:size(p,1),
               p(rgb,:) = hex2rgb(hexcolors(rgb));
            end;

            f = figure; colormap(p);
            acr_pcolormapr4_nocoast((cmipdhwduration/2).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[4.5 28.5],['MERRA2 (' gcmtitle ' Adjusted) Mean Annual Max Duration of Hotspots>1C ' thisscennice ' ' num2str(sliceyrs(1)) '-' num2str(sliceyrs(end))],['maximum consecutive weeks with >1C anomaly from baseline']);
            hold on;
            h1 = borders('countries','color','k','linewidth',1);
            hold on;
            plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',1.4);
            print(f,'-dpng',[thisgcm '_MERRA2_HotSpots_MaxDuration_' thisscen '_' num2str(sliceyrs(1)) '-' num2str(sliceyrs(end)) '.png']);
            print(f,'-depsc',[thisgcm '_MERRA2_HotSpots_MaxDuration_' thisscen '_' num2str(sliceyrs(1)) '-' num2str(sliceyrs(end)) '.eps']);
            close all

            f = figure; colormap(p);
            acr_pcolormapr4_nocoast((durationdiffmap/2).*cmipmask,sublat,sublon,[slat nlat],[wlon elon],[4.5 28.5],['MERRA2 (' gcmtitle ' Adjusted) Change in Mean Annual Max Duration of Hotspots>1C ' thisscennice ' ' num2str(sliceyrs(1)) '-' num2str(sliceyrs(end))],['change from baseline in maximum consecutive weeks with >1C anomaly']);
            hold on;
            h1 = borders('countries','color','k','linewidth',1);
            hold on;
            plotm([regionslat regionnlat regionnlat regionslat regionslat],[regionwlon regionwlon regionelon regionelon regionwlon],'-','color','k','linewidth',1.4);
            print(f,'-dpng',[thisgcm '_MERRA2_HotSpots_DeltaMaxDuration_' thisscen '_' num2str(sliceyrs(1)) '-' num2str(sliceyrs(end)) '.png']);
            print(f,'-depsc',[thisgcm '_MERRA2_HotSpots_DeltaMaxDuration_' thisscen '_' num2str(sliceyrs(1)) '-' num2str(sliceyrs(end)) '.eps']);
            close all

         end;

         save([thisgcm '_MERRA2_HotSpots_MaxDuration_' thisscen '_' num2str(sliceyrs(1)) '-' num2str(sliceyrs(end)) '.mat'],'cmipdhwduration');
         save([thisgcm '_MERRA2_HotSpots_DeltaMaxDuration_' thisscen '_' num2str(sliceyrs(1)) '-' num2str(sliceyrs(end)) '.mat'],'durationdiffmap');

      end;  % slice
   end;     % scen

end;     %%% xgcm
