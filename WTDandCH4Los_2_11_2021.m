%% Load Data
clear all;close all;
clc;

load('USLos2020_2019_calcresp.mat')
%% Change -9999s to NaNs
USLosHH2019.Properties.VariableNames = USLosHH2020.Properties.VariableNames;
USLos = [USLosHH2019;USLosHH2020];

for i= 3:width(USLos)
    USLos.(i)(USLos.(i)==-9999) = NaN;
end

% indexing
% yr = year(USLos.TIMESTAMP_END);
% ind = find(ismember(yr,2019));
% USLos = USLos(ind,:);
%% Make daily averages
NEE = USLos.NEE_F;
CH4 = USLos.FCH4_RF_filled;% HH data
daily_CH4 = dailyaverage2(USLos.FCH4_RF_filled);
dailyCH4 = daily_CH4(:,1);% Daily data

WTD = USLos.WTD;% HH data
daily_WTD = dailyaverage2(USLos.WTD);
dailyWTD = daily_WTD(:,1);% Daily data

GPP = USLos.GPP_F;% HH data
daily_GPP = dailyaverage2(USLos.GPP_F);
dailyGPP = daily_GPP(:,1);% Daily data

% Eliminating missing data from both datasets to ensure we are comparing
% data collected at the same times
dailyWTD(isnan(dailyCH4))=NaN;
dailyCH4(isnan(dailyWTD))=NaN;

% [R,P] = corrcoef(dailyWTD(~isnan(dailyWTD)),dailyCH4(~isnan(dailyCH4)))

WTD(isnan(CH4))=NaN;
CH4(isnan(WTD))=NaN;

[R,P] = corrcoef(WTD(~isnan(WTD)),CH4(~isnan(CH4)))
%% Calculate the cumulative sum of FCH4 for each year
CH4grams = gramconvertnmol(CH4(year(USLos.TIMESTAMP_END)==2019));
cumulativeCH4 = cumsum(CH4grams,'omitnan');
cumulativeCH4(end)
stderror = std(cumulativeCH4)/sqrt(length(cumulativeCH4))

CH4grams = gramconvertnmol(CH4(year(USLos.TIMESTAMP_END)==2020));
cumulativeCH4 = cumsum(CH4grams,'omitnan');
cumulativeCH4(end)
stderror = std(cumulativeCH4)/sqrt(length(cumulativeCH4))

NEEgrams = gramconvert(NEE(year(USLos.TIMESTAMP_END)==2019));
cumulativeNEE = cumsum(NEEgrams,'omitnan');
cumulativeNEE(end)
stderror = std(cumulativeNEE)/sqrt(length(cumulativeNEE))

NEEgrams = gramconvert(NEE(year(USLos.TIMESTAMP_END)==2020));
cumulativeNEE = cumsum(NEEgrams,'omitnan');
cumulativeNEE(end)
stderror = std(cumulativeNEE)/sqrt(length(cumulativeNEE))
%% Lag Analysis of WTD and FCH4
Z = [];

for k = 1:90
    B = NaN(length(dailyWTD),1);
    B((k+1):end) = dailyWTD(1:(end-k));% shift WTD up
    Z(:,k) = B;
end

Z = [dailyWTD Z(:,1:89)];
% Remove any missing data from both datasets
TF = isnan(Z);
TF2 = isnan(dailyCH4);
location = [];

for i = 1:size(TF,1)
    for j = 1:90
        if (TF(i,j)==1) || (TF2(i)==1)
            location = [location;[i j]];
        end
    end
end

A = repmat(dailyCH4,1,90);

for i = 1:length(location)
    Z(location(i,1),location(i,2)) = NaN;
    A(location(i,1),location(i,2)) = NaN;
end

% Calculate cross correlations

for i = 1:90
    dayWTD = Z(:,i);
    dayCH4 = A(:,i);
    [rho(i),pval(i)] = corr(dayWTD(~isnan(dayWTD)),dayCH4(~isnan(dayCH4)));
end

% Plot 
lag = 0:89;
w = lag;
v = rho;
f = fit(w',v','smoothingspline','SmoothingParam',1);

figure()
plot(f,w',v')
xlabel('time lag, days')
ylabel('correlation coefficient')
title('Sign. Lag corr. of WTD & FCH_4 US-Los 2019 & 2020')
%% Plot 
dates = USLos.TIMESTAMP_END(48*daily_CH4(:,2));
GPPdates = USLos.TIMESTAMP_END(48*daily_GPP(:,2));

figure()
yyaxis left
plot(dates,dailyCH4,'LineWidth',3)
ylabel('FCH_4 (\etamol CH_4 m^-^2s^-^1)')
yyaxis right
plot(dates,dailyWTD,'LineWidth',3)
ylabel('WTD')
title('FCH_4 and WTD US-Los 2019 & 2020')

%% Replace missing data with a random number from the standard normal distribution
aveCH4 = nanmean(dailyCH4);
aveWTD = nanmean(dailyWTD);
aveGPP = nanmean(dailyGPP);

for i = 1:length(dailyCH4)
    if isnan(dailyCH4(i))||isnan(dailyWTD(i))
        dailyCH4(i) = randn*aveCH4;
        dailyWTD(i) = randn*aveWTD;
    end
end

for i = 1:length(dailyCH4)
    if isnan(dailyCH4(i))||isnan(dailyGPP(i))
        dailyCH4(i) = randn*aveCH4;
        dailyGPP(i) = randn*aveGPP;
    end
end
%% Plot wavelets
figure()
wcoherence(dailyWTD(year(dates)==2019),dailyCH4(year(dates)==2019),days(1))
title('US-Los FCH_4& WTD Wavelet Coherence with Daily Data 2019')

figure()
wcoherence(dailyWTD(year(dates)==2020),dailyCH4(year(dates)==2020),days(1))
title('US-Los FCH_4& WTD Wavelet Coherence with Daily Data 2020')

figure()
wcoherence(dailyGPP(year(GPPdates)==2019),dailyCH4(year(GPPdates)==2019),days(1))
title('US-Los FCH_4& GPP Wavelet Coherence with Daily Data 2019')

figure()
wcoherence(dailyGPP(year(GPPdates)==2020),dailyCH4(year(GPPdates)==2020),days(1))
title('US-Los FCH_4& GPP Wavelet Coherence with Daily Data 2020')
%% Both years
[wcoh,wcs,period,coi] = wcoherence(dailyGPP',dailyCH4',days(1),'phasedisplaythreshold',0.75);

figure()
wcoherence(dailyGPP,dailyCH4,days(1));
title('US-Los FCH_4 & GPP Wavelet Coherence with Daily Data')
CMap =colormap(linspecer);
xticklabels('manual')
xt = xticks;
xticklabels(datestr(dates(xt+1)));
xtickangle(45)
xlabel('')
%% Repeat
figure
period = days(period);
coi = days(coi);
h = pcolor(dates,log2(period),wcoh);
h.EdgeColor = 'none';
ax = gca;
ytick=round(pow2(ax.YTick),1);
ax.YTickLabel=ytick;
% ax.XLabel.String='Days';
ax.YLabel.String='Period (days)';
ax.Title.String = 'US-Los FCH_4 & GPP Wavelet Coherence with Daily Data';
hcol = colorbar;
hcol.Label.String = 'Magnitude-Squared Coherence';
hold on;
plot(ax,dates,log2(coi),'w--','linewidth',2)

% CMap = colormap(hsv);
% map = customcolormap(linspace(0,1,21), {'#ffffff','#f2ffff','#e6ffff','#d9ffff',...
%  	  '#ccffff','#d9ffff','#00ffff','#00f2f7','#00e6f0','#00d9e8','#00cce0',...
%     '#00bfd9','#00a6c9','#008cba','#0080b2','#0066a3',...
%     '#00599c','#00408c','#003385','#001975','#000d6e'});
% CMap = colormap(map);
CMap =colormap(linspecer); 

%% Repeat with half-hourly data
aveCH4 = nanmean(CH4);
aveWTD = nanmean(WTD);
aveGPP = nanmean(GPP);

for i = 1:length(CH4)
    if isnan(CH4(i))||isnan(WTD(i))
        CH4(i) = randn*aveCH4;
        WTD(i) = randn*aveWTD;
    end
end

for i = 1:length(CH4)
    if isnan(CH4(i))||isnan(GPP(i))
        CH4(i) = randn*aveCH4;
        GPP(i) = randn*aveGPP;
    end
end
%% Plot wavelets
figure()
wcoherence(WTD(year(USLos.TIMESTAMP_END)==2019),CH4(year(USLos.TIMESTAMP_END)==2019),days(1/48))
title('US-Los FCH_4& WTD Wavelet Coherence with HH Data 2019')

figure()
wcoherence(WTD(year(USLos.TIMESTAMP_END)==2020),CH4(year(USLos.TIMESTAMP_END)==2020),days(1/48))
title('US-Los FCH_4& WTD Wavelet Coherence with HH Data 2020')

figure()
wcoherence(GPP(year(USLos.TIMESTAMP_END)==2019),CH4(year(USLos.TIMESTAMP_END)==2019),days(1/48))
title('US-Los FCH_4& GPP Wavelet Coherence with HH Data 2019')

figure()
wcoherence(GPP(year(USLos.TIMESTAMP_END)==2020),CH4(year(USLos.TIMESTAMP_END)==2020),days(1/48))
title('US-Los FCH_4& GPP Wavelet Coherence with HH Data 2020')