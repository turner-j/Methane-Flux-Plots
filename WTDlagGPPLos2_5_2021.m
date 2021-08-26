%% Clear workspace, clear variables, and clear figures
clear all;close all;clc;
load('USLos2020_2019_calcresp.mat')
load('woodruffprecip.mat')
load('dischargevars.mat')
%% Replace placeholders for missing values with NaN
USLosHH2019.Properties.VariableNames = USLosHH2020.Properties.VariableNames;
USLos = [USLosHH2019;USLosHH2020];

for i= 3:width(USLos)
    USLos.(i)(USLos.(i)==-9999) = NaN;
end

%% Sort data by the half hour for whichever month or months you want
precip = woodruffprecipS3.Precipitationin(2:end,:)*2.54;% convert in to cm

dis = USGSLos;
dis.discharge = (dis.discharge_cfs).*0.028316847;% convert cfs to m3/s
yr2 = year(dis.datetime);
mon2 = month(dis.datetime);

%% Make daily averages
mo = month(USLos.TIMESTAMP_END);
monthlyCH4 = accumarray(mo,USLos.FCH4_RF_filled,[],@nanmean);

daily_CH4 = dailyaverage2(USLos.FCH4_RF_filled);
dailyCH4 = daily_CH4(:,1);
CH4dates = USLos.TIMESTAMP_END(daily_CH4(:,2).*48);

daily_Tairf = dailyaverage2(USLos.TA_F_1_1_1);
dailyTairf = daily_Tairf(:,1);

daily_WTD = dailyaverage2(USLos.WTD);
dailyWTD = daily_WTD(:,1);

daily_GPP = dailyaverage2(USLos.GPP_F);
dailyGPP = daily_GPP(:,1);

daily_VPD = dailyaverage2(USLos.VPD_F_1_1_1);
dailyVPD = daily_VPD(:,1);
%% Lag Analysis of WTD and FCH4 2019
yr = year(CH4dates);
Z = [];
dailyGPP1 = dailyGPP(year(CH4dates)==2019);
dailyWTD1 = dailyWTD(year(CH4dates)==2019);

%For WTD
for k = 1:200
    B = NaN(length(dailyWTD1),1);
    B((k+1):end) = dailyWTD1(1:(end-k));% shift WTD down
    Z(:,k) = B;
end

Z = [dailyWTD1 Z(:,1:200)];
% Remove any missing data from both datasets
TF = isnan(Z);
TF2 = isnan(dailyGPP1);
location = [];

for i = 1:size(TF,1)
    for j = 1:200
        if (TF(i,j)==1) || (TF2(i)==1)
            location = [location;[i j]];
        end
    end
end

A = repmat(dailyGPP1,1,200);

for i = 1:length(location)
    Z(location(i,1),location(i,2)) = NaN;
    A(location(i,1),location(i,2)) = NaN;
end

% Calculate cross correlations

for i = 1:200
    dayWTD = Z(:,i);
    dayGPP = A(:,i);
    [rho3(i),pval3(i)] = corr(dayWTD(~isnan(dayWTD)),dayGPP(~isnan(dayGPP)));
end
%% For 2020
Z = [];
dailyGPP = dailyGPP(year(CH4dates)==2020);
dailyWTD = dailyWTD(year(CH4dates)==2020);

for k = 1:200
    B = NaN(length(dailyWTD),1);
    B((k+1):end) = dailyWTD(1:(end-k));% shift WTD down
    Z(:,k) = B;
end

Z = [dailyWTD Z(:,1:200)];
% Remove any missing data from both datasets
TF = isnan(Z);
TF2 = isnan(dailyGPP);
location = [];

for i = 1:size(TF,1)
    for j = 1:200
        if (TF(i,j)==1) || (TF2(i)==1)
            location = [location;[i j]];
        end
    end
end

A = repmat(dailyGPP,1,200);

for i = 1:length(location)
    Z(location(i,1),location(i,2)) = NaN;
    A(location(i,1),location(i,2)) = NaN;
end

% Calculate cross correlations

for i = 1:200
    dayWTD = Z(:,i);
    dayGPP = A(:,i);
    [rho(i),pval(i)] = corr(dayWTD(~isnan(dayWTD)),dayGPP(~isnan(dayGPP)));
end

%% Plot 
lag2 = 1:200;

figure()
subplot(1,2,1)
hold on
plot(lag2,rho3,'o','MarkerEdgeColor','k','MarkerSize',11)
plot(lag2(pval3<0.05 & rho3>0),rho3(pval3<0.05 & rho3>0),'o','MarkerFaceColor','b','MarkerEdgeColor','k','MarkerSize',11);
hold off
xlabel('time lag (days)')
ylabel('correlation coefficient')
% title('Lag Corr of WTD & GPP US-Los 2019')
title('US-Los 2019')
set(gca,'FontSize',17)
ylim([0 1])
xlim([0 200])

subplot(1,2,2)
hold on
plot(lag2,rho,'o','MarkerEdgeColor','k','MarkerSize',11)
plot(lag2(pval<0.05 & rho>0),rho(pval<0.05 & rho>0),'o','MarkerFaceColor','b','MarkerEdgeColor','k','MarkerSize',11);
hold off
xlabel('time lag (days)')
ylabel('correlation coefficient')
% title('Lag Corr of WTD & GPP US-Los 2020')
title('US-Los 2020')
set(gca,'FontSize',17)
ylim([0 1])
xlim([0 200])