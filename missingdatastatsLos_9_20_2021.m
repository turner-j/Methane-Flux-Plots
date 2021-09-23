%% Load Data
clear all;close all
clc;

load('USLos2020_2019_calcresp.mat')
%% Replace missing placeholder values with NaNs
USLosHH2019.Properties.VariableNames = USLosHH2020.Properties.VariableNames;
USLos = [USLosHH2019;USLosHH2020];

for i= 3:width(USLos)
    USLos.(i)(USLos.(i)==-9999) = NaN;
end
%% Calculate Average Monthly FCH4 and NEE
x = month(USLos.TIMESTAMP_END);

monthlyWTD = accumarray(month(USLos.TIMESTAMP_END),USLos.WTD,[],@nanmean);% This is for Los 

monthlyNEE = accumarray(x,USLos.NEE_F,[],@nanmean);
monthlyFCH4 = accumarray(x,USLos.FCH4_RF_filled,[],@nanmean);

monthlyNEEnums = accumarray(x,~isnan(USLos.NEE_F),[]);
monthlyFCH4nums = accumarray(x,~isnan(USLos.FCH4_RF_filled),[]);

monthlyNEEstd = accumarray(x,USLos.NEE_F,[],@nanstd);
monthlyFCH4std = accumarray(x,USLos.FCH4_RF_filled,[],@nanstd);

err = monthlyNEEstd./sqrt(monthlyNEEnums);
err2 = monthlyFCH4std./sqrt(monthlyFCH4nums);

FCH4variance = accumarray(x,USLos.FCH4_RF_filled,[],@nanvar);
NEEvariance = accumarray(x,USLos.NEE_F,[],@nanvar);

monthlytempLos = accumarray(month(USLos.TIMESTAMP_END),USLos.TA_F_1_1_1,[],@nanmean);
%% Scatterplot
figure()
h = scatter(monthlyWTD,monthlyFCH4,'filled');% third variable represents marker size
xlabel('WTD (m)')
ylabel('FCH4 (\etamol CH_4 m^-^2s^-^1)')

% Add a colormap to represent variance
CMap = colormap(parula);
set(h,'CData',monthlytempLos);
c = colorbar;
c.Label.String = 'Monthly Average Temperature (C)';
caxis([min(monthlytempLos) max(monthlytempLos )]);
%% Create Bar Plots with Error Bars
months = 1:12;

figure()
hBar = bar(months,monthlyNEE,'hist');                
hold on
er = errorbar(months,monthlyNEE,err);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off
title('Bar Plot of Monthly NEE Los 2019-20')
xlabel('Month')
ylabel('CO_2 flux (\mumol CO_2 m^-^2s^-^1)')

% Add a colormap to represent variance
CMap = colormap(spring);
set(hBar,'CData',NEEvariance,'CDataMapping','scaled');
c = colorbar;
c.Label.String = 'Variance';
caxis([min(NEEvariance) max(NEEvariance)]);

figure()
hBar2 = bar(months,monthlyFCH4,'hist');                
hold on
er = errorbar(months,monthlyFCH4,err2');    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off
title('Bar Plot of Monthly FCH_4 Los 2019-20')
xlabel('Month')
ylabel('CH_4 flux (\etamol m^-^2s^-^1)')

% Add a colormap to represent variance
CMap = colormap(winter);
set(hBar2,'CData',FCH4variance,'CDataMapping','scaled');
c = colorbar;
c.Label.String = 'Variance';
caxis([min(FCH4variance) max(FCH4variance)]);
%% Calculate Missing Data
missingNEE19 = sum(isnan(USLos.NEE(year(USLos.TIMESTAMP_END)==2019)));
NEE19 = sum(~isnan(USLos.NEE_F(year(USLos.TIMESTAMP_END)==2019)));
missingtotalNEE19 = (missingNEE19./length(USLos.NEE_F(year(USLos.TIMESTAMP_END)==2019)))*100;
totalNEE19 = (NEE19./length(USLos.NEE_F(year(USLos.TIMESTAMP_END)==2019)))*100;

missingNEE20 = sum(isnan(USLos.NEE(year(USLos.TIMESTAMP_END)==2020)));
NEE20 = sum(~isnan(USLos.NEE_F(year(USLos.TIMESTAMP_END)==2020)));
missingtotalNEE20 = (missingNEE20./length(USLos.NEE_F(year(USLos.TIMESTAMP_END)==2020)))*100;
totalNEE20 = (NEE20./length(USLos.NEE_F(year(USLos.TIMESTAMP_END)==2020)))*100;

mo = month(USLos.TIMESTAMP_END);
missingNEE = sum(isnan(USLos.NEE(ismember(mo,4:10))));
missingtotalNEEseason = (missingNEE./length(USLos.NEE_F(ismember(mo,4:10))))*100;

missingNEE = sum(isnan(USLos.NEE(~ismember(mo,4:10))));
missingtotalNEEoffseason  = (missingNEE./length(USLos.NEE_F(~ismember(mo,4:10))))*100;
%%
missingFCH419 = accumarray(x,isnan(USLos.FCH4_1_1_1(year(USLos.TIMESTAMP_END)==2019)),[]);
monthlyFCH419 = accumarray(x,~isnan(USLos.FCH4_1_1_1(year(USLos.TIMESTAMP_END)==2019)),[]);
missingtotalFCH419 = (missingFCH419./monthlyFCH419)*100;
TF = 100*sum(isnan(USLos.FCH4_1_1_1(year(USLos.TIMESTAMP_END)==2019))/length(USLos.FCH4_1_1_1(year(USLos.TIMESTAMP_END)==2019)))

yr = 2019*ones(12,1);

x2 = month(USLos.TIMESTAMP_END(year(USLos.TIMESTAMP_END)==2020));
missingNEE20 = accumarray(x2,isnan(USLos.NEE_F(year(USLos.TIMESTAMP_END)==2020)),[]);
monthlyNEE20 = accumarray(x2,~isnan(USLos.NEE_F(year(USLos.TIMESTAMP_END)==2020)),[]);
missingtotalNEE20 = (missingNEE20./monthlyNEE20)*100;

missingFCH420 = accumarray(x2,isnan(USLos.FCH4_1_1_1(year(USLos.TIMESTAMP_END)==2020)),[]);
monthlyFCH420 = accumarray(x2,~isnan(USLos.FCH4_1_1_1(year(USLos.TIMESTAMP_END)==2020)),[]);
missingtotalFCH420 = (missingFCH420./monthlyFCH420)*100;
TF = 100*sum(isnan(USLos.FCH4_1_1_1(year(USLos.TIMESTAMP_END)==2020))/length(USLos.FCH4_1_1_1(year(USLos.TIMESTAMP_END)==2020)))