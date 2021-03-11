%% Load Data
clear all;close all
clc;

load('USALQ2020_2019_calcresp.mat')
%% Replace missing placeholder values with NaNs
USALQ = [USALQHH2019; USALQHH2020]; 

for i= 3:width(USALQ)
    USALQ.(i)(USALQ.(i)==-9999) = NaN;
end
%% Calculate Average Monthly FCH4 and NEE
x = month(USALQ.TIMESTAMP_END);

monthlyNEE = accumarray(x,USALQ.NEE_F,[],@nanmean);
monthlyFCH4 = accumarray(x,USALQ.FCH4_RF_filled,[],@nanmean);

monthlyNEEnums = accumarray(x,~isnan(USALQ.NEE_F),[]);
monthlyFCH4nums = accumarray(x,~isnan(USALQ.FCH4_RF_filled),[]);

monthlyNEEstd = accumarray(x,USALQ.NEE_F,[],@nanstd);
monthlyFCH4std = accumarray(x,USALQ.FCH4_RF_filled,[],@nanstd);

err = monthlyNEEstd./sqrt(monthlyNEEnums);
err2 = monthlyFCH4std./sqrt(monthlyFCH4nums);

FCH4variance = accumarray(x,USALQ.FCH4_RF_filled,[],@nanvar);
NEEvariance = accumarray(x,USALQ.NEE_F,[],@nanvar);

nanmean(USALQ.FCH4_RF_filled)
%% Create Bar Plots with Error Bars
months = 1:12;

figure()
hBar = bar(months,monthlyNEE,'hist');                
hold on
er = errorbar(months,monthlyNEE,err);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off
title('Bar Plot of Monthly NEE ALQ 2019-20')
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
title('Bar Plot of Monthly FCH_4 ALQ 2019-20')
xlabel('Month')
ylabel('CH_4 flux (\etamol m^-^2s^-^1)')

% Add a colormap to represent variance
CMap = colormap(winter);
set(hBar2,'CData',FCH4variance,'CDataMapping','scaled');
c = colorbar;
c.Label.String = 'Variance';
caxis([min(FCH4variance) max(FCH4variance)]);
%% Calculate Missing Data Each Month
TF = 100*sum(isnan(USALQ.FCH4)/length(USALQ.FCH4))
TF = 100*sum(isnan(USALQ.FCH4(ismember(x,4:10))))/length(USALQ.FCH4(ismember(x,4:10)))
TF = 100*sum(isnan(USALQ.FCH4(~ismember(x,4:10))))/length(USALQ.FCH4(~ismember(x,4:10)))
TF = 100*sum(isnan(USALQ.FCH4_RF_filled)/length(USALQ.FCH4_RF_filled))

x = month(USALQ.TIMESTAMP_END(year(USALQ.TIMESTAMP_END)==2019));
missingNEE19 = accumarray(x,isnan(USALQ.NEE(year(USALQ.TIMESTAMP_END)==2019)),[]);
monthlyNEE19 = accumarray(x,~isnan(USALQ.NEE(year(USALQ.TIMESTAMP_END)==2019)),[]);
missingtotalNEE19 = (missingNEE19./monthlyNEE19)*100;

missingFCH419 = accumarray(x,isnan(USALQ.FCH4(year(USALQ.TIMESTAMP_END)==2019)),[]);
monthlyFCH419 = accumarray(x,~isnan(USALQ.FCH4(year(USALQ.TIMESTAMP_END)==2019)),[]);
missingtotalFCH419 = (missingFCH419./monthlyFCH419)*100;
TF = 100*sum(isnan(USALQ.FCH4(year(USALQ.TIMESTAMP_END)==2019))/length(USALQ.FCH4(year(USALQ.TIMESTAMP_END)==2019)))

yr = 2019*ones(12,1);

x2 = month(USALQ.TIMESTAMP_END(year(USALQ.TIMESTAMP_END)==2020));
missingNEE20 = accumarray(x2,isnan(USALQ.NEE(year(USALQ.TIMESTAMP_END)==2020)),[]);
monthlyNEE20 = accumarray(x2,~isnan(USALQ.NEE(year(USALQ.TIMESTAMP_END)==2020)),[]);
missingtotalNEE20 = (missingNEE20./monthlyNEE20)*100;

missingFCH420 = accumarray(x2,isnan(USALQ.FCH4(year(USALQ.TIMESTAMP_END)==2020)),[]);
monthlyFCH420 = accumarray(x2,~isnan(USALQ.FCH4(year(USALQ.TIMESTAMP_END)==2020)),[]);
missingtotalFCH420 = (missingFCH420./monthlyFCH420)*100;
TF = 100*sum(isnan(USALQ.FCH4(year(USALQ.TIMESTAMP_END)==2020))/length(USALQ.FCH4(year(USALQ.TIMESTAMP_END)==2020)))

yr2 = 2020*ones(length(unique(x2)),1);

T = table([yr;yr2], [unique(x);unique(x2)] ,[missingtotalNEE19;missingtotalNEE20] ,[missingtotalFCH419;missingtotalFCH420]);
T.Properties.VariableNames = {'Year','Month','Missing NEE','Missing FCH4'};
