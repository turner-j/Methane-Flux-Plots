%% Clear workspace, clear variables, and clear figures
clear all;close all;clc;
load('USLos2020_2019_calcresp.mat')
load('USALQ2020_2019_calcresp.mat')
load('woodruffprecip.mat')
load('dischargevars.mat')
load('aveFCH4.mat')
%% Replace placeholders for missing values with NaN
USALQ = [USALQHH2019; USALQHH2020]; 

for i= 3:width(USALQ)
    USALQ.(i)(USALQ.(i)==-9999) = NaN;
end

%% Make daily averages
CH4 = USALQ.FCH4_RF_filled;% HH data
NEE = USALQ.NEE_F;

daily_CH4 = dailyaverage2(USALQ.FCH4_RF_filled);
dailyCH4 = daily_CH4(:,1);% Daily data

GPP = USALQ.GPP_F;% HH data
daily_GPP = dailyaverage2(USALQ.GPP_F);
dailyGPP = daily_GPP(:,1);% Daily data

dates = USALQ.TIMESTAMP_END(48*daily_CH4(:,2));
GPPdates = USALQ.TIMESTAMP_END(48*daily_GPP(:,2));

%% Calculate the cumulative sum of FCH4 for each year
CH4grams = gramconvertnmol(CH4(year(USALQ.TIMESTAMP_END)==2019));
cumulativeCH4 = cumsum(CH4grams,'omitnan');
cumulativeCH4(end)
stderror = std(cumulativeCH4)/sqrt(length(cumulativeCH4))

CH4grams = gramconvertnmol(CH4(year(USALQ.TIMESTAMP_END)==2020));
cumulativeCH4 = cumsum(CH4grams,'omitnan');
cumulativeCH4(end)
stderror = std(cumulativeCH4)/sqrt(length(cumulativeCH4))

NEEgrams = gramconvert(NEE(year(USALQ.TIMESTAMP_END)==2019));
cumulativeNEE = cumsum(NEEgrams,'omitnan');
cumulativeNEE(end)
stderror = std(cumulativeNEE)/sqrt(length(cumulativeNEE))

NEEgrams = gramconvert(NEE(year(USALQ.TIMESTAMP_END)==2020));
cumulativeNEE = cumsum(NEEgrams,'omitnan');
cumulativeNEE(end)
stderror = std(cumulativeNEE)/sqrt(length(cumulativeNEE))
