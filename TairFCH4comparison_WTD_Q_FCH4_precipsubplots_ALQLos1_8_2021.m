%% Clear workspace, clear variables, and clear figures
clear all;close all;clc;
load('USLos2020_2019_calcresp.mat')
load('USALQ2020_2019_calcresp.mat')
load('woodruffprecip.mat')   
load('dischargevars.mat')
load('aveFCH4.mat')
%% Replace placeholders for missing values with NaN
USLosHH2019.Properties.VariableNames = USLosHH2020.Properties.VariableNames;
USLos = [USLosHH2019;USLosHH2020];  

for i= 3:width(USLos)
    USLos.(i)(USLos.(i)==-9999) = NaN;
end

USALQ = [USALQHH2019; USALQHH2020]; 

for i= 3:width(USALQ)
    USALQ.(i)(USALQ.(i)==-9999) = NaN;
end

%% Calculate daily averages
% US-Los
daily_CH4Los = dailyaverage2(USLos.FCH4_RF_filled);
dailyCH4Los = daily_CH4Los(:,1);
datesLos = USLos.TIMESTAMP_END(48*daily_CH4Los(:,2));

daily_WTDLos = dailyaverage2(USLos.WTD);
dailyWTDLos = daily_WTDLos(:,1);

daily_GPPLos = dailyaverage2(USLos.GPP_F);
dailyGPPLos = daily_GPPLos(:,1);

dailyNEELos = dailyaverage2(USLos.NEE_F);
dailyNEELos = dailyNEELos(:,1);

dailyTaLos = dailyaverage2(USLos.TA_F_1_1_1);
dailyTaLos = dailyTaLos(:,1);

daily_Rg = dailyaverage2(USLos.SW_IN_1_1_1);
dailyRg = daily_Rg(:,1);

daily_VPD = dailyaverage2(USLos.VPD_F_1_1_1);
dailyVPD = daily_VPD(:,1);
%% US-ALQ
dailyALQ = dailyaverage2(USALQ.FCH4_RF_filled);
dailyFCH4 = dailyALQ(:,1);

dates = dailyALQ(:,2);
dates = USALQ.TIMESTAMP_END(dates.*48);

daily_GPPALQ = dailyaverage2(USALQ.GPP_F);
dailyGPPALQ = daily_GPPALQ(:,1);

dailyNEEALQ = dailyaverage2(USALQ.NEE_F);
dailyNEEALQ = dailyNEEALQ(:,1);

dailyTaALQ = dailyaverage2(USALQ.TA_F);
dailyTaALQ = dailyTaALQ(:,1);

daily_RgALQ = dailyaverage2(USALQ.SW_IN);
dailyRgALQ = daily_RgALQ(:,1);

daily_VPDALQ = dailyaverage2(USALQ.VPD);
dailyVPDALQ = daily_VPDALQ(:,1);
%% Index for the same dates
precip = woodruffprecipS3.Precipitationin*2.54;

Losdischarge = USGSLos.discharge_cfs*0.028316847;% Convert USLos discharge from cfs to m3s (site no 5357335)

ALQdischarge = USGSALQ.discharge_cfs*0.028316847;% Convert USALQ discharge from cfs to m3s-1 (site no 5357205)
%% Plotting
fig = figure();
right_color = [0.6510    0.6510    0.6510];
left_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

subplot(5,1,1)
yyaxis left
bar(dates,precip)% Convert from in to cm
ylabel('Precip (cm)')
yyaxis right
plot(datesLos,dailyWTDLos,'-','LineWidth',2)
ylabel('WTD (m)')
set(gca,'FontSize',17)
xlim([dates(1) dates(end)])

left_color = [0.3020 0.7451 0.9333];
right_color = [0.8784    0.4588    0.7176];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

subplot(5,1,2)
yyaxis left
% US-Los
plot(dates,Losdischarge,'^')
ylabel('Q (m^3s^{-1})')
yyaxis right
% US-ALQ
plot(dates,ALQdischarge,'o')
ylabel('Q (m^3s^{-1})')
xlim([dates(1) dates(end)])
set(gca,'FontSize',17)
legend('US-Los','US-ALQ','Location','Northwest');

right_color = [0.9294    0.2235    0.5059];
left_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

subplot(5,1,3)
yyaxis left
hold on
plot(dates,dailyTaALQ,'-','Color',[0.8784    0.4588    0.7176],'LineWidth',2)
plot(datesLos,dailyTaLos,'-','Color',[0.3020 0.7451 0.9333],'LineWidth',2)
xlim([dates(1) dates(end)])
hold off
ylabel('T_{air} (^{\circ}C)')
box on
set(gca,'FontSize',17)
yyaxis right
plot(dates,USGSALQ.Twater_mean,'-','Color',[0.9294    0.2235    0.5059],'LineWidth',2)
ylabel('T_{water} (^{\circ}C)')
lgnd = legend('US-ALQ','US-Los','US-ALQ','Location','Northwest','NumColumns',2)
set(lgnd,'color','none');

subplot(5,1,4)
hold on
plot(dates,dailyGPPALQ,'+','Color',[0.8784    0.4588    0.7176])
xlim([dates(1) dates(end)])
plot(datesLos,dailyGPPLos,'s','Color',[0.3020 0.7451 0.9333])
hold off
ylabel('GPP')
box on
set(gca,'FontSize',17)
lgnd = legend('US-ALQ','US-Los','Location','Northwest')
set(lgnd,'color','none');

A = repmat(aveFCH4,2,1);

subplot(5,1,5)
hold on
plot(datesLos,dailyCH4Los,'-','Color',[0.3020 0.7451 0.9333],'LineWidth',2)
plot(dates,A(1:length(dates)),'Color', [  0.5333    0.8039    0.9216 0.7],'LineWidth',2)
plot(dates,dailyFCH4,'-','Color',[0.8784    0.4588    0.7176],'LineWidth',2)
xlim([dates(1) dates(end)])
hold off
ylabel('FCH_4')% (\etamol CH_4 m^-^2s^-^1)
lgnd = legend('US-Los','US-Los 2014-18','US-ALQ','Location','Northwest');
box on
set(gca,'FontSize',17)
set(lgnd,'color','none');
%% Statistics
nanvar(USLos.WTD(year(USLos.TIMESTAMP_END)==2019))
nanvar(USLos.WTD(year(USLos.TIMESTAMP_END)==2020))

nanvar(Losdischarge(year(datesLos)==2019))
nanvar(Losdischarge(year(datesLos)==2020))

nanvar(ALQdischarge(year(dates)==2019))
nanvar(ALQdischarge(year(dates)==2020))

range = [min(Losdischarge(year(datesLos)==2019)) max(Losdischarge(year(datesLos)==2019))]
range = [min(Losdischarge(year(datesLos)==2020)) max(Losdischarge(year(datesLos)==2020))]
range = [min(ALQdischarge(year(datesLos)==2019)) max(ALQdischarge(year(datesLos)==2019))]
range = [min(ALQdischarge(year(datesLos)==2020)) max(ALQdischarge(year(datesLos)==2020))]

C = cov(precip,dailyWTDLos,'omitrows')
C = cov(precip,USGSALQ.Twater_mean,'omitrows')
C = cov(precip,dailyTaLos,'omitrows')
C = cov(precip,dailyTaALQ,'omitrows')
C = cov(precip,Losdischarge,'omitrows')
C = cov(precip,ALQdischarge,'omitrows')

[h,p] = ttest2(dailyRg,dailyRgALQ)
[h,p] = ttest2(Losdischarge,ALQdischarge)
[h,p] = ttest2(dailyTaLos,dailyTaALQ)
[h,p] = ttest2(dailyVPD,dailyVPDALQ)

[R,P] = corrcoef(Losdischarge(~isnan(Losdischarge)&~isnan(ALQdischarge)),ALQdischarge(~isnan(Losdischarge)&~isnan(ALQdischarge)))
[R,P] = corrcoef(ALQdischarge(~isnan(ALQdischarge)&~isnan(dailyWTDLos)),dailyWTDLos(~isnan(ALQdischarge)&~isnan(dailyWTDLos)))
[R,P] = corrcoef(dailyTaLos(~isnan(dailyTaLos)&~isnan(dailyTaALQ)),dailyTaALQ(~isnan(dailyTaLos)&~isnan(dailyTaALQ)))
[R,P] = corrcoef(dailyVPD(~isnan(dailyVPD)&~isnan(dailyVPDALQ)),dailyVPDALQ(~isnan(dailyVPD)&~isnan(dailyVPDALQ)))

nanmean(dailyCH4Los)
nanmean(dailyFCH4)
nanmean(dailyCH4Los(year(datesLos)==2019))
nanmean(dailyCH4Los(year(datesLos)==2020))
nanmean(dailyFCH4(year(dates)==2019))
nanmean(dailyFCH4(year(dates)==2020))

max(gramconvertnmolday(dailyCH4Los))
max(gramconvertnmolday(dailyFCH4))
min(gramconvertnmolday(dailyCH4Los))
min(gramconvertnmolday(dailyFCH4))