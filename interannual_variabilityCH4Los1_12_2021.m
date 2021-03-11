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
hhCH4 = nanmean(USLos.FCH4_RF_filled)
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
%% Heuristic Model
for i = 1:length(dailyCH4)
    if isnan(dailyCH4(i))||isnan(dailyTairf(i))
        dailyTairf(i) = NaN;
        dailyCH4(i) = NaN;
    end
end

x = dailyTairf(~isnan(dailyTairf));
y = dailyCH4(~isnan(dailyCH4));
modeldates = USLos.TIMESTAMP_END(daily_Tairf(:,2)*48);
modeldates = modeldates(~isnan(dailyTairf)&~isnan(dailyCH4));

modelfun = @(b,x)(b(1).*exp((b(2).*x))); % This is eqn 3 in Rinne et al.
beta0 = [.1;.1];
[beta,r,J,COVB,mse] = nlinfit(x,y,modelfun,beta0); % perform the nonlinear regression of FCH4 and Tair
%% Normalizing daily FCH4 to remove dominant temp effect on emission
% F = (beta(1)+(beta(2).*x));% Utilize this line if there was a linear fit
F = (beta(1).*(exp(beta(2).*x)));% nonlinear best fit
Fn =  y./(F);% Fn represents temp-normalized FCH4

%% Analyze relationship between temp-normalized FCH4 and dis w linear eqn
T=(-5:1:25);

SSE = sum((y-F).^2);% sum of squared estimate of errors
SST = sum((y-mean(y)).^2);% total sum of squares
R2T = 1-SSE./SST;% R-squared
[R,P]= corrcoef(y,beta(1).*exp(beta(2).*x))

figure()
subplot(2,2,[1,2])
plot(x,y,'o','MarkerFaceColor',[0.3020 0.7451 0.9333],'MarkerEdgeColor',[0.3020 0.7451 0.9333])
hold on
% plot(T,beta(1)+(beta(2).*T),'-k','LineWidth',3)% use this line for linear
% fit
plot(T,beta(1).*exp(beta(2).*T),'-','Color',[0.0980    0.3176    0.4118],'LineWidth',3)% exponential fit
hold off
ylabel('FCH_4 [\etamol m^{-2} s^{-1}]')
title(['Discharge and T_{air} versus FCH_4 US-Los 2019-2020'])% make sure this is the right year or manually insert both years
xlim([-5 25])
xlabel('T_{air} {\circ}(C)')% Daily average air temp
formatSpec = '%.2f';
text(-3,30,[{'R^2 = '} num2str(R2T,formatSpec)],'FontSize',17)
set(gca,'FontSize',17)

x= linspace(0,10);

subplot(2,2,3)
Y = round(10*dis.discharge((yr2==2019)&~isnan(dis.discharge)));
Fn19 = accumarray(Y,Fn((yr2==2019)&(~isnan(dis.discharge))),[],@nanmean);
% make sure not to include NaNs in polyfit
plot(unique(Y)/10,Fn19(unique(Y)),'o','MarkerFaceColor',[0.3020 0.7451 0.9333],...
    'MarkerEdgeColor',[0.3020 0.7451 0.9333])
hold on
p=polyfit(unique(Y)/10,Fn19(unique(Y)),1);
f = polyval(p,x);
% plot(x,f,'-','Color',[0.0980    0.3176    0.4118],'LineWidth',3)% Linear fit
hold off
xlabel('Discharge (^{m^3}/_{s})')
ylabel('T_{air} FCH_4')
newdis = unique(Y)/10;
Fb = polyval(p,newdis); 
SSEb = sum((Fn19(unique(Y))-Fb).^2);% sum of squared estimate of errors
SSTb = sum((Fn19(unique(Y))-nanmean(Fn19(unique(Y)))).^2);% total sum of squares
R2Tb = 1-SSEb./SSTb;% R-squared
formatSpec = '%.2f';
% text(0.5,4,[{'R^2 = '} num2str(R2Tb,formatSpec)],'FontSize',17)
set(gca,'FontSize',17)
ylim([-4 6])
title('2019')
[R,P] = corrcoef(polyval(p,unique(Y)/10),Fn19(unique(Y)))
R2 = (R).^2

subplot(2,2,4)
% make sure not to include NaNs in polyfit
hold on
Y = round(10*dis.discharge((yr2==2020)));
Y = Y(~isnan(dailyCH4(yr2==2020)));
Fn20 = Fn(year(modeldates)==2020);
Fn20 = accumarray(Y(~isnan(Y)),Fn20(~isnan(Y)),[],@nanmean);
plot(unique(Y(~isnan(Y)))/10,Fn20(unique(Y(~isnan(Y)))),'o','MarkerFaceColor',[0.3020 0.7451 0.9333],...
    'MarkerEdgeColor',[0.3020 0.7451 0.9333])
p=polyfit(unique(Y(~isnan(Y)))/10,Fn20(unique(Y(~isnan(Y)))),1);
f = polyval(p,x); 
plot(x,f,'-','Color',[0.0980    0.3176    0.4118],'LineWidth',3)% Linear fit
hold off
ylim([-4 6])
newdis = unique(Y(~isnan(Y)))/10;
Fb = polyval(p,newdis); 
SSEb = sum((Fn20(unique(Y(~isnan(Y))))-Fb).^2);% sum of squared estimate of errors
SSTb = sum((Fn20(unique(Y(~isnan(Y))))-nanmean(Fn20(unique(Y(~isnan(Y)))))).^2);% total sum of squares
R2Tb = 1-SSEb./SSTb;% R-squared
[R,P] = corrcoef(polyval(p,unique(Y(~isnan(Y)))/10),Fn20(unique(Y(~isnan(Y)))))
formatSpec = '%.2f';
text(0.5,4,[{'R^2 = '} num2str(R2Tb,formatSpec)],'FontSize',17)
set(gca,'FontSize',17)
title('2020')

%% Lag Analysis
Z = [];
dailyGPP19 = dailyGPP(~isnan(dailyTairf)&~isnan(dailyCH4));
dailyGPP19 = dailyGPP19(year(CH4dates)==2019);
Fn19 = Fn(year(CH4dates)==2019);

for k = 1:90
    B = NaN(length(dailyGPP19),1);
    B((k+1):end) = dailyGPP19(1:(end-k));% SHIFT GPP DOWN
    Z(:,k) = B;
end

Z = [dailyGPP19 Z(:,1:89)];% GPP
% Remove any missing data from both datasets
TF = isnan(Z);
TF2 = isnan(Fn19);
location = [];

for i = 1:size(TF,1)
    for j = 1:90
        if (TF(i,j)==1) || (TF2(i)==1)
            location = [location;[i j]];
        end
    end
end

A = repmat(Fn19,1,90);

for i = 1:length(location)
    Z(location(i,1),location(i,2)) = NaN;
    A(location(i,1),location(i,2)) = NaN;
end

% Calculate cross correlations

for i = 1:90
    TaCH4 = A(:,i);
    GPP = Z(:,i);
    [rho(i),pval(i)] = corr(GPP(~isnan(GPP)),TaCH4(~isnan(TaCH4)));
end
%% Repeat Lag Analysis with just GPP and FCH4
Z = [];
dailyCH419 = dailyCH4(~isnan(dailyTairf)&~isnan(dailyCH4));
dailyCH419 = dailyCH419(year(CH4dates)==2019);

for k = 1:90
    B = NaN(length(dailyGPP19),1);
    B((k+1):end) = dailyGPP19(1:(end-k));% shift GPP down
    Z(:,k) = B;
end

Z = [dailyGPP19 Z(:,1:89)];
% Remove any missing data from both datasets
TF = isnan(Z);
TF2 = isnan(dailyCH419);
location = [];

for i = 1:size(TF,1)
    for j = 1:90
        if (TF(i,j)==1) || (TF2(i)==1)
            location = [location;[i j]];
        end
    end
end

A = repmat(dailyCH419,1,90);

for i = 1:length(location)
    Z(location(i,1),location(i,2)) = NaN;
    A(location(i,1),location(i,2)) = NaN;
end

% Calculate cross correlations

for i = 1:90
    CH4 = A(:,i);
    GPP = Z(:,i);
    [rho2(i),pval2(i)] = corr(GPP(~isnan(GPP)),CH4(~isnan(CH4)));
end
%% Rename variables
lag = 0:89;
% GPP and Ta-CH4
w = lag(pval<0.05 & rho>0);
v = rho(pval<0.05 & rho>0);
% GPP and FCH4
Los2019lag = lag(pval2<0.05 & rho2>0);
Los2019lagcorr = rho2(pval2<0.05 & rho2>0);
%% Repeat for 2020
Z = [];
dailyGPP20 = dailyGPP(~isnan(dailyTairf)&~isnan(dailyCH4)&year(CH4dates)==2020);
Fn20 = Fn(~isnan(dailyTairf)&~isnan(dailyCH4)&year(CH4dates)==2020);

for k = 1:90
    B = NaN(length(dailyGPP20),1);
    B((k+1):end) = dailyGPP20(1:(end-k));% SHIFT GPP DOWN
    Z(:,k) = B;
end

Z = [dailyGPP20 Z(:,1:89)];% GPP
% Remove any missing data from both datasets
TF = isnan(Z);
TF2 = isnan(Fn20);
location = [];

for i = 1:size(TF,1)
    for j = 1:90
        if (TF(i,j)==1) || (TF2(i)==1)
            location = [location;[i j]];
        end
    end
end

A = repmat(Fn20,1,90);

for i = 1:length(location)
    Z(location(i,1),location(i,2)) = NaN;
    A(location(i,1),location(i,2)) = NaN;
end

% Calculate cross correlations

for i = 1:90
    TaCH4 = A(:,i);
    GPP = Z(:,i);
    [rho(i),pval(i)] = corr(GPP(~isnan(GPP)),TaCH4(~isnan(TaCH4)));
end
%% Repeat Lag Analysis with just GPP and FCH4
Z = [];
dailyCH420 = dailyCH4(~isnan(dailyTairf)&~isnan(dailyCH4)&year(CH4dates)==2020);

for k = 1:90
    B = NaN(length(dailyGPP20),1);
    B((k+1):end) = dailyGPP20(1:(end-k));% shift GPP down
    Z(:,k) = B;
end

Z = [dailyGPP20 Z(:,1:89)];
% Remove any missing data from both datasets
TF = isnan(Z);
TF2 = isnan(dailyCH420);
location = [];

for i = 1:size(TF,1)
    for j = 1:90
        if (TF(i,j)==1) || (TF2(i)==1)
            location = [location;[i j]];
        end
    end
end

A = repmat(dailyCH420,1,90);

for i = 1:length(location)
    Z(location(i,1),location(i,2)) = NaN;
    A(location(i,1),location(i,2)) = NaN;
end

% Calculate cross correlations

for i = 1:90
    CH4 = A(:,i);
    GPP = Z(:,i);
    [rho2(i),pval2(i)] = corr(GPP(~isnan(GPP)),CH4(~isnan(CH4)));
end
%% Rename variables
lag = 0:89;
% GPP and Ta-CH4
w2 = lag(pval<0.05 & rho>0);
v2 = rho(pval<0.05 & rho>0);
% GPP and FCH4
Los2020lag = lag(pval2<0.05 & rho2>0);
Los2020lagcorr = rho2(pval2<0.05 & rho2>0);
%% Plot
figure()
subplot(1,2,1)
hold on
plot(Los2019lag',Los2019lagcorr','--','Color',[0.3020 0.7451 0.9333],'LineWidth',2)% original
plot(w',v','-','Color',[0.0980    0.3176    0.4118],'LineWidth',2)% normalized
hold off
xlabel('time lag, days')
ylabel('correlation coefficient')
title('GPP lag influence FCH_4 Los 2019')
ylim([0.1 1])
set(gca,'FontSize',17)

subplot(1,2,2)
hold on
plot(Los2020lag',Los2020lagcorr','--','Color',[0.3020 0.7451 0.9333],'LineWidth',2)% original
plot(w2',v2','-','Color',[0.0980    0.3176    0.4118],'LineWidth',2)% normalized
hold off
title('GPP lag influence FCH_4 Los 2020')
legend('original','normalized')
set(gca,'FontSize',17)
ylim([0.1 1])
%%
figure()
[c,lags] = xcorr(dailyGPP19,dailyCH419,70);
stem(lags,c)
xlabel('time lag, days')  
ylabel('correlation')
title('Lag corr. of GPP, CH_4 residuals')

figure()
[c,lags] = xcorr(dailyGPP19,Fn19,70);
stem(lags,c)
xlabel('time lag, days')  
ylabel('correlation')
title('Lag corr. of GPP, CH_4 residuals')

%% Lag Analysis of WTD and FCH4
Z = [];
dailyWTD19 = dailyWTD(~isnan(dailyTairf)&~isnan(dailyCH4)&~isnan(dailyWTD)&~isnan(dailyGPP)&(year(CH4dates)==2019));
dailyGPP19 = dailyGPP(~isnan(dailyTairf)&~isnan(dailyCH4)&~isnan(dailyWTD)&~isnan(dailyGPP)&(year(CH4dates)==2019));
dailyCH419 = dailyCH4(~isnan(dailyTairf)&~isnan(dailyCH4)&~isnan(dailyWTD)&~isnan(dailyGPP)&(year(CH4dates)==2019));
Fn19 = Fn(~isnan(dailyTairf)&~isnan(dailyCH4)&~isnan(dailyWTD)&~isnan(dailyGPP)&(year(CH4dates)==2019));

%For GPP
for k = 1:200
    B = NaN(length(dailyWTD19),1);
    B((k+1):end) = dailyWTD19(1:(end-k));% shift WTD down
    Z(:,k) = B;
end

Z = [dailyWTD19 Z(:,1:200)];
% Remove any missing data from both datasets
TF = isnan(Z);
TF2 = isnan(dailyGPP19);
location = [];

for i = 1:size(TF,1)
    for j = 1:200
        if (TF(i,j)==1) || (TF2(i)==1)
            location = [location;[i j]];
        end
    end
end

A = repmat(dailyGPP19,1,200);

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
%% For CH4
for k = 1:90
    B = NaN(length(dailyWTD19),1);
    B((k+1):end) = dailyWTD19(1:(end-k));% shift WTD down
    Z(:,k) = B;
end

Z = [dailyWTD19 Z(:,1:90)];
% Remove any missing data from both datasets
TF = isnan(Z);
TF2 = isnan(dailyCH419);
location = [];

for i = 1:size(TF,1)
    for j = 1:90
        if (TF(i,j)==1) || (TF2(i)==1)
            location = [location;[i j]];
        end
    end
end

A = repmat(dailyCH419,1,90);

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
%% Lag Analysis of WTD and Ta-FCH4
Z = [];

for k = 1:90
    B = NaN(length(dailyWTD19),1);
    B((k+1):end) = dailyWTD19(1:(end-k));% shift WTD down
    Z(:,k) = B;
end

Z = [dailyWTD19 Z(:,1:89)];
% Remove any missing data from both datasets
TF = isnan(Z);
TF2 = isnan(Fn19);
location = [];

for i = 1:size(TF,1)
    for j = 1:90
        if (TF(i,j)==1) || (TF2(i)==1)
            location = [location;[i j]];
        end
    end
end

A = repmat(Fn19,1,90);

for i = 1:length(location)
    Z(location(i,1),location(i,2)) = NaN;
    A(location(i,1),location(i,2)) = NaN;
end

% Calculate cross correlations

for i = 1:90
    dayWTD = Z(:,i);
    dayCH4 = A(:,i);
    [rho2(i),pval2(i)] = corr(dayWTD(~isnan(dayWTD)),dayCH4(~isnan(dayCH4)));
end
%% Plot 
lag = 1:90;
% test = 1:200;

figure()
subplot(1,2,1)
hold on
plot(lag,rho,'o','MarkerEdgeColor','k','MarkerSize',11)
plot(lag(pval<0.05 & rho>0),rho(pval<0.05 & rho>0),'o','MarkerSize',11,'MarkerFaceColor',[0.3020 0.7451 0.9333],'MarkerEdgeColor','k');
hold off
xlabel('lag (days)')
ylabel('corr. coeff.')
title('Lag Corr of WTD & FCH_4 US-Los 2019')
set(gca,'FontSize',17)
ylim([0 0.3])
xlim([0 100])

subplot(1,2,2)
hold on
plot(lag,rho2,'o','MarkerEdgeColor','k','MarkerSize',11)
plot(lag(pval2<0.05 & rho2>0),rho2(pval2<0.05 & rho2>0),'o','MarkerSize',11,'MarkerFaceColor',[0.3020 0.7451 0.9333],'MarkerEdgeColor','k');
hold off
xlabel('lag (days)')
ylabel('corr. coeff.')
title('Lag Corr of WTD & T_{air}-FCH_4 US-Los 2019')
set(gca,'FontSize',17)
ylim([0 0.3])
xlim([0 100])