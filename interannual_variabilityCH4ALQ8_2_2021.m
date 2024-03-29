%% Clear workspace, clear variables, and clear figures
clear all;close all;clc;
load('USALQ2020_2019_calcresp.mat')
load('woodruffprecip.mat')
load('dischargevars.mat')
%% Replace placeholders for missing values with NaN
USALQ = [USALQHH2019; USALQHH2020]; 

for i= 3:width(USALQ)
    USALQ.(i)(USALQ.(i)==-9999) = NaN;
end
%% Sort data by the half hour for whichever month or months you want
precip = woodruffprecipS3.Precipitationin*2.54;% convert in to cm

dis = USGSALQ;
dis.discharge = (dis.discharge_cfs).*0.028316847;% convert cfs to m3/s
yr2 = year(dis.datetime);
mon2 = month(dis.datetime);

%% Make daily averages
hhCH4 = nanmean(USALQ.FCH4_RF_filled)
mo = month(USALQ.TIMESTAMP_END);
monthlyCH4 = accumarray(mo,USALQ.FCH4_RF_filled,[],@nanmean);

daily_CH4 = dailyaverage2(USALQ.FCH4_RF_filled);
dailyCH4 = daily_CH4(:,1);
CH4dates = USALQ.TIMESTAMP_END(daily_CH4(:,2).*48);

daily_Tairf = dailyaverage2(USALQ.TA_F);
dailyTairf = daily_Tairf(:,1);

daily_GPP = dailyaverage2(USALQ.GPP_F);
dailyGPP = daily_GPP(:,1);

daily_VPD = dailyaverage2(USALQ.VPD);
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
modeldates = USALQ.TIMESTAMP_END(daily_Tairf(:,2)*48);
modeldates = modeldates(~isnan(dailyTairf)&~isnan(dailyCH4));

modelfun = @(b,x)(b(1).*exp((b(2).*x))); % This is eqn 3 in Rinne et al.
beta0 = [.1;.1];
[beta,r,J,COVB,mse] = nlinfit(x,y,modelfun,beta0); % perform the nonlinear regression of FCH4 and Tair
%% Normalizing daily FCH4 to remove dominant temp effect on emission
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
plot(x,y,'o','MarkerFaceColor',[0.8784    0.4588    0.7176],'MarkerEdgeColor',[0.8784    0.4588    0.7176])
hold on
plot(T,beta(1).*exp(beta(2).*T),'-','Color',[ 0.4902    0.2588    0.4000],'LineWidth',3)% exponential fit
hold off
ylabel('FCH_4 [\etamol m^{-2} s^{-1}]')
% title(['Discharge and T_{air} versus FCH_4 US-ALQ 2019-2020'])% make sure this is the right year or manually insert both years
title('US-ALQ')
xlim([-5 25])
xlabel('T_{air} {\circ}(C)')% Daily average air temp
formatSpec = '%.2f';
text(-3,30,[{'R^2 = '} num2str(R2T,formatSpec)],'FontSize',17)
set(gca,'FontSize',17)
xl = xlim;
yl = ylim;
xt = 0.05 * (xl(2)-xl(1)) + xl(1)
yt = 0.90 * (yl(2)-yl(1)) + yl(1)
caption = sprintf('y = %.2f*exp(%.2f*x)', beta(1), beta(2));
text(xt, yt, caption, 'FontSize', 16, 'Color', 'k');

x= linspace(0,.5);

subplot(2,2,3)
Y = round(100*dis.discharge((yr2==2019)&~isnan(dis.discharge)));
Fn19 = accumarray(Y,Fn((yr2==2019)&(~isnan(dis.discharge))),[],@nanmean);
% make sure not to include NaNs in polyfit
plot(unique(Y)/100,Fn19(unique(Y)),'o','MarkerFaceColor',[0.8784    0.4588    0.7176],'MarkerEdgeColor',[0.8784    0.4588    0.7176])
hold on
p=polyfit(unique(Y)/100,Fn19(unique(Y)),1);
f = polyval(p,x);
plot(x,f,'-','Color',[ 0.4902    0.2588    0.4000],'LineWidth',3)% Linear fit
hold off
xlabel('Discharge (^{m^3}/_{s})')
ylabel('T_{air} FCH_4')
newdis = unique(Y)/100;
Fb = polyval(p,newdis); 
SSEb = sum((Fn19(unique(Y))-Fb).^2);% sum of squared estimate of errors
SSTb = sum((Fn19(unique(Y))-nanmean(Fn19(unique(Y)))).^2);% total sum of squares
R2Tb = 1-SSEb./SSTb;% R-squared
formatSpec = '%.2f';
text(0.05,0.5,[{'R^2 = '} num2str(R2Tb,formatSpec)],'FontSize',17)
set(gca,'FontSize',17)
xlim([0.0 0.5])
% title('2019')
[R,P] = corrcoef(polyval(p,unique(Y)/100),Fn19(unique(Y)))
R2 = (R).^2
ylim([-0.5 2])
title('2019')

xl = xlim;
yl = ylim;
xt = 0.05 * (xl(2)-xl(1)) + xl(1)
yt = 0.90 * (yl(2)-yl(1)) + yl(1)
caption = sprintf('y = %.2f*x + %.2f',p(1),p(2));
text(xt, yt, caption, 'FontSize', 16, 'Color', 'k');

subplot(2,2,4)
% make sure not to include NaNs in polyfit
hold on
Y = round(100*dis.discharge((yr2==2020)));
Y = Y(~isnan(dailyCH4(yr2==2020)));
Fn20 = Fn(year(modeldates)==2020);
Fn20 = accumarray(Y(~isnan(Y)),Fn(~isnan(Y)),[],@nanmean);
plot(unique(Y(~isnan(Y)))/100,Fn20(unique(Y(~isnan(Y)))),'o','MarkerFaceColor',[0.8784    0.4588    0.7176],'MarkerEdgeColor',[0.8784    0.4588    0.7176])
hold off
xlim([0.0 0.5])
newdis = unique(Y(~isnan(Y)))/100;
set(gca,'FontSize',17)
% title('2020')
ylim([-0.5 2])
xlabel('Discharge (^{m^3}/_{s})')
ylabel('T_{air} FCH_4')
title('2020')
%% Lag Analysis
Z = [];
dailyGPP19 = dailyGPP(~isnan(dailyTairf)&~isnan(dailyCH4)&year(CH4dates)==2019);
CH4dates19 = CH4dates(year(CH4dates)==2019);

idx = ismember(modeldates,CH4dates19,'rows');
modeldates19 = modeldates(idx);
Fn19 = Fn(idx);

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
dailyGPP19 = dailyGPP(~isnan(dailyTairf)&~isnan(dailyCH4)&year(CH4dates)==2019);
dailyCH419 = dailyCH4(~isnan(dailyTairf)&~isnan(dailyCH4)&year(CH4dates)==2019);

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
ALQ2019lag = lag(pval2<0.05 & rho2>0);
ALQ2019lagcorr = rho2(pval2<0.05 & rho2>0);
%% Repeat for 2020
Z = [];
dailyGPP20 = dailyGPP(~isnan(dailyTairf)&~isnan(dailyCH4)&year(CH4dates)==2020);
% Fn20 = Fn(~isnan(dailyTairf)&~isnan(dailyCH4)&year(CH4dates)==2020);
Fn20 = Fn(year(modeldates)==2020);

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
dailyGPP20 = dailyGPP(~isnan(dailyTairf)&~isnan(dailyCH4)&year(CH4dates)==2020);
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
ALQ2020lag = lag(pval2<0.05 & rho2>0);
ALQ2020lagcorr = rho2(pval2<0.05 & rho2>0);
%% Plot
figure()
subplot(1,2,1)
hold on
plot(ALQ2019lag',ALQ2019lagcorr','--','Color',[0.8784    0.4588    0.7176],'LineWidth',2)% original
plot(w',v','-','Color',[ 0.4902    0.2588    0.4000],'LineWidth',2)% normalized
hold off
xlabel('time lag (days)')
ylabel('correlation coefficient')
% title('GPP lag influence FCH_4 ALQ 2019')
title('US-ALQ 2019')
ylim([0.1 1])
set(gca,'FontSize',17)

subplot(1,2,2)
hold on
plot(ALQ2020lag',ALQ2020lagcorr','--','Color',[0.8784    0.4588    0.7176],'LineWidth',2)% original
plot(w2',v2','-','Color',[ 0.4902    0.2588    0.4000],'LineWidth',2)% normalized
hold off
% title('GPP lag influence FCH_4 ALQ 2020')
title('US-ALQ 2020')
legend('original','normalized')
set(gca,'FontSize',17)
ylim([0.1 1])
xlabel('time lag (days)')
ylabel('correlation coefficient')
%% Water temperature
dis = dis(51:end,:);

figure()
subplot(1,2,1)
hold on

x = dis.Twater_mean(ismember(month(dis.datetime),4:10));
disdates = dis.datetime(ismember(month(dis.datetime),4:10));
y = Fn(ismember(month(modeldates),4:10));
modeldates2 = modeldates(ismember(month(modeldates),4:10));
idx = ismember(modeldates2,disdates,'rows');
x = x(idx);

plot(x,y,'o','Color',[0.8784    0.4588    0.7176])
coefficients = polyfit(x,y,1);
xFit = linspace(0,25,100);
yFit = polyval(coefficients,xFit);
plot(xFit,yFit,'-','Color',[ 0.4902    0.2588    0.4000],'LineWidth',2);
[R,P] = corrcoef(x,y,'Rows','complete')
hold off
xlabel('T_{stream} (^{\circ} C)')
ylabel('T_{air}-FCH_4')
% title('T_{water} vs. T_{air}-FCH_4 US-ALQ Apr-Oct')
xlim([0 25])
ylim([-0.5 4.5])
set(gca,'FontSize',17)
%%
subplot(1,2,2)
hold on
x = dis.Twater_mean(~ismember(month(dis.datetime),4:10));
disdates = dis.datetime(~ismember(month(dis.datetime),4:10));
y = Fn(~ismember(month(modeldates),4:10));
modeldates2 = modeldates(~ismember(month(modeldates),4:10));
idx = ismember(modeldates2,disdates,'rows');
x = x(idx);

plot(x,y,'o','Color',[0.8784    0.4588    0.7176])
coefficients = polyfit(x,y,1);
xFit = linspace(0,25,100);
yFit = polyval(coefficients,xFit);
plot(xFit,yFit,'-','Color',[ 0.4902    0.2588    0.4000],'LineWidth',2);
[R,P] = corrcoef(x,y,'Rows','complete');
hold off
% title('T_{water} vs. T_{air}-FCH_4 US-ALQ, Not Apr-Oct')
xlim([0 5])
ylim([-0.5 4.5])
set(gca,'FontSize',17)
xlabel('T_{stream} (^{\circ} C)')
ylabel('T_{air}-FCH_4')

nanmean(dis.Twater_mean(ismember(year(dis.datetime),2019)))
nanmean(dis.Twater_mean(ismember(year(dis.datetime),2020)))