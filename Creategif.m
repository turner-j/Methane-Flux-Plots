% Clear workspace and command window
clear all;clc;close all;
% Load data
load('USALQ19smallrange_withCH4.mat')
load('USALQ2020_smallrange_withCH4.mat')
%% ALQ 2019
Table = USALQ19smallrange_withCH4;

for i= 2:width(Table)
    Table.(i)(Table.(i)==-9999) = NaN;
end

Table2 = USALQ2020_smallrange_withCH4;

for i= 3:width(Table2)
    Table2.(i)(Table2.(i)==-9999) = NaN;
end
%% Make the Scatter Wind Roses for LOST CREEK
% Eliminate unrealistic data
WS = [Table.WS;Table2.WS];
WD = [Table.WD;Table2.WD];

WS(WS<0) = NaN;
WD(WD<0) = NaN;

timestamp = [Table.TIMESTAMP_END;Table2.TIMESTAMP_END];
CH4 = [Table.FCH4_RF_filled;Table2.FCH4_RF_filled];
%% Set limits and variable names
limU = [min(WS),max(WS)];

name_U = 'U (m/s)'; % #4 name of variable U
name_IU = 'CH4 (nmol/m2/s)'; % #4 name of variable IU
%% Plot
starting = 0;
% Plot the first frame:
ScatterWindRose(WD((1+(starting*144)):(144*(1+starting))),WS((1+(starting*144)):(144*(1+starting))),'Ylim',limU,'labelY',name_U,'labelZ',name_IU,'Z',CH4((1+(starting*144)):(144*(1+starting))));
title([datestr(timestamp(starting+1),'mmm dd yyyy')])
%% Split up data into even numbers of 3-day sets
remainder = mod(size(WS,1),144);
ending = size(WS,1)-remainder;
repetitions = ending/144;
%% Plot all figures
for starting = 1:(repetitions-1)
    newCH4 = CH4((1+(starting*144)):(144*(1+starting)));
    newWD = WD((1+(starting*144)):(144*(1+starting)));
    newWS = WS((1+(starting*144)):(144*(1+starting)));
    if (sum(isnan(newWD))>130)||(sum(isnan(newWS))>130)
        continue
    end
    fig = figure;
    ScatterWindRose(newWD,newWS,'Ylim',limU,'labelY',name_U,'labelZ',name_IU,'Z',newCH4);
    title([datestr(timestamp(starting*144+1),'mmm dd yyyy')])
    drawnow
    frame = getframe(fig);
    im{starting} = frame2im(frame);
end

close;
%%
filename = 'testAnimated.gif'; % Specify the output file name
for starting = 1:(repetitions-1)
    newCH4 = CH4((1+(starting*144)):(144*(1+starting)));
    newWD = WD((1+(starting*144)):(144*(1+starting)));
    newWS = WS((1+(starting*144)):(144*(1+starting)));
    if (sum(isnan(newWD))>130)||(sum(isnan(newWS))>130)
        continue
    end
    [A,map] = rgb2ind(im{starting},256);
    if starting == 1
        if (sum(isnan(newWD))>130)||(sum(isnan(newWS))>130)
            continue
        end
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.2);
    else
        if (sum(isnan(newWD))>130)||(sum(isnan(newWS))>130)
            continue
        end
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.2);
    end
end