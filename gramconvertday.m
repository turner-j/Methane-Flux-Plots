function [gram] = gramconvertday(umol)
%GRAMCONVERTDAY Converts umol/m2/s to units of gC/m2/DAY
%   Works for array of CO2 data that needs to be converted to cumulative
%   NEE.
gram = umol.*(1E-6).*(12).*(86400);

end

