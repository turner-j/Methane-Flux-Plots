function [gram] = gramconvertnmolday(nmol)
%GRAMCONVERT Converts nmol/m2/s CH4 to units of gC/m2/half hour CO2
%   Works for array of CH4 data
gram = nmol.*(1E-9).*(12).*(86400);

end

