function [f_d] = estimate_pdf_innovations(ens,obs,min,max,dx)
% Given ens(size can vary) and obs (a scalar), estimate the pdf of the
% innovation
% input::
% ens         : input samples (draw from f_X)
% obs         : input observation value
% min, max, dx: parameters to determine the x-axis of the histogram
% output::
% f_x1_x2     : the pdf value of f_{x1-x2}
% 2022/01/07

np = length(ens);        % number of available ensemble members

% parameter for the histogram (pdf)
len      = (max-min)/dx + 1;     % the x-axis length of the histogram
half_len = (max-min)/(2*dx);     % half length of the histogram
center   = min:dx:max;             % the positions the pdf is evaluated
edges    = min-dx/2:dx:max+dx/2;   % the two edges of the positions

% innovation
d = obs - ens;

% histogram of these differences:
ct = histcounts(d,edges,'Normalization','probability');
f_d = ct/dx; % pdf f_X estimated by finite samples

end