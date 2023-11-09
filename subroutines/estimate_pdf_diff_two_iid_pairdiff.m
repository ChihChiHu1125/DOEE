function [f_x1_x2] = estimate_pdf_diff_two_iid_pairdiff(x,min,max,dx,option,n_rnd)
% Given X1,X2~i.i.d~f_X, where f_X is given by finite samples, estimate the
% pdf of f_{X1-X2}
% input::
% x           : input samples (draw from f_X)
% min, max, dx: parameters to determine the x-axis of the histogram
% option      : option == 1 calculate all the pairwise differences
%               option == 2 randomly select n_rnd differences
%               n_rnd only used when option == 2
% output::
% f_x1_x2     : the pdf value of f_{x1-x2}
% 2022/01/06

np = length(x);        % number of available ensemble members

% parameter for the histogram (pdf)
len      = (max-min)/dx + 1;     % the x-axis length of the histogram
half_len = (max-min)/(2*dx);     % half length of the histogram
center   = min:dx:max;             % the positions the pdf is evaluated
edges    = min-dx/2:dx:max+dx/2;   % the two edges of the positions

% option 1 for calculating the differences=================================
% calculate the pairwise differences between samples
if option == 1  
tmp_eta = 0;
ct = 1;
for j=1:np
    for k=1:np
        if k>j
            tmp_eta(ct) = x(j) - x(k);
            ct = ct+1;
        end
    end
end
tmp_eta(ct:2*(ct-1)) = -tmp_eta(1:ct-1);
end             
%==========================================================================

% option 2 for calculating the differences=================================
% randomly pick 100 differences
if option ==2
    
if ~ exist('n_rnd','var')
    n_rnd = 100; % default setting
end

for i=1:n_rnd
    repeat = 1;
    while repeat == 1
        pick1 = randi(np);
        pick2 = randi(np);
        if pick1==pick2
            repeat = 1;
        else
            tmp_eta(i) = x(pick1) - x(pick2);
            repeat = 0;
        end
    end
end
tmp_eta(n_rnd+1:2*n_rnd) = -tmp_eta(1:n_rnd);
end
%==========================================================================

% histogram of these differences:
ct = histcounts(tmp_eta,edges,'Normalization','probability');
f_x1_x2 = ct/dx; % pdf f_X estimated by finite samples

end