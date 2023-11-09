function [f_x1_x2] = estimate_pdf_diff_two_iid_randompick(x,min,max,dx,random_seed)
% Given X1,X2~i.i.d~f_X, where f_X is given by finite samples, estimate the
% pdf of f_{X1-X2}
% by randomly picking a member as truth and estimate the histogram
% input::
% x           : input samples (draw from f_X)
% min, max, dx: parameters to determine the x-axis of the histogram
% output::
% f_x1_x2     : the pdf value of f_{x1-x2}
% 2022/01/12

rng(random_seed)

np = length(x);        % number of available ensemble members
truth_mem = randi(np); % truth
% truth_mem = 1;

% tmp_eta = zeros(1,2*(np-1));
tmp_eta = zeros(1,np-1);
tmp_ctr = 1; % counter
for i=1:np
    if i~=truth_mem
        tmp_eta(tmp_ctr) = x(truth_mem)-x(i);
        tmp_ctr = tmp_ctr + 1;
    end
end

% tmp_eta(tmp_ctr+1:2*tmp_ctr) = -tmp_eta(1:tmp_ctr); % force the pdf to be symmetric

% parameter for the histogram (pdf)
len      = (max-min)/dx + 1;     % the x-axis length of the histogram
half_len = (max-min)/(2*dx);     % half length of the histogram
center   = min:dx:max;             % the positions the pdf is evaluated
edges    = min-dx/2:dx:max+dx/2;   % the two edges of the positions

% histogram of these differences:
ct = histcounts(tmp_eta,edges,'Normalization','probability');
ct_sym = zeros(1,len);

ct_sym(1:half_len)        = (ct(1:half_len) + ct(end:-1:half_len+2))/2;
ct_sym(end:-1:half_len+2) = (ct(1:half_len) + ct(end:-1:half_len+2))/2;
ct_sym(half_len+1)        = ct(half_len+1);

% f_x1_x2 = ct/dx; % pdf f_X estimated by finite samples
f_x1_x2 = ct_sym/dx;
end