function [f_x1_x2] = estimate_pdf_diff_randompick_from_same_cloud(x, x_cloud, y_cloud, min,max,dx,option,random_seed)
% Given X1,X2~i.i.d~f_X, where f_X is given by finite samples, estimate the
% pdf of f_{X1-X2}
% instead of randomly picking a member as the truth, we pick the "truth member" that
% has similar cloud amount compared with the observation
% input::
% x           : input samples (brightness temperature) (draw from f_X)
% x_cloud     : input samples (cloud amount)
% y_cloud     : observation cloud amount
% min, max, dx: parameters to determine the x-axis of the histogram
% option      : =1 use a randomly picked member
%               =2 use average of the pool member
% output::
% f_x1_x2     : the pdf value of f_{x1-x2}
% 2022/06/26

% parameter for the histogram (pdf)
len      = (max-min)/dx + 1;     % the x-axis length of the histogram
half_len = (max-min)/(2*dx);     % half length of the histogram
center   = min:dx:max;             % the positions the pdf is evaluated
edges    = min-dx/2:dx:max+dx/2;   % the two edges of the positions
rng(random_seed)

np = length(x);        % number of available ensemble members

% we pick the truth member that has similar cloud amount compared to the
% observation cloud amount

% continuous method
%{
pick_cri = 0.05; % only pick the member that has cloud amount [y_cloud-pick_cri, y_cloud+pick_cri]

% which ensemble member has cloud fall into [y_cloud-pick_cri,
% y_cloud+pick_cri]:
select_pool = find( (x_cloud >= y_cloud - pick_cri) + (x_cloud <= y_cloud + pick_cri) ==2 ); 
pool_size   = numel(select_pool);   
%}

% block method
cat     = [0:0.1:1]; % define the categories

tmp_cat = find(y_cloud-cat>0,1,'last');
cloud_lower_limit = cat(tmp_cat);
cloud_upper_limit = cat(tmp_cat+1);

select_pool = find( (x_cloud>=cloud_lower_limit)+(x_cloud<=cloud_upper_limit) ==2 );
pool_size   = numel(select_pool);   

%=========================================================================
if option==1 % only choose one member to be the truth:

% random pick one member in the select pool:
truth_mem = select_pool(randi(pool_size)); % truth member

% tmp_eta = zeros(1,2*(np-1));
tmp_eta = zeros(1,np-1);
tmp_ctr = 1; % counter
for i=1:np
    if i~=truth_mem
        tmp_eta(tmp_ctr) = x(truth_mem)-x(i);
        tmp_ctr = tmp_ctr + 1;
    end
end

% histogram of these differences:
ct = histcounts(tmp_eta,edges,'Normalization','probability');
ct_sym = zeros(1,len);

ct_sym(1:half_len)        = (ct(1:half_len) + ct(end:-1:half_len+2))/2;
ct_sym(end:-1:half_len+2) = (ct(1:half_len) + ct(end:-1:half_len+2))/2;
ct_sym(half_len+1)        = ct(half_len+1);

f_x1_x2 = ct/dx; % pdf f_X estimated by finite samples

%=========================================================================
elseif option == 2 % average over all possible truth member:

f_x1_x2 = zeros(1,len);

for mem = 1:pool_size
truth_mem = select_pool(mem); % truth member

% tmp_eta = zeros(1,2*(np-1));
tmp_eta = zeros(1,np-1);
tmp_ctr = 1; % counter
for i=1:np
    if i~=truth_mem
        tmp_eta(tmp_ctr) = x(truth_mem)-x(i);
        tmp_ctr = tmp_ctr + 1;
    end
end

% histogram of these differences:
ct = histcounts(tmp_eta,edges,'Normalization','probability');
f_x1_x2 = f_x1_x2 + ct/dx/pool_size; % pdf f_X estimated by finite samples
end

end