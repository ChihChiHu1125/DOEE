% Main script for estimating the observation error
% for synthetic data
% run generate_synthetic_data.m first to generate the input for this script
% 2022/08/29

%% Step 3 - obtain the pdf of the "difference between ensemble members" and "innovations"

% parameters for the histogram
bmin = -45;
bmax =  45;
dx  =   0.25;
len      = (bmax-bmin)/dx + 1;     % the x-axis length of the histogram
half_len = (bmax-bmin)/(2*dx);     % half length of the histogram
center   = bmin:dx:bmax;           % the positions the pdf is evaluated
edges    = bmin-dx/2:dx:bmax+dx/2; % the two edges of the positions

% parameters for estimating the pdf of X1-X2:
% if we use the subroutine "estimate_pdf_diff_two_iid_pairdiff"; if not,
% below parameters do not matter!
% option = 1: calculate all pairwise differences
% option = 2: randomly select some pairwise difference
option   = 1;                    % parameters for estimate_pdf_diff_two_iid_pairdiff
n_rnd    = 100;                  % parameters for estimate_pdf_diff_two_iid_pairdiff

f_x1_x2   = zeros(1, len);         % the pdf of "difference between ensemble members"
f_d       = zeros(1, len);         % the pdf of "tb innovation"
f_truth   = zeros(1, len);         % the pdf of "tb innovation"
f_x1_x2_pair = zeros(1, len);

tic
random_seed = 1125; % used if use estimate_pdf_diff_two_iid_randompick method

for i=1:num_data
    tmp     = estimate_pdf_diff_two_iid_pairdiff(ens(i,:),bmin,bmax,dx,1,n_rnd);
%     tmp     = estimate_pdf_diff_two_iid_randompick(ens(i,:),bmin,bmax,dx,random_seed);
%     tmp     = estimate_pdf_diff_randompick_from_same_cloud(ens_input{i},cloud_ens_input{i},cloud_obs_input(i),bmin,bmax,dx,2,random_seed);
    f_x1_x2 = f_x1_x2 + tmp/num_data;

%     tmp     = estimate_pdf_diff_two_iid_pairdiff(ens(i,:),bmin,bmax,dx,1,n_rnd);
%     f_x1_x2_pair = f_x1_x2_pair + tmp/num_data;
    
    tmp_d   = estimate_pdf_innovations(ens(i,:),obs(i),bmin,bmax,dx);
    f_d     = f_d + tmp_d/num_data;

    tmp_d   = estimate_pdf_innovations(ens(i,:),truth(i),bmin,bmax,dx);
    f_truth = f_truth + tmp_d/num_data;
%     tmp_cloud_d = estimate_pdf_innovations(cloud_ens_input{i},cloud_obs_input(i),cmin, cmax,cdx);
%     f_cloud_d = f_cloud_d + tmp_cloud_d/num_data;
    
    if mod(i,1000) == 1
        disp(['complete ', num2str(i/num_data*100),' % ...'])
    end
end
toc
%}