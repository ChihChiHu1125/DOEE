% generate synthetic data (for idealized experiment)
% 2022/08/29

rng('default') % random seed

% parameters:
num_data = 10000;   % number of observations
Ne       = 100;   % number of ensemble members

% definitions of the main variables:
truth   = zeros(num_data,1);  % true state
obs     = zeros(num_data,1);  % observation
ens     = zeros(num_data,Ne); % the ensemble member

% errors
eps     = zeros(num_data,1);  % observation error
eps_est = zeros(num_data,1);  % a crude way to estimate of observation error (observation - pseudo truth)
eta     = zeros(num_data,Ne); % real background error (for each ensemble member)
eta_est = zeros(num_data,Ne); % estimated background error
ens_per = zeros(num_data,Ne); % ensemble perturbation

% generate synthetic data:
for t=1:num_data
    %% part I: the distribution of the ensemble/truth:
  % draw the truth and the ensemble from the same distribution
    
    % non-symmetric distribution
    
    % gamma distribution 
%     shape = 2; % shape parameter
    shape = 2+(rand(1)-0.5); % shape parameter

%     scale = 2; % scale parameter
    scale = 2+(rand(1)-0.5); % scale parameter

    draw  = gamrnd(shape,scale,Ne+1,1);
    ens(t,:) = draw(1:Ne);
%     ens(t,:) = draw(1:Ne)+1.5*randn(Ne,1); % over-dispersive ensemble
    truth(t) = draw(end);
%     truth(t) = draw(end)+1.5*randn(1,1); % under-dispersive ensemble

    
    % beta distribution B(alpha = 2, beta = 5)
    %{
    alpha = 2;
    beta  = 5;
    draw  = 5*(betarnd(alpha, beta, Ne+1, 1)-0.5);
    ens(t,:)  = draw(1:Ne);
    
%     draw  = 5*(betarnd(alpha, beta, Ne+1, 1)-0.5);
    truth(t) = draw(Ne+1);
    %}
    
    % parameters for skew-normal distribution
    %{
    mu   = 0;
    sig  = 1;
    skew = 1;
    kurt = 4;
    draw = pearsrnd(mu,sig,skew,kurt,[Ne+1,1]);
    xp(t,:)  = draw(1:Ne);
    xtrue(t) = draw(Ne+1);
    %}
    
    % Normal distribution
    %{
    std  = 1;
%     std = 3*abs(randn(1,1));
%     std  = raylrnd(1);
    men = 0;
%    men = randn(1,1);
    draw = std*randn(1,Ne+1) + men; % total Ne+1 draw
%    xp(t,:)  = 1.1*draw(1:Ne);
    ens(t,:)  = draw(1:Ne);
    truth(t)  = draw(Ne+1);
%     xtrue(t) = std*1.2*randn(1,1) + men;
    %}
    
    %% part II: the distribution of the observation error:
    u = rand(1,1);
    
    eps(t) = 3*randn(1,1);

    % biased Gaussian
    %{
%     eps(t) = 2 + 2*randn(1,1); % positive bias
    eps(t) = -2 + 2*randn(1,1); % negative bias
    %}
    
    % bi-Gaussian mixture
    
    if u>1/2
        eps(t) = -4+1*randn(1,1);
    else
        eps(t) =  4+1*randn(1,1);
    end
    %}
    
    % multi-Gaussian mixture
    %{
    if u>2/3
        eps(t) = 2+0.5*randn(1,1);
    elseif u<=2/3 && u>1/3
        eps(t) = 0+0.5*randn(1,1);
    else
        eps(t) = -2+0.5*randn(1,1);
    end
    %}
    
    % student t distribution
    %{
    mu     = 5; % degree of freedom in student t distribution 
    eps(t) = trnd(mu,1);
    %}
    
    % Cauchy Lorenz distribution
    %{
    b = 1; % half width
    eps(t) = b*tan(pi*(rand(1,1)-1/2));
    %}
    
    % mixed distribution with Gaussian on the right, Cauchy on the left
    %{
    b = 2; % half width for the Cauchy
    if u>1/2
        eps(t) = abs(randn(1,1));
    else
        eps(t) = -abs(b*tan(pi*(rand(1,1)-1/2)));
    end
    %}
    
    % Gamma distribution
%     eps(t) = gamrnd(2,2,1,1);

    


    %% generate the observation:
    obs(t) = truth(t) + eps(t);
     
    %% add model error (optional):
    %{
    model_err_mag = 0.5;
%     model_err = model_err_mag*randn(1,1);
%     model_err = 0.5; 
%     xp(t,:) = xp(t,:) + model_err;
    
%     eps_est(t) = y(t)-xp(t,kk); % one crude way to estimate of observation error: observation - "pseudo truth"
    %}
end