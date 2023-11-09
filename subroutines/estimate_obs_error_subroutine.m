function [f_o, ax_b_error,D]=estimate_obs_error_subroutine(f_d, f_x1_x2, min,max,dx, var_cov_mat)
% input:
% f_d    : pdf of the innovation
% f_x1_x2: pdf of the differences between the ensemble members
% min,max,dx : the parameters to define the pdf of f_d and f_x1_x2
% var_cov_mat: the variance of the "prior" in the cost-function when 
%              solving the minimization problem (or the strength of the penalty term)
%              larger value means less constraint of the smoothness of the resulting pdf
%
% You can plot the result using the below 
% plot(x_output, y_output)
% revised - 2022/01/10

% parameters:
len      = (max-min)/dx + 1;     % the x-axis length of the histogram
half_len = (max-min)/(2*dx);     % half length of the histogram
center   = min:dx:max;           % the positions the pdf is evaluated
edges    = min-dx/2:dx:max+dx/2; % the two edges of the positions

% extended pdf of f_x1_x2:
f_x1_x2_extended = zeros(1,4*half_len+1);
f_x1_x2_extended(half_len+1:half_len+len) = f_x1_x2;

% generate the matrix for background error
eta_matrix = zeros(len, len);
for i=1:len
    eta_matrix(:,i) = f_x1_x2_extended(len-i+1:1:2*len-i);
end

%% turn the problem into an optimization problem (see matlab documentation for detailed usage of quadprog)

% constructing the prior constrain on the derivative:
% the "penalty term" is formulated as: (Dx)^t (S)^{-1} (Dx)

% derivative matrix (size: n-1*n)

% [-1 1] derivative
tmp_mat = diag(ones(len-1,1),1) - diag(ones(len,1));
D       = tmp_mat(1:len-1,1:len);

% [-1 0 1] derivative
% tmp_mat = diag(ones(len-2,1),2) - diag(ones(len,1));
% D       = tmp_mat(1:len-2,1:len);

% "error covariance matrix" S for first-order derivative 
% (note: dimension = (n-1)*(n-1))
% if for [-1 0 1] derivative, dimension = (n-2)*(n-2)

dim     = len-1;
% dim     = len-2; 
cov_mat = zeros(dim);
r_influ = 2; % tuning parameter

for i=1:4*r_influ
    cov_mat = cov_mat + exp(-i^2/r_influ^2)* ( diag(ones(dim-i,1),i) + diag(ones(dim-i,1),-i));
end

cov_mat     = cov_mat + diag(ones(dim,1)); % diagonal element
cov_mat     = var_cov_mat* cov_mat;
% cond_num    = -8;
% inv_cov_mat = inv_SVD(cov_mat, cond_num, 1);
inv_cov_mat = inv(cov_mat);

% the cost-function (see MATLAB documentation for definitions)
H = 2*( (dx^2)*(eta_matrix')*eta_matrix + D'*inv_cov_mat*D); % the costfunction in the minimization problem
% H = 2*( (dx^2)*(eta_matrix')*eta_matrix ); % the costfunction in the minimization problem
f = -2*dx*(eta_matrix')*f_d;
A = -eye(len);
b = zeros(len,1);
Aeq = ones(1,len);
beq = 1/dx;
% x0  = zeros(len,1);
x0  = ones(len,1)/dx;
options = optimset('Display', 'off','LargeScale', 'off','MaxIter',1000); % 1000 maximum # of iterations

f_o = quadprog(H,f,A,b,Aeq,beq,[],[],x0,options);
% f_o = quadprog(H,f,A,b,[],[],[],[],x0,options);
% f_o = quadprog(H,f,[],[],Aeq,beq,[],[],x0,options);

% verification:
predict = eta_matrix*f_o;
predict = predict/sum(predict)/dx;
ax_b_error = norm(predict-f_d);


end