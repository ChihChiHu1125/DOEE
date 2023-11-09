function [A_inv] = inv_SVD(A,cond_num,option);
% compute the inverse of a matrix A using the given condition number
% option = 1: calculate the inverse of A
% option = 2: calculate the square root of A

[u,sig,v] = svd(A);
sig_val   = diag(sig);
sig_inv   = sig; % a matrix for storge

if option == 1
for i=1:length(sig_val)
    if sig_val(i) < (10^(cond_num))*sig_val(1)
        sig_inv(i,i) = 0; % to make the condition number < cond_num (zero out too small singular value)
    else
        sig_inv(i,i) = 1/sig_val(i);
    end
end

% calculate the inverse of C by SVD:
A_inv = u*sig_inv*v'; 
end % end option 1============================

if option ==2
for i=1:length(sig_val)
    if sig_val(i) < (10^(cond_num))*sig_val(1)
        sig_inv(i,i) = 0; % to make the condition number < cond_num (zero out too small singular value)
    else
        sig_inv(i,i) = sqrt(sig_val(i));
    end
end

% calculate the inverse of C by SVD:
A_inv = u*sig_inv*v';
end % end option 2============================


end