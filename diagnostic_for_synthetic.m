% diagnostics of the results from 
% estimate_obs_error_main.m (only for synthetic data)
% 2022/01/16

%% Step 4 - estimate the observation error f_o, given f_d and f_x1_x2

% sensitivity test (determine an optimal alpha first):
%{
list_var_cov_mat = [ 1e-2 2.5e-2 5e-2 1e-1 2.5e-1 5e-1 ...
                     1e0  2.5e0  5e0  1e1  2.5e1  5e1  ...
                     1e2  2.5e2  5e2  1e3  2.5e3  5e3  ...
                     1e4  2.5e4  5e4  1e5  2.5e5  5e5  ...
                     1e6  2.5e6  5e6  1e7  2.5e7  5e7];
ntest = length(list_var_cov_mat);
error = zeros(1,ntest);
stnd  = zeros(1,ntest);
meann = zeros(1,ntest);
foiqr = zeros(1,ntest);
noise = zeros(1,ntest);

cal_iqr = 1; % whether you want to calcaulte iqr

for i=1:ntest
var_cov_mat = list_var_cov_mat(i);
[f_o, ax_b_error] = estimate_obs_error_subroutine(f_d', f_x1_x2', bmin, bmax, dx, var_cov_mat);
noise(i) = sqrt( mean(((f_o(2:end) - f_o(1:end-1))/dx).^2 ) );
error(i) = ax_b_error;
meann(i)  = sum(center.*f_o'*dx);
second_moment = sum(center.^2 .*f_o'*dx);
stnd(i) = sqrt(second_moment - meann(i)^2);

if cal_iqr == 1
    first  = 0.25*sum(f_o);
    second = 0.75*sum(f_o);
    
    for j=1:len-1
        low_val  = sum(f_o(1:j));
        high_val = sum(f_o(1:j+1));
        if ((low_val < first)&&(first < high_val))
            first_quantile = 0.5*(center(j) + center(j+1));
        elseif ((low_val < second)&&(second < high_val))
            second_quantile = 0.5*(center(j) + center(j+1));
            foiqr(i) = second_quantile - first_quantile;
            break
        end
    end
end

end


figure;
yyaxis left
semilogx(list_var_cov_mat, error,'linewidth',1.5)
axis([1e-2 1e7 0 max(error)*1.1])
ylabel('$\| \mathbf{A} f_{\varepsilon^{o}} - f_{D} \|$','interpreter','latex','fontsize',14)

yyaxis right
% plot(list_var_cov_mat, stnd,'linewidth',1.5)
plot(list_var_cov_mat, 50*noise,'linewidth',1.5)
hold on
plot(list_var_cov_mat, foiqr,'-.','linewidth',1.5)
hold off
axis([1e-2 1e7 0 max([stnd, foiqr])*1.1])
ylabel('noise level (solid) and IQR (dashed) of $f_{\epsilon_{o}}$','interpreter','latex','fontsize',14)
grid on

set(gca,'fontsize',11)
xlabel('$\alpha$','interpreter','latex','fontsize',14)
title(['hist = [',num2str(bmin),':',num2str(dx),':',num2str(bmax),'], sample size = ',num2str(num_data)],'interpreter','latex','fontsize',13)

%}

% plot result

var_cov_mat = 500;
[f_o, ax_b_error,inv_cov_mat] = estimate_obs_error_subroutine(f_d', f_x1_x2', bmin, bmax, dx, var_cov_mat);

% calculate the mean, variance, skewness of the distribution
first_moment  = sum(center.*f_o'*dx);
second_moment = sum(center.^2 .*f_o'*dx);

f_o_mean = first_moment; 
f_o_std  = sqrt(second_moment - f_o_mean^2);
f_o_mode = center(find(f_o==max(f_o)));

mid_val = 0.5*sum(f_o);
for i=1:len
    low_val  = sum(f_o(1:i));
    high_val = sum(f_o(1:i+1));
    if ((low_val < mid_val)&&(mid_val < high_val))
        f_o_med = 0.5*(center(i) + center(i+1));
        break
    end
end

% histogram for the true observation errors
ct = histcounts(eps,edges,'Normalization','probability');
f_o_truth = ct/dx;

figure;
% plot(center, f_o,'linewidth',1.5)
plot(center, f_o,'linewidth',2.5)
% semilogy(center, f_o,'linewidth',2.5)

hold on
% plot(center, f_o_truth,'k:','linewidth',1.5)
plot(center, f_o_truth,'k-.','linewidth',2.5)
% semilogy(center, f_o_truth,'k-.','linewidth',2.5)

hold off
legend('estimated','truth','fontsize',20)
% axis([-8 8 0 0.4])
axis([-10 10 0 0.3])
% axis([-12 12 0 0.4])
% axis([-3 20 0 0.3])

grid on
set(gca,'fontsize',16)
ylabel('pdf','interpreter','latex','fontsize',20)
% ylabel('$ \mathbf{f}_{\varepsilon^{o}}$','interpreter','latex','fontsize',16)
% title(['(ens negative bias) $ \mathbf{f}_{\varepsilon^{o}}$ hist = [',num2str(bmin),':',num2str(dx),':',num2str(bmax),'] $\alpha = $',num2str(var_cov_mat),' size = ',num2str(num_data)], ...
%        'interpreter','latex','fontsize',14)
% title(['$ \mathbf{f}_{\varepsilon^{o}}$ hist = [',num2str(bmin),':',num2str(dx),':',num2str(bmax),'] $\alpha = $',num2str(var_cov_mat),' size = ',num2str(num_data)], ...
%        'interpreter','latex','fontsize',14)
% title(['$ \mathbf{f}_{\varepsilon^{o}}$  ($\alpha = $',num2str(var_cov_mat),')'],'interpreter','latex','fontsize',22)
% title(['$ \mathbf{f}_{\varepsilon^{o}}$ :  $ N(2,2^2) $'],'interpreter','latex','fontsize',22)
% title(['$ \mathbf{f}_{\varepsilon^{o}}$ :  $ N(-2,2^2) $'],'interpreter','latex','fontsize',22)
% title(['$ \mathbf{f}_{\varepsilon^{o}}$ (bi-modal)'],'interpreter','latex','fontsize',22)
% title(['$ \mathbf{f}_{\varepsilon^{o}}$ (skewed)'],'interpreter','latex','fontsize',22)
% title(['$ \mathbf{f}_{\varepsilon^{o}}$ (under-dispersive ensemble)'],'interpreter','latex','fontsize',22)
% title(['$ \mathbf{f}_{\varepsilon^{o}}$ (ens negative bias)'],'interpreter','latex','fontsize',22)

% 2023/10/10: to address the comments from reviewer
% title('correlated penalty ($r_{corr} = 3$)','interpreter','latex')
title('[-1 1] derivative','interpreter','latex')

% title(['$ \mathbf{f}_{\varepsilon^{o}}$ [pick mem 1] ($N_e =$',num2str(Ne),' $N_s=$',num2str(num_data),')'], ...
%        'interpreter','latex','fontsize',20)
% title(['$ \mathbf{f}_{\varepsilon^{o}}$ [pick mem010] ($N_e =$',num2str(Ne),' $N_s=$',num2str(num_data),')'], ...
%        'interpreter','latex','fontsize',14)
% title(['$ \mathbf{f}_{\varepsilon^{o}}$ [random pick] ($N_e =$',num2str(Ne),' $N_s=$',num2str(num_data),')'], ...
%        'interpreter','latex','fontsize',20)
% title(['$ \mathbf{f}_{\varepsilon^{o}}$ [ens difference] ($N_e =$',num2str(Ne),' $N_s=$',num2str(num_data),')'], ...
%        'interpreter','latex','fontsize',20)
%}

% plot observation error vs f_x1_x2 vs innovation
%{
var_cov_mat = 10;
[f_o, ax_b_error] = estimate_obs_error_subroutine(f_d', f_x1_x2', bmin,bmax,dx, var_cov_mat);

% extended pdf of f_x1_x2:
f_x1_x2_extended = zeros(1,4*half_len+1);
f_x1_x2_extended(half_len+1:half_len+len) = f_x1_x2;

% generate the matrix for background error
eta_matrix = zeros(len, len);
for i=1:len
    eta_matrix(:,i) = f_x1_x2_extended(len-i+1:1:2*len-i);
end

Ax = eta_matrix*f_o;
Ax_norm = Ax*dx;
% figure;
plot(center, f_o,'color',[0 114 189]/255,'linewidth',2.5)
hold on
plot(center, f_x1_x2,'color',[255 130 0]/255,'linewidth',2.5)
plot(center, f_d,'color',[204 0 0]/255,'linewidth',2.5)
plot(center, Ax_norm,'--','color',[130 0 0]/255,'linewidth',2.5)

hold off
% legend('observation error','ensemble difference','true innovation','estimated innovation', ...
%        'fontsize',17,'location','northwest')

grid on
set(gca,'fontsize',16)
% axis([-15 15 0 0.4])
axis([-15 15 0 0.3])
% axis([-20 20 0 0.3])

ylabel('pdf','interpreter','latex','fontsize',20)

% title(['hist = [',num2str(bmin),':',num2str(dx),':',num2str(bmax),']  $\alpha = $',num2str(var_cov_mat),' size = ',num2str(num_data)],...
%       'interpreter','latex','fontsize',14)
% title(['$\alpha = $',num2str(var_cov_mat)],'interpreter','latex','fontsize',22)
% title(['$ \mathbf{f}_{\varepsilon^{o}}$ (bi-modal)'],'interpreter','latex','fontsize',22)
% title(['$ \mathbf{f}_{\varepsilon^{o}}$ (skewed)'],'interpreter','latex','fontsize',22)
% title(['$ \mathbf{f}_{\varepsilon^{o}}$ (bi-modal)'],'interpreter','latex','fontsize',22)
% title(['$ \mathbf{f}_{\varepsilon^{o}}$ :  $ N(2,2^2) $'],'interpreter','latex','fontsize',22)
title(['$ \mathbf{f}_{\varepsilon^{o}}$ :  $ N(-2,2^2) $'],'interpreter','latex','fontsize',22)

%}

% figure for the graphical abstract
%{
var_cov_mat = 100;
[f_o, ax_b_error] = estimate_obs_error_subroutine(f_d', f_x1_x2', bmin,bmax,dx, var_cov_mat);

ct = histcounts(eps,edges,'Normalization','probability');
f_o_truth = ct/dx;

% extended pdf of f_x1_x2:
f_x1_x2_extended = zeros(1,4*half_len+1);
f_x1_x2_extended(half_len+1:half_len+len) = f_x1_x2;

% generate the matrix for background error
eta_matrix = zeros(len, len);
for i=1:len
    eta_matrix(:,i) = f_x1_x2_extended(len-i+1:1:2*len-i);
end

Ax = eta_matrix*f_o;
Ax_norm = Ax*dx;
% figure;
plot(center, f_o_truth,'k-.','linewidth',2.5)
hold on
plot(center, f_o,'color',[0 114 189]/255,'linewidth',2.5)
plot(center, f_x1_x2,'color',[255 130 0]/255,'linewidth',2.5)
plot(center, f_d,'color',[204 0 0]/255,'linewidth',2.5)

hold off
legend('true obs error','estimated obs error','ensemble difference','innovation', ...
       'fontsize',17,'location','northwest')

grid on
set(gca,'fontsize',16)
axis([-20 20 0 0.36])

ylabel('pdf','interpreter','latex','fontsize',20)
%}

%% Step 5 - sensitivity test (test the validity of randomly choosing one member as the truth)
%{
figure;
plot(center, f_x1_x2,'linewidth',2.5)
hold on
% plot(center, f_x1_x2_pair,'linewidth',1.5)
plot(center, f_truth,'k:','linewidth',2.5)
hold off
% legend('estimate (random pick)','estimate (all pairs)','truth','fontsize',13)
legend('estimate','truth','fontsize',20)
axis([-15 15 0 0.2])

grid on
set(gca,'fontsize',16)
ylabel('pdf','interpreter','latex','fontsize',20)
% title(['$f_{H(X^{truth})-H(X^{b})}$ [pick mem 1] ($N_e =$',num2str(Ne),' $N_s=$',num2str(num_data),')'],...
%      'interpreter','latex','fontsize',17.5)
% title(['$f_{H(X^{truth})-H(X)}$ [pick mem010] ($N_e =$',num2str(Ne),' $N_s=$',num2str(num_data),')'],...
%      'interpreter','latex','fontsize',14)
% title(['$f_{H(X^{truth})-H(X^{b})}$ [random pick] ($N_e =$',num2str(Ne),' $N_s=$',num2str(num_data),')'],...
%      'interpreter','latex','fontsize',17.5)
title(['$f_{H(X^{truth})-H(X^{b})}$ [ens difference] ($N_e =$',num2str(Ne),' $N_s=$',num2str(num_data),')'],...
     'interpreter','latex','fontsize',17.5)
%}