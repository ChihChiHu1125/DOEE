function [ind_pixel,bounds] = sort_by_cloud_amount(cloud_amount, inc, lower_bound, upper_bound);
% sort subroutine based on user-defined cloud amount
% new subroutine - 2022/01/05
% TB: brightness temperature data to be sorted
% cloud amount: cloud amount proxy [0,1] (sometimes can be slightly smaller than 0)
% inc: the inc to discretize the cloud amount. e.g., inc = 0.1, then cloud
% amount will be 0-0.1, 0.1-0.2, 0.2-0.3,...

% make sure if the lower bound, inc are consistent:
if mod(upper_bound-lower_bound, inc)~=0
    error('[Error from SUBROUTINE: sort_by_cloud_amount] Please make sure lower, upper bound, increment are consistent')
end

% parameters
bounds    = lower_bound:inc:upper_bound;
num_bin = length(bounds)-1;

% the output is the indices of pixels that fall in the bin defined above
ind_pixel = cell(num_bin,1);

for i=1:num_bin
    ind_pixel{i} = find( ((cloud_amount>=bounds(i)) + (cloud_amount<bounds(i+1))) ==2 );
end

end