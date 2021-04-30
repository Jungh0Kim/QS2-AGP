function [U_val] = gp_U_eval(gp_mean, gp_sigma)
% Learning function 'U' (introduced in AK-MCS, Echard et al., 2011 )
% Jungho Kim
U_val = abs(gp_mean)./gp_sigma;
end % function end
