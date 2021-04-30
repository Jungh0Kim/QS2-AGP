function [f, g] = Bayopt_func(x)
% Learning function for adaptive training of quantile surrogates
% Jungho Kim
global model1 model2 bt comp_id f_opt_sv hyp_optim_o meanfunc...
    covfunc likfunc d_doe Objval  iter

[mu_GPu_g1, s2_GPu_g1] = gpm_m(model1,x);
std_GPu_g1 = sqrt(s2_GPu_g1);
[mu_GPu_g2, s2_GPu_g2] = gpm_m(model2,x);
std_GPu_g2 = sqrt(s2_GPu_g2);

% Find function component that influences the system failure domain
[~,comp_id] = min([abs(mu_GPu_g1./std_GPu_g1) abs(mu_GPu_g2./std_GPu_g2)],[],2);

del_tol = 5e-4;
if iter==1
    penalty_func = 1;
elseif iter > 1
    % Penalty function
    dist_f = abs(f_opt_sv(end-1,:) - f_opt_sv(end,1))./abs(f_opt_sv(end,1));
    if dist_f > del_tol
        [mu_o, ~] = gp(hyp_optim_o, @infGaussLik, meanfunc, covfunc, likfunc, d_doe, Objval, x);
        penalty_func = abs(mu_o - f_opt_sv(end,1))./abs(f_opt_sv(end,1));
    else % dist_f < del_tol
        penalty_func = 1;
    end
end
if comp_id==1
    Up_val(:,1) = abs(mu_GPu_g1 - bt.*std_GPu_g1).*penalty_func;
elseif comp_id==2
    Up_val(:,1) = abs(mu_GPu_g2 - bt.*std_GPu_g2).*penalty_func;
end
f = Up_val;

if nargout > 1
    g = [];
end
end % function end
