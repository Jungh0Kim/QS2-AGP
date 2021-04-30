function [c,ceq,GC,GCeq] = prob_func_fmc(x)
% G quantile (Gp) and its gradient function (dGp/dtheta)
% Jungho Kim
global bt model1 model2

% There are no equality constrints
ceq = []; GCeq = [];

[mu_GPu_g1, s2_GPu_g1, l_calc1, Ldf_calc1] = gpm_m(model1,x);
std_GPu_g1 = sqrt(s2_GPu_g1);
[mu_GPu_g2, s2_GPu_g2, l_calc2, Ldf_calc2] = gpm_m(model2,x);
std_GPu_g2 = sqrt(s2_GPu_g2);
   
g_p1 = mu_GPu_g1 - bt.*std_GPu_g1;
g_p2 = mu_GPu_g2 - bt.*std_GPu_g2;

grad_calc_mu = calcdf_mu(model1, x, l_calc1);
grad_calc_s2 = Ldf_calc1 - 2.*mu_GPu_g1.*grad_calc_mu;
grad_gp1 = grad_calc_mu - bt./2./sqrt(s2_GPu_g1).*grad_calc_s2;
grad_calc_mu = calcdf_mu(model2, x, l_calc2);
grad_calc_s2 = Ldf_calc2 - 2.*mu_GPu_g2.*grad_calc_mu;
grad_gp2 = grad_calc_mu - bt./2./sqrt(s2_GPu_g2).*grad_calc_s2;

c = -[g_p1; g_p2]; GC = -[grad_gp1; grad_gp2];
c = c'; GC = GC';
end % function end
