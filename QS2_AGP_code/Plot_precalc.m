% Jungho Kim
% pre-calculation for surface plot (only required for visualizing(plotting))

global meanfunc covfunc likfunc model1 model2 bt

[mu_GPu_g1, s2_GPu_g1] = gpm_m(model1,d_plot);
std_GPu_g1 = sqrt(s2_GPu_g1);
[mu_GPu_g2, s2_GPu_g2] = gpm_m(model2,d_plot);
std_GPu_g2 = sqrt(s2_GPu_g2);
g_p1 = mu_GPu_g1 - bt.*std_GPu_g1; % estimated quantile surface
g_p2 = mu_GPu_g2 - bt.*std_GPu_g2;
[mu_GP_g1,s2_GP_g1] = gp(hyp_optim_g1, @infGaussLik, meanfunc, covfunc, likfunc, x_doe1, G_val1, d_plot);
[mu_GP_g2,s2_GP_g2] = gp(hyp_optim_g2, @infGaussLik, meanfunc, covfunc, likfunc, x_doe2, G_val2, d_plot);
mu_GP_g1 = reshape(mu_GP_g1 ,size(xs1)); % estimated limit-state surface
mu_GP_g2 = reshape(mu_GP_g2 ,size(xs1));
