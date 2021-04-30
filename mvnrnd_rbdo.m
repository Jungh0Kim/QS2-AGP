function [MCsample_Gp] = mvnrnd_rbdo(d_sample,N_mc)
% Generate multivariate normal random variables
% Jungho Kim
% d denotes design parameters (mean value) here
global stdx nd
% size(d_sample) = (1,nd) 
% When non-normal variables are employed, generate multi-variate samples
% by Nataf transformation from standard normal space U (also required when
% the optimizer updates the mean values during optim process)
u_sample_MCS = mvnrnd(zeros(1,nd),diag(ones(1,nd)).^2,N_mc);
for i=1:nd
    MCsample_Gp(:,i) = u_sample_MCS(:,i).*stdx(1,i) + d_sample(1,i); 
end
end % function end
