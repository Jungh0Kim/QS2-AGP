function [model] = Modelstruct(hyp_gpml,x_DoE,g_DoE)
global stdx
% [model] = Modelstruct(hyp_gpml,x_DoE,g_DoE)
% Function to make model struct of GP with uncertain input
% Inputs
%   hyp_gpml    hyperparameters of GP
%   x_DoE       DoE input
%   g_DoE       DoE ouput
% Outputs
%   model       model with R and alpha fields
%      .x       training inputs, N-by-D
%      .y       training targets, N-by-E
%      .seard   squared exp kernel hyperparameters
%      .lsipn   log input noise standard deviation, D-by-1
%      .R        cholesky of training cov matrix, N-by-N-by E
%      .alpha   [K + sn2I]^-1 * y
% Andrew McHutchon, July 2012, some modifications are made by Jungho Kim

model.seard = [hyp_gpml.cov'; hyp_gpml.lik];
model.lsipn = log(stdx)';
model.x = x_DoE;
model.y = g_DoE;

[N E] = size(model.y);
if ~isfield(model,'R')
    XmX2 = bsxfun(@minus,permute(model.x,[1,3,2]),permute(model.x,[3,1,2])).^2;
    K = calcK(model.seard, XmX2);                        % N-by-N-by-E
    sn2I = bsxfun(@times,permute(exp(2*model.seard(end,:)),[1,3,2]),eye(N));% N-by-N-by-E
    Kn = K + sn2I;
    model.R = chol3(Kn);
end
if ~isfield(model,'alpha')
    model.alpha = findAlpha(model.R,model.y);
end

end % function end
