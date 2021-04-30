function df = calcdf_mu(model,xs,l_calc)
% df = calcdf_mu(model,x, xs,l_calc)
% Function to calculate the derivative of the exact posterior mean of a GP about
% set of test points (GP with uncertain inputs)
% Inputs
%  model          model struct 
%      .x       training inputs, N-by-D
%      .y       training targets, N-by-E
%      .seard   squared exp kernel hyperparameters
%      .lsipn   log input noise standard deviation, D-by-1
%      .R        cholesky of training cov matrix, N-by-N-by E
%      .alpha   Kn^-1 * y, N-by-E
%   x            training inputs matrix, N-by-D
%   xs           a set of test points at which to evaluate the derivative
% l_calc         l matrix, N-by-Ns
%
% Output
%   df            matrix of derivatives, N-by-D
%
% Andrew McHutchon, July 2012, some modifications are made by Jungho Kim

x = model.x;
[N D] = size(x); E = size(model.seard,2);

% iell2 = exp(-2*lhyp.seard(1:D,:));                    % squared lengthscales D-by-E
sx_vec = exp(2*model.lsipn(1:D,:));
ell_vec = exp(2*model.seard(1:D,:));
iellsx2 = (sx_vec + ell_vec).^-1;   % squared lengthscales + input noise variance D-by-E

% Test set
Ns = size(xs,1); df = zeros(Ns,E,D);
% Find the derivative of the covariance function
% model.x = x; model.y = model.y;
% [d,Ks] = calcKtest(model,xs);

for i=1:E
    Xmx = bsxfun(@minus,permute(x,[1,3,2]),permute(xs,[3,1,2])); % N-by-Ns-by-D
    %         XmxiLam = bsxfun(@times,Xmx,permute(iell2(:,i),[3,2,1])); % N-by-Ns-by-D
    XmxiLam = bsxfun(@times,Xmx,permute(iellsx2(:,i),[3,2,1])); % N-by-Ns-by-D
    %         dKtsdx = bsxfun(@times,XmxiLam,Ks(:,:,i));                  % N-by-Ns-by-D
    dKtsdx = bsxfun(@times,XmxiLam,l_calc(:,:,i));                  % N-by-Ns-by-D
    % Compute derivative
    df(:,i,:) = etprod('12',dKtsdx,'412',model.alpha(:,i),'4');  % Ns-by-E-by-D
end
df = reshape(df,[Ns,D]);

end % function end
