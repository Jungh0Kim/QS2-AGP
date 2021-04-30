function [M, S, l_calc, Ldf_calc] = gpm_m(model, m)
% Function to make NIGP predictions for Gaussian distributed test points.
% gpm.m in NIGP by Andrew McHutchon, some modifications are made by Jungho Kim

input = model.x; target = model.y; alpha = model.alpha;
s = diag(exp(2*model.lsipn));     % input cov matrix is diag of input noise, Sigma.x matrix
[N D] = size(input);                                  % Training dimensions
[Ns Ds] = size(m);                                        % Test dimensions
if Ds ~= D; error('Test point(s) supplied wrongly'); end
E = size(target,2);

M = zeros(Ns,E); S = zeros(Ns,E);
l_calc = zeros(N,Ns); del_k = zeros(N,D);
del_Qexp1 = zeros(N,N,D); del_Qexp2 = zeros(N,N,D);
Ldf_calc = zeros(Ns,E,D);
iella = exp(-model.seard(1:D,:))'; % E-by-D, 1/lambda
sf2a = exp(2*model.seard(D+1,:));

for n = 1:Ns
    inp = bsxfun(@minus,input,m(n,:)); % N-by-D, (x-theta)
    for i=1:E
        L = diag(iella(i,:));                                      % D-by-D, L = 1/lmabda
        iLam = L.^2;                                               % D-by-D, Lambda matrix (L^2)
        in = inp*L;         % N-by-D, training data minus m divided by ells, (x-theta)*L^2
        B = L*s*L+eye(D);                                          % D-by-D, Sigma.x/lmabda^2 diag matrix
        iB = B\eye(D);      % D-by-D
        t = in*iB;          % in*inv(B) - N-by-D, O(D^3 x E)
        l = exp(-sum(in.*t,2)/2);  % N-by-1
        lb = l.*alpha(:,i);    % N-by-1
        c = exp(2*model.seard(D+1,i))/sqrt(det(B));
  
        M(n,i) = c*sum(lb); % Ns-by-E
        l_calc(:,n) = c.*l; % N-by-Ns
        
        k = 2*model.seard(D+1,i)-sum((inp*L).*(inp*L),2)/2;       % N-by-1, O(ND x E), log(sig.f^2) - (x-theta)^2*L^4/2
        ii = bsxfun(@times,inp,iella(i,:).^2);               % N-by-D, O(ND x E), (x-theta)/lambda^2
        
        C = 2*s*iLam+eye(D);                     % D-by-D, O(D^2 x E/2 x E), 2*Sigma.x*L^2 + I
        t = 1/sqrt(det(C));                        % scalar, O(D^3 x E^2/2)
        Qexp = exp(bsxfun(@plus,k,k') + maha(ii,-ii,C\s/2)); % N-by-N, O(N^2 xE^2/2) or O(D^3 xE^2/2),
        
%         del_k_1D = -sum(2.*inp*L.^2,2)/2;                       % must be N-by-D                         
%         del_Qexp1_1D = bsxfun(@plus,del_k_1D,del_k_1D');  % N-by-N-by-D 
%         del_Qexp2_1D = 4.*(C\s/2).*iella(i,:).^2.*bsxfun(@times,ones(N),2.*ii); % N-by-N-by-D
%         del_Qexp = -del_Qexp1_1D - del_Qexp2_1D;
               
        A = alpha(:,i)*alpha(:,i)';                 % N-by-N, O(N^2 xE^2/2)
        A = A - solve_chol(model.R(:,:,i),eye(N)); % incorporate model uncertainty
        Aq = A.*Qexp;                                % N-by-N, O(N^2 xE^2/2)
        
        S(n,i) = sf2a(i) + t*sum(sum(Aq)) - M(n,i)^2;        % O(N^2 xE^2/2)
        if nargout > 2   % gradient calculation
            Csc = (C\s/2).*iella(i,:).^2;    % D-by-D, E-by-D
            for Dn=1:D
                del_k(:,Dn) = -sum(2.*inp(:,Dn)*L(Dn,Dn).^2,2)/2;     % must be N-by-D
                del_Qexp1(:,:,Dn) = bsxfun(@plus,del_k(:,Dn),del_k(:,Dn)');  % N-by-N-by-D
                del_Qexp2(:,:,Dn) = 4.*Csc(Dn,Dn).*bsxfun(@times,ones(N),2.*ii(:,Dn));  % N-by-N-by-D 
                del_Qexp = - del_Qexp1 - del_Qexp2; % should be N-by-N-by-D
                Ldf_calc(n,:,Dn) = t*sum(sum(Aq.*del_Qexp(:,:,Dn)));
            end
        end
    end
end

if nargout > 2 % gradient calculation
    Ldf_calc = reshape(Ldf_calc,[Ns,D]);
end

S = bsxfun(@plus,S,exp(2*model.seard(end,:)));  % Add output noise

end % function end
