function [g] = prob_func(x,k)
% Probabilistic constraints (Gx) in RBDO
% The high-fidelity computational simulation results, e.g. finite element
% analysis, can be applied. Interface must be implemented here.
% Jungho Kim
    if k == 1
        g = -x(:,1).*sin(4.*x(:,1)) - 1.1.*x(:,2).*sin(2.*x(:,2));
    elseif k == 2
        g = x(:,1) + x(:,2) - 3;
    end
end % function end
