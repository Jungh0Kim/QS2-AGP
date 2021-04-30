function [f,fg] = obj_func(x)
% Cost (objective) function we want to minimize
% Jungho Kim
f = (x(:,1)-3.7).^2 + (x(:,2)-4).^2;
if nargout > 1 % gradient required
    fg = [2.*(x(:,1)-3.7), 2.*(x(:,2)-4)];
end
end % function end
