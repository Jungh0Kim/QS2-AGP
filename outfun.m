function stop = outfun(x,optimValues,state)
% save(or plot) the progress of optimization
% Jungho Kim
global history
stop = false;
    switch state
       case 'init'
%            hold on
       case 'iter'
           % Concatenate current point and objective function
           % value with history. x must be a row vector.
           history.fval = [history.fval; optimValues.fval];
           history.x = [history.x; x];
%            plot(x(1),x(2),'o');
%            % Label points with iteration number.
%            % Add .15 to x(1) to separate label from plotted 'o'
%            text(x(1)+.15,x(2),num2str(optimValues.iteration));
       case 'done'
%            hold off
       otherwise
    end
end % function end
