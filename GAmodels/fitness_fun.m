function [J,Jv] = fitness_fun(x,t,N,y0)

    [T,Jv] = ode45(@fun, t, y0);
    
    %differential equations to ode45
    function dY = fun(t,y)
        
        dydt = zeros(3,1);
        a = x(1);
        b = x(2);
        dydt(1) = -a*y(1)*y(2)/N; %-a*S*I
        dydt(2) =  a*y(1)*y(2)/N - b*y(2); % a*S*I
        dydt(3) =  b*y(2);        % b*I
       dY = dydt;
       
    end

J = Jv(:,2);

end