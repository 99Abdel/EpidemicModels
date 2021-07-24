function [J,Jv] = fitness_funSEIR(x,t,N,y0)

    [T,Jv] = ode45(@fun, t, y0);
    
    %differential equations to ode45
    function dY = fun(t,y)
        
        dydt = zeros(4,1);
        a = x(1);
        b = x(2);
        c = x(3);
        dydt(1) = -a*y(1)*y(3)/N;          %-a*S*I/N
        dydt(2) =  a*y(1)*y(3)/N - b*y(2); % a*S*I/N - b*E
        dydt(3) =  b*y(2) - c*y(3);        % b*E - c*I
        dydt(4) =  c*y(3);                 % c*I
        dY = dydt;
       
    end

J = Jv(:,3);

end