clear all
tmax = 25;
N = 1e4;
dt = 0.01;
tspan = 0:dt:tmax;
y0 = [N 1 0];
[x,y] = ode45(@fun, tspan, y0);
plot(x,y(:,2),'r-','Linewidth', 2)
hold on

%% extra to comapre with ODE

t = 0:dt:tmax;
Nt = length(t);
I0 = y0(2);
I = zeros(1,Nt);
R = zeros(1,Nt);

a = 1.8;
b = 0.4;

S = zeros(1,Nt);
I(1) = I0;

for i = 1:Nt-1
    
%     I(i) = N*((I0*exp(b*t(i)))/(I0*exp(b*t(i)) + N - I0));
%     S(i) = N-I(i);
    S(i) = N-I(i)-R(i); %total susceptible people in this day
    
    dI = a*I(i)*S(i)/N -b*I(i);  % rate of change per day of infection(dI = dI/dt)
    
    I(i+1) = I(i) + dI*dt;  %total infected people in the day
    
    dR = b*I(i);  %rate of change per day of recovery ()
    
    R(i+1) = R(i) + dR*dt; %total removed people in the day
    
end 

plot(t,I, 'b-')




%% series of function for the ODE solver

function dydt = fun(t,y)
    
    N = 1e4;
    dydt = zeros(3,1);
    a = 1.8;
    b = 0.4;
    dydt(1) = -a*y(1)*y(2)/N; %-a*S*I
    dydt(2) =  a*y(1)*y(2)/N - b*y(2); % a*S*I
    dydt(3) =  b*y(2);        % b*I
    
end
