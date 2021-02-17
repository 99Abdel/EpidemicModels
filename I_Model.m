% Easy I model
%This model studies the grouth of an infection only 
%dependent on the variable I

%% Set Parameters
I0 = 3; % initial number of infected people
a = 0.8; % coefficient in days^-1 (speed of the infection)

tmax = 100; % number of days to plot

dt = 0.1; %size of time steps in weeks

Imax = 10.6e6;%Max number infected per graph in millions

%% Initialize Vectors

t = 0:dt:tmax; % vettore di tempo
Nt = length(t); % number of time steps
I = zeros(1,Nt); % initialize the vector of infections that i will plot in the graph,with zeros
I(1) = I0; % first element of the vector initialized
%% calculations 

for i = 1:Nt-1
    
    dI = a*I(i); % rate of change per day
    I(i+1) = I(i) + dI*dt ;
    
end 
%% Plots 

plot(t,I, '-r', 'LineWidth', 2)
axis([0 tmax 0 Imax])
grid on
grid minor
xlabel('Time(days)')
ylabel('Number of infected')
title('Number of infection vs Time')
