%SI model
%This model studies the grouth of an infection
%dependent on the variable I and also the parameter S

%% Set Parameters
I0 = 0.001; % portion of infected people
a = 1.0; % coefficient in days^-1 (speed of the infection)

tmax = 10; % number of days to plot

dt = 0.01; %size of time steps in weeks

Imax = 1.1;%Max number infected per graph in millions

plotCase = 3; %the graph to plot 1=S, 2=I, 3=All;

%% Initialize Vectors
t = 0:dt:tmax; % vettore di tempo
Nt = length(t); % number of time steps
S = zeros(1,Nt); %initialize the vector of supsectable people
I = zeros(1,Nt); % initialize the vector of infections that i will plot in the graph,with zeros
I(1) = I0; % first element of the vector infected initialized
S(1) = 1-I0; %first element of the vector supsectable initialized
%% calculations

for i = 1:Nt-1
    
    S(i) = 1-I(i);
    dI = a*I(i)*S(i);  % rate of change per day
    I(i+1) = I(i) + dI*dt;
      
end
S(Nt) = 1 - I(Nt);

%% Plots

switch plotCase
    case 1
        plot(t, S, '-b', 'LineWidth', 2)
        axis([0 tmax 0 Imax])
        grid on
        grid minor
        xlabel('Time(days)')
        ylabel('proportion of susceptible')
        title('proportion of susceptible vs time')
        
        
    case 2
        
        plot(t,I, '-r', 'LineWidth', 2)
        axis([0 tmax 0 Imax])
        grid on
        grid minor
        xlabel('Time(days)')
        ylabel('proportion of infected')
        title('proportion of infection vs Time')
    case 3
        
        plot(t,I, '-r',...
            t,S,'-b','LineWidth', 2)
        axis([0 tmax 0 Imax])
        grid on
        grid minor
        xlabel('Time(days)')
        ylabel('proportion of I and S')
        title('proportion of I and S vs Time')
        
end
