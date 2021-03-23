%SIR model
%This model studies the grouth of an infection
%dependent on the variable I, the parameter S and R.


%% Set Parameters
I0 = 1; % portion of infected people

a = 0.2; % coefficient in days^-1 (speed of the infection)

b = 0.05; % coefficient in days^-1 (speed of recovery)

tmax = 122; % number of days to plot

N = 4e6;

dt = 0.01; %size of time steps in weeks

Imax = N + 1e5;%Max number infected per graph in millions

plotCase = 4; %the graph to plot 1=S, 2=I, 3=R, 4=All;


%% Initialize Vectors
t = 0:dt:tmax; % vettore di tempo
Nt = length(t); % number of time steps
S = zeros(1,Nt); %initialize the vector of supsectable people
I = zeros(1,Nt); % initialize the vector of infections that i will plot in the graph,with zeros
R = zeros(1,Nt); % initialize the removed vector
I(1) = I0; % first element of the vector infected initialized
S(1) = N-I0; %first element of the vector supsectable initialized


%% calculations

for i = 1:Nt-1
    
    S(i) = N-I(i)-R(i); %total susceptible people in this day
    
    dI = a*I(i)*S(i)/N -b*I(i);  % rate of change per day of infection(dI = dI/dt)
    
    I(i+1) = I(i) + dI*dt;  %total infected people in the day
    
    dR = b*I(i);  %rate of change per day of recovery ()
    
    R(i+1) = R(i) + dR*dt; %total removed people in the day
      
end
S(Nt) = N - I(Nt) - R(Nt);

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
        plot(t,R, '-g', 'LineWidth', 2)
        axis([0 tmax 0 Imax])
        grid on
        grid minor
        xlabel('Time(days)')
        ylabel('proportion of Removed')
        title('proportion of Removed vs Time')
        
    
    case 4
        
        plot(t,S,'-b',...
            t,I, '-r',...
            t,R,'-g','LineWidth', 2)
        axis([0 tmax 0 Imax])
        grid on
        grid minor
        xlabel('Time(days)')
        ylabel('proportion of I S R')
        title('proportion of I S R vs Time')
        legend('S','I','R')
end
