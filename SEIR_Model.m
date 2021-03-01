%% Desription
% SEIR model with vaccine 
% This model studies the grouth of an infection dependent on the variables S, E, I and R with the help of a
% vaccine to stop the epidemic (in function of deelay time and rate of vaccination)
% the core parameters are the a, b and c which describes the epidemic evolution speed;
% a rappresents the exposition rate in other words the speed in which each susceptible individual becomes exposed.
% b describes the infection rate, the speed of becoming infectious for each
% exposed individual.
% c this is the removal rate, healing (dead) speed.
% This parameters are used in the following differential equations

% in my last program:
% I have decided to use directly an adimensional number of population, so
% S + E + I + R = 1.
% it is only a semplification other wise i would have had to devide each
% equation with the parameter N (total number of population)

% IN THIS CASE I do USE S + E + I + R = N, and i have devided later oe equation by N; 

%                                aIS/N        bE        cI
% this is the block diagram: S  ------->  E  ---->  I  ---->  R

%{ dS/dt = -aIS/N  (minus because people are extracted not added)
%{ dE/dt = aIS/N - bE 
%{ dI/dt = bE - cI
%{ dR/dt = cI 

% this equations are easily taken from the block diagram, everithing wich
% enters to our variable has a +, and the outgoing parameters has a -.
% we need to consider that this is an invariant time system.

%% Import Data
data = xlsread('MyData.xlsx', 'Foglio1', 'B1:B123');
D = data(:,1); %copy data into array
data1 = xlsread('MyData.xlsx', 'Foglio1', 'B123:B355');
D2 = data1(:,1); %copy data into array


%% First wave parapeters
% (with this exact parameter i was abdle to describe the italian situation 
% in the first period) first wave. only the growth because it has taken
% place a control during the epidemic.
% I0 = 2; a1 = 1;
% b = 0.805; c = 0.61;  t2 = 15;
% data taken from D;
% for a better envelope I0 = 700; and a1 = 0.8; t2 = 27;
%% Second wave data
% (with this exact parameter i describe the italian situation 
% in the second period) second wave only growth.
% I0 = 190; a1 = 0.688;
% b = 0.805; c = 0.61;  t2 = 143;
% data taken from D2;

%% Set Parameters

N = 60e6;   % total number of population (if 1 is already adimensional)

I0 = 190;     % infected people (in the beginning) (if less than one is a portion of population)

a = 0.688;    % S to E coefficient (days^-1) (exposition rate)
b = 0.805;  % E to I coefficient (days^-1) (infection rate)
c = 0.61;   % I to R coefficient (days^-1) (removal rate)


tDelay = 20; % delay time of the vaccine
dv = 0.0;    % vaccination rate in proportion of population per day dv/dt
tmax = 300;  % number of days to plot
dt = .01;   % size of time steps in days
Imax = 8e4;  % Max number of the graph (ordinata max), because i want to plot certain range

plotCase = 5; % the graph to plot 1= S, 2= E, 3= I, 4 = R , 5= All, 6 = real vs calculated I;


%% Initialize Vectors
t = 0:dt:tmax;   % time vector
Nt = length(t);  % number of time steps

S = zeros(1,Nt); % Supsectable vector (initialization)
E = zeros(1,Nt); % Exposed vector     (initialization)
I = zeros(1,Nt); % Infected vector    (initialization)
R = zeros(1,Nt); % Removed vector     (initialization)

I(1) = I0;   % first element of the Infected vector 


%% calculations

for i = 1:Nt-1
    
    S(i) = N - E(i) - I(i) - R(i);  %total susceptible people in this day
    
    dE = a*I(i)*S(i)/N - b*E(i);    %Exposition rate of change (dE/dt)
    
    E(i+1) = E(i) + dE*dt;          %Exposed people in this day
    
    dI = b*E(i) - c*I(i);       %Infection rate of change (dI/dt)
    I(i+1) = I(i) + dI*dt;      %Infected people in this day
    
    if S(i) > 0 && t(i) >= tDelay  %if not all are vaccinated
        dR = b*I(i) + dv;          %rate of change per day of recovery
        
    else %if all are vaccinated
        dR = c*I(i);
    end
    
    R(i+1) = R(i) + dR*dt;     %total removed people in the day
end

S(Nt) = N - E(Nt) - I(Nt) - R(Nt); %last S value

%% Plots

switch plotCase
    
    case 1   % S vs t
        plot(t, S, '-b', 'LineWidth', 2)
        axis([0 tmax 0 Imax])
        grid on
        grid minor
        xlabel('Time(days)')
        ylabel('proportion of Susceptible')
        title('proportion of Susceptible vs time')
    
    case 2   % E vs t
        plot(t, E, '-g', 'LineWidth', 2)
        axis([0 tmax 0 Imax])
        grid on
        grid minor
        xlabel('Time(days)')
        ylabel('proportion of Exposed')
        title('proportion of Exposed vs time')
        
        
    case 3   % I vs t
        plot(t,I, '-r', 'LineWidth', 2)
        axis([0 tmax 0 Imax])
        grid on
        grid minor
        xlabel('Time(days)')
        ylabel('proportion of Infected')
        title('proportion of Infection vs Time')
        
    case 4   % R vs t
        plot(t,R, '-m', 'LineWidth', 2)
        axis([0 tmax 0 Imax])
        grid on
        grid minor
        xlabel('Time(days)')
        ylabel('proportion of Removed')
        title('proportion of Removed vs Time')
        
        
    case 5   % SEIR vs t
        plot(t,S,'-b',...
            t,E,'-g',...
            t,I, '-r',...
            t,R,'-m','LineWidth', 2)
        axis([0 tmax 0 Imax])
        grid on
        grid minor
        xlabel('Time(days)')
        ylabel('proportion of S E I R')
        title('proportion of S E I R vs Time')
        legend('S','E','I','R')
    
end
