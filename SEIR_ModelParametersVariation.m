%% Desription
% SEIR model with vaccine 
% This model studies the grouth of an infection dependent on the variables S, E, I and R with the help of a
% vaccine to stop the epidemic (in function of deelay time and rate of vaccination)
% the core parameters are the a, b and c which describes the epidemic evolution speed;
% a rappresents the exposition rate in other words the speed in which each susceptible individual becomes exposed.
% b describes the infection rate, the speed of becoming infectious for each
% exposed individual.
% c this is the removal rate, healing (dead) speed.
% This parameters are used in the following differential equations:

% S + E + I + R = N total number of population
% equation with the parameter N (total number of population)

%                                aIS/N        bE        cI
% this is the block diagram: S  ------->  E  ---->  I  ---->  R

%{ dS/dt = -aIS  (minus because people are extracted not added)
%{ dE/dt = aIS - bE
%{ dI/dt = bE - cI
%{ dR/dt = cI

% this equations are easily taken from the block diagram, everithing wich
% enters to our variable has a +, and the outgoing parameters has a -.
% we need to consider that this is an variant time system.

%% First wave parapeters
% (with this exact parameter i describe the italian situation 
% in the first period) first wave.
% I0 = 2; 
% a1 = 2; 
% a2 = 0.565;
% b = 0.805; 
% c = 0.61; 
% t2 = 15;
% data taken from D;
% for a better envelope I0 = 700; a1 = 0.8; t2 = 27;

%% Second wave data
% (with this exact parameter i describe the italian situation 
% in the second period) second wave.
% I0 = 190; 
% a1 = 0,688;
% a2 = 0.565;
% b = 0.805; 
% c = 0.61; 
% t2 = 143;
% data taken from D2;

%% Total Wave Data
% a3 = 0,63
% t3 = 174;
% 

%% Import Data
data = xlsread('MyData.xlsx', 'Foglio1', 'B1:B123'); %first wave
D = data(:,1); %copy data into array
data1 = xlsread('MyData.xlsx', 'Foglio1', 'B123:B355'); %se3cond wave
D2 = data1(:,1); %copy data into array

totalData = xlsread('MyData.xlsx', 'Foglio1', 'B1:B355'); %both waves total graph
D3 = totalData(:,1);


%% Set Parameters

I0 = 700; % infected people (in the beginning)

a1 = 0.8;   % S to E coefficient (days^-1) (exposition rate) before intevention
a2 = 0.565;   % S to E coefficient (days^-1) (exposition rate) after intervention
a3 = 0.75;    % another change

a = a1;    % S to E coefficient (days^-1) (exposition rate) in the equation

N = 60e6; % total number of population (italy)

b = 0.805; % E to I coefficient (days^-1) (infection rate)
c = 0.61; % I to R coefficient (days^-1) (removal rate)

t2 = 27;     % Time passed until the intervention phase (days)
t3 = 174;     % third time...

tVac = 20;   % delay time of the vaccine (days)
dv = 0;    % vaccination rate in proportion of population per day dv/dt

tmax = 400;   % number of days to plot
dt = 0.01;   % size of time steps in days
Imax = 50000;  % Max number of the graph (ordinata max), because i want to plot certain range

plotCase = 3; % the graph to plot 1= S, 2= E, 3= I, 4 = R , 5= All;


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
    
         if t(i) > t2 && t(i) < t3
            a = a2;
        elseif t(i) > t3
                a = a3;
         
        end
    
    dE = a*I(i)*S(i)/N - b*E(i);      %Exposition rate of change (dE/dt)
    E(i+1) = E(i) + dE*dt;          %Exposed people in this day
    
    dI = b*E(i) - c*I(i);       %Infection rate of change (dI/dt)
    I(i+1) = I(i) + dI*dt;      %Infected people in this day
    
    if S(i) > 0 && t(i) >= tVac  %if not all are vaccinated
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
        ylabel('number of Susceptible')
        title('number of Susceptible vs time')
    
    case 2   % E vs t
        plot(t, E, '-g', 'LineWidth', 2)
        axis([0 tmax 0 Imax])
        grid on
        grid minor
        xlabel('Time(days)')
        ylabel('number of Exposed')
        title('number of Exposed vs time')
        
        
    case 3   % I vs t
        plot(t,I, '-r', 'LineWidth', 2)
        axis([0 tmax 0 Imax])
        grid on
        grid minor
        xlabel('Time(days)')
        ylabel('number of Infected')
        title('number of Infection vs Time')
        
    case 4   % R vs t
        plot(t,R, '-m', 'LineWidth', 2)
        axis([0 tmax 0 Imax])
        grid on
        grid minor
        xlabel('Time(days)')
        ylabel('number of Removed')
        title('number of Removed vs Time')
        
        
    case 5   % SEIR vs t
        plot(t,S,'-b',...
            t,E,'-g',...
            t,I, '-r',...
            t,R,'-m','LineWidth', 2)
        axis([0 tmax 0 Imax])
        grid on
        grid minor
        xlabel('Time(days)')
        ylabel('number of S E I R')
        title('number of S E I R vs Time')
        legend('S','E','I','R')
        
end

function b = myfunction(c,d)

    b = c+d;

end