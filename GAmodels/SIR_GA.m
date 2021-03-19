clear all

% uncomment the data which you want to use.
tmax = 25;
N = 1e4;
dt = 0.5;
%% first  data curve

t = 0:dt:tmax;
% data = xlsread('MyData.xlsx', 'Foglio1', 'B1:B100'); %infected people per day

%% second data curve 
% Tmax = 55;
% t = 1:dt:Tmax;
% data = xlsread('MyData.xlsx', 'Foglio1', 'B211:B265');

%% Loading data into arrays
% I = data(:,1)'; %copy data into array
% I0 = I(1);
% SIRT = [1000 20 0
%     950 50 30
%     800 100 20
%     600 300 60
%     400 330 100
%     100 100 40
%     80 30 90];
% SIR = SIRT';
%% setting conditions for the ga
A = []; 
b = [];
Aeq = [];
beq = [];
lb = [0 0];
ub = [5 5];
variables = 2;
%% defining the fitness function and the distance fitness function which will be inside the ga
load data.mat
fun = @(x) fitness_fun(x,t,N);
objFun = @(x) norm(fun(x) - y);

%% solution given by the ga 

%opts = optimoptions('ga', 'PopulationSize',N, 'InitialPopulationMatrix',randi(1E+4,N,variables)*1E-3, 'MaxGenerations',50, 'PlotFcn',@gaplotbestf, 'PlotInterval',1);

[x,fval] = ga(objFun,variables,A,b,Aeq,beq,lb,ub,[],[]);

%% graph of the solution compared withj the real one
axes();
plot(t,y(:,2), 'b+');
hold on
J = fitness_fun(x,t,N);
plot(t, J(:,2), 'r-');
%legend({'Data points', 'Fitted Curve'})
