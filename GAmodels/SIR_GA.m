clear all

%prima salita [1:28]
%prima discesa [28:150]
%seconda salita [200:265]
%seconda discesa [265:359]
%terza salita [359:383]

%% first  data curve
t1 = 359; t2 = 383;
range1 = "B"+t1+":B"+t2;
range2 = "C"+t1+":C"+t2;
range3 = "D"+t1+":D"+t2;

data = xlsread('MyData.xlsx', 'Foglio1', range1);
B = data(:,1); %copy data into array
data = xlsread('MyData.xlsx', 'Foglio1', range2);
C = data(:,1); %copy data into array
data = xlsread('MyData.xlsx', 'Foglio1', range3);
D = data(:,1); %copy data into array

FIRST = 1;
y = [B C D];
y0 = [B(FIRST),C(FIRST),D(FIRST)];
%load data.mat

% uncomment the data which you want to use.
tmax = t2-t1;
N = 60e6;
dt = 1;

%N = 1e4;
%tmax = 51;
%dt = 0.5;

t = 0:dt:tmax;

%% second data curve 
% Tmax = 55;
% t = 1:dt:Tmax;
% data = xlsread('MyData.xlsx', 'Foglio1', 'B211:B265');

%% Loading data into arrays
% I = data(:,1)'; %copy data into array
% I0 = I(1);

%% setting conditions for the ga
A = []; 
b = [];
Aeq = [];
beq = [];
lb = [0 0];
ub = [1 1];
variables = 2;
%% defining the fitness function and the distance fitness function which will be inside the ga

fun = @(x) fitness_fun(x,t,N,y0);
objFun = @(x) norm(fun(x) - y(:,2));

%% solution given by the ga 
pop = 280;
maxGen = 30;
opts = optimoptions('ga', 'PopulationSize',pop, 'MaxGenerations',maxGen, 'PlotFcn',@gaplotbestf, 'PlotInterval',1);

[x,fval] = ga(objFun,variables,A,b,Aeq,beq,lb,ub,[],[],opts);

J = fitness_fun(x,t,N,y0);

%% graph of the solution compared withj the real one

tiledlayout(2, 2);
nexttile
plot(t,y(:,1), 'b+');
hold on
%plot(t, J(:,1), 'm-');

nexttile([2 1])
plot(t,y(:,2), 'b+');
hold on
plot(t, J, 'r-');
legend({'Data points', 'Fitted Curve'})

nexttile
plot(t,y(:,3), 'b+');
%plot(t, J(:,3), 'g-');




