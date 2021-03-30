clear all

%prima salita [2:28]
%prima discesa [28:170]
%seconda salita [170:265]
%seconda discesa [265:359]
%terza salita [359:383]


%% setting conditions for the ga
A = [];
b = [];
Aeq = [];
beq = [];
lb = [0 0];
ub = [2 1];
variables = 2;

parametri = zeros(5,2);
F = zeros;

%% Work station

for i = 1:5
    
    t3 = [2 28 170 265 343 383];
    range1 = "B"+t3(i)+":B"+t3(i+1);
    range2 = "C"+t3(i)+":C"+t3(i+1);
    range3 = "D"+t3(i)+":D"+t3(i+1);
    
    data = xlsread('MyData.xlsx', 'Foglio1', range1);
    B = data(:,1); %copy data into array
    data = xlsread('MyData.xlsx', 'Foglio1', range2);
    C = data(:,1); %copy data into array
    data = xlsread('MyData.xlsx', 'Foglio1', range3);
    D = data(:,1); %copy data into array
    
    FIRST = 1;
    y = [B C D];
    y0 = [B(FIRST),C(FIRST),D(FIRST)];
    
    tmax = t3(i+1)-t3(i);
    N = 60e6;
    dt = 1;
    t = 0:dt:tmax;
    
    %% defining the fitness function and the distance fitness function which will be inside the ga
    
    fun = @(x) fitness_fun(x,t,N,y0);
    objFun = @(x) norm(fun(x) - y(:,2));
    
    %% solution given by the ga
    pop = 280;
    maxGen = 50;
    opts = optimoptions('ga', 'PopulationSize',pop, 'TolFun',1e-5,'MaxGenerations',maxGen, 'PlotFcn',@gaplotbestf, 'PlotInterval',1);
    
    [x,fval] = ga(objFun,variables,A,b,Aeq,beq,lb,ub,[],[],opts);
    J = fitness_fun(x,t,N,y0);
    
    parametri(i,1) = x(1);
    parametri(i,2) = x(2);
    
    for j = 1:size(J)-1
        F(size(F) + 1) = J(j);
    end
end

%% data to plot later 
t1 = 2; t2 = 383;
tmax = t2-t1;
tg = 0:dt:tmax;

range1 = "B"+t1+":B"+t2;
range2 = "C"+t1+":C"+t2;
range3 = "D"+t1+":D"+t2;

data = xlsread('MyData.xlsx', 'Foglio1', range1);
S = data(:,1); %copy data into array
data = xlsread('MyData.xlsx', 'Foglio1', range2);
I = data(:,1); %copy data into array
data = xlsread('MyData.xlsx', 'Foglio1', range3);
R = data(:,1); %copy data into array

M = [S I R];

%% graph of the solution compared withj the real one
tiledlayout(2, 2);
nexttile
plot(tg,M(:,1), 'b+');
%hold on
%plot(tg, J(:,1), 'm-');

nexttile([1 2])
plot(tg,M(:,2), 'b+');
hold on
plot(tg,F, 'r-');
legend({'Data points', 'Fitted Curve'})

nexttile
plot(tg,M(:,3), 'b+');
%plot(t, J(:,3), 'g-');




