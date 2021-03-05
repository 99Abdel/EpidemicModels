clear all

% uncomment the data which you want to use.

%% first  data curve
xx = 1:1:27;
data = xlsread('MyData.xlsx', 'Foglio1', 'B1:B28');

%% second data curve 
% xx = 1:1:55;
% data = xlsread('MyData.xlsx', 'Foglio1', 'B211:B265');

%% Loading data into arrays
y = data(:,1); %copy data into array
yy = y';
I0 = yy(1);

%% setting conditions for the ga
A = []; 
b = [];
Aeq = [];
beq = [];
lb = [0];
ub = [1];

%% defining the fitness function and the distance fitness function which will be inside the ga
fun = @(a) (I0*exp(a(1)*xx));
objFun = @(a) norm(fun(a)-yy);

%% solution given by the ga 
sol = ga(objFun, 1,A,b,Aeq,beq,lb,ub);

%% graph of the solution compared withj the real one
axes();
plot(xx, yy, 'b+');
hold on
plot(xx, fun(sol), 'r-');
legend({'Data points', 'Fitted Curve'})


