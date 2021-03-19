clear all

% uncomment the data which you want to use.

%% first  data curve
t = 1:1:27;
data = xlsread('MyData.xlsx', 'Foglio1', 'B1:B28'); %infected people per day

%% second data curve 
% t = 1:1:55;
% data = xlsread('MyData.xlsx', 'Foglio1', 'B211:B265');

%% Loading data into arrays
y = data(:,1); %copy data into array
I = y';  % vector of infected people per day
I0 = I(1);

%% setting conditions for the ga
A = []; 
b = [];
Aeq = [];
beq = [];
lb = [0];
ub = [1];

%% defining the fitness function and the distance fitness function which will be inside the ga
fun = @(a) myfun(a,I0,t);
objFun = @(a) norm(fun(a)-I);

%% solution given by the ga 
sol = ga(objFun, 1,A,b,Aeq,beq,lb,ub);

%% graph of the solution compared withj the real one
axes();
plot(t, I, 'b+');
hold on
plot(t, fun(sol), 'r-');
legend({'Data points', 'Fitted Curve'})

function y = myfun(a,I0,t)

y = (I0*exp(a(1)*t));

end

