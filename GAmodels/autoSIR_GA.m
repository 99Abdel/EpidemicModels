clear all

%prima salita [2:28]
%prima discesa [28:170]
%seconda salita [170:265]
%seconda discesa [265:359]
%terza salita [359:383]

t3 = [2 28 170 265 294 311 343 383]; % vettore tempi degli interventi del governo
n = length(t3);
dt = 1;
tmax = t3(n)-t3(1);
t = 0:dt:tmax;

% settagio iniziale con i dati

N = 60e6;
FIRST = 1;
y = letturaExcel(t3(1),t3(n));
D = y;
y0 = [y(FIRST,1),y(FIRST,2),y(FIRST,3)];

%% setting conditions for the ga
A = [];
b = [];
Aeq = [];
beq = [];
lb = [0 0];
ub = [2 1];
variables = 2;

pop = 50;
maxGen = 30;

%preparo i vettori dove salvare i risultati del GA.
c = n-1;
parametri = zeros(c,variables);
FVAL_TOT = zeros(c,1);


%% defining the fitness function and the distance fitness function which will be inside the ga

fun = @(x) autofitness_fun(x,t,N,y0);
objFun = @(x) norm(fun(x) - y(:,2));

opts = optimoptions('ga', 'PopulationSize',pop, 'TolFun',1e-5,'MaxGenerations',maxGen, 'PlotFcn',@gaplotbestf, 'PlotInterval',1);

%% solution given by the ga

[x,fval] = ga(objFun,variables,A,b,Aeq,beq,lb,ub,[],[],opts);

[J,Jv] = fitness_fun(x,t,N,y0);


%% graph of the solution compared with the real one
figure
tiledlayout(2, 2);
nexttile
plot(t,D(:,1), 'b+');
hold on
plot(t, Jv(:,1), 'm-');

nexttile([1 2])
plot(t,D(:,2), 'b+');
hold on
plot(t, Jv(:,2), 'r-');
legend({'Data points', 'Fitted Curve'})

nexttile
plot(t,D(:,3), 'b+');
hold on
plot(t, Jv(:,3), 'g-');


%% funzione lettura foglio excel

function y = letturaExcel(t1,t2)

range1 = "B"+t1+":B"+t2;
range2 = "C"+t1+":C"+t2;
range3 = "D"+t1+":D"+t2;

data = xlsread('MyData.xlsx', 'Foglio1', range1);
S = data(:,1); %copy data into array
data = xlsread('MyData.xlsx', 'Foglio1', range2);
I = data(:,1); %copy data into array
data = xlsread('MyData.xlsx', 'Foglio1', range3);
R = data(:,1); %copy data into array

y = [S I R];

end




% for i = 1:c
%     
%     parametri(i,1) = x(1);
%     parametri(i,2) = x(2);
%     FVAL_TOT(i,1) = fval;
%     
%     %% cloning the pieces of the vector J into one only
%     
%     if i ==1
%         for j = 1:length(J)
%             pos = j+FIRST;
%             F(pos) = J(j);
%         end
%         pos = pos+1;
%     else
%         for j = 1:length(J)-1
%             pos = j+FIRST-i+1;
%             F(pos) = J(j+1);
%         end
%     end
%     
%     FIRST = FIRST + length(J);
%     y0 = [Jv(length(J),1),Jv(length(J),2),Jv(length(J),3)];
%     
% end