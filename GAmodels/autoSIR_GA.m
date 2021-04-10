clear all

%prima salita [2:28]
%prima discesa [28:170]
%seconda salita [170:265]
%seconda discesa [265:359]
%terza salita [359:383]

% t3 = [2 28 170 265 294 311 343 383]; % vettore tempi degli interventi del governo

% settagio iniziale con i dati
t1 = 2; t2 = 383;
N = 59999771;
FIRST = 1;
y = letturaExcel(t1,t2);
y0 = [y(FIRST,1),y(FIRST,2),y(FIRST,3)];

t3 = [t1 timeChanges(t1,t2,y(:,2))];

n = length(t3);
tmax = t3(n)-t3(1);
%% setting conditions for the ga
A = [];
b = [];
Aeq = [];
beq = [];
lb = [0 0];
ub = [5 2];
variables = 2;

%preparo i vettori dove salvare i risultati del GA.
c = n-1;
parametri = zeros(c,variables);
FVAL_TOT = zeros(c,1);
R0 = zeros(c,1);
%condizioni per l'ingresso al ciclo for
F = zeros(tmax,3);
FIRST = 0;

%% Work station

for i = 1:c
    
    tmax = t3(i+1)-t3(i);
    dt = 1;
    t = 0:dt:tmax;
    y = letturaExcel(t3(i),t3(i+1));
    %% defining the fitness function and the distance fitness function which will be inside the ga
    
    fun = @(x) fitness_fun(x,t,N,y0);
    objFun = @(x) norm(fun(x) - y(:,2));
    
    %% solution given by the ga
    pop = 380;
    maxGen = 100;
    opts = optimoptions('ga', 'PopulationSize',pop, 'TolFun',1e-5,'MaxGenerations',maxGen, 'PlotFcn',@gaplotbestf, 'PlotInterval',1);
    
    [x,fval] = ga(objFun,variables,A,b,Aeq,beq,lb,ub,[],[],opts);
    [J,Jv] = fitness_fun(x,t,N,y0);
    
    parametri(i,1) = x(1);
    parametri(i,2) = x(2);
    R0(i,1) = x(1)/x(2);
    FVAL_TOT(i,1) = fval;
    
%% cloning the pieces of the vector J into one only
    
    if i ==1
        for j = 1:length(J)
            pos = j+FIRST;
            F(pos,1) = Jv(j,1);
            F(pos,2) = Jv(j,2);
            F(pos,3) = Jv(j,3);
        end
        pos = pos+1;
    else
        for j = 1:length(J)-1
            pos = j+FIRST-i+1;
            F(pos,1) = Jv(j+1,1);
            F(pos,2) = Jv(j+1,2);
            F(pos,3) = Jv(j+1,3);
        end
   end
    
    FIRST = FIRST + length(J);
    y0 = [Jv(length(J),1),Jv(length(J),2),Jv(length(J),3)];
    
end

%% data to plot later 
t1 = t3(1); t2 = t3(n);
tmax = t2-t1;
tg = 0:dt:tmax;

D = letturaExcel(t1,t2);


%% graph of the solution compared with the real one
figure
tiledlayout(2, 2);
nexttile
plot(tg,D(:,1), 'b+');
hold on
plot(F(:,1), 'm-');

nexttile([1 2])
plot(tg,D(:,2), 'b+');
hold on
plot(F(:,2), 'r-');
legend({'Data points', 'Fitted Curve'})

nexttile
plot(tg,D(:,3), 'b+');
hold on
plot(F(:,3), 'g-');


