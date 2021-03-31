t1 = 2; t2 = 380;
M = letturaExcel(t1,t2);
figure
tiledlayout(2, 2);
nexttile
plot(M, 'b+');

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