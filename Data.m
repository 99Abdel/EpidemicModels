data = xlsread('MyData.xlsx', 'Foglio1', 'B1:B123');

D = data(:,1)

t = 1:1:122

plot(t,D)
