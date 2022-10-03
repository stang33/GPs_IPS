ds = [1 1.5 2 2.5 3 3.5  4 4.5 5 5.5 6 6.5 7 7.5 8 ];% 22 24 29 30];
length(ds)

k = 4;%pts
M = 5; %boxes - 1, edge case is 


rder = EMM(ds, k, M)
