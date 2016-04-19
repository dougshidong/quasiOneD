clc;clear;
N = 100;
x = linspace(0, 10, N)';
y = sin(x);
noise = y + 0.2 * rand(length(y), 1);
y = y + noise;
plot(x,y,'-o')



k = 3;
m = N + k + 1;
L = m + 1;
T = [zeros(k+1, 1); linspace(0.2, 9.8, m - 2 *(k+1))' ;10 * ones(k+1, 1)]
for i = 1 : N - 1
    for j = 1 : L
        A(i,j) = getbij(x(i), j, k-1, T);
    end
end
A