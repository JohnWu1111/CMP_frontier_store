%给定任意单粒子态求时间演化
J = 1;
U = 0.1;
% global dt L phi0

L = 8;
dt = 0.01;

H1 = zeros(L,L);
for i = 1:L-1
    H1(i,i+1) = -J;
    H1(i+1,i) = -J;
    H1(i,i) = -U;
end

H1(L,1) = -J;
H1(1,L) = -J;
H1(L,L) = -U;

[V,D] = eig(H1);

H2 = zeros(L,L);
k = zeros(1,L);
temp = zeros(1,L);

for i = 1:L
    k(i) = 2*pi*i/L;
    temp(i) = -2*J*cos(k(i))-U;
end

for i = 1:L
    for j = i:L
        if (temp(i)>temp(j))
            tem = temp(i);
            temp(i) = temp(j);
            temp(j) = tem;
            tem = k(i);
            k(i) = k(j);
            k(j) = tem;
        end
    end
end

% syms tt
t = 0:dt:10;
% G = zeros(L,L);

% phi0 = [1 0 1 0 1 0 1 0]';
phi0 = [1 0 1 0 1 0 1 0]'/sqrt(4);
coeff = inv(V)*phi0;
coefft = zeros(L,1);

len = length(t);
num = zeros(L,len);

% figure;
% bart = bar(abs(num(:,1)));
% axis([0 9 0 3])

for n = 1:len
    for i = 1:L
        coefft(i) = coeff(i)*exp(-1i*D(i,i)*t(n));
    end
    phit = V*coefft;

    G0 = zeros(L,L);
    for i = 1:L
        trans = zeros(L,L);
        trans(i,i) = 1; 
        G0(i,i) = phit'*trans*phit;
    end
    
    for i = 1:L
        num(i,n) = G0(i,i);
    end
    
%     set(bart,'Ydata',abs(num(:,n)));
%     drawnow
%     pause(0.05)
end

figure;
plot(t,num(1,:));
