%给定任意单粒子态求时间演化
J = 1;
U = 0;
% global dt L phi0

L = 8;
dt = 0.1;

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

% syms tt
t = 0:dt:10;
% G = zeros(L,L);

% phi0 = [1 0 1 0 1 0 1 0]';
phi0 = [1 0 1 0 1 0 1 0];

len = length(t);
num = zeros(L,len);

% figure;
% bart = bar(abs(num(:,1)));
% axis([0 9 0 3])

G = zeros(L,L);
tar = zeros(len,1);
for i = 1:L
    G(i,i) = phi0(i);
    num(i,1) = G(i,i);
    tar(1) = tar(1) + (-1)^i*num(i,1);
end

trans = expm(-1i*H1*dt);
for n = 2:len   
    G = trans*G*trans';
    
    for i = 1:L
        num(i,n) = G(i,i);
        tar(n) = tar(n) + (-1)^i*num(i,n);
    end
    
%     set(bart,'Ydata',abs(num(:,n)));
%     drawnow
%     pause(0.05)
end


figure;
plot(t,tar);
