J = 1;
U = 0.1;
global dt L phi0

L = 8;
dt = 0.01;

H1 = zeros(L,L);
for i = 1:L-1
    H1(i,i+1) = -J;
    H1(i+1,i) = -J;
    H1(i,i) = -U;
end

H1(L,1) = J;
H1(1,L) = J;
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
% G = zeros(L,L);
G0 = zeros(L,L);
% phi0 = [1 0 1 0 1 0 1 0]';
phi0 = [1 1 1 0 1 0 1 0]'./sqrt(5);
for i = 1:L
    for j = 1:L
        trans = zeros(L,L);
        trans(i,j) = 1; 
        G0(i,j) = phi0'*trans*phi0;
    end
end

t = 0:dt:10;
G = runge(t,G0);

function y = f(x,~,m,n)
    global phi0
    temp = 1i*phi0'*(x(mod(m,8)+1,n)+x(mod(m-2,8)+1,n))*phi0;  
    y = temp;
end

function y = runge(t,y0) 
    global dt L
    n = length(t);
    y = zeros(L,L,n);
    y(:,:,1) = y0;
    
    for i = 1:n-1
        for j = 1:L
            for k = 1:L
                c1 = f(y(:,:,i),t(i),j,k);
                c2 = f(y(:,:,i)+c1.*dt/2,t(i)+dt/2,j,k);
                c3 = f(y(:,:,i)+c2.*dt/2,t(i)+dt/2,j,k);
                c4 = f(y(:,:,i)+c3.*dt,t(i)+dt,j,k);
                y(j,k,i+1) = y(j,k,i) + dt*(c1+2*c2+2*c3+c4)/6;
            end
        end
    end
end