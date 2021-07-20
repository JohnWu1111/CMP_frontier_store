% calculation of (-1)^i*n_i, with time-dependent pertubation
tic;
J = 1;
U = 0.1;
L = [10 20 40];
dt = 0.01;
t = 0:dt:200;
lt = length(t);
G0 = [1 -1; -1 1]./2;
tar = zeros(length(L),lt);

ww = 1;
dd = 0.1;

for m = 1:length(L)
    k = zeros(1,L(m)/2);
    Hk = zeros(1,L(m)/2);

    for i = 1:(L(m)/2)
        k(i) = 2*pi*i/L(m) - pi/2;
        Hk(i) = 2*J*cos(k(i));
    end
    
    for n = 1:lt
        for i = 1:(L(m)/2)
            coeff = dd*cos(ww*t(n));
            H = [-Hk(i) coeff; coeff Hk(i)];
            Tev = expm(-1i*H*t(n));
            G = Tev'*G0*Tev;
            tar(m,n) = tar(m,n) + G(1,2) + G(2,1);
        end    
    end

%     plot(t,real(tar(m,:))./(L(m)/2))
    [w,tarw] = Fourier(t,real(tar(m,:))./(L(m)/2));
    tarw = tarw./sum(tarw);
    plot(w(1:200),tarw(1:200))
    hold on;
end
% plot(t,-besselj(0,4*J*t))
[w,tarw] = Fourier(t,-besselj(0,4*J*t));
tarw = tarw./sum(tarw);
plot(w(1:200),tarw(1:200))
legend('L=10','L=20','L=40','L=\infty')
toc;

function [omega,y] = Fourier(t,x)
    len = length(t);
    T = t(end) - t(1);
%     dt = T/(len-1);
%     Omega = 2*pi/dt;
    domega = 2*pi/T;
    omega0 = 0;
    omega = zeros(len,1);
    y = zeros(len,1);
    for i = 1:len
        omega(i) = (i-1)*domega + omega0;
        for j = 1:len
            y(i) = y(i) + exp(-1i*omega(i)*t(j))*x(j);
        end 
        y(i) = abs(y(i));
    end
end