J = 1;
U = 0.1;
L = [10 20 40];
dt = 0.01;
t = 0:dt:20;
lt = length(t);
G0 = [1 -1; -1 1]./2;
tar = zeros(length(L),lt);

for m = 1:length(L)
    k = zeros(1,L(m)/2);
    Hk = zeros(1,L(m)/2);
    Tev = zeros(2,2,L(m)/2);

    for i = 1:(L(m)/2)
        k(i) = 2*pi*i/L(m) - pi/2;
        Hk(i) = 2*J*cos(k(i));
    end
    
    for i = 1:(L(m)/2)
        H = [-Hk(i) 0; 0 Hk(i)];
        Tev(:,:,i) = expm(-1i*H*dt);
        tar(m,1) = tar(m,1) + G0(1,2) + G0(2,1);
    end
    Tevt = Tev;
    
    for n = 2:lt
        for i = 1:(L(m)/2)
            H = [-Hk(i) 0; 0 Hk(i)];
            Tevt(:,:,i) = Tevt(:,:,i)*Tev(:,:,i);
            G = Tevt(:,:,i)'*G0*Tevt(:,:,i);
            tar(m,n) = tar(m,n) + G(1,2) + G(2,1);
        end    
    end

    plot(t,real(tar(m,:))./(L(m)/2))
    hold on;
end
plot(t,-besselj(0,4*J*t))
legend('L=10','L=20','L=40','L=\infty')