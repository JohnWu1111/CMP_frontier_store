J = 1;
U = 0.01;
global L

L = 8;
dt = 0.01;
len = 2^L;
t = 0:dt:10;
omega = 0.1;

I = eye(2);
K = [1 0; 0 -1];
C = [0 0; 1 0];
C_dag = [0 1; 0 0];

phi_pos = [1 0 1 0 1 0 1 0 1 0];
phi_s = zeros(2,2);
phi_s(:,1) = [1 0]';
phi_s(:,2) = [0 1]';

phi0 = kron(phi_s(:,phi_pos(1)+1),phi_s(:,phi_pos(2)+1));
for i = 3:L
    phi0 = kron(phi0,phi_s(:,phi_pos(i)+1));
end
phit = phi0;

lt = length(t);
num = zeros(L,lt);

for n = 1:lt
    H = zeros(len,len);
    % pos=1单独赋值
    H1 = C_dag; 
    H1 = kron(H1,C);
    for j = 3:L
        H1 = kron(H1,I);
    end
    H = H -J.*(H1 + H1');

    for i = 2:L-1
        H1 = I;
        for j = 2:i-1
            H1 = kron(H1,I);
        end
        H1 = kron(H1,C_dag);
        H1 = kron(H1,C);
        for j = i+2:L
            H1 = kron(H1,I);
        end
        H = H -J.*(H1 + H1');
    end

    % pos=L单独赋值
    H1 = C_dag; 
    for j = 2:L-1
        H1 = kron(H1,K);
    end
    H1 = kron(H1,C);
    H = H -J.*(H1 + H1');

    % pos=1单独赋值
    H2 = C_dag*C;
    for j = 2:L
        H2 = kron(H2,I);
    end
    H = H - (-1)^i*cos(omega*dt*n)*U.*H2;

    for i = 2:L
        H2 = I;
        for j = 2:i-1
            H2 = kron(H2,I);
        end
        H2 = kron(H2,C_dag*C);
        for j = i+1:L
            H2 = kron(H2,I);
        end
        H = H - (-1)^i*cos(omega*dt*n)*U.*H2;
    end

    [V,D] = eig(H);

    coeff = inv(V)*phit;
    coefft = zeros(len,1);

    % figure;
    % bart = bar(abs(num(:,1)));
    % axis([0 9 0 3])
    ten = zeros(L,len);
    for i = 1:len
        ten(:,i) = tentotwo(i);
    end

    for i = 1:len
        coefft(i) = coeff(i)*exp(-1i*D(i,i)*dt);
    end
    phit = V*coefft;

    G0 = zeros(len,1);
    
    for i = 1:len
        G0(i) = conj(phit(i))*phit(i);
    end
    
    for i = 1:len        
        num(:,n) = num(:,n) + ten(:,i)*G0(i);
    end

%     set(bart,'Ydata',abs(num(:,n)));
%     drawnow
%     pause(0.05)
end

for i = 1:L
    figure;
    plot(t,num(i,:));
    xlabel('time')
    str = strcat('density of number at x=',num2str(i));
    ylabel(str)
    str = strcat('L=',num2str(L),',pos=',num2str(i));
    title(str)
    str = strcat('L=',num2str(L),'_pos=',num2str(i));
    fname = [str,'.png '];
    saveas(gcf, fname, 'png')
end

function y = tentotwo(n)
    global L
    y = zeros(L,1);
    m = n-1;
    for i = L:-1:1
        y(i) = mod(m,2);
        m = (m - y(i))/2;
    end
end
