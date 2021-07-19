% 1D repulsive spinless Ferimions, half-filling

% Method1: ED in full space
J = 1;
U = 0.1;
L = 8;
len = 2^L;

H = zeros(len,len);
I = eye(2);
K = [1 0; 0 -1];
C_dag = [0 0; 1 0];
C = [0 1; 0 0];
nn = [-1/2 0; 0 1/2];

% pos=1单独赋值
H1 = C_dag; 
H1 = kron(H1,C);
H2 = nn;
H2 = kron(H2,nn);
for j = 3:L
    H1 = kron(H1,I);
    H2 = kron(H2,I);
end
H = H -J.*(H1 + H1')./2 + U.*H2;

for i = 2:L-1
    H1 = I;
    H2 = I;
    for j = 2:i-1
        H1 = kron(H1,I);
        H2 = kron(H2,I);
    end
    H1 = kron(H1,C_dag);
    H1 = kron(H1,C);
    H2 = kron(H2,nn);
    H2 = kron(H2,nn);
    for j = i+2:L
        H1 = kron(H1,I);
        H2 = kron(H2,I);
    end
    H = H -J.*(H1 + H1')./2 + U.*H2;
end

% pos=L单独赋值
H1 = C_dag;
H2 = nn;
for j = 2:L-1
    H1 = kron(H1,K);
    H2 = kron(H2,I);
end
H1 = kron(H1,C);
H2 = kron(H2,nn);
H = H -J.*(H1 + H1')./2 + U.*H2;

e1 = eig(H);

e1(1)


% Method 2: ED in N = L/2 sub-space

len_s = factorial(L)/(factorial(L/2)^2);
H3 = zeros(len_s,len_s);
store1 = zeros(L,len_s);
store2 = zeros(len,1);
count = 1;
for i = 1:len
    temp = tentotwo(i,L);
    temp_num = sum(temp);
    if temp_num == L/2
        store1(:,count) = temp;
        store2(i) = count;
        count = count +1;
    end  
end

for k = 1:len_s
    phi_k = store1(:,k);
    for i = 1:L-1
        j = i+1;
        if phi_k(i) == phi_k(j)
            H3(k,k) = H3(k,k) + U/4;
        else
            H3(k,k) = H3(k,k) - U/4;
            phi_b = phi_k;
            phi_b(j) = phi_k(i);
            phi_b(i) = phi_k(j);
            b = twototen(phi_b);
            b = store2(b);
            H3(k,b) = H3(k,b) - J/2;
        end
    end
    if phi_k(1) == phi_k(8)
        H3(k,k) = H3(k,k) + U/4;
    else
        H3(k,k) = H3(k,k) - U/4;
        phi_b = phi_k;
        phi_b(1) = phi_k(8);
        phi_b(8) = phi_k(1);
        b = twototen(phi_b);
        b = store2(b);
        H3(k,b) = H3(k,b) + J/2;
    end
end

e2 = eig(H3);
e2(1)

check = zeros(len,1);
for i = 1:len
    if twototen(tentotwo(i,L)) == i
        check(i) = 1;
    end
end

function y = twototen(phi)
    len = length(phi);
    y = 0;
    for i = 1:len
        y = 2*y+phi(i);
    end
    y = y+1;
end

function y = tentotwo(n,L)
    y = zeros(L,1);
    m = n-1;
    for i = L:-1:1
        y(i) = mod(m,2);
        m = (m - y(i))/2;
    end
end
