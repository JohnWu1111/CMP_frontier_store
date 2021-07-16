%对于单自由费米子，分别从本征态、能量、关联函数的角度证明哈密顿量对角化等价于傅里叶变换。
J = 1;
U = 0.1;
L = 8;

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

G = zeros(L,L);
for i = 1:L
    for j = 1:L
        trans = zeros(L,L);
        trans(i,j) = 1;
        G(i,j) = V(:,1)'*trans*V(:,1);
    end
end

for i = 1:L
    H2(i,i) = temp(i);
end

V2 = zeros(L,L);
for i = 1:L
    for j = 1:L
        V2(i,j) = exp(-1i*k(i)*j);
    end
end

diff = D - H2;