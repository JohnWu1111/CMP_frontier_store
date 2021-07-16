clear all

clc

M = rand(10000,10000);            % 生成一个随机矩阵

tic

[A1,B1] = eig(M);                    % 求该随机矩阵的特征值和特征向量

t1=toc

 

tic

M = single(M);                     % 将数据转换为单精度型

M = gpuArray(M);                % 将数据从CPU中搬到GPU

[A2,B2] = eig(M);                 % 求特征值和特征向量

A2 = gather(A2);                 % 将数据从GPU中搬到CPU

t2 = toc