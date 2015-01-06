close all, clear all;
TwoD = 1;
if(TwoD)
    LLM=dlmread('LLM_s.txt');
    W=dlmread('V_s.txt');
else
    LLM=dlmread('data3D/LLM.txt');
    W=dlmread('data3D/V.txt');
end

size(LLM);
size(W);

spy(LLM);

figure(2);
spy(W);

M = tril(LLM)*triu(LLM);
x=M\W;
g= W' * x;
figure(3);
spy(g),

y= tril(LLM)\W;
g2= y'*y;
figure(4);
spy(g2);



