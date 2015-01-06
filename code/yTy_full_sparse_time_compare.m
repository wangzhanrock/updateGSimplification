format long;
close all, clear all;
TwoD = 1;
if(TwoD)
    LLM=dlmread('LLM_s.txt');
    W=dlmread('V_s.txt');
else
    LLM=dlmread('data3D/LLM.txt');
    W=dlmread('data3D/V.txt');
end

yFull= tril(LLM)\W;
ySparse = sparse(yFull),

isequal(yFull, ySparse),

t = cputime;
    gFull= yFull'*yFull;
eFull = cputime -t,


t = cputime;
    gSparse= ySparse'*ySparse;
eSparse = cputime -t,

isequal(gSparse, gFull),
