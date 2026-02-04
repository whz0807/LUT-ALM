%% leading one detector
%% N: the bit width of the adder/subtractor
%% e: the bit width of the integer part
%% f: the bit width of the fractional part

function [k,q] = LOD(x,N,e,f)
    k = floor(log2(x));
    q = x - 2.^k;
