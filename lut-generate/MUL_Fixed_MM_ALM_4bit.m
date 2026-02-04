%% An accurate fixed number multiplier for multiple precision
%% Inputs:
%% e: the bit width of input
%% A,B: input operands of the multiplier/divier
%% ctrl_bits: control the width of multiplier, ctrl_bits = 1'b1 for 8bits, ctrl_bits = 1'b0 for 4 bits
%% Outputs:
%% mul_div: approximate output product/quotient in decimal

function [mul_div] = MUL_Fixed_MM_ALM_4bit(e, A, B, var, m_d, v)

len = length(A);

f = e-1;

A_z = zeros(1,len);
A_z(A == 0) = 1;
B_z = zeros(1,len);
B_z(B == 0) = 1;

mul_z = zeros(1,len);
mul_z(A_z==1 | B_z==1) = 1;


[k1,p1] = LOD(A,e+1,e,0);
[k2,p2] = LOD(B,e+1,e,0);

q1 = p1 ./ (2.^k1);
q2 = p2 ./ (2.^k2);
q1 = floor(q1 * (2^f)) / (2^f);
q2 = floor(q2 * (2^f)) / (2^f);
k1(A == 0) = 0;
k2(B == 0) = 0;
q1(A == 0) = 0;
q2(B == 0) = 0;

a_val = k1;
a_val(q1 == 0) = k1(q1 == 0) + var(1);
a_val(q1 == (1/8)) = k1(q1 == (1/8)) + var(2);
a_val(q1 == (2/8)) = k1(q1 == (2/8)) + var(3);
a_val(q1 == (3/8)) = k1(q1 == (3/8)) + var(4);
a_val(q1 == (4/8)) = k1(q1 == (4/8)) + var(5);
a_val(q1 == (5/8)) = k1(q1 == (5/8)) + var(6);
a_val(q1 == (6/8)) = k1(q1 == (6/8)) + var(7);
a_val(q1 == (7/8)) = k1(q1 == (7/8)) + var(8);
b_val = k2;
b_val(q2 == 0) = k2(q2 == 0) + var(1);
b_val(q2 == (1/8)) = k2(q2 == (1/8)) + var(2);
b_val(q2 == (2/8)) = k2(q2 == (2/8)) + var(3);
b_val(q2 == (3/8)) = k2(q2 == (3/8)) + var(4);
b_val(q2 == (4/8)) = k2(q2 == (4/8)) + var(5);
b_val(q2 == (5/8)) = k2(q2 == (5/8)) + var(6);
b_val(q2 == (6/8)) = k2(q2 == (6/8)) + var(7);
b_val(q2 == (7/8)) = k2(q2 == (7/8)) + var(8);
a_val(A == 0) = 0;
b_val(B == 0) = 0;

if m_d == 1
    mul_val = a_val - b_val;
else
    mul_val = a_val + b_val;
end

mul_val = floor(mul_val .* (2.^(v))) ./ (2.^(v));
mul_int = floor(mul_val);
mul_frac = mul_val - mul_int;

mul_div = floor(not(mul_z).*(1 + mul_frac).*2.^(mul_int));

