%% the code for generating the optimized LUTs for MM-ALM
%% MM-ALM is the logarithmic multiplier presented in the paper "LUT-ALM: Trading Off Accuracy and Power for Approximate Logarithmic Multipliers via LUT Optimization"

clear all;

e=4;
m_d = 0;

[X,Y] = meshgrid(0:(2^(e)-1), 0:(2^(e)-1));
X_T = reshape(X, 1, []);
Y_T = reshape(Y, 1, []);
A = X_T;
B = Y_T;
MUL_acc = A .* B;

t = 3;
v = 4;
str1 = strcat('MM_ALM_LUT_',num2str(v-1),'.mat');
str2 = strcat('MM_ALM_MED_',num2str(v-1),'.mat');
lut = load(str1);
med = load(str2);
var_n_1 = lut.var_b;
min_MED = med.min_MED;
var_b = var_n_1;

for i = 0:(t^8-1)
    bin_i = dec2base(i,t,8);
    dec_i = zeros(1,8);
    dec_i(1) = double(bin_i(1)) - double('0');
    dec_i(2) = double(bin_i(2)) - double('0');
    dec_i(3) = double(bin_i(3)) - double('0');
    dec_i(4) = double(bin_i(4)) - double('0');
    dec_i(5) = double(bin_i(5)) - double('0');
    dec_i(6) = double(bin_i(6)) - double('0');
    dec_i(7) = double(bin_i(7)) - double('0');
    dec_i(8) = double(bin_i(8)) - double('0');

    var = zeros(1,8);
    var = var + dec_i - 1;

    var_n = var_n_1 + var .* 2.^(-v);
    var_n(var_n < 0) = 0; 
    var_n(var_n > (1-2^(-v))) = (1-2^(-v)); 

    % transform A and B to the logarithmic numbers
    len = length(A);

    [mul] = MUL_Fixed_MM_ALM_4bit(e, A, B, var_n, 0, v);

    ED = mul-MUL_acc;
    MED = mean(abs(ED));

    if(MED < min_MED)
        min_MED = MED;
        var_b = var_n;
    end

end

save(strcat('MM_ALM_LUT_',num2str(v),'.mat'),'var_b');
save(strcat('MM_ALM_MED_',num2str(v),'.mat'),'min_MED');



