%% the code for generating the optimized LUTs for 8-bit LUT-ALM 
%% LUT-ALM is the logarithmic multiplier presented in the paper "LUT-ALM: Trading Off Accuracy and Power for Approximate Logarithmic Multipliers via LUT Optimization"

clear
clc

e=8; % the bit width of the fixed-point multiplication
f=e-1; % the fractional bits of the fixed-point multiplication
v_h = 3; v_l = 4; v=v_h+v_l; % the widths of LUT_{v_h} and LUT_{v_l}
% for 8-bit LUT-ALM, we set n_h=n_l=3
n_h = 3; n_l = 3; d_h = 2^n_h; d_l = 2^n_l; % the sizes of LUT_{v_h} and LUT_{v_l}
mode = 2; % the mode of the multiplier, 0: LUT-ALM_0, 1: LUT-ALM_1, 2: LUT-ALM_2
t = 4; % the approximate bits in inexact adders

% initialize the lut
lut_h = (0:1:(2^(n_h)-1)) / (2^(n_h)) ; 
lut_l = (0:1:(2^(n_l)-1)) / (2^(n_h+n_l)); 
lut_h = floor(lut_h * 2^v) / 2^v;
lut_l = floor(lut_l * 2^v) / 2^v;

init_lut_h = lut_h;
init_lut_l = lut_l;
init_lut_h_final = init_lut_h;
init_lut_l_final= init_lut_l;

MED = 1e18;
MRED = 99999;
MED_be = 0;

max_iter = 10000;  % 最大迭代次数

% the input operands for 8-bit LUT-ALM
[X,Y] = meshgrid(-2^(e-1):(2^(e-1)-1), -2^(e-1):(2^(e-1)-1));
X_T = reshape(X, 1, []);
Y_T = reshape(Y, 1, []);
A = X_T;
B = Y_T;
MUL_acc = A .* B;
len = length(A);

for i = 1:max_iter
    fprintf('i %d has changed.\n', i);
    if(MED_be == MED)
        save(['LUT_ALM_MED_e' num2str(e) '_vh' num2str(v_h) '_vl' num2str(v_l) '_model' num2str(mode) '.mat'], 'MED')
        save(['LUT_ALM_luth_e' num2str(e) '_vh' num2str(v_h) '_vl' num2str(v_l) '_model' num2str(mode) '.mat'], 'init_lut_h_final')
        save(['LUT_ALM_lutl_e' num2str(e) '_vh' num2str(v_h) '_vl' num2str(v_l) '_model' num2str(mode) '.mat'], 'init_lut_l_final')
        break;
    else
        MED_be = MED;
        for num = 1:1:2^n_h
            max_lut_l = max(init_lut_l);
            for j=0:2^(-v):(1-2^(-v)-max_lut_l)
                tmp = init_lut_h(num);
                a = floor(j * 2^v_h) / 2^v_h;
                init_lut_h(num)=a;

                lut_h = reshape(init_lut_h, [1, d_h]);
                lut_l = reshape(init_lut_l, [d_l, 1]);
                % 利用广播机制相加
                lut_total = lut_h + lut_l;

                mul = MUL_Fixed_LUT_ALM(mode, e, A, B, lut_total, n_h, n_l, v, t);

                ED_h = mul-MUL_acc;
                RED_h = zeros(1,len);
                RED_h(MUL_acc~=0) = ED_h(MUL_acc~=0)./MUL_acc(MUL_acc~=0);
                MED_h = mean(abs(ED_h(:)));
                MRED_h = mean(abs(RED_h));
                
                % computr the difference of MED
                delta_MED_h = MED_h - MED;
                if((delta_MED_h < 0))
                    MED = MED_h;
                    MRED = MRED_h;
                    init_lut_h_final = init_lut_h;
                    init_lut_l_final = init_lut_l;
                    fprintf('v_h: Num %d has changed.\n', num);
                    fprintf('MED %d has changed.\n', MED);
                    fprintf('MRED %d has changed.\n', MRED);
                else
                    init_lut_h(num) = tmp;
                    init_lut_h_final = init_lut_h;
                    init_lut_l_final= init_lut_l;
                end

            end
        end

        for num = 1:1:2^n_l
            max_lut_h = max(init_lut_h);
            for j=0:2^(-v):(1-2^(-v)-max_lut_h)
                tmp = init_lut_l(num);
                a = floor(j * (2^v)) / (2^v);
                a = floor(a * (2^v)) / (2^v) - floor(a * 2^v_h) / 2^v_h;
                init_lut_l(num)=a;

                lut_h = reshape(init_lut_h, [1, d_h]);
                lut_l = reshape(init_lut_l, [d_l, 1]);
                % 利用广播机制相加
                lut_total = lut_h + lut_l;

                mul = MUL_Fixed_LUT_ALM(mode, e, A, B, lut_total, n_h, n_l, v, t);

                ED_l = mul-MUL_acc;
                RED_l = zeros(1,len);
                RED_l(MUL_acc~=0) = ED_l(MUL_acc~=0)./MUL_acc(MUL_acc~=0);
                MED_l = mean(abs(ED_l(:)));
                MRED_l = mean(abs(RED_l));
                
                % computr the difference of MED
                delta_MED_l = MED_l - MED;
                if((delta_MED_l < 0))
                    MED = MED_l;
                    init_lut_h_final = init_lut_h;
                    init_lut_l_final= init_lut_l;
                    fprintf('v_l: Num %d has changed.\n', num);
                    fprintf('MED %d has changed.\n', MED);
                    fprintf('MRED %d has changed.\n', MRED);
                else
                    init_lut_l(num) = tmp;
                    init_lut_h_final = init_lut_h;
                    init_lut_l_final= init_lut_l;
                end

            end
        end
        if(i == max_iter)
            save(['LUT_ALM_MED_e' num2str(e) '_vh' num2str(v_h) '_vl' num2str(v_l) '_model' num2str(mode) '.mat'], 'MED')
            save(['LUT_ALM_luth_e' num2str(e) '_vh' num2str(v_h) '_vl' num2str(v_l) '_model' num2str(mode) '.mat'], 'init_lut_h_final')
            save(['LUT_ALM_lutl_e' num2str(e) '_vh' num2str(v_h) '_vl' num2str(v_l) '_model' num2str(mode) '.mat'], 'init_lut_l_final')
        end
    end
end
fprintf('Iteration ends');