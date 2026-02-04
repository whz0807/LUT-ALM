function mul = MUL_Fixed_LUT_ALM_0(e, A, B, lut, n_h, n_l, f_o, app_k)

    len = length(A);
    len_h = 2^n_h;
    len_l = 2^n_l;
    
    f = zeros(1,len);
    f(1,:) = e - 1;
    
    A_z = zeros(1,len);
    A_z(A == 0) = 1;
    A_s = zeros(1,len);
    A_s(A < 0) = 1;
    B_z = zeros(1,len);
    B_z(B == 0) = 1;
    B_s = zeros(1,len);
    B_s(B < 0) = 1;
    
    A_abs = abs(A);
    B_abs = abs(B);
    
    [k1,p1] = LOD(A_abs,e+1,e,0);
    [k2,p2] = LOD(B_abs,e+1,e,0);
    
    k1(A_abs == 0) = 0;
    p1(A_abs == 0) = 0;
    k2(B_abs == 0) = 0;
    p2(B_abs == 0) = 0;
    
    q1 = p1 ./ (2.^k1);
    q2 = p2 ./ (2.^k2);
    q1 = floor(q1 .* (2.^f)) ./ (2.^f);
    q2 = floor(q2 .* (2.^f)) ./ (2.^f);
    
    q1_int = floor(q1 .* (2.^f));
    q2_int = floor(q2 .* (2.^f));
    p1_h = floor(q1_int ./ (2.^(f-n_h)));
    p1_l = floor(q1_int ./ (2.^(f-n_h-n_l))) - p1_h .* 2.^n_l;
    p2_h = floor(q2_int ./ (2.^(f-n_h)));
    p2_l = floor(q2_int ./ (2.^(f-n_h-n_l))) - p2_h .* 2.^n_l;

    mul_s = zeros(1,len);
    mul_s(A_s ~= B_s) = 1;
    mul_z = zeros(1,len);
    mul_z(A_z==1 | B_z==1) = 1;
    
    a_val = zeros(1,len);
    b_val = zeros(1,len);
    for i = 0:len_h-1
        for j = 0:len_l-1
            a_val((p1_h == i) & (p1_l == j)) = k1((p1_h == i) & (p1_l == j)) + lut(j+1,i+1);
            b_val((p2_h == i) & (p2_l == j)) = k2((p2_h == i) & (p2_l == j)) + lut(j+1,i+1);
        end
    end
    
    mul_val = a_val + b_val;
    
    mul_val = floor(mul_val .* (2.^f_o)) ./ (2.^f_o);
    mul_int = floor(mul_val);
    mul_frac = mul_val - mul_int;
    
    exp_tmp = 1 + mul_frac;
    exp_tmp = floor(exp_tmp.*(2.^f_o)) ./ (2.^f_o);
    
    mul = floor(not(mul_z).*exp_tmp.*2.^(mul_int));
    mul(mul_s == 1) = mul(mul_s == 1) * (-1);
    
    
    