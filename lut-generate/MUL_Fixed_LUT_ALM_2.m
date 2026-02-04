function mul = MUL_Fixed_LUT_ALM_2(e, A, B, lut, n_h, n_l, f_o, app_k)

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
    
    a_0 = A - floor(A ./ 2) .* 2;
    a_1 = floor(A ./ 2) - floor(A ./ 4) .* 2;
    b_0 = B - floor(B ./ 2) .* 2;
    b_1 = floor(B ./ 2) - floor(B ./ 4) .* 2;
    a_t = floor(A ./ 4);
    b_t = floor(B ./ 4);
    at_abs = abs(a_t);
    bt_abs = abs(b_t);
    at_abs(A < 0) = at_abs(A < 0) - 1;
    bt_abs(B < 0) = bt_abs(B < 0) - 1;
    at_abs = at_abs .* 4;
    bt_abs = bt_abs .* 4;
    a_abs_1 = zeros(1,len);
    a_abs_1(A_s == 0) =  a_1(A_s == 0);
    a_abs_1(A_s == 1 & a_0 == 0) = 1;
    a_abs_1(A_s == 1 & a_1 == 0) = 1;
    b_abs_1 = zeros(1,len);
    b_abs_1(B_s == 0) =  b_1(B_s == 0);
    b_abs_1(B_s == 1 & b_0 == 0) = 1;
    b_abs_1(B_s == 1 & b_1 == 0) = 1;

    a_abs_0 = zeros(1,len);
    a_abs_0(A_s == 0) = a_0(A_s == 0);
    a_abs_0(A_s == 1 & a_1 == 0) = 1;
    a_abs_0(A_s == 1 & a_0 == 1) = 1;
    b_abs_0 = zeros(1,len);
    b_abs_0(B_s == 0) = b_0(B_s == 0);
    b_abs_0(B_s == 1 & b_1 == 0) = 1;
    b_abs_0(B_s == 1 & b_0 == 1) = 1;

    A_abs = at_abs + a_abs_1 .* 2 + a_abs_0;
    B_abs = bt_abs + b_abs_1 .* 2 + b_abs_0;

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
    
    x1 = a_val;
    x2 = b_val;
    e_x1 = floor(log2(abs(floor(x1)))) + 2;
    e_x2 = floor(log2(abs(floor(x2)))) + 2;
    e_x1(x1 == 0) = 1;
    e_x2(x2 == 0) = 1;
    e_x = e_x1;
    e_x(e_x2 > e_x1) = e_x2(e_x2 > e_x1);
    all_result = LOA(x1,x2,e+f_o+1,e_x,f_o,app_k);
    mul_val = all_result;
    
    mul_val = floor(mul_val .* (2.^f_o)) ./ (2.^f_o);
    mul_int = floor(mul_val);
    mul_frac = mul_val - mul_int;

    exp_tmp = 1 + mul_frac;
    exp_tmp = floor(exp_tmp.*(2.^f_o)) ./ (2.^f_o);


    mul = exp_tmp.*2.^(mul_int);

    d_0 = mul - floor(mul ./ 2) .* 2;
    d_1 = floor(mul ./ 2) - floor(mul ./ 4) .* 2;
    d_t = floor(mul ./ 4);
    dt_abs = d_t;
    dt_abs(mul_s == 1) = -1*d_t(mul_s == 1) -1;
    dt_abs = dt_abs .* 4;
    d_abs_1 = zeros(1,len);
    d_abs_1(mul_s == 0) =  d_1(mul_s == 0);
    d_abs_1(mul_s == 1 & d_0 == 0) = 1;
    d_abs_1(mul_s == 1 & d_1 == 0) = 1;
    d_abs_0 = zeros(1,len);
    d_abs_0(mul_s == 0) = d_0(mul_s == 0);
    d_abs_0(mul_s == 1 & d_1 == 0) = 1;
    d_abs_0(mul_s == 1 & d_0 == 1) = 1;
    mul = dt_abs + d_abs_1 .* 2 + d_abs_0;

    mul(mul_z == 1) = 0;
    
    
    