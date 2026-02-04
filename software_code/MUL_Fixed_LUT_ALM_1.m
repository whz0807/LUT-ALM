function mul = MUL_Fixed_LUT_ALM_1(e, A, B, n_h, n_l, v_h, v_l, app_k)

    f_o = v_h+v_l;

    len = length(A);
    len_h = 2^n_h;
    len_l = 2^n_l;
    
    f = zeros(1,len);
    f(1,:) = e - 1;

    lut_h = load(['LUT_ALM_luth_e' num2str(e) '_vh' num2str(v_h) '_vl' num2str(v_l) '_model1.mat']);
    lut_h = lut_h.init_lut_h_final;
    lut_l = load(['LUT_ALM_lutl_e' num2str(e) '_vh' num2str(v_h) '_vl' num2str(v_l) '_model1.mat']);
    lut_l = lut_l.init_lut_l_final;

    % 将每个查找表重构成四维数组（每个数组在一个维度上）
    lut_h_2d = reshape(lut_h, [1, len_h]);
    lut_l_2d = reshape(lut_l, [len_l, 1]);
    lut_total = lut_h_2d + lut_l_2d;
    
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
    
    % 将多维索引转换为线性索引
    max_dims = [len_h, len_l];
    % 计算p1和p2的线性索引
    p1_linear = p1_h * max_dims(2) + p1_l + 1;
    p2_linear = p2_h * max_dims(2) + p2_l + 1;
    % 直接使用线性索引赋值
    a_val = k1 + lut_total(p1_linear);
    b_val = k2 + lut_total(p2_linear);
    
    mul_val = a_val + b_val;
    
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
    
    
    