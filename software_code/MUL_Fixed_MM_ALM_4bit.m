function mul = MUL_Fixed_MM_ALM_4bit(e, A, B, v)
    len = length(A);
    
    f = zeros(1,len);
    f(1,:) = e - 1;

    lut = load(['MM_ALM_LUT_',num2str(v),'.mat']);
    lut = lut.var_b;

    A_z = zeros(1,len);
    A_z(A == 0) = 1;
    A_s = zeros(1,len);
    A_s(A < 0) = 1;
    B_z = zeros(1,len);
    B_z(B == 0) = 1;
    B_s = zeros(1,len);
    B_s(B < 0) = 1;

    mul_s = zeros(1,len);
    mul_s(A_s ~= B_s) = 1;
    mul_z = zeros(1,len);
    mul_z(A_z==1 | B_z==1) = 1;

    A_abs = abs(A);
    B_abs = abs(B);

    [k1,p1] = LOD(A_abs,e+1,e,0);
    [k2,p2] = LOD(B_abs,e+1,e,0);

    q1 = p1 ./ (2.^k1);
    q2 = p2 ./ (2.^k2);
    q1 = floor(q1 .* (2.^f)) ./ (2.^f);
    q2 = floor(q2 .* (2.^f)) ./ (2.^f);
    k1(A == 0) = 0;
    k2(B == 0) = 0;
    q1(A == 0) = 0;
    q2(B == 0) = 0;

    q1_linear = q1 .* (2.^f) + 1;
    q2_linear = q2 .* (2.^f) + 1;
    % 直接使用线性索引赋值
    a_val = k1 + lut(q1_linear);
    b_val = k2 + lut(q2_linear);
    a_val(A == 0) = 0;
    b_val(B == 0) = 0;
    
    mul_val = a_val + b_val;
    
    mul_val = floor(mul_val .* (2.^(v))) ./ (2.^(v));
    mul_int = floor(mul_val);
    mul_frac = mul_val - mul_int;
    
    exp_tmp = 1 + mul_frac;
    exp_tmp = floor(exp_tmp.*(2.^v)) ./ (2.^v);
    
    mul = floor(not(mul_z).*exp_tmp.*2.^(mul_int));
    
    mul = (-1).^mul_s.*mul;

    
    
    