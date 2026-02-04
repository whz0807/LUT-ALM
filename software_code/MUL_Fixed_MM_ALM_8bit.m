function mul = MUL_Fixed_MM_ALM_8bit(e_in, A, B, v_h, v_l)
    len = length(A);
    e = e_in / 2;
    
    f = zeros(1,len);
    f(1,:) = e - 1;

    lut_h = load(['MM_ALM_LUT_' num2str(v_h) '.mat']);
    lut_h = lut_h.var_b;
    lut_l  = load(['MM_ALM_LUT_' num2str(v_l) '.mat']);
    lut_l = lut_l.var_b;

    A_s = zeros(1,len);
    A_s(A < 0) = 1;
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


    AH_abs = floor(A_abs ./ (2.^e));
    AL_abs = A_abs - AH_abs .* (2.^e);
    BH_abs = floor(B_abs ./ (2.^e));
    BL_abs = B_abs - BH_abs .* (2.^e);
    
    
    AH_z = zeros(1,len);
    AH_z(AH_abs == 0) = 1;
    BH_z = zeros(1,len);
    BH_z(BH_abs == 0) = 1;
    AL_z = zeros(1,len);
    AL_z(AL_abs == 0) = 1;
    BL_z = zeros(1,len);
    BL_z(BL_abs == 0) = 1;
    
    [k_ah,p_ah] = LOD(AH_abs,e+1,e,0);
    [k_bh,p_bh] = LOD(BH_abs,e+1,e,0);
    [k_al,p_al] = LOD(AL_abs,e+1,e,0);
    [k_bl,p_bl] = LOD(BL_abs,e+1,e,0);
    
    x_ah = p_ah ./ (2.^k_ah);
    x_bh = p_bh ./ (2.^k_bh);
    x_al = p_al ./ (2.^k_al);
    x_bl = p_bl ./ (2.^k_bl);
    x_ah = floor(x_ah .* (2.^f)) ./ (2.^f);
    x_bh = floor(x_bh .* (2.^f)) ./ (2.^f);
    x_al = floor(x_al .* (2.^f)) ./ (2.^f);
    x_bl = floor(x_bl .* (2.^f)) ./ (2.^f);
    k_ah(AH_abs == 0) = 0;
    x_ah(AH_abs == 0) = 0;
    k_bh(BH_abs == 0) = 0;
    x_bh(BH_abs == 0) = 0;
    k_al(AL_abs == 0) = 0;
    x_al(AL_abs == 0) = 0;
    k_bl(BL_abs == 0) = 0;
    x_bl(BL_abs == 0) = 0;

    xah_linear = x_ah .* (2.^f) + 1;
    xal_linear = x_al .* (2.^f) + 1;
    xbh_linear = x_bh .* (2.^f) + 1;
    xbl_linear = x_bl .* (2.^f) + 1;

    x_ah_log = lut_h(xah_linear);
    x_bh_log = lut_h(xbh_linear);
    x_al_log = lut_l(xal_linear);
    x_bl_log = lut_l(xbl_linear);
    x_ah_log = floor(x_ah_log .* (2.^(v_h))) ./ (2.^(v_h));
    x_bh_log = floor(x_bh_log .* (2.^(v_h))) ./ (2.^(v_h));
    x_al_log = floor(x_al_log .* (2.^(v_l))) ./ (2.^(v_l));
    x_bl_log = floor(x_bl_log .* (2.^(v_l))) ./ (2.^(v_l));
    
    ah_val = k_ah + x_ah_log;
    bh_val = k_bh + x_bh_log;
    al_val = k_al + x_al_log;
    bl_val = k_bl + x_bl_log;
    
    L_hh = ah_val + bh_val;
    L_ll = al_val + bl_val;
    L_hl = ah_val + bl_val;
    L_lh = al_val + bh_val;

    max_v = max(v_h,v_l);
    L_hl = floor(L_hl .* (2.^(max_v))) ./ (2.^(max_v));
    L_lh = floor(L_lh .* (2.^(max_v))) ./ (2.^(max_v));
    L_hh = floor(L_hh .* (2.^(v_h))) ./ (2.^(v_h));
    L_ll = floor(L_ll .* (2.^(v_l))) ./ (2.^(v_l));
    
    charac_hh = floor(L_hh);
    mant_hh = L_hh - charac_hh;
    mant_hh = floor(mant_hh.*(2.^(v_h))) ./ (2.^(v_h));
    charac_ll = floor(L_ll);
    mant_ll = L_ll - charac_ll;
    mant_ll = floor(mant_ll.*(2.^(v_l))) ./ (2.^(v_l));
    charac_hl = floor(L_hl);
    mant_hl = L_hl - charac_hl;
    mant_hl = floor(mant_hl.*(2.^(max_v))) ./ (2.^(max_v));
    charac_lh = floor(L_lh);
    mant_lh = L_lh - charac_lh;
    mant_lh = floor(mant_lh.*(2.^(max_v))) ./ (2.^(max_v));
    
    D_hh = floor((1 + mant_hh) .* 2.^(charac_hh));
    D_ll = floor((1 + mant_ll) .* 2.^(charac_ll));
    D_hl = floor((1 + mant_hl) .* 2.^(charac_hl));
    D_lh = floor((1 + mant_lh) .* 2.^(charac_lh));
    
    D_hh((AH_z == 1) | (BH_z == 1)) = 0;
    D_ll((AL_z == 1) | (BL_z == 1)) = 0;
    D_hl((AH_z == 1) | (BL_z == 1)) = 0;
    D_lh((AL_z == 1) | (BH_z == 1)) = 0;
    D_mid = D_hl + D_lh;
    D = D_hh .* (2.^(2.*e)) + D_ll + D_mid .* (2.^e);

    mul_s = zeros(1,len);
    mul_s(A_s ~= B_s) = 1;

    d_0 = D - floor(D ./ 2) .* 2;
    d_1 = floor(D ./ 2) - floor(D ./ 4) .* 2;
    d_t = floor(D ./ 4);
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
    