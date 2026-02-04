function mul = MUL_Fixed_LUT_ALM(mode, e, A, B, n_h, n_l, v_h, v_l, t)
    func_table = {
        @MUL_Fixed_LUT_ALM_0
        @MUL_Fixed_LUT_ALM_1
        @MUL_Fixed_LUT_ALM_2
    };

    if mode < 0 || mode >= numel(func_table)
        error('Unsupported mode: %d', mode);
    end

    mul  = func_table{mode+1}(e, A, B, n_h, n_l, v_h, v_l, t);
end
