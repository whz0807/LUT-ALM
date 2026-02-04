%% lower-part-or adder: sum[n-1:0]=LOA(a[n-1:0],b[n-1:0],k)
%% describtion: sum[k-1:0]=a[k-1:0] | b[k-1:0], cin=a[k-1] & b[k-1], sum[n-1:k]=FA(a[n-1:k],b[n-1:k],cin)
%% an approximate adder to calculate A+B
%% N: the bit width of the adder/subtractor
%% e: the bit width of the integer part
%% f: the bit width of the fractional part
%% k: the bit width of approximate calculation

function y = LOA(A,B,N,e_in,f_in,k)

    len_A = length(A);
    len_f = length(f_in);
    len_e = length(e_in);
    f = ones(1,len_A);
    e = ones(1,len_A);
    if(len_f == len_A)
        f = f_in;
    else
        f = f * f_in;
    end

    if(len_e == len_A)
        e = e_in;
    else
        e = e * e_in;
    end

    atemp = floor(A .* (2.^(f)));
    btemp = floor(B .* (2.^(f)));
    
    %% deal with A/B < 0 
    atemp(A < 0) = atemp(A < 0) + 2.^(e(A < 0)+f(A < 0));
    btemp(B < 0) = btemp(B < 0) + 2.^(e(B < 0)+f(B < 0));

    atemp_acc = floor(atemp ./ (2.^k)) .* (2.^k);
    atemp_app = atemp - atemp_acc;
    btemp_acc = floor(btemp ./ (2.^k)) .* (2.^k);
    btemp_app = btemp - btemp_acc;
    
    if k <= 0
        Cin = zeros(1,len_A);
    else
        Cin = and(bitget(atemp_app,k),bitget(btemp_app,k)).*(2.^(k));
    end
    result_app = bitor(atemp_app,btemp_app);
    result_acc = Cin + atemp_acc + btemp_acc;
    Cout = result_acc + result_app;
    Cout(A < 0) = Cout(A < 0) - 2.^(e(A < 0) + f(A < 0));
    Cout(B < 0) = Cout(B < 0) - 2.^(e(B < 0)+f(B < 0));

    y = Cout ./ (2.^f);



