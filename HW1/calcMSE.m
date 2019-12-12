%%calculate MSE for uniform quantization
function [MSE,rep,dec] = calcMSE(PDF,b)
    rep_num=2^b;
    delta=256/rep_num;
    MSE=0;
    dec = [];
    rep = [];
    val = 1:1:256;
    for i=1 :rep_num
        dec = [dec (i-1)*delta+1];
        Q_x = sum(val((i-1)*delta+1 : i*delta))/delta;
        rep = [rep Q_x];
        for j=(i-1)*delta+1 : i*delta
            MSE= MSE+(PDF(j)*(j-Q_x).^2);
        end
    end
    dec = [dec i*delta];
end