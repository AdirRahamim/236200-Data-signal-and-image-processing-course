function [dec_f,rep,MSE] = Max_Loyd(PDF,dec,epsilon)
    rep = zeros(1,size(dec,2)-1);
    MSE_old = 0;
    MSE_curr = 0;
    val = 1:1:256;
    while(1)
        %calculating optimal representation levels
        for i = 1 : size(rep,2)
            rep(i)=sum(val(dec(i):dec(i+1)).*PDF(dec(i):dec(i+1))');
            rep(i) = rep(i)/sum(PDF(dec(i):dec(i+1)));
        end
        %calculate optimal decesion levels
        for i =2 : size(dec,2)-1
            dec(i) = round((rep(i)+rep(i-1))/2);
        end
        
        %calculat new MSE and comapre to old one
        for i=1 : size(rep,2)
           for k=dec(i) : dec(i+1) 
               MSE_curr=MSE_curr+PDF(k).*(k-rep(i)).^2;
           end
        end
        if abs(MSE_old-MSE_curr)<epsilon
            break;
        end
        MSE_old = MSE_curr;
        MSE_curr = 0;
    end
    MSE = MSE_curr;
    dec_f = dec;
end
