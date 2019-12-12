%% read the image
clear;clc;
A = imread('img.png');
A = rgb2gray(A);
%imshow(A,[0,255]);
PDF=imhist(A)/(256^2); 
epsilon = 1*10^10;
rng(42);
%% Q2 - Uniform quantization
figure(1);
hold on;
figure(2);
hold on;
for b=1 : 8
   [MSE,rep,dec]=calcMSE(PDF,b);
   figure(1);
   scatter(b,MSE,'k','filled');
    
   figure(2);
   for i=1 : size(dec,2)
      scatter(b,dec(i),'r');
   end
   for i=1 : size(rep,2)
      scatter(b,rep(i),'b') 
   end
end

figure(1);
title('MSE as function of bit budget b');
xlabel('b');
ylabel('MSE');
hold off;

figure(2);
title('decision and representation levels for representive b values');
xlabel('b');
ylabel('representive and decision levels');
hold off;

%% Max_Lloyd-Q4
figure(3);
hold on;
figure(4);
hold on;
for b = 1 : 8
    dec = [1];
    delta = 256/(2^b);
    for i = 2 : (2^b)
        dec = [dec delta*(i-1)];
    end
    dec = [dec 256];
   [dec,rep,MSE] = Max_Loyd(PDF,dec,epsilon);
   
   figure(3);
   scatter(b,MSE,'k','filled');
   
   figure(4);
   for i=1 : size(dec,2)
      scatter(b,dec(i),'r');
   end
   for i=1 : size(rep,2)
      scatter(b,rep(i),'b') 
   end
end

figure(3);
title('Q4 - MSE values');
xlabel('b');
ylabel('MSE');
hold off;

figure(4);
title('Q4-decision and representation levels for representive b values');
xlabel('b');
ylabel('representive and decision levels');
hold off;

%% Q5 - random quantizations

figure(5);
hold on;
figure(6);
hold on;
for b = 1 : 8
    dec = sort(ceil((255-2).*rand(1,(2^b)-1)+2),'ascend');
    dec = [1 dec 256];
   [dec,rep,MSE] = Max_Loyd(PDF,dec,epsilon);
   
   figure(5);
   scatter(b,MSE,'k','filled');
   
   figure(6);
   for i=1 : size(dec,2)
      scatter(b,dec(i),'r');
   end
   for i=1 : size(rep,2)
      scatter(b,rep(i),'b') 
   end
end

figure(5);
title('Q5 - MSE values');
xlabel('b');
ylabel('MSE');
hold off;

figure(6);
title('Q5-decision and representation levels for representive b values');
xlabel('b');
ylabel('representive and decision levels');
hold off;
