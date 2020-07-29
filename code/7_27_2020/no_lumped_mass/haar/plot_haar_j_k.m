clc;
clear;

syms x

j_index = [0,1,2];
k_index = [0,1,2];


basis = [];


for i = 1:length(k_index)     

    haar_j = 1;
    haar_k = k_index(i);

    haar_one_begin = haar_k/(2^haar_j);
    haar_one_end = (haar_k+1/2)/(2^haar_j);
    haar_minus_one_begin = (haar_k+1/2)/(2^haar_j);
    haar_minus_one_end = (haar_k+1)/(2^haar_j);
    haar_factor = 2^(haar_j/2);

    y = piecewise( haar_one_begin<x<haar_one_end,haar_factor,...
        haar_minus_one_begin<x<haar_minus_one_end,-1*haar_factor);

    y1 = piecewise( 0<x<0.25,1.4142,0.25<x<0.5,-1.4142);
    y2 = piecewise( 0.5<x<0.75,1.4142,0.75<x<1,-1.4142);
    y3 = piecewise( 1<x<1.25,1.4142,1.25<x<1.5,-1.4142);
    
%     fplot(y)
%     xlim([0 3])
%     xlabel('x axis');
%     ylabel('y axis');
%     
%     %figure
%     %fplot(y1)
%     %xlim([0 3])
%     
%     hold on;

end



y1 = piecewise( 0<x<0.25,1.4142,0.25<x<0.5,-1.4142);
y2 = piecewise( 0.5<x<0.75,1.4142,0.75<x<1,-1.4142);
y3 = piecewise( 1<x<1.25,1.4142,1.25<x<1.5,-1.4142);

figure
fplot(y1); hold on;
fplot(y2); hold on;
fplot(y3); hold on;
xlim([0 3])
xlabel('x axis');
ylabel('y axis');
