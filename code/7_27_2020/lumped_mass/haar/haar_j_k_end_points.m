clc;
clear;

j_index = [0,1,2,3,4,5,6];
k_index = [0,1,2,3,4];



haar_j_k_plus_end_points_matrix = zeros(length(j_index),length(k_index));
haar_j_k_minus_end_points_matrix = zeros(length(j_index),length(k_index));


for row=1:length(j_index)
    for column = 1:length(k_index)
        
        haar_j = j_index(row);
        haar_k = k_index(column);
        
        plus_end = (haar_k+1/2)/(2^haar_j);
        minus_end = (haar_k+1)/(2^haar_j);
        
        haar_j_k_plus_end_points_matrix(row,column) = plus_end;
        haar_j_k_minus_end_points_matrix(row,column) = minus_end;
        
    end
end


haar_j_k_plus_end_points_matrix
haar_j_k_minus_end_points_matrix

jk_pair = [0,0;1,1;,2,2;1,0;2,0;2,1]

for row=1:6
    haar_j = jk_pair(row,1)
    haar_k = jk_pair(row,2)
    
    plus_begin = haar_k/(2^haar_j)
    plus_end = (haar_k+1/2)/(2^haar_j)
    minus_begin = (haar_k+1/2)/(2^haar_j)
    minus_end = (haar_k+1)/(2^haar_j)
    
end
