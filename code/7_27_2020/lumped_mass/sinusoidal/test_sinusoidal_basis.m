clc;
clear;

basis_rank = 5;
record_integral = zeros(basis_rank,basis_rank);

for row = 1:basis_rank
    for column = 1:basis_rank
        fun = @(x) sqrt(2)*sin (2*pi*x*row).*sqrt(2).*sin(2*pi*x*column);
        integral_result = integral(fun,0,1);
        
        record_integral(row,column) = integral_result;
        
        
    end
end

record_integral


x  = linspace(0, 1, 100).';
plot_basis = [];
for index = 1:basis_rank
    U0 = sin(2*pi*x*index);
    plot_basis = [plot_basis,U0];
end

plot(plot_basis)