% Oliviery Amadeo Eide Rusli
% 13123103
% General Version of Fourier Series

clc; clear; close all;

% Define parameters
a = -2;
b = 2;
T = 4;
N_values = [1, 2, 4, 8, 16, 32, 64, 128, 194];

% Define the piecewise function
f_piecewise = @(x) (-x/2).*((x >= -2) & (x < 0)) + (2*x - (x.^2)/2).*((x >= 0) & (x <= 2));

% Define x values
x = linspace(a, b, 1000);
x_func = linspace(a, b, 1000);

pdf_filename = 'Fourier_Series_Results.pdf';
figures = [];
fid = fopen('Fourier_Coefficients.txt', 'w');

for N = N_values
    % Compute Fourier Series approximation
    a0 = (2 / T) * trapz(x, f_piecewise(x));
    sum_series = a0 / 2;
    
    fprintf(fid, 'For N = %d:\n', N);
    fprintf(fid, 'a0 = %.8f\n', a0);
    
    for n = 1:N
        an = (2 / T) * trapz(x, f_piecewise(x) .* cos(2 * pi * n * x / T));
        bn = (2 / T) * trapz(x, f_piecewise(x) .* sin(2 * pi * n * x / T));
        fprintf(fid, 'an %d = %.8f    bn %d = %.8f\n', n, an, n, bn);
        sum_series = sum_series + an * cos(2 * pi * n * x / T) + bn * sin(2 * pi * n * x / T);
    end
    
    % Create new figure for each N
    fig = figure;
    figures = [figures, fig];
    clf;
    subplot(2,1,1);
    plot(x_func, f_piecewise(x_func), 'k--', 'LineWidth', 2); hold on;
    plot(x_func + T, f_piecewise(x_func), 'k--', 'LineWidth', 2);
    plot(x_func - T, f_piecewise(x_func), 'k--', 'LineWidth', 2);
    
    % Plot the Fourier series approximation
    y_fourier = sum_series;
    plot(x, y_fourier, 'r', 'LineWidth', 2);
    plot(x + T, y_fourier, 'r', 'LineWidth', 2);
    plot(x - T, y_fourier, 'r', 'LineWidth', 2);
    
    % Labels and title
    title(sprintf('Fourier Series Approximation (N = %d)', N));
    xlabel('x');
    ylabel('y');
    grid on;
    legend('Piecewise Function', 'Fourier Approximation');
    
    % Compute error
    y_actual = f_piecewise(x);
    error = abs(y_fourier - y_actual);
    
    % Plot the error
    subplot(2,1,2);
    plot(x, error, 'b', 'LineWidth', 2); hold on;
    plot(x - T, error, 'b', 'LineWidth', 2);
    plot(x + T, error, 'b', 'LineWidth', 2);
    
    title(sprintf('Error Between Fourier Approximation and Actual Function (N = %d)', N));
    xlabel('x');
    ylabel('Error');
    grid on;
    
    % Compute total error metrics
    total_error = sum(error) / length(error);
    error_integral = trapz(x, error);
    fprintf(fid, 'Error (integral) = %.8f\n', error_integral);
    fprintf(fid, 'Max error for one x value = %.8f\n\n', max(error));
end

fclose(fid);

% Save all figures to separate PDF files
for i = 1:length(figures)
    saveas(figures(i), sprintf('Fourier_Series_N%d.pdf', N_values(i)), 'pdf');
end

% Convert coefficient text file to PDF
system('pandoc Fourier_Coefficients.txt -o Fourier_Coefficients.pdf');