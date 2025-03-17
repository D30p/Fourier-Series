% Oliviery Amadeo Eide Rusli
% 13123103
% Code to solve Odd and Even Function

clc; clear; close all;

% Define parameters
a = 0;
b = 3;
T = 3;
N_values = [1, 2, 4, 8, 16, 32, 64, 128];
N = 128;
ecos = zeros(1, N);
esin = zeros(1, N);


% Define the piecewise function
f_piecewise = @(x) (x.^3 - 5*x.^2 + 5*x + 1) .* ((x > 0) & (x < 3));

% Define x values
x = linspace(-b, b, 1000); % Extend x to make function even and odd

fid = fopen('Fourier_Coefficients.txt', 'w');

for N = N_values
    % Compute Fourier Series approximation
    a0 = (1 / T) * trapz(x, f_piecewise(x)); % Ensure even symmetry for cosine terms
    sum_series = a0;
    sum_sins = 0;
    
    fprintf(fid, 'For N = %d:\n', N);
    fprintf(fid, 'a0 = %.8f\n', a0);
    
    for n = 1:N
        an = (2 / T) * trapz(x, f_piecewise(x) .* cos(pi * n * x / T)); % Even function for cosine
        bn = (2 / T) * trapz(x, f_piecewise(x) .* sin(pi * n * x / T)); % Keep sine as is for odd function
        fprintf(fid, 'an %d = %.8f    bn %d = %.8f\n', n, an, n, bn);
        
        sum_series = sum_series + an * cos(pi * n * x / T);
        sum_sins = sum_sins + bn * sin(pi * n * x / T);
        y_actual = f_piecewise(x);
        errorcos = abs(sum_series - y_actual);
        errorsin = abs(sum_sins - y_actual);
        error_integral_cos = trapz(x, errorcos);
        error_integral_sin = trapz(x, errorsin);
        ecos(n) = error_integral_cos;
        esin(n) = error_integral_sin;
    end
    
    % Create new figure for each N
    fig = figure('Visible', 'off', 'Position', [100, 100, 1200, 800]);
    subplot(2,1,1);
    plot(x, f_piecewise(x), 'k--', 'LineWidth', 2); hold on;
    plot(x(x>0 & x<3), sum_series(x>0 & x<3), 'r', 'LineWidth', 2);
    plot(x(x>0 & x<3), sum_sins(x>0 & x<3), 'b', 'LineWidth', 2);
    plot(x(x>0 & x<3) - 2.*T, sum_series(x>0 & x<3), 'LineWidth', 2, 'Color', [0.4 0 0 0.4]);
    plot(x(x>0 & x<3) + 2.*T, sum_series(x>0 & x<3), 'LineWidth', 2, 'Color', [0.4 0 0 0.4]);
    plot(x(x>0 & x<3) - 2.*T, sum_sins(x>0 & x<3), 'LineWidth', 2, 'Color', [0 0 0.4 0.4]);
    plot(x(x>0 & x<3) + 2.*T, sum_sins(x>0 & x<3), 'LineWidth', 2, 'Color', [0 0 0.4 0.4]);
    plot(x(x>-3 & x<0), sum_series(x>-3 & x<0), 'LineWidth', 2, 'Color', [0.4 0 0 0.4]);
    plot(x(x>-3 & x<0), sum_sins(x>-3 & x<0), 'LineWidth', 2, 'Color', [0 0 0.4 0.4]);
    plot(x(x>-3 & x<0)+2.*T, sum_series(x>-3 & x<0), 'LineWidth', 2, 'Color', [0.4 0 0 0.4]);
    plot(x(x>-3 & x<0)+2.*T, sum_sins(x>-3 & x<0), 'LineWidth', 2, 'Color', [0 0 0.4 0.4]);
    plot(x(x>-3 & x<0)-2.*T, sum_series(x>-3 & x<0), 'LineWidth', 2, 'Color', [0.4 0 0 0.4]);
    plot(x(x>-3 & x<0)-2.*T, sum_sins(x>-3 & x<0), 'LineWidth', 2, 'Color', [0 0 0.4 0.4]);
    
    
    
    % Labels and title
    title(sprintf('Fourier Series Approximation (N = %d)', N));
    xlabel('x');
    ylabel('y');
    grid on;
    legend('Piecewise Function (Real)', 'Fourier Cosine Approximation (Even)', 'Fourier Sine Approximation (Odd)');
    
    % Compute error
    y_actual = f_piecewise(x);
    errorcos = abs(sum_series - y_actual);
    errorsin = abs(sum_sins - f_piecewise(x));
    
    % Plot the error
    subplot(2,1,2);
    plot(x(x>0 & x<3), errorcos(x>0 & x<3), 'b', 'LineWidth', 2); hold on;
    plot(x(x>0 & x<3), errorsin(x>0 & x<3), 'r', 'LineWidth', 2);
    
    title(sprintf('Error Between Fourier Approximation and Actual Function (N = %d)', N));
    xlabel('x');
    ylabel('Error');
    legend('Error for cosine series', 'Error for sine series');
    grid on;
    
    % Compute total error metrics
    total_error_cos = sum(errorcos) / length(errorcos);
    total_error_sin = sum(errorsin) / length(errorsin);
    error_integral_cos = trapz(x, errorcos);
    error_integral_sin = trapz(x, errorsin);
    ecos(n) = error_integral_cos;
    esin(n) = error_integral_sin;
    fprintf(fid, 'Error cos (integral) = %.8f\n', error_integral_cos);
    fprintf(fid, 'Max error for one x value (cos)= %.8f\n', max(errorcos));
    fprintf(fid, 'Error sin (integral) = %.8f\n', error_integral_sin);
    fprintf(fid, 'Max error for one x value (sin)= %.8f\n\n', max(errorsin));

    %Plot total integral error to n
    
    % Save figure to PDF
    set(fig, 'PaperPositionMode', 'auto');
    set(fig, 'PaperOrientation', 'landscape');
    set(fig, 'PaperUnits', 'inches', 'PaperSize', [35, 16]);
    print(fig, sprintf('Fourier_Series_N%d.pdf', N), '-dpdf', '-bestfit');
    close(fig);
end

fclose(fid);

% Convert coefficient text file to PDF
system('pandoc Fourier_Coefficients.txt -o Fourier_Coefficients.pdf');

% Graph Error to n
fig_even = figure('Visible', 'off');
plot(1:N, ecos, 'r-');
title('Error Even Function Terhadap Jumlah Suku (n)');
xlabel('n');
ylabel('error');
grid on;
saveas(fig_even, 'Error_Even_Function.png');
close(fig_even);

fig_odd = figure('Visible', 'off');
plot(1:N, esin, 'b-');
title('Error Odd Function Terhadap Jumlah Suku (n)');
xlabel('n');
ylabel('error');
grid on;
saveas(fig_odd, 'Error_Odd_Function.png');
close(fig_odd);