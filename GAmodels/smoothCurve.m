%% funzione per ammorbidire la curva sperimetnale e trovare i pti di discontinuità.

function [interpolatedX,interpolatedY] = smoothCurve(t,I)

fontSize = 16;

% Find the coefficients.
n = 25; %numero coefficienti.

coeffs = polyfit(t, I, n);
%plot(t, I, 'ro', 'MarkerSize', n);

% Make a finer sampling so we can see what it
% does in between the training points.
interpolatedX = linspace(min(t), max(t), length(t));
interpolatedY = polyval(coeffs, interpolatedX);

end 

% % Plot the interpolated points.
% hold on;
% plot(interpolatedX, interpolatedY, 'b-', 'LineWidth', 3);
% grid on;
% title('Interpolating Polynomial', 'FontSize', fontSize);
% xlabel('t', 'FontSize', fontSize);
% ylabel('I', 'FontSize', fontSize);
% % Enlarge figure to full screen.
% set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
