clc;
clearvars;
close all;
%imtool close all;  % Close all imtool figures.
workspace;
format longg;
format compact;
fontSize = 16;
numberOfElements = 1001;
middleElement = ceil(numberOfElements/2);
x = linspace(-10,10, numberOfElements);
halfWindowWidth = 150;
% Make rectangular pulse.
rect = zeros(1, numberOfElements);
rect(middleElement-halfWindowWidth: middleElement+halfWindowWidth) = 1;
subplot(3, 1, 1);
plot(x, rect,'b-', 'LineWidth', 2);
ylim([0 2]);
grid 'on';
% Enlarge figure to full screen.
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
% Make sinc.
subplot(3, 1, 2);
s = sin(pi*x/2)./(pi*x/2);
plot(x, s, 'b-', 'LineWidth', 2);
grid 'on';
% Convolve
y = conv(s, rect,'full');
predictedLength = length(s)+length(rect)-1
newXAxis = linspace((-10 - halfWindowWidth), (+10 + halfWindowWidth), length(y)) ;
subplot(3, 1, 3);
plot(newXAxis, y, 'b-', 'LineWidth', 2);
grid 'on';