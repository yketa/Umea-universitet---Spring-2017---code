% This script is designed to display the results of the 'simulation'
% executable.

T = readtable('results.csv'); % results table
len = size(T); % length of table T

% nEllipsoids = (len(2) - 8)/18; % number of ellipsoids
nEllipsoids = 4;
N = 100; % 1 every N points will be displayed

figure;
hold on;
for i = 0:nEllipsoids -1
    eval(strcat('scatter3(T.x',num2str(i),'(1:N:end),T.y',num2str(i),'(1:N:end),T.y',num2str(i),'(1:N:end),''DisplayName'',''particle',num2str(i),''');')); % 3D trajectory of ellipsoid i
end
legend('-DynamicLegend'); % displaying legend

savefig('results.fig'); % saving figure
