clear


load test_6U.mat
y_avg = y_max';
%
newf = 0:1e-6:.5;

for i = 1:80
    y_avg(i) = 0 ;
end

skip = 1;
newf = newf(1:skip:20000);
y_avg = y_avg(1:skip:20000);

newf = [newf(1:1000),newf(1001:50:end)];
y_avg = [y_avg(1:1000),y_avg(1001:50:end)];

for i = 1:length(y_avg)
    y_avg(i) = round(y_avg(i),3,'significant');
end

% anything after 1000 can sample way way slower 

figure
hold on 
% plot(newf,y_avg,'k')
a = area(newf,y_avg);
a.FaceColor = 'k';
xlim([1e-4 1000])
% h=fill([1.65e-4,1.85e-4,1.85e-4,1.65e-4],[0,0,1,1],'red');
h.EdgeAlpha = 0.0;
%     h.FaceAlpha=0.3;
% plot([1.7e-4,1.7e-4],[0,1],'r')
h=fill([10,1000,1000,10],[0,0,1,1],'b');
    h.FaceAlpha=0.3;
    h.EdgeAlpha = 0.0;
legend('Maximum Normalized Magnitude','Reaction Wheel Control/Jitter')
% title('Monte Carlo Spectral Analysis of Disturbance Torques (1000 Trials)')
xlabel('Frequency (hz)')
ylabel('Maximum Normalized Magnitude')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off 
% saveas(gcf,'testplot.png')
% matlab2tikz('mc_frequencies_2.tex')


%% legend test 

clear

figure
hold on 
plot(1:10,randn(10,10))
legend('uea','location','best')
hold off 
matlab2tikz('legend_test.tex')

