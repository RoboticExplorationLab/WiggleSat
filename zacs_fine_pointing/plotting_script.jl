using JLD2, MATLAB


@load "mc_pointing.jld2" W_telescope_range mean_line lower_bound upper_bound


# /Users/kevintracy/devel/Low-Thrust-TrajOpt/matlab2tikz/src

mat"
figure
hold on
plot($W_telescope_range,$mean_line,'b','linewidth',2)
plot($W_telescope_range,$lower_bound,'r','linewidth',2)
plot($W_telescope_range,$upper_bound,'r','linewidth',2)
xlabel('Sensing Error (arcseconds)')
ylabel('RMS Body Pointing Error (arcseconds)')
xlim([1e-4,10])
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
legend('Mean','2 sigma Bound','Location','SouthEast')
hold off
addpath('/Users/kevintracy/devel/Low-Thrust-TrajOpt/matlab2tikz/src')
matlab2tikz('ponting_mc.tikz')
%saveas(gcf,'rms_bounds_pointing.png')
"
