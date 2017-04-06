function h_rose = rose_with_stats(theta, n_bins, with_stats)

if nargin < 3
    with_stats = false;
end
[t_out,r_out] = rose(theta,n_bins);
r_max = max(r_out)*1.1;

if with_stats
    % Make an invisible line to push polar axis limit out to allow mean and
    % stddev lines inside plot
    t_temp = 0:0.01:2*pi;
    P_temp = polar(t_temp, 1.1*r_max * ones(size(t_temp)));
    set(P_temp, 'Visible', 'off')
else
    h_rose = polar(t_out,r_out,'b');
end

if with_stats
    hold on;
    addpath('CircStat2012a/')
    h_rose = polar(t_out,r_out,'b');
    theta_mean = circ_mean(theta);
    theta_std = circ_std(theta);
    p = circ_otest(theta);
    line_colour = 'k';
    polar([theta_mean,theta_mean],[0,r_max],line_colour);
    n_points = round(100*theta_std);
    arc_min = theta_mean-theta_std;
    arc_max = theta_mean+theta_std;
    arc_range = linspace(arc_min,arc_max,n_points);
    r_range = r_max*ones(1,n_points);
    polar(arc_range,r_range,line_colour);
    
    % end lines
    polar([arc_min,arc_min],[0.95*r_max,1.05*r_max],line_colour);
    polar([arc_max,arc_max],[0.95*r_max,1.05*r_max],line_colour);
    
    p_string = sprintf('Omnibus test\np = %.2g',p);
    text(r_max*1.2,r_max*1.2,p_string);
    
    hold off
end