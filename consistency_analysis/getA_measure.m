function [A_measure, scaled_A_measure] = getA_measure(x0,x1)
    [p,h,stats]=ranksum(x0,x1);
    A_measure=(stats.ranksum/length(x0)-(length(x0)+1)/2)/length(x1);
    scaled_A_measure=0.5+abs(0.5-A_measure);
end