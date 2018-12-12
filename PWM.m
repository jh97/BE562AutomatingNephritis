function pwm = PWM(Im);
%{
    Get Progressive Weight Mean curve.
    Args:
        Im - Gray-scale image

    Returns:
        pwm - Data values on pwm curve
%}
% Get histogram values
hist_handle = histogram(Im, 0:256);
title('GrayScale Distr');
bin_edges = hist_handle.BinEdges;
bin_values = hist_handle.Values;
% PWM curve
pwm = [];
for i = 0:length(bin_values)-1
    top = 0;
    bottom = 0;
    for j = 0:i
        x = bin_edges(j+1);
        w = bin_values(j+1);
        top = top + x*w;
        bottom = bottom + w;
    end
    
    if (bottom == 0)
        newEl = 0;
    else
        newEl = top/bottom;
    end
    pwm = [pwm newEl];
end