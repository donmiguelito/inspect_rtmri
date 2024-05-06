function [xx_ fx_] = spline_smooth( xx,fx,res, vb )
%
%   FUNCTION:
%   smooth 1-d input function by fitting a spline
%
%   USAGE:
%   [fx_ xx_] = spline_smooth( xx,fx, vb );
%
%   INPUTS:
%   xx [1 x n]:     input range (function x-vector)
%   fx [1 x n]:     input function (curve to smooth)
%   res (int):      resolution: factor by which to oversample input function
%   vb (int):       plot verbosity: 0: plot nothing;
%                                   1: plot curves
%
%   OUTPUTS:
%   xx_ [1 x res*n]:	output range (interpolated x-vector)
%   fx_ [1 x res*n]:	output function (smoothed, oversampled curve)
%

    nfx = length(fx);
    xx_	= linspace(xx(1),xx(end),res*nfx);
    fx_	= spline(xx,fx,xx_);
    
    if (vb),
        figure;
        hold on;
        plot( xx, fx,  'bo-' );
        plot( xx_,fx_, 'r--' );
    end

end %main function
