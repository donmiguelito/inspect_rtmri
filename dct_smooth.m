function fx_ = dct_smooth( fx,thr, vb )
%
%   FUNCTION:
%   smooth 1-d input function by eliminating least significant DCT
%   components
%
%   USAGE:
%   fx_ = dct_smooth( fx,thr, vb );
%
%   INPUTS:
%   fx [1 x n]:     input function (curve to smooth)
%   thr:            DCT smoothing threshold - specify either:
%                   [0..1]: percentage of contribution of DCT components
%                   (int):	explicit number of DCT components to use
%   vb (int):       plot verbosity: 0: work silently
%                                   1: report number of components, plot nothing
%                                   2: plot smoothed curve
%                                   3: plot iDCT, components & threshold
%
%   EXAMPLE:
%   fx_	= dct_smooth( fx,0.95, 1 );     reconstruct signal from DCT components contributing 95% of curvature
%   fx_	= dct_smooth( fx,5,    1 );     reconstruct signal from 5 most significant DCT components
%
%   AUTHOR:
%   M.Proctor (2010)
%
        
    ifx	= dct(fx);
    nfx = length(fx);
    
    [XX,ix] = sort(abs(ifx));
    ix = fliplr(ix);
    
    if (thr >= 1)
        thrsh = XX(end-thr-1);
        ifx( abs(ifx)<thrsh ) = 0;
        ndct = thr;
    else
        i = 1;
        while (norm([ifx(ix(1:i)) zeros(1,100-i)])/norm(ifx)<thr)
            i = i+1;
        end
        thrsh = XX(end-i-1);
        ifx( abs(ifx)<thrsh ) = 0;
        ndct = i;
    end
    fx_  = idct(ifx);
    
    if (vb),
        fprintf( '   Using %d most significant DCT components (of %d)\n\n', ndct,max(ix) );
        if (vb>1),
            figure;	hold on;
            plot( fx, 'b-'  );
            plot( fx_,'k--' );
            if (vb>2),
                figure;
                plot( ifx,'b' );
                figure;	hold on; 
                plot( XX, 'b' );
                plot( thrsh*ones(1,nfx),'r:' );
            end
        end
    end

end %main function
