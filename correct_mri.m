function [vid, ROI] = correct_mri( fn_source, fn_out, ff, ROI, vb, displ )
%
%  CORRECT MRI: correct for intensity bias in rtMRI video.
%
%  [vid, roi] = correct_mri( fn_source, fn_out, ff, ROI, vb, displ )
%
%  INPUTS
%    fn_source (string)	full filename (inclduing path) of source AVI file
%    fn_out (string)	full filename of AVI file to write adjusted video
%    ff (1xn int)       indices of frames to use in mean image
%                       ff = []: create mean image from all frames in video
%    ROI (2xn int)      set of pts defining region in which to apply mean image mask
%                       0:  don't mask any part of image
%                       []: manually select mask region limits with mouse
%    vb (flag)        	activity verbosity:	0: work silently
%                                           1: report progress
%    displ (flag)       display verbosity:	0: don't display any image
%                                           1: display mean image
%
%  OUTPUTS
%    vid (4-D double)	corrected video frames
%    ROI (2xn int)      set of pts defining region in which image mask was applied
%
%  EXAMPLE
%    [vid, roi] = correct_mri( 'avi_/lac07072013_11_54_01.avi', 'avi/lac07072013_11_54_01.avi', [], [],  1, 1 );
%    [vid, roi] = correct_mri( 'avi_/lac07072013_11_54_01.avi', 'avi/lac07072013_11_54_01.avi', [], roi, 1, 1 );
%
%  AUTHOR:	mike.i.proctor@gmail.com
%           Intensity correction algorithm: Adam Lammert (2010)
%  CREATED:	26-aug-2014
%

    % specify image, display & correction settings
    wbi         = 7;            % waitbar interval
    bfm         = 1.7;          % mask brightening factor
    bfi         = 1.0;          % image brightening factor
    sensloc1	= [20 20];      % coil1 location for intensity correction
    sensloc2	= [50 20];      % coil2 location for intensity correction
    mask        = [];           % initialize ROI mask
    mask_fg     = [];           % initialize image mask
    
    % fetch video & metadata; calculate range of frames to analyze
    v	= VideoReader( fn_source );
    vid	= read(v);
    h	= v.Height;
    w	= v.Width;
    len	= v.Duration;
    fps	= v.FrameRate;
    if isempty(ff)
        ff	= 1:v.NumberOfFrames;
    end
    nf	= length(ff);
    
    % report video metadata (and display middle frame from original video)
    if (vb)
        fprintf('\n    %d frames (%d x %d px) found in <%s>\n', v.NumberOfFrames,w,h,fn_source );
        fprintf(  '    %0.00f sec of video @ %0.00f f.p.s\n', len,fps );
        fprintf(  '    Analyzing %d frames: %d to %d\n', nf,ff(1),ff(end) );
    end
    
    % reshape image sequence -> matrix of row vectors (1 row = 1 frame)
    if (vb), hWB = waitbar(0,'Correcting MRImage intensity: fetching image frames ...'); end;
    vlen	= h*w;
    rowim	= zeros( nf,vlen );
    for f = ff
        if ((vb) && ~mod(f,wbi)), waitbar(f/nf); figure(hWB); end;	% report progress
        rowim(f,:) = reshape(vid(:,:,1,f),1,vlen);
    end
    if (vb), close(hWB); end;

    % calculate mean image
    ImSize = sqrt(size(rowim,2));
    if ImSize ~= round(ImSize)
        error('Images must be square!');
    end
    immean	= uint8(reshape(mean(rowim,1),ImSize,ImSize));
    
    % Get the distance maps (D1 and D2) from each coil
    % and calculate combined, quantized distance, D3
    D1 = zeros(ImSize,ImSize);
    D2 = zeros(ImSize,ImSize);
    for i = 1:ImSize
        for j = 1:ImSize
            D1(i,j) = sqrt(sum((sensloc1-[i j]).^2));
            D2(i,j) = sqrt(sum((sensloc2-[i j]).^2));
        end
    end
    D3 = round(D1+D2);

    % Get average intensity at each quantized distance
    X = min(min(D3)):max(max(D3));
    I = zeros(1,length(X));
    for i = 1:length(X)
        I(i) = mean(immean(D3==X(i)));
    end

    % Assume intensities monotonically increasing
    I2 = zeros(1,length(I));
    I2(end) = I(end);
    for i = length(I)-1:-1:1
        if I(i) > I2(i+1)
            I2(i) = I(i);
        else
            I2(i) = I2(i+1);
        end
    end

    % build correction map (Mask) & correct video sequence
    Mask = zeros(ImSize,ImSize);
    for i = 1:length(X)
        Mask(D3==X(i)) = I2(i);
    end
    if (vb), hWB = waitbar(0,'Correcting individual frame intensities ...'); end;
    Mask2 = reshape(Mask,1,ImSize^2);
    M2 = zeros(size(rowim));
    for f = 1:size(rowim,1)
        if ((vb) && ~mod(f,wbi)), waitbar(f/nf); figure(hWB); end;	% report progress
        M2(f,:) = rowim(f,:)./Mask2;
    end
    if (vb), close(hWB); end;

    M2 = M2-min(min(M2));
    M2 = M2./max(max(M2));
    M2 = M2.*255;

    % recalculate (and display) mean image from intensity-corrected video data
    immean	= uint8(reshape(mean(M2,1),ImSize,ImSize));
    if (displ)
        new_tractplot;
        imagesc(immean)
    end

    % show mean frame with corrected intensity to locate region limits for mask
    if isempty(ROI)
        new_tractplot;	imagesc( immean );
        [BW,roiX,roiY]	= roipoly;
        ROI	= floor([roiX roiY]);
    elseif (length(ROI) > 1)
        BW	= poly2mask( ROI(:,1),ROI(:,2), h,w );
    end
    if (length(ROI) > 1)
        mask_fg	= uint8(~BW);
        mask	= immean.*uint8(BW);
        mask	= uint8(bfm .* mask);
    end
    
    % recreate movie from intensity-corrected image row matrix 
    if (vb), hWB = waitbar(0,'Reconstructing intensity-corrected AVI frames ...'); end;
    vid	= zeros( h,w,1,nf);
    for f = ff
        if ((vb) && ~mod(f,wbi)), waitbar(f/nf); figure(hWB); end;	% report progress
        im = uint8(reshape(M2(f,:),h,w));
        if (length(ROI) > 1)
            im = im .* mask_fg;
            im = uint8(bfi .* im);
            im = im + mask;
        end
        vid(:,:,1,f) = im;
    end
    if (vb), close(hWB); end;

    % write corrected video to output file
    v_out = VideoWriter(fn_out);
    v_out.FrameRate = fps;
    v_out.Quality	= 100;
    open(v_out);
    writeVideo(v_out,uint8(vid));
    close(v_out);

    % report frame and sample correspondances
    if (vb)
        fprintf('    Output video written to <%s>\n\n', fn_out );
    end    
        
end %of file

