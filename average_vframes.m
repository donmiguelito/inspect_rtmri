function immean = average_vframes( fn_source, fn_out, ff, imadj, vb, displ )
%
%  AVERAGE VIDEO FRAMES: Mean frames in rtMRI video.
%
%  immean = average_vframes( fn_source, fn_out, ff, vb, displ )
%
%    fn_source (string)	full filename (inclduing path) of source AVI file
%    fn_out (string)	full filename of file to write mean image
%    ff (1xn int)       indices of frames to use in mean image
%                       ff = []: create mean image from all frames in video
%    v (struct)         video object containing enhanced MRI data
%    imadj (flag)       image adjust:       0: don't adjust image
%                                           1: correct image intensity
%    vb (flag)        	activity verbosity:	0: work silently
%                                           1: report progress
%    displ (flag)       display verbosity:	0: don't display any image
%                                           1: display mean image
%
%    eg. im = average_vframes( 'lac07072013_11_54_01.avi', 'lac07072013_11_54_01_ave.jpg', 1:25, 0, 1, 1 );
%        im = average_vframes( 'lac07072013_11_54_01.avi', 'lac07072013_11_54_01_ave.jpg', [], 1, 1, 1 );
%
%  AUTHOR:	mike.i.proctor@gmail.com
%           Intensity correction algorithm: Adam Lammert (2010)
%  CREATED:	26-aug-2014
%

    % specify image, display & correction settings
    wbi         = 7;            % waitbar interval
    sensloc1	= [20 20];      % coil1 location for intensity correction
    sensloc2	= [50 20];      % coil2 location for intensity correction
    
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
    if (displ)
        new_tractplot;
        fmid	= ff(floor(nf/2));
        imagesc( vid(:,:,1,fmid) );
    end
    
    % reshape image sequence -> matrix of row vectors (1 row = 1 frame)
    if (vb)     % show progress bar if verbose mode
        hWB	= waitbar(0,'Correcting MRImage intensity: fetching image frames ...');
        jf	= get( get(hWB,'JavaFrame'), 'FigurePanelContainer');
        jf.getComponent(0).getRootPane.getTopLevelAncestor.setAlwaysOnTop(1);
    end
    vlen	= h*w;
    rowim	= zeros( nf,vlen );
    for f = ff
        if ((vb) && ~mod(f,wbi)), waitbar(f/nf); end;	% report progress
        rowim(f,:) = reshape(vid(:,:,1,f),1,vlen);
    end
    if (vb), close(hWB); end;

    % calculate mean image
    ImSize = sqrt(size(rowim,2));
    if ImSize ~= round(ImSize)
        error('Images must be square!');
    end
    immean	= uint8(reshape(mean(rowim,1),ImSize,ImSize));
    
    % correct image intensity, if flagged
    if (imadj)
        
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
        if (vb)     % show progress bar if verbose mode
            hWB	= waitbar(0,'Correcting individual frame intensities ...');
            jf	= get( get(hWB,'JavaFrame'), 'FigurePanelContainer');
            jf.getComponent(0).getRootPane.getTopLevelAncestor.setAlwaysOnTop(1);
        end
        Mask2 = reshape(Mask,1,ImSize^2);
        M2 = zeros(size(rowim));
        for f = 1:size(rowim,1)
            if ((vb) && ~mod(f,wbi)), waitbar(f/nf); end;	% report progress
            M2(f,:) = rowim(f,:)./Mask2;
        end
        if (vb), close(hWB); end;

        M2 = M2-min(min(M2));
        M2 = M2./max(max(M2));
        M2 = M2.*255;
        
        % recalculate mean image from intensity-corrected video data
        immean	= uint8(reshape(mean(M2,1),ImSize,ImSize));

    end

    % save and return (corrected) mean image
    immean	= imresize(immean, 5);
    imwrite( immean, fn_out, 'jpg' );
    if (displ)
        new_tractplot;
        imagesc(immean)
    end
    if (vb)
        fprintf('    Mean image written to <%s>\n\n', fn_out );
    end
        
end %of file

