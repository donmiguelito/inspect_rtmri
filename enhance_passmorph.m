function vid = enhance_passmorph( fn_source, fn_out, offset, nmean, vb, displ )
%
%  ENHANCE PASSIVE MORPHOLOGY: Adjust SNR in regions of MRI corresponding
%  to passive structures for improved analysis in inspect_drmi GUI.
%
%  v = enhance_passmorph( fn_source, fn_out, offset, nmean )
%
%    fn_source (string)	full filename (inclduing path) of source AVI file
%    fn_out (string)	full filename (inclduing path) of cropped AVI file
%    offset (int)       enhancement region [px] from top of image frame
%                       (if offset = [], load image for manual location)
%    nmean (int)        number of frames to use in mean image enhancement
%                       (if nf = [], create mean image from 1% of total no. of frames in video)
%    v (struct)         video object containing enhanced MRI data
%    vb (flag)        	activity verbosity:	0: work silently
%                                           1: report progress
%    displ (flag)       display verbosity:	0: don't display any image
%                                           1: display only composite image
%                                           2: display each image frame
%
%    eg. v = enhance_passmorph( 'lac07072013_11_54_01.avi', 'lac07072013_11_54_01_en.avi', 10, 10, 1, 1 );
%        v = enhance_passmorph( 'lac07072013_11_54_01.avi', 'lac07072013_11_54_01_en.avi', [], [], 1, 1 );
%
%  AUTHOR:	mike.i.proctor@gmail.com
%  CREATED:	21-aug-2014
%  
%  

    % specify image contrast & brightness levels & display settings
    vloc	= 250;              % vertical location of window for manual cropping
    hloc	= 450;              % horizonatl location of window for manual cropping
    imscale	= 10;               % image scale factor for manual crop frame placement
    cfactor	= 0.74;             % contrast level of output JPG images
    bfactor	= 0.07;             % background brightness factor of output JPG images
    imcont	= cfactor * 255;
    imbrite	= bfactor * 255;
    
    % fetch video & metadata
    v	= VideoReader( fn_source );
    vid	= read(v);
    h	= v.Height;
    w	= v.Width;
    nf	= v.NumberOfFrames;
    len	= v.Duration;
    fps	= v.FrameRate;

    % report frame and sample correspondances
    if (vb)
        fprintf('\n    %d frames found in <%s>\n', nf,fn_source );
        fprintf(  '    %0.00f sec of video @ %d f.p.s\n', len,fps );
        fprintf(  '    Original frame size: %d x %d px\n\n', w,h );
    end

    % show frame from original video to locate crop frame
    if isempty(offset)
        hf = figure; 
        set( hf, 'position', [hloc vloc w*imscale h*imscale] );
        imagesc( vid(:,:,1,1) );
        set( gcf,'Colormap',gray );
        [x,y]	= ginput(1);
        offset	= floor(y);
        close(hf)
    end

    % calculate mean frame over subset of frames from whole video
    if isempty(nmean)
        nmean = floor(nf/100);
    end
    if (nmean<1)
        nmean = 1;
    end
    f_mean	= vid(:,:,1,1);
    for f = 2:nmean
        f_mean	= imfuse(f_mean,vid(:,:,1,f));
        %f_mean	= uint8( uint8(f_mean) + uint8(vid(:,:,1,f)) );
    end
    %f_mean	= uint8(f_mean/nmean);
    f_mean	= (f_mean(:,:,1) + f_mean(:,:,3) + f_mean(:,:,3))/3;
    f_mean	= (f_mean(:,:,1));
    if (displ)
        new_tractplot;
        imagesc(f_mean);
        colormap(gray);
    end

    % create output video: top region replaced with mean frame data
    mask	= repmat(f_mean, [1 1 1 nf]);
    vid( 1:offset-1, :, :, : ) = mask( 1:offset-1, :, :, : );

    % display middle frame from output video
    if (displ)
        h	= new_tractplot; 
        imagesc( vid(:,:,1,floor(nf/2)) );
        set( h,'Colormap',gray );
    end

    % write cropped video to output file
    v_out = VideoWriter(fn_out);
    v_out.FrameRate = fps;
    v_out.Quality	= 100;
    open(v_out);
    writeVideo(v_out,vid);
    close(v_out);

    % report frame and sample correspondances
    if (vb)
        fprintf('    Output video written to <%s>\n', fn_out );
        fprintf('    Top %d rows replaced with mean image data\n', offset );
        fprintf('    Mean image constructed from %d frames\n\n', nmean );
    end
        
end %of file

