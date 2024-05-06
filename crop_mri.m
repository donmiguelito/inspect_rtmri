function vid_ = crop_mri( fn_source, fn_out, offset )
%
%  CROP MRI: Prepare MRI for analysis in inspect_drmi GUI.
%
%  v = crop_mri( fn_source, fn_out, offset )
%
%    fn_source (string)	full filename (inclduing path) of source AVI file
%    fn_out (string)	full filename (inclduing path) of cropped AVI file
%    offset (int [1x2])	location of top left corner of crop frame (px)
%                       if offset = [], load image for manual crop location
%    v (struct)         video object containing cropped MRI data
%
%    eg. v = crop_mri( 'SW_realtime_vowels_trial1.avi', 'trial1_68px.avi', [20 30] );
%        v = crop_mri( 'SW_realtime_vowels_trial1.avi', 'trial1_68px.avi', [] );
%
%  AUTHOR:	mike.i.proctor@gmail.com
%  CREATED:	17-aug-2014
%  
%  

    % declare constants
    vloc	= 150;                  % vertical location of window for manual cropping
    hloc	= 550;                  % horizonatl location of window for manual cropping
    imscale	= 6;                    % image scale factor for manual crop frame placement
    im_w	= 68;                   % frame width of target video (px)
    im_h	= 68;                   % frame height of target video (px)
    crop_ok = '';
    
    % parse inputs
    if (nargin ~= 3)
        disp('   Specify filenames (with path) of source and output AVI files + crop window location');
        disp('   ');
        vid_ = [];
        return;
    end

    % fetch video & metadata
    v	= VideoReader( fn_source );
    vid	= read(v);
    h	= v.Height;
    w	= v.Width;
    nf	= v.NumberOfFrames;
    len	= v.Duration;
    fps	= v.FrameRate;

    % report frame and sample correspondances
	fprintf('\n    %d frames found in <%s>\n', nf,fn_source );
	fprintf(  '    %0.00f sec of video @ %d f.p.s\n', len,fps );
	fprintf(  '    Original frame size: %d x %d px\n\n', w,h );

    % show frame from original video to locate crop frame
    if isempty(offset)
        hf = figure; 
        set( hf, 'position', [hloc vloc w*imscale h*imscale] );
        while ~(strcmp( crop_ok, 'y' ))
            imagesc( vid(:,:,1,1) );
            set( gcf,'Colormap',gray );
            [x,y] = ginput(1);
            x	= floor(x);
            y	= floor(y);
            hl	= line( [x x x+im_w x+im_w x], [y y+im_h y+im_h y y] );
            set( hl, 'Color','y', 'LineWidth',2, 'LineStyle','--' );
            str	= sprintf( '    Crop location (x = %d, y = %d) OK (y/n)?', x,y );
            crop_ok	= input( str,'s' );
        end
        offset = [x y];
        close(hf)
    end

    % crop video frames
	xoff	= offset(2);
	yoff	= offset(1);
    vid_	= vid( xoff:(xoff+im_w-1), yoff:(yoff+im_h-1), :, : );
    hf      = figure; 
    imagesc( vid_(:,:,1,1) );
    set( gcf,'Colormap',gray );
    axis equal; axis tight;

    % write cropped video to output file
    v_out = VideoWriter(fn_out);
    v_out.FrameRate = fps;
    v_out.Quality	= 100;
    open(v_out);
    writeVideo(v_out,vid_);
    close(v_out);


end %of file

