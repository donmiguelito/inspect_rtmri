function dat = annotate_dmri( x,grid )
%
%  INSPECT DMRI: GUI for annotating vocal tract contours on MRI images.
%  (see Narayanan et al. 2004 for details of acquisition and protocols).
%
%  Assumes existence of synchronised audio and video files which share the 
%  base filename 'token', in dedicated subdirectories of working directory.
%
%  dat = annotate_dmri( mri_video )
%
%

% declare constants
mypi	= 3.141593;
deg2rad	= mypi/180;
rad2deg	= 1/deg2rad;
gui_wd	= 700;
gui_ht	= 900;

% global configuration flags
invert_video	= 0;        % set flag to invert video on some OS/codec combinations

% Interpanel spacing
hoff = 0.008;  voff = 0.008; % inter-panel spacing

% Video panel height
hvid = 600/gui_ht;                 % panel height
wvid = 4/5;

FOV     = 100;                                 % MRI field of view (%)
wd_im	= 200;                                 % width (mm) of MR Image at 100% FOV
wd_vid	= wvid*gui_wd;          % width (px) of video display window
wbi     = 5;                                   % waitbar interval
cm256	= gray(256);        % default colormap
cm_sg	= cm256.^5;         % spectrogram color map
cm_corr	= jet(64);          % correlation analysis color map
dcr     = 5;                % size of cross at displayed points
dpx     = 15;               % size of cross indicating correlation pixel location
colwin	= [1   1   0.7];	% bg color of zoomed inspection window
colwlim	= [1   0.8 0.4];	% color of zoomed window borders
colgrn	= [0   0.7 0  ];	% color of green pushbutton controls
colcyc	= {'red','blue','green','black','yellow'};

% declare initial parameters
sig_wwd	= 1;                % default zoomed audio window width (seconds)
min_wwd	= 0.01;             % min width of zoomed audio window (seconds)
sg_pree	= 1;                % spectrogram signal pre-emphasis enable
sg_win	= 20;               % spectrogram window width (msec)
sg_step	= 5;                % spectrogram inter-window step length (msec)
sg_Flim	= 5000;             % spectrogram frequency upper plot limit
sg_nfft	= 512;              % spectrogram number of FFT bins
sg_nf	= 4;                % spectrogram number of formants to track
gwidth	= 60;               % initial width of Ohman analysis grid (px)
gint	= 8;                % intervals between parallel gridlines (px)
CoilL1	= [23 1];           % assumed location of 1st MRI sensor coil
CoilL2	= [45 1];           % assumed location of 1st MRI sensor coil

% supress warning associated with java code managing waitbar
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');

% check input
if isstruct(x)
    old_dat = x;
    tok = x.tok;
elseif ischar(x)
    old_dat = [];
    tok = x;
else
    fprintf('\n   Input argument must be a string or structure ...\n');
    dat = [];
    return;
end
    
% create main GUI + global data structure
disp(' ');
hGUI	= figure( 'Name',['Review dMRI frames: ' tok] ); axis off;
set( hGUI, 'menubar','none', 'toolbar','none', 'Resize','off');
dat     = guihandles(hGUI);
guidata(hGUI,dat);
dat.hGUI_id = hGUI;
dat.tok	= tok;

% configure GUI color scheme
bgc	= get( hGUI, 'Color' );
set( hGUI, 'DefaultUIControlBackgroundColor',bgc );	% default control color to that of background
bgd = get(0, 'FactoryUIControlBackgroundColor');    % interactive controls default to darker bgc

% construct filenames/paths
home	= pwd();        % base directory containing audio and video subdirs
dir_vid	= 'avi';        % directory containing video
dir_ann = 'annotations';% directory with annotated images 
fn_vid	= [ tok '.avi' ];
dat.ffn_vid	= fullfile( home, dir_vid, fn_vid );

% inialize global configuration flags
dat.status.reload	= 0;
dat.status.FOV      = FOV;
dat.status.saved	= 0;
dat.status.exit     = 0;
dat.status.radanal	= 5;
if isstruct(old_dat)
    dat.status.reload	= 1;
end

% initialize data structures
dat.aud	= [];           % audio signal properties and display parameters
dat.sg	= [];           % spectral properties and display parameters
dat.vid	= [];           % MRI video signal & properties
dat.tvd	= [];           % tracked MRI video
dat.cvd	= [];           % intensity-corrected MRI video
dat.lvd	= [];           % low pass-filtered MRI video
dat.vt	= [];           % vocal tract delimitation data
dat.seg	= [];           % segmental analysis data
dat.time	= 0;
dat.fnum	= 1;
dat.fstep =1;

dat.vid.info	= [];   dat.tvd.info	= [];
dat.vid.avi     = [];   dat.tvd.avi     = [];
dat.cvd.avi     = [];   dat.lvd.avi     = [];
dat.vid.avi_	= [];
dat.vid.fint	= [0 0];
dat.vid.scale	= 1;
dat.vid.dims	= 0;
%dat.vid.invert	= invert_video;
dat.wwd = sig_wwd;
dat.annotation.show = false;
dat.annotation.labels = {'Epiglottis','Tongue','Subling. Cavity','Lower Lip',...
    'Chin','Lower Chin','Pharyngeal Wall','Upper back head','Back head','Low back head',...
    'Hard Palate','Velum', 'Label 2','Nose','Upper Lip' };
dat.annotated_part_id = 1;
dat.annotation.point_handles = [];
dat.edit_annotations = false;

dat.casy.loaded = false;

% create main GUI
ssz	= get(0,'ScreenSize');
gui_lf	= (ssz(3)-gui_wd)/2;
gui_bt  = (ssz(4)-gui_ht)/2;
gui_pos	= [gui_lf gui_bt gui_wd gui_ht];
set( hGUI, 'Position', gui_pos, 'toolbar','none' ); axis off;
set( hGUI, 'Resize','off');

% GUI font sizes
fsx = 11;  fsh = 10;  fsl = 9;  fss = 8;  fst = 7;  fsm = 6;

% calculate dimension of GUI panels
hana = 0.33;
haud = hvid-hana-voff;
hclo = 0.275*hana;
hexp = hana-hclo-voff;
hlab = 0.275*hana;
hcor = hana-hlab-voff;
baud = hana+2*voff;             % locate bottom edge of panels
bexp = hclo+2*voff;
bcor = hlab+2*voff;
wcor = 0.08;
wexp = wcor;
laud = wvid+2*hoff;             % locate left edge of panels
waud = 1-wvid-3*hoff;
wana = waud-wexp-wcor-2*hoff;
lcor = laud+wana+hoff;
lexp = lcor+wcor+hoff;

% create GUI panels
pVID = uipanel(	'Parent',hGUI, 'BackgroundColor',bgc, 'Title','Video',	...
                'FontSize',fsx, 'Position',[hoff 1-voff-hvid wvid hvid]);
pANN = uipanel( 'Parent',hGUI, 'BackgroundColor',bgc, 'Title','Annotation',...
                'FontSize',fsx, 'Position',[2*hoff+wvid 1-voff-hvid 1-2*hoff-wvid hvid]);
pCASY = uipanel( 'Parent',hGUI, 'BackgroundColor',bgc, 'Title','CASY',...
                'FontSize',fsx, 'Position',[hoff voff 1-2*hoff 1-2*voff-hvid]);
            
% calculate gui panel dimensions in pixels 
hpvid = hvid * gui_ht;	% panel heights
wpvid = wvid * gui_wd;	% panel widths
              
% fetch original video from AVI file
if ( exist( dat.ffn_vid, 'file' ) )
    dat.vid = read_dmri_video( dat.ffn_vid );
    dat.vid.dur = dat.vid.info.NumFrames/dat.vid.info.FramesPerSecond;
else
    error( '   Can''t find video file <%s>\n', dat.ffn_vid );
end

% find properties of audio & video streams:
%	2005a: 88 x 88 px  @ 23.78 fps
%	2005b: 68 x 68 px  @ 20.84 fps
%	2009a: 68 x 68 px  @ 22.41 fps
%	2010a: 68 x 68 px  @ 23.80 fps
%	2010b: 84 x 84 px  @ 23.80 fps
%	2010b: 68 x 68 px  @ 33.18 fps
nfv	= 0; frv = 0; vlv = 0;
if ~isempty(dat.vid.info)
    nfv	= dat.vid.info.NumFrames;
    frv	= dat.vid.info.FramesPerSecond;
    vlv	= nfv/frv;
end

% create video display axes
hvax = 8/9;
lvax = 0.08;    bvax = 1-hvax;
wvax = 1-2*lvax;	
dat.vid.pos	= [lvax bvax wvax hvax];

% image crop frame
dat.vid.frame.pos	= [0 0];
dat.vid.frame.size	= [0 0];

fVID = axes( 'Parent',pVID, 'Position',dat.vid.pos, 'Visible','on', 'FontSize',fss );
dat.vid.fVID = fVID;

% video time and frame selection controls
% create Navigaton panel
pNAV = uipanel(	'Parent',pVID, 'BackgroundColor',bgc, ...
                'FontSize',fsx, 'Position',[hoff voff 1-2*hoff 1-hvax]);

slMin	= 1;	slMax = dat.vid.info.NumFrames;
slRng	= (slMax-slMin);
slStpSm = slMin/slRng;	slStpLg = 5*slStpSm;
nav_width = 1-2*hoff;
nav_hght = 1-2*voff;
pbFRLT	= uicontrol( 'Parent',pNAV, 'Style','pushbutton', 'String','<', ...
                     'FontSize',fsm, 'ForegroundColor','black', 'Units','normalized',...
                     'Position',[hoff 0.5 nav_width/10 nav_hght/2-hoff],'Tag','FRLT',...
                     'Callback',{@update_frame_CB});
lbFRST	= uicontrol( 'Parent',pNAV, 'Style','Edit', 'String',num2str(dat.fstep), 'Tag', 'FRST',...
                     'Units','normalized','FontSize',fsl, ...
                     'Position',[hoff+nav_width/10 0.5 nav_width/10 nav_hght/2-hoff],...
                     'BackgroundColor',bgd,'Callback',{@update_framestep_CB});
pbFRRT	= uicontrol( 'Parent',pNAV, 'Style','pushbutton', 'String','>', ...
                     'FontSize',fsm, 'ForegroundColor','black', 'Units','normalized',...
                     'Position',[hoff+2*nav_width/10 0.5 nav_width/10 nav_hght/2-hoff],...
                     'Tag','FRRT','Callback',{@update_frame_CB});
slVINC	= uicontrol( 'Parent',pNAV, 'Style','Slider', 'BackgroundColor',bgd, ...
                     'Value',dat.fnum, 'Min',slMin, 'Max',slMax, 'SliderStep',[slStpSm slStpLg],...
                     'Units','normalized','Tag','FRSL', ...                 
                     'Position',[0.3+hoff 0.5+2*voff nav_width/2+nav_width/5-hoff nav_hght/3],	'Callback',{@update_frame_CB} );
lbVFRAM	= uicontrol( 'Parent',pNAV, 'Style','Edit', 'String','Frame No.:', 'Enable','Inactive', ...
                     'Units','normalized','FontSize',fsl, 'Position',[hoff voff nav_width/6 nav_hght/2] );
txVFRAM	= uicontrol( 'Parent',pNAV, 'Style','Edit', 'String',num2str(dat.fnum,'%d'), ...
                     'FontSize',fsl, 'Units','normalized','Position',[hoff+nav_width/6 voff nav_width/6 nav_hght/2], ...
                     'BackgroundColor',bgd,'Callback',{@edit_frameno_CB} );
lbVFTOT	= uicontrol( 'Parent',pNAV, 'Style','Edit', 'String',[' of ' num2str(slMax)], ...
                     'FontSize',fsl, 'Units','normalized','Position',[hoff+nav_width/3 voff nav_width/6 nav_hght/2], ...
                     'Enable','Inactive', 'HorizontalAlignment','left' );
lbVTIME	= uicontrol( 'Parent',pNAV, 'Style','Edit', 'String','Time (sec):', 'Enable','Inactive', ...
                     'Units','normalized','FontSize',fsl, 'Position',[hoff+nav_width/2 voff nav_width/6 nav_hght/2] );
txVTIME	= uicontrol( 'Parent',pNAV, 'Style','Edit', 'String',num2str(dat.time,'%.2f'), ...
                     'FontSize',fsl, 'Units','normalized','Position',[hoff+2*nav_width/3 voff nav_width/6 nav_hght/2], ...
                     'BackgroundColor',bgd,	'Callback',{@edit_time_CB} );
lbVTTOT	= uicontrol( 'Parent',pNAV, 'Style','Edit', 'String',[' of ' num2str(dat.vid.dur,'%.2f')], ...
                     'FontSize',fsl, 'Units','normalized','Position',[hoff+5*nav_width/6 voff nav_width/6 nav_hght/2], ...
                     'Enable','Inactive', 'HorizontalAlignment','left' );

                 
% Annotation controls        
pbLoad	= uicontrol( 'Parent',pANN, 'Style','pushbutton', 'String','Load', ...
                     'FontSize',fsh, 'ForegroundColor','black', 'Units','normalized',...
                     'Position',[hoff 1-voff-0.1 1-2*hoff 0.1],'Tag','load_button',...
                     'Callback',{@load_annotations_CB});
pbLoad	= uicontrol( 'Parent',pANN, 'Style','pushbutton', 'String','Save', ...
                     'FontSize',fsh, 'ForegroundColor','black', 'Units','normalized',...
                     'Position',[hoff 1-2*voff-0.2 1-2*hoff 0.1],'Tag','load_button',...
                     'Callback',{@save_annotations_CB});
puAnn = uicontrol( 'Parent',pANN, 'Style','popup', 'String',dat.annotation.labels, ...
                     'FontSize',fsh, 'ForegroundColor','black', 'Units','normalized',...
                     'Position',[hoff 1-3*voff-0.3 1-2*hoff 0.1],'Tag','popup_menu',...
                     'Callback',{@annotate_part_CB});
pbLoad	= uicontrol( 'Parent',pANN, 'Style','togglebutton', 'String','Edit', ...
                     'FontSize',fsh, 'ForegroundColor','black', 'Units','normalized',...
                     'Position',[hoff 1-4*voff-0.4 1-2*hoff 0.1],'Tag','edit_button',...
                     'Callback',{@activate_edit_mode_CB});
dat.annotation.popup_handle = puAnn;
                 
% CASY controls
pbLoad_CASY = uicontrol( 'Parent',pCASY, 'Style','togglebutton', 'String','Show', ...
                     'FontSize',fsh, 'ForegroundColor','black', 'Units','normalized',...
                     'Position',[hoff 1-voff-0.15 1/8 0.15],'Tag','show_casy_button',...
                     'Callback',{@show_casy_CB},'Visible','Off');
             
% rescale original video
if ( isstruct(dat.vid.avi) )
	nf = dat.vid.info.NumFrames;
    num_px	= dat.vid.info.Width;
    scale	= wd_vid/num_px;
    pixres	= wd_im * FOV/100 / num_px;
    dat.status.pixres = pixres;
    hWB     = waitbar(0,['Rescaling untracked AVI: processing ' num2str(nf) ' frames ...']);
    vid_rs	= repmat(struct('cdata', [], 'colormap', []), 1, nf);
    jf      = get( get(hWB,'JavaFrame'), 'FigurePanelContainer');
    jf.getComponent(0).getRootPane.getTopLevelAncestor.setAlwaysOnTop(1);
    for s = 1:nf
        if ~mod(s,wbi), waitbar(s/nf); end;	% report progress
        im	= dat.vid.avi(s).cdata;
        im	= im(:,:,1);
        map	= dat.vid.avi(s).colormap;
        im_ = imresize( im, scale );
        %if (~strncmp(computer(),'MAC',3) || (invert_video))
        if (invert_video)
            vid_rs(s).cdata     = flipud(im_);
        else
            vid_rs(s).cdata     = im_;
        end
        vid_rs(s).colormap  = map;
    end
    close(hWB);
    dat.vid.avi_	= dat.vid.avi;
    dat.vid.avi     = vid_rs;
    dat.vid.scale	= scale;
    dat.vid.dims	= size( dat.vid.avi(1).cdata );
    if (length(dat.vid.avi)>nf),	dat.vid.avi(nf+1:end)  = [];	end;
    if (length(dat.vid.avi_)>nf),	dat.vid.avi_(nf+1:end) = [];	end;
    % centre and resize crop frame to scaled data
    fdims	= floor( 0.8*dat.vid.dims );
    fcenter	= floor( dat.vid.dims/2 );
    dat.vid.frame.size	= fdims;
    dat.vid.frame.pos	= fcenter;
end

% reinitialize null data structure arrays now that signal lengths are known
nf = dat.vid.info.NumFrames;

% initialize data structure representing vocal tract parameters for each image frame
dat.vt(nf).imap	= [];	% intensity map
dat.vt(nf).Mxy	= [];   % coodinates of minima nodes in graph
dat.vt(nf).Mwt	= [];   % weightings of minima nodes in graph
dat.vt(nf).Ixy	= [];   % coodinates of inner inflection nodes in graph
dat.vt(nf).Iwt	= [];   % weightings of inner inflection nodes in graph
dat.vt(nf).Oxy	= [];   % coodinates of inner inflection nodes in graph
dat.vt(nf).Owt	= [];   % weightings of inner inflection nodes in graph
dat.vt(nf).cen	= [];   % tract center coordinates
dat.vt(nf).len	= [];   % tract length
dat.vt(nf).pts	= [];   % coordinates of tissue boundary points
dat.vt(nf).spts	= [];   % coordinates of smoothed boundary points
dat.vt(nf).man	= [];   % coordinates of hand-segmented boundary points
dat.vt(nf).aut	= [];   % coordinates of auto-segmented boundary points
dat.vt(nf).iix	= [];   % intensity indices
dat.vt(nf).int	= [];   % intensity of lf & rt tissue boundary thresholds
dat.vt(nf).AF	= [];	% tract area function
dat.vt(nf).AFx	= [];   % x-index of tract area function
dat.vt(nf).ngl	= [];	% number of gridlines used in frame
dat.vt(nf).dx	= [];	% horizontal image displacement from reference frame
dat.vt(nf).dy	= [];	% vertical image displacement from reference frame
dat.vt(nf).mod	= [];	% flag indicating that frame has been manually edited

% initialize image displacement for all frames in vocal tract
for fno = 1:nf
    dat.vt(fno).dx	= 0;	% horizonatal shift
    dat.vt(fno).dy	= 0;	% vertical shift
    dat.vt(fno).mod	= 0;	% no frames initially edited
end

% initialize audio & video centre frame at 2 sec
set( fVID, 'Visible','On' );
dat.status.showmloc = 0;
dat.status.mloc	= [40 40];
dat.status.ptpicked = 0;
set( hGUI, 'WindowButtonMotionFcn',{@image_hover_CB} );
guidata(hGUI,dat);

if isempty(dat.seg)
    update_time(2);
else
    lbl = list_labels_CB();
    lab = char(lbl(1));
    ff	= dat.seg.(lab).fint;
    fprintf('   Going to segment [%s]: frames %d to %d\n', lab,ff(1),ff(2) );
    set(txLABEL,'String',lab);
    guidata(hGUI,dat);
    update_frame_lims( ff );
end


% % don't assign output until all data structures have been updated
% waitfor( pbSAVE,'String','Saved' );
% if ( (nargout) && (dat.status.exit) )
%     dat	= guidata(hGUI);
% end
% 
% waitfor( pbEXIT,'String','Exiting' );
% if (dat.status.exit)
%     disp('   Exiting GUI.' );
%     close(hGUI);
% end
% disp(' ');



%------------------------------------
% subfunctions and callbacks
%------------------------------------
    function update_time( t )
        dat	= guidata(hGUI);
        dt = 1/dat.vid.info.FramesPerSecond;	% min timestep = length of video frame
        if (t >= dt) 
            fv	= floor( t * dat.vid.info.FramesPerSecond );
            dat.time	= t;
            dat.fnum	= fv;
            guidata(hGUI,dat);
        else
            fprintf('   Time must be specified in range [%0.3f .. %0.3f]\n', dt, dat.aud.dur-dt );
        end
        update_windows();
        refresh_GUI();
    end %update_time()

    function update_frame( f )
        dat	= guidata(hGUI);
        if (f >= 1) && (f <= dat.vid.info.NumFrames)
            t	= f/dat.vid.info.FramesPerSecond;
            dat.time	= t;
            dat.fnum	= f;
            set( slVINC,  'Value', f );
            set( txVFRAM, 'String',num2str(f) );
            set( txVTIME, 'String',num2str(t,'%.2f') );
            guidata(hGUI,dat);
        else 
            fprintf('   Frame number must be specified in range [1..%d]\n', dat.vid.info.NumFrames );
        end 
        update_windows();
        refresh_GUI();
    end %update_frame()

    function update_frame_lims( fint )
        dat	= guidata(hGUI);
        dat.fnum	 = round(mean(fint));
        dat.aud.wwd  = diff(fint/dat.vid.info.FramesPerSecond);
        guidata(hGUI,dat);
        update_frame(dat.fnum);
    end %update_frame_lims()

    function update_windows()
        dat	= guidata(hGUI);
        int	= [dat.time-dat.wwd/2 dat.time+dat.wwd/2];
        ff	= round( int .* dat.vid.info.FramesPerSecond );
        try
            if ff(1)<1,	ff(1)=1; end;
            if ff(2)>dat.vid.info.NumFrames, ff(2)=dat.vid.info.NumFrames; end;
        catch ME
            warning(ME.message);
        end
        dat.vid.fint	= ff;
        guidata(hGUI,dat);
    end %update_windows()

    function refresh_GUI( s,e )
        dat	= guidata(hGUI);
        t	= dat.time;

        % return focus to main GUI in case there are subplots active
        set(0,'CurrentFigure',dat.hGUI_id);
        
        im	= dat.vid.avi(dat.fnum).cdata;
        map	= dat.vid.avi(dat.fnum).colormap;
        h = subimage( im );
        set(h,'HitTest','off');

        % In case there is annotation information, show the vocal tract
        % contours
        if dat.annotation.show
            scale	= dat.vid.scale;
            if size(dat.annotation.contours(dat.fnum).points,1)>=1
                hold on
                x_coords = dat.annotation.contours(dat.fnum).points(:,1);
                y_coords = dat.annotation.contours(dat.fnum).points(:,2);
                n_points = size(x_coords,1);
                
                point_handles = zeros(n_points,1);
                ids = dat.annotation.contours(dat.fnum).label_ids(:);
                inactive_ids = ids~=dat.annotated_part_id;
                point_handles(inactive_ids) = plot(scale*x_coords(inactive_ids),...
                    scale*y_coords(inactive_ids),'.r','MarkerSize',8);
                active_ids = ids==dat.annotated_part_id;
                point_handles(active_ids) = plot(scale*x_coords(active_ids),...
                    scale*y_coords(active_ids),'.g','MarkerSize',8);
                
                for c = 1:length(point_handles)
                    set(point_handles(c),'HitTest','off');
                end
                
                %             h = plot(scale*x_coords,scale*y_coords);
                %             set(h, 'HitTest', 'off');
                dat.annotation.contours(dat.fnum).point_handles = point_handles;
                hold off
            end
        end
        
        if dat.edit_annotations
            set(gca, 'ButtonDownFcn', @edit_annotations);
        end
                
        % Show the CASY model
        if dat.casy.loaded && dat.casy.show
            hold on
            show_casy(dat.casy.parameters);
            hold off
        end

        
        % update video image axes
        imdim = size(im);
        imwd  = imdim(1);
        %axis off;


        guidata(hGUI,dat);
    end %refresh_GUI()

    % edit frame limits of audio signal display window
    function edit_wflim_CB( s,e )
        dat	= guidata(hGUI);
        flf	= round(str2double( get(txWFRML,'String')) );
        if flf<1, flf=1; end;
        frt	= round(str2double( get(txWFRMR,'String')) );
        if frt<=flf, frt=flf+diff(dat.vid.fint); end;
        if frt>dat.vid.info.NumFrames, frt = dat.vid.info.NumFrames; end;
        update_frame_lims( [flf frt] );
        guidata(hGUI,dat);
    end

    % export segment of video as AVI movie
    function export_video_CB( s,e )
        dat	= guidata(hGUI);
        fdn	= fullfile( home, dir_cpt, dir_mov );
        if ~(exist(fdn,'dir')), mkdir(fdn); end;
        lbl	= get(txLABEL,'String');
        fa	= dat.vid.fint(1);
        fb	= dat.vid.fint(2);
        nf	= fb-fa;
        fn	= [tok '_' lbl '_f' num2str(fa) '-f' num2str(fb) '.avi'];
        ffn	= fullfile( fdn, fn );
        int	= floor(fa:fb);
        if get( cbTMODE,'Value' )	% use tracked Video
            vobj = VideoReader(dat.ffn_trk);
        else                        % use original video 
            vobj = VideoReader(dat.ffn_vid);
        end
        hh	 = vobj.Height;
        ww	 = vobj.Width;
        mov(1:nf) = struct('cdata', zeros(hh,ww,3,'uint8'), 'colormap',[] );
        for f = 1:nf
            mov(f).cdata = read(vobj,fa+f);
        end
        fps = str2double(get(txVRate,'String'));
        fprintf( '   Writing video segment to file   <%s> (@ %0.1f fps)\n', fn,fps );
        movie2avi( mov,ffn, 'fps',fps, 'compression','none', 'quality',100 );
        %movie2avi( mov,ffn, 'fps',fps, 'compression','RLE', 'quality',100 );
    end %export_video_CB

    function export_audio_CB( s,e )
        dat = guidata(hGUI);
        fdn	= fullfile( home, dir_cpt, dir_wav );
        if ~(exist(fdn,'dir')), mkdir(fdn); end;
        lbl	= get(txLABEL,'String');
        fa	= dat.vid.fint(1);
        fb	= dat.vid.fint(2);
        fn	= [tok '_' lbl '_f' num2str(fa) '-f' num2str(fb) '.wav'];
        ffn	= fullfile( fdn, fn );
        disp([ '   Writing audio segment to file   <' fn '>' ]);
        sa	= dat.aud.sint(1);
        sb	= dat.aud.sint(2);
        int = floor(sa:sb);
        wavwrite( dat.aud.sig(int), dat.aud.Fs, 16, ffn );
    end %export_audio_CB

    % export all JPG images in specified range
    function export_frames_CB( s,e )
        dat	= guidata(hGUI);
        ff  = dat.vid.fint;
        fa  = min(ff);	fb = max(ff);
        lbl	= get(txLABEL,'String');
        dn	= [ dir_cpt '/' tok '_' lbl '_f' num2str(fa) ];
        disp([ '   Writing frames ' num2str(fa) ':' num2str(fb) ' to directory: <' dn '>' ]);
        jpgint = str2double(get(txFSTEP,'String'));
        for fr = fa:jpgint:fb
            export_frame( fr, lbl, dn, 1 );
        end
    end %export_frames_CB

    % export current frame
    function export_frame_CB( s,e )
        dat	= guidata(hGUI);
        lbl	= get(txLABEL,'String');
        %export_frame( dat.fnum, lbl, dir_cpt, 1 );
        export_gca( dat.fnum, lbl, dir_cpt, 1 );
    end %export_frame_CB

    % export JPG image of single AVI frame
    function export_frame( f, label, dir, vb )
        dat	= guidata(hGUI);
        fdn	= fullfile( home, dir_cpt, dir_jpg );
        if ~(exist(fdn,'dir')), mkdir(fdn); end;
        if get( cbTMODE,'Value' )       % Use tracked Video frame
            im	= dat.tvd.avi(f).cdata; im_ = im;
        else 
            if get( cbICORR,'Value' )	% Use intensity-corrected Video frame
                im	= dat.cvd.avi(f).cdata; im_ = im;
            else
                im	= dat.vid.avi(f).cdata; im_ = im;
            end
        end
        fn	= [ tok '_' label '_f' num2str(f) ];
        
        if get( cbSHOWF,'Value' )       % crop image if frame currently activated
            fn	= [ fn '_crp' ];
            clf = dat.vid.frame.pts.x(1);	crt = dat.vid.frame.pts.x(2);
            ctp = dat.vid.frame.pts.y(1);	cbt = dat.vid.frame.pts.y(3);
            % adjust limits if frame extends beyond image matrix limits
            cxa	= max([clf 1]);	cxb	= min([crt size(im_,2)]);
            cya	= max([ctp 1]);	cyb	= min([cbt size(im_,1)]);
            climx = round(cxa:cxb); climy = round(cya:cyb);
            im_ = im( climy, climx );
        end
        fn	= [ fn '.jpg' ];
        ffn	= fullfile( fdn, fn );
        if (vb),   disp([ '   Writing frame ' num2str(f) ' to image file: <' fn '>' ]); end;
        imwrite( im_, ffn, 'jpg' );
    end %export_frame

    % export JPG image of single AVI frame
    function export_gca( f, label, dir, vb )
        dat	= guidata(hGUI);
        fdn	= fullfile( home, dir_cpt, dir_jpg );
        if ~(exist(fdn,'dir')), mkdir(fdn); end;
        fr	= getframe(fVID);
        im	= fr.cdata; im_ = im;
        fn	= [ tok '_' label '_f' num2str(f) ];
        if get( cbSHOWF,'Value' )       % crop image if frame currently activated
            fn	= [ fn '_crp' ];
            clf = dat.vid.frame.pts.x(1);	crt = dat.vid.frame.pts.x(2);
            ctp = dat.vid.frame.pts.y(1);	cbt = dat.vid.frame.pts.y(3);
            % adjust limits if frame extends beyond image matrix limits
            cxa	= max([clf 1]);	cxb	= min([crt size(im_,2)]);
            cya	= max([ctp 1]);	cyb	= min([cbt size(im_,1)]);
            climx = round(cxa+1:cxb); climy = round(cya+2:cyb);
            im_ = im( climy, climx );
        end
        fn	= [ fn '.jpg' ];
        ffn	= fullfile( fdn, fn );
        if (vb),   disp([ '   Writing frame ' num2str(f) ' to image file: <' fn '>' ]); end;
        imwrite( im_, ffn, 'jpg' );
    end %export_gca

    % load annotation callback
    function load_annotations_CB ( s, e )
        dat = guidata(hGUI);
        [ann_file, pathname] = uigetfile('.mat','Please select the annotation file');
        contour_info = load(fullfile(pathname,ann_file));
        
        if isfield(contour_info, 'data') 
            if isfield(contour_info.data,'frame')
                if isfield(contour_info.data.frame{1},'modelEvolution')
                    contours = contour_info.data.frame;
                    n_frames = length(contour_info.data.frame);
                    n_pixels = dat.vid.info.Width;
                    dat.annotation.contours = struct('points',cell(1,n_frames),...
                        'label_ids',cell(1,n_frames), 'point_handles', cell(1,n_frames));                    
                    first_contour = contours{1}.modelEvolution{end};
                    seg_points = zeros(3,1);
                    for seg_counter = 1:3
                        seg = first_contour.segment{seg_counter};
                        seg_points(seg_counter) = size(seg.v,1);
                    end
                    n_points = sum(seg_points);
                    
                    for counter = 1:n_frames
                        frame_points = zeros(n_points,2);
                        frame_label_ids = zeros(n_points,1);
                        ind_start = 1;
                        n_unique_ids = 0;
                        
                        for seg_counter = 1:3 
                            ind_end = ind_start + seg_points(seg_counter) - 1;
                            x_coords = contours{counter}.modelEvolution{end}.segment{seg_counter}.v(:,1);
                            y_coords = contours{counter}.modelEvolution{end}.segment{seg_counter}.v(:,2);
                            x_coords = x_coords + 0.5*n_pixels;
                            y_coords = n_pixels-(y_coords + 0.5*n_pixels);

                            frame_points(ind_start:ind_end,:) = [x_coords y_coords];
                            label_ids = contours{counter}.modelEvolution{end}.segment{seg_counter}.i;

                            frame_label_ids(ind_start:ind_end,:) = label_ids+n_unique_ids;
                            n_unique_ids = n_unique_ids + length(unique(label_ids));
                            ind_start = ind_end+1;
                        end
                        dat.annotation.contours(counter).points = frame_points;
                        dat.annotation.contours(counter).label_ids = frame_label_ids;
                    end
                end
            end
            dat.annotation.dir = pathname;
            dat.annotation.file = ann_file;
            dat.annotation.show = true;
            guidata(hGUI,dat)
            refresh_GUI();
        elseif isfield(contour_info, 'contours')
            dat.annotation.contours = contour_info.contours;
            dat.annotation.dir = pathname;
            dat.annotation.file = ann_file;
            dat.annotation.show = true;
            if isfield(contour_info, 'labels')
                set(dat.annotation.popup_handle,'String',contour_info.labels);
            end
            guidata(hGUI,dat)
            refresh_GUI();            
        else
            error('Unrecognized data format');
        end
    end

    function save_annotations_CB(s, e)
        dat = guidata(hGUI);
        avi_fname = dat.ffn_vid;
        [pth, bname] = fileparts(avi_fname);
        [ann_file, pathname] = uiputfile(sprintf('%s_contours.mat', bname),'Please select the file where the contours will be saved');
        contours = dat.annotation.contours;
        labels = dat.annotation.labels;
        save(fullfile(pathname, ann_file), 'contours', 'labels');
    end

    function activate_edit_mode_CB( s, e)
        button_state = get(s, 'Value');
        dat = guidata(hGUI);
               
        if button_state == get(s, 'Max')
            dat.edit_annotations = true;
        elseif button_state == get(s, 'Min')
            dat.edit_annotations = false;
        end                        

        guidata(hGUI,dat);
        refresh_GUI();
    end      

    function edit_annotations ( s, e)
        dat = guidata(hGUI);
        scale	= dat.vid.scale;
        n_frames = dat.vid.info.NumFrames;

        cur_pnt = get(gca, 'CurrentPoint');
        coords = (1/scale)*cur_pnt(2, 1:2);
        
        mouse_side = get(gcf, 'SelectionType');
        
        if isfield(dat.annotation, 'contours')
            all_points = dat.annotation.contours(dat.fnum).points;
            ids = dat.annotation.contours(dat.fnum).label_ids;
            point_handles = dat.annotation.contours(dat.fnum).point_handles;
            active_ids = ids==dat.annotated_part_id;
        else
            dat.annotation.contours = struct('points',cell(1,n_frames),...
                'label_ids',cell(1,n_frames), 'point_handles', cell(1,n_frames));                    
            all_points = [];
            ids = [];
            point_handles = [];
            active_ids = [];
        end
                 
        if strcmp(mouse_side,'normal')
            p_handle = plot(coords(1), coords(2), 'g.');
            if isempty(all_points) 
                all_points(1,:) = coords;
                ids = dat.annotated_part_id;
                point_handles = p_handle;
            else
                all_points(end+1,:) = coords;
                ids(end+1) = dat.annotated_part_id;
                point_handles(end+1) = p_handle;
            end
        elseif strcmp(mouse_side, 'alt')
            if ~isempty(all_points)
                active_points = all_points(active_ids,:);
                if ~isempty(active_points)
                    ind = dsearchn(active_points, coords);                    
                    active_inds = find(active_ids);
                    ids(active_inds(ind)) = nan;
                end                        
            end
        end           
        point_ids = find(~isnan(ids));
        dat.annotation.contours(dat.fnum).points = all_points(point_ids,:);
        dat.annotation.contours(dat.fnum).point_handles = ...
            point_handles(point_ids);
        dat.annotation.contours(dat.fnum).label_ids = ids(point_ids);
        dat.annotation.show = true;
        guidata(hGUI, dat);
        refresh_GUI();
    end

    function plot_vt_outline( to_be_plotted, params)
        dat = guidata(hGUI);
        
        scale = dat.vid.scale;
        n_pixels = dat.vid.info.Width*scale;
        % plot outlines
        outline_scale = 1.6;
        outline_trans_x = -0.5*n_pixels-50;
        outline_trans_y = -0.5*n_pixels-50;
        u_outline_x = outline_scale*(n_pixels - (to_be_plotted.UpperOutline(:,1)))+outline_trans_x;
        u_outline_y = outline_scale*(n_pixels - to_be_plotted.UpperOutline(:,2))+outline_trans_y;
        
        plot( u_outline_x, u_outline_y, 'y-', 'LineWidth', 1.5 )
        % upper outline: red,  star  markers
        
        % subsecond plots will not destroy the pevious ones
        %                     set( allAxes.VTAxes, 'NextPlot', 'add' )
        hold on
        if ~isreal( to_be_plotted.BottomOutline )
            to_be_plotted.BottomOutline
        end
        b_outline_x = outline_scale*(n_pixels - (to_be_plotted.BottomOutline(:,1)))+outline_trans_x;
        b_outline_y = outline_scale*(n_pixels - to_be_plotted.BottomOutline(:,2))+outline_trans_y;        
        plot( b_outline_x, b_outline_y, 'y-','LineWidth', 1.5)
        % lower outline: cyan, cross markers
        
        
        if false
            for i = 1:size(gridLines,1)
                temp = squeeze(gridLines(i,:,:));
                plot( temp(1,:), temp(2,:), 'g', 'Parent', allAxes.VTAxes )
            end
        end
        
        if false
            global midPoints
            
            if ~isempty( midPoints )
                plot( midPoints(1,:), midPoints(2,:), 'o', 'Color', [1 0.65 0 ], ...
                    'MarkerFaceColor', [1 0.65 0 ], 'Parent', allAxes.VTAxes  )
            end
        end
        
        if false
            global crossSection
            
            if ~isempty( crossSection )
                delta = 0.5* crossSection';
                a = midPoints - delta;
                b = midPoints + delta;
                
                for i = 1:size( midPoints, 2 )
                    plot( [a(1,i) b(1,i)], [a(2,i) b(2,i)], '-', 'Color', [1 0.65 0 ], 'Parent', allAxes.VTAxes  )
                end
            end
            
            for i = GeoPts
                plot( i.coords(1), i.coords(2), 'go', 'MarkerSize',12, ...
                    'MarkerFaceColor','g', ...
                    'ButtonDownFcn', ['starts ' 'GEO_' i.Tags], 'Parent', allAxes.VTAxes )
            end
            
            for i = KeyPts
                plot( i.coords(1), i.coords(2), 'r+', 'MarkerSize',12, ...
                    'ButtonDownFcn', ['starts ' 'KEY_' i.Tags], 'Parent', allAxes.VTAxes )
            end
            
            global TUN_GridCenter
            
            plot(  TUN_GridCenter(1),  TUN_GridCenter(2), 'ms', 'MarkerSize', 6, ...
                'MarkerFaceColor','m', ...
                'ButtonDownFcn', [ 'starts ' 'TUN_GridCenter'], 'Parent', allAxes.VTAxes )
        end
        
    end

    

    % show CASY
    function show_casy( casy_parameters )        
        params = eall_casy(casy_parameters);
        to_be_plotted = outline_casy(params);
        plot_vt_outline(to_be_plotted, params)
        
    end

    % show CASY callback
    function show_casy_CB ( s, e )
        button_state = get(s, 'Value');
        dat = guidata(hGUI);
        
        % Check whether this is the first time that the articulatory model
        % is displayed, which means that we need to load all its parameters
        if ~dat.casy.loaded
            casy_parameters = init_casy();
            dat.casy.parameters = casy_parameters;
            dat.casy.loaded = true;
        end
        
        if button_state == get(s, 'Max')
            dat.casy.show = true;
        elseif button_state == get(s, 'Min')
            dat.casy.show = false;
        end                        

        guidata(hGUI,dat);
        refresh_GUI();
    end

    
    % choose a part to be active and annotated
    function annotate_part_CB ( s, e )
        dat = guidata(hGUI);
        dat.annotated_part_id =  get(s,'Value');
        guidata(hGUI,dat);
        refresh_GUI();
    end

    % update location after mouse click
    function mouse_down_CB( s,e )
        dat	= guidata(hGUI);
        obj	= gco();
        if (obj)
            
            obp = get( obj,'Parent');
            
            % mouse activity on audio windows:
            if (obj == fAUD) || (obj == fAUDw) || (obp == fAUD) || (obp == fAUDw)
                set( s,'WindowButtonUpFcn',@mouse_up_CB );
                pt	= get(gca,'CurrentPoint');
                t1	= pt(1,1);
                
            % mouse activity on sgram:
            elseif (obj == fFFT) || (obp == fFFT)
                if (obj == fFFT)
                    pt	= get(gco,'CurrentPoint');
                elseif (obp == fFFT)
                    pt	= get(gca,'CurrentPoint');
                end %if
                t = pt(1,1);
                dat.sg.tfmt = t;
                dat.sg.fmt	= 1;
                guidata(hGUI,dat);
                refresh_GUI();
                fetch_formants( dat.sg.tfmt, 1 );
                
            % mouse activity on video display:
            elseif (obp == fVID)
                [pt,Ip] = get_pixel_pt;
                if (dat.status.showmloc)    % if in analysis mode, fetch pt
                    dat.status.ptpicked = dat.status.ptpicked+1;
                    guidata(hGUI,dat);
                    fetch_correl(pt);
                else                        % otherwise just display pt
                    if ( get(cbSHVAX,'Value') && get(cbSHVMM,'Value') )	% in mm
                        fprintf('   x = %0.1fmm\ty = %0.1fmm\t\tIp = %d\n', pt(1),pt(2),Ip );
                    else                                                % or px
                        fprintf('   x = %d\ty = %d\t\tIp = %d\n', pt(1),pt(2),Ip );
                    end
                    guidata(hGUI,dat);
                    refresh_GUI();
                end
            end %if
            
        end %if
        
        function mouse_up_CB( s,e )
            pt2 = get(gca,'CurrentPoint');
            obj	= gco();
            obp = get( obj,'Parent');
            if ( ~sum(pt2(:,3)) & (obj ~= fFFT) & (obp ~= fFFT) )
                t2	= pt2(1,1);
                if (t2 == t1)	% user has selected single point in time
                    update_time( t1 );
                else            % user has dragged mouse over interval to measure duration
                    if (t2 < t1)
                        tt = t1; t1 = t2; t2 = tt;
                    end
                    tc = mean([t1 t2]);
                    tw = diff([t1 t2]);
                    dat.aud.wwd = tw;
                    set(txWINWD,'String',num2str(tw,'%.1f'));
                    guidata(hGUI,dat);
                    update_time( tc );
                end %if
            end %if
        end %mouse_up_CB

    end %mouse_down_CB

    % mouse click on gridline: display intensity profile and minima
    function gridline_sel( s,e )
        dat	= guidata(hGUI);
        gl	= get(s,'UserData');
        fr	= dat.fnum;
        fn	= ['Intensity profile: GL ' num2str(gl) ', Fr ' num2str(fr)];
        hh	= findobj('Name',fn);
        if isempty( hh )
            fpos = get(hGUI,'Position');
            ht	 = fpos(4)*0.75; wd = ht;	fpos(3) = wd;  fpos(4) = ht;
            fpos = fpos + [-wd-9 ht/2 0 0];
            if ~isempty( dat.vt(fr).imap )
                hPLT = figure( 'Name',fn );
                set( hPLT, 'Position',fpos, 'ToolBar','none', 'MenuBar','none' ); hold on;
                set( hPLT, 'HitTest','on', 'ButtonDownFcn',{@iplot_click,fr,gl} );
                xlabel('Position on Gridline (pixel no.)');
                ylabel('Pixel Intensity');
                imap = dat.vt(fr).imap(gl).imap;
                imin = dat.vt(fr).imap(gl).minix;
                imax = dat.vt(fr).imap(gl).maxix;
                hPRO = plot( 1:length(imap),imap,'b-', 'LineWidth',2 );
                plot( imin,imap(imin),'k.', 'MarkerSize',15 );
                plot( imax,imap(imax),'r.', 'MarkerSize',15 );
                axis tight; yy = ylim; ylim([0 yy(2)]);
                set( hPRO, 'HitTest','on', 'ButtonDownFcn',{@intensityprofile_click,gl} );
                if ~isempty( dat.vt(fr).iix )
                    if (gl <= length(dat.vt(fr).iix))
                        iix	 = dat.vt(dat.fnum).iix(gl);
                        int	 = dat.vt(dat.fnum).int(gl);
                        if (~isempty( iix.lf ) && ~isempty( int.lf ))
                            plot( iix.lf,int.lf,'go',   'MarkerSize',13 );
                        end
                        if (~isempty( iix.rt ) && ~isempty( int.rt ))
                            plot( iix.rt+1,int.rt,'go', 'MarkerSize',13 );
                        end
                        if ~isempty( int.lf )
                            line( xlim,[int.lf int.lf], 'Color','g', 'LineStyle','--', 'LineWidth',2 );
                        elseif ~isempty( int.rt )
                            line( xlim,[int.rt int.rt], 'Color','g', 'LineStyle','--', 'LineWidth',2 );
                        end
                    end
                end
            else
                nmin = get_grid_profiles(fr);
                gridline_sel( s,e );
            end
            if get( cbTMODE,'Value' )	% Use tracked tissue boundaries
                fn	= ['Tracked Video Intensity profile: GL ' num2str(gl) ', Fr ' num2str(fr)];
                hPL2 = figure( 'Name',fn );
                fpos = fpos + [0 -ht/3 0 0];
                set( hPL2, 'Position',fpos, 'ToolBar','none', 'MenuBar','none' ); hold on;
                xlabel('Position on Gridline (pixel no.)');
                ylabel('Pixel Intensity');
                glx	= dat.vid.grid.pts(gl).xx;
                gly	= dat.vid.grid.pts(gl).yy;
                pxr	= dat.tvd.avi(dat.fnum).cdata(:,:,1);
                pxg	= dat.tvd.avi(dat.fnum).cdata(:,:,2);
                pix = 1:length(glx);
                for p = pix
                    rmap(p) = pxr( glx(p),gly(p) );
                    gmap(p) = pxg( glx(p),gly(p) );
                end;
                plot( pix,rmap,'r-', 'LineWidth',2 ); hold on
                plot( pix,gmap,'g-', 'LineWidth',2 );
            end
        else % intensity profile plot already displayed: show location on profile
            pt	= round( get(gca,'CurrentPoint') );	pt = pt(1,1:2);
            xx	= dat.vid.grid.pts(gl).xx;
            ix	= find( xx == pt(1) );
            if length(ix)>1
                yy	= dat.vid.grid.pts(gl).yy;
                ix	= find( yy == pt(2) );
            end;
            imap = dat.vt(fr).imap(gl).imap;
            plot( get(hh,'CurrentAxes'), ix,imap(ix), 'r+', 'MarkerSize',13 );
        end
    end %gridline_sel

    % mouse click on intensity profile plot: show corresponding pixel
    function intensityprofile_click( s,e,gix )
        dat	= guidata(hGUI);
        pt	= round( get(gca,'CurrentPoint') );
        xix	= pt(1,1);
        row	= dat.vid.grid.pts(gix).yy;
        col	= dat.vid.grid.pts(gix).xx;
        cx	= col(xix); cy = row(xix);
        line( [cx-dcr cx+dcr], [cy cy], [0 0], 'Color','r', 'Parent',fVID );
        line( [cx cx], [cy-dcr cy+dcr], [0 0], 'Color','r', 'Parent',fVID );
    end %intensityprofile_click


    % mouse click on intensity profile plot: plot difference function
    function iplot_click( s,e, f,gl )
        dat	= guidata(hGUI);
        fn	= ['Intensity difference: GL ' num2str(gl) ', Fr ' num2str(f)];
        fpos = get(hGUI,'Position');
        ht	 = fpos(4)*0.75; wd = ht;	fpos(3) = wd;  fpos(4) = ht;
        fpos = fpos + [-wd-9 ht/3 0 0];
        hDPL = figure( 'Name',fn );
        set( hDPL, 'Position',fpos, 'ToolBar','none', 'MenuBar','none' ); hold on;
        xlabel('Position on Gridline (px)');
        ylabel('dIntensity');
        imap = dat.vt(f).imap(gl).imap;
        pix  = 2:length(imap);
        dint = diff(imap);
        plot( pix,0, 'k--', 'LineWidth',1 ); hold on;	% x-axis
        hDIF = plot( pix,dint, 'b-',  'LineWidth',2 );
        set( hDIF, 'HitTest','on', 'ButtonDownFcn',{@dintprofile_click, f,gl,s} );
        axis tight;
    end %iplot_click

    % click on intensity diff function plot: show corresponding px and inflection
    function dintprofile_click( s,e, f,gl,hIPlot )
        dat	= guidata(hGUI);
        dc  = 2;
        hAx = get(hIPlot,'CurrentAxes');
        pt	= round( get(gca,'CurrentPoint') );
        xix	= pt(1,1);
        % show corresponding pixel on video image on main GUI
        row	= dat.vid.grid.pts(gl).yy;
        col	= dat.vid.grid.pts(gl).xx;
        cx	= col(xix); cy = row(xix);
        if get( cbSHIFT,'Value' )
            cx	= cx + dat.vt(f).dx;
            cy	= cy + dat.vt(f).dy;
        end
        line( [cx-dcr cx+dcr], [cy cy], [0 0], 'Color','r', 'Parent',fVID );
        line( [cx cx], [cy-dcr cy+dcr], [0 0], 'Color','r', 'Parent',fVID );
        % show corresponding inflection pt on intensity plot
        imap = dat.vt(f).imap(gl).imap;
        cy	 = imap(xix);
        line( [xix-dc xix+dc], [cy cy], [0 0], 'Color','r', 'Parent',hAx );
        line( [xix xix], [cy-dc cy+dc], [0 0], 'Color','r', 'Parent',hAx );
    end %dintprofile_click

    % update the frame step used for video display
    function update_framestep_CB( s, e)
        dat = guidata(hGUI);
        new_fstep = round(str2double(get(s,'String')));
        set(s,'String',num2str(new_fstep));
        dat.fstep = new_fstep;
        guidata(hGUI, dat);
    end

    % update position after GUI control increment/decrement
    function update_frame_CB( s,e )
        dat = guidata(hGUI);
        ctag = get(s,'Tag');
        switch ctag
            case 'FRSL'
                fr = round( get(slVINC,'Value') );
            case 'FRLT'
                fr = max(1,dat.fnum-dat.fstep);
            case 'FRRT'
                fr = min(dat.vid.info.NumFrames,dat.fnum+dat.fstep);
        end
        update_frame( fr );
    end %update_frame_CB

    % move to left frame of currently displayed interval
    function frame_left_CB( s,e )
        dat = guidata(hGUI);
        update_frame( dat.vid.fint(1) );
    end %frame_left_CB

    % move to right frame of currently displayed interval
    function frame_right_CB( s,e )
        dat = guidata(hGUI);
        update_frame( dat.vid.fint(2) );
    end %frame_right_CB

    function edit_frameno_CB( s,e )
        dat = guidata(hGUI);
        fs	= get(txVFRAM,'String');
        f	= round(str2double(fs));
        update_frame(f);
    end %edit_frameno_CB

    function edit_time_CB( s,e )
        dat = guidata(hGUI);
        t	= str2double( get(txVTIME,'String') );
        update_time(t);
    end %edit_time_CB

    % play audio/video segments
    function play_aud_CB( s,e )
        dat = guidata(hGUI);
        int = dat.aud.sint(1):dat.aud.sint(2);
        soundsc( dat.aud.sig(int), dat.aud.Fs );
    end %play_aud_CB

    function play_vid_CB( s,e )
        dat	= guidata(hGUI);
        fps = str2double(get(txVRate,'String'));
        if ( (fps<1) || (fps>100) )
            fps = dat.vid.info.FramesPerSecond;
            fprintf('   Select a playback framerate between 1 and 100. (Original video rate = %0.2f fps)\n', fps );
        end 
        int	= dat.vid.fint(1):dat.vid.fint(2);
        if get( cbTMODE,'Value' )	% Show tracked Video
            mov = dat.tvd.avi(int);
        else 
            if get( cbICORR,'Value' )	%  Show intensity-corrected video
                mov_ = dat.cvd.avi(int);
            else
                mov_ = dat.vid.avi(int);
            end
            mov = mov_;
            for i=1:length(int)
                if strcmp(dat.status.OS,'lin')
                    mov(i) = im2frame( mov_(i).cdata, gray );
                else
                    mov(i) = im2frame( flipud(mov_(i).cdata), gray );
                end
            end
        end
        movie( fVID, mov, 1, fps, [0 0 0 0] );
    end %play_vid_CB

    function fmt = fetch_formants( t, vb ) % fetch formant vector at time t (secs)
        dat = guidata(hGUI);
        Fs	= dat.aud.Fs;
        fmt	= fetch_fmts( dat.aud.sig, Fs, t, 1 );
        if (vb)
            fprintf('   Time %2.2f\t%5.1f\t%5.1f\t%5.1f\t%5.1f Hz\n', t,fmt(1),fmt(2),fmt(3),fmt(4) );
        end
    end %fetch_formants

	% display crop frame on video window
    function toggle_frame_CB( s,e )
        dat = guidata(hGUI);
        if get( cbSHOWF,'Value' )
            set_frame();
        end
        refresh_GUI();
    end %toggle_frame_CB

    % define crop frame
    function set_frame()
        dat	= guidata(hGUI);
        lf	= dat.vid.frame.pos(1) - dat.vid.frame.size(1)/2;
        rt	= dat.vid.frame.pos(1) + dat.vid.frame.size(1)/2;
        tp	= dat.vid.frame.pos(2) - dat.vid.frame.size(2)/2;
        bt	= dat.vid.frame.pos(2) + dat.vid.frame.size(2)/2;
        % define frame borders
        dat.vid.frame.pts.x	= [lf rt rt lf lf];
        dat.vid.frame.pts.y	= [tp tp bt bt tp];
        guidata(hGUI,dat);
    end %set_frame 

	% display analysis grid on video window
    function toggle_grid_CB( s,e )
        dat = guidata(hGUI);
        if get( cbSHOWG,'Value' )
            if ~dat.vid.grid.nlines
                place_grid();
            end
        end
        refresh_GUI();
    end %toggle_grid_CB

    function toggle_vidaxes_CB( s,e )
        if get( cbSHVAX,'Value' )
            set( lbSCVAX, 'Visible','on' );
            set( lbSHVMM, 'Visible','on' );
            set( cbSCVAX, 'Visible','on' );
            set( cbSHVMM, 'Visible','on', 'Enable','On' );
            if get( cbSHVMM,'Value' )
                set( lbSHFOV, 'Visible','on' );
                set( txSHFOV, 'Visible','on', 'Enable','On' );
                set( cbSCVAX, 'Value',0 );
                set( cbSCVAX, 'Enable','Inactive' );
            else
                set( lbSHFOV, 'Visible','off' );
                set( txSHFOV, 'Visible','off', 'Enable','Inactive' );
                set( cbSCVAX, 'Enable','On' );
            end
        else
            set( lbSCVAX, 'Visible','off' );
            set( cbSCVAX, 'Visible','off', 'Enable','Inactive' );
            set( lbSHVMM, 'Visible','off' );
            set( cbSHVMM, 'Visible','off', 'Enable','Inactive' );
            set( lbSHFOV, 'Visible','off' );
            set( txSHFOV, 'Visible','off', 'Enable','Inactive' );
        end
        refresh_GUI();
    end %toggle_vidaxes_CB

	% overlay Ohman analysis grid on video frame
    function place_grid( s,e )
        dat	= guidata(hGUI);
        if (dat.vid.grid.nlines)
            fprintf('\n   An analysis grid has already been configured ... \n');
            action	= input('   Clear grid and start again? (y/n):','s');
            if strcmpi(action,'n') || strcmpi(action,'no'),
                fprintf('\n   Leaving existing analysis grid unchanged. \n\n' );
                return;
            end;
        end;
        % clear any existing grid
        dat.vid.grid.nlines	= 0;
        dat.vid.grid.ends	= [];
        set( cbSHOWG,'Value',1 )
        guidata(hGUI,dat);
        refresh_GUI();
        fprintf( '\n   Select glottis ... ' );
        dat.vid.grid.glot	= ginput(1);
        fprintf( '\n   Select highest point on palate ... ' );
        [xpal,ypal]         = ginput(1);
        dat.vid.grid.mpal	= round([xpal ypal]);
        dat.vid.grid.rad1	= round( dat.vid.grid.glot(1) - xpal );
        set( txGRAD, 'String', num2str(round(dat.vid.grid.rad1)) );
        fprintf( '\n   Select point on alveolar ridge ... ' );
        dat.vid.grid.dent	= round( ginput(1) );
        fprintf( '\n   Select mid-labial limit of tract ... ' );
        dat.vid.grid.mlab	= round( ginput(1) );
        % locate 1st (lingual) origin: directly below mid-palatal pt
        dy	= dat.vid.grid.wid/6;
        cx	= round( xpal );
        cy	= round( ypal + dat.vid.grid.rad1 + dy );
        dat.vid.grid.or1 = [cx cy];
        set( txOR1X, 'String',num2str(cx) );
        set( txOR1Y, 'String',num2str(cy) );
        % locate 2nd origin: extend line from OR1 through alveolar ridge
        tht	= atan( (cy-dat.vid.grid.dent(2))/(cx-dat.vid.grid.dent(1)) );
        dat.vid.grid.theta = tht;
        dR	= dat.vid.grid.wid/4;
        RR	= dat.vid.grid.rad1 + dat.vid.grid.wid/2 +dR;
        ax	= round( cx - RR*cos(tht) );
        ay	= round( cy - RR*sin(tht) );
        dat.vid.grid.or2 = [ax ay];
        set( txOR2X, 'String',num2str(ax) );
        set( txOR2Y, 'String',num2str(ay) );
        dat.vid.grid.rad2  = norm( dat.vid.grid.or1 - dat.vid.grid.or2 ) - dat.vid.grid.rad1;
        fprintf( '\n   All anatomical landmarks acquired ... setting analysis grid\n\n' );
        guidata(hGUI,dat);
        set_grid(1);
        refresh_GUI();
    end %place_grid

    % extract AF from frame callback
    function extract_AF_frame_CB( s,e )
        dat	= guidata(hGUI);
        extract_AF_frame( dat.fnum, 1 );
    end %extract_AF_frame_CB

    % extract AF from frame
    function extract_AF_frame( f, vb )
        dat	= guidata(hGUI);
        if isempty( dat.vt(f).imap )
            find_tissue_bnds( f );
            dat.vt(f).spts	  = [];
            dat.vt(f).spts.lf = [];
            refresh_GUI();
        end
        pts	= dat.vt(f).pts;
        cn_ = dat.vid.grid.glot;
        ngl	= dat.vt(f).ngl;
        AF	= zeros(1,ngl);
        len	= 0;
        AFx	= AF;
        for gl = 1:ngl
            cen	= dat.vt(f).cen(gl).xy;
            if get( cbSHIFT,'Value' )
                cen	= cen + [dat.vt(f).dx dat.vt(f).dy];
            end
            len	= len + norm( cn_-cen )*pixres/dat.vid.scale;
            ap	= norm( pts(gl).lf - pts(gl).rt )*pixres/dat.vid.scale;
            AF(gl)  = ap;
            AFx(gl) = len;
            if (vb>1), fprintf( '   Gridline %d:\tVT len = %3.0f\taperture = %2.0f\n', gl,len,ap ); end;
            cn_ = cen;
        end
        dat.vt(f).AF	= AF;
        dat.vt(f).AFx	= AFx;
        dat.vt(f).VTA	= AF(1:end-1) * diff(AFx)';
        fprintf('   Frame %d:\tVT Area = %0.1f\n', f,dat.vt(f).VTA );
        if (vb>0),
            hAF = figure;
            fig_pos = get(hGUI,'Position') + [0 gui_ht+29 0 0];
            fig_pos(3) = gui_wd/2;  fig_pos(4) = gui_wd/9;
            set( hAF, 'Name',['Vocal Tract Area Function: Frame ' num2str(f)], 'Position', fig_pos );
            set( hAF, 'ToolBar','none', 'MenuBar','none', 'WindowButtonDownFcn',{@AF_mouse_down_CB, f} );
            stairs( AFx,AF );
            %axis('equal');
            axis tight; yy = ylim; ylim([0 yy(2)]);
            %xlabel('Distance from Glottis (mm)');
            ylabel('Aperture (mm)');
        end
        guidata(hGUI,dat);
    end %extract_AF_frame

    % handle mouse activity on area function plot
    function AF_mouse_down_CB( s,e, fnum )
        dat	= guidata(hGUI);
        hAX	= gca();
        pt	= get( hAX,'CurrentPoint' );
        d	= pt(1,1);
        AFx	= dat.vt(fnum).AFx;
        gix	= find( AFx >= d, 1,'first' ) - 1;
        ap	= dat.vt(fnum).AF(gix);
        fprintf('   Dist from glottis = %0.1f:\tVocal Tract aperture = %0.1f\n', d,ap );
        % indicate which gridline this section of AF was calculated at
        grd = dat.vid.grid.ends(gix);
        set(0,'CurrentFigure',hGUI);
        set(hGUI,'CurrentAxes',fVID);
        line( grd.x, grd.y, [0 0], 'Color','r', 'Linewidth',dat.vid.grid.linewd );
    end %corr_mouse_down_CB

    % save segment callback
    function save_segment_CB( s,e, mode )
        dat	= guidata(hGUI);
        f	= dat.fnum;
        pts	= dat.vt(f).pts;
        if strcmp( mode,'manual' )
            dat.vt(f).man = pts;
        else
            dat.vt(f).aut = pts;
        end
        guidata(hGUI,dat);
    end %save_segment_CB

    % compare area function callback
    function compare_AF_CB( s,e )
        dat	= guidata(hGUI);
        f	= dat.fnum;
        if isempty( [dat.vt(f).man] )
            fprintf( '   No manual boundary data for frame %d:\t segment and save\n', f );
        elseif isempty( [dat.vt(f).aut] )
            fprintf( '   No auto boundary data for frame %d:\t segment and save\n', f );
        else
            compare_AF( f, 1 );
        end
    end %extract_AF_frame_CB

    % compare area function
    function compare_AF( f, vb )
        dat	= guidata(hGUI);
        man	= dat.vt(f).man;
        aut	= dat.vt(f).aut;
        cn_ = dat.vid.grid.glot;
        ngl	= dat.vt(f).ngl;
        AFm	= zeros(1,ngl);
        AFa	= zeros(1,ngl);
        len	= 0;
        AFx = AFm;
        for gl = 1:ngl
            cen	= dat.vt(f).cen(gl).xy;
            len	= len + norm( cn_-cen )*pixres/dat.vid.scale;
            AFm(gl)	= norm( man(gl).lf - man(gl).rt )*pixres/dat.vid.scale;
            AFa(gl)	= norm( aut(gl).lf - aut(gl).rt )*pixres/dat.vid.scale;
            AFx(gl) = len;
            if (vb>1), fprintf( '   AFm(%d) = %3.0f:\tAFa(%d) = %0.0f\n', gl,AFm(gl),gl,AFa(gl) ); end;
            cn_ = cen;
        end
        dat.vt(f).AFm = AFm;
        dat.vt(f).AFa = AFa;
        dat.vt(f).AFx = AFx;
        err = AFm-AFa;
        if (vb>0),
            hAF = figure;
            fig_pos = get(hGUI,'Position') + [0 gui_ht+29 0 0];
            fig_pos(3) = gui_wd/2;  fig_pos(4) = gui_wd/9;
            set( hAF, 'Name',['Vocal Tract Area Functions: Frame ' num2str(f)], 'Position', fig_pos );
            set( hAF, 'ToolBar','none', 'MenuBar','none', 'WindowButtonDownFcn',{@AF_mouse_down_CB, f} );
            stairs( AFx,AFm, 'b-' ); hold on
            stairs( AFx,AFa, 'r-' );
            stairs( AFx,err, 'k:' );
            %axis('equal');
            axis tight; yy = ylim; ylim([0 yy(2)]);
            %xlabel('Distance from Glottis (mm)');
            ylabel('Aperture (mm)');
        end
        guidata(hGUI,dat);
    end %compare_AF

    % plot lingual movement callback
    function plot_lingualmovt_CB( s,e )
        dat	= guidata(hGUI);
        seg = get(txLABEL,'String');
        if isfield(dat.seg, seg)
            fint  = dat.seg.(seg).fint;
            fprintf('   Plotting lingual movement for segment [%s]: frames %d to %d\n', seg,fint(1),fint(2) );
            plot_lingualmovt( fint, 1 );
        else
            fprintf('   Can find lingual data for specified segment: %s\n', seg );
        end
    end %plot_lingualmovt_CB

    % plot lingual movement over selected interval
    function plot_lingualmovt( fint, vb )
        dat	= guidata(hGUI);
        if (vb>0),
            hLM = figure; hold on;
            set( hLM, 'Name',['Lingual Movement: Frames ' num2str(fint(1)) ':' num2str(fint(2))] );
            set( hLM, 'ToolBar','none', 'MenuBar','none', 'WindowButtonDownFcn',{@LM_mouse_down_CB,fint} );
        end
        tgl = dat.vid.grid.tng;     % range of lingual gridlines
        ff	= fint(1):fint(2);
        nf	= fint(2)-fint(1)+1;
        fm	= round(mean(fint));
        im	= dat.vid.avi(fm).cdata;
        map	= dat.vid.avi(fm).colormap;
        imshow( im,map );
        cix = 1;
        %colspec	= 'winter';
        colspec	= 'autumn';
        cmap	= eval(['flipud( ' colspec '(' num2str(nf) '))']);
        for f = ff
            if ~isempty( dat.vt(f).pts )
                tongue	= [dat.vt(f).pts(tgl).lf];
                line( tongue(1:2:end),tongue(2:2:end), 'Color',cmap(cix,:) );
            end
            cix = cix+1;
        end
        axis('equal');
    end %plot_lingualmovt

    % get palate trace callback
    function edit_boundaries_CB( s,e )
        dat	= guidata(hGUI);
        set( cbSHOWB, 'Value',0 );
        refresh_GUI();
        edit_boundaries( dat.fnum );
    end %edit_boundaries_CB


    % go to specified segment, or specify frame limits of segment if it doesn't already exist
    function goto_frame_CB( s,e )
        dat	= guidata(hGUI);
        lab = get(txLABEL,'String');
        if ( isempty(lab) || ~isempty(regexpi(lab,'\W')) )
            fprintf('   Segment labels must be non-empty alphanumeric strings ... \n\n' );
        else
            if isfield(dat.seg, lab)
                ff	= dat.seg.(lab).fint;
                fprintf('   Going to segment [%s]: frames %d to %d\n', lab,ff(1),ff(2) );
                update_frame_lims( ff );
            else
                f1 = str2double( get( txWFRML,'String' ) );
                f2 = str2double( get( txWFRMR,'String' ) );
                dat.seg.(lab).fint = [f1 f2];
                fprintf('   Frames %d to %d labelled: [%s]\n', f1,f2,lab );
            end
        end
        guidata(hGUI,dat);
    end %goto_frame_CB

    % list all tagged segments
    function lbl = list_labels_CB( s,e )
        dat	= guidata(hGUI);
        if isempty(dat.seg)
            fprintf('\n   No segment label(s) defined.\n' );
        else
            lbl	= fieldnames(dat.seg);
            nlab = length(lbl);
            fprintf('\n   %d segment label(s) defined:\n', nlab );
            for lno = 1:nlab
                lab	= lbl{lno};
                ff	= dat.seg.(lab).fint;
                fprintf('   Segment [%s]: frames %d to %d\n', lab,ff(1),ff(2) );
            end
            fprintf('\n\n');
        end
    end %list_labels_CB

    % clear segment label and associated data structures
    function clear_label_CB( s,e )
        dat	= guidata(hGUI);
        lab = get(txLABEL,'String');
        if isfield(dat.seg, lab)
            prompt	= sprintf( '   Delete segment [%s] and associated data structures [y/n]?  ', lab );
            action	= input( prompt, 's');
            if strcmpi( action,'y' ) || strcmpi( action,'yes' )
                ff	= dat.seg.(lab).fint;
                fprintf('   Clearing segment label [%s]: frames %d to %d\n', lab,ff(1),ff(2) );
                dat.seg = rmfield(dat.seg,lab);
                set( txLABEL, 'String','' );
            end
        else
            fprintf('   Cannot find segment with label [%s]\n', lab );
        end
        guidata(hGUI,dat);
    end %clear_label_CB

    % extract AFs across segment
    function analyze_segment( s,e )
        vb = 1;
        dat = guidata(hGUI);
        seg = get(txLABEL,'String');
        if (regexp(seg, '[\s\[\]]')),
            fprintf('\n   Enter alphanumeric segment name (*%s) and try again.', seg );
            return;
        elseif isfield(dat.seg, seg),
            %fprintf('\n   Data for segment (%s) already exisits ... \n', seg );
            %action	= input('   Overwrite? (y/n):','s');
            %if strcmpi(action,'n') || strcmpi(action,'no'),
            %    fprintf('\n   Leaving data for segment (%s) unchanged. \n\n', seg );
            %    return;
            %else
                dat.seg.(seg) = [];
            %end;
        end;
        fa	= dat.vid.fint(1);
        fb	= dat.vid.fint(2);
        fprintf('\n   Fetching vocal tract boundries for segment [%s]:\n', seg );
        for fn = fa:fb
            extract_AF_frame( fn, 0 );
            if (vb), fprintf( '   Frame %d:\tVT len = %3.0f\n', fn,dat.vt(fn).len ); end;
        end
        if (vb), fprintf( '   \n' ); end;
        dat.seg.(seg).fint = [fa fb];
        guidata(hGUI,dat);
    end %analyze_segment

    % set crop frame size
    function set_frm_size_CB( s,e )
        dat	= guidata(hGUI);
        dat.vid.frame.size(1) = str2double(get(txFRMW, 'String'));
        dat.vid.frame.size(2) = str2double(get(txFRMH, 'String'));
        guidata(hGUI,dat);
        set_frame();
        refresh_GUI();
    end %set_frm_size_CB
    % set crop frame position
    function set_frm_pos_CB( s,e )
        dat	= guidata(hGUI);
        dat.vid.frame.pos(1) = str2double(get(txFOFFx, 'String'));
        dat.vid.frame.pos(2) = str2double(get(txFOFFy, 'String'));
        guidata(hGUI,dat);
        set_frame();
        refresh_GUI();
    end %set_frm_pos_CB
    % shift crop frame size
    function alt_frm_size_CB( s,e, dx,dy )
        dat	= guidata(hGUI);
        dat.vid.frame.size(1) = str2double(get(txFRMW, 'String'))+dx;
        dat.vid.frame.size(2) = str2double(get(txFRMH, 'String'))+dy;
        set( txFRMW, 'String',num2str(dat.vid.frame.size(1)) );
        set( txFRMH, 'String',num2str(dat.vid.frame.size(2)) );
        guidata(hGUI,dat);
        set_frame();
        refresh_GUI();
    end %alt_frm_size_CB
    % shift crop frame position
    function alt_frm_pos_CB( s,e, dx,dy )
        dat	= guidata(hGUI);
        dat.vid.frame.pos(1) = str2double(get(txFOFFx, 'String'))+dx;
        dat.vid.frame.pos(2) = str2double(get(txFOFFy, 'String'))+dy;
        set( txFOFFx, 'String',num2str(dat.vid.frame.pos(1)) );
        set( txFOFFy, 'String',num2str(dat.vid.frame.pos(2)) );
        guidata(hGUI,dat);
        set_frame();
        refresh_GUI();
    end %alt_frm_pos_CB

    % set Ohman grid radius
    function set_grad_CB( s,e )
        dat	= guidata(hGUI);
        r1	= str2double(get(txGRAD,'String'));
        r2	= norm( dat.vid.grid.or1-dat.vid.grid.or2 ) - r1;
        dat.vid.grid.rad1 = round(r1);
        dat.vid.grid.rad2 = round(r2);
        guidata(hGUI,dat);
        set_grid(1);
        refresh_GUI();
    end %set_grad_CB
    % increment/decrement Ohman grid radius
    function alt_grad_CB( s,e, dr )
        dat	= guidata(hGUI);
        r1	= str2double(get(txGRAD,'String'))+dr;
        r2	= norm( dat.vid.grid.or1-dat.vid.grid.or2 ) - r1;
        dat.vid.grid.rad1 = round(r1);
        dat.vid.grid.rad2 = round(r2);
        set( txGRAD, 'String',num2str(dat.vid.grid.rad1) );
        guidata(hGUI,dat);
        set_grid(1);
        refresh_GUI();
    end %alt_grad_CB
    % set Ohman grid width
    function set_gwid_CB( s,e )
        dat	= guidata(hGUI);
        wid = str2double(get(txGWID,'String'));
        dat.vid.grid.wid = round(wid);
        guidata(hGUI,dat);
        set_grid(1);
        refresh_GUI();
    end %set_grad_CB
    % increment/decrement Ohman grid width
    function alt_gwid_CB( s,e, dr )
        dat	= guidata(hGUI);
        wid	= str2double(get(txGWID,'String'))+dr;
        dat.vid.grid.wid = round(wid);
        set( txGWID, 'String',num2str(dat.vid.grid.wid) );
        guidata(hGUI,dat);
        set_grid(1);
        refresh_GUI();
    end %alt_grad_CB
    % set Ohman grid primary origin location
    function set_gor1_CB( s,e )
        dat	= guidata(hGUI);
        cx	= str2double(get(txOR1X,'String'));
        cy	= str2double(get(txOR1Y,'String'));
        dat.vid.grid.or1	= round([cx cy]);
        dat.vid.grid.rad2	= norm( dat.vid.grid.or1-dat.vid.grid.or2 ) - dat.vid.grid.rad1;
        dat.vid.grid.theta	= atan( (cy-dat.vid.grid.or2(2)) / (cx-dat.vid.grid.or2(1)) );
        guidata(hGUI,dat);
        set_grid(1);
        refresh_GUI();
    end %set_gor1_CB
    % increment/decrement Ohman grid primary origin location
    function alt_gor1_CB( s,e, dx,dy )
        dat	= guidata(hGUI);
        cx	= round( str2double(get(txOR1X,'String')) + dx );
        cy	= round( str2double(get(txOR1Y,'String')) + dy );
        set( txOR1X, 'String',num2str(cx) );
        set( txOR1Y, 'String',num2str(cy) );
        dat.vid.grid.or1	= [cx cy];
        dat.vid.grid.rad2	= norm( dat.vid.grid.or1-dat.vid.grid.or2 ) - dat.vid.grid.rad1;
        dat.vid.grid.theta	= atan( (cy-dat.vid.grid.or2(2)) / (cx-dat.vid.grid.or2(1)) );
        guidata(hGUI,dat);
        set_grid(1);
        refresh_GUI();
    end %alt_gor1_CB
    % set Ohman grid secondary origin location
    function set_gor2_CB( s,e )
        dat	= guidata(hGUI);
        cx	= round( str2double(get(txOR2X,'String')) );
        cy	= round( str2double(get(txOR2Y,'String')) );
        dat.vid.grid.or2	= [cx cy];
        dat.vid.grid.rad2	= norm( dat.vid.grid.or1-dat.vid.grid.or2 ) - dat.vid.grid.rad1;
        dat.vid.grid.theta	= atan( (dat.vid.grid.or1(2)-cy) / (dat.vid.grid.or1(1)-cx) );
        guidata(hGUI,dat);
        set_grid(1);
        refresh_GUI();
    end %set_gor2_CB
    % increment/decrement Ohman grid secondary origin location
    function alt_gor2_CB( s,e, dx,dy )
        dat	= guidata(hGUI);
        cx	= round( str2double(get(txOR2X,'String')) + dx );
        cy	= round( str2double(get(txOR2Y,'String')) + dy );
        set( txOR2X, 'String',num2str(cx) );
        set( txOR2Y, 'String',num2str(cy) );
        dat.vid.grid.or2	= [cx cy];
        dat.vid.grid.rad2	= norm( dat.vid.grid.or1-dat.vid.grid.or2 ) - dat.vid.grid.rad1;
        dat.vid.grid.theta	= atan( (dat.vid.grid.or1(2)-cy) / (dat.vid.grid.or1(1)-cx) );
        guidata(hGUI,dat);
        set_grid(1);
        refresh_GUI();
    end %alt_gor2_CB
    % set Ohman gridline spacing
    function set_gspc_CB( s,e )
        dat	= guidata(hGUI);
        gi	= str2double(get(txGDSPC,'String'));
        dat.vid.grid.gint = round(gi);
        guidata(hGUI,dat);
        reset_vt_data();
        set_grid(1);
        refresh_GUI();
    end %set_gspc_CB
    % change Ohman grid linewidth
    function set_glinewd_CB( s,e )
        dat	= guidata(hGUI);
        lwd	= str2double(get(txGDLNW,'String'));
        dat.vid.grid.linewd = round(lwd);
        guidata(hGUI,dat);
        refresh_GUI();
    end %alt_gang_CB

    % set reference palate gridline limits
    function set_refphagl_CB( s,e )
        dat	= guidata(hGUI);
        gl1	= str2double(get(txPHA1,'String'));
        gl2	= str2double(get(txPHA2,'String'));
        dat.vid.grid.pha  = gl1:gl2;
        for gl = 1:dat.vid.grid.nlines
            dat.vid.grid.pharynx.pts(gl).rt	= [];
        end
        guidata(hGUI,dat);
        set( cbSHOWP, 'Value',1 );
        refresh_GUI();
    end
    % set reference palate gridline limits
    function set_refpalgl_CB( s,e )
        dat	= guidata(hGUI);
        gl1	= str2double(get(txPAL1,'String'));
        gl2	= str2double(get(txPAL2,'String'));
        dat.vid.grid.pal  = gl1:gl2;
        for gl = 1:dat.vid.grid.nlines
            dat.vid.grid.palate.pts(gl).rt	= [];
        end
        guidata(hGUI,dat);
        set( cbSHOWR, 'Value',1 );
        refresh_GUI();
    end
    % set lingual gridline limits
    function set_reflingl_CB( s,e )
        dat	= guidata(hGUI);
        gl1	= str2double(get(txTNG1,'String'));
        gl2	= str2double(get(txTNG2,'String'));
        dat.vid.grid.tng  = gl1:gl2;
        for gl = 1:dat.vid.grid.nlines
            dat.vt(dat.fnum).spts(gl).lf  = [];
        end
        for f = 1:dat.vid.info.NumFrames
            dat.vt(f).spts = [];
        end
        guidata(hGUI,dat);
        smooth_tongue_CB();
        %set( cbSHOWS, 'Value',1 );
        %refresh_GUI();
    end

    % update MRI FOV parameter
    function update_FOV_CB( s,e )
        dat	= guidata(hGUI);
        dat.status.FOV      = str2double(get(txSHFOV,'String'));
        dat.status.pixres	= wd_im * (dat.status.FOV/100) / dat.vid.info.Width;
        guidata(hGUI,dat);
        refresh_GUI();
    end %update_FOV_CB

	% toggle pixel intensity correction mode
    function tog_icorr_CB( s,e )
        dat = guidata(hGUI);
        if get( cbICORR,'Value' )
            if ( length(dat.cvd.avi) < length(dat.vid.avi)  )
                pixintensity_correct();
            end
        end
        refresh_GUI();
    end
	% update pixel intensity correction
    function pixintensity_correct( s,e )
        dat = guidata(hGUI);
        fprintf('   Correcting for uneveness in MRI video intensity\n');
        fprintf('   Assuming Coil 1 at [%d,%d], Coil 2 [%d,%d]\n', CoilL1(1),CoilL1(2), CoilL2(1),CoilL2(2) );
        dat.cvd.avi = mri_intensitycorrect( dat.vid.avi, CoilL1, CoilL2 );
        set( cbICORR, 'Value',1 );
        guidata(hGUI,dat);
        refresh_GUI();
    end

	% toggle video mode: tracked or normal
    function tog_vmode_CB( s,e )
        dat = guidata(hGUI);
        if get( cbTMODE,'Value' )
            if ( length(dat.tvd.avi) < 1  )
                fprintf('   Tracked video data not available\n');
                set( cbTMODE,'Value',0 )
            end
        end
        guidata(hGUI,dat);
        refresh_GUI();
    end

	% toggle filtered video mode
    function lpf_video_CB( s,e )
        dat = guidata(hGUI);
        if ( length(dat.lvd.avi) < length(dat.vid.avi)  )
            lpf_video();
        end
        refresh_GUI();
    end
	% filter video using discrete cosine transformation
    function lpf_video( s,e )
        dat = guidata(hGUI);
        fprintf('   Filtering video using discrete cosine transformation\n');
        dctthr = str2double(get(txVLPF,'String'));
        dat.lvd.avi = mri_dctlpf( dat.vid.avi, dctthr );
        guidata(hGUI,dat);
        refresh_GUI();
    end
    % set lpf/dct threshold parameter
    function set_dctthresh_CB( s,e )
        dat	= guidata(hGUI);
        dctthr	= str2double(get(txVLPF,'String'));
        dctdat	= dct2(dat.vid.avi(dat.fnum).cdata);
        dctmin	= min(min(dctdat));
        dctmax	= max(max(dctdat));
        if (dctthr > dctmin) && (dctthr < dctmax)
            lpf_video();
        else
            fprintf('\n   Specify a DCT LPF threshold in the range [%0f .. %0f]\n\n', dctmin,dctmax );
        end
    end %set_thresh_CB

    % set tract boundary detection threshold parameter
    function set_ithresh_CB( s,e )
        dat	= guidata(hGUI);
        ithrsh = str2double(get(txITHRS,'String'));
        if (ithrsh >= 0) && (ithrsh <= 1)
            dat.vid.grid.ithrsh = ithrsh;
            guidata(hGUI,dat);
            find_tissue_bnds( dat.fnum );
            refresh_GUI();
        else
            fprintf('\n   Specify a Tissue Boundary intensity threshold in the range [0..1]\n\n' );
        end
    end %set_thresh_CB

    % set tract centreline detection graph edge weighting parameter
    function set_diwtcl_CB( s,e )
        dat	= guidata(hGUI);
        diwtcl = str2double(get(txDIWTC,'String'));
        if (diwtcl >= 0) && (diwtcl <= 1)
            dat.vid.grid.diwtcl = diwtcl;
            guidata(hGUI,dat);
            find_tissue_bnds( dat.fnum );
            refresh_GUI();
        else
            fprintf('\n   Specify a Distance-Intensity Weighting factor in the range [0..1]\n\n' );
        end
    end %set_diwtcl_CB

    % set inner tract boundary graph edge weighting parameter
    function set_gamma_CB( s,e )
        dat	= guidata(hGUI);
        diwtib = str2double(get(txGAMMA,'String'));
        if (diwtib >= 0) && (diwtib < 1)
            dat.vid.grid.diwtib = diwtib;
            guidata(hGUI,dat);
            find_tissue_bnds( dat.fnum );
            refresh_GUI();
        else
            fprintf('\n   Specify a Distance-Intensity Weighting factor 0 <= wt < 1\n\n' );
        end
    end %set_gamma_CB

    % set boundary smoothing dct threshold parameter
    function set_bndsmthresh_CB( s,e )
        dat	= guidata(hGUI);
        dctth1 = 0;	dctth2 = 100;
        dctthr = str2double(get(txDCTTH,'String'));
        if ( (dctthr >= 0) && (dctthr <= 100) )
            if ~isempty(dat.vid.grid.tng)
                glines	= dat.vid.grid.tng;
                dctthr	= str2double(get(txDCTTH,'String'));
                smooth_ibnd( dat.fnum, dctthr, glines, 1 );
            end
            set( cbSHOWS, 'Value',1 );
            refresh_GUI();
        else
            fprintf('\n   Specify a DCT LPF threshold in the range [%0.1f .. %0.1f]\n\n', dctth1,dctth2 );
        end
    end %set_thresh_CB

    % set pixel intensity analysis radius
    function set_intanalrad_CB( s,e )
        dat	= guidata(hGUI);
        irad = str2double(get(txANRAD,'String'));
        rmax = floor(max(dat.vid.dims)/2);
        if (irad > 0) && (irad < rmax)
            dat.status.radanal = irad;
            guidata(hGUI,dat);
            refresh_GUI();
        else
            fprintf('\n   Specify an analysis radius in the range [1..%d]\n\n', rmax );
        end
    end %set_thresh_CB


    % intialize correlation analysis
    function init_correl_CB( s,e, cdur )
        dat = guidata(hGUI);
        % ensure images displayed as uninterpolated
        set( cbVMODE, 'Value',1 );
        guidata(hGUI,dat);
        refresh_GUI();
        % prep analysis frame interval
        fa	= dat.vid.fint(1);
        fb	= dat.vid.fint(2);
        nf	= fb-fa;
        
        vobj = VideoReader(dat.ffn_vid);
        hh	 = vobj.Height;
        ww	 = vobj.Width;
        avi(1:nf) = struct('cdata', zeros(hh,ww,3,'uint8'), 'colormap',[] );
        for i = 1:nf
            avi(i).cdata = read(vobj,ff(1)+i);
        end
        
        h	= dat.vid.info.Height;
        w	= dat.vid.info.Width;
        hh	= 68;	ww = 68;
        offsetH = 0;
        if hh<h
            offsetH = floor((h-hh)/2);
        end
        offsetW = 0;
        if ww<w
            offsetW = floor((w-ww)/2);
        end
        FRAMES = zeros(nf,hh*ww);
        for fr = 1:nf
            im = double( avi(fr).cdata(offsetH+1:h-offsetH,offsetW+1:w-offsetW,1) );
            im = reshape( im,1,hh*ww );
            FRAMES(fr,:) = im;
        end
        nT	= str2double(get(txCORTR,'String'));
        dat.corr.nT	= nT;                   % number of analysis traces
        dat.corr.px(nT).coord = [];         % coordinates of selected pixel
        dat.corr.tl = dat.aud.tint(1);      % left edge of analysis window (sec)
        dat.status.ptpicked = 0;
        dat.status.showmloc = 1;
        dat.status.cduranal = cdur;
        set( txSHPOX, 'Visible','On' );
        set( txSHPOY, 'Visible','On' );
        set(hGUI,'CurrentAxes',fVID);       % plot mean image
        guidata(hGUI,dat);
        im	= reshape( mean(FRAMES,1), hh,ww );
        [ im_, map_ ] = imresize( im, cm_corr, scale );
        subimage( im_,map_ ); axis off;
        gpos = get( hGUI, 'Position' );
        set(0, 'pointerlocation', [gpos(1)+0.24*gui_wd gpos(2)+0.65*gui_ht]);
        set(hGUI, 'Pointer','crosshair');
        guidata(hGUI,dat);
    end %init_correl_CB

    % correlation analysis: fetch selected pt
    function fetch_correl(pt)
        dat = guidata(hGUI);
        npt	= dat.status.ptpicked;
        pts	= dat.corr.nT;
        dat.corr.px(npt).coord = pt;
        guidata(hGUI,dat);
        fprintf('   Selected analysis point %d of %d: x=%d y=%d\n', npt,pts, pt(1),pt(2) );
        if (npt == pts)
            dat.status.showmloc = 0;
            set( txSHPOX, 'Visible','Off' );
            set( txSHPOY, 'Visible','Off' );
            set(hGUI, 'Pointer','arrow');
            guidata(hGUI,dat);
            analyze_correl();
        end
    end %fetch_correl

    % correlation analysis
    function analyze_correl()
        dat = guidata(hGUI);
        fa	= dat.vid.fint(1);
        fb	= dat.vid.fint(2);
        nf	= fb-fa;
        
        vobj = VideoReader(dat.ffn_vid);
        h	 = vobj.Height;
        w	 = vobj.Width;
        hh	= 68;	ww = 68;
        avi(1:nf) = struct('cdata', zeros(hh,ww,3,'uint8'), 'colormap',[] );
        for i = 1:nf
            avi(i).cdata = read(vobj,fa+i);
        end
        
        offsetH = 0;
        if hh<h
            offsetH = floor((h-hh)/2);
        end
        offsetW = 0;
        if ww<w
            offsetW = floor((w-ww)/2);
        end
        FRAMES = zeros(nf,hh*ww);
        for fr = 1:nf
            im = double( avi(fr).cdata(offsetH+1:h-offsetH,offsetW+1:w-offsetW,1) );
            im = reshape( im,1,hh*ww );
            FRAMES(fr,:) = im;
        end
        C       = corrcoef(FRAMES);
        POINTS	= [];
        nT      = dat.corr.nT;
        for tr = 1:nT
            pt	= dat.corr.px(tr).coord;
            p	= round((ww*(round(pt(1))-1)) + round(pt(2)));
            POINTS = [POINTS; p];
        end
        % image regions
        hCG = figure; hold on;
        fig_pos     = get(hGUI,'Position');
        fig_pos(2)	= fig_pos(2)+gui_ht+25;
        fig_pos(3)	= 0.8*gui_wd;           % correlation figure width
        if ( (nT>1) && get( cbCORTR,'Value' ) )
            fig_pos(4) = nT*(0.215*gui_ht);	% multi-row figure height
        else
            fig_pos(4) = 0.275*gui_ht;      % single row figure height
        end;
        set( hCG, 'Name','Correlation Analysis', 'Position', fig_pos );
        set( hCG, 'ToolBar','none'  );
        %set( hCG, 'MenuBar','none'  );
        set( hCG, 'WindowButtonDownFcn',{@corr_mouse_down_CB} );
        for tr = 1:nT
            p1R = C(POINTS(tr),:);
            if ( (nT>1) && get( cbCORTR,'Value' ) )
                [l b w h] = calc_subplot_dims( nT,7, tr,1, 0.02,0.05, 0.05,0.05 );
            else
                [l b w h] = calc_subplot_dims( 1,7,  1,tr, 0.02,0.2,  0.02,0 );
            end;
            subplot('position',[l b w h]);
            im = reshape(p1R,hh,ww);
            imagesc( im );
            axis('square'); axis('tight'); axis('off');
            pt	= dat.corr.px(tr).coord;
            xx	= pt(1); yy = pt(2);
            txt	= sprintf( 'R%d (x=%d y=%d)\n', tr,xx,yy );
            hTitle = title( txt );
            posTtl = get(hTitle,'Position') + [0 15 0];
            ctxt = colcyc{tr};
            set( hTitle, 'Color',ctxt, 'Position',posTtl );
            line( [min(xlim) xx-dpx],[yy yy],[0 0], 'Color','k' );
            line( [xx+dpx max(xlim)],[yy yy],[0 0], 'Color','k' );
            line( [xx xx],[min(ylim) yy-dpx],[0 0], 'Color','k' );
            line( [xx xx],[yy+dpx max(ylim)],[0 0], 'Color','k' );
        end
         % plot traces
        for tr = 1:nT
            if ( (nT>1) && get( cbCORTR,'Value' ) )
                switch(nT)
                   case 2
                      dy = 0.1;
                   case 3
                      dy = 0.075;
                   case 4
                      dy = 0.05;
                end
                [l b w h] = calc_subplot_dims( nT,7, tr,2:7,     0.02,dy,  0.05,dy );
            else
                [l b w h] = calc_subplot_dims( 1,7,  1,(nT+1):7, 0.05,0.2,   0.05,0   );
            end;
            hPLT = subplot('position',[l b w h]); hold on;
            set( hPLT,'UserData','CATrace' );
            p1R = C(POINTS(tr),:);
            thr = str2double(get(txCORTH,'String'));
            pp	= ( p1R >= thr );
            pLT	= mean(double(FRAMES(:,pp)),2);
            %pLT	= (pLT-min(pLT))./range(pLT);
            pLT	= (pLT-min(pLT));
            cpl	= colcyc{tr};
            plot( pLT, 'color',cpl, 'linewidth',1 ); hold on;
            %plot( pLT, 'b-o', 'linewidth',2 ); hold on;
            axis('tight');
            dat.corr.tr(tr).pp	= pLT;
            xmax = size(FRAMES,1);
            xlim([1 xmax]); %ylim([0 1]);
            if (tr==nT),
                set( hPLT, 'XTick',[1 median(xlim) xmax] );
                dt	= diff(dat.aud.tint);
                l1	= sprintf('%d',0);
                l2	= sprintf('%0.2f',dt/2);
                l3	= sprintf('%0.2f',dt);
                set( hPLT, 'XTickLabel',{l1,l2,l3} );
                %xlabel( 'Duration (sec)' );
            else
                set( hPLT, 'XTick',[] );
            end;
            txt1 = sprintf( 'Corr > %0.2f: ', thr );
            txt2 = sprintf( '%0.2f to %0.2f sec', dat.aud.tint(1), dat.aud.tint(2) );
            txt3 = sprintf( ' (Frames %d to %d)', dat.vid.fint(1), dat.vid.fint(2) );
            if (tr==1),
                hTitle = title( [txt1 txt2 txt3] );
                posTtl = get(hTitle,'Position') + [0 -0.05 0];
                set( hTitle, 'Position',posTtl );
            end;
            % calculate duration of correlated pixel intensity peak if flagged
            if (dat.status.cduranal)
                Ith = 0.25;
                xx	= find( pLT >= Ith );
                t1	= xx(1); t2	= xx(end); 
                dur	= (t2-t1)/dat.vid.info.FramesPerSecond;
                txt	= sprintf( 'D = %0.2f sec', dur );
                text( xmax-2.5,0.875,txt, 'FontSize',9 );
                plot( t1,pLT(t1),'ko', 'MarkerSize',10 );
                plot( t2,pLT(t2),'ko', 'MarkerSize',10 );
                % calculate area under plot
                xx	= linspace(0,dt,xmax);
                Ac	= trapz(xx,pLT);
                txt	= sprintf( 'A = %0.2f', Ac );
                text( xmax-2.5,0.775,txt, 'FontSize',9 );
                % store analysis metrics
                seg = get(txLABEL,'String');
                if isfield(dat.seg, seg),
                    ff	= dat.seg.(seg).fint;
                    fprintf('\n   Writing analysis data for segment ''%s'':\n', seg );
                    dat.seg.(seg).Icurve	= pLT;
                    dat.seg.(seg).tint      = [dat.aud.tint(1) dat.aud.tint(2)];
                    dat.seg.(seg).threshint	= [t1 t2];
                    dat.seg.(seg).cen_pixel	= pt;
                    dat.seg.(seg).Cth       = thr;
                    dat.seg.(seg).Ith       = Ith;
                    dat.seg.(seg).duration	= dur;
                    dat.seg.(seg).peak_area	= Ac;
                    fprintf('     Frame interval:         [%0d %0d]\n', ff(1),ff(2) );
                    fprintf('     Time interval:          [%0.2f %0.2f]\n', dat.aud.tint(1),dat.aud.tint(2) );
                    fprintf('     Analysis center pixel:  [%d %d]\n',   pt(1),pt(2)  );
                    fprintf('     Correlation threshold:  %0.2f\n',     thr );
                    fprintf('     Intensity threshold:    %0.2f\n',     Ith );
                    fprintf('     Duration:               %0.2f sec\n', dur );
                    fprintf('     Area under curve:       %0.2f\n\n',   Ac  );
                end;
            end
            %fps	= dat.vid.info.FramesPerSecond;
            %t1	= ceil(dat.aud.tint(1));
            %t2	= floor(dat.aud.tint(2));
            %set( hPLT, 'XTick',((t1-dat.aud.tint(1))*fps):fps:t2*fps );
            %set( hPLT, 'XTickLabel',t1:t2 );
        end
        guidata(hGUI,dat);
        refresh_GUI();
    end %analyze_correl

    % handle mouse activity on correlation analysis figure
    function corr_mouse_down_CB( s,e )
        dat	= guidata(hGUI);
        set( s,'WindowButtonUpFcn',@corr_mouse_up_CB);
        hAX	= gca();
        pt	= get( hAX,'CurrentPoint' );
        t1	= pt(1,1)/dat.vid.info.FramesPerSecond + dat.corr.tl;
        function corr_mouse_up_CB( s,e )
            pt2	= get( hAX,'CurrentPoint' );
            t2	= pt2(1,1)/dat.vid.info.FramesPerSecond + dat.corr.tl;
            if (t2 == t1)	% user has selected single point in time
                fprintf('\n   Selected time: %0.2f\n',t1);
                dat.aud.wwd = 0.5;	% zoom into half second interval for inspection
                set(txWINWD,'String',num2str(dat.aud.wwd,'%.1f'));
                guidata(hGUI,dat);
                update_time( t1 );
                update_correl_win( hAX );
            else            % user has dragged mouse over interval to measure duration
                if (t2 < t1)
                    tt = t1; t1 = t2; t2 = tt;
                end
                tc = mean([t1 t2]);
                tw = diff([t1 t2]);
                dat.aud.wwd = tw;
                set(txWINWD,'String',num2str(tw,'%.1f'));
                guidata(hGUI,dat);
                fprintf( '\n   Start time:  %3.2f',   t1 );
                fprintf( '\n   End time:    %3.2f',   t2 );
                fprintf( '\n   Duration:    %3.2f\n', tw );
                update_time( tc );
                update_correl_win( hAX );
            end %if
        end %corr_mouse_up_CB
    end %corr_mouse_down_CB

    % handle mouse activity on correlation analysis figure
    function  update_correl_win( hAX )
        xx	= (dat.aud.tint - dat.corr.tl)*dat.vid.info.FramesPerSecond;
        yy	= get( hAX,'YLim' );
        % clear any objects indicating previous windowed region on subplot
        hl	= findobj( 'UserData','WinCen' ); delete(hl);
        hr	= findobj( 'UserData','WinLim' ); delete(hr);
        % indicate new windowed region on subplot
        %hl	= line( [mean(xx) mean(xx)],yy,[0 0], 'Color','k','LineStyle','--', 'Parent',hAX );
        %hr	= rectangle( 'Position',[min(xx),min(yy),diff(xx),diff(yy)], 'FaceColor',colwin, ...
        %                 'EdgeColor',colwlim, 'EraseMode','xor', 'Parent',hAX );
        set( hl,'UserData','WinCen' );
        set( hr,'UserData','WinLim' );
    end %update_correl_win

    % correlation frequency analysis
    function analyze_cofreq_CB( s,e )
        dat = guidata(hGUI);
        if isfield( dat,'corr')
            Fs	= dat.vid.info.FramesPerSecond;
            NFT = 2048;
            hFG = figure; hold on;
            set( hFG, 'toolbar','none' );
            title('Correlation Analysis Trace Spectra');
            xlabel('Frequency (Hz)');	ylabel('|Y(f)|');
            for tr = 1:dat.corr.nT;
                Y	= fft(dat.corr.tr(tr).pp,NFT);
                f	= Fs/2*linspace(0,1,NFT/2+1);
                plot( f,2*abs(Y(1:NFT/2+1)), 'color',colcyc{tr} );
                %plot( f,2*abs(Y(1:NFT/2+1)), 'color','b' );
                xlim([0 6]);
            end
        else
            fprintf('\n   No correlation analysis data found ... \n');
        end
    end %analyze_cofreq_CB

    % mean local pixel intensity analysis
    function intensity_anal( s,e )
        dat = guidata(hGUI);
        lab = get(txLABEL,'String');
        if ~(isempty(dat.seg))
            if isfield(dat.seg, lab)
                ff	= dat.seg.(lab).fint(1):dat.seg.(lab).fint(2);
                pt	= dat.seg.(lab).cen_pixel;
                r	= dat.status.radanal;
                fprintf('\n   Segment [%s]: frames %d to %d\n', lab,ff(1),ff(end) );
                fprintf(  '   Center pixel:    [%d %d]\n',      pt(1),pt(2) );
                fprintf(  '   Analysis radius = %d\n\n',        r );
                hIA = figure; hold on;
                set( hIA, 'toolbar','none' );
                title(sprintf('Intensity in Vicinity: ''%s'' (x=%d; y=%d; r=%d)', lab,pt(1),pt(2),dat.status.radanal));
                xlabel('Frame');  ylabel('Mean Pixel Intensity');
                Ipp	= zeros(1,length(ff));
                for f = ff
                    Ipp(f)	= calc_local_px_intensity( f,pt,r );
                end
                plot( ff,Ipp(ff) );
                dat.seg.(lab).Ipp = Ipp;
                guidata(hGUI,dat);
            end
        end
    end %intensity_anal

    % calculate mean local pixel intensity in circle of radius r 
    % around point pt on frame f
    function It = calc_local_px_intensity( f,pt,r )
        It	= 0;
        Ii	= 0;
        for xx = pt(1)-r:pt(1)+r
            for yy = pt(2)-r:pt(2)+r
                if (norm([xx yy]-[pt(1) pt(2)]) < r)
                    %Ip	= double(fr(xx,yy));
                    Ip	= uint32(get_pix_pt( f,[xx yy] ));
                    It	= uint32(It + Ip);
                    Ii	= Ii + 1;
                end
            end
        end
        It	= uint32(round(It./Ii));
    end %calc_local_px_intensity

    function [pt,Ip] = get_pixel_pt( s,e )
        dat = guidata(hGUI);
        pt	= get(gca,'CurrentPoint');
        pt	= pt(1,1:2);
        if ( (get(cbSHVAX,'Value')&& ~(get(cbSCVAX,'Value'))&& ~(get(cbVMODE,'Value'))) || dat.status.showmloc )
            pt	= pt/dat.vid.scale;
        end
        im	 = dat.vid.avi(dat.fnum).cdata;
        if ( (get(cbSHVAX,'Value')&&~(get(cbSCVAX,'Value'))) || get(cbVMODE,'Value') )
            im	 = dat.vid.avi_(dat.fnum).cdata;
        end
        imsz = size(im);
        imsz = imsz(1:2);
        if ( get(cbSHVAX,'Value') && get(cbSHVMM,'Value') )	% fetch coordinates in mm
            if ( (pt>0) & (pt<imsz) )
                Ip	= im( round(pt(2)),round(pt(1)) );
                pt	= pt * dat.status.pixres;
            end
        else
            if ( (pt>0) & (pt<imsz) )
                pt	= ceil(pt);
                Ip	= im( pt(2),pt(1) );
            end
        end
        if ( (pt>0) & (pt<imsz) )
            dat.status.mloc	= pt;
            posx = sprintf('%d',pt(1));
            posy = sprintf('%d',pt(2));
            set( txSHPOX, 'String',posx );
            set( txSHPOY, 'String',posy );
            guidata(hGUI,dat);
        end
    end %get_pixel_pt

    function Ip = get_pix_pt( f,pt )
        dat = guidata(hGUI);
        if ( (get(cbSHVAX,'Value')&& ~(get(cbSCVAX,'Value'))&& ~(get(cbVMODE,'Value'))) || dat.status.showmloc )
            pt	= pt/dat.vid.scale;
        end
        im	= dat.vid.avi(f).cdata;
        if ( (get(cbSHVAX,'Value')&&~(get(cbSCVAX,'Value'))) || get(cbVMODE,'Value') )
            im	= dat.vid.avi_(f).cdata;
        end
        imsz = size(dat.vid.avi_(f).cdata);
        if ( get(cbSHVAX,'Value') && get(cbSHVMM,'Value') )	% fetch coordinates in mm
            if ( (pt>0) & (pt<imsz) )
                Ip	= im( round(pt(2)),round(pt(1)) );
                pt	= pt * dat.status.pixres;
            end
        else
            if ( (pt>0) & (pt<imsz) )
                pt	= ceil(pt);
                Ip	= im( pt(2),pt(1) );
            end
        end
    end %get_pix_pt

    % video reading wrapper function to hide platform and matlab-version
    % specific details 
    function vid = read_dmri_video( file_name )
        vid.info = aviinfo( file_name );
        comp_type = computer();
        if ~strncmp(comp_type,'MAC',3)
            vobj = VideoReader(dat.ffn_vid);
            hh	 = vobj.Height;
            ww	 = vobj.Width;
            nf	 = vobj.NumberOfFrames;
            avi(1:nf) = struct('cdata', zeros(hh,ww,3,'uint8'), 'colormap',[] );
            for fr = 1:nf
                avi(fr).cdata = read(vobj,fr);
            end
        else
            avi = aviread(file_name); 
        end
        vid.avi	 = avi;
        if ( isempty(vid.info) )
            error( '   Can''t load video from file <%s>\n', file_name );
        end
    end


    % reset vocal tract data structures
    function reset_vt_data( s,e )
        dat	= guidata(hGUI);
        nf = dat.vid.info.NumFrames;
        for f = 1:nf
            dat.vt(f).imap      = [];	% intensity map
            dat.vt(f).Mxy       = [];   % coodinates of minima nodes in graph
            dat.vt(f).Mwt       = [];   % weightings of minima nodes in graph
            dat.vt(f).cen       = [];   % tract center coordinates
            dat.vt(f).len       = [];   % tract length
            dat.vt(f).pts       = [];   % coordinates of tissue boundaries
            dat.vt(f).spts      = [];   % coordinates of smoothed tissue boundaries
            dat.vt(f).spts.lf	= [];
            dat.vt(f).man       = [];
            dat.vt(f).aut       = [];
            dat.vt(f).iix       = [];   % intensity indices
            dat.vt(f).int       = [];   % intensity of lf & rt tissue boundary thresholds
            dat.vt(f).AF        = [];	% tract area function
            dat.vt(f).AFx       = [];   % x-index of tract area function
            dat.vt(f).ngl       = [];	% number of gridlines used in this frame
            dat.vt(f).dx        = [];	% horizontal image displacement from reference frame
            dat.vt(f).dy        = [];	% vertical image displacement from reference frame
        end
        guidata(hGUI,dat);
        refresh_GUI();
    end %reset_vt_data

    % reset reference boundary data structures
    function reset_refbnd_data( s,e )
        dat	= guidata(hGUI);
        dat.vid.grid.pal        = [];
        dat.vid.grid.palate     = [];
        dat.vid.grid.palfnum	= 0;
        dat.vid.grid.pha        = [];
        dat.vid.grid.pharynx	= [];
        dat.vid.grid.phafnum	= 0;
        dat.vid.grid.tng        = [];
        guidata(hGUI,dat);
        refresh_GUI();
    end %reset_refbnd_data

    % exit functions
    function save_GUI_CB( s,e )
        dat = guidata(hGUI);
        dat.status.saved = 1;
        disp('   Saving all data structures.' );
        set( pbSAVE,'String','Saved' );
        guidata(hGUI,dat);
    end %save_GUI_CB

    function close_GUI_CB( s,e )
        dat = guidata(hGUI);
        if (~dat.status.saved)	% clean up big variables to save memory
            disp('   Saving metadata structures only.' );
            dat.vid.avi  = [];
            dat.vid.avi_ = [];
            dat.tvd.avi	 = [];
            dat.cvd.avi	 = [];
            dat.lvd.avi	 = [];
            set( pbSAVE,'String','Saved' );
            dat.status.exit = 1;
        end
        set( pbEXIT,'String','Exiting' );
        guidata(hGUI,dat);
    end %close_GUI_CB

    function image_hover_CB( s,e )
        dat = guidata(hGUI);
        if ( dat.status.showmloc )	% if mouse over pixel correlation image
            get_pixel_pt( s,e );
        end
    end %image_hover_CB


end %of main function
