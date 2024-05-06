function dat = inspect_rtmri( x, grid )
%
%  INSPECT rtMRI: GUI for synchronizing, inspecting and analyzing paired WAV
%  & AVI files containing dynamic MRI data acquired at USC.
%  (see Narayanan et al. 2004 for details of acquisition and protocols).
%
%  Assumes existence of synchronised audio and video files which share the 
%  base filename 'token', in dedicated subdirectories of working directory.
%
%  dat = inspect_rtmri( token, (grid) )
%
%    eg. inspect_rtmri( 'paLam' )
%
%        - load audio file <paLam.wav> from subdir <./wav/>
%        - check for video file <paLam.avi> in subdir <./avi/>
%        - check for tracked video file <paLam.avi> in subdir <./avi_trak/>
%        - display synchronized audio & spectral & video data for analysis
%
%    eg. dat = inspect_rtmri( 'paLam' );         % load data from file
%    eg. dat = inspect_rtmri( tamil1 );          % reload data structure <tamil1>
%    eg. dat = inspect_rtmri( 'paLam', gd );     % (use analysis grid 'gd')
%


% video/OS configuration
invert_video	= 0;        % set flag to invert video on some OS/codec combinations
%if isunix
%    invert_video	= 1;	% set flag to invert video on some OS/codec combinations
%end
avoid_wmark     = 1;        % set flag to configure specific frame crop parameters

% audio filter configuration
truncate_burst  = 1;        % set flag to truncate initial noise burst
len_burst       = 0.55;     % length of initial noise burst (secs)
atten_burst     = 0.1;      % attenuation factor to apply to initial noise burst

% detect Matlab version to configure Signal Processing Toolbox toolset
verMatlab	=  version('-release');
relMatlab	=  str2num(verMatlab(1:4));
fprintf( '\n   Using Matlab version %s\n', verMatlab );

% declare constants
mypi	= 3.141593;
deg2rad	= mypi/180;
rad2deg	= 1/deg2rad;
FOV     = 100;              % MRI field of view (%)
wd_im	= 200;              % width (mm) of MR Image at 100% FOV
wd_vid	= 340;              % width (px) of video display window
wbi     = 5;                % waitbar interval
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
set( hGUI, 'menubar','none', 'toolbar','none', 'Resize','off', 'WindowButtonDownFcn',{@mouse_down_CB} );
dat     = guihandles(hGUI);
guidata(hGUI,dat);
dat.tok	= tok;

% configure GUI color scheme
bgc	= get( hGUI, 'Color' );
set( hGUI, 'DefaultUIControlBackgroundColor',bgc );	% default control color to that of background
bgd = get(0, 'FactoryUIControlBackgroundColor');    % interactive controls default to darker bgc

% construct filenames/paths
home	= pwd();        % base directory containing audio and video subdirs
dir_aud	= 'wav';        % directory containing audio
dir_vid	= 'avi';        % directory containing video
dir_trk	= 'avi_trak';	% directory containing tracked video
dir_cpt	= 'capture';	% directory in which to write captured frames
dir_out	= 'jpg';        % jpg output directory
dir_mov	= 'avi';        % video output directory
dir_wav	= 'wav';        % audio output directory
dir_jpg	= 'jpg';        % image output directory
fn_aud	= [ tok '.wav' ];
fn_vid	= [ tok '.avi' ];
dat.ffn_aud	= fullfile( home, dir_aud, fn_aud );
dat.ffn_vid	= fullfile( home, dir_vid, fn_vid );
dat.ffn_trk	= fullfile( home, dir_trk, fn_vid );

% inialize global configuration flags
dat.status.reload	= 0;
dat.status.FOV      = FOV;
dat.status.saved	= 0;
dat.status.exit     = 0;
dat.status.radanal	= 5;
if isstruct(old_dat)
    dat.status.reload	= 1;
    cropframe	= old_dat.vid.frame;
end
% detect OS
if	strcmp(getenv('OS'),'Windows_NT')
	dat.status.OS = 'win';
else
	dat.status.OS = 'lin';
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

dat.vid.info	= [];   dat.tvd.info	= [];
dat.vid.avi     = [];   dat.tvd.avi     = [];
dat.cvd.avi     = [];   dat.lvd.avi     = [];
dat.vid.avi_	= [];
dat.vid.fint	= [0 0];
dat.vid.scale	= 1;
dat.vid.dims	= 0;
dat.vid.invert	= invert_video;


% create main GUI
ssz	= get(0,'ScreenSize');
gui_wd	= 975;
gui_ht	= 535;
gui_lf	= (ssz(3)-gui_wd)/2;
gui_bt  = (ssz(4)-gui_ht)/2;
gui_pos	= [gui_lf gui_bt gui_wd gui_ht];
set( hGUI, 'Position', gui_pos, 'toolbar','none' ); axis off;
set( hGUI, 'Resize','off',  'WindowButtonDownFcn',{@mouse_down_CB} );

% GUI font sizes
fsx = 11;  fsh = 10;  fsl = 9;  fss = 8;  fst = 7;  fsm = 6;

% calculate dimension of GUI panels
hoff = 0.008;	voff = 0.008;   % inter-panel spacing
hvid = 1-2*voff;                % panel heights 
hana = 0.33;
haud = hvid-hana-voff;
hclo = 0.275*hana;
hexp = hana-hclo-voff;
hlab = 0.275*hana;
hcor = hana-hlab-voff;
baud = hana+2*voff;             % locate bottom edge of panels
bexp = hclo+2*voff;
bcor = hlab+2*voff;
wvid = 0.44;                    % panel widths
wcor = 0.08;
wexp = wcor;
laud = wvid+2*hoff;             % locate left edge of panels
waud = 1-wvid-3*hoff;
wana = waud-wexp-wcor-2*hoff;
lcor = laud+wana+hoff;
lexp = lcor+wcor+hoff;

% create GUI panels
pVID = uipanel(	'Parent',hGUI, 'BackgroundColor',bgc, 'Title','Video',	...
                'FontSize',fsx, 'Position',[hoff voff wvid hvid],       ...
                'ButtonDownFcn',{@play_vid_CB}  );
pAUD = uipanel(	'Parent',hGUI, 'BackgroundColor',bgc, 'Title','Audio',	...
                'FontSize',fsx, 'Position',[laud baud waud haud],       ...
                'ButtonDownFcn',{@play_aud_CB}  );
pANA = uipanel(	'Parent',hGUI, 'BackgroundColor',bgc, 'Title','Analyze Tract',	...
                'FontSize',fsh, 'Position',[laud voff wana hana]  );
pCOR = uipanel(	'Parent',hGUI, 'BackgroundColor',bgc, 'Title','CorAnal', ...
                'FontSize',fsh, 'Position',[lcor bcor wcor hcor]  );
pLAB = uipanel(	'Parent',hGUI, 'BackgroundColor',bgc,	...
                'FontSize',fsh, 'Position',[lcor voff wcor hlab]  );
pEXP = uipanel(	'Parent',hGUI, 'BackgroundColor',bgc, 'Title','Export',	...
                'FontSize',fsh, 'Position',[lexp bexp wexp hexp]  );
pCLO = uipanel(	'Parent',hGUI, 'BackgroundColor',bgc,	...
                'FontSize',fsh, 'Position',[lexp voff wexp hclo]  );

% calculate gui panel dimensions in pixels 
hpvid = hvid * gui_ht;	% panel heights
hpaud = haud * gui_ht;
hpana = hana * gui_ht;
hpclo = hclo * gui_ht;
wpvid = wvid * gui_wd;	% panel widths
wpaud = waud * gui_wd;
wpexp = wexp * gui_wd;
wpana =	wana * gui_wd;

% gui control dimensions and placement coordinates 
dx1	= 1; dy1 = 2;       % inter-element displacements 
dx2	= 4; dy2 = 3;
dxo	= 5; dyo = 4;       % displacement from panel edges

% width & height of standard analysis panel controls:
% arrange (pb) elements on a grid of 6 rows x 4 cols
wp1 = (wpana-2*dxo)/4 - dx2;    % std pb width
hp1 = (hpana-2*dyo)/6 - 5;      % std pb height
hp2 = (hpana-2*dyo)/7 - 2;      % thinner pb height

hcb = 18;               % checkbox height
hl1 = hcb-dy1;          % label height

wcb = hcb;              % checkbox width
wt1 = 0.400*wp1;        % textbox width
wlc = wp1-wcb+dx1;      % label widths
wlt = wp1-wt1+dx1;

le1 = dxo;              % left of 1st element in frame
lp2 = le1+wp1+dx2;      % left of nth pb-width element in frame
lp3 = lp2+wp1+dx2;
lp4 = lp3+wp1+dx2;
lc1 = le1+wlc;          % left of nth std checkbox
lc2 = lp2+wlc;
lc3 = lp3+wlc;
lc4 = lp4+wlc;
lt1 = le1+wlt;          % left of nth std textbox
lt2 = lp2+wlt;
lt3 = lp3+wlt;
lt4 = lp4+wlt;

bo1 = dyo;              % bottom of lowest element in frame
bp2 = bo1+hp1+dy2;      % bottom of nth pb-height element
bp3 = bp2+hp1+dy2;
bp4 = bp3+hp1+dy2;
bp5 = bp4+hp1+dy2;
bp6 = bp5+hp1+dy2;
bl2 = bo1+hl1+1;      % bottom of nth lbl-height element
bl3 = bl2+hl1+1;
bl4 = bl3+hl1+1;
bl5 = bl4+hl1+1;
bl6 = bl5+hl1+1;
bl7 = bl6+hl1+1;

                
                
% read audio, get length & sampling rate
if ( exist(dat.ffn_aud,'file') )
    fprintf( '   Using audio file:     <%s>\n', dat.ffn_aud );
    if exist('audioread')
        [y, Fs]	= audioread( dat.ffn_aud );
    else
        [y, Fs]	= wavread( dat.ffn_aud );
    end
    if (truncate_burst)
        y(1:Fs*len_burst) = y(1:Fs*len_burst).*atten_burst;
    end
    dat.aud.sig	= y;
    dat.aud.Fs	= Fs;
    dat.aud.len	= length(dat.aud.sig);
    dat.aud.dur	= dat.aud.len/dat.aud.Fs;
else
    fprintf( '   Can''t find audio file <%s.wav> in subdirectory <%s> ...\n', tok, dir_aud );
    action	= input( '   Continue without audio [y/n]?  ', 's') ;
    if strcmpi( action,'n' ) || strcmpi( action,'no' )
        fprintf( '   \n\n' );
        close(hGUI);
        return;
    else	% use 4 sec dummy audio signal until length of video is known ...
        dat.ffn_aud = '';
        dat.aud.Fs	= 20000;
        dat.aud.dur	= 4;
        dat.aud.len	= round(dat.aud.dur * dat.aud.Fs);
        dat.aud.sig	= zeros(1,dat.aud.len);
    end
end

% specify configuration of audio and sg subplots in audio frame
lwav = 0.10;    wwav = 0.79;    % location of waveform plot
bwav = 0.82;    hwav = 0.17;
bwin = 0.59;                    % location of waveform window plot
bsgm = 0.11;	hsgm = 0.41;    % location of spectrogram plot

% display complete audio signal
fAUD = subplot( 4,1,1, 'Parent',pAUD );
dat.aud.tvec = (1:dat.aud.len)./dat.aud.Fs;
plot( dat.aud.tvec, dat.aud.sig );
xlim([0 dat.aud.dur]);
set( fAUD, 'XMinorTick','on', 'FontSize',fss );
set( fAUD, 'Position',[lwav bwav wwav hwav] )

% display window of audio signal
fAUDw = subplot( 4,1,2, 'Parent',pAUD );
plot( dat.aud.tvec, dat.aud.sig );
dat.aud.wwd	 = sig_wwd;                             % audio window width
dat.aud.tint = [0 sig_wwd];                         % audio interval (seconds)
dat.aud.sint = floor([0 sig_wwd].*dat.aud.Fs)+1;    % audio interval (samples)
xlim([0 sig_wwd]);
set( fAUDw, 'Color',colwin, 'FontSize',fss );
set( fAUDw, 'Position',[lwav bwin wwav hwav] )
%xlabel('Time (sec)'); ylabel('Amplitude');

% audio GUI controls: adjust center and width of waveform displays
dx = 8;	 dy = 2;
lac = 0.905 * wpaud;        % left edge of all audio & spectrogram uicontrols
wac = wpaud-lac-dx;
wat = 0.65*wac;	 wap = wac-wat;	 lat = lac+wap;
hac = 0.0520 * hpaud;
bac = bwin*hpaud - 7*dy;    % bottom of windowed waveform plot
ba2 = bac+hac+dy;	ba3 = ba2+hac+dy;
bwf = bwav*hpaud - 9*dy;    % bottom of main waveform plot
bw2 = bwf+hac+dy;	bw3 = bw2+hac+dy;
pbSELAL	= uicontrol( 'Parent',pAUD, 'Style','pushbutton', 'String','All', ...
                     'FontSize',fss, 'Callback',{@aud_selall_CB}, ...
                     'Position',[lac bw3 wac hac] );
pbWFRML	= uicontrol( 'Parent',pAUD, 'Style','pushbutton', 'String','L', ...
                     'FontSize',fst, 'Callback',{@frame_left_CB}, ...
                     'Position',[lac bw2 wap hac] );
txWFRML	= uicontrol( 'Parent',pAUD, 'Style','Edit', 'String',num2str(dat.vid.fint(1),'%d'), ...
                     'FontSize',fst, 'Callback',{@edit_wflim_CB}, 'BackgroundColor',bgd, ...
                     'Position',[lat bw2 wat hac], 'HorizontalAlignment','right' );
pbWFRMR	= uicontrol( 'Parent',pAUD, 'Style','pushbutton', 'String','R', ...
                     'FontSize',fst, 'Callback',{@frame_right_CB}, ...
                     'Position',[lac bwf wap hac] );
txWFRMR	= uicontrol( 'Parent',pAUD, 'Style','Edit', 'String',num2str(dat.vid.fint(2),'%d'), ...
                     'FontSize',fst, 'Callback',{@edit_wflim_CB}, 'BackgroundColor',bgd, ...
                     'Position',[lat bwf wat hac], 'HorizontalAlignment','right' );
pbWWINC	= uicontrol( 'Parent',pAUD, 'Style','pushbutton', 'String','< >', ...
                     'FontSize',fss, 'Callback',{@alt_winwd_CB, +0.1}, ...
                     'Position',[lac ba3 wac hac]  );
txWINWD	= uicontrol( 'Parent',pAUD, 'Style','Edit', 'String',num2str(dat.aud.wwd,'%.2f'), ...
                     'FontSize',fss, 'Callback',{@edit_winwd_CB}, 'BackgroundColor',bgd, ...
                     'Position',[lac ba2 wac hac] );
pbWWDEC	= uicontrol( 'Parent',pAUD, 'Style','pushbutton', 'String','> <', ...
                     'FontSize',fss, 'Callback',{@alt_winwd_CB, -0.1}, ...
                     'Position',[lac bac wac hac]  );

% compute spectrogram of windowed audio signal
fFFT = subplot( 4,1,[3 4], 'Parent',pAUD );
set( fFFT, 'Visible','off', 'FontSize',fss );
set( fFFT, 'Position',[lwav bsgm wwav hsgm] )
if ~strcmp( dat.ffn_aud,'' )
    dat.sg.ft	= 0;            % formant tracking flag
    dat.sg.fmt	= 0;            % formant detection flag
    dat.sg.pree	= sg_pree;      % signal pre-emphasis flag
    dat.sg.win	= sg_win;       % spectrogram window width (msec)
    dat.sg.step	= sg_step;      % spectrogram window step length (msec)
    dat.sg.Flim	= sg_Flim;      % spectrogram maximum frequency
    dat.sg.nfft	= sg_nfft;      % spectrogram number of FFT bins
    dat.sg.nf	= sg_nf;        % spectrogram number of formants to track
    dat.sg.fwin	= 3*sg_win;     % width of formant tracking window
    dat.sg.fstp	= 3*sg_step;	% width of formant tracking step
    [SG,F,t]	= calc_sg( dat.aud.sig, dat.aud.Fs, 0, sg_pree, sg_win, sg_step, sg_Flim, sg_nfft );
    dat.sg.t	= t;            % spectrogram time value vector
    dat.sg.f	= F;            % spectrogram frequency value vector
    dat.sg.SG	= SG;           % spectrogram value matrix

    % spectrogram GUI controls: adjust spectrogram parameters
    hsg = 0.0400 * hpaud;       % height of spectrogram control text box
    bsg = bsgm*hpaud - 5*dy;	% bottom of spectrogram plot
                        bf1 = bsg+hsg-1;
    bs2 = bf1+hsg+dy;   bf2 = bs2+hsg-1;
    bs3 = bf2+hsg+dy;   bf3 = bs3+hsg-1;
    bs4 = bf3+hsg+dy;   bf4 = bs4+hsg-1;
    bs5 = bf4+hsg+dy;   bf5 = bs5+hsg-1;
    lbSGWIN	= uicontrol( 'Parent',pAUD, 'Style','Edit', 'String','Win', ...
                         'Enable','Inactive', ...
                         'FontSize',fst, 'Position',[lac bf5 wac hsg] );
    txSGWIN	= uicontrol( 'Parent',pAUD, 'Style','Edit', 'String',num2str(dat.sg.win,'%1d'), ...
                         'BackgroundColor',bgd,	'Callback',{@edit_sgwin_CB}, ...
                         'FontSize',fst, 'Position',[lac bs5 wac hsg] );
    lbSGSTP	= uicontrol( 'Parent',pAUD, 'Style','Edit', 'String','Step', ...
                         'Enable','Inactive', ...
                         'FontSize',fst, 'Position',[lac bf4 wac hsg] );
    txSGSTP	= uicontrol( 'Parent',pAUD, 'Style','Edit', 'String',num2str(dat.sg.step,'%1d'), ...
                         'BackgroundColor',bgd,	'Callback',{@edit_sgstep_CB}, ...
                         'FontSize',fst, 'Position',[lac bs4 wac hsg] );
    lbSGFMX	= uicontrol( 'Parent',pAUD, 'Style','Edit', 'String','Fmax', ...
                         'Enable','Inactive', ...
                         'FontSize',fst, 'Position',[lac bf3 wac hsg] );
    txSGFMX	= uicontrol( 'Parent',pAUD, 'Style','Edit', 'String',num2str(dat.sg.Flim,'%1d'), ...
                         'BackgroundColor',bgd,	'Callback',{@edit_sgfmax_CB}, ...
                         'FontSize',fst, 'Position',[lac bs3 wac hsg] );
    lbFTWIN	= uicontrol( 'Parent',pAUD, 'Style','Edit', 'String','FWin', ...
                         'Enable','Inactive', ...
                         'FontSize',fst, 'Position',[lac bf2 wac hsg] );
    txFTWIN	= uicontrol( 'Parent',pAUD, 'Style','Edit', 'String',num2str(dat.sg.fwin,'%1d'), ...
                         'BackgroundColor',bgd,	'Callback',{@edit_ftwin_CB}, ...
                         'FontSize',fst, 'Position',[lac bs2 wac hsg] );
    lbFSTEP	= uicontrol( 'Parent',pAUD, 'Style','Edit', 'String','FStep', ...
                         'Enable','Inactive', ...
                         'FontSize',fst, 'Position',[lac bf1 wac hsg] );
    txFSTEP	= uicontrol( 'Parent',pAUD, 'Style','Edit', 'String',num2str(dat.sg.fstp,'%1d'), ...
                         'BackgroundColor',bgd,	'Callback',{@edit_ftstep_CB}, ...
                         'FontSize',fst, 'Position',[lac bsg wac hsg] );
    bfm = 2*dy;
    lpc = lac-wcb-2*dx;	 lpr = lpc-wac;
    pbFTTOG	= uicontrol( 'Parent',pAUD, 'Style','pushbutton', 'String','Fmt', ...
                         'FontSize',fst, 'Position',[lac bfm wac hcb],	'Callback',{@ft_onoff_CB} );
    cbPETOG	= uicontrol( 'Parent',pAUD, 'Style','checkbox', 'Value',1, ...
                         'FontSize',fst, 'Position',[lpc bfm wcb hcb],	'Callback',{@preemph_onoff_CB} );
    lbPETOG	= uicontrol( 'Parent',pAUD, 'Style','Edit', 'String','PreE', 'Enable','Inactive', ...
                         'FontSize',fst, 'Position',[lpr bfm wac hl1]);
end
                 
% fetch original video from AVI file
if ( exist( dat.ffn_vid, 'file' ) )
    fprintf( '   Using video:          <%s>\n', dat.ffn_vid );
    vobj = VideoReader(dat.ffn_vid);
    hh	 = vobj.Height;
    ww	 = vobj.Width;
    %nf	 = vobj.NumberOfFrames;
    nf	 = vobj.Duration * vobj.FrameRate;
    avi(1:nf) = struct('cdata', zeros(hh,ww,3,'uint8'), 'colormap',[] );
    for f = 1:nf
        avi(f).cdata = read(vobj,f);
    end
    dat.vid.avi	 = avi;
    dat.vid.info = get(vobj);
    dat.vid.info.NumberOfFrames = nf;
    if ( isempty(dat.vid.info) )
        fprintf( '   Can''t load video from file <%s>\n', dat.ffn_vid );
        return;
    end
else
    fprintf( '   Can''t find video file <%s>\n', dat.ffn_vid );
    return;
end

% fetch tracked video from AVI file
if ( exist( dat.ffn_trk, 'file' ) )
    fprintf( '   Using tracked video:  <%s>\n', dat.ffn_trk );

    vobj = VideoReader(dat.ffn_trk);
    hh	 = vobj.Height;
    ww	 = vobj.Width;
    %nf	 = vobj.NumberOfFrames;
    nf	 = vobj.Duration * vobj.FrameRate;
    avi(1:nf) = struct('cdata', zeros(hh,ww,3,'uint8'), 'colormap',[] );
    for f = 1:nf
        avi(f).cdata = read(vobj,f);
    end
    dat.tvd.avi	 = avi;
    dat.tvd.info = get( vobj );
    dat.tvd.info.NumberOfFrames = nf;
    if ( isempty(dat.tvd.info) )
        fprintf( '   Can''t load video from file <%s>\n', dat.ffn_trk );
        return;
    end
else
    fprintf( '   No tracked-video found ... \n' );
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
    nfv	= dat.vid.info.NumberOfFrames;
    frv	= dat.vid.info.FrameRate;
    vlv	= nfv/frv;
end
% construct dummy audio signal of correct length if no WAV file found ....
if strcmp( dat.ffn_aud,'' )
    dat.aud.Fs	= 20000;
    dat.aud.dur	= dat.vid.info.NumberOfFrames/dat.vid.info.FrameRate;
    dat.aud.len	= round(dat.aud.dur * dat.aud.Fs);
    dat.aud.sig	= zeros(1,dat.aud.len);
    dat.aud.tvec = (1:dat.aud.len)./dat.aud.Fs;
    fprintf('\n   Using %0.2f seconds of null audio in lieu of matching WAV file', dat.aud.dur );
else
    fprintf('\n   Found %0.2f seconds of audio         (%d samples at %d Hz)', dat.aud.dur, dat.aud.len, dat.aud.Fs );
end
fprintf('\n   Found %0.2f seconds of video         (%d frames  at %0.2f fps)\n', vlv,   nfv,   frv  );
if ~isempty(dat.tvd.info)
    nft	= dat.tvd.info.NumberOfFrames;
    frt	= dat.tvd.info.FrameRate;
    vlt	= nft/frt;
    fprintf(  '   Found %0.2f seconds of tracked video (%d frames  at %0.2f fps)\n', vlt, nft, frt );
end
disp(' ');

% create video display axes
lvax = 0.110;	bvax = 0.250;
wvax = 0.805;	hvax = wvax;
dat.vid.pos	= [lvax bvax wvax hvax];

if (nargin>1)
    dat.vid.grid            = grid;     % fetch grid parameters specified in input argument
    old_dat.vid.grid        = grid;     % and overwrite previously specified grid in reloaded data
else                                    % ... or set initial params for analysis grid if not specified
    dat.vid.grid.or1        = [0 0];    % alveolar origin (px)
    dat.vid.grid.or2        = [0 0];    % lingual origin (px)
    dat.vid.grid.mpal       = [0 0];    % mid palatal location (px)   
    dat.vid.grid.dent       = [0 0];    % alveolar location (px)
    dat.vid.grid.glot       = [0 0];    % glottal location (px)
    dat.vid.grid.mlab       = [0 0];    % mid-labial point (px)
    dat.vid.grid.rad1       = 0;        % semipolar grid central radius (px)
    dat.vid.grid.rad2       = 0;        % alveolar subgrid central radius (px)
    dat.vid.grid.theta      = 0;        % angle of ray joining alveolar & lingual origins (rad)
    dat.vid.grid.wid        = gwidth;	% width of analysis gridlines (px)
    dat.vid.grid.gint       = gint;     % interval between gridlines (px)
    dat.vid.grid.nlines     = 0;        % number of gridlines
    dat.vid.grid.linewd     = 1;        % thickness of gridlines (pt)
    dat.vid.grid.ends       = [];       % endpoint coords of gridlines (px)
    dat.vid.grid.txt        = [];       % coords of gridline labels (px)
    dat.vid.grid.pts        = [];       % points defining gridlines (px)
    dat.vid.grid.palfnum    = 0;        % index of palatal reference frame (initially undefined)
    dat.vid.grid.pal        = [];       % range of gridlines defining palate
    dat.vid.grid.palate     = [];       % set of points defining palatal boundary
    dat.vid.grid.phafnum    = 0;        % index of pharyngeal reference frame (initially undefined)
    dat.vid.grid.pha        = [];       % range of gridlines defining pharynx
    dat.vid.grid.pharynx	= [];       % set of points defining pharyngeal boundary
    dat.vid.grid.tng        = [];       % range of gridlines defining tongue
    dat.vid.grid.shifted    = 0;        % global flag: interfame corelation shift algorithm has been applied
    dat.vid.grid.ithrsh     = 0.55;     % default intensity threshold
    dat.vid.grid.diwtcl     = 0.15;     % default distance/intensity threshold
    dat.vid.grid.diwtib     = 0.15;     % default distance/intensity threshold
end

fVID = axes( 'Parent',pVID, 'Position',dat.vid.pos, 'Visible','Off', 'FontSize',fss );


% video playback frame rate control
wvr	= 0.07*wpvid;   lvr	= wpvid-wvr-6;
hvr	= 1.25*hl1;     bvl	= hpvid-hvr-26;     bvr = bvl-hvr+1;
vr	= num2str(dat.vid.info.FrameRate,'%0.1f');                 
lbVRate	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String','fps',	...
                     'FontSize',fss,      ...
                     'Position',[lvr bvl wvr hvr],	'Enable','Inactive' );
txVRate	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String',vr,                            ...
                     'FontSize',fss, 'HorizontalAlignment','right', 'BackgroundColor',bgd,	...
                     'Position',[lvr bvr wvr hvr] );
                 
% video axis checkbox controls
wla = wvr-wcb+5;    lac	= lvr+wla-2;
ba1	= bvr-hl1-9;    ba2	= ba1-hl1;	ba3	= ba2-hl1;
ba4	= ba3-hvr-5;	ba5	= ba4-hvr;  bab = 0.305*hpvid;
lbSHVAX	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String','Ax',	...
                     'FontSize',fsm, 'HorizontalAlignment','left', ...
                     'Position',[lvr ba1 wla hl1],	'Enable','Inactive' );
cbSHVAX	= uicontrol( 'Parent',pVID, 'Style','checkbox', 'Value',0, ...
                     'Position',[lac ba1 wcb hcb],	'Callback', {@toggle_vidaxes_CB}  );
lbSCVAX	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String','Sc',	...
                     'FontSize',fsm, 'HorizontalAlignment','left',	'Visible','off', ...
                     'Position',[lvr ba2 wla hl1],	'Enable','Inactive' );
cbSCVAX	= uicontrol( 'Parent',pVID, 'Style','checkbox', 'Value',0,	'Visible','off', ...
                     'Position',[lac ba2 wcb hcb],	'Callback', {@toggle_vidaxes_CB} );
lbSHVMM	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String','mm',	...
                     'FontSize',fsm, 'HorizontalAlignment','left',	'Visible','off', ...
                     'Position',[lvr ba3 wla hl1],	'Enable','Inactive' );
cbSHVMM	= uicontrol( 'Parent',pVID, 'Style','checkbox', 'Value',0,	'Visible','off', ...
                     'Position',[lac ba3 wcb hcb],	'Callback', {@toggle_vidaxes_CB} );
lbSHFOV	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String','FOV',	...
                     'FontSize',fst, 'HorizontalAlignment','left',	'Visible','off', ...
                     'Position',[lvr ba4 wvr hvr],	'Enable','Inactive' );
txSHFOV	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String', num2str( dat.status.FOV,'%d' ), ...
                     'FontSize',fst, 'BackgroundColor',bgd,         'Visible','off', ...
                     'Position',[lvr ba5 wvr hvr],	'Callback', {@update_FOV_CB} );
pbFIXPT	= uicontrol( 'Parent',pVID, 'Style','pushbutton', 'String','Fix',	...
                     'FontSize',fss, 'HorizontalAlignment','left',	'Visible','off', ...
                     'Position',[lvr bab wvr hvr],	'Callback', {@edit_boundaries_CB}  );

% pixel location display controls
bax = bab+hvr+2;
txSHPOX	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String', num2str( 0,'%d' ), ...
                     'FontSize',fst, 'BackgroundColor',bgd,	'Visible','off', ...
                     'Position',[lvr bax wvr hvr] );
txSHPOY	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String', num2str( 0,'%d' ), ...
                     'FontSize',fst, 'BackgroundColor',bgd,	'Visible','off', ...
                     'Position',[lvr bab wvr hvr] );
                 
                 
% video control placement parameters
boff = 5;
lvc = 10;
txtb = 7;	txtw = 60;	txth = 24;  lblw = 80;  lbw2 = 60; 
sldw = 406;	sldh = 18;  sldb = txth+13;
flft = lblw+lvc;   ftlf = flft+txtw;
vlf2 = ftlf+lbw2+boff;   
tlft = lblw+vlf2;   ttlf = tlft+txtw;

% video time and frame selection controls
slMin	= 1;	slMax = dat.vid.info.NumberOfFrames;
slRng	= (slMax-slMin);
slStpSm = slMin/slRng;	slStpLg = 5*slStpSm;
slVINC	= uicontrol( 'Parent',pVID, 'Style','Slider', 'BackgroundColor',bgd, ...
                     'Value',dat.fnum, 'Min',slMin, 'Max',slMax, 'SliderStep',[slStpSm slStpLg], ...
                                     'Position',[lvc sldb sldw sldh],	'Callback',{@update_frame_CB} );
lbVFRAM	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String','Frame No.:', 'Enable','Inactive', ...
                     'FontSize',fsl, 'Position',[lvc txtb lblw txth] );
txVFRAM	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String',num2str(dat.fnum,'%d'), ...
                     'FontSize',fsl, 'Position',[flft txtb txtw txth], ...
                     'BackgroundColor',bgd,	'Callback',{@edit_frameno_CB} );
lbVFTOT	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String',[' of ' num2str(slMax)], ...
                     'FontSize',fsl, 'Position',[ftlf txtb lbw2 txth], ...
                     'Enable','Inactive', 'HorizontalAlignment','left' );
lbVTIME	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String','Time (sec):', 'Enable','Inactive', ...
                     'FontSize',fsl, 'Position',[vlf2 txtb lblw txth] );
txVTIME	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String',num2str(dat.time,'%.3f'), ...
                     'FontSize',fsl, 'Position',[tlft txtb txtw txth], ...
                     'BackgroundColor',bgd,	'Callback',{@edit_time_CB} );
lbVTTOT	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String',[' of ' num2str(dat.aud.dur,'%.2f')], ...
                     'FontSize',fsl, 'Position',[ttlf txtb lbw2 txth], ...
                     'Enable','Inactive', 'HorizontalAlignment','left' );
                 

% calculate locations of video display uicontrols
bvc = 0.120*hpvid;      % bottom location of video display controls
lcf = lvc;              % left location of crop frame controls
lgp = 0.215*wpvid;      % left location of grid placement controls
llv = 0.500*wpvid;      % left location of video display cbs
                 
% calculate locations of crop-frame controls
wtc	= 0.058*wpvid;      % widths of crop-frame placement controls
wpc	= 0.450*wtc;
htc = 0.048*hpvid;      % heights of crop-frame placement controls
hpc	= 0.450*htc;
lf1 = lcf+wpc;          % left locations of crop-frame controls
lf2 = lf1+wtc;
lf3 = lf2+wtc;
lfc = lf3-1;
wlb	= lf3-lcf;
bf2 = bvc+hpc;          % bottom locations of crop-frame controls
bf3 = bf2+hpc;
bf4 = bf3+hpc;
bf5 = bf4+hpc+dy2;
bf6 = bf5+hpc+dy2;

% place crop-frame controls on GUI
lbSHOWF	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String','Crop Frame', ...
                     'FontSize',fst, 'HorizontalAlignment','left', ...
                     'Position',[lcf bf5 wlb hl1],	'Enable','Inactive' );
cbSHOWF	= uicontrol( 'Parent',pVID, 'Style','checkbox', 'Value',0, ...
                     'Position',[lfc bf5 wcb hcb],	'Callback', {@toggle_frame_CB} );
lbSHOFN	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String','Frame No.', 'Visible','off', ...
                     'FontSize',fst, 'HorizontalAlignment','left', ...
                     'Position',[lcf bf6 wlb hl1],	'Enable','Inactive' );
cbSHOFN	= uicontrol( 'Parent',pVID, 'Style','checkbox', 'Value',0,  'Visible','off', ...
                     'Position',[lfc bf6 wcb hcb],	'Callback', {@refresh_GUI} );
txFOFFx	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String','', ...
                     'FontSize',fst, 'BackgroundColor',bgd, ...
                     'Position',[lf1 bf3 wtc htc],	'Callback', {@set_frm_pos_CB} );
pbFRLF	= uicontrol( 'Parent',pVID, 'Style','pushbutton', 'String','<', ...
                     'FontSize',fsm, 'ForegroundColor','black', ...
                     'Position',[lcf bf4 wpc hpc],	'Callback', {@alt_frm_pos_CB,-10, 0} );
pbFRRT	= uicontrol( 'Parent',pVID, 'Style','pushbutton', 'String','>', ...
                     'FontSize',fsm, 'ForegroundColor','black', ...
                     'Position',[lcf bf3 wpc hpc],	'Callback', {@alt_frm_pos_CB, 10, 0} );
txFOFFy	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String','', ...
                     'FontSize',fst, 'BackgroundColor',bgd, ...
                     'Position',[lf2 bf3 wtc htc],	'Callback', {@set_frm_pos_CB} );
pbFRUP	= uicontrol( 'Parent',pVID, 'Style','pushbutton', 'String','^', ...
                     'FontSize',fsm, 'ForegroundColor','black', ...
                     'Position',[lf3 bf4 wpc hpc],	'Callback', {@alt_frm_pos_CB, 0,-10} );
pbFRDN	= uicontrol( 'Parent',pVID, 'Style','pushbutton', 'String','v', ...
                     'FontSize',fsm, 'ForegroundColor','black', ...
                     'Position',[lf3 bf3 wpc hpc],	'Callback', {@alt_frm_pos_CB, 0, 10} );
txFRMW	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String','', ...
                     'FontSize',fst, 'BackgroundColor',bgd, ...
                     'Position',[lf1 bvc wtc htc],	'Callback', {@set_frm_size_CB} );
pbFINCW	= uicontrol( 'Parent',pVID, 'Style','pushbutton', 'String','<', ...
                     'FontSize',fsm, 'ForegroundColor','black', ...
                     'Position',[lcf bf2 wpc hpc],	'Callback', {@alt_frm_size_CB,-10, 0} );
pbFDECW	= uicontrol( 'Parent',pVID, 'Style','pushbutton', 'String','>', ...
                     'FontSize',fsm, 'ForegroundColor','black', ...
                     'Position',[lcf bvc wpc hpc],	'Callback', {@alt_frm_size_CB, 10, 0} );
txFRMH	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String','', ...
                     'FontSize',fst, 'BackgroundColor',bgd, ...
                     'Position',[lf2 bvc wtc htc],	'Callback', {@set_frm_size_CB} );
pbFINCH	= uicontrol( 'Parent',pVID, 'Style','pushbutton', 'String','^', ...
                     'FontSize',fsm, 'ForegroundColor','black', ...
                     'Position',[lf3 bf2 wpc hpc],	'Callback', {@alt_frm_size_CB, 0, 10} );
pbFDECH	= uicontrol( 'Parent',pVID, 'Style','pushbutton', 'String','v', ...
                     'FontSize',fsm, 'ForegroundColor','black', ...
                     'Position',[lf3 bvc wpc hpc],	'Callback', {@alt_frm_size_CB, 0,-10} );

                 
% calculate locations of analysis grid controls
lg1 = lgp+wpc;      % left locations of grid controls
lg2 = lg1+wtc;
lg3 = lg2+wtc;
lg4 = lg3+wpc+dx2;
lg5 = lg4+wtc;
wlb	= lg5-lgp;
lgs = lg5-wcb;
lgl = lgs-wcb;
lcb = lg5-1;
bg2 = bvc+hpc;      % bottom locations of grid controls
bg3 = bg2+hpc;
bg4 = bg3+hpc;

% place analysis grid controls on GUI
lbSHOWG	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String','Analysis Grid',	...
                     'FontSize',fst, 'HorizontalAlignment','left', ...
                     'Position',[lgp bf5 wlb hl1],	'Enable','Inactive' );
txGDLNW	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String', num2str(dat.vid.grid.linewd), ...
                     'FontSize',fst, 'BackgroundColor',bgd, ...
                     'Position',[lgs bf5 wcb hl1],	'Callback', {@set_glinewd_CB} );
txGDSPC	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String', num2str(dat.vid.grid.gint), ...
                     'FontSize',fst, 'BackgroundColor',bgd, ...
                     'Position',[lgl bf5 wcb hl1],	'Callback', {@set_gspc_CB} );
cbSHOWG	= uicontrol( 'Parent',pVID, 'Style','checkbox', 'Value',0, ...
                     'Position',[lcb bf5 wcb hcb],	'Callback', {@toggle_grid_CB} );

txOR1X	= uicontrol( 'Parent',pVID, 'Style','Edit',	'String', num2str(dat.vid.grid.or1(1)), ...
                     'FontSize',fss, 'BackgroundColor',bgd, ...
                     'Position',[lg1 bg3 wtc htc],	'Callback', {@set_gor1_CB} );
pbOR1XD	= uicontrol( 'Parent',pVID, 'Style','pushbutton', 'String','<', ...
                     'FontSize',fsm, 'ForegroundColor','black', ...
                     'Position',[lgp bg4 wpc hpc],	'Callback', {@alt_gor1_CB,-5, 0} );
pbOR1XI	= uicontrol( 'Parent',pVID, 'Style','pushbutton', 'String','>', ...
                     'FontSize',fsm, 'ForegroundColor','black', ...
                     'Position',[lgp bg3 wpc hpc],	'Callback', {@alt_gor1_CB,+5, 0} );
txOR1Y	= uicontrol( 'Parent',pVID, 'Style','Edit',	'String', num2str(dat.vid.grid.or1(2)), ...
                     'FontSize',fss, 'BackgroundColor',bgd, ...
                     'Position',[lg2 bg3 wtc htc],	'Callback', {@set_gor1_CB} );
pbOR1YD	= uicontrol( 'Parent',pVID, 'Style','pushbutton', 'String','^', ...
                     'FontSize',fsm, 'ForegroundColor','black', ...
                     'Position',[lg3 bg4 wpc hpc],	'Callback', {@alt_gor1_CB, 0,-5} );
pbOR1YI	= uicontrol( 'Parent',pVID, 'Style','pushbutton', 'String','v', ...
                     'FontSize',fsm, 'ForegroundColor','black', ...
                     'Position',[lg3 bg3 wpc hpc],	'Callback', {@alt_gor1_CB, 0,+5} );
txOR2X	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String', num2str(dat.vid.grid.or2(1)), ...
                     'FontSize',fss, 'BackgroundColor',bgd, ...
                     'Position',[lg1 bvc wtc htc],	'Callback', {@set_gor2_CB} );
pbOR2XD	= uicontrol( 'Parent',pVID, 'Style','pushbutton', 'String','<', ...
                     'FontSize',fsm, 'ForegroundColor','black',  ...
                     'Position',[lgp bg2 wpc hpc],	'Callback', {@alt_gor2_CB,-5, 0} );
pbOR2XI	= uicontrol( 'Parent',pVID, 'Style','pushbutton', 'String','>', ...
                     'FontSize',fsm, 'ForegroundColor','black',  ...
                     'Position',[lgp bvc wpc hpc],	'Callback', {@alt_gor2_CB,+5, 0} );
txOR2Y	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String', num2str(dat.vid.grid.or2(2)), ...
                     'FontSize',fss, 'BackgroundColor',bgd, ...
                     'Position',[lg2 bvc wtc htc],	'Callback', {@set_gor2_CB} );
pbOR2YD	= uicontrol( 'Parent',pVID, 'Style','pushbutton', 'String','^', ...
                     'FontSize',fsm, 'ForegroundColor','black', ...
                     'Position',[lg3 bg2 wpc hpc],	'Callback', {@alt_gor2_CB, 0,-5} );
pbOR2YI	= uicontrol( 'Parent',pVID, 'Style','pushbutton', 'String','v', ...
                     'FontSize',fsm, 'ForegroundColor','black', ...
                     'Position',[lg3 bvc wpc hpc],	'Callback', {@alt_gor2_CB, 0,+5} );

txGRAD	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String', num2str(dat.vid.grid.rad1,'%d'), ...
                     'FontSize',fss, 'BackgroundColor',bgd, ...
                     'Position',[lg4 bg3 wtc htc],	'Callback', {@set_grad_CB} );
txGWID	= uicontrol( 'Parent',pVID, 'Style','Edit',	'String', num2str(dat.vid.grid.wid,'%d'), ...
                     'FontSize',fss, 'BackgroundColor',bgd, ...
                     'Position',[lg4 bvc wtc htc],	'Callback', {@set_gwid_CB} );
pbGRADI	= uicontrol( 'Parent',pVID, 'Style','pushbutton', 'String','^', ...
                     'FontSize',fsm, 'ForegroundColor','black', ...
                     'Position',[lg5 bg4 wpc hpc],	'Callback', {@alt_grad_CB,+5} );
pbGRADD	= uicontrol( 'Parent',pVID, 'Style','pushbutton', 'String','v', ...
                     'FontSize',fsm, 'ForegroundColor','black', ...
                     'Position',[lg5 bg3 wpc hpc],	'Callback', {@alt_grad_CB,-5} );
pbGWIDI	= uicontrol( 'Parent',pVID, 'Style','pushbutton', 'String','^', ...
                     'FontSize',fsm, 'ForegroundColor','black', ...
                     'Position',[lg5 bg2 wpc hpc],	'Callback', {@alt_gwid_CB,+5} );
pbGWIDD	= uicontrol( 'Parent',pVID, 'Style','pushbutton', 'String','v', ...
                     'FontSize',fsm, 'ForegroundColor','black', ...
                     'Position',[lg5 bvc wpc hpc],	'Callback', {@alt_gwid_CB,-5} );

                 
% calculate locations of video display checkbox controls
wlv = 0.115*wpvid+1;	% width of video display cb labels
lcv = llv+wlv-1;	
bv2 = bvc+hl1;          % bottom of nth video control cb
bv3 = bv2+hl1;
bv4 = bv3+hl1;

% place video display control checkboxes
%   leftmost column
lbSHOWC	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String','VT Center',	...
                     'FontSize',fst, 'HorizontalAlignment','left', ...
                     'Position',[llv bv4 wlv hl1],	'Enable','Inactive' );
cbSHOWC	= uicontrol( 'Parent',pVID, 'Style','checkbox', 'Value',0, ...
                     'Position',[lcv bv4 wcb hcb],	'Callback', {@refresh_GUI} );
lbSHOWB	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String','Boundaries',	...
                     'FontSize',fst, 'HorizontalAlignment','left', ...
                     'Position',[llv bv3 wlv hl1],	'Enable','Inactive' );
cbSHOWB	= uicontrol( 'Parent',pVID, 'Style','checkbox', 'Value',0, ...
                     'Position',[lcv bv3 wcb hcb],	'Callback', {@refresh_GUI} );
lbSHOWP	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String','Pharynx',	...
                     'FontSize',fst, 'HorizontalAlignment','left', ...
                     'Position',[llv bv2 wlv hl1],	'Enable','Inactive' );
cbSHOWP	= uicontrol( 'Parent',pVID, 'Style','checkbox', 'Value',0, ...
                     'Position',[lcv bv2 wcb hcb],	'Callback', {@refresh_GUI} );
lbSHOWR	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String','Palate',	...
                     'FontSize',fst, 'HorizontalAlignment','left', ...
                     'Position',[llv bvc wlv hl1],	'Enable','Inactive' );
cbSHOWR	= uicontrol( 'Parent',pVID, 'Style','checkbox', 'Value',0, ...
                     'Position',[lcv bvc wcb hcb],	'Callback', {@refresh_GUI} );
%   2nd column of video display control checkboxes
lv2 = lcv+wcb+2;	lcv = lv2+wlv;
wtx = wcb+5;        ltx = lcv+wcb-wtx-2;
lbSHCGR	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String','Cen Graph',	...
                     'FontSize',fst, 'HorizontalAlignment','left', ...
                     'Position',[lv2 bv4 wlv hl1],	'Enable','Inactive' );
cbSHCGR	= uicontrol( 'Parent',pVID, 'Style','checkbox', 'Value',0, ...
                     'Position',[lcv bv4 wcb hcb],	'Callback', {@refresh_GUI} );
lbSMTHB	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String','Smoothed',	...
                     'FontSize',fst, 'HorizontalAlignment','left', ...
                     'Position',[lv2 bv3 wlv hl1],	'Enable','Inactive' );
cbSHOWS	= uicontrol( 'Parent',pVID, 'Style','checkbox', 'Value',0, ...
                     'Position',[lcv bv3 wcb hcb],	'Callback', {@refresh_GUI} );
lbSHCIR	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String','Int Circle',	...
                     'FontSize',fst, 'HorizontalAlignment','left', ...
                     'Position',[lv2 bv2 wlv hl1],	'Enable','Inactive' );
cbSHCIR	= uicontrol( 'Parent',pVID, 'Style','checkbox', 'Value',0, ...
                     'Position',[lcv bv2 wcb hcb],	'Callback', {@refresh_GUI} );
lbANRAD	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String','rAnalysis', ...
                     'FontSize',fst, 'ForegroundColor','black', 'HorizontalAlignment','left', ...
                     'Position',[lv2 bvc wlv hl1],	'Enable','Inactive' );
txANRAD	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String',num2str(dat.status.radanal,'%d'), ...
                     'FontSize',fst, 'BackgroundColor',bgd, ...
                     'Position',[ltx bvc wtx hl1],	'Callback', {@set_intanalrad_CB} );
%   3rd column of video display control checkboxes
lv3 = lcv+wcb+2;	lcv = lv3+wlv;	wpb = wlv+wcb;
lbVMODE	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String','Uninterp',	...
                     'FontSize',fst, 'HorizontalAlignment','left', ...
                     'Position',[lv3 bv4 wlv hl1],   'Enable','Inactive' );
cbVMODE	= uicontrol( 'Parent',pVID, 'Style','checkbox', 'Value',0, ...
                     'Position',[lcv bv4 wcb hcb],	'Callback', {@refresh_GUI} );
lbTMODE	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String','Tracked',	...
                     'FontSize',fst, 'HorizontalAlignment','left', ...
                     'Position',[lv3 bv3 wlv hl1],   'Enable','Inactive' );
cbTMODE	= uicontrol( 'Parent',pVID, 'Style','checkbox', 'Value',0, ...
                     'Position',[lcv bv3 wcb hcb],	'Callback', {@tog_vmode_CB} );
lbICORR	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String','Corrected',	...
                     'FontSize',fst, 'HorizontalAlignment','left', ...
                     'Position',[lv3 bv2 wlv hl1],   'Enable','Inactive' );
cbICORR	= uicontrol( 'Parent',pVID, 'Style','checkbox', 'Value',0, ...
                     'Position',[lcv bv2 wcb hcb],	'Callback', {@tog_icorr_CB} );
lbSHLPF	= uicontrol( 'Parent',pVID, 'Style','Edit', 'String','LPF Vid',	...
                     'FontSize',fst, 'HorizontalAlignment','left', ...
                     'Position',[lv3 bvc wlv hl1],	'Enable','Inactive' );
cbSHLPF	= uicontrol( 'Parent',pVID, 'Style','checkbox', 'Value',0, ...
                     'Position',[lcv bvc wcb hcb],	'Callback', {@refresh_GUI} );


% place analysis frame uicontrols
%   leftmost column
lcb = lt1-wcb+2;	
pbPGRID	= uicontrol( 'Parent',pANA, 'Style','pushbutton', 'String','Place Grid', ...
                     'FontSize',fss, 'ForegroundColor',colgrn, ...
                     'Position',[le1 bp6 wp1 hp1],	'Callback', {@place_grid} );
pbGETAF	= uicontrol( 'Parent',pANA, 'Style','pushbutton', 'String','Segment Frame', ...
                     'FontSize',fss, 'ForegroundColor',colgrn, 'HorizontalAlignment','left', ...
                     'Position',[le1 bp5 wp1 hp1],	'Callback', {@find_tissue_bnds_CB} );
pbAFSEG	= uicontrol( 'Parent',pANA, 'Style','pushbutton', 'String','Segment Intvl', ...
                     'FontSize',fss, 'ForegroundColor',colgrn, ...
                     'Position',[le1 bp4 wp1 hp1],	'Callback', {@analyze_segment} );
lbUSERF	= uicontrol( 'Parent',pANA, 'Style','Edit', 'String','Use Ref Bnds',	...
                     'FontSize',fst, 'HorizontalAlignment','left', ...
                     'Position',[le1 bl4 wlc hl1],	'Enable','Inactive' );
cbUSERF	= uicontrol( 'Parent',pANA, 'Style','checkbox', 'Value',0, ...
                     'Position',[lc1 bl4 wcb hcb],	'Callback', {@refresh_GUI} );
lbDIWTC	= uicontrol( 'Parent',pANA, 'Style','Edit', 'String','Alpha', ...
                     'FontSize',fst, 'ForegroundColor','black', 'HorizontalAlignment','left', ...
                     'Position',[le1 bl3 wlt hl1],	'Enable','Inactive' );
txDIWTC	= uicontrol( 'Parent',pANA, 'Style','Edit', 'String',num2str(dat.vid.grid.diwtcl,'%1.2f'), ...
                     'FontSize',fst, 'BackgroundColor',bgd, ...
                     'Position',[lt1 bl3 wt1 hl1],	'Callback', {@set_diwtcl_CB} );
lbITHRS	= uicontrol( 'Parent',pANA, 'Style','Edit', 'String','Beta', ...
                     'FontSize',fst, 'ForegroundColor','black', 'HorizontalAlignment','left', ...
                     'Position',[le1 bl2 wlt hl1],	'Enable','Inactive' );
txITHRS	= uicontrol( 'Parent',pANA, 'Style','Edit', 'String',num2str(dat.vid.grid.ithrsh,'%1.2f'), ...
                     'FontSize',fst, 'BackgroundColor',bgd, ...
                     'Position',[lt1 bl2 wt1 hl1],	'Callback', {@set_ithresh_CB} );
lbGAMMA	= uicontrol( 'Parent',pANA, 'Style','Edit', 'String','Gamma', ...
                     'FontSize',fst, 'ForegroundColor','black', 'HorizontalAlignment','left', ...
                     'Position',[le1 bo1 wlt hl1],	'Enable','Inactive' );
cbGAMMA	= uicontrol( 'Parent',pANA, 'Style','checkbox', 'Value',0, ...
                     'Position',[lcb bo1 wcb hcb],	'Callback', {@refresh_GUI} );
txGAMMA	= uicontrol( 'Parent',pANA, 'Style','Edit', 'String',num2str(dat.vid.grid.diwtib,'%1.2f'), ...
                     'FontSize',fst, 'BackgroundColor',bgd, ...
                     'Position',[lt1 bo1 wt1 hl1],	'Callback', {@set_gamma_CB} );
                 
%	tract analysis controls: col 2
wgt = 0.25*wp1; wgl = wp1-2*wgt;
lg1 = lp2+wgl;	lg2 = lg1+wgt;
pbPAREF	= uicontrol( 'Parent',pANA, 'Style','pushbutton', 'String','Edit Palate', ...
                     'FontSize',fss, 'ForegroundColor',colgrn, 'HorizontalAlignment','left', ...
                     'Position',[lp2 bp6 wp1 hp1],	'Callback', {@edit_palate_CB}	);
pbAVEPH	= uicontrol( 'Parent',pANA, 'Style','pushbutton', 'String','Mean Pharynx', ...
                     'FontSize',fss, 'ForegroundColor',colgrn, 'HorizontalAlignment','left', ...
                     'Position',[lp2 bp5 wp1 hp1],	'Callback', {@mean_pharynx_CB}	);
pbPHREF	= uicontrol( 'Parent',pANA, 'Style','pushbutton', 'String','Edit Pharynx', ...
                     'FontSize',fss, 'ForegroundColor',colgrn, 'HorizontalAlignment','left', ...
                     'Position',[lp2 bp4 wp1 hp1],	'Callback', {@edit_pharynx_CB}	);
pbSHOWT	= uicontrol( 'Parent',pANA, 'Style','pushbutton', 'String','Lingual Movt', ...
                     'FontSize',fss, 'ForegroundColor',colgrn, ...
                     'Position',[lp2 bp3 wp1 hp1],	'Callback', {@plot_lingualmovt_CB} );
lbPALIM	= uicontrol( 'Parent',pANA, 'Style','Edit',       'String','Palate',	...
                     'FontSize',fst, 'HorizontalAlignment','left', ...
                     'Position',[lp2 bl3 wgl hl1],   'Enable','Inactive' );
txPAL1	= uicontrol( 'Parent',pANA, 'Style','Edit',       'String',num2str(0,'%d'), ...
                     'FontSize',fst, 'BackgroundColor',bgd, 'HorizontalAlignment','right', ...
                     'Position',[lg1 bl3 wgt hl1],	'Callback', {@set_refpalgl_CB}	);
txPAL2	= uicontrol( 'Parent',pANA, 'Style','Edit',       'String',num2str(0,'%d'), ...
                     'FontSize',fst, 'BackgroundColor',bgd, 'HorizontalAlignment','right', ...
                     'Position',[lg2 bl3 wgt hl1],	'Callback', {@set_refpalgl_CB}	);
lbPHLIM	= uicontrol( 'Parent',pANA, 'Style','Edit',       'String','Pharynx',	...
                     'FontSize',fst, 'HorizontalAlignment','left', ...
                     'Position',[lp2 bl2 wgl hl1],   'Enable','Inactive' );
txPHA1	= uicontrol( 'Parent',pANA, 'Style','Edit',       'String',num2str(0,'%d'), ...
                     'FontSize',fst, 'BackgroundColor',bgd, 'HorizontalAlignment','right', ...
                     'Position',[lg1 bl2 wgt hl1],	'Callback', {@set_refphagl_CB}	);
txPHA2	= uicontrol( 'Parent',pANA, 'Style','Edit',       'String',num2str(0,'%d'), ...
                     'FontSize',fst, 'BackgroundColor',bgd, 'HorizontalAlignment','right', ...
                     'Position',[lg2 bl2 wgt hl1],	'Callback', {@set_refphagl_CB}	);
lbSHIFT	= uicontrol( 'Parent',pANA, 'Style','Edit',       'String','Shift Refs.',	...
                     'FontSize',fst, 'HorizontalAlignment','left', ...
                     'Position',[lp2 bo1 wlc hl1],	'Enable','Inactive' );
cbSHIFT	= uicontrol( 'Parent',pANA, 'Style','checkbox', 'Value',0, ...
                     'Position',[lc2 bo1 wcb hcb],	'Callback', {@shift_tract_bnds_CB} );
                 
%	tract analysis controls: col 3
lg1 = lp3+wgl;	lg2 = lg1+wgt;
pbAFFRM	= uicontrol( 'Parent',pANA, 'Style','pushbutton', 'String','Area Function', ...
                     'FontSize',fss, 'ForegroundColor',colgrn, ...
                     'Position',[lp3 bp6 wp1 hp1],	'Callback', {@extract_AF_frame_CB} );
pbDCTSM	= uicontrol( 'Parent',pANA, 'Style','pushbutton', 'String','Smooth Tongue', ...
                     'FontSize',fss, 'ForegroundColor',colgrn, ...
                     'Position',[lp3 bp5 wp1 hp1],	'Callback', {@smooth_tongue_CB} );
pbSAVEM	= uicontrol( 'Parent',pANA, 'Style','pushbutton', 'String','Manu', ...
                     'FontSize',fss, 'ForegroundColor',colgrn, ...
                     'Position',[lp3 bp4 wgl hp1],	'Callback', {@save_segment_CB,'manual'} );
pbSAVEA	= uicontrol( 'Parent',pANA, 'Style','pushbutton', 'String','Auto', ...
                     'FontSize',fss, 'ForegroundColor',colgrn, ...
                     'Position',[lg1 bp4 wgl hp1],	'Callback', {@save_segment_CB,'auto'} );
pbCMPAF	= uicontrol( 'Parent',pANA, 'Style','pushbutton', 'String','Compare AF', ...
                     'FontSize',fss, 'ForegroundColor',colgrn, ...
                     'Position',[lp3 bp3 wp1 hp1],	'Callback', {@compare_AF_CB} );
lbUSERF	= uicontrol( 'Parent',pANA, 'Style','Edit',       'String','Use Smooth',	...
                     'FontSize',fst, 'HorizontalAlignment','left', ...
                     'Position',[lp3 bl3 wlc hl1],	'Enable','Inactive' );
cbUSESM	= uicontrol( 'Parent',pANA, 'Style','checkbox', 'Value',0, ...
                     'Position',[lc3 bl3 wcb hcb],	'Callback', {@refresh_GUI} );
lbTLIMS	= uicontrol( 'Parent',pANA, 'Style','Edit',       'String','Tongue',	...
                     'FontSize',fst, 'HorizontalAlignment','left', ...
                     'Position',[lp3 bl2 wgl hl1],   'Enable','Inactive' );
txTNG1	= uicontrol( 'Parent',pANA, 'Style','Edit',       'String',num2str(0,'%d'), ...
                     'FontSize',fst, 'BackgroundColor',bgd, 'HorizontalAlignment','right', ...
                     'Position',[lg1 bl2 wgt hl1],	'Callback', {@set_reflingl_CB}	);
txTNG2	= uicontrol( 'Parent',pANA, 'Style','Edit',       'String',num2str(0,'%d'), ...
                     'FontSize',fst, 'BackgroundColor',bgd, 'HorizontalAlignment','right', ...
                     'Position',[lg2 bl2 wgt hl1],	'Callback', {@set_reflingl_CB}	);
lbDCTTH	= uicontrol( 'Parent',pANA, 'Style','Edit',       'String','DCT LPF',	...
                     'FontSize',fst, 'HorizontalAlignment','left', ...
                     'Position',[lp3 bo1 wlt hl1],   'Enable','Inactive' );
txDCTTH	= uicontrol( 'Parent',pANA, 'Style','Edit',       'String',num2str(0.60,'%0.2f'), ...
                     'FontSize',fst, 'BackgroundColor',bgd, 'HorizontalAlignment','right', ...
                     'Position',[lt3 bo1 wt1 hl1],	'Callback', {@set_bndsmthresh_CB}	);
                 
%	tract analysis controls: col 4
pbAFFRM	= uicontrol( 'Parent',pANA, 'Style','pushbutton', 'String','Correct Video', ...
                     'FontSize',fss, 'ForegroundColor',colgrn, ...
                     'Position',[lp4 bp6 wp1 hp1],	'Callback', {@pixintensity_correct} );
pbAFFRM	= uicontrol( 'Parent',pANA, 'Style','pushbutton', 'String','Filter Video', ...
                     'FontSize',fss, 'ForegroundColor',colgrn, ...
                     'Position',[lp4 bp5 wp1 hp1],	'Callback', {@lpf_video_CB} );
lbVLPF	= uicontrol( 'Parent',pANA, 'Style','Edit',       'String','LPF Param',	...
                     'FontSize',fst, 'HorizontalAlignment','left', ...
                     'Position',[lp4 bl2 wlt hl1],   'Enable','Inactive' );
txVLPF	= uicontrol( 'Parent',pANA, 'Style','Edit',       'String',num2str(50,'%d'), ...
                     'FontSize',fst, 'BackgroundColor',bgd, 'HorizontalAlignment','right', ...
                     'Position',[lt4 bl2 wt1 hl1],	'Callback', {@set_dctthresh_CB} );
lbVLPF	= uicontrol( 'Parent',pANA, 'Style','Edit',       'String','MAF Wind',	...
                     'FontSize',fst, 'HorizontalAlignment','left', ...
                     'Position',[lp4 bo1 wlt hl1],   'Enable','Inactive' );
txVLPF	= uicontrol( 'Parent',pANA, 'Style','Edit',       'String',num2str(50,'%d'), ...
                     'FontSize',fst, 'BackgroundColor',bgd, 'HorizontalAlignment','right', ...
                     'Position',[lt4 bo1 wt1 hl1],	'Callback', {@set_dctthresh_CB} );

                 
% correlation analysis controls
bc3 = bl2+hl1+dy1;
bc4 = bc3+hp2+dy1/2;
bc5 = bc4+hp2+dy1/2;
wpe	= wpexp-2*le1-dx1;
wte = 0.435*wpe;	wle = wpe-wte;  lte = le1+wle;
wp2 = 0.50*wpe;     lp2 = le1+wp2;
wla = 0.35*wpe;     wtb = wpe-wcb-wla;
lcb = le1+wla;      ltb = lcb+wcb;
pbANCOR	= uicontrol( 'Parent',pCOR, 'Style','pushbutton', 'String','Trace', ...
                     'FontSize',fss, 'ForegroundColor','blue', 'HorizontalAlignment','left', ...
                     'Position',[le1 bc5 wpe hp2],	'Callback', {@init_correl_CB,0} );
pbANDUR	= uicontrol( 'Parent',pCOR, 'Style','pushbutton', 'String','CDur', ...
                     'FontSize',fss, 'ForegroundColor','blue', 'HorizontalAlignment','left', ...
                     'Position',[le1 bc4 wp2 hp2],	'Callback', {@init_correl_CB,1} );
pbIANAL	= uicontrol( 'Parent',pCOR, 'Style','pushbutton', 'String','CInt', ...
                     'FontSize',fss, 'ForegroundColor','blue', 'HorizontalAlignment','left', ...
                     'Position',[lp2 bc4 wp2 hp2],	'Callback', {@intensity_anal} );
pbANCOF	= uicontrol( 'Parent',pCOR, 'Style','pushbutton', 'String','Frequency', ...
                     'FontSize',fss, 'ForegroundColor','blue', 'HorizontalAlignment','left', ...
                     'Position',[le1 bc3 wpe hp2],	'Callback', {@analyze_cofreq_CB} );
lbCORTR	= uicontrol( 'Parent',pCOR, 'Style','Edit', 'String','nTr', 'Enable','Inactive',	 ...
                     'FontSize',fst, 'ForegroundColor','black', 'HorizontalAlignment','left', ...
                     'Position',[le1 bl2 wt1 hl1] );
cbCORTR	= uicontrol( 'Parent',pCOR, 'Style','checkbox', 'Value',0, ...
                     'Position',[lcb bl2 wcb hcb],	'Callback', {@refresh_GUI} );
txCORTR	= uicontrol( 'Parent',pCOR, 'Style','Edit', 'String',num2str(1,'%1d'), ...
                     'FontSize',fst, 'BackgroundColor',bgd, ...
                     'Position',[ltb bl2 wtb hl1] );
lbCORTR	= uicontrol( 'Parent',pCOR, 'Style','Edit', 'String','Correl', 'Enable','Inactive', ...
                     'FontSize',fst, 'ForegroundColor','black', 'HorizontalAlignment','left', ...
                     'Position',[le1 bo1 wle hl1] );
txCORTH	= uicontrol( 'Parent',pCOR, 'Style','Edit', 'String',num2str(0.8,'%0.1f'), ...
                     'FontSize',fst, 'BackgroundColor',bgd, ...
                     'Position',[lte bo1 wte hl1] );
        
% labelling controls
hpe = (hpclo-2*bo1-2*dy1)/2;
be2 = bo1+hpe+dy1;
wte = 0.625*wpe;	wle = wpe-wte;	lte = le1+wle;
wpb = 0.50*wpe;     lpb = le1+wpb;
pbCLLAB	= uicontrol( 'Parent',pLAB, 'Style','pushbutton', 'String','Clear', ...
                     'FontSize',fss, 'ForegroundColor','yellow', ...
                     'Position',[le1 be2 wpb hpe],	'Callback', {@clear_label_CB} );
pbLILAB	= uicontrol( 'Parent',pLAB, 'Style','pushbutton', 'String','List', ...
                     'FontSize',fss, 'ForegroundColor','yellow', ...
                     'Position',[lpb be2 wpb hpe],	'Callback', {@list_labels_CB} );
lbLABEL	= uicontrol( 'Parent',pLAB, 'Style','Edit', 'String','Lab:', 'Enable','Inactive', ...
                     'FontSize',fss, 'ForegroundColor','black', 'HorizontalAlignment','left', ...
                     'Position',[le1 bo1 wle hpe]	);
txLABEL	= uicontrol( 'Parent',pLAB, 'Style','Edit', 'String','lab1', ...
                     'FontSize',fss, 'BackgroundColor',bgd, ...
                     'Position',[lte bo1 wte hpe],	'Callback', {@goto_frame_CB} );

% export controls
wpb = 0.65*wpe;	 wtb = wpe-wpb;  ltb = le1+wpb;
pbEXPAU	= uicontrol( 'Parent',pEXP, 'Style','pushbutton', 'String','Audio', ...
                     'FontSize',fss, 'ForegroundColor','red', ...
                     'Position',[le1 bp4 wpe hp1],	'Callback', {@export_audio_CB} );
pbEXPVD	= uicontrol( 'Parent',pEXP, 'Style','pushbutton', 'String','Video', ...
                     'FontSize',fss, 'ForegroundColor','red', ...
                     'Position',[le1 bp3 wpe hp1],	'Callback', {@export_video_CB} );
pbEXPFM	= uicontrol( 'Parent',pEXP, 'Style','pushbutton', 'String','Frame', ...
                     'FontSize',fss, 'ForegroundColor','red', ...
                     'Position',[le1 bp2 wpe hp1],  'Callback', {@export_frame_CB} );
pbEXPJS	= uicontrol( 'Parent',pEXP, 'Style','pushbutton', 'String','Seq', ...
                     'FontSize',fss, 'ForegroundColor','red', ...
                     'Position',[le1 bo1 wpb hp1],  'Callback', {@export_frames_CB} );
step = 5;
if isfield(dat,'sg')
    if isfield(dat.sg,'step')
        step = dat.sg.step;
    end
end
txSGSTP	= uicontrol( 'Parent',pEXP, 'Style','Edit', 'String',num2str(step,'%1d'), ...
                     'FontSize',fst, 'BackgroundColor',bgd, ...
                     'Position',[ltb bo1 wtb hp1],	'Callback', {@edit_sgstep_CB} );

% exit controls
pbSAVE	= uicontrol( 'Parent',pCLO, 'Style','pushbutton', 'String','Save', ...
                     'FontSize',fsh,  'ForegroundColor','black', ...
                     'Position',[le1 be2 wpe hpe],	'Callback', {@save_GUI_CB} );
pbEXIT	= uicontrol( 'Parent',pCLO, 'Style','pushbutton', 'String','Exit', ...
                     'FontSize',fsh,  'ForegroundColor','black', ...
                     'Position',[le1 bo1 wpe hpe],	'Callback', {@close_GUI_CB} );


                 
% rescale original video
if ( isstruct(dat.vid.avi) )
	nf = dat.vid.info.NumberOfFrames;
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
        if (dat.vid.invert)
            vid_rs(s).cdata	= flipud(im_);
        else
            vid_rs(s).cdata	= im_;
        end %if
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
    if (avoid_wmark)
        fdims	= [280 280];
        fcenter	= [165 185];
    else
        fdims	= floor( 0.8*dat.vid.dims );
        fcenter	= floor( dat.vid.dims/2 );
    end
    if (dat.status.reload)
        dat.vid.frame	= cropframe;
    else
        dat.vid.frame.size	= fdims;
        dat.vid.frame.pos	= fcenter;
    end
    set(txFOFFx,'String',num2str(dat.vid.frame.pos(1), '%d'));
    set(txFOFFy,'String',num2str(dat.vid.frame.pos(2), '%d'));
    set(txFRMW, 'String',num2str(dat.vid.frame.size(1),'%d'));
    set(txFRMH, 'String',num2str(dat.vid.frame.size(2),'%d'));
end

% reinitialize null data structure arrays now that signal lengths are known
nf = dat.vid.info.NumberOfFrames;

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
if (dat.status.reload)	% reload old data and update GUI params
    new_dat	= dat;
    dat     = old_dat;
    if isempty(dat.aud)
        dat.aud	= new_dat.aud;
    end
    if isempty(dat.sg)
        dat.sg	= new_dat.sg;
    end
    dat.ffn_aud     = new_dat.ffn_aud;
    dat.ffn_vid     = new_dat.ffn_vid;
    dat.ffn_trk     = new_dat.ffn_trk;
    dat.vid         = new_dat.vid;
    dat.vid.grid	= old_dat.vid.grid;
	dat.status.OS	= new_dat.status.OS;
    set(txWINWD,'String',num2str(dat.aud.wwd,'%.1f'));
    set(txFOFFx,'String',num2str(dat.vid.frame.pos(1), '%d'));
    set(txFOFFy,'String',num2str(dat.vid.frame.pos(2), '%d'));
    set(txFRMW, 'String',num2str(dat.vid.frame.size(1),'%d'));
    set(txFRMH, 'String',num2str(dat.vid.frame.size(2),'%d'));
    if ~(isempty(dat.vid.grid.tng))
        set(txTNG1, 'String',num2str(dat.vid.grid.tng(1),  '%d'));
        set(txTNG2, 'String',num2str(dat.vid.grid.tng(end),'%d'));
    end
    if ~(isempty(dat.vid.grid.pal))
        set(txPAL1, 'String',num2str(dat.vid.grid.pal(1),  '%d'));
        set(txPAL2, 'String',num2str(dat.vid.grid.pal(end),'%d'));
    end
    if ~(isempty(dat.vid.grid.pha))
        set(txPHA1, 'String',num2str(dat.vid.grid.pha(1),  '%d'));
        set(txPHA2, 'String',num2str(dat.vid.grid.pha(end),'%d'));
    end
    if (dat.vid.grid.nlines)
        set( txOR1X,'String', num2str(dat.vid.grid.or1(1),'%d') );
        set( txOR1Y,'String', num2str(dat.vid.grid.or1(2),'%d') );
        set( txOR2X,'String', num2str(dat.vid.grid.or2(1),'%d') );
        set( txOR2Y,'String', num2str(dat.vid.grid.or2(2),'%d') );
        set( txGRAD,'String', num2str(dat.vid.grid.rad1,  '%d') );
        set( txGWID,'String', num2str(dat.vid.grid.wid,   '%d') );
    end
end
if (nargin>1)
    dat.vid.grid	= grid;
end
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


% don't assign output until all data structures have been updated
waitfor( pbSAVE,'String','Saved' );
if ( (nargout) && (dat.status.exit) )
    dat	= guidata(hGUI);
end

waitfor( pbEXIT,'String','Exiting' );
if (dat.status.exit)
    disp('   Exiting GUI.' );
    close(hGUI);
end
disp(' ');



%------------------------------------
% subfunctions and callbacks
%------------------------------------
    function update_time( t )
        dat	= guidata(hGUI);
        dt = 1/dat.vid.info.FrameRate;	% min timestep = length of video frame
        if (t >= dt) && (t < dat.aud.dur-dt)
            fv	= floor( t * dat.vid.info.FrameRate );
            dat.time	= t;
            dat.fnum	= fv;
            set( slVINC,  'Value', fv );
            set( txVFRAM, 'String',num2str(fv) );
            set( txVTIME, 'String',num2str(t,'%.3f') );
            guidata(hGUI,dat);
        else
            fprintf('   Time must be specified in range [%0.3f .. %0.3f]\n', dt, dat.aud.dur-dt );
        end
        update_windows();
        refresh_GUI();
    end %update_time()

    function update_frame( f )
        dat	= guidata(hGUI);
        if (f >= 1) && (f <= dat.vid.info.NumberOfFrames)
            t	= f/dat.vid.info.FrameRate;
            dat.time	= t;
            dat.fnum	= f;
            set( slVINC,  'Value', f );
            set( txVFRAM, 'String',num2str(f) );
            set( txVTIME, 'String',num2str(t,'%.3f') );
            guidata(hGUI,dat);
        else 
            fprintf('   Frame number must be specified in range [1..%d]\n', dat.vid.info.NumberOfFrames );
        end 
        update_windows();
        refresh_GUI();
    end %update_frame()

    function update_frame_lims( fint )
        dat	= guidata(hGUI);
        dat.fnum	 = round(mean(fint));
        dat.aud.wwd  = diff(fint/dat.vid.info.FrameRate);
        guidata(hGUI,dat);
        update_frame(dat.fnum);
    end %update_frame_lims()

    function update_windows()
        dat	= guidata(hGUI);
        w	= dat.aud.wwd;
        int	= [dat.time-w/2 dat.time+w/2];
        ss	= floor( int .* dat.aud.Fs );
        if ss(1)<1
            ss(1) = 1;
            int	  = ss/dat.aud.Fs;
            dat.aud.wwd = diff(int);
        end
        if ss(2)>dat.aud.len,
            ss(2) = dat.aud.len;
            int	  = [dat.aud.dur-w dat.aud.dur];
            ss    = floor( int .* dat.aud.Fs );
            dat.aud.wwd = diff(int);
        end
        set( txWINWD, 'String',num2str(dat.aud.wwd,'%.2f') );
        ff	= round( int .* dat.vid.info.FrameRate );
        if ff(1)<1,	ff(1)=1; end;
        if ff(2)>dat.vid.info.NumberOfFrames, ff(2)=dat.vid.info.NumberOfFrames; end;
        dat.aud.tint	= int;
        dat.aud.sint	= ss;
        dat.vid.fint	= ff;
        set( txWFRML, 'String',num2str(ff(1)) );
        set( txWFRMR, 'String',num2str(ff(2)) );
        guidata(hGUI,dat);
    end %update_windows()

    function refresh_GUI( s,e )
        dat	= guidata(hGUI);
        t	= dat.time;
        ta	= dat.aud.tint(1);
        tb	= dat.aud.tint(2);

        % return focus to main GUI in case there are subplots active
        set(0,'CurrentFigure',1);
        
        % update main audio signal display
        cla(fAUD);
        plot( fAUD, dat.aud.tvec, dat.aud.sig );
        xlim( fAUD, [0 dat.aud.dur]);
        set( fAUD, 'XMinorTick','on' );
        yy	= get(fAUD,'YLim');
        line( [t t],yy, 'Color','k','LineStyle','--', 'Parent',fAUD )
        %rectangle(	'Position',[ta,min(yy),tb-ta,diff(yy)], 'FaceColor',colwin, ...
        %            'EdgeColor',colwlim, 'EraseMode','xor', 'Parent',fAUD );
        hWin = patch( [ta,tb,tb,ta], [min(yy),min(yy),max(yy),max(yy)], colwin );
        set( hWin, 'EdgeColor',colwlim, 'FaceAlpha',0.5, 'Parent',fAUD );
 
        % update windowed audio display
        cla(fAUDw);
        plot( fAUDw, dat.aud.tvec, dat.aud.sig );
        xlim( fAUDw, dat.aud.tint );
        line( [t t],yy,	'Color','k','LineStyle','--', 'Parent',fAUDw );
        set( fAUDw, 'Color',colwin );

        % update windowed video display
        set(hGUI,'CurrentAxes',fVID); cla(fVID);
        if get( cbTMODE,'Value' )           % Show tracked video
            im	= dat.tvd.avi(dat.fnum).cdata;
            map	= dat.tvd.avi(dat.fnum).colormap;
            imshow( im,map );
        else
            if get( cbICORR,'Value' )       % Show intensity-corrected video
                im	= dat.cvd.avi(dat.fnum).cdata;
                map	= dat.cvd.avi(dat.fnum).colormap;
                imshow( im,map );
            elseif get( cbSHLPF,'Value' )	% Show dct/lpfiltered video
                im	= dat.lvd.avi(dat.fnum).cdata;
                map	= dat.lvd.avi(dat.fnum).colormap;
                subimage( im,map );
            elseif get( cbVMODE,'Value' )	% Show uninterpolated (unscaled) video
                im	= dat.vid.avi_(dat.fnum).cdata;
                map	= dat.vid.avi_(dat.fnum).colormap;
                subimage( im );
            else                            % Show original video
                im	= dat.vid.avi(dat.fnum).cdata;
                map	= dat.vid.avi(dat.fnum).colormap;
                %subimage( im,map );
                subimage( im );
            end
        end

        % update video image axes
        imdim = size(im);
        imwd  = imdim(1);
        if get( cbSHVAX,'Value' )           % display pixel coords on axes if flagged
            if get( cbSCVAX,'Value' )       % display pixel coords as scaled
                tik	= 0:50:imwd;
                if get( cbVMODE,'Value' )
                    tik	= 0:10:imwd;
                end
                lab	= cellstr(num2str(tik'));
            else                            % display original pixel coords
                tik	= 0:10:dat.vid.info.Width;
                lab = cellstr(num2str(tik'));
                if ~(get( cbVMODE,'Value' ))
                    tik	= tik * dat.vid.scale;
                end
            end %if
            if get( cbSHVMM,'Value' )       % display distances in mm
                tik	= 0:50:(wd_im*dat.status.FOV/100);
                lab = cellstr(num2str(tik'));
                tik	= tik/dat.status.pixres;
                if ~(get( cbVMODE,'Value' ))
                    tik	= tik * dat.vid.scale;
                end
            end %if
            axis on;
            axis([1 imdim(1) 1 imdim(2)]);
            set( gca, 'XTick',tik, 'YTick',tik );
            set( gca, 'XTickLabel',lab, 'YTickLabel',lab );
        else
            axis off;
        end %if

        % update video analysis graphics: frame, grid, data points
        if get( cbSHOWF,'Value' )           % superimpose crop frame
            frm = dat.vid.frame.pts;
            line( [frm.x], [frm.y], 'Color','y','LineStyle','--' );
        end %if
        if get( cbSHOFN,'Value' )           % superimpose crop frame
            fn	 = sprintf( 'f%d', dat.fnum );
            fn_x = round(dat.vid.frame.pts.x(1) + 10);
            fn_y = round(dat.vid.frame.pts.y(3) - 15);
            text( fn_x,fn_y, fn, 'Color','y', 'FontSize',16, 'FontName','liberation sans' );
        end %if
        if get( cbSHOWG,'Value' )           % superimpose analysis grid
            cx = dat.vid.grid.glot(1);      % place cross at glottis
            cy = dat.vid.grid.glot(2);
            line( [cx-dcr cx+dcr-1], [cy cy], [0 0], 'Color','b', 'Linewidth',dat.vid.grid.linewd );
            line( [cx cx], [cy-dcr cy+dcr-1], [0 0], 'Color','b', 'Linewidth',dat.vid.grid.linewd );
            cx = dat.vid.grid.or1(1);      % place cross at lingual origin
            cy = dat.vid.grid.or1(2);
            line( [cx-dcr cx+dcr-1], [cy cy], [0 0], 'Color','b', 'Linewidth',dat.vid.grid.linewd );
            line( [cx cx], [cy-dcr cy+dcr-1], [0 0], 'Color','b', 'Linewidth',dat.vid.grid.linewd );
            cx = dat.vid.grid.or2(1);      % place cross at alveolar origin
            cy = dat.vid.grid.or2(2);
            line( [cx-dcr cx+dcr-1], [cy cy], [0 0], 'Color','b', 'Linewidth',dat.vid.grid.linewd );
            line( [cx cx], [cy-dcr cy+dcr-1], [0 0], 'Color','b', 'Linewidth',dat.vid.grid.linewd );
            for gl = 1:dat.vid.grid.nlines	% draw line across tract at each grid interval
                txt = dat.vid.grid.txt(gl);
                grd = dat.vid.grid.ends(gl);
                if get( cbSHIFT,'Value' )
                    grd.x = grd.x + dat.vt(dat.fnum).dx;
                    grd.y = grd.y + dat.vt(dat.fnum).dy;
                end
                hGL = line( grd.x, grd.y, [0 0], 'Color','b', 'Linewidth',dat.vid.grid.linewd );                
                set( hGL, 'UserData',gl, 'HitTest','on', 'ButtonDownFcn',@gridline_sel );
                if ~mod(gl,5)
                    text( txt.x,txt.y, num2str(gl), 'Color','b', 'Linewidth',dat.vid.grid.linewd );
                end
            end
        end %if
        if get( cbSHOWC,'Value' )           % draw a cross at each point on VT centerline
            if ~isempty( dat.vt(dat.fnum).cen )
                cc = [dat.vt(dat.fnum).cen.xy];
                line( cc(1:2:end),cc(2:2:end), 'LineStyle','none', ...
                      'Marker','+', 'MarkerEdgeColor','y', 'MarkerSize',dat.vid.grid.linewd );
            end
        end %if
        if get( cbSHCGR,'Value' )           % draw nodes and edges of graph constructed over intensity minima
            if ~isempty( dat.vt(dat.fnum).imap )
                nodes = dat.vt(dat.fnum).Mxy;
                edges = dat.vt(dat.fnum).Mwt;
                for r = 1:size(edges,1)
                    for c = 1:size(edges,2)
                        if isfinite(edges(r,c))
                            line( [nodes(r,1) nodes(c,1)],[nodes(r,2) nodes(c,2)], 'Color','y', 'MarkerSize',10 );
                        end
                    end
                end
            end
        end %if
        if get( cbSHCIR,'Value' )           % draw intensity-detection circle
            x	= dat.status.mloc(1);
            y	= dat.status.mloc(2);
            r	= dat.status.radanal;
            %f	= dat.vid.avi(dat.fnum).cdata;
            It	= calc_local_px_intensity( dat.fnum,[x y],r );
            fprintf('   Px = [%d %d]  f%d  Imean = %d\n', x,y,dat.fnum,It );
            rectangle(	'Position',[x-r,y-r,2*r,2*r], 'FaceColor',colwin, ...
                        'Curvature',[1,1], 'EdgeColor',colwlim, ...
                        'EraseMode','xor', 'Parent',fVID );
        end %if
        if get( cbSHOWB,'Value' )           % superimpose tissue boundaries on displayed frame
            if ~isempty( dat.vt(dat.fnum).pts )
                glot  = dat.vid.grid.glot;
                inner = [dat.vt(dat.fnum).pts.lf];
                outer = [dat.vt(dat.fnum).pts.rt];
                ibnd  = [glot inner];
                obnd  = [glot outer ibnd(end-1:end)];
                line( ibnd(1:2:end),ibnd(2:2:end), 'Color','g' );
                line( obnd(1:2:end),obnd(2:2:end), 'Color','g' );
                set( pbFIXPT, 'Visible','On' )
            else
                set( pbFIXPT, 'Visible','Off' )
            end
        else
            set( pbFIXPT, 'Visible','Off' )
        end %if
        if get( cbSHOWS,'Value' )           % superimpose smoothed tissue boundaries on displayed frame
            if ~isempty( dat.vt(dat.fnum).spts )
                bnd	= [dat.vt(dat.fnum).spts.lf];
                line( bnd(1:2:end),bnd(2:2:end), 'Color','w' );
            end
            gl1 = dat.vid.grid.tng(1);
            gl2 = dat.vid.grid.tng(end);
            tx1 = dat.vid.grid.txt(gl1);
            tx2 = dat.vid.grid.txt(gl2);
            text( tx1.x,tx1.y, num2str(gl1), 'Color','w', 'Linewidth',dat.vid.grid.linewd );
            text( tx2.x,tx2.y, num2str(gl2), 'Color','w', 'Linewidth',dat.vid.grid.linewd );
        end %if
        if get( cbSHOWR,'Value' )           % superimpose palatal reference boundary on displayed frame
            if isempty( dat.vid.grid.palate )
                fprintf('   No palatal boundary defined: first place grid and extract AF ...\n');
            else
                pal	 = [dat.vid.grid.palate.pts.rt];
                palx = pal(1:2:end);
                paly = pal(2:2:end);
                if get( cbSHIFT,'Value' )
                    palx  = palx + dat.vt(dat.fnum).dx;
                    paly  = paly + dat.vt(dat.fnum).dy;
                end
                line( palx,paly, 'Color','r', 'LineStyle','-' );
                gl1 = dat.vid.grid.pal(1);
                gl2 = dat.vid.grid.pal(end);
                tx1 = dat.vid.grid.txt(gl1);
                tx2 = dat.vid.grid.txt(gl2);
                text( tx1.x,tx1.y, num2str(gl1), 'Color','r', 'Linewidth',dat.vid.grid.linewd );
                text( tx2.x,tx2.y, num2str(gl2), 'Color','r', 'Linewidth',dat.vid.grid.linewd );
            end
        end %if
        if get( cbSHOWP,'Value' )           % superimpose pharyngeal reference boundary on displayed frame
            if isempty( dat.vid.grid.pharynx )
                fprintf('   No pharyngeal boundary defined: first place grid and extract AF ...\n');
            else
                pha	 = [dat.vid.grid.pharynx.pts.rt];
                phax = pha(1:2:end);
                phay = pha(2:2:end);
                line( phax,phay, 'Color','r', 'LineStyle','-' );
                gl1 = dat.vid.grid.pha(1);
                gl2 = dat.vid.grid.pha(end);
                tx1 = dat.vid.grid.txt(gl1);
                tx2 = dat.vid.grid.txt(gl2);
                text( tx1.x,tx1.y, num2str(gl1), 'Color','r', 'Linewidth',dat.vid.grid.linewd );
                text( tx2.x,tx2.y, num2str(gl2), 'Color','r', 'Linewidth',dat.vid.grid.linewd );
            end
        end %if

        % update spectrogram
        if ~strcmp( dat.ffn_aud,'' )
            
            set(hGUI,'CurrentAxes',fFFT); cla(fFFT);
            sg      = dat.sg.SG;
            ttot	= size(sg,2);
            tmin	= ta/max(dat.sg.t)*ttot;	if (tmin<1), tmin = 1; end;
            tmax	= tb/max(dat.sg.t)*ttot;
            tvals	= ceil( tmin:tmax );
            sg_win	= sg(:,tvals);
            clims	= [min(min(sg_win)) max(max(sg_win))];
            imagesc( dat.aud.tint, [1 dat.sg.Flim], sg_win, clims );
            xlim( dat.aud.tint ); ylim( [1 dat.sg.Flim] );
            set(gca, 'ydir', 'normal', 'box','on');
            colormap(fFFT,cm_sg); %freezeColors;
            xlabel('Time (sec)');
            ylabel('Frequency (Hz)');
            
            % update formant display data
            if (dat.sg.fmt)
                yy	= get(fFFT,'YLim');
                line( [dat.sg.tfmt dat.sg.tfmt],yy,[0 0], 'Color','y' );
            end
            if (dat.sg.ft)
                sig = dat.aud.sig(dat.aud.sint(1):dat.aud.sint(2));
                [ft,tt] = fmt_track( sig, dat.aud.Fs, ta,tb, dat.sg.Flim,dat.sg.nfft, ...
                                     dat.sg.nf, dat.sg.fwin,dat.sg.fstp, 0 );
                hold on;  plot( tt, ft, 'r.', 'MarkerSize',5 );
                dat.sg.ffreqs = ft;
                dat.sg.ftimes = tt;
            end %if            
        end

        % update any companion windows
        hAX	= findobj( findobj( 'Name','Correlation Analysis' ), 'UserData','CATrace' );
        if ~isempty(hAX),
            update_correl_win( hAX );
        end;
        guidata(hGUI,dat);
    end %refresh_GUI()

    % enable/disable formant tracking on spectrogram
    function ft_onoff_CB( s,e )
        dat = guidata(hGUI);
        if (dat.sg.ft)
            dat.sg.ft = 0;
            set( pbFTTOG, 'String','Fmt' );
        else
            dat.sg.ft = 1;
            set( pbFTTOG, 'String','Off' );
        end
        guidata(hGUI,dat);
        refresh_GUI();
    end %ft_onoff_CB()

	% enable/disable signal pre-emphasis on spectrogram
    function preemph_onoff_CB( s,e )
        dat = guidata(hGUI);
        if get( cbPETOG,'Value' )
           dat.sg.pree = 1;
        else
           dat.sg.pree = 0;
        end
        [dat.sg.SG,dat.sg.f,dat.sg.t] = calc_sg( dat.aud.sig, dat.aud.Fs, ...
         0, dat.sg.pree, dat.sg.win, dat.sg.step, dat.sg.Flim, dat.sg.nfft );
        guidata(hGUI,dat);
        refresh_GUI();
    end %preemph_onoff_CB

    % edit frame limits of audio signal display window
    function edit_wflim_CB( s,e )
        dat	= guidata(hGUI);
        flf	= round(str2double( get(txWFRML,'String')) );
        if flf<1, flf=1; end;
        frt	= round(str2double( get(txWFRMR,'String')) );
        if frt<=flf, frt=flf+diff(dat.vid.fint); end;
        if frt>dat.vid.info.NumberOfFrames, frt = dat.vid.info.NumberOfFrames; end;
        update_frame_lims( [flf frt] );
        guidata(hGUI,dat);
    end

    % edit (time) width of audio signal display window edit_wflim_CB
    function edit_winwd_CB( s,e )
        dat	= guidata(hGUI);
        wwd = str2double(get(txWINWD,'String'));
        if (wwd >= min_wwd) && (wwd <= dat.aud.dur)
            dat.aud.wwd = wwd;
            guidata(hGUI,dat);
            update_windows();
            refresh_GUI();
        else 
            fprintf('   Window width must be specified in range [%0.3f .. %0.3f]\n', min_wwd, dat.aud.dur );
        end 
    end %edit_winwd_CB
	% inc/decrement audio signal display window
    function alt_winwd_CB( s,e, dwd )
        dat	= guidata(hGUI);
        wwd = dat.aud.wwd + dwd;
        if (wwd < min_wwd)
            fprintf('   Window width cannot be decreased below %0.2f sec\n', min_wwd );
        elseif (wwd > dat.aud.dur)
            fprintf('   Window width cannot be increased beyond %0.2f sec\n', dat.aud.dur );
        else 
            set(txWINWD,'String',num2str(dat.aud.wwd,'%.1f'));
            dat.aud.wwd = wwd;
            guidata(hGUI,dat);
            update_windows();
            refresh_GUI();
        end 
    end %alt_winwd_CB
    % select entire audio signal
    function aud_selall_CB( s,e )
        dat	= guidata(hGUI);
        dat.aud.wwd = dat.aud.dur;
        set(txWINWD,'String',num2str(dat.aud.wwd,'%.1f'));
        guidata(hGUI,dat);
        update_time(dat.aud.dur/2);
    end %aud_selall_CB

	% edit width of spectrogram time window
    function edit_sgwin_CB( s,e )
        dat	= guidata(hGUI);
        dat.sg.win = str2double(get(txSGWIN,'String'));
        [dat.sg.SG,dat.sg.f,dat.sg.t] = calc_sg( dat.aud.sig, dat.aud.Fs, ...
         0, dat.sg.pree, dat.sg.win, dat.sg.step, dat.sg.Flim, dat.sg.nfft );
        guidata(hGUI,dat);
        update_windows();
        refresh_GUI();
    end
	% edit spectrogram time step
    function edit_sgstep_CB( s,e )
        dat	= guidata(hGUI);
        dat.sg.step = str2double(get(txSGSTP,'String'));
        [dat.sg.SG,dat.sg.f,dat.sg.t] = calc_sg( dat.aud.sig, dat.aud.Fs, ...
         0, dat.sg.pree, dat.sg.win, dat.sg.step, dat.sg.Flim, dat.sg.nfft );
        guidata(hGUI,dat);
        update_windows();
        refresh_GUI();
    end
	% edit spectrogram maximum frequency
    function edit_sgfmax_CB( s,e )
        dat	= guidata(hGUI);
        dat.sg.Flim = str2double(get(txSGFMX,'String'));
        [dat.sg.SG,dat.sg.f,dat.sg.t] = calc_sg( dat.aud.sig, dat.aud.Fs, ...
         0, dat.sg.pree, dat.sg.win, dat.sg.step, dat.sg.Flim, dat.sg.nfft );
        guidata(hGUI,dat);
        update_windows();
        refresh_GUI();
    end
	% edit formant tracking window
    function edit_ftwin_CB( s,e )
        dat	= guidata(hGUI);
        dat.sg.fwin = str2double(get(txFTWIN,'String'));
        [dat.sg.SG,dat.sg.f,dat.sg.t] = calc_sg( dat.aud.sig, dat.aud.Fs, ...
         0, dat.sg.pree, dat.sg.win, dat.sg.step, dat.sg.Flim, dat.sg.nfft );
        guidata(hGUI,dat);
        update_windows();
        refresh_GUI();
    end
	% edit formant tracking step size
    function edit_ftstep_CB( s,e )
        dat	= guidata(hGUI);
        dat.sg.fstp = str2double(get(txFSTEP,'String'));
        [dat.sg.SG,dat.sg.f,dat.sg.t] = calc_sg( dat.aud.sig, dat.aud.Fs, ...
         0, dat.sg.pree, dat.sg.win, dat.sg.step, dat.sg.Flim, dat.sg.nfft );
        guidata(hGUI,dat);
        update_windows();
        refresh_GUI();
    end %edit_ftstep_CB()

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
        fprintf( '   ...Done!\n' );   
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
        fprintf( '   ...Done!\n' );   
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
        % write .eps too
        %fn	= [ fn '.eps' ];
        %ffn	= fullfile( fdn, fn );
        %imwrite( im_, ffn, 'eps' );
        %print(gcf,'new_image.eps','-depsc2','-r300');
    end %export_gca

    % update location after mouse click
    function mouse_down_CB( s,e )
        dat	= guidata(hGUI);
        obj	= gco();
        if ~isempty(obj)
            
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

    % update position after GUI control increment/decrement
    function update_frame_CB( s,e )
        dat = guidata(hGUI);
        fr	= round( get(slVINC,'Value') );
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
        if (fps<1)
            fps_orig = dat.vid.info.FrameRate;
            fprintf('   Select a playback framerate > 1. (Original video rate = %0.2f fps)\n', fps_orig );
        else
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
                    if ((relMatlab < 2014) && ~strcmp(dat.status.OS,'lin'))
                        mov(i) = im2frame( flipud(mov_(i).cdata), gray(256) );
                    else
                        mov(i) = im2frame( mov_(i).cdata, gray(256) );
                    end
                end
            end
            movie( fVID, mov, 1, fps, [0 0 0 0] );
        end 
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
            set( lbSHOFN, 'Visible','on' );
            set( cbSHOFN, 'Visible','on' );
        else
            set( lbSHOFN, 'Visible','off' );
            set( cbSHOFN, 'Value',0 );
            set( cbSHOFN, 'Visible','off' );
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

	% define lingual analysis grid
    function set_grid(reset)
        dat	= guidata(hGUI);
        tof	= 15;	% grid numbering offset from edge of gridlines
        dtx	= dat.vid.grid.gint/2;
        cx	= dat.vid.grid.or1(1);
        cy	= dat.vid.grid.or1(2);
        r	= dat.vid.grid.rad1 - dat.vid.grid.wid/2;
        R	= dat.vid.grid.rad1 + dat.vid.grid.wid/2;
        % horizontal grid through pharynx
        ix	= 0;
        dy	= dat.vid.grid.gint;
        gy	= dat.vid.grid.glot(2);
        for yy = gy-dy:-dy:cy
            ix	= ix+1;
            dat.vid.grid.ends(ix).x	= [cx+r cx+R];
            dat.vid.grid.ends(ix).y	= [yy yy];
            dat.vid.grid.txt(ix).x	= cx+R+tof;
            dat.vid.grid.txt(ix).y	= yy+dtx;
        end;
        % fetch or set lower limit of pharyngeal reference contour
        if (reset)
            phagl1	= floor(ix/2);
            tongl1	= phagl1+1;
            set( txPHA1, 'String',num2str(phagl1) );
            set( txTNG1, 'String',num2str(tongl1) );
        else
            phagl1	= str2double( get(txPHA1,'String') );
            tongl1	= str2double( get(txTNG1,'String') );
        end
        % construct polar grid, extending radially until alveolar ridge
        dy	= yy - cy;
        dth	= (gint - dy) * deg2rad;
        tht	= mypi - dat.vid.grid.theta;      % angle at which rays change from lingual to alveolar
        gth	= dat.vid.grid.gint * deg2rad;	% angle between rays in semipolar grid
        pal1 = 0;
        ix_vel = ix;
        for th = dth:gth:tht
            ix = ix+1;
            x1 = cx+r*cos(th);	y1 = cy-r*sin(th);
            x2 = cx+R*cos(th);	y2 = cy-R*sin(th);
            dat.vid.grid.ends(ix).x	= [x1 x2];
            dat.vid.grid.ends(ix).y	= [y1 y2];
            dat.vid.grid.txt(ix).x	= cx+(R+tof)*cos(th)-dtx;
            dat.vid.grid.txt(ix).y	= cy-(R+tof)*sin(th);
            if (th < 20 * deg2rad)          % initial estimate of end of oro-pharynx:
                ix_vel	= ix;               % ix of last gridline before grid crosses velum
            end
            if ( (x2 < dat.vid.grid.mpal(1)) && ~pal1 )
                pal1	= ix-1;             % ix of first gridline beyond vertical
            end
        end
        % fetch or set upper limit of pharyngeal reference contour
        if (reset)
            phagl2	= ix_vel;
            set( txPHA2, 'String',num2str(phagl2) );
        else
            phagl2	= str2double( get(txPHA2,'String') );
        end
        dat.vid.grid.pha = phagl1:phagl2;
        % fetch or set limits of palatal reference contour
        if (reset)
            palgl1	= pal1;
            palgl2	= ix+4;
            tongl2	= palgl2;
            set( txPAL1, 'String',num2str(palgl1) );
            set( txPAL2, 'String',num2str(palgl2) );
            set( txTNG2, 'String',num2str(tongl2) );
        else
            palgl1	= str2double( get(txPAL1,'String') );
            palgl2	= str2double( get(txPAL2,'String') );
            tongl2	= str2double( get(txTNG2,'String') );
        end
        dat.vid.grid.pal = palgl1:palgl2;
        dat.vid.grid.tng = tongl1:tongl2;
        % reverse polar grid to vertical
        tx	= dat.vid.grid.or2(1);
        ty	= dat.vid.grid.or2(2);
        dth = gth - (tht - th);
        th1 = dat.vid.grid.theta + dth;
        r	= dat.vid.grid.rad2 - dat.vid.grid.wid/2;
        R	= dat.vid.grid.rad2 + dat.vid.grid.wid/2;
        for th = th1:gth:mypi/2
            ix = ix+1;
            x1 = tx+r*cos(th);	y1 = ty+r*sin(th);
            x2 = tx+R*cos(th);	y2 = ty+R*sin(th);
            dat.vid.grid.ends(ix).x	= [x2 x1];
            dat.vid.grid.ends(ix).y	= [y2 y1];
            dat.vid.grid.txt(ix).x	= tx+(R+tof)*cos(th)-dtx;
            dat.vid.grid.txt(ix).y	= ty+(R+tof)*sin(th);
        end
        % vertical grid through oral cavity
        dth = gth - (mypi/2 - th);
        x0	= tx - dth*rad2deg;
        xmn	= dat.vid.grid.mlab - 4*gint;
        if (xmn	< 0), xmn = gint; end;
        ym	= dat.vid.grid.ends(ix).y(1);
        dx	= dat.vid.grid.gint;
        for xx = x0:-dx:xmn;
            ix = ix+1;
            dat.vid.grid.ends(ix).x	= [xx xx];
            dat.vid.grid.ends(ix).y	= [ym ty+r];
            dat.vid.grid.txt(ix).x	= xx-dtx;
            dat.vid.grid.txt(ix).y	= ym+tof;
        end        
        dat.vid.grid.nlines = ix;
        dat.vid.grid.palate.pts(ix).rt	= [];
        dat.vid.grid.pharynx.pts(ix).rt	= [];
        guidata(hGUI,dat);
        reset_vt_data();
    end %set_grid

	% find tissue boundaries on current frame
    function find_tissue_bnds_CB( s,e )
        dat	= guidata(hGUI);
        f	= dat.fnum;
        %if ( get(cbUSESM,'Value') && ~isempty([dat.vt(f).spts.lf]) )
        if get( cbUSESM,'Value' )
            find_tissue_bnds( f );
            smooth_tongue_CB();
        end
        find_tissue_bnds( f );
        refresh_GUI();
    end %find_tissue_bnds_CB

	% find tissue boundaries for specified frame
    function find_tissue_bnds( f )        
        vb = 0;
        dat = guidata(hGUI);
        if ~dat.vid.grid.nlines
            place_grid();
        end
        nmin = get_grid_profiles(f);
        imap = dat.vt(f).imap;
        a	 = dat.vid.grid.diwtcl;	% weight blending factor: distance vs px intensity
        
        % initialize graph for calculating optimal path through tract
        glot = dat.vid.grid.glot;
        if get( cbSHIFT,'Value' )
            glot = glot+[dat.vt(f).dx dat.vt(f).dy];
        end
        Mxy = glot;                 % Mxy: pixel coordinates of intensity minima
        Mwt = Inf*ones(nmin,nmin);	% Mwt: weighted edge matrix constraining search path
        
        % traverse analysis grid and create graph of pixel intensity minima
        ix	= 2;                    % node index
        lnodes	= 1;                % indices of nodes on previous gridline
        ngl	= dat.vt(f).ngl;
        for gl = 1:ngl
            row     = dat.vid.grid.pts(gl).yy;
            col     = dat.vid.grid.pts(gl).xx;
            ixmin	= imap(gl).minix;
            nmin	= length(ixmin);
            minxy	= [col(ixmin); row(ixmin)]';
            if (vb>2), for r=1:size(minxy,1), plot(fVID, minxy(r,1),minxy(r,2),'r+'); end; end;
            Mxy     = [Mxy; minxy];             % insert minima coordinates into graph
            mnodes	= ix:(ix+nmin-1);           % update graph edge weights
            intense	= imap(gl).imap;
            wti     = intense(ixmin);           % weighting due to pixel intensity
            for n = lnodes
                wtd = zeros(1,nmin);
                for nn = 1:nmin                 % calc dist between nodes on adjacent gridlines
                    m = mnodes(nn);
                    wtd(nn) = norm( Mxy(n,:) - Mxy(m,:) );
                end
                wt	= wti.*a + wtd.*(1-a);      % combined wt: internode distance & px intensity
                Mwt(n,mnodes) = wt;             % wt between nodes on adjacent gridlines
                Mwt(mnodes,n) = wt';            % edge matrix is symmetric
            end
            ix      = ix+nmin;
            lnodes	= mnodes;
            dat.vt(f).imap(gl).nodes = lnodes;	% store node numbers assigned to minima
        end
        n   = lnodes(1);
        dat.vt(f).Mxy = Mxy;
        dat.vt(f).Mwt = sparse(Mwt);
        
        % now use graph to find optimal centerline path through tract
        [d pred] = dijkstra(Mwt,1);             % weighted dist glottis -> all other nodes
        path = [];
        while (n ~= 1),
            path = [n path];
            n = pred(n);
        end
        gl = 1;
        vtlen = 0;
        lastn = 1;
        for n = path
            p1	= Mxy(lastn,:);
            p2	= Mxy(n,:);
            mdist = norm( p1 - p2 );
            if (vb>2), fprintf( '   Gridline %d:\tn%d\t[%3.0f %3.0f]\tvtlen = %0.0f\n', gl,n,p1(1),p1(2),vtlen*pixres/dat.vid.scale ); end;
            dat.vt(f).cen(gl).xy = p2;
            dat.vt(f).imap(gl).cnode = n;       % store node selected as vt centre
            nodes = dat.vt(f).imap(gl).nodes;
            cix = find( nodes == n );
            pix = dat.vt(f).imap(gl).minix(cix);
            dat.vt(f).imap(gl).cix = pix;       % ... and pixel index of centre node on gridline
            vtlen = vtlen+mdist;
            lastn = n;
            gl = gl+1;
        end
        %set( cbSHOWC, 'Value',1 );
        if (vb>1), fprintf( '   Gridline %d:\tn%d\t[%3.0f %3.0f]\tvtlen = %0.0f\n', gl,n,p1(1),p1(2),vtlen*pixres/dat.vid.scale ); end;
        
        dat.vt(f).len = vtlen*pixres/dat.vid.scale;
        if (vb>1), fprintf( '\n   Frame %d:\ttract length = %0.1f mm\n', f,dat.vt(f).len ); end;
        
        % traverse tract centreline and find soft tissue boundaries
        %nlab = 3;   % first treat only pre-labial boundaries
        nlab = 0;   % first treat only pre-labial boundaries
        imap = dat.vt(f).imap;
        for gl = 1:ngl-nlab
            n	= imap(gl).cnode;
            nix	= ( imap(gl).nodes==n );
            iix	= imap(gl).minix( nix );                % pixel index of vt centre node on intensity fn
            int = imap(gl).imap( iix );                 % intensity of vt centre node
            rmx	= max( imap(gl).imap(iix+1:end) );      % max intensity in right half of intensity fn
            lmx	= max( imap(gl).imap(1:iix-1) );        % max intensity in left half of intensity fn
            imx	= min( [lmx rmx] );
            itb = int + (imx-int)*dat.vid.grid.ithrsh;  % threshold set to specified ratio of dynamic range
            rit	= imap(gl).imap(iix+1:end);             % right half of intensity function
            rix	= find( rit<itb )+iix;                  % pixel ix of region with intensity < specified threshold
            if isempty(rix)
                [mn,mnix] = min(rit);
                [mx,mxix] = max(rit);
                if (mx<itb)
                    rix	= mxix+iix;
                    dat.vt(f).int(gl).rt = mx;
                else %(mn>itb)
                    rix	= mnix+iix;
                    dat.vt(f).int(gl).rt = mn;
                end
            else
                if isempty(find(diff(rix)>1, 1))
                    rix	= rix( end );
                else
                    rix	= rix( diff(rix)>1 );
                    rix	= rix(1);
                end
                dat.vt(f).int(gl).rt = itb;
            end
            lit	= imap(gl).imap(1:iix-1);
            lix	= find( lit<itb );
            if isempty(lix)
                [mn,mnix] = min(lit);
                [mx,mxix] = max(lit);
                if (mx<itb)
                    lix	= mxix;
                    dat.vt(f).int(gl).lf = mx;
                else %(mn>itb)
                    lix	= mnix;
                    dat.vt(f).int(gl).lf = mn;
                end
            else
                if isempty(find(diff(lix)>1, 1))
                    lix	= lix( 1 );
                else
                    ix	= find( diff(lix)>1 );
                    lix	= lix(ix+1);
                    lix	= lix(end);
                end
                dat.vt(f).int(gl).lf = itb;
            end
            row	= dat.vid.grid.pts(gl).yy;
            col	= dat.vid.grid.pts(gl).xx;
            pts	= [col; row]';
            dat.vt(f).iix(gl).lf = lix;
            dat.vt(f).iix(gl).rt = rix;
            % if using reference tissue boundaries, fetch reference pt if palatal or pharyngeal gridline
            if get( cbUSERF,'Value' )
                if ( ismember(gl,dat.vid.grid.pal)    && ~isempty([dat.vid.grid.palate.pts(gl).rt]) )
                    pt	= dat.vid.grid.palate.pts(gl).rt;
                    dat.vt(f).pts(gl).rt = pt;
                elseif (ismember(gl,dat.vid.grid.pha) && ~isempty([dat.vid.grid.pharynx.pts(gl).rt]))
                    pt	= dat.vid.grid.pharynx.pts(gl).rt;
                    dat.vt(f).pts(gl).rt = pt;
                else
                    dat.vt(f).pts(gl).rt = pts(rix,:);
                end
            else	% otherwise use thresholded pt for this gridline
                dat.vt(f).pts(gl).rt = pts(rix,:);
            end
            % if using smoothed tissue boundaries, fetch reference pt if lingual gridline
            if get( cbUSESM,'Value' )
                if (ismember(gl,dat.vid.grid.tng) && ~isempty([dat.vt(f).spts.lf]))
                    pt	= dat.vt(f).spts(gl).lf;
                    dat.vt(f).pts(gl).lf = pt;
                else
                    dat.vt(f).pts(gl).lf = pts(lix,:);
                end
            else	% otherwise use thresholded pt for this gridline
                dat.vt(f).pts(gl).lf = pts(lix,:);
            end
            rr	= norm( dat.vid.grid.or1 - dat.vt(f).pts(gl).rt );    % distance to outer bnd pt
            rl	= norm( dat.vid.grid.or1 - dat.vt(f).pts(gl).lf );    % distance to inner bnd pt
            if (rl>rr)
                dat.vt(f).pts(gl).lf = dat.vt(f).pts(gl).rt;
            end
            if (vb), fprintf( '   Gridline %d:\ttb intensity %0.0f\n', gl,itb ); end;
        end

        % now find labial tissue boundaries in last nlab gridlines
        for gl = (ngl-nlab+1):ngl
            n	= imap(gl).cnode;
            nix	= ( imap(gl).nodes==n );
            iix	= imap(gl).minix( nix );                % pixel index of vt centre node on intensity fn
            int = imap(gl).imap( iix );                 % intensity of vt centre node
            rmx	= max( imap(gl).imap(iix+1:end) );      % max intensity in right half of intensity fn
            lmx	= max( imap(gl).imap(1:iix-1) );        % max intensity in left half of intensity fn
            imn	= min( [lmx rmx] );
            imx	= max( [lmx rmx] );
            itb = int + (imx-int)*dat.vid.grid.ithrsh;  % threshold set to specified ratio of dynamic range
            if (imn > itb )
                % if lower lip intensity > ithr% of ulip, use same algorithm as for 
                % all pre-labial gridlines, but weight llip threshold independently
                itb = int + (imx-int)*dat.vid.grid.ithrsh;  % threshold set to specified ratio of dynamic range
                rit	= imap(gl).imap(iix+1:end);             % right half of intensity function
                rix	= find( rit<itb )+iix;                  % pixel ix of region with intensity < specified threshold
                if isempty(rix)
                    [mn,mnix] = min(rit);
                    [mx,mxix] = max(rit);
                    if (mx<itb)
                        rix	= mxix+iix;
                        dat.vt(f).int(gl).rt = mx;
                    else %(mn>itb)
                        rix	= mnix+iix;
                        dat.vt(f).int(gl).rt = mn;
                    end
                else
                    if isempty(find(diff(rix)>1, 1))
                        rix	= rix( end );
                    else
                        rix	= rix( diff(rix)>1 );
                        rix	= rix(1);
                    end
                    dat.vt(f).int(gl).rt = itb;
                end
                itb = int + (imn-int)*dat.vid.grid.ithrsh;  % llip threshold: specified ratio of llip dynamic range
                lit	= imap(gl).imap(1:iix-1);
                lix	= find( lit<itb );
                if isempty(lix)
                    [mn,mnix] = min(lit);
                    [mx,mxix] = max(lit);
                    if (mx<itb)
                        lix	= mxix;
                        dat.vt(f).int(gl).lf = mx;
                    else %(mn>itb)
                        lix	= mnix;
                        dat.vt(f).int(gl).lf = mn;
                    end
                else
                    if isempty(find(diff(lix)>1, 1))
                        lix	= lix( 1 );
                    else
                        ix	= find( diff(lix)>1 );
                        lix	= lix(ix+1);
                        lix	= lix(end);
                    end
                    dat.vt(f).int(gl).lf = itb;
                end
                row	= dat.vid.grid.pts(gl).yy;
                col	= dat.vid.grid.pts(gl).xx;
                pts	= [col; row]';
                dat.vt(f).iix(gl).lf = lix;
                dat.vt(f).pts(gl).lf = pts(lix,:);
                dat.vt(f).iix(gl).rt = rix;
                dat.vt(f).pts(gl).rt = pts(rix,:);
            else
                % if lower lip intensity < ithr% of ulip, only select an ulip
                % tissue boundary
                itb = imx*dat.vid.grid.ithrsh;              % threshold set to specified ratio of dynamic range
                rit	= imap(gl).imap(iix+1:end);             % right half of intensity function
                rix	= find( rit<itb )+iix;                  % pixel ix of region with intensity < specified threshold
                if isempty(rix)
                    [mn,mnix] = min(rit);
                    [mx,mxix] = max(rit);
                    if (mx<itb)
                        rix	= mxix+iix;
                        dat.vt(f).int(gl).rt = mx;
                    else %(mn>itb)
                        rix	= mnix+iix;
                        dat.vt(f).int(gl).rt = mn;
                    end
                else
                    if isempty(find(diff(rix)>1, 1))
                        rix	= rix( end );
                    else
                        rix	= rix( diff(rix)>1 );
                        rix	= rix(1);
                    end
                    dat.vt(f).int(gl).rt = itb;
                end
                row	= dat.vid.grid.pts(gl).yy;
                col	= dat.vid.grid.pts(gl).xx;
                pts	= [col; row]';
                dat.vt(f).iix(gl).rt = rix;
                dat.vt(f).pts(gl).rt = pts(rix,:);
                if (imn > 0.4*itb )
                    lit	= imap(gl).imap(1:iix-1);
                    lix	= find( lit == max(lit),1,'first' );
                    dat.vt(f).iix(gl).lf = lix;
                    dat.vt(f).pts(gl).lf = pts(lix,:);
                else
                    dat.vt(f).iix(gl).lf = [];
                    dat.vt(f).pts(gl).lf = dat.vt(f).pts(gl-1).lf;
                end
            end
        end

        set( cbSHOWB, 'Value',1 );
        guidata(hGUI,dat);
        
    end %find_tissue_bnds

	% extract intensity profiles across gridlines for specified image frame
    function nmin = get_grid_profiles( fnum )
        dat = guidata(hGUI);
        vb = 1;     % debugging verbosity
        if     get( cbICORR,'Value' )	% Use intensity-corrected image if flagged
            im	= dat.cvd.avi(fnum).cdata;
        elseif get( cbSHLPF,'Value' )	% Use low-pass filtered image if flagged
            im	= dat.lvd.avi(fnum).cdata;
        else
            im	= dat.vid.avi(fnum).cdata;
        end
        % traverse whole grid and find all minima
        ngl	= dat.vid.grid.nlines;
        for gl = 1:ngl
            grd	 = dat.vid.grid.ends(gl);
            %npts = floor( norm( [grd.x(1) grd.y(1)]-[grd.x(2) grd.y(2)] ) );
            npts = dat.vid.grid.wid;
            col	 = round( linspace( grd.x(1),grd.x(2),npts ) );
            row	 = round( linspace( grd.y(1),grd.y(2),npts ) );
            imap = zeros(1,npts);
            ix	 = 1:npts;
            for pix = ix
                imap(pix) = im( row(pix), col(pix) );
            end
            dat.vid.grid.pts(gl).xx = col;
            dat.vid.grid.pts(gl).yy = row;
            dat.vt(fnum).imap(gl).imap = imap;
            %imap_  = spline(ix,imap, ix);
            [ymax,imax,ymin,imin] = extrema( imap );
            imin = imin( imin>1 & imin<npts );  % eliminate minima on gridline endpoints
            dat.vt(fnum).imap(gl).minix	 = imin;
            dat.vt(fnum).imap(gl).minima = imap(imin);
            dat.vt(fnum).imap(gl).maxix  = imax;
            dat.vt(fnum).imap(gl).maxima = ymax;
            % find inflection pts on difference function
            dint  = diff(imap);
            dat.vt(fnum).imap(gl).dint      = dint;
            [ymax,imax,ymin,imin] = extrema( dint );
            imin = imin( imin>1 & imin<length(dint) );	% eliminate endpoint extrema
            imax = imax( imax>1 & imax<length(dint) );
            dat.vt(fnum).imap(gl).dminix	= imin;
            dat.vt(fnum).imap(gl).dminima	= dint(imin);
            dat.vt(fnum).imap(gl).dmaxix	= imax;
            dat.vt(fnum).imap(gl).dmaxima	= dint(imax);
        end
        % estimate mean intensities of background and soft tissue pixels
        % and use these values to calculate tissue transition threshold
        int_bg	= mean([dat.vt(fnum).imap(ngl-3:ngl).maxima]);
        int_st	= mean([dat.vt(fnum).imap(1:ngl-6).maxima]);
        int_th	= 0.7*(int_st-int_bg);
        if (vb>1), fprintf( '   Mean maximum background pixel intensity = %0.0f\n', int_bg ); end;
        if (vb>1), fprintf( '   Mean maximum soft tissue intensity      = %0.0f\n', int_st ); end;
        if (vb>1), fprintf( '   Labial tissue intensity theshold set to:  %0.0f\n', int_th ); end;
        % now backtrack through oral region of grid to locate lips
        for gl = dat.vid.grid.nlines:-1:1
            mpi	= mean([dat.vt(fnum).imap(gl).maxima]);
            if (vb>2), fprintf( '   Gridline %d:\tmean tissue intensity = %0.0f\n', gl,mpi ); end;
            if mpi > int_th     % 1st gline where mean max intensity > background noise
                break;
            end
        end
        dat.vt(fnum).ngl = gl;
        if (vb>1), fprintf( '   Labial limit detected at gridline %d\n\n', gl ); end;
        guidata(hGUI,dat);
        nmin = length( [dat.vt(fnum).imap.minix] );
    end %get_grid_profiles

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
        if get( cbSHCGR,'Value' )       % if intensity minima graph displayed, edit that
            edit_graph( dat.fnum );
        elseif get( cbSHOWC,'Value' )	% if tract centerline displayed, edit that
            edit_cline( dat.fnum );
        else                            % otherwise, edit boundaries
            set( cbSHOWB, 'Value',0 );
            refresh_GUI();
            edit_boundaries( dat.fnum );
        end %if
    end %edit_boundaries_CB

    % manually edit centerline for frame f
    function edit_cline( f )
        dat	= guidata(hGUI);
        cc	= [dat.vt(dat.fnum).cen.xy];
        cc	= reshape(cc,2,[])'
        set(hGUI,'CurrentAxes',fVID)
        xl	= get(fVID, 'XLim');
        % move pt. closest to mouse click to pt. selected by user
        [x,y]	= ginput(1);
        pt      = ceil([x y]);
        while ( (pt > xl(1)) & (pt < xl(2)) )
            pd	= Inf;
            % find node closest to pt selected
            for ci	= 1:length(cc)
                cxy	= cc(ci,:);
                dp	= abs(norm(pt-cxy));
                if (dp<pd),     % find gridline which pt lies on
                    cix	= ci;
                    pd	= dp;
                end
            end
            cxy	= cc(cix,:);
            dat.vt(dat.fnum).cen(cix).xy = pt;	% update point on centerline 
            %set( hBndI, 'Visible','Off');
            %hOld  = line( oldbi(gls*2-1), oldbi(gls*2),  'Color','g', 'LineStyle',':' );
            %hBndI = line( newpts(gls*2-1),newpts(gls*2), 'Color','r', 'LineStyle','-' );
            fprintf('   Frame %d gline %d:\t[%0.0f,%0.0f] -> [%0.0f,%0.0f]\n', f,cix, cxy(1),cxy(2), pt(1),pt(2) );
            [x,y]	= ginput(1);
            pt      = ceil([x y]);
            %set( hOld, 'Visible','Off');
        end
        dat.vt(f).mod = 1;      % flag frame as having been edited
        guidata(hGUI,dat);
        refresh_GUI();
    end %edit_cline

    % manually edit graph for frame f
    function edit_graph( f )
        dat	= guidata(hGUI);
        nodes = dat.vt(dat.fnum).Mxy;
        edges = dat.vt(dat.fnum).Mwt;
        set(hGUI,'CurrentAxes',fVID)
        xl	= get(fVID, 'XLim');
        % move pt. closest to mouse click to pt. selected by user
        [x,y]	= ginput(1);
        pt      = ceil([x y]);
        while ( (pt > xl(1)) & (pt < xl(2)) )
            pd	= Inf;
            % find node closest to pt selected
            for ni = 1:length(nodes)
                n	= nodes(ni,:);
                dp	= abs(norm(pt-n));
                if (dp<pd),     % find gridline which pt lies on
                    nix	= ni;
                    pd	= dp;
                end
            end
            n	= nodes(nix,:);
            dat.vt(f).Mxy(nix,:) = [];	% delete node
            dat.vt(f).Mwt(nix,:) = [];	% delete all edges leading to deleted node
            dat.vt(f).Mwt(:,nix) = [];	% delete all edges leading from deleted node
            %set( hBndI, 'Visible','Off');
            %hOld  = line( oldbi(gls*2-1), oldbi(gls*2),  'Color','g', 'LineStyle',':' );
            %hBndI = line( newpts(gls*2-1),newpts(gls*2), 'Color','r', 'LineStyle','-' );
            fprintf('   Frame %d Node %d:\t[%0.0f,%0.0f] -> []\n', f,nix, n(1),n(2) );
            [x,y]	= ginput(1);
            pt      = ceil([x y]);
            %set( hOld, 'Visible','Off');
        end
        dat.vt(f).mod = 1;      % flag frame as having been edited
        guidata(hGUI,dat);
        refresh_GUI();
    end %edit_graph

    % manually edit boundaries for frame f
    function edit_boundaries( f )
        dat	= guidata(hGUI);
        if ( isempty( dat.vt(f).imap ) )
            find_tissue_bnds( f );
        end
        set(hGUI,'CurrentAxes',fVID)
        xl	= get(fVID, 'XLim');
        % show current tissue boundaries
        gls     = 1:dat.vt(f).ngl;
        pts     = dat.vt(f).pts;
        oldbi	= [pts.lf];
        oldbo	= [pts.rt];
        glot	= dat.vid.grid.glot;
        ibnd	= [glot oldbi];
        obnd	= [glot oldbo ibnd(end-1:end)];
        hBndI	= line( ibnd(gls*2-1),ibnd(gls*2), 'Color','r', 'LineStyle','-' );
        hBndO	= line( obnd(gls*2-1),obnd(gls*2), 'Color','r', 'LineStyle','-' );
        % move pt. closest to mouse click to pt. selected by user
        [x,y,mb] = ginput(1);
        inner	 = (mb == 1);   % left mouse button: inner bnd; rt button: outer bnd
        pt       = ceil([x y]);
        while ( (pt > xl(1)) & (pt < xl(2)) )
            pd	= Inf;
            for gl = gls        % find boundary pt closest to pt selected
                gg	= dat.vid.grid.ends(gl);
                p1	= [gg.x(1) gg.y(1)];
                p2	= [gg.x(2) gg.y(2)];
            	dp	= min(abs(norm(pt-p1)),abs(norm(pt-p2)));
                if (dp<pd),     % find gridline which pt lies on
                    gix	= gl;
                    pd	= dp;
                end;
            end
            % get angle of closest gridline to point selected
            gl	= dat.vid.grid.ends(gix);
            x1	= gl.x(1);	y1	= gl.y(1);                      % [x1 y1]: inner gridline endpt
            x2	= gl.x(2);	y2	= gl.y(2);                      % [x2 y2]: outer gridline endpt
            tht	= real(acos( (x2-x1)/norm([x2 y2]-[x1 y1]) ));	% angle of gridline
            d	= norm(pt-[x1 y1]);                             % d(inner gl endpt-> new bnd pt)
            if (inner)	% find closest pt on inner boundary
                old	= pts(gix).lf;
                dx	= d*cos(tht);	dy = d*sin(tht);
                x	= x1 + dx;      y  = y1 - dy;
                dat.vt(f).pts(gix).lf = round([x y]);
                set( hBndI, 'Visible','Off');
                newpts	= [dat.vt(f).pts.lf];
                hOld  = line( oldbi(gls*2-1), oldbi(gls*2),  'Color','g', 'LineStyle',':' );
                hBndI = line( newpts(gls*2-1),newpts(gls*2), 'Color','r', 'LineStyle','-' );
            else        % find closest pt on outer boundary
                old	= pts(gix).rt;
                dx	= d*cos(tht);	dy = d*sin(tht);
                x	= x1 + dx;      y  = y1 - dy;
                dat.vt(f).pts(gix).rt = round([x y]);
                set( hBndO, 'Visible','Off');
                newpts	= [dat.vt(f).pts.rt];
                hOld  = line( oldbo(gls*2-1), oldbo(gls*2),  'Color','g', 'LineStyle',':' );
                hBndO = line( newpts(gls*2-1),newpts(gls*2), 'Color','r', 'LineStyle','-' );
            end;
            fprintf('   Frame %d GLine %d:\t[%0.0f,%0.0f] -> [%0.0f,%0.0f]\n', f,gix, old(1),old(2), pt(1),pt(2) );
            [x,y,mb] = ginput(1);
            inner	 = (mb == 1);
            pt       = ceil([x y]);
            set( hOld, 'Visible','Off');
        end
        dat.vt(f).mod = 1;      % flag frame as having been edited
        set( cbSHOWB, 'Value',1 );
        guidata(hGUI,dat);
        refresh_GUI();
    end %edit_boundaries

    % get palate trace callback
    function edit_palate_CB( s,e )
        dat	= guidata(hGUI);
        set( cbSHOWR, 'Value',0 );  % turn off current boundary displays
        set( cbSHOWB, 'Value',0 );
        refresh_GUI();
        edit_palate( dat.fnum );
    end
    % get palate trace from reference frame and allow editing with mouse
    function edit_palate( f )
        dat	= guidata(hGUI);
        gls	= dat.vid.grid.pal;     % extent of palate (range of gridlines)
        if ( isempty(dat.vt(f).imap) )
            find_tissue_bnds( f );
        end
        pts	= dat.vt(f).pts;
        dat.vid.grid.palfnum = f;
        if ( get( cbUSERF,'Value' ) && ~isempty([dat.vid.grid.palate.pts.rt]) )
            % use reference boundary if flagged and ref pts exist
            palpts	= dat.vid.grid.palate.pts;
            for pgl = dat.vid.grid.pal
                pts(pgl).rt	= palpts(pgl).rt;
            end
        end
        oldpts	= round([pts(gls(1):gls(end)).rt]);
        set(hGUI,'CurrentAxes',fVID)
        xl	= get(fVID, 'XLim');
        % show current palatal tissue boundaries
        hPal	= line( oldpts(1:2:end),oldpts(2:2:end), 'Color','r', 'LineStyle','-' );
        % move pt. closest to mouse click to pt. selected by user
        pt	= ceil( ginput(1) );
        while ( (pt > xl(1)) & (pt < xl(2)) )
            pd	= Inf;
            for gl = gls
                gg	= dat.vid.grid.ends(gl);
                p1	= [gg.x(1) gg.y(1)];
                p2	= [gg.x(2) gg.y(2)];
            	dp	= min(abs(norm(pt-p1)),abs(norm(pt-p2)));
                if (dp<pd),     % find closest point to point selected
                    pd	= dp;
                    gix	= gl;
                end;
            end
            % get angle of closest gridline to point selected
            gl	= dat.vid.grid.ends(gix);
            x1	= gl.x(1);	y1	= gl.y(1);                      % [x1 y1]: inner gridline endpt
            x2	= gl.x(2);	y2	= gl.y(2);                      % [x2 y2]: outer gridline endpt
            tht	= real(acos( (x2-x1)/norm([x2 y2]-[x1 y1]) ));	% angle of gridline
            d	= norm(pt-[x1 y1]);                             % d(inner gl endpt-> new bnd pt)
            % project new palatal pt onto closest gridline
            old	= pts(gix).rt;
            dx	= d*cos(tht);	dy = d*sin(tht);
            x	= x1 + dx;      y  = y1 - dy;
            dat.vt(f).pts(gix).rt = [x y];
            set( hPal, 'Visible','Off');
            newpts	= [dat.vt(f).pts.rt];
            hOld = line( oldpts(1:2:end),oldpts(2:2:end), 'Color','g', 'LineStyle',':' );
            hPal = line( newpts(gls*2-1),newpts(gls*2),   'Color','r', 'LineStyle','-' );
            fprintf('   Gridline %d:\t[%0.0f,%0.0f] -> [%0.0f,%0.0f]\n', gix, old(1),old(2), pt(1),pt(2) );
            pt = ceil( ginput(1) );
            set( hOld, 'Visible','Off');
        end
        % store reference palate for access from all frames
        for gl = gls
            dat.vid.grid.palate.pts(gl).rt = dat.vt(dat.fnum).pts(gl).rt;
        end
        set( cbSHOWR, 'Value',1 );
        guidata(hGUI,dat);
        refresh_GUI();
    end %edit_palate

    % get pharyngeal trace callback
    function edit_pharynx_CB( s,e )
        dat	= guidata(hGUI);
        set( cbSHOWP, 'Value',0 );  % turn off current boundary displays
        set( cbSHOWB, 'Value',0 );
        refresh_GUI();
        edit_pharynx( dat.fnum );
    end
    % get pharyngeal trace and allow editing with mouse
    function edit_pharynx( f )
        dat	= guidata(hGUI);
        gls	= dat.vid.grid.pha;     % extent of pharynx (range of gridlines)
        if ( isempty(dat.vt(f).imap) )
            find_tissue_bnds( f );
        end
        pts	= dat.vt(f).pts;
        dat.vid.grid.phafnum = f;
        if ( get( cbUSERF,'Value' ) && ~isempty([dat.vid.grid.pharynx.pts.rt]) )
            % use reference boundary if flagged and ref pts exist
            phapts	= dat.vid.grid.pharynx.pts;
            for pgl = dat.vid.grid.pha
                pts(pgl).rt	= phapts(pgl).rt;
            end
        end
        oldpts	= [pts(gls(1):gls(end)).rt];
        set(hGUI,'CurrentAxes',fVID)
        xl	= get(fVID, 'XLim');
        % show current pharyngeal tissue boundaries
        hPha	= line( oldpts(1:2:end),oldpts(2:2:end), 'Color','r', 'LineStyle','-' );
        % move pt. closest to mouse click to pt. selected by user
        pt	= ceil( ginput(1) );
        while ( (pt > xl(1)) & (pt < xl(2)) )
            pd	= Inf;
            for gl = gls
                gg	= dat.vid.grid.ends(gl);
                p1	= [gg.x(1) gg.y(1)];
                p2	= [gg.x(2) gg.y(2)];
            	dp	= min(abs(norm(pt-p1)),abs(norm(pt-p2)));
                if (dp<pd),     % find closest point to point selected
                    pd	= dp;
                    gix	= gl;
                end;
            end
            % get angle of closest gridline to point selected
            gl	= dat.vid.grid.ends(gix);
            x1	= gl.x(1);	y1	= gl.y(1);                      % [x1 y1]: inner gridline endpt
            x2	= gl.x(2);	y2	= gl.y(2);                      % [x2 y2]: outer gridline endpt
            tht	= real(acos( (x2-x1)/norm([x2 y2]-[x1 y1]) ));	% angle of gridline
            d	= norm(pt-[x1 y1]);                             % d(inner gl endpt-> new bnd pt)
            % project new palatal pt onto closest gridline
            old	= pts(gix).rt;
            dx	= d*cos(tht);	dy = d*sin(tht);
            x	= x1 + dx;      y  = y1 - dy;
            dat.vt(f).pts(gix).rt = [x y];
            set( hPha, 'Visible','Off');
            newpts	= [dat.vt(f).pts.rt];
            hOld = line( oldpts(1:2:end),oldpts(2:2:end), 'Color','g', 'LineStyle',':' );
            hPha = line( newpts(gls*2-1),newpts(gls*2),   'Color','r', 'LineStyle','-' );
            fprintf('   Gridline %d:\t[%0.0f,%0.0f] -> [%0.0f,%0.0f]\n', gix, old(1),old(2), pt(1),pt(2) );
            pt = ceil( ginput(1) );
            set( hOld, 'Visible','Off');
        end
        % store reference pharynx for access from all frames
        for gl = gls
            dat.vid.grid.pharynx.pts(gl).rt = dat.vt(f).pts(gl).rt;
        end
        set( cbSHOWP, 'Value',1 );
        guidata(hGUI,dat);
        refresh_GUI();
    end %edit_pharynx

	% toggle averaged pharynx mode
    function mean_pharynx_CB( s,e )
        dat = guidata(hGUI);
        set( cbUSESM, 'Value',1 );
        refresh_GUI();
        mean_pharynx();
    end
	% time-filter video using moving average filter
    function  mean_pharynx( s,e )
        dat = guidata(hGUI);
        pha = dat.vid.grid.pha;
        fsp	= str2double(get(txFSTEP,'String'));
        nf	= dat.vid.info.NumberOfFrames;
        fix = 1:fsp:nf;
        fprintf('\n   Averaging theshold values for pharyngeal gridlines %d to %d', pha(1),pha(end) );
        fprintf('\n   Locating every %dth pharyngeal boundary ... ', fsp );
        % show progress bar
        hWB	= waitbar(0,'Averaging theshold values for pharyngeal gridlines.');
        jf	= get( get(hWB,'JavaFrame'), 'FigurePanelContainer');
        jf.getComponent(0).getRootPane.getTopLevelAncestor.setAlwaysOnTop(1);
        for f = fix
            extract_AF_frame( f, 0 );
            waitbar(f/nf);	% report progress
        end
        waitbar(1);	close(hWB);
        for p = pha
            allphar(p).xy = zeros(length(fix),2);
        end
        i = 1;
        for f = fix
            pts	= dat.vt(f).pts;
            for p = pha
                allphar(p).xy(i,:)	= pts(p).rt;
            end
            i = i+1;
        end
        for p = pha
            dat.vid.grid.pharynx.pts(p).rt	= mean(allphar(p).xy);
        end
        fprintf('done.\n\n');
        dat.vid.grid.phafnum = nf+1;    % flag meaned pharynx by setting reference frame > no. frames
        guidata(hGUI,dat);
        set( cbSHOWP, 'Value',1 );
        refresh_GUI();
    end

    % shift tract boundaries callback
    function shift_tract_bnds_CB( s,e )
        dat	= guidata(hGUI);
        if (dat.vid.grid.palfnum)
            if get( cbSHIFT,'Value' )
                if ~(dat.vid.grid.shifted)
                    shift_tract_bnds();
                    dat.vid.grid.shifted = 1;	% flag that shift algorithm has already been applied
                    guidata(hGUI,dat);
                else
                    fprintf('\n   Shifting frame %d [%d %d] px to align with palatal reference frame (%d).\n\n', dat.fnum, dat.vt(dat.fnum).dx, dat.vt(dat.fnum).dy, dat.vid.grid.palfnum);
                end
                refresh_GUI();
            end
        else
            fprintf('\n   No reference frame defined: adjust palate and try again ...\n');
            set( cbSHIFT, 'Value',0 );
            guidata(hGUI,dat);
        end
    end %shift_tract_bnds_CB

    % shift tract boundaries callback
    function shift_tract_bnds( s,e )
        dat	= guidata(hGUI);
        ref	= dat.vid.grid.palfnum;
        fprintf('\n   Shifting boundaries to align with reference frame (%d) ... ', ref);
        M1	= dat.vid.avi(ref).cdata;
        ff  = 39;	% interframe spacing for checking overall variance
        nf	= dat.vid.info.NumberOfFrames;
        stp	= floor(nf/ff);
        j	= 1;
        for i = 1:stp:nf
            M(j,:,:) = dat.vid.avi(i).cdata;
            j = j+1;
        end
        M       = double(M);
        VARM	= var(M,1);
        tmp(:,:)= VARM(1,:,:);
        VARM	= tmp;
        [xs1,ys1,xs2,ys2] = vtShift(M,VARM);
        for i = 1:dat.vid.info.NumberOfFrames
            M2              = dat.vid.avi(i).cdata;
            [xshift,yshift]	= find_shift(M1,M2,VARM,xs1,ys1,xs2,ys2);
            dat.vt(i).dx	= xshift;
            dat.vt(i).dy	= yshift;
        end
        fprintf('\n   Done: shifted this frame [%d %d] pixels.\n\n', dat.vt(dat.fnum).dx, dat.vt(dat.fnum).dy);
        guidata(hGUI,dat);
    end %shift_tract_bnds

    % smooth tract boundaries callback
    function smooth_tongue_CB( s,e)
        dat	= guidata(hGUI);
        if ( isempty( dat.vt(dat.fnum).imap ) )
            find_tissue_bnds(dat.fnum);
        end
        if ~isempty(dat.vid.grid.tng)
            glines	= dat.vid.grid.tng;
            dctthr	= str2double(get(txDCTTH,'String'));
            smooth_ibnd( dat.fnum, dctthr, glines, 1 );
            set( cbSHOWS, 'Value',1 );
            %set( cbSHOWB, 'Value',0 );
        end
        refresh_GUI();
    end %smooth_tongue_CB

    % smooth inner tract boundary
    function smooth_ibnd( f, dctthr, glines, vb )
        dat	= guidata(hGUI);
        pts	= dat.vt(f).pts;
        ngl	= length(glines);
        Di	= zeros(1,ngl);
        gix	= 1;
        for g = glines
            gline	= dat.vid.grid.ends(g);
            orig	= [gline.x(1) gline.y(1)];
            Di(gix)	= norm( pts(g).lf - orig) * pixres/dat.vid.scale;
            if (vb>1), fprintf( '   Gridline %d:\td = %3.0f\n', g, Di(gix) ); end;
            gix     = gix+1;
        end
        bnd_ = dct(Di);
        % Find how many DCT coefficients represent dctthr% of the energy in a sequence
        dctthr = 1-(dctthr/100);
        orderedcomponents = 1;
        if (orderedcomponents)
            [XX,ix] = sort(abs(bnd_)); ix = fliplr(ix);
            i = 1;
            while (norm([bnd_(ix(1:i)) zeros(1,100-i)])/norm(bnd_)<dctthr)
                i = i+1;
            end
            thrsh = XX(end-i-1);
            fprintf( '   Using %d most significant DCT components (of %d)\n', i,max(ix) );
            if (vb>2),
                figure;  plot( bnd_,'b' );
                figure;	 plot( XX,'b' );
                hold on; plot( thrsh*ones(1,ngl),'r:' );
            end
        else
            XX = abs(bnd_);
            dct_lpf = cumsum(XX)/sum(XX);
            n_comps = length( find( dct_lpf <= dctthr) );
            thrsh	= XX(n_comps);
        end
        bnd_( abs(bnd_)<thrsh ) = 0;
        bnd  = idct(bnd_);
        if (vb>1),
            hTB = figure;
            fig_pos = get(hGUI,'Position') + [0 gui_ht+29 0 0];
            fig_pos(3) = gui_wd/2;  fig_pos(4) = gui_wd/9;
            set( hTB, 'Name',['VT Tissue Boundary Contours: Frame ' num2str(f)], 'Position', fig_pos );
            set( hTB, 'ToolBar','none', 'MenuBar','none', 'WindowButtonDownFcn',{@TB_mouse_down_CB, f} );
            plot( glines,Di, 'b-'  ); hold on;
            plot( glines,bnd,'k--' );
            axis tight;	yy=ylim; ylim([yy(1)-5 yy(2)+5]);
            xlabel('Gridline');	ylabel('Dist from Origin (mm)');
        end
        fprintf( '\n' );
        gix	= 1;
        for g = glines
            gl	= dat.vid.grid.ends(g);
            x1	= gl.x(1);      y1	= gl.y(1);                  % p1 (origin): innermost grid pt
            xo	= pts(g).rt(1); yo	= pts(g).rt(2);             % original outer boundary pt
            th	= real(acos( (xo-x1)/norm([xo yo]-[x1 y1]) ));	% angle of gridline
            di	= bnd(gix) / pixres*dat.vid.scale;              % d -> inner bnd
            dix	= di*cos(th);	diy	= di*sin(th);
            xi_	= x1 + dix;     yi_	= y1 - diy;
            dat.vt(f).spts(g).lf = [xi_ yi_];
            if (vb>1), fprintf( '   Gridline %d:\ttheta = %3.0f:\tpi = [%3.0f %3.0f]\n', g,th*rad2deg,xi_,yi_ ); end;
            gix	= gix+1;
        end
        guidata(hGUI,dat);
    end %smooth_ibnd

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
        for f = 1:dat.vid.info.NumberOfFrames
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
                dur	= (t2-t1)/dat.vid.info.FrameRate;
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
            %fps	= dat.vid.info.FrameRate;
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
        t1	= pt(1,1)/dat.vid.info.FrameRate + dat.corr.tl;
        function corr_mouse_up_CB( s,e )
            pt2	= get( hAX,'CurrentPoint' );
            t2	= pt2(1,1)/dat.vid.info.FrameRate + dat.corr.tl;
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
        xx	= (dat.aud.tint - dat.corr.tl)*dat.vid.info.FrameRate;
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
            Fs	= dat.vid.info.FrameRate;
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

    % reset vocal tract data structures
    function reset_vt_data( s,e )
        dat	= guidata(hGUI);
        nf = dat.vid.info.NumberOfFrames;
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
