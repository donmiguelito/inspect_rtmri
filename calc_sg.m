function [SG,f,t] = calc_sg( s, Fs, vb, varargin )
%
% CALC_SG: calculate and display spectrogram from signal
% Useage: [SG,f,t] = calc_sg( s, Fs, vb, (preemph),(winlen),(winstep),(Fmax),(nfft) )
%
% Input Arguments:
%   s (1xn float):      signal vector
%   Fs (int):           sampling frequency
%   vb (bool):          1: plot spectrogram; 0: don't plot sg
%
% Optional Input Arguments (defaults):
%   preemph (bool):     1: pre-emphasize signal; 0: no pre-emphasis (1)
%   winlen (float):     length (msec) of time analysis window (Fs/100 msec)
%   steplen (float):	length (msec) of inter-window increment (winlen/2)
%   Fmax (int):         maximum frequency to display (Fs/2)
%   nfft (int):         number of FFT bins to use (512)
%
% Output Arguments:
%   SG (TxF float):     spectrogram matrix
%   f (1xF float):      vector of F frequency values at which SG calculated
%   t (1xT float):      vector of T time values at which SG calculated
%
%  eg.
%	[SG,f,t] = calc_sg( s,Fs,0 );       Use default spectrogram params
%	[SG,f,t] = calc_sg( s,Fs,0,0,20 );	No preemph, 20msec windows
%

% parse inputs, or assign defaults for unspecified values
if ( nargin>=3 && nargin<=8 )
    preemph	= 1;
    winlen	= 10;
    steplen	= 5;
    Fnyq	= floor(Fs/2);
    Fmax	= Fnyq;
    nfft	= 512;
    if ( nargin >= 4 )
        preemph	= varargin{1};
        if ( nargin >= 5 )
            winlen = varargin{2};
            if ( nargin >= 6 )
                steplen = varargin{3};
                if ( nargin >= 7 )
                    Fmax = varargin{4};
                    if ( nargin == 8 )
                        nfft = varargin{5};
                    end
                end
            end
        end
    end
else
    eval('help calc_sg');
    sg=[]; f=[]; t=[];
	return;
end

    [r c] = size(s);
    if c>r, s = s'; end;    % ensure signal is a column vector
    
	% unless preemphasis disabled, take 1st difference of signal
    if (preemph)
        s = diff(s);
        %s = real( filter([1 preemph], 1, s) );
    end;

    nSamps	= length(s);                % length of signal (samples)
    sigLen	= nSamps/Fs;                % length of signal (seconds)
    winwd	= floor(winlen*Fs/1000);	% window size (samples)
    winwd	= winwd + mod(winwd,2);		% ensure window size is even
    winstep	= steplen * Fs/1000;		% overlap (fractional samples)
    nFrames = round(nSamps/winstep);
    w       = hanning(winwd);

    if (vb > 1)
        disp(' ');
        disp(['   Input signal: ' num2str(nSamps) ' samples @ '  num2str(Fs) 'Hz = ' num2str(sigLen) 'sec' ] );
        disp(['   Using ' num2str(nFrames) ' frames, length ' num2str(winlen) 'msec, step ' num2str(steplen) 'msec' ]);
        disp(['   Using ' num2str(nfft) ' frequency bins, 0 to ' num2str(Fnyq) 'Hz' ]);
        disp(' ');
    end
    
	% initialize spectrogram
    fratio	= Fmax/Fnyq;                                % find maximum FFT bin
    binmax	= floor(fratio*nfft);
    fvals	= 1:binmax;
    sg      = zeros(binmax, nFrames);                   % initialize SG for speed
    sx      = winwd/2 + 1;                              % fractional sample index
    s       = [zeros(winwd/2,1); s; zeros(winwd,1)];    % pad signal

    % create spectrogram
    for fi = 1:nFrames,
        si = round(sx);
        pf = abs(fft(w .* s(si:si+winwd-1),nfft*2));
        sg(:,fi) = pf(2:binmax+1);                      % drop DC and values>Fmax
        sx = sx + winstep;
    end;
    
	% create/condition output variables
    sg	= filter(ones(3,1)/3,1,abs(sg),[],2);           % clean holes
    sg_ = sg./max(max(sg));                             % normalize
    sg_	= ones(size(sg_))-sg_;                          % invert
    SG	= uint8(255*(sg_));
    f	= linspace(1,Fmax,binmax);
    t	= linspace(0,sigLen,nFrames);

	% display spectrogram if flagged
    if (vb)
        figure;
        imagesc( t, f, SG( fvals,:) );
        colormap(gray(256).^5);
        xlabel('msecs');
        ylabel('kHz');
        set(gca, 'ydir', 'normal', 'box','on');
        set(gca, 'YTickLabel',{'0','2','4','6','8','10','12','14','16','18','20'}, ...
        'YTick',[0 2000 4000 6000 8000 10000 12000 14000 16000 18000 20000]);
    end
    
end
