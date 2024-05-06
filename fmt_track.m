function [fmts, times]	= fmt_track( sig,Fs, beg,fin, Fmax,Fres,nf, win,step, vb )
%
% [fmts, times] = fmt_track( sig,Fs, beg,fin, Fmax,Fres, win, vb )
%
% FMT TRACK:	Identifies formant tracks in input speech signal,
%               returns matrix of formant frequencies at each frame.
%
% Input Arguments:
%   sig (1xn float):	speech signal
%   Fs (int):   		sampling frequency
%   beg (float):   		start time (secs) of segment
%   fin (float):   		end time (secs) of segment
%	Fmax (int):         maximum frequency of all formants	(try 5000)
%	Fres (int):         frequency resolution (nfft)         (try 512)
%	nf (int):           number of formants to track         (try 4)
%   win (float):   		width of analysis window (msec)     (try 10)
%   step (float):       length of analysis step (msec)      (try 2)
%   vb (bool):   		verbosity:	0: work silently
%                                   1: plot formant tracks
%
% Output Arguments:
%   fmts (FxT float):	row F: list of frequencies for formant N
%                       col M: set of N formant frequencies at time T
%   times (1xT float):	vector: times (secs) at which formant vals calculated
%
% Other parameters which can be tweaked:                (default)
%	MaxBWFreq:      maximum bandwidth of each formant	(600)
%	PreEmphFact:	pre-emphasis factor                 (-0.98)
%	Pinit:          initial order chosen for LPC model	(15)
% -------------------------------------------------------------------------
% Speech signal first segmented, pre-emphasized, and windowed.
% LPC Model applied to segments, roots calculated from coefFs, sorted by
% increasing frequency, extreme poles (f<0, f>max, BW>max) eliminated.
% If no. poles insufficient, segments remodelled using larger order LPC.
% Sorted poles grouped with regards to no. poles required.
% If more than one group, each group evaluated in terms of trajectory
% continuity.
%
% Original author Esfandiar Zavarehei	(5-Nov-2003)
% Modified by Michael Proctor           (19-Dec-2008)
%

% specify analysis parameters
PreEmphFact	= -0.98;
nf          = 4;
Pinit       = 15;
MaxBWFreq	= 500;

% segment, pre-emphasize signal; calculate times vector
siglen      = length(sig);
seglen      = floor( win*Fs/1000 );
steplen     = floor( step*Fs/1000 );
nsegs       = floor( fix(siglen - seglen + steplen) / steplen );
sig         = real( filter([1 PreEmphFact], 1, sig) );
times       = linspace( beg, fin, nsegs );

% initialize loop
BegPtr      = 1;
FCFreq      = [];
FCBW        = [];
hamwin      = hamming(seglen);
AddP        = 1;

    for n = 1:nsegs   
        segment     = sig(BegPtr:(BegPtr+seglen-1));            % segment signal    
        segment     = segment.*hamwin;                          % window signal
        segment     = [segment; zeros((Fres-seglen),1)];        % adjust freq resolution
        Candidates	= FormantCand(segment);
        FCseg       = angle(Candidates)*Fs/(2*pi);
        FCBW        = abs(log(abs(Candidates))*Fs/pi);
        FCseg       = FCseg(find(FCBW<MaxBWFreq));
        FCseg       = FCseg(find(FCseg<Fmax));    
        while ( nf > length(FCseg) )
            Candidates	= FormantCand(segment,AddP);
            FCseg       = angle(Candidates)*Fs/(2*pi);
            FCBW        = abs( log(abs(Candidates))*Fs/pi );
            FCseg       = FCseg(find(FCBW<MaxBWFreq));
            FCseg       = FCseg(find(FCseg<Fmax));
            AddP        = AddP + 1;
        end

        if (n>1)
            Formants = distFormant( FCseg, nf, lastfmt );
        else
            Formants = distFormant( FCseg, nf);
        end
        ForMat      = Formants.FM;
        chance      = 1./Formants.dst;
        [a b]       = max(chance);
        lastfmt     = ForMat(:,b);
        fmts(:,n)   = lastfmt;    

        AddP	= 1;
        BegPtr	= BegPtr + steplen;	% move to next analysis window

    end %for n

    if (vb)
        figure; hold on;
        for f = 1:nf
            plot( times, fmts(f,:), 'r.' );
        end
        xlabel(['Time (sec)']);
        xlim([beg fin]);
        ylabel(['Frequency (Hz)']);
        ylim([0 Fmax]);
    end %if vb


    %----- subfunction FormantCand() ------
    % Return formant candidates for a signal segment
    function cand = FormantCand( seg, incP )

        if (nargin<2)
            incP = 0;
        end
        P       = Pinit + incP;     % LPC order
        lpccof  = lpc( seg, P );
        rts     = roots(lpccof);

        ph      = angle(rts);
        amp     = abs(rts);
        ph_rng	= 0.000001;
        ph_idx	= find( ph>ph_rng & ph<(pi-ph_rng) );
        ph      = ph(ph_idx);
        amp     = amp(ph_idx);

        [ph ix]	= sort(ph);
        amp     = amp(ix);
        cand	= amp.*exp(j*ph);

    end % FormantCand()

    %----- subfunction distFormant() ------
    function dist = distFormant( Poles, nf, LastChosenFormant )
        % Calculate distance from last formant to all candidates.
        % Return distances in row vector.

        fmnts	= (nchoosek(Poles,nf))';
        nCands	= size(fmnts,2);
        if (nargin==3)
            dst	= sum( (fmnts-LastChosenFormant(:,ones(1,nCands))).^2 );
        elseif (nargin==2)
            dst = ones( 1, size(fmnts,2) );
        else
            error('wrong number of inputs')
        end
        dist.FM	 = fmnts;
        dist.dst = dst;

    end % distFormant()

end % fmt_track()
