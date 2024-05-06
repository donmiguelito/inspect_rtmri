function fmts	= fetch_fmts( sig, Fs, t, vb )
%
% fmts = fetch_fmts( sig, Fs, t, vb )
%
% FETCH FMTS:	Identifies formants in input speech signal at specified time,
%               returns vector of formant frequencies.
%
% Input Arguments:
%   sig (1xn float):	speech signal
%   Fs (int):   		sampling frequency
%   t (float):   		time (secs) at which to seek formant freqs
%   vb (bool):   		verbosity:	0: work silently
%                                   1: plot formant tracks
%
% Output Arguments:
%   fmts (1xN float):	set of N formant frequencies at time t
%
% Variables which can be tweaked:
%	MaxFmtFreq:     maximum frequency of all formants
%	MaxBWFreq:      maximum bandwidth of each formant
%	SegLen:         number of samples in window around time t 
%	NumFormants:    number of formants to be chosen
%	PreEmphFact:	pre-emphasis factor
%	P:              order of LPC model
% -------------------------------------------------------------------------
% Signal segment extracted with midpoint at time t, pre-emphasized, hamming windowed.
% LPC Model applied to segments, roots calculated from coefFs, sorted by
% increasing frequency, extreme poles (f<0, f>max, BW>max) eliminated.
% If no. poles insufficient, segment remodelled using larger order LPC.
%

% Declare constants
WinRatio	= 0.0250;       % Window width (percentage of sampling frequency)
SegLen      = WinRatio*Fs;
PreEmphFact	= -0.98;
NumFormants	= 6;
MaxFmtFreq	= 5000;
MaxBWFreq	= 600;
FreqRes     = 1024;

% segment, pre-emphasize and window signal
siglen      = length(sig);
sig         = real( filter([1 PreEmphFact], 1, sig) );
seg_lft     = floor( t*Fs-SegLen/2 );
seg_rht     = seg_lft + SegLen-1;
segment     = sig( seg_lft:seg_rht );
hamwin      = hamming(SegLen);
segment     = segment.*hamwin;	% plot(segment);

% initialize search variables
FCBW	= [];
poles	= 1;

    % find formants
    segment	= [segment; zeros((FreqRes-SegLen),1)];     % adjust freq resolution
    cands	= FormantCand(segment);
    FCseg	= angle(cands)*Fs/(2*pi);
    FCBW	= abs(log(abs(cands))*Fs/pi);
    FCseg	= FCseg(find(FCBW<MaxBWFreq));
    FCseg	= FCseg(find(FCseg<MaxFmtFreq));    
    while ( NumFormants > length(FCseg) )
        cands	= FormantCand( segment, poles );
        FCseg	= angle(cands)*Fs/(2*pi);
        FCBW	= abs( log(abs(cands))*Fs/pi );
        FCseg	= FCseg(find(FCBW<MaxBWFreq));
        FCseg	= FCseg(find(FCseg<MaxFmtFreq));
        poles	= poles + 1;
    end
    fmnts	= distFormant( FCseg, NumFormants);
    ForMat	= fmnts.FM;
    prob	= 1./fmnts.dst;
    [a b]	= max(prob);
    fmts	= ForMat(:,b)';


%----- subfunction FormantCand() ------
% Return formant candidates for a signal segment
function cand = FormantCand( seg, incP )
    
    if (nargin<2)
        incP = 0;
    end
    P       = 15 + incP;	% LPC order
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
function dist = distFormant( Poles, NumFormants, LastChosenFormant )
    % Calculate distance from last formant to all candidates.
    % Return distances in row vector.

    fmnts	= (nchoosek(Poles,NumFormants))';
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
