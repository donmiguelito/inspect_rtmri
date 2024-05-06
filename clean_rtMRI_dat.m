function dat = clean_rtMRI_dat( dat )
%
%	clean_rtMRI_dat( dat );
%
%   strip large audio and video components from rtMRI 
%   data structure to reduce memory and disc space
%

    dat.aud         = [];
    dat.sg          = [];
    dat.vid.avi     = [];
    dat.vid.avi_	= [];
    dat.cvd.avi     = [];


end %of main function
