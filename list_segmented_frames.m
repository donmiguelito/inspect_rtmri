function ff = list_segmented_frames( dat )
%
%   FUNCTION:
%   list all frames in specified data structure containing
%   tissue segmentation data
%
%   USAGE:
%   ff = list_segmented_frames( dat );
%
%   INPUTS:
%   dat (string):       name of data structure to update
%
%   EXAMPLE:
%   ff = list_segmented_frames( le2 );
%

    ff	= [];
    nf	= length( dat.vt );
    for f = 1:nf
        if ~isempty( dat.vt(f).pts)
            ff	= [ff f];
        end
    end

end %of main function
