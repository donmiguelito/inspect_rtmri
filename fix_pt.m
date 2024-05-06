function fix_pt( dat,gl,bnd, newpt )
%
%   FUNCTION:
%   replaces all points on specified gridline in segmented frames 
%   in specified data structure with new coordinate pair [newpt]
%
%   USAGE:
%   fix_pts( dat, gl, bnd, newpt );
%
%   INPUTS:
%   dat (string):       name of data structure to update
%   gl (int):           gridline on which to update coordinate
%   bnd (string):       'i'nner or 'o'uter: boundary to update
%   newpt [int int]:	[x y] updated coordinate pair
%
%   EXAMPLE:
%   fix_pt( 'le2',19,'o', [] );
%

    txt	= [ 'length(' char(dat) '.vt)' ];
    nf	= evalin( 'base', txt );
    for f = 1:nf
        txt	= [ 'isempty(' dat '.vt(' num2str(f) ').pts)' ];
        if ~evalin( 'base', txt )
            if strcmp( bnd, 'i' )
                txt	= [ dat '.vt(' num2str(f) ').pts(' num2str(gl) ').lf = [' num2str(newpt(1)) ' ' num2str(newpt(2)) '];' ];
            else
                txt	= [ dat '.vt(' num2str(f) ').pts(' num2str(gl) ').rt = [' num2str(newpt(1)) ' ' num2str(newpt(2)) '];' ];
            end
            evalin( 'base', txt );
        end
    end

end %of main function
