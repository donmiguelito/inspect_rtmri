function go = shift_grid( gi, off )
%
%   FUNCTION:
%   shifts all coordinates and data points in input grid 
%   structure, by specified x- and y- offsets
%
%   USAGE:
%   grid_ = shift_grid( grid, offset );
%
%   INPUTS:
%   grid (struct):      grid data structure to update
%   offset [int int]:	[x y] shift offset
%
%   EXAMPLE:
%   grid_aa2_ = shift_grid( grid_aa2, [-2 1] );
%

    go	= gi;
    go.or1	= go.or1  + off;
    go.or2	= go.or2  + off;
    go.mpal	= go.mpal + off;
    go.dent	= go.dent + off;
    go.glot	= go.glot + off;
    go.mlab	= go.mlab + off;
    
    for i = 1:length(go.ends)
        go.ends(i).x = go.ends(i).x + off(1);
        go.ends(i).y = go.ends(i).y + off(2);
    end

    for i = 1:length(go.pts)
        go.pts(i).xx = go.pts(i).xx + off(1);
        go.pts(i).yy = go.pts(i).yy + off(2);
        if ~isempty(go.palate.pts(i).rt)
            go.palate.pts(i).rt	 = go.palate.pts(i).rt + off;
        end
        if ~isempty(go.pharynx.pts(i).rt)
            go.pharynx.pts(i).rt = go.pharynx.pts(i).rt + off;
        end
    end

    
end %of main function

