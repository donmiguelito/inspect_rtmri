function fix_pts( gl,bnd, newpt )
%
%   FUNCTION:
%   searches all data structures currently loaded into workspace
%   replaces any pair of coordinates [oldpt] on specified gridline 
%   with new coordinate pair [newpt]
%
%   USAGE:
%   fix_pts( gl, bnd, newpt );
%
%   INPUTS:
%   gl (int):           gridline on which to update coordinate
%   bnd (string):       'i'nner or 'o'uter: boundary to update
%   newpt [int int]:	[x y] updated coordinate pair
%
%   EXAMPLE:
%   fix_pts( 19, 'o', [190 115] );
%

    x	= newpt(1); 
    y	= newpt(2);
    
    dat	= evalin('base','whos');
    for i = 1:length(dat)
        if strcmp(dat(i).class,'struct')
            nv	= dat(i).name;
            txt = [ 'isfield( ' nv ', ''vt'' );' ];
            if evalin( 'base', ['isfield( ' nv ', ''vt'' );'] )
                fprintf( '    Updating point [%d %d] on gridline %d in structure <%s>\n', x,y, gl, nv  );
                txt	= [ 'fix_pt( ''' nv ''', ' num2str(gl) ', ''' bnd ''', [' num2str(x) ' ' num2str(y) '] );' ];
                evalin( 'base', txt );
            end
        end
    end
    fprintf( '    \n'  );

end %of main function
