function set = find_token( lab, varargin )
%
%   FUNCTION:
%   searches all data structures currently loaded into workspace
%   to find dataset containing specified segmentation label 'lab'
%
%   USAGE:
%   set = find_token( 'lab', (vb) );
%
%   INPUTS:
%   lab (char):     specified segmentation label
%   vb (bool):      verbosity: 0: operate silently
%
%   EXAMPLE:
%   set = find_token( 'bacco' );
%   set = find_token( 'bacco',0 );
%

    if (nargin<2)
        vb = 1;
    else
        vb = varargin{1};
    end;

    dat     = evalin('base','whos');
    ndat	= length(dat);
    set     = {};
    if (vb), fprintf('\n'); end;
    for i = 1:ndat
        if strcmp(dat(i).class,'struct')
            nv	= dat(i).name;
            if evalin('base', ['isfield( ' nv ',''seg'' )'])
                if evalin('base', ['isfield( ' nv '.seg, ''' lab ''' )'])
                    if (vb), fprintf('   Segmental label ''%s'' found in data structure ''%s''\n', lab,nv ); end;
                    set	= union(set, nv);
                end
            end
        end
    end
    if isempty(set)
        if (vb), fprintf('   Segmental label ''%s'' not found in any data structures in workspace.\n', lab ); end;
    end
    if (vb), fprintf('\n\n'); end;

end %of main function
