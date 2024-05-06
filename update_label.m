function set = update_label( lab_old,lab_new, varargin )
%
%   FUNCTION:
%   searches all data sstructures currently loaded into workspace
%   to find dataset containing specified segmentation label 'lab_old'
%   and replaces those labels witn new string 'lab_new'
%
%   USAGE:
%   set = update_label( lab_old,lab_new, (vb) );
%
%   INPUTS:
%   lab_old (char): current segmentation label
%   lab_new (char): replacement segmentation label
%   vb (bool):      verbosity: 0: operate silently
%
%   EXAMPLE:
%   set = update_label( 't2','t' );
%   set = update_label( 'tamman2','tamman',0 );
%

    if (nargin<3)
        vb = 1;
    else
        vb = varargin{2};
    end;

    dat     = evalin('base','whos');
    ndat	= length(dat);
    set     = {};
    if (vb), fprintf('\n'); end;
    for i = 1:ndat
        if strcmp(dat(i).class,'struct')
            nv	= dat(i).name;
            if evalin('base', ['isfield( ' nv ',''seg'' )'])
                if evalin('base', ['isfield( ' nv '.seg, ''' lab_old ''' )'])
                    evalin('base', [ nv '.seg.(''' lab_new ''' ) = ' nv '.seg.(''' lab_old ''' )' ]);
                    evalin('base', [ nv '.seg = rmfield(' nv '.seg,''' lab_old ''' )' ]);
                    if (vb), fprintf('   Segmental label ''%s -> %s'' updated in data structure ''%s''\n', lab_old,lab_new,nv ); end;
                    set	= union(set, nv);
                end
            end
        end
    end
    if isempty(set)
        if (vb), fprintf('   Segmental label ''%s'' not found in any data structures in workspace.\n', lab_old ); end;
    end
    if (vb), fprintf('\n\n'); end;

end %of main function
