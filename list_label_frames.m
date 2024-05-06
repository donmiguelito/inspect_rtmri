function segs = list_label_frames(dat,s)
%
%   FUNCTION:
%   list all segment labels in specified data structure, or
%   all segment labels in all data structures in current workspace,
%   along with the start and end frames of each labelled sequence.
%
%   USAGE:
%   segs = list_label_frames();
%   segs = list_label_frames(dat);
%   segs = list_label_frames(dat,1);
%
%   INPUTS:
%   dat (string/struct):    optional input argument specifing name of data
%                           structure to search for labels
%   s (bool):               optional input flag: return field list sorted
%
%

    segs = [];
    if (nargin)
        if ischar(dat)
            if evalin('base', ['isstruct(' dat ')'])
                if evalin('base', ['isfield( ' dat ',''seg'' )'])
                    segs = evalin('base', ['fieldnames(' dat '.seg)']);
                end
                for s = 1:length(segs)
                    segname	= char(segs(s));
                    fint	= evalin( 'base', [ dat '.seg.(''' segname ''').fint'] );
                    fprintf('    %8s:\t%3d\t%3d\n', segname, fint(1),fint(2) );
                end
            else
                fprintf('\n   Specified variable ''%s'' is not a structure.\n\n', dat );
            end
        else
            if isstruct(dat)
                if isfield( dat,'seg' )
                    segs = fieldnames( dat.seg);
                end
                for s = 1:length(segs)
                    segname	= char(segs(s));
                    fint	= dat.seg.(segname).fint;
                    fprintf('    %8s:\t%3d\t%3d\n', segname, fint(1),fint(2) );
                end
            else
                fprintf('\n   Input variable ''%s'' is not a structure.\n\n', inputname(1) );
            end
        end

    else
        dat     = evalin('base','whos');
        ndat	= length(dat);
        fprintf('\n    %8s\t%8s\t%3s\t%3s\n', 'Struct','Label','Start','End' );
        fprintf(  '    %8s\t%8s\t%3s\t%3s\n', '------','-----','-----','---' );
        for i = 1:ndat'
            if strcmp(dat(i).class,'struct')
                nv	= dat(i).name;
                if evalin('base', ['isfield(' nv ',''seg'')'])
                    if ~evalin('base', ['isempty(' nv '.seg'')'])
                        toks = evalin('base', ['fieldnames(' nv '.seg)']);
                        for t = 1:length(toks)
                            segname	= char(toks(t));
                            fint	= evalin( 'base', [ nv '.seg.(''' segname ''').fint'] );
                            fprintf('    %8s\t%8s:\t%3d\t%3d\n', nv, segname, fint(1),fint(2) );
                        end
                        segs = [segs; toks];
                    end
                end
            end
        end
        fprintf( '\n' );
    end
    
    if (nargin>1)
        if (s)
            segs = sort(unique(segs));
        end
    end


end %of main function
