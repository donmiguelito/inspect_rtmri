function segs = list_labels(dat,s)
%
%   FUNCTION:
%   list all segment labels in specified data structure, or
%   all segment labels in all data structures in current workspace
%
%   USAGE:
%   segs = list_labels();
%   segs = list_labels(dat);
%   segs = list_labels(dat,1);
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
            else
                fprintf('   Specified variable ''%s'' is not a structure.\n', dat );
            end
        else
            if isstruct(dat)
                if isfield( dat,'seg' )
                    segs = fieldnames(dat.seg);
                end
            else
                fprintf('   Input variable ''%s'' is not a structure.\n', inputname(1) );
            end
        end
    else
        dat     = evalin('base','whos');
        ndat	= length(dat);
        for i = 1:ndat'
            if strcmp(dat(i).class,'struct')
                nv	= dat(i).name;
                if evalin('base', ['isfield(' nv ',''seg'')'])
                    if ~evalin('base', ['isempty(' nv '.seg'')'])
                        toks = evalin('base', ['fieldnames(' nv '.seg)']);
                        segs = [segs; toks];
                    end
                end
            end
        end
    end
    
    if (nargin>1)
        if (s)
            segs = sort(unique(segs));
        end
    end


end %of main function
