function update_all_paths( newpath )
%
%   FUNCTION:
%   update path to source media files in all data structures in workspace
%
%   USAGE:
%   update_paths( newpath );
%
%   ARGUMENTS:
%   newpath (string):	new (local) path to media files
%
%   EXAMPLE:
%   update_all_paths( '/home/mproctor/data_rtMRI/geminates/italian/ip1/' );
%

    dat     = evalin('base','whos');
    ndat	= length(dat);
    set     = {};
    
    for i = 1:ndat
        if strcmp(dat(i).class,'struct')
            nv	= dat(i).name;
            if evalin('base', ['isfield( ' nv ',''ffn_aud'' )'])
                set	= union(set, nv);
            end
        end
    end
    
    if isempty(set)
        fprintf('   No updatable data structures found in workspace.\n' );
    else
        for i = 1:length(set)
            nm	= char(set(i));
            evalstr	= [ nm ' = update_path( ' nm ',''' newpath ''' );' ];
            fprintf('   %s\n',evalstr);
            evalin( 'base', evalstr );
        end;
    end

        
end %of main function
