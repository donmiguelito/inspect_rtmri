function update_data()
%
%   FUNCTION:
%   append specified data structure to all structures currently loaded into workspace
%
%

    struc1	= 'status';
    newnm1	= 'radanal';
    defval	= 4;

    dat     = evalin('base','whos');
    ndat	= length(dat);
    fprintf('\n');
    for i = 1:ndat
        if strcmp(dat(i).class,'struct')
            nv	= dat(i).name;
            if evalin('base', ['isfield(' nv ', ''' struc1 ''' )'])
                if ~evalin('base', ['isfield(' nv '.' struc1 ',''' newnm1 ''' )'])
                    fieldstr = [ struc1 '.' newnm1 ];
                    fprintf('   Adding field ''%s'' to data structure ''%s''\n', fieldstr,nv );
                    evalin('base', [ nv '.' fieldstr '=' num2str(defval) ';']);
                end
            end
        end
    end
    fprintf('\n\n');

end %of main function
