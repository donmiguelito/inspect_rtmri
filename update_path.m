function dat_ = update_path( dat, newpath )
%
%   FUNCTION:
%   update path to source media files in analysis data structures
%
%   USAGE:
%   dat_  = update_path( dat, newpath );
%
%   ARGUMENTS:
%   dat (struct):       input  MRI analysis data structure
%   newpath (string):   new path to media files
%   dat_ (struct):      output MRI analysis data structure
%
%   EXAMPLE:
%   gem1_ = update_path( gem1, '/home/mproctor/data_rtMRI/geminates/italian/ip1/' );
%

    dat_ = dat;
    
    if ~(isstruct(dat))
        fprintf('   Input variable ''%s'' is not a structure.\n', inputname(1) );
    else
        if ~(isfield( dat,'ffn_aud' ))
            fprintf('   Can''t find path to source audio in input structure %s.\n', inputname(1) );
        else
            audpath	= dat.ffn_aud;
            ixstr	= regexp( audpath,'wav', 'start' );
            ffn_aud = [newpath audpath(ixstr:end)];
            dat_.ffn_aud = regexprep( ffn_aud,'[\/\\]',filesep );
            fprintf('\n   Original audio source:  %s\n', audpath      );
            fprintf(  '   Updated audio source:   %s\n', dat_.ffn_aud );

            vidpath	= dat.ffn_vid;
            ixstr	= regexp( vidpath,'avi', 'start' );
            ffn_vid = [newpath vidpath(ixstr:end)];
            dat_.ffn_vid = regexprep( ffn_vid,'[\/\\]',filesep );
            fprintf('\n   Original video source:  %s\n', vidpath      );
            fprintf(  '   Updated video source:   %s\n', dat_.ffn_vid );

            trkpath	= dat.ffn_trk;
            ixstr	= regexp( trkpath,'avi', 'start' );
            ffn_trk = [newpath trkpath(ixstr:end)];
            dat_.ffn_trk = regexprep( ffn_trk,'[\/\\]',filesep );
            fprintf('\n   Original tracked video: %s\n', trkpath      );
            fprintf(  '   Updated tracked video:  %s\n', dat_.ffn_trk );

        end
    end
    fprintf( '   \n\n' );

end %of main function
