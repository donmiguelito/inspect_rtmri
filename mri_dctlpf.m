function mov_ = mri_dctlpf( mov, cof )
%
% MRI_DCT_LPF: reduce video noise using discrete cosine transformation
%
% Useage:	mov_ =  mri_dctlpf( mov, cof )
%
% Inputs:	mov:	MRI data in Matlab mov format
%           cof:	Cut-off frequency
%
% Output:	mov_:	filtered MRI data in Matlab mov format
%
% (Michael Proctor 12-04-2010)
%

    % declare consts
    wbi     = 3;       % waitbar interval

    % fetch movie dimensions
    movlen	= size(mov,2);

    % progress bar
    hWB     = waitbar(0,'Applying DCT transform to each frame ...');
    jf      = get( get(hWB,'JavaFrame'), 'FigurePanelContainer');
    jf.getComponent(0).getRootPane.getTopLevelAncestor.setAlwaysOnTop(1);

    % apply DCT transform to each frame, truncate values<threshold, reconstruct frame
    mov_ = struct('cdata',{},'colormap',{});
    for f = 1:movlen
        if ~mod(f,wbi), waitbar(f/movlen); end;	% report progress
        im	= mov(f).cdata;
        im_	= dct2(im);
        im_(abs(im_)<cof) = 0;
        mov_(f).cdata	 = idct2(im_);
        mov_(f).colormap = mov(f).colormap;
    end
    close(hWB);
    
end

