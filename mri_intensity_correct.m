function mov_ = mri_intensity_correct( v, SensLoc1, SensLoc2 )
%
% MRI INTENSITY CORRECT: correct for intensity bias in MRImage
%
% Useage:	mov_ = mri_intensity_correct(mov, SensLoc1, SensLoc2)
%
% Inputs:	mov:        MRI data in Matlab mov format
%           SensLoc1:   Location of 1st MRI coil (default [23 1])
%           SensLoc2:   Location of 2nd MRI coil (default [45 1])
%
% Output:	mov_:	corrected MRI data in Matlab mov format
%
% (Adam Lammert 2010; modified by Michael Proctor)
%
%	eg. v_ = mri_intensity_correct( v, [20 20], [50 20] );
%

    % declare consts
    wbi     = 7;       % waitbar interval

    % fetch video & metadata
    vid     = read(v);
    movlen	= v.NumberOfFrames;
    frameht = v.Height;
    framewd = v.Width;
    veclen	= frameht*framewd;

    % first progress bar
    hWB     = waitbar(0,'Correcting MRImage intensity: fetching image frames ...');
    jf      = get( get(hWB,'JavaFrame'), 'FigurePanelContainer');
    jf.getComponent(0).getRootPane.getTopLevelAncestor.setAlwaysOnTop(1);
    
    % reshape image matrix -> rox vector
    M = zeros(movlen,veclen);
    for f = 1:movlen
        if ~mod(f,wbi), waitbar(f/movlen); end;	% report progress
        M(f,:) = reshape(vid(:,:,1,1),1,veclen);
    end
    close(hWB);

    % Check size of images
    ImSize = sqrt(size(M,2));
    if ImSize ~= round(ImSize)
        error('The images must be square!');
    end

    % Obtain the mean image
    Mean = reshape(mean(M,1),ImSize,ImSize);

    % Get the distance maps (D1 and D2) from each coil
    D1 = zeros(ImSize,ImSize);
    D2 = zeros(ImSize,ImSize);
    for itor = 1:ImSize
        for jtor = 1:ImSize
            D1(itor,jtor) = sqrt(sum((SensLoc1-[itor jtor]).^2));
            D2(itor,jtor) = sqrt(sum((SensLoc2-[itor jtor]).^2));
        end
    end

    % Combined, quantized distance, D3
    D3 = round(D1+D2);

    % Get the average intensity at each quantized distance
    X = min(min(D3)):max(max(D3));
    I = zeros(1,length(X));
    for itor = 1:length(X)
        I(itor) = mean(Mean(D3==X(itor)));
    end

    % Assume that intensities are monotonically increasing
    I2 = zeros(1,length(I));
    I2(end) = I(end);
    for itor = length(I)-1:-1:1
        if I(itor) > I2(itor+1)
            I2(itor) = I(itor);
        else
            I2(itor) = I2(itor+1);
        end
    end

    % Build the correction map (Mask)
    Mask = zeros(ImSize,ImSize);
    for itor = 1:length(X)
        %Mask(find(D3==X(itor))) = I2(itor);
        Mask(D3==X(itor)) = I2(itor);
    end

    % second progress bar
    hWB     = waitbar(0,'Correcting individual frame intensities ...');
    jf      = get( get(hWB,'JavaFrame'), 'FigurePanelContainer');
    jf.getComponent(0).getRootPane.getTopLevelAncestor.setAlwaysOnTop(1);
    
    % Correct the sequence, generating the output
    Mask2 = reshape(Mask,1,ImSize^2);
    M2 = zeros(size(M));
    for f = 1:size(M,1)
        if ~mod(f,wbi), waitbar(f/movlen); end;	% report progress
        M2(f,:) = M(f,:)./Mask2;
    end
    close(hWB);

    M2 = M2-min(min(M2));
    M2 = M2./max(max(M2));
    M2 = M2.*255;
    
    % third progress bar
    hWB     = waitbar(0,'Reconstructing intensity-corrected AVI frames ...');
    jf      = get( get(hWB,'JavaFrame'), 'FigurePanelContainer');
    jf.getComponent(0).getRootPane.getTopLevelAncestor.setAlwaysOnTop(1);

    mov_ = zeros( frameht,framewd,1,movlen);
    for f = 1:movlen
        if ~mod(f,wbi), waitbar(f/movlen); end;	% report progress
        im = uint8(reshape(M2(f,:),frameht,framewd));
        mov_(:,:,1,f) = im;
    end
    close(hWB);
    
end

