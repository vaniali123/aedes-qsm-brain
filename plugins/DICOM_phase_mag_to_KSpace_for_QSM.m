function dicom_phase_mag_swi_to_kspace(DATA,AddInfo,savefname)
%DICOM_PHASE_MAG_SWI_TO_KSPACE Re-processes input into k-space

% Changes:
%
% 04-06/2024. B0 orientation calculation. Olli Nykänen
%    - include calculation of B0 orientation vector from
%    ImagePositionPatient / ImageOrientationPatient dicom tags
%
% 08/2024. Include Pha/Mag/Re/Im separation & support for GE/Philips. Mikko Nissi
%    - check for manufacturer
%    - separation of volumus to Re/Im
%    - fix for GE +1 -1 +1 -1 +1 -1 Re/Im mystery
%
% 06/09/2024. Updated switch case. Aliisa Paulus
%    - Added 'Siemens Healthineers' to Siemens case
%    - Included an 'otherwise' case to treat other manufacturers as Siemens

% (c) Mikko Nissi 07/2013, niss@cmrr.umn.edu

% Input check
if nargin<3
    savefname='';
end

info = DATA{1}.HDR.dicominfo;
count=length(info); % number of dicom entries, at least for old fashion GE/Philips..

% Check manufacturer
switch lower(info(1).Manufacturer)
    case 'philips'
        disp('Philips data. Trying stuff.');
        % Philips:
        % Private_2005_1011 --> M / P / I / R

        % Find M / P / I / R  mag/phase/re/im info
        if(isfield(info,'Private_2005_1011'))
            PMRItmp1={};
            [PMRItmp1{1:count,1}]=deal(info.Private_2005_1011);
            PMRI_vector=cell2mat(PMRItmp1); % put slice locations into single var
            PMRI_vector=unique(PMRI_vector);
        end

        % find echo times
        if(isfield(info,'EchoTime'))
            echotmp1={};
            [echotmp1{1:count,1}]=deal(info.EchoTime);
            te_vector=cell2mat(echotmp1); % put slice locations into single var
            te_vector=unique(te_vector);
        end

        % Now we assume the following:
        %  - data is SORTED into volumes, i.e. ordered in this fashion
        %  - # of echoes volumes of 0 (Mag)
        %  - # of echoes volumes of 1 (Pha)
        %  - # of echoes volumes of 2 (Re)
        %  - # of echoes volumes of 3 (Im)
        
        if size(DATA{1}.FTDATA,4) == length(te_vector)*length(PMRI_vector)
            % assumption fits
            % let's go for Re and Im
            index_for_Re=find(PMRI_vector=='R')-1;
            index_for_Im=find(PMRI_vector=='I')-1;
            
            Re=double(DATA{1}.FTDATA(:,:,:,index_for_Re*length(te_vector)+[1:length(te_vector)]));  % i.e. pick the third set of volumes
            Im=double(DATA{1}.FTDATA(:,:,:,index_for_Im*length(te_vector)+[1:length(te_vector)]));  % i.e. pick the fourth set of volumes

          
        else
            % Something does not add up. complain and stop
            disp('Could not figure out which volume is Re/Im. Giving up. Please figure this out and fix the code for Philips!');
            return;
        end

        cplx=Re+1i*Im;

        mag=abs(cplx);
        phase=angle(cplx);



    case 'ge medical systems'
        disp('GE data. Trying stuff.');
        % IF GE DATA, look into:
        %  info.Private_0043_102f:  0 = data1, 1=data2, 2=data3,
        %  3=data4, where data1-4 are likely: mag, phase, real, imag. Or something else.
        % Data should be split to these four categories (or maybe fewer) in
        % any case.

        % Find 0-1-2-3 mag/phase/re/im info
        if (isfield(info,'Private_0043_102f'))
            PMRItmp1={};
            [PMRItmp1{1:count,1}]=deal(info.Private_0043_102f);
            PMRI_vector=cell2mat(PMRItmp1); % put slice locations into single var
            PMRI_vector=unique(PMRI_vector);
        end

        % find echo times
        if(isfield(info,'EchoTime'))
            echotmp1={};
            [echotmp1{1:count,1}]=deal(info.EchoTime);
            te_vector=cell2mat(echotmp1); % put slice locations into single var
            te_vector=unique(te_vector);
        end

        % Now we assume the following:
        %  - data is SORTED into volumes, i.e. ordered in this fashion
        %  - # of echoes volumes of 0 (Mag)
        %  - # of echoes volumes of 1 (Pha)
        %  - # of echoes volumes of 2 (Re)
        %  - # of echoes volumes of 3 (Im)
        
        if size(DATA{1}.FTDATA,4) == length(te_vector)*length(PMRI_vector)
            % assumption fits
            % let's go for Re and Im
            index_for_Re=find(PMRI_vector==2)-1;
            index_for_Im=find(PMRI_vector==3)-1;
            
            Re=double(DATA{1}.FTDATA(:,:,:,index_for_Re*length(te_vector)+[1:length(te_vector)]));  % i.e. pick the third set of volumes
            Im=double(DATA{1}.FTDATA(:,:,:,index_for_Im*length(te_vector)+[1:length(te_vector)]));  % i.e. pick the fourth set of volumes

            % Then multiply the Re & Im with -1 1 -1 1 -1 1 -1 1 (or 1 -1 1 -1 1 -1 ...) vector to correct them
            aa=ones(1,size(DATA{1}.FTDATA,3));
            aa=reshape(aa,[numel(aa)/2 2]);
            aa(:,2)=aa(:,2)*-1;
            aa=aa';
            aa=aa(:);
            aa=permute(aa,[3 2 1]); % I'm sure there must be a simpler way of creating this vector..
            aa=repmat(aa,[1 1 1 length(te_vector)]);

            % multiply with the -1 1 -1 1 ... vector:
            Re=Re.*aa;
            Im=Im.*aa;
          
        else
            % Something does not add up. complain and stop
            disp('Could not figure out which volume is Re/Im. Giving up. Please figure this out and fix the code for GE!');
            return;
        end

        cplx=Re+1i*Im;

        mag=abs(cplx);
        phase=angle(cplx);





    case {'siemens', 'siemens healthineers'}
        % all good as is?
        disp('Siemens data. Trying stuff.');

        % Assumptions
        % - 3-D input data, first volume is Mag, second is
        %   Phase. Or multiple volumes if multiecho sequence was used, first half
        %   corresponding to Mags, latter half to Phases.
        % - Data coming from d_loader
        % - Mag scale is random @ 14bits
        % - Phase scale is -4096 ... +4095  == -pi .. +pi

        % Separate data and rescale
        mag=double(DATA{1}.FTDATA(:,:,:,1:end/2));
        mag=mag/max(mag(:));

        phase=double(DATA{1}.FTDATA(:,:,:,end/2+1:end));

        phase= (-phase/4096)*pi;
        % NOTE -- MINUS ADDED TO SWAP THE SIGN OF THE PHASE !!!!!
        %         THIS SEEMS MORE CONSISTENT WITH VARIAN, BUT MAY NOT BE
        %         CORRECT!!!!
        %   See: Haacke et al, JOURNAL OF MAGNETIC RESONANCE IMAGING 32:561?576 (2010)

        % Generate complex image
        cplx=mag.*exp(1i*phase);


    otherwise
        % Handle any other manufacturer
        disp(['Other manufacturer: ' info(1).Manufacturer '. Treating as Siemens data.']);

        % Separate data and rescale
        mag=double(DATA{1}.FTDATA(:,:,:,1:end/2));
        mag=mag/max(mag(:));

        phase=double(DATA{1}.FTDATA(:,:,:,end/2+1:end));

        phase= (-phase/4096)*pi;

        % Generate complex image
        cplx=mag.*exp(1i*phase);
end



% Now that mag / phase / cplx images should be set, work through the rest:


% Generate kspace(s), FFT'ing volume-by-volume
kspace=zeros(size(cplx));
for ii=1:size(cplx,4)
    kspace(:,:,:,ii)=ifftshift(ifftn(ifftshift(cplx(:,:,:,ii))));
end

% Re-generate complex image(s) from k-space, FFT'ing volume-by-volume
for ii=1:size(cplx,4)
    cplx(:,:,:,ii)=fftshift(fftn(fftshift(kspace(:,:,:,ii))));
end

% Stuff these to Aedes-format DATA structure
DATA{1}.KSPACE=kspace;
DATA{1}.FTDATA=abs(cplx);

% remove extra dicomheaders
if (isfield(DATA{1}.HDR,'dicominfo'))
    DATA{1}.HDR.dicominfo=DATA{1}.HDR.dicominfo(1,:);
end


if isfield(info,'PerFrameFunctionalGroupsSequence')
    % enhanced DICOM

    info2 = info(1).PerFrameFunctionalGroupsSequence; % pixel measures
    % these info should be the ssame in each dicom file
    voxel_size(1,1) = single(info2.Item_1.PixelMeasuresSequence.Item_1.PixelSpacing(1));
    voxel_size(2,1) = single(info2.Item_1.PixelMeasuresSequence.Item_1.PixelSpacing(2));
    voxel_size(3,1) = single(info2.Item_1.PixelMeasuresSequence.Item_1.SliceThickness);

    % matrix size needs to correspond the size of the data!
    matrix_size(1) = size(mag,1);
    matrix_size(2) = size(mag,2);
    matrix_size(3) = size(mag,3);
    nechoes = size(mag,4);

    % echo times
    for ii=1:nechoes
        TEs(ii) = info(ii).PerFrameFunctionalGroupsSequence.Item_1.MREchoSequence.Item_1.EffectiveEchoTime;
    end

    % B0 direction
    FF = fieldnames(info2);
    for ii=1:matrix_size(3)
        apu = getfield(info2,FF{ii});
        locs(:,ii) = apu.PlanePositionSequence.Item_1.ImagePositionPatient;
    end

    % % third is the slice direction in head scans at least, grab that here
    locs3 = locs(3,:);
    [mmax mxind] = max(locs3);
    [mmin mnind] = min(locs3);
    maxLoc = locs(:,mxind); % not sure if these are good, but lets try
    minLoc = locs(:,mnind);

    % search the direction
    Affine2D = reshape(info2.Item_1.PlaneOrientationSequence.Item_1.ImageOrientationPatient,[3 2]);
    Affine3D = [Affine2D (maxLoc-minLoc)/( (matrix_size(3)-1)*voxel_size(3))];
    B0_dir = Affine3D\[0 0 1]';

else
    % old dicom
    
    % these info should be the ssame in each dicom file
    voxel_size(1,1) = single(info(1).PixelSpacing(1));
    voxel_size(2,1) = single(info(1).PixelSpacing(2));
    voxel_size(3,1) = single(info(1).SliceThickness);

    % matrix size needs to correspond the size of the data!
    matrix_size(1) = size(mag,1);
    matrix_size(2) = size(mag,2);
    matrix_size(3) = size(mag,3);

    % grab corner(?) points for each slice and echo times as well
    locs = [info.ImagePositionPatient];
    TEs = [info.EchoTime];

    % % third is the slice direction in head scans at least, grab that here
    locs3 = locs(3,:);
    [mmax mxind] = max(locs3);
    [mmin mnind] = min(locs3);
    maxLoc = locs(:,mxind); % not sure if these are good, but lets try
    minLoc = locs(:,mnind);

    % search the direction
    Affine2D = reshape(info(1).ImageOrientationPatient,[3 2]);
    Affine3D = [Affine2D (maxLoc-minLoc)/( (matrix_size(3)-1)*voxel_size(3))];
    B0_dir = Affine3D\[0 0 1]';

   
end



% round the ridiculously small values from B0_dir
B0_dir(abs(B0_dir)<1e-3) = 0;

% save B0 direction into the data to copy it into the qsm computation
% plugin
DATA{1}.B0_dir = B0_dir;

% also save TEs for fallback
TEs = unique(TEs);
DATA{1}.TE = TEs;
DATA{1}.voxel_size = voxel_size;

% Open in Aedes
aedes(DATA);

end