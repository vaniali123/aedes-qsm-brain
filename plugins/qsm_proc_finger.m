function f = qsm_proc_finger(DATA, ROI, AddInfo, bgfr, inversion, thold, lambda)
% QSM Processing Function
% clear all;clc;close all;
% phase unwrapping and background field removal using MEDI toolbox
% functions and aedes

% Further modified by Mikko Nissi / 2020 for
%  - dicom header processing
%  - assumes bipolar readout gradients used: process even and odd echoes
%  separately: combine at RDF or CHI stage
%  - Using ppm complex fit
%
% modified by Olli Nykänen 2016.05.01 based on MEDI toolbox and Mikko
% Nissi's aedes plugins
% 
% Some modifications 30/05/2024
% grabbing TE for enhanced dicom and old dicom, and bruker in ms as the
% code later on sivides it by 1000

% this version incorporates multiple possib ilities, dunno how it works from
% start to end with any data, as this is done KYS data in mind
g = 2*pi*42.58;  % gyromagnetic ratio
polarity = 'unipolar'; % can be 'unipolar' or 'bipolar' or 'single' (for single echo) Depends on the acquisition, not something to tune afterwards!

% select QSM methodologyae
unwrapping = 'laplacian'; % can be 'other' or 'laplacian' 

% Default parameters
if nargin < 7 || isempty(lambda)
    lambda = 2000; % Default lambda for MEDI
end

if nargin < 6 || isempty(thold)
    thold = 0.2; % Default threshold for V-SHARP
end

if nargin < 4
    bgfr = 'pdf'; % can be 'lbv', 'vsharp', 'pdf', others if implemented background field correction
    inversion = 'sdi'; % can be 'sdi', 'medi', 'star' and 'ndi', other if implemented 
end

% changed DICOM combination script to compute B0_dir, retrieve it from the
% DATA if using Siemens' DICOM images (e.q. KYS data) 
if isfield(DATA{1},'B0_dir')
    B0_dir = DATA{1}.B0_dir;
else % field does not exist, set B0-direction manually
    B0_dir = [0 0 1];   
end

% display intermadiate results such as unwrapped phase on/off
dsp_intermediate_results = 0; %1 avaa ikkunat

if isfield(DATA{1},'PROCPAR') && ~isempty(DATA{1}.PROCPAR)
    % grab procpar
    p=DATA{1}.PROCPAR;
    B0=p.B0/1e4;  % B0 in Tesla
    if polarity == "single"
        TE=p.te;
    else
        TE = p.TE;
    end
    bruker=0;
    varian=1;
    dicom=0;
elseif isfield(DATA{1}.HDR,'FileHeader') && isfield(DATA{1}.HDR.FileHeader,'acqp') 
    % Bruker data
    bruker=1;
    dicom=0;
    varian=0;
    p=DATA{1}.HDR.FileHeader.acqp;
    TE=p.ACQ_echo_time;
%     TE=TE(AddInfo.CurrentVol);
    B0=p.BF1/42.58;
    if length(p.ACQ_echo_time)>size(DATA{1}.FTDATA,4)
        % assume this is ppm complex fit, get the echo spacing
        % echo spacing
         TE=p.ACQ_echo_time;
         esp=mean(diff(TE));
        TE=esp;
    end
    mtx=p.ACQ_size./[2 1 1];   % first dim 2x, since it is complex
    px=mean(p.ACQ_fov./mtx);    % AVG PIXEL SIZE!!! NOTE NOTE!!
    
    
elseif isfield(DATA{1}.HDR,'dicominfo')
    p=DATA{1}.HDR.dicominfo;
    if AddInfo.CurrentVol>1
        % if more than one volume, choose the correct dicominfo range
        plen=length(p)/size(DATA{1}.FTDATA,4);
        tt=[(AddInfo.CurrentVol-1)*plen+1,(AddInfo.CurrentVol)*plen];
        p=p(tt(1):tt(2));       
    end
    if size(DATA{1}.FTDATA,4)>1
        % assume this is multiecho acquisition
        % also assume dim-3 is the slice dimension
        %  nstack=size(DATA{1}.FTDATA,3);
        % alternative approach for slice dim size:
        if isfield(p,'PerFrameFunctionalGroupsSequence')
            % new dicom, not sure id this does anything as I got the TE
            % already in the combining phase
            te_vector = fieldgetter(p,1:length(p),'PerFrameFunctionalGroupsSequence.Item_1.MREchoSequence.Item_1.EffectiveEchoTime');
            TE = unique(te_vector);
        else
            % OG DICOM
            nstack=length(p)/size(DATA{1}.FTDATA,4);
            TE=vertcat(p(1:nstack:end).EchoTime);
        end
        % indices to even and odd TEs:
        ndx_o=[1:2:length(TE)-1];
        ndx_e=[2:2:length(TE)];
        
        TEo=TE(ndx_o);
        TEe=TE(ndx_e);
    else
        TE=p(1).EchoTime/1000;   % convert to s
    end
    B0=p(1).MagneticFieldStrength;
    disp(sprintf('DICOM/origin, TE=%1.5f, B0=%1.5f',TE,B0));
    bruker=0; 
    dicom=1;
    varian=0;
else
    % fallback: TE=15ms, B0=9.4T
    TE=15e-3; % assume TE if no header..
    TE=4.05e-3;  % now for bruker/multiechodata
    B0=9.4;
    bruker=0;
    dicom=0;
    varian=0;
end

    disp('Re-calculating cplx image');

if ndims(DATA{1}.KSPACE)<4
    im=fftshift(fftn(fftshift(DATA{1}.KSPACE)));
    if ~exist(DATA{1}.KSPACE)
        im = DATA{1}.FTDATA;
    end
else
    % calculate volume-by-volume
    multivol=1;
    im=zeros(size(DATA{1}.KSPACE));
    for jj=1:size(DATA{1}.KSPACE,5)
        for ii=1:size(DATA{1}.KSPACE,4)
            im(:,:,:,ii,jj)=fftshift(fftn(fftshift(DATA{1}.KSPACE(:,:,:,ii,jj))));
        end
    end
end

% fix orientation:  (this seems to depend on aedes prefs -- perhaps
if ~bruker&&(strcmpi(getpref('Aedes','ReadfidOrientImages'),'on'))&&~strcmpi(DATA{1}.DataFormat,'dcm') && getpref('Aedes','VnmrUseOldReadFcn')
    % preference is set, re-orient and NOT dicom data
    keyboard
    im=aedes_rot3d(im,2,3);
else
    % don't re-rotate data..
    % if DICOM
    % if Bruker
end

nvol=size(DATA{1}.FTDATA,4);  % number of volumes
nvol2=size(DATA{1}.FTDATA,5);  % number of secondary volumes
if nvol2>1
    disp('WARNING --- MULTIPLES OF VOLUMES --- ASSUMING AVERAGES (TE may not be correct)');
end

% grab fallback TE if such is found, something goes wrong for picking TE
% data in the above for DICOMS
if isfield(DATA{1},'TE')
    TE = DATA{1}.TE;
end

% Use ROI in the first volume as the mask
Mask=ROI(1).voxels{1}(:,:,:,1);
% keyboard
% unnecessary?
iField = im;

% if mask was not provided, generate automatically, but with bad quality 
if (~exist('Mask','var'))                     
    Mask = genMask(iField, voxel_size);
end

% compute noise level if not pre computed
if (~exist('noise_level','var'))
    noise_level = calfieldnoise(iField, Mask);
end

% Generate the Magnitude image
iMag = sqrt(sum(abs(iField).^2,4));

% scale with noise level
iField = iField/noise_level;

% Remove the last two volumes due to bad data
iField = iField(:,:,:,1:end-2);

% automatically define matrix and voxel sizes
matrix_size = size(iMag);

if isfield(DATA{1},'voxel_size') 
    % as voxel size is already needed in the B0 dir computation when combining dicoms, why not just save it 
    % the grand idea being that there would not be a need to redig this info
    % for dicoms
    voxel_size = DATA{1}.voxel_size;
else

    if varian
        dimension = 10*[DATA{1}.PROCPAR.lpe2 DATA{1}.PROCPAR.lro DATA{1}.PROCPAR.lpe2];
    elseif bruker
        dimension = 10*p.ACQ_fov;
    elseif dicom
        dimension=matrix_size.*[DATA{1}.HDR.dicominfo(1).PixelSpacing(1),DATA{1}.HDR.dicominfo(1).PixelSpacing(1),DATA{1}.HDR.dicominfo(1).SliceThickness];
    else
        % something went wrong
        keyboard
    end

    voxel_size = dimension./matrix_size;
end

% readout gradient polarity, usually all the data should be unipolar or
% single echo
switch polarity
    case 'bipolar'
        fprintf(1,'\nCalculating ppm complex fit for bipolar echoes...\n');
        [iFreq_raw,N_std] = Fit_ppm_complex_bipolar(iField);      
    
    case 'unipolar'
        fprintf(1,'\nCalculating ppm complex fit for unipolar echoes...\n');
        [iFreq_raw,N_std] = Fit_ppm_complex(iField);  

    case 'single'
        fprintf(1,'\nUsing signle echo data...\n');
        iFreq_raw = angle(iField);
        N_std = 1/(Mask+0.01); % just something could be just ones

    otherwise
        fprintf(1,'\ndont do thatn');
        keyboard
end

% unwrapping method selection
switch unwrapping
    case 'laplacian'
        fprintf(1,'\nPerforming Laplacian unwrapping...\n');
        iFreq = unwrapLaplacian(iFreq_raw, matrix_size,voxel_size); 
    case 'other'
        fprintf(1,'\nPerforming Region growing unwrapping...\n');
        iFreq = unwrapPhase(iMag, iFreq_raw, matrix_size);

    otherwise
        fprintf(1,'\ndont try that either\n');
        keyboard
end

% Mask the unwrapped map, as the part outside the Mask is not wanted in
% processing
iFreq2 = iFreq.*Mask;

%intermediate results
if dsp_intermediate_results==1
    aedes(angle(iField));
    aedes(iFreq_raw);
    aedes(iFreq2);   
end

% figure out delta-TE
if length(TE) > 2
    % seems to work most of the time
    delta_TE = (TE(3)-TE(2))*10^-3; 

elseif length(TE) == 2 % double echo data
    delta_TE = (TE(2)-TE(1)); % minnesota data in seconds
    
else % single echo!
    delta_TE = TE;
end

% background field removal method selection
switch bgfr
    case 'lbv'
        fprintf(1,'\nPerforming background field removal using LBV\n');
        tol = 0.01;
        depth = -1;
        peel = 1;
        N1 = 2000;
        N2 = 2000;
        N3 = 2000;
        RDF = LBV(iFreq,Mask,matrix_size,voxel_size,tol,depth,peel,N1,N2,N3);
        % note that LBV does not output eroded mask

    case 'pdf'
        fprintf(1,'\nPerforming background field removal using PDF\n');
        [RDF,new_mask] = PDF(iFreq, N_std, Mask,matrix_size,voxel_size, B0_dir);
        Mask = new_mask; % use eroded mask for the rest of the processing

    case 'vsharp'
        fprintf(1,'\nPerforming background field removal using V-SHARP with threshold %.2f\n', thold');
        % thold = 0.2; % if using vsharp, different values should be tested
        ksize = [2 9];

        [RDF,~,new_mask] = vsharp(iFreq,Mask,matrix_size,voxel_size,thold,ksize);
        Mask = new_mask;

    otherwise
        fprintf(1,'\nNot implemented\n');
        keyboard

end

% scale the regional difference field to ppm
RDF = RDF./(B0*g*delta_TE); % -1 for some data!!!!
if bruker
%     keyboard
    RDF = -RDF;
end

% Stuff into aedes:
DATA{1}.FTDATA = RDF;
DATA{1}.STD = N_std;
DATA{1}.magnitudi = iMag;
DATA{1}.bdir = B0_dir;
DATA{1}.msize = matrix_size;
DATA{1}.vsize = voxel_size;

if dsp_intermediate_results==1
    aedes(RDF);
end 
% return

% Dipole inversion method selection
switch inversion
    case 'medi'

        % regularization parameter for MEDI depends on the target, around
        % 2000 seems to be okay for human brains at 3T, test to find
        % optimal
        %lambda = 2000;
        %fprintf(1,'\nCalculating dipole inversion using MEDI...\n');

        fprintf(1,'\nCalculating dipole inversion using MEDI with lambda=%d...\n', lambda);
       
      
        %%%% run MEDI %%%%%
        QSM = MEDI_on(lambda, RDF, N_std, iMag, Mask, matrix_size, voxel_size, delta_TE, B0_dir);
        
    case 'sdi'
        fprintf(1,'\nCalculating dipole inversion using SDI...\n');
        % simple algorithm but streaky, streaks are difficult to suppress
        QSM = susceptibility_deconvolution_sim(RDF,Mask,B0_dir,matrix_size,voxel_size);
   
    case 'ndi'
        fprintf(1,'\nCalculating dipole inversion using NDI...\n');
        % ndi might need further tuning within the code, seems to be more
        % streaky than MEDI
        QSM = ndi(RDF,iMag,Mask,g,B0,delta_TE,B0_dir,matrix_size,voxel_size);
            
    case 'star'

        % streaking artefact reduction QSM but using MEDI instead of iLSQR,
        % the lambdas should be set for both runs off the algorithm and
        % carefully!
        fprintf(1,'\nCalculating dipole inversion using STAR...\n');

        %%%%% run 1 use low lambda to get only the high suceptibility
        %%%%% sources
        lambda = 100;
        QSM = MEDI_on(lambda, RDF, N_std, iMag, Mask, matrix_size, voxel_size, delta_TE, B0_dir);
        QSM = QSM.*Mask;

        %%%% run 2
        % second run with less regularization to get the details in STAR
        
        % STAR step - remove the effect of high susceptinbility sources
        % prior to the inversion with less regularization
        D = dipole_kernel(matrix_size, voxel_size, B0_dir);
        RDF = RDF-real(ifftn(D .* fftn(QSM)));
        
        % lambda should be around the same as the optmal MEDI
        lambda = 3000;
        QSM2 = MEDI_on(lambda, RDF, N_std, iMag, Mask, matrix_size, voxel_size, delta_TE, B0_dir);
        QSM2 = QSM2.*Mask;

        % add the results of both runs together to get a results with small
        % streakst and good overall quality, supposedly at least
        QSM = QSM+QSM2;

    otherwise
        fprintf(1,'\nnot implemented\n');
        keyboard
end

% Stuff things into aedes, the idea being that the QSM comutation could be
% repeated using saved variables starting from the RDF
DATA{1}.RDF=RDF;
DATA{1}.FTDATA=QSM;
DATA{1}.STD = N_std;
DATA{1}.magnitudi = iMag;
DATA{1}.bdir = B0_dir;
DATA{1}.msize = matrix_size;
DATA{1}.vsize = voxel_size;
DATA{1}.dTE = delta_TE;

%% Return or open data
if nargout==1
    % return if requested
    f=DATA;
else
    % or open
    aedes(DATA);
end