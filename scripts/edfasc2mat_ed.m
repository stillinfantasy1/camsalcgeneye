%----------------------------------------------------------------------------
%EDFASC2MAT: Converts eyelink ASC files to matlab .mat files
%
%Reads ASC files converted from EDF files using EDF2ASC for the eyelink system
%And create matlab MAT files that contain the following:
%
%header: information copied from the ascii header
%rec{1..n}: indexed structure for each recording in the ASC file with:
%   rec{1..n}.eye: LEFT, RIGHT, or BOTH. (if BOTH, Xpos2, Ypos2, and PS2 fields appear)
%   rec{1..n}.Xpos1: X position for eye 1 (indicated by eye field)
%   rec{1..n}.Ypos1: Y position for eye 1
%   rec{1..n}.PS1: Pupil Size for eye 1
%   rec{1..n}.Xpos2: X position for eye 2 (if binocular recording)
%   rec{1..n}.Ypos2: Y position for eye 2
%   rec{1..n}.PS2: Pupil Size for eye 2
%   rec{1..n}.valid: boolean indicating if any data point is missing across channels
%   rec{1..n}.samplerate: sample rate of recording in Hz
%   rec{1..n}.TTLtimestamps: Time stamps of (scanner) TTL pulses
%
%EJH 2011
%----------------------------------------------------------------------------

function edfasc2mat_ed(arg_ascfile)

% Before the loop over nlines, initialize On and Off

On.cc = [];
On.sr = [];
On.b = [];
Off.cc = [];
Off.sr = [];
Off.b = [];

%Get file name
if nargin<1
    clear
    %Get the ASC file from the GUI
    [ascfile,ascpath]=uigetfile('*.asc;*.ASC', ...
    'Please specify your ASCII-file (converted from EDF)');
else
    if exist(arg_ascfile)
        %Decompose the input
        [ascpath,ascfile,ascext] = fileparts(arg_ascfile);
        ascfile=[ascfile,ascext];
        %If no path was specified, check which file is requested
        if numel(ascpath) == 0
            [ascpath,ascfile,ascext] = fileparts(which(arg_ascfile));
            ascfile=[ascfile,ascext];
        end
    else
       error(['Cannot find a file named ', arg_ascfile]);
    end
    
end



%Set data channel names
dcns={'Xpos1','Ypos1','PS1','Xpos2','Ypos2','PS2'};

%Open and read in the ASCII file
disp('Reading the ASCII file...')
fid=fopen(fullfile(ascpath,ascfile));
dat=textscan(fid,'%s%s%s%s%s%s%s%s%s%s%s%s','delimiter',char(9)); %Read maximum number of columns of 12
fclose(fid);


disp('Searching events and recordings...')

%Convert to character arrays
for i=1:7,chdat{i}=char(dat{i});,end;

%Get number of lines
nlines=numel(dat{1});

%Find START, END lines indicating start of recording
startrec = find(chdat{1}(:,1)=='S'&chdat{1}(:,2)=='T'&chdat{1}(:,3)=='A'&chdat{1}(:,4)=='R'&chdat{1}(:,5)=='T'&chdat{1}(:,6)==' ');
endrec = find(chdat{1}(:,1)=='E'&chdat{1}(:,2)=='N'&chdat{1}(:,3)=='D'&chdat{1}(:,4)==' ');

%Find all scanner TTLs
scannerTTLs = find(chdat{1}(:,1)=='M'&chdat{1}(:,2)=='S'&chdat{1}(:,3)=='G'&chdat{1}(:,4)==' ');

%EVG MOD: Flag start and end of trials by type
ci=1;si=1;bi=1; ci2=1;si2=1;bi2=1;
for j = 1:nlines
    % look for trial start
    idx = strfind(char(dat{2}(j,:)),'stimOn');
    if ~isempty(idx)
        % determine trial type
        ccT=strfind(char(dat{2}(j,:)),'cc');
        if ~isempty(ccT)
            On.cc(ci) = j;
            ci=ci+1;
        else
            srT=strfind(char(dat{2}(j,:)),'s-r');
            if ~isempty(srT)
                On.sr(si) = j;
                si=si+1;
            else
                bT=strfind(char(dat{2}(j,:)),'baseline');
                if ~isempty(bT)
                    On.b(bi) = j;
                    bi=bi+1;
                end
            end
        end
    end
    
    % look for trial end
    idx = strfind(char(dat{2}(j,:)),'stimOff');
    if ~isempty(idx)
        % determine trial type
        ccT2=strfind(char(dat{2}(j,:)),'cc');
        if ~isempty(ccT2)
            Off.cc(ci2) = j;
            ci2=ci2+1;
        else
            srT2=strfind(char(dat{2}(j,:)),'s-r');
            if ~isempty(srT2)
                Off.sr(si2) = j;
                si2=si2+1;
            else
                bT2=strfind(char(dat{2}(j,:)),'baseline');
                if ~isempty(bT2)
                    Off.b(bi2) = j;
                    bi2=bi2+1;
                end
            end
        end
    end
    
end

%Error if no start of recording can be found
if numel(startrec)==0
    error(['Could not find START label of any recording in ',ascfile]);
end

%If last recording did not end, take the end of file as end
if numel(endrec)<numel(startrec)
    endrec(numel(endrec)+1)=nlines+1;
    %Error if it still doesn't match
    if numel(endrec)<numel(startrec)
        error(['Found START label of new recording before END label of previous recording in ',ascfile]);
    end
end

%Scan the first 100 lines searching header information by finding lines starting with **
for i=1:100
    tcline=strtrim(dat{1}{i});
    if numel(tcline)>16 && strcmp(tcline(4:17),'CONVERTED FROM'),header.converted_from = strtrim(tcline(19:end));,end
    if numel(tcline)>6 && strcmp(tcline(4:7),'DATE'),header.date = strtrim(tcline(10:end));,end
    if numel(tcline)>6 && strcmp(tcline(4:7),'TYPE'),header.type = strtrim(tcline(10:end));,end
    if numel(tcline)>9 && strcmp(tcline(4:10),'VERSION'),header.version = strtrim(tcline(13:end));,end
    if numel(tcline)>8 && strcmp(tcline(4:9),'SOURCE'),header.source = strtrim(tcline(12:end));,end
    if numel(tcline)>9 && strcmp(tcline(4:10),'EYELINK'),header.eyelink = strtrim(tcline(3:end));,end
    if numel(tcline)>9 && strcmp(tcline(4:10),'CAMERA:'),header.camera = strtrim(tcline(12:end));,end
    if numel(tcline)>16 && strcmp(tcline(4:17),'SERIAL NUMBER:'),header.serial_number = strtrim(tcline(19:end));,end
    if numel(tcline)>16 && strcmp(tcline(4:17),'CAMERA_CONFIG:'),header.camera_config = strtrim(tcline(19:end));,end
end

%Find lines containing samples (starting with a number)
samplines=double(chdat{1}(:,1));
samplines=find(samplines>47 & samplines<58);


%Loop over recordings to convert the data
for crec=1:numel(startrec)

    disp(['Converting recording ',num2str(crec),'...']);
    %Determine the sample lines for the current recording
    csamplines = samplines(samplines>startrec(crec) & samplines<endrec(crec));
    
    %Check determine which eye(s) is recorded by searching string after
    %"SAMPLES' line within 20 first (header) lines of recording
    rec{crec}.eye=char(dat{3}(find(ismember(dat{1}(startrec(1):startrec(1)+20),'SAMPLES'))+startrec(crec)-1));

    %Try to read the sample rate
    try
        rec{crec}.samplerate=str2num(char(dat{5}(find(ismember(dat{1}(startrec(1):startrec(1)+20),'SAMPLES'))+startrec(crec)-1)));
    catch
        disp('Warning: Could not determine sample rate from ASC file. Please check.');
    end
    
    %If this variable says BOTH, use binocular recording. Otherwise, assume monocular
    ndcs=3;
    if strcmp(rec{crec}.eye,'BOTH'),ndcs=5;,end %Number of data channels is 5 for binocular
    
    %Copy the time variable
    disp('... channels time and valid');
    rec{crec}.time=str2num(chdat{1}(csamplines,:));
    rec{crec}.valid=ones(numel(rec{crec}.time),1);
    
    %Set missing values to zero and copy the data channels
    for cdc=1:ndcs;
        disp(['... channel ',dcns{cdc}]);
        cdcd=chdat{cdc+1}(csamplines,:);
        %Lines starting with '. ' are considered missing values and set to
        %zero in the channel. The valid vector is also set to zero.
        rec{crec}.valid(find(double(cdcd(:,1))==46 & double(cdcd(:,2))==32))=0;
        cdcd(find(double(cdcd(:,1))==46 & double(cdcd(:,2))==32),1)=char(48);
        %Copy the data
        eval(['rec{crec}.',dcns{cdc},' = ','str2num(cdcd);']);
        eval(['ccl=length(rec{crec}.',dcns{cdc},');']);
        if ccl~=length(rec{crec}.time)
            error('Check for number of samples failed.');
        end
    end
    
    %Find the scanner TTL pulses for this recording
    cscannerTTLs = scannerTTLs(scannerTTLs>startrec(crec) & scannerTTLs<endrec(crec));
    
    %Loop over TTL pulses
    for cTTL=1:length(cscannerTTLs)
        %Determine the time stamps of the TTL pulses
        ctimest=dat{2}(cscannerTTLs(cTTL));
        ctimest=textscan(char(ctimest),'%d');
        rec{1}.TTLtimestamps(cTTL)=ctimest{1};
    end
    
    %EVG MOD: Find trial on/offsets, timestamps
    types={'cc', 'sr', 'b'};
    for ty = 1:length(types)
        cOn.(types{ty}) = On.(types{ty})(On.(types{ty})>startrec(crec) & On.(types{ty})<endrec(crec));
        cOff.(types{ty}) = Off.(types{ty})(Off.(types{ty})>startrec(crec) & Off.(types{ty})<endrec(crec));
    end
    
    for ty = 1:length(types)
        
        for cTy = 1:length(cOn.(types{ty}))
            ontimest=dat{2}(cOn.(types{ty})(cTy));
            ontimest=textscan(char(ontimest),'%d');
            rec{1}.On.(types{ty})(cTy)=ontimest{1};
            
            offtimest=dat{2}(cOff.(types{ty})(cTy));
            offtimest=textscan(char(offtimest),'%d');
            rec{1}.Off.(types{ty})(cTy)=offtimest{1};
            
        end
    end
end

%Save the .mat file
[irr,matfile,irr]=fileparts(ascfile);
matfile=[matfile,'.mat'];
disp(['Saving converted data to ', fullfile(ascpath,matfile)]);
save(fullfile(ascpath,matfile),'rec','header');



