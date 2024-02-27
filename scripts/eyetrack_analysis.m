subjIDs=2:10;
% subjIDs = [2 3 7 8 9];
   
meanEye = zeros(4,6,2,length(subjIDs));

dbstop if error
for subj = 1:length(subjIDs)
    subjNum = subjIDs(subj);
    
    disp(['working on: ' num2str(subjNum)])
    
types={'cc' 'sr' 'b'};

dir = '/Volumes/phelpslab2/Lizzie/MMS_Ctrl_Eye/Data/behav/';

cd (dir)

%% Determine accurate trials
nBlocks=24;
nTrials=24;
typeCol = 6;
accCol = 18;
blockCol = 2;
rtCol = 15;

runTrials={1:4; 5:8; 9:12; 13:16; 17:20; 21:24};

fid = fopen(['STPsr_cc_eye.' num2str(subjNum) 'dat.txt']);
[importSearch] = textscan(fid, '%d %d %d %d %s %s %d %d %d %d %4.0f %4.0f %d %d %4.0f %d %d %d %d %*d %*d %f %f %f %f %f',nBlocks*nTrials,...
    'delimiter', '\t', 'Headerlines', 1);
fclose(fid);

for run = 1:length(runTrials)
    trials=find(ismember(importSearch{blockCol},runTrials{run}));
    ccTrials=find(strcmp(importSearch{typeCol}(trials),'cc'));
    srTrials=find(strcmp(importSearch{typeCol}(trials),'s-r'));
    bTrials=find(strcmp(importSearch{typeCol}(trials),'baseline'));
    
    acc.cc{run} = find(importSearch{accCol}(trials(ccTrials))==1);
    acc.sr{run} = find(importSearch{accCol}(trials(srTrials))==1);
    acc.b{run} = find(importSearch{accCol}(trials(bTrials))==1);
    
    rt.cc{run} = importSearch{rtCol}(trials(ccTrials(acc.cc{run})));
    rt.sr{run} = importSearch{rtCol}(trials(srTrials(acc.sr{run})));
    rt.b{run} = importSearch{rtCol}(trials(bTrials(acc.b{run})));
end


%% Read in eye tracking data
% terminal (tsch): edf2asc
cd ../eye/

clear vars trialInfo 
% run=1;
for run = 1:6

    clear vars On Off msgs
    
% modified from edfasc2mat
fid=fopen(['st' num2str(subjNum) '_' num2str(run) '.asc']);
dat=textscan(fid,'%s%s%s%s%s%s%s%s%s%s%s%s','delimiter',char(9)); %Read maximum number of columns of 12
fclose(fid);

for i=1:7,chdat{i}=char(dat{i});end;

%Get number of lines
nlines=numel(dat{1});

% Find all "MSG" instances
mi=1;
for j = 1:nlines
    idx = strfind(char(dat{1}(j,:)),'MSG');
    if ~isempty(idx)
        msgs(mi) = j;
        mi=mi+1;
    end
end

% Flag start of trials by type
oi=1; ci=1; si=1; bi=1; 
for m = 1:length(msgs)
    onIdx = strfind(char(dat{2}(msgs(m),:)),'stimOn');
    
    if ~isempty(onIdx)
        
        onStore(oi) = m;
        oi=oi+1;
        
        % determine trial type
        ccT=strfind(char(dat{2}(msgs(m),:)),'cc');
        if ~isempty(ccT)
            On.cc(ci) = msgs(m);
            ci=ci+1;
        else
            srT=strfind(char(dat{2}(msgs(m),:)),'s-r');
            if ~isempty(srT)
                On.sr(si) = msgs(m);
                si=si+1;
            else
                bT=strfind(char(dat{2}(msgs(m),:)),'baseline');
                if ~isempty(bT)
                    On.b(bi) = msgs(m);
                    bi=bi+1;
                end
            end
        end
    end
end

% Find trials without end
timeout = 0;

ti=1;
for t = 2:length(onStore)
    diff = onStore(t)-onStore(t-1);
    if diff == 1
        timeout(ti) = onStore(t-1);
        ti=ti+1;
    end
end

% Flag end of trials by type
ci2=1; si2=1; bi2=1; n=1;

for m = 1:length(msgs)
    
    
    if ~ismember(m, timeout)    
       
        % look for trial end
        offIdx = strfind(char(dat{2}(msgs(m),:)),'stimOff');
        
        if ~isempty(offIdx)
            off_store(n)=m;
            n=n+1;
            
            % determine trial type
            ccT2=strfind(char(dat{2}(msgs(m),:)),'cc');
            if ~isempty(ccT2)
                Off.cc(ci2) = msgs(m);
                ci2=ci2+1;
            else
                srT2=strfind(char(dat{2}(msgs(m),:)),'s-r');
                if ~isempty(srT2)
                    Off.sr(si2) = msgs(m);
                    si2=si2+1;
                else
                    bT2=strfind(char(dat{2}(msgs(m),:)),'baseline');
                    if ~isempty(bT2)
                        Off.b(bi2) = msgs(m);
                        bi2=bi2+1;
                    end
                end
            end
            
        end
    
    
    else
        
         % auto-end if no response
        %                     disp(['timeout: ' num2str(msgs(onStore(oi-2))) ': ' num2str(msgs(onStore(oi-1)))])
        if strfind(char(dat{2}(msgs(m))),'cc')
            Off.cc(ci2) = msgs(m)+1000;
            ci2=ci2+1;
        elseif strfind(char(dat{2}(msgs(m))),'s-r')
            Off.sr(si2) = msgs(m)+1000;
            si2=si2+1;
        elseif strfind(char(dat{2}(msgs(m))),'baseline')
            Off.b(bi2) = msgs(m)+1000;
            bi2=bi2+1;
        end
        
    end
        
end


% Count events within each accurate trial
for ty = 1:length(types)
    %         for tr = 1:length(On.(types{ty}))
    for ta = 1:length(acc.(types{ty}){run})
        tr = acc.(types{ty}){run}(ta);
        tOn = On.(types{ty})(tr);
        tOff = Off.(types{ty})(tr);
        nBlink = length(find(chdat{1}(tOn:tOff,1)=='S'&chdat{1}(tOn:tOff,2)=='B'&...
            chdat{1}(tOn:tOff,3)=='L'&chdat{1}(tOn:tOff,4)=='I'&...
            chdat{1}(tOn:tOff,5)=='N'&chdat{1}(tOn:tOff,6)=='K'));
        nFix = length(find(chdat{1}(tOn:tOff,1)=='S'&chdat{1}(tOn:tOff,2)=='F'&...
            chdat{1}(tOn:tOff,3)=='I'&chdat{1}(tOn:tOff,4)=='X'&...
            chdat{1}(tOn:tOff,5)==' '));
        nSacc = length(find(chdat{1}(tOn:tOff,1)=='S'&chdat{1}(tOn:tOff,2)=='S'&...
            chdat{1}(tOn:tOff,3)=='A'&chdat{1}(tOn:tOff,4)=='C'&...
            chdat{1}(tOn:tOff,5)=='C'));
        trialInfo.(types{ty}){run}(ta,1)=nFix;
        trialInfo.(types{ty}){run}(ta,2)=nSacc;
        trialInfo.(types{ty}){run}(ta,3)=nBlink;
    end
end

end

%% Look at performance over time
for meas=1:2
    for ty = 1:length(types)
        meanEye(ty,:,meas,subj)=cat(1,mean(trialInfo.(types{ty}){1}(:,meas)),mean(trialInfo.(types{ty}){2}(:,meas)),...
            mean(trialInfo.(types{ty}){3}(:,meas)),mean(trialInfo.(types{ty}){4}(:,meas)),...
            mean(trialInfo.(types{ty}){5}(:,meas)),mean(trialInfo.(types{ty}){6}(:,meas)));
        meanRT(ty,:,subj) = cat(1,mean(rt.(types{ty}){1}),mean(rt.(types{ty}){2}),...
            mean(rt.(types{ty}){3}), mean(rt.(types{ty}){4}), mean(rt.(types{ty}){5}),...
            mean(rt.(types{ty}){6}));
%         meanErr(ty,:)=cat(1,std(trialInfo.(types{ty}){1}(:,meas))/sqrt(length(trialInfo.(types{ty}){1})),...
%             std(trialInfo.(types{ty}){2}(:,meas))/sqrt(length(trialInfo.(types{ty}){2})),...
%             std(trialInfo.(types{ty}){3}(:,meas))/sqrt(length(trialInfo.(types{ty}){3})),...
%             std(trialInfo.(types{ty}){4}(:,meas))/sqrt(length(trialInfo.(types{ty}){4})),...
%             std(trialInfo.(types{ty}){5}(:,meas))/sqrt(length(trialInfo.(types{ty}){5})),...
%             std(trialInfo.(types{ty}){6}(:,meas))/sqrt(length(trialInfo.(types{ty}){6})));
    end
    
%     figure, clf
%     errorBarSeries(meanDat',meanErr')
%     legend(types)
%     if meas==1
%         title('Eye movement per run: fixations')
%     elseif meas==2
%         title('Eye movement per run: saccades')
%     end
end

%% correlate eye movement with rt

if subjNum ~= 6
allRT = cat(1,rt.cc{1:6},rt.sr{1:6},rt.b{1:6});
allEye = cat(1,trialInfo.cc{1:6},trialInfo.sr{1:6},trialInfo.b{1:6});
else
    allRT = cat(1,rt.cc{1:5},rt.sr{1:5},rt.b{1:5});
allEye = cat(1,trialInfo.cc{1:5},trialInfo.sr{1:5},trialInfo.b{1:5});
end
for meas = 1:2
    eye_rt_corr(meas,subj) = corr(allRT,allEye(:,meas));
end
for t = 1:length(types)
for run = 1:6
allRT_store.(types{t})(run,subj) = mean(rt.(types{t}){run});
end
end

allEye_store{subj} = allEye;

end

%% plot eye movements over time
addpath(genpath('/Applications/matlab7.14/toolbox/stats/stats/'))
measName={'Fixations', 'Saccades'};
learnNum=[2 3 7 8 9];
learn = find(ismember(subjIDs, learnNum));
denom=sqrt(length(learn));

% learn=1:5;
figure, clf

fixdat=cat(1,squeeze(mean(meanEye(3,4:6,1,learn)))',squeeze(mean(meanEye(1,4:6,1,learn)))',...
    squeeze(mean(meanEye(2,4:6,1,learn)))');
saccdat=cat(1,squeeze(mean(meanEye(3,4:6,2,learn)))',...
    squeeze(mean(meanEye(1,4:6,2,learn)))',squeeze(mean(meanEye(2,4:6,2,learn)))');
rtdat=cat(1,mean(allRT_store.b(4:6,learn)),mean(allRT_store.cc(4:6,learn)),...
    mean(allRT_store.sr(4:6,learn)));
for k = 1:size(rtdat,2)
    rtsplit(:,:,k)=rtdat(:,k)/1000;
    fixsplit(:,:,k)=fixdat(:,k);
    saccsplit(:,:,k)=saccdat(:,k);
end

sechalf=cat(2,fixsplit,saccsplit,rtsplit);

for meas=1:2
disp(measName{meas})
[h,p, ~, stats]=ttest(mean(squeeze(meanEye(3,4:6,meas,learn))),mean(squeeze(meanEye(1,4:6,meas,learn))));
disp(['B v CC: p = ' num2str(p)])
[h,p]=ttest(mean(squeeze(meanEye(3,4:6,meas,learn))),mean(squeeze(meanEye(2,4:6,meas,learn))));
disp(['B v SR: p = ' num2str(p)])
[h,p]=ttest(mean(squeeze(meanEye(1,4:6,meas,learn))),mean(squeeze(meanEye(2,4:6,meas,learn))));
disp(['CC v SR: p = ' num2str(p)])
end
disp('RT')
[h,p]=ttest(mean(allRT_store.b(4:6,learn)),mean(allRT_store.cc(4:6,learn)));
disp(['B v CC: p = ' num2str(p)])
[h,p]=ttest(mean(allRT_store.b(4:6,learn)),mean(allRT_store.sr(4:6,learn)));
disp(['B v SR: p = ' num2str(p)])
[h,p]=ttest(mean(allRT_store.sr(4:6,learn)),mean(allRT_store.cc(4:6,learn)));
disp(['SR v CC: p = ' num2str(p)])

cd ~/Dropbox/NYU/lab/MMS/fMRI/experiment/
sechalf_mean = mean(sechalf,3);
h=bar(sechalf_mean');
set(gca,'XTickLabel',{'No. Fixation', 'No. Saccades', 'RT'})
h_ylab = ylabel('Performance (varied)');
h_xlab = xlabel('Measure');
h_title = title('Performance: Second Half');
% set(gca,'YLim',[0 40])
subj_scatter(h,sechalf_mean,sechalf,length(learn))

figure, clf

for meas = 1:2
subplot(1,3,meas)
ccEye=squeeze(meanEye(1,:,meas,learn));
srEye=squeeze(meanEye(2,:,meas,learn));
bEye=squeeze(meanEye(3,:,meas,learn));

disp(measName{meas})
for epoch = 1:6
        disp('CC v B')
        [h,p]=ttest(ccEye(epoch,:),bEye(epoch,:));
        disp(['E' num2str(epoch) ': p = ' num2str(p)]);
        
        disp('SR v B')
        [h,p]=ttest(srEye(epoch,:),bEye(epoch,:));
        disp(['E' num2str(epoch) ': p = ' num2str(p)]);
        
        disp('SR v CC')
        [h,p]=ttest(srEye(epoch,:),ccEye(epoch,:));
        disp(['E' num2str(epoch) ': p = ' num2str(p)]);
end

for s = 1:length(learn)
    dat_sub(:,:,s) = cat(2,bEye(:,s),ccEye(:,s),srEye(:,s))';
end
dat_epoch = mean(dat_sub,3);
h=bar(dat_epoch');
subj_scatter(h,dat_epoch,dat_sub,length(learn))
legend('No Cue', 'CC', 'SR')
xlabel('Epoch')
ylabel(['No. ' measName{meas}])
title(measName{meas})
end 

subplot(1,3,3)
for s = 1:length(learn)
rt_sub(:,:,s)=cat(2,mean(allRT_store.b(:,learn(s)),2)/1000,mean(allRT_store.cc(:,learn(s)),2)/1000,...
    mean(allRT_store.sr(:,learn(s)),2)/1000)';
end
rt_mean = mean(rt_sub,3);
h=bar(rt_mean')
subj_scatter(h,rt_mean,rt_sub,length(learn))
xlabel('Epoch')
ylabel('RT (sec)')
title('RT')

disp('RT')
for epoch = 1:6
        disp('CC v B')
        [h,p]=ttest(allRT_store.cc(epoch,learn),allRT_store.b(epoch,learn));
        disp(['E' num2str(epoch) ': p = ' num2str(p)]);
        
        disp('SR v B')
        [h,p]=ttest(allRT_store.sr(epoch,learn),allRT_store.b(epoch,learn));
        disp(['E' num2str(epoch) ': p = ' num2str(p)]);
        
        disp('SR v CC')
        [h,p]=ttest(allRT_store.sr(epoch,learn),allRT_store.cc(epoch,learn));
        disp(['E' num2str(epoch) ': p = ' num2str(p)]);
end

xlabel('Epoch')
ylabel('RT')