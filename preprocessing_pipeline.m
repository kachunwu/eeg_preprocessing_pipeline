%% Define Batch parameters
data_path = '';
output_dir = '';
eeglab_dir = '';
erp_subj = [];
filetype = '.vhdr';

%% Batch Preprocessing Pipeline
for nsubj = 1:length(erp_subj)
    
    % Define Variables
    subj_path = fullfile(output_dir,sprintf('ss%04d',erp_subj(nsubj)));
    if ~exist(subj_path, 'dir')
       mkdir(subj_path)
    end
    
    raw_name = sprintf('ss%04d_raw',erp_subj(nsubj));
    raw_path = fullfile(subj_path,strcat(raw_name,'.set'));
    
%     % Load file
%     [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
%     EEG = pop_loadset(strcat(raw_name,'.set'), fullfile(data_path,sprintf('ss%04d',erp_subj(nsubj))));
%     [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname',raw_name,'gui','off'); 
%     EEG = eeg_checkset( EEG );
    
    % Load .vhdr file
    filelist = dir(sub_dir);
    nfile = length(filelist);
    for i = 3:nfile
        if contains(filelist(i).name,filetype)
            raw = filelist(i).name;
            set_name = sprintf('ss%04d_raw',subj(nsubj));
            break
        end
        continue
    end
    
    EEG = pop_loadbv(sub_dir, raw);
    
    % Assign channel location
    EEG = pop_chanedit(EEG, 'lookup',fullfile(eeglab_dir, 'plugins\\dipfit4.3\\standard_BEM\\elec\\standard_1005.elc'));
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );
    
%     % Remove EOG channel for CUHK data
%     try
%         EEG = pop_select( EEG, 'nochannel',{'VEOGU' 'VEOGL' 'HEOGL' 'HEOGR'});
%     catch
%         warning('EOG does not exist in this dataset.');
%     end    
%     [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'overwrite','on','gui','off'); 
%     EEG = eeg_checkset( EEG );
%     
    % Resampling to 250 Hz
    EEG = pop_resample( EEG, 250);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',raw_name,'savenew',raw_path,'overwrite','on','gui','off'); 
    
    % High-pass filter at 0.01 Hz
    filter_name = char(strcat(raw_name,"_fil"));
    filter_path = fullfile(subj_path,strcat(filter_name,'.set'));
    EEG = pop_eegfiltnew(EEG, 'locutoff',0.01,'plotfreqz',0);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',filter_name,'savenew',filter_path,'gui','off'); 
    EEG = eeg_checkset( EEG );
    
    % CleanLine for removing line noise
%     cleanline_name = char(strcat(EEG.setname,"_cl"));
%     cleanline_path = fullfile(subj_path,strcat(cleanline_name,'.set'));
%     EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',[1:EEG.nbchan] ,'computepower',1,'linefreqs',50,'newversion',0,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',0,'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',4,'winstep',1);
%     [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, s,'setname',cleanline_name,'savenew',cleanline_path,'gui','off'); 
%     EEG = eeg_checkset( EEG );
    
    % Remove bad channel
    badchannel_name = char(strcat(EEG.setname,"_br"));
    badchannel_path = fullfile(subj_path,strcat(badchannel_name,'.set'));
    EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion',5,'ChannelCriterion',0.8,'LineNoiseCriterion',4,'Highpass','off','BurstCriterion','off','WindowCriterion','off','BurstRejection','off','Distance','Euclidian');
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',badchannel_name,'savenew',badchannel_path,'gui','off'); 
    EEG = eeg_checkset( EEG );                   
    
    % Interpolate bad channel
    interpolation_name = char(strcat(EEG.setname,"_ip"));
    interpolation_path = fullfile(subj_path,strcat(interpolation_name,'.set'));
    EEG = pop_interp(EEG, ALLEEG(1).chanlocs, 'spherical');
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',interpolation_name,'savenew',interpolation_path,'gui','off'); 
    EEG = eeg_checkset( EEG );
    
    % Re-reference
    reref_name = char(strcat(EEG.setname,"_reref"));
    reref_path = fullfile(subj_path,strcat(reref_name,'.set'));
    EEG = pop_chanedit(EEG, 'insert',1,'changefield',{1,'labels','CPz'},'lookup',fullfile(eeglab_dir, 'plugins\\dipfit4.3\\standard_BEM\\elec\\standard_1005.elc'));
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );
    if erp_subj(nsubj) < 2000
        EEG = pop_reref( EEG, [26 36] ,'refloc',struct('labels',{'CPz'},'type',{''},'theta',{179.5329},'radius',{0.14139},'X',{-47.318},'Y',{-0.3858},'Z',{99.432},'sph_theta',{-179.5329},'sph_phi',{64.5503},'sph_radius',{110.1175},'urchan',{1},'ref',{''},'datachan',{0}));
    else
        EEG = pop_reref( EEG, [13 19] ,'refloc',struct('labels',{'CPz'},'type',{''},'theta',{179.5329},'radius',{0.14139},'X',{-47.318},'Y',{-0.3858},'Z',{99.432},'sph_theta',{-179.5329},'sph_phi',{64.5503},'sph_radius',{110.1175},'urchan',{1},'ref',{''},'datachan',{0}));
    end
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',reref_name,'savenew',reref_path,'gui','off');

%     Segmentation
%     epoch_name = char(strcat(EEG.setname,"_epoch"));
%     epoch_path = fullfile(subj_path,strcat(epoch_name,'.set'));
%     EEG = pop_epoch( EEG, {  's13'  's23'  's33'  's43'  's53'  's57'  's58'  's63'  's67'  's68'  's73'  's77'  's78'  's83'  's87'  's88'  }, [-0.2           6], 'newname', epoch_name, 'epochinfo', 'yes');
%     [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'savenew',epoch_path,'gui','off'); 
%     EEG = eeg_checkset( EEG );          
%     
%     % ICA and ICLabel
%     ica_name = char(strcat(EEG.setname,"_ica"));
%     ica_path = fullfile(subj_path,strcat(ica_name,'.set'));
%     
%     % runica (not recommended)
%     % EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on');
%     
%     % run amica
%     outdir = fullfile(subj_path,'amicaout');
%     if ~exist(outdir, 'dir')
%        mkdir(outdir)
%     end
% 
%     aminca_rank = sum(eig(cov(double(EEG.data'))) > 1E-6);
%     runamica15_pk(EEG.data, 'num_chans', EEG.nbchan, 'num_models',1, 'outdir',outdir, 'pcakeep', aminca_rank, 'do_reject', 1, 'numrej', 15, 'rejsig', 3, 'rejint', 1);
%     EEG = pop_loadmodout(EEG,outdir);
%     EEG = pop_iclabel(EEG, 'default');
%     EEG = pop_icflag(EEG, [NaN NaN;0.9 1;0.9 1;NaN NaN;NaN NaN;NaN NaN;NaN NaN]);
%     [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',ica_name,'savenew',ica_path,'gui','off');
%     EEG = eeg_checkset( EEG );
%     
    % Update GUI
    eeglab redraw;

end

%% Create Study

% eeglab
% 
% for nsubj = 1:length(erp_subj)
% 
%     % Load dataset
%     subj_path = fullfile(output_dir,sprintf('ss%04d',erp_subj(nsubj)));
%     loadName = sprintf('ss%04d_raw_fil_cl_br_ip_reref_epoch_ica_regressor.set',erp_subj(nsubj));
%     EEG = pop_loadset('filename', loadName, 'filepath', subj_path, 'loadmode', 'info');
%     
%     % Enter EEG.subjuct.
%     EEG.subject = sprintf('ss%04d',erp_subj(nsubj)); 
%  
%     % Enter EEG.group.
%     EEG.group = 1; 
%  
%     % Store the current EEG to ALLEEG.
%     [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 0);
%     
% end
% 
% eeglab redraw;
