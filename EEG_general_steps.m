%% Script pour étapes finales

addpath('/Users/charlotte.piette/Desktop/EEG_Analysis')

%% Importation of data: 
cd('/Users/charlotte.piette/Desktop/EEG_Analysis/DATA/DBS_4')
filename='/Users/charlotte.piette/Desktop/EEG_Analysis/DATA/DBS_4/9 DBS OFF ttt ON 5_8.txt';
[full_data_good,pulse_light] = importdata_EEG(filename);
save('bloc_5_8_no_filtre.mat','full_data_good','pulse_light')



% Find specific files
rootpath=strcat('/Volumes/LaCie/EEG_Analysis/DATA_EEG/TEMOINS/',Temoins{k},'/Decodage');
cd(rootpath)
theFiles = fullfile(rootpath, '/*random_sequences*.mat');
files = dir(theFiles) ; 
load(files(1).name); 

%% Rejected trials 

clear
mvt_possible = {'index','poing','marionnette'}; 
condition='Dopa_On';

seuil_rejection = 80; 

for f=1:length(mvt_possible)
    
    clearvars -except f mvt_possible condition DBS seuil_rejection

    mvt = mvt_possible{f} ;
    load(strcat(condition,'_',mvt,'_EMG_based_full_cortical_ERP_filtered.mat'))
    %load(strcat(condition,'_',mvt,'_random_sequences_full_cortical_ERP_filtered.mat'))
        
    % Delete rejected trials with artefactual points
    rejected_electro_trial=[];
    full_cortical_ERP_artefact_free=full_cortical_ERP;
    full_cortical_ERP_filtered_artefact_free=full_cortical_ERP_filtered;
    for m=[4:8,9,11,12] % most relevant electrodes
        for k=1:size(full_cortical_ERP_filtered{1},1)
            if max(abs(full_cortical_ERP_filtered{m}(k,:)))>=seuil_rejection;
                rejected_electro_trial=[rejected_electro_trial,k];
            end
        end
    end
    rejected_electro_trial=unique(rejected_electro_trial); 
    disp(rejected_electro_trial)
    
    if ~isempty(rejected_electro_trial)
        for k=[1:15] %size(full_cortical_ERP,2)
            full_cortical_ERP_artefact_free{k}(rejected_electro_trial,:)=[];
        end
        
        for k=[1:15]
            full_cortical_ERP_filtered_artefact_free{k}(rejected_electro_trial,:)=[];
        end
        
        start_mov_EEG_index_final(rejected_electro_trial) = [];
        bloc_belonging_index_final(rejected_electro_trial) = [];
        %start_mov_random_index_final(rejected_electro_trial) = [];
    end
    
    %         figure();
    %         for j=1:25
    %         subplot(5,5,j)
    %             plot(full_cortical_ERP_filtered_artefact_free{6}(j,:))
    %         end
    %         figure();
    %         for j=26:50
    %         subplot(5,5,j-25)
    %             plot(full_cortical_ERP_filtered_artefact_free{6}(j,:))
    %         end
    
    
    label={'Fp1','Fp2','F3','F4','Fz','T3','T4','C3','C4','Cz','P3','P4'};
    h=figure(); m=4;
    for j=1:12
        subplot(4,3,j); hold on; 
        plot(mean(full_cortical_ERP_filtered_artefact_free{1,m}));
        plot([1024 1024],[-5 5],'k-')
        plot([308 308],[-5 5],'k-')
        m=m+1;
        title(label{j})
    end
    sgtitle(strcat(condition,' EEG electrodes for mvt ',mvt))
    savefig(h,strcat(condition,'_',mvt,'_EEG_EMG_based_electrode_summary'),'compact')
    %savefig(h,strcat(condition,'_',mvt,'_random_sequences_electrode_summary'),'compact')
    
    save(strcat(condition,'_',mvt,'_EMG_based_full_cortical_ERP_filtered'),'full_cortical_ERP','full_cortical_ERP_filtered','full_cortical_ERP_artefact_free','full_cortical_ERP_filtered_artefact_free','video_baseline_points','video_post_onset_points','start_mov_EEG_index_final','rejected_electro_trial','seuil_rejection')
    %save(strcat(condition,'_',mvt,'_random_sequences_full_cortical_ERP_filtered'),'full_cortical_ERP','full_cortical_ERP_filtered','full_cortical_ERP_artefact_free','full_cortical_ERP_filtered_artefact_free','video_baseline_points','video_post_onset_points','start_mov_random_index_final','rejected_electro_trial','seuil_rejection')
    
end



%% Extract ERs frequency bands for decoding (based on average spectrogram from 6 trials) 

temoins_list = {'Temoin_6','Temoin_Laurent_Bis','Temoin_Jean','Temoin_Isabelle', ...
    'Temoin_Francois','Temoin_7','Temoin_8', ...
    'Temoin_9','Temoin_10','Temoin_Sabrina'};

U=1;
for Z=2:length(temoins_list)
    
    U=U+1; 
    
    addpath('/Volumes/DBS/DBS/EEG_Analysis/GOOD_FILES')
    mvt_possible = {'index','poing','marionnette'};
    condition=temoins_list{U} %'Temoin_6';
    savepath = strcat('/Volumes/DBS/DBS/EEG_Analysis/DATA_EEG/TEMOINS/',condition,'/Average_Permutations_Frequency');
    globalpath = strcat('/Volumes/DBS/DBS/EEG_Analysis/DATA_EEG/TEMOINS/',condition);
    
    cd(globalpath)
    mkdir('Average_Permutations_Frequency')
    
    Fs          = 1024;
    samples_average = 6;% trials averaged per group
    nPermut     = 25;
    nChannels   = 15;
    use_power   = true;     % true -> |Hilbert|^2 (power); false -> |Hilbert| (amplitude)
    
    % band definitions (Hz); "delta cycle" taken as 1-4 Hz
    bandNames  = {'delta','theta','alpha','beta','gamma'};
    bandRanges = [ 0.5  4;
        4  8;
        8 12;
        13 30;
        30 48];
    filtOrder  = 4;               % Butterworth order (per band)
    
    % Pre-design one zero-phase band-pass filter per band
    nyq = Fs/2;
    bpFilt = cell(numel(bandNames),1);
    for bi = 1:numel(bandNames)
        [bb,aa] = butter(filtOrder, bandRanges(bi,:)/nyq, 'bandpass');
        bpFilt{bi} = {bb, aa};
    end
    
    for f = 1:length(mvt_possible)
        
        clearvars -except mvt_possible condition savepath globalpath f ...
            Fs samples_average nPermut nChannels use_power ...
            bandNames bandRanges bpFilt temoins_list U
        
        cd(globalpath)
        mvt = mvt_possible{f};
        load(strcat('Temoin_',mvt,'_EMG_based_full_cortical_ERP_filtered.mat'))
        
        source = full_cortical_ERP_filtered_artefact_free;
        
        nTrials = size(source{1}, 1);
        number_average_trials = floor(nTrials / samples_average);
        
        for Permutations = 1:nPermut
            
            permutation_list = randperm(nTrials);
            
            % init output: env.<band>{channel} = (groups x time)
            env = struct();
            for bi = 1:numel(bandNames)
                env.(bandNames{bi}) = cell(1, nChannels);
            end
            
            for k = 1:nChannels
                for bi = 1:numel(bandNames)
                    bb = bpFilt{bi}{1};  aa = bpFilt{bi}{2};
                    for p = 1:number_average_trials
                        idx    = permutation_list(1+(p-1)*samples_average : p*samples_average);
                        trials = source{k}(idx, :);                  % 6 x time
                        env.(bandNames{bi}){k}(p,:) = band_envelope(trials, bb, aa, use_power);
                    end
                end
            end
            
            cd(savepath)
            frequency_envelopes = env;
            save(strcat(condition,'_',mvt,'_average_frequency_envelopes_', ...
                num2str(Permutations),'.mat'), 'frequency_envelopes', 'bandNames', 'bandRanges', 'Fs');
        end
        
        disp(strcat('Done mvt = ', mvt, ' | groups = ', num2str(number_average_trials), ...
            ' | trials = ', num2str(nTrials)))
    end
    
end



patient_list = {'DBS_1','DBS_2','DBS_3','DBS_4','DBS_5','DBS_6','DBS_7'};

U=0;
for Z=1:length(patient_list)
    
    U=U+1; 
    
    addpath('/Volumes/DBS/DBS/EEG_Analysis/GOOD_FILES')
    mvt_possible = {'index','poing'};
    patient = patient_list{U};
    condition='DBS_On_Dopa_On';
    savepath = strcat('/Volumes/DBS/DBS/EEG_Analysis/DATA_EEG/DBS/',patient,'/',condition,'/Average_Permutations_Frequency');
    globalpath = strcat('/Volumes/DBS/DBS/EEG_Analysis/DATA_EEG/DBS/',patient,'/',condition);
    
    cd(globalpath)
    mkdir('Average_Permutations_Frequency')
    
    Fs          = 1024;
    samples_average = 6;% trials averaged per group
    nPermut     = 25;
    nChannels   = 15;
    use_power   = true;     % true -> |Hilbert|^2 (power); false -> |Hilbert| (amplitude)
    
    % band definitions (Hz); "delta cycle" taken as 1-4 Hz
    bandNames  = {'delta','theta','alpha','beta','gamma'};
    bandRanges = [ 0.5  4;
        4  8;
        8 12;
        13 30;
        30 48];
    filtOrder  = 4;               % Butterworth order (per band)
    
    % Pre-design one zero-phase band-pass filter per band
    nyq = Fs/2;
    bpFilt = cell(numel(bandNames),1);
    for bi = 1:numel(bandNames)
        [bb,aa] = butter(filtOrder, bandRanges(bi,:)/nyq, 'bandpass');
        bpFilt{bi} = {bb, aa};
    end
    
    for f = 1:length(mvt_possible)
        
        clearvars -except mvt_possible condition savepath globalpath f ...
            Fs samples_average nPermut nChannels use_power ...
            bandNames bandRanges bpFilt patient condition patient_list U
        
        cd(globalpath)
        mvt = mvt_possible{f};
        load(strcat(condition,'_',mvt,'_EMG_based_full_cortical_ERP_filtered.mat'))
        
        source = full_cortical_ERP_filtered_artefact_free;
        
        nTrials = size(source{1}, 1);
        number_average_trials = floor(nTrials / samples_average);
        
        for Permutations = 1:nPermut
            
            permutation_list = randperm(nTrials);
            
            % init output: env.<band>{channel} = (groups x time)
            env = struct();
            for bi = 1:numel(bandNames)
                env.(bandNames{bi}) = cell(1, nChannels);
            end
            
            for k = 1:nChannels
                for bi = 1:numel(bandNames)
                    bb = bpFilt{bi}{1};  aa = bpFilt{bi}{2};
                    for p = 1:number_average_trials
                        idx    = permutation_list(1+(p-1)*samples_average : p*samples_average);
                        trials = source{k}(idx, :);                  % 6 x time
                        env.(bandNames{bi}){k}(p,:) = band_envelope(trials, bb, aa, use_power);
                    end
                end
            end
            
            cd(savepath)
            frequency_envelopes = env;
            save(strcat(condition,'_',mvt,'_average_frequency_envelopes_', ...
                num2str(Permutations),'.mat'), 'frequency_envelopes', 'bandNames', 'bandRanges', 'Fs');
        end
        
        disp(strcat('Done mvt = ', mvt, ' | groups = ', num2str(number_average_trials), ...
            ' | trials = ', num2str(nTrials)))
    end
    
end






%% Plot average cortical ERPs for each movement: 
label={'Fp1','Fp2','F3','F4','Fz','T3','T4','C3','C4','Cz','P3','P4'};

temoins_list = {'Temoin_Laurent_Bis','Temoin_Jean','Temoin_Isabelle', ...
    'Temoin_Francois','Temoin_6','Temoin_7','Temoin_8', ...
    'Temoin_9','Temoin_10','Temoin_Sabrina'};
side_movement= {'droite','droite','droite','droite','droite','gauche','droite','droite','droite','droite'}

average_patient={}
U=0;
for Z=1:length(temoins_list)
    
    U=U+1; 
    
    addpath('/Volumes/DBS/DBS/EEG_Analysis/GOOD_FILES')
    mvt_possible = {'index','poing','marionnette'};
    condition=temoins_list{U} 
    globalpath = strcat('/Volumes/DBS/DBS/EEG_Analysis/DATA_EEG/TEMOINS/',condition);
    
    cd(globalpath)
    
    Fs          = 1024;

    average_mvt = []; 
    
    for f = 1:length(mvt_possible)
                
        cd(globalpath)
        mvt = mvt_possible{f};
        load(strcat('Temoin_',mvt,'_EMG_based_full_cortical_ERP_filtered_50Hz.mat'))
        
        if strcmp(side_movement{U}, 'gauche')
            average_mvt(f,:) = median(full_cortical_ERP_filtered_artefact_free{1,12});
        else
            average_mvt(f,:) = median(full_cortical_ERP_filtered_artefact_free{1,11});
        end
        
    end
    
    average_patient{1}(U,:) = decimate(average_mvt(1,304:1023),3); 
    average_patient{2}(U,:) = decimate(average_mvt(2,304:1023),3); 
    average_patient{3}(U,:) = decimate(average_mvt(3,304:1023),3); 
    
end

time =linspace(0,703,240);
figure(); for k=1:length(temoins_list)
    subplot(5,2,k); hold on; 
plot(time,average_patient{1}(k,:),'b-')
plot(time,average_patient{2}(k,:),'g-')
plot(time,average_patient{3}(k,:),'r-')
%legend('index','poing','marionnette')
title(temoins_list{k})
xlabel('Time (ms)')
ylabel('Voltage (uV)')
end


patient_list = {'DBS_1','DBS_2','DBS_3','DBS_4','DBS_5','DBS_6','DBS_7'};
side_movement={'gauche','droite','gauche','droite','droite','droite','droite'}
average_patient={}
U=0;
for Z=1:length(patient_list)
    
    U=U+1; 
    
    addpath('/Volumes/DBS/DBS/EEG_Analysis/GOOD_FILES')
    mvt_possible = {'index','poing'};
    patient = patient_list{U};
    condition='DBS_On_Dopa_On';
    globalpath = strcat('/Volumes/DBS/DBS/EEG_Analysis/DATA_EEG/DBS/',patient,'/',condition);
       cd(globalpath)
    
    Fs          = 1024;

    average_mvt = []; 
    
    for f = 1:length(mvt_possible)
                
        cd(globalpath)
        mvt = mvt_possible{f};
        load(strcat(condition,'_',mvt,'_EMG_based_full_cortical_ERP_filtered.mat'))
        
        if strcmp(side_movement{U}, 'gauche')
            average_mvt(f,:) = median(full_cortical_ERP_filtered_artefact_free{1,7}); %F4
        else
            average_mvt(f,:) = median(full_cortical_ERP_filtered_artefact_free{1,6}); %F3
        end
        
    end
    
    average_patient{1}(U,:) = decimate(average_mvt(1,304:1023),3); 
    average_patient{2}(U,:) = decimate(average_mvt(2,304:1023),3); 
    
end

time =linspace(0,703,240);
figure(); for k=1:length(patient_list)
    subplot(4,2,k); hold on; 
plot(time,average_patient{1}(k,:),'b-')
plot(time,average_patient{2}(k,:),'g-')
%legend('index','poing','marionnette')
title(patient_list{k})
xlabel('Time (ms)')
ylabel('Voltage (uV)')
sgtitle('DBS ON')
end


%% Extract average cortical ERPs

mvt_possible = {'index','poing','marionnette'}; 
condition='Dopa_On';
mkdir('Average_Permutations')
mkdir('Decodage')
savepath = '/Volumes/LaCie/EEG_Analysis/DATA_EEG/PARK/Park_10/On/Average_Permutations';
globalpath = '/Volumes/LaCie/EEG_Analysis/DATA_EEG/PARK/Park_10/On';

for f=1:length(mvt_possible)
    
    clearvars -except mvt_possible condition EMG_based savepath globalpath f
    
    cd(globalpath) 
    
    mvt = mvt_possible{f};
    load(strcat(condition,'_',mvt,'_EMG_based_full_cortical_ERP_filtered'))

        for Permutations=1:25
            
            % Average trials together
            permutation_list = randperm(size(full_cortical_ERP_filtered_artefact_free{1},1));
            samples_average = 6;
            number_average_trials = floor(size(full_cortical_ERP_filtered_artefact_free{1},1)/samples_average);
            
            average_cortical_ERP_filtered={}; average_cortical_ERP={}; average_cortical_ERP_pre_mvt={}; average_cortical_ERP_pre_mvt={}; 
            for k=1:15%size(full_cortical_ERP_filtered_artefact_free,2)
                for p=1:number_average_trials
                    index = permutation_list(1+(p-1)*samples_average:p*samples_average);
                    average_cortical_ERP_filtered{k}(p,:) = mean(full_cortical_ERP_filtered_artefact_free{k}(permutation_list(1+(p-1)*samples_average:p*samples_average),:),1);
                    average_cortical_ERP{k}(p,:) = mean(full_cortical_ERP_artefact_free{k}(permutation_list(1+(p-1)*samples_average:p*samples_average),:),1);

                    start_pre=308 ; end_pre = 1024; 
                    average_cortical_ERP_filtered_pre_mvt{k}(p,:) = mean(full_cortical_ERP_filtered_artefact_free{k}(permutation_list(1+(p-1)*samples_average:p*samples_average),start_pre:end_pre),1);
                    average_cortical_ERP_pre_mvt{k}(p,:) = mean(full_cortical_ERP_artefact_free{k}(permutation_list(1+(p-1)*samples_average:p*samples_average),start_pre:end_pre),1);
                end
            end
            
             cd(savepath)
             save(strcat(condition,'_',mvt,'_average_cortical_EMG_based_ERP_filtered_',num2str(Permutations),'.mat'),'average_cortical_ERP','average_cortical_ERP_filtered','average_cortical_ERP_filtered_pre_mvt','average_cortical_ERP_pre_mvt')   
            
        end
       
        disp(strcat('Size for mvt =',mvt))
        disp(num2str(size(full_cortical_ERP{1,4},1)))
        disp(num2str(size(full_cortical_ERP_artefact_free{1,4},1)))
        disp(num2str(size(average_cortical_ERP_filtered{1,4},1)))

end






%% Add a jitter of on average 50 ms (from -50 to 50 ms max, uniform distribution)
mvt_possible = {'index','poing','marionnette'}; 
condition='Temoin';
savepath = '/Volumes/DBS/EEG_Analysis/DATA_EEG/TEMOINS/Temoin_6';
globalpath = '/Volumes/DBS/EEG_Analysis/DATA_EEG/TEMOINS/Temoin_6';
mkdir('Jittered_Averaged_Permutations')

sampling_rate = 1024; 

for f=1:length(mvt_possible)
    
    clearvars -except mvt_possible condition EMG_based savepath globalpath sampling_rate f
    
    cd(globalpath) 
    
    mvt = mvt_possible{f};
    load(strcat(condition,'_',mvt,'_EMG_based_full_cortical_ERP_filtered'))

    jitter_50 = randi([0 100],1,size(full_cortical_ERP_filtered_artefact_free{1,1},1)) -50; % both negatives or positives jitter
    jitter_100 = randi([0 200],1,size(full_cortical_ERP_filtered_artefact_free{1,1},1)) -100; % both negatives or positives jitter
    jitter_200 = randi([0 400],1,size(full_cortical_ERP_filtered_artefact_free{1,1},1)) -200; % both negatives or positives jitter
    total_vector_length = size(full_cortical_ERP_filtered_artefact_free{1,1},2);
    
    jittered_full_cortical_ERP_filtered_artefact_free={}; 
    for m=1:size(full_cortical_ERP_filtered_artefact_free{1,1},1)
        
        for k=1:15%size(full_cortical_ERP_filtered_artefact_free,2)
            jitter_trial = ceil(jitter_50(m)*sampling_rate/1000); 
            if jitter_trial>=0 
                jittered_full_cortical_ERP_filtered_artefact_free{1,k}(m,1:jitter_trial) = NaN*ones(jitter_trial,1);
                jittered_full_cortical_ERP_filtered_artefact_free{1,k}(m,jitter_trial+1:total_vector_length) = full_cortical_ERP_filtered_artefact_free{1,k}(m,1:total_vector_length-jitter_trial);
            elseif jitter_trial<0
                jittered_full_cortical_ERP_filtered_artefact_free{1,k}(m,1:total_vector_length-abs(jitter_trial)+1) = full_cortical_ERP_filtered_artefact_free{1,k}(m,abs(jitter_trial):total_vector_length);
                jittered_full_cortical_ERP_filtered_artefact_free{1,k}(m,total_vector_length-abs(jitter_trial)+2:total_vector_length) = NaN*ones(abs(jitter_trial)-1,1);
            end    
        end
    end
    jittered_full_cortical_ERP_filtered_artefact_free_50 = jittered_full_cortical_ERP_filtered_artefact_free; 
    
   jittered_full_cortical_ERP_filtered_artefact_free={}; 
   for m=1:size(full_cortical_ERP_filtered_artefact_free{1,1},1)        
        for k=1:15%size(full_cortical_ERP_filtered_artefact_free,2)
            jitter_trial = ceil(jitter_100(m)*sampling_rate/1000); 
            if jitter_trial>=0 
                jittered_full_cortical_ERP_filtered_artefact_free{1,k}(m,1:jitter_trial) = NaN*ones(jitter_trial,1);
                jittered_full_cortical_ERP_filtered_artefact_free{1,k}(m,jitter_trial+1:total_vector_length) = full_cortical_ERP_filtered_artefact_free{1,k}(m,1:total_vector_length-jitter_trial);
            elseif jitter_trial<0
                jittered_full_cortical_ERP_filtered_artefact_free{1,k}(m,1:total_vector_length-abs(jitter_trial)+1) = full_cortical_ERP_filtered_artefact_free{1,k}(m,abs(jitter_trial):total_vector_length);
                jittered_full_cortical_ERP_filtered_artefact_free{1,k}(m,total_vector_length-abs(jitter_trial)+2:total_vector_length) = NaN*ones(abs(jitter_trial)-1,1);
            end    
        end
    end
    jittered_full_cortical_ERP_filtered_artefact_free_100 = jittered_full_cortical_ERP_filtered_artefact_free; 

   jittered_full_cortical_ERP_filtered_artefact_free={}; 
   for m=1:size(full_cortical_ERP_filtered_artefact_free{1,1},1)        
        for k=1:15%size(full_cortical_ERP_filtered_artefact_free,2)
            jitter_trial = ceil(jitter_200(m)*sampling_rate/1000); 
            if jitter_trial>=0 
                jittered_full_cortical_ERP_filtered_artefact_free{1,k}(m,1:jitter_trial) = NaN*ones(jitter_trial,1);
                jittered_full_cortical_ERP_filtered_artefact_free{1,k}(m,jitter_trial+1:total_vector_length) = full_cortical_ERP_filtered_artefact_free{1,k}(m,1:total_vector_length-jitter_trial);
            elseif jitter_trial<0
                jittered_full_cortical_ERP_filtered_artefact_free{1,k}(m,1:total_vector_length-abs(jitter_trial)+1) = full_cortical_ERP_filtered_artefact_free{1,k}(m,abs(jitter_trial):total_vector_length);
                jittered_full_cortical_ERP_filtered_artefact_free{1,k}(m,total_vector_length-abs(jitter_trial)+2:total_vector_length) = NaN*ones(abs(jitter_trial)-1,1);
            end    
        end
    end
    jittered_full_cortical_ERP_filtered_artefact_free_200 = jittered_full_cortical_ERP_filtered_artefact_free; 
    
    
    save(strcat(condition,'_',mvt,'_jittered_full_cortical_ERP_filtered.mat'),'jittered_full_cortical_ERP_filtered_artefact_free_50','jittered_full_cortical_ERP_filtered_artefact_free_100','jittered_full_cortical_ERP_filtered_artefact_free_200')
    
    
    for Permutations=1:25
            
            % Average trials together
            permutation_list = randperm(size(full_cortical_ERP_filtered_artefact_free{1},1));
            samples_average = 6;
            number_average_trials = floor(size(full_cortical_ERP_filtered_artefact_free{1},1)/samples_average);
            
            average_cortical_ERP_filtered_jittered_50={}; average_cortical_ERP_filtered_jittered_100={}; average_cortical_ERP_filtered_jittered_200={};
            for k=1:15%size(full_cortical_ERP_filtered_artefact_free,2)
                for p=1:number_average_trials
                    index = permutation_list(1+(p-1)*samples_average:p*samples_average);
                    average_cortical_ERP_filtered_jittered_50{k}(p,:) = mean(jittered_full_cortical_ERP_filtered_artefact_free_50{k}(permutation_list(1+(p-1)*samples_average:p*samples_average),:),1);
                    average_cortical_ERP_filtered_jittered_100{k}(p,:) = mean(jittered_full_cortical_ERP_filtered_artefact_free_100{k}(permutation_list(1+(p-1)*samples_average:p*samples_average),:),1);
                    average_cortical_ERP_filtered_jittered_200{k}(p,:) = mean(jittered_full_cortical_ERP_filtered_artefact_free_200{k}(permutation_list(1+(p-1)*samples_average:p*samples_average),:),1);
                end
            end
            
             cd(savepath)
             save(strcat(condition,'_',mvt,'_average_cortical_EMG_based_ERP_filtered_jittered_',num2str(Permutations),'.mat'),'average_cortical_ERP_filtered_jittered_50','average_cortical_ERP_filtered_jittered_100','average_cortical_ERP_filtered_jittered_200')   
            
    end
end
