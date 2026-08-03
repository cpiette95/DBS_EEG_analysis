clear

addpath(genpath('/Volumes/DBS/DBS/EEG_Analysis/OUTILS_ANALYSE'))
 
label = {'Fp1','Fp2','F3','F4','Fz','T3','T4','C3','C4','Cz','P3','P4'};
 
patient_list  = {'Temoin_Laurent_Bis','Temoin_Jean','Temoin_Isabelle','Temoin_Francois', ...
                 'Temoin_6','Temoin_7','Temoin_8','Temoin_9','Temoin_10','Temoin_Sabrina'};
side_movement = {'droite','droite','droite','droite','droite','gauche','droite','droite','droite','droite'};
mvt_names     = {'index','poing'};
 
% -------- analysis parameters --------
fs     = 1024;     % sampling rate (Hz)
dec    = 1;        % 1 = no decimation, 3 = decimate 3x
width  = 3;        % wavelet cycles
gwidth = 3;        % wavelet gaussian width (in std)
freqtarget = 4:40; % desired freqs; auto-clipped to what each epoch supports
 
%With width=7, gwidth=3, the prep window (~0.70 s) only supports ~10 Hz and up, 
%and even then the edges are wavelet-contaminated within roughly gwidth·(width/(2πf)) seconds — for a 0.7 s epoch at low frequencies that's most of the window, so trust the central time bins. Lowering width (e.g. 4–5) buys you lower frequencies at the cost of frequency resolution.

% accumulators (cells indexed {patient, movement})
pow_all_mean  = cell(numel(patient_list),2);
%pow_prep_mean = cell(numel(patient_list),2);
freq_all = cell(numel(patient_list),2);  time_all = cell(numel(patient_list),2);
%freq_prep= cell(numel(patient_list),2);  time_prep= cell(numel(patient_list),2);
 
for k =1:length(patient_list)
 
    rootpath = strcat('/Volumes/DBS/DBS/EEG_Analysis/DATA_EEG/TEMOINS/',patient_list{k});
    cd(rootpath)
 
    % contralateral frontal electrode (index into your data cell, per your rule)
    if strcmp(side_movement{k},'gauche')
        t = 12;
    else
        t = 11;
    end
 
    for m = 1:2
        load(strcat('Temoin_',mvt_names{m},'_EMG_based_full_cortical_ERP_filtered.mat'))
 
        number_trials = size(full_cortical_ERP_filtered_artefact_free{t},1);
 
        acc_all = [];  acc_prep = [];   % power accumulators
 
        for r = 1:number_trials
 
            x_all  = full_cortical_ERP_filtered_artefact_free{t}(r,:);
            x_prep = full_cortical_ERP_filtered_artefact_free{t}(r,304:1024);
 
            if dec > 1
                sig_all  = decimate(x_all , dec);
                sig_prep = decimate(x_prep, dec);
                fs_use   = fs/dec;
            else
                sig_all  = x_all;   sig_prep = x_prep;   fs_use = fs;
            end
 
            sig_all  = sig_all(:).';      % 1 x nsample (channel in row)
            sig_prep = sig_prep(:).';
 
            time_vec_all  = (0:numel(sig_all )-1)/fs_use;   % time axis at the ACTUAL rate
            time_vec_prep = (0:numel(sig_prep)-1)/fs_use;
 
            % lowest frequency that physically fits in each epoch
            fmin_all  = ceil(gwidth*width/(pi*(numel(sig_all )/fs_use)));
            %fmin_prep = ceil(gwidth*width/(pi*(numel(sig_prep)/fs_use)));
            foi_all   = freqtarget(freqtarget >= max(1,fmin_all));
            %foi_prep  = freqtarget(freqtarget >= max(1,fmin_prep));
 
            [spec_all , f_all , to_all ] = ft_specest_wavelet(sig_all , time_vec_all , ...
                'width',width,'gwidth',gwidth,'freqoi',foi_all);
            %[spec_prep, f_prep, to_prep] = ft_specest_wavelet(sig_prep, time_vec_prep, ...
            %    'width',width,'gwidth',gwidth,'freqoi',foi_prep);
 
            % power = |Fourier coeff|^2 ; output is [nchan x nfreq x ntime] -> squeeze to nfreq x ntime
            pow_all  = squeeze(abs(spec_all ).^2);
            %pow_prep = squeeze(abs(spec_prep).^2);
 
            % accumulate across trials (init on first trial)
            if isempty(acc_all)
                acc_all  = pow_all;   %acc_prep = pow_prep;
            else
                acc_all  = acc_all  + pow_all;
                %acc_prep = acc_prep + pow_prep;
            end
        end
 
        % average power across trials
        pow_all_mean{k,m}  = acc_all  ./ number_trials;
        %pow_prep_mean{k,m} = acc_prep ./ number_trials;
        freq_all{k,m} = f_all;   time_all{k,m} = to_all;
        %freq_prep{k,m}= f_prep;  time_prep{k,m}= to_prep;
    end
end
 
cd('/Volumes/DBS/DBS/EEG_Analysis/GOOD_FILES/revisions')
save('Temoins_Spectrogram_Celectrode.mat','pow_all_mean','freq_all','time_all','patient_list','mvt_names')


      

%% ---------------- plotting: one figure per patient, 2 mvts x (all / prep) ----------------
for k = 1:length(patient_list)
    figure('Name',patient_list{k},'Color','w');
    for m = 1:2
        % full epoch
        subplot(2,2,(m-1)*2+1);
        power_normalized = pow_all_mean{k,m}./nanmean(pow_all_mean{k,m},2);
        imagesc(time_all{k,m}, freq_all{k,m}, power_normalized); %_all_mean{k,m});
        axis xy; colorbar;
        xlabel('Time (s)'); ylabel('Frequency (Hz)');
        title(sprintf('%s — full', mvt_names{m}));
 
        % prep window
        subplot(2,2,(m-1)*2+2);
        imagesc(time_prep{k,m}, freq_prep{k,m}, pow_prep_mean{k,m});
        axis xy; colorbar;
        xlabel('Time (s)'); ylabel('Frequency (Hz)');
        title(sprintf('%s — prep', mvt_names{m}));
    end
    sgtitle(strrep(patient_list{k},'_','\_'));
end


figure('Name','Index tapping');
for k=1:length(patient_list)
    subplot(5,2,k); hold on ;
    power_normalized = pow_all_mean{k,m}./nanmean(pow_all_mean{k,m},2);
    imagesc(time_all{k,m}, freq_all{k,m}, power_normalized); %_all_mean{k,m});
     axis xy; 
    ylim([2 40])
    plot([1 1], [2 40], 'k-', 'LineWidth',1.5); % mvt onset at 1 s
    plot([0.3 0.3], [2 40], 'r-', 'LineWidth',1.5); % preparatory period
    colorbar; 
    xlabel('Time (s)'); ylabel('Frequency (Hz)');
    title(sprintf('%s — full', patient_list{k}));
end

figure("Name",'Fist clenching');
for k=1:length(patient_list)
    subplot(5,2,k); 
    power_normalized = pow_all_mean{k,m}./nanmean(pow_all_mean{k,m},2);
    imagesc(time_all{k,m}, freq_all{k,m}, power_normalized); %_all_mean{k,m});
    axis xy; hold on; 
    ylim([2 40])
    plot([1 1], [2 40], 'k-', 'LineWidth',1.5)
    plot([0.3 0.3], [2 40], 'r-', 'LineWidth',1.5); % preparatory period
    colorbar;
    xlabel('Time (s)'); ylabel('Frequency (Hz)');
    title(sprintf('%s — full', patient_list{k}));
end