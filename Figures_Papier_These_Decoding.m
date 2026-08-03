%% Figure papier % Decoding
clear

condition='Temoin'
patient_list = {'Temoin_Laurent_Bis','Temoin_Jean','Temoin_Isabelle','Temoin_Francois','Temoin_6','Temoin_7','Temoin_8','Temoin_9','Temoin_10','Temoin_Sabrina'};
side_movement= {'droite','droite','droite','droite','droite','gauche','droite','droite','droite','droite'};
minimal_number_per_patient_permut=[37,28,30,23,33,30,18,32,29,33] ;

index_final = 12;

for k=1:length(patient_list)
    
    rootpath=strcat('/Volumes/DBS/EEG_Analysis/DATA_EEG/TEMOINS/',patient_list{k},'/Decodage');
    
    cd(rootpath)
    
    for z=1:index_final
        
        if z==1
            saving_name = 'decoding_2mvts_individual_trials_C_contra_decimate_window2_100trials';
        elseif z==2
            saving_name = 'decoding_2mvts_permutation_trials_C_contra_decimate_window2_20trials';
        elseif z==3
            saving_name = 'decoding_2mvts_individual_trials_Fp_F_C_P_contra_decimate_window2_100trials';
        elseif z==4
            saving_name = 'decoding_2mvts_individual_trials_Fp_F_C_P_ipsi_decimate_window2_100trials';
        elseif z==5
            saving_name = 'decoding_2mvts_permutation_trials_Fp_F_C_P_decimate_window2_20trials'; %contra
        elseif z==6
            saving_name = 'decoding_2mvts_individual_trials_all_decimate_window2_100trials';
        elseif z==7
            saving_name = 'decoding_2mvts_individual_trials_all_decimate_window2_50trials';
        elseif z==8
            saving_name = 'decoding_2mvts_permutation_trials_all_decimate_window2_20trials';
        elseif z==9
            saving_name = 'decoding_2mvts_permutation_trials_Fp_F_C_P_contra_decimate_window2_perpatient';
        elseif z==10
            saving_name = 'decoding_2mvts_permutation_trials_Fp_F_C_P_contra_decimate_window2_jittered_50ms_perpatient';
        elseif z==11
            saving_name = 'decoding_2mvts_permutation_trials_Fp_F_C_P_contra_decimate_window2_jittered_100ms_perpatient';
        elseif z==12
            saving_name = 'decoding_2mvts_permutation_trials_Fp_F_C_P_contra_decimate_window2_jittered_200ms_perpatient';
        end
        
        load(strcat(rootpath,'/',saving_name,'.mat'))
        variable = score_algo_LDA;
        score_mean(k,z) = mean(variable(:));
        score_std(k,z) = std(variable(:));
        
    end
        
end

for z=1:8
        [r,p]=signrank(score_mean(:,z)-50)
        pval(z)=r;
end


% Jittered decoding for temoins
figure(); hold on ;
%bar(mean(score_mean),'EdgeColor','k','FaceColor',[0.9 0.9 0.9])
for k=1:length(patient_list)
    plot(1:4,[score_mean(k,9) score_mean(k,10) score_mean(k,11) score_mean(k,12)],'o-','color',[0.8 0.8 0.8])
end
errorbar(1:4,mean(score_mean(:,9:12)),std(score_mean(:,9:12),0,1)/sqrt(length(patient_list)),'ks','Linestyle','-','LineWidth',2)
plot(1:4,50*ones(1,4),'k-.')
xlim([0.5 4+0.5])
set(gca,'xtick',1:4,'xticklabel',{'No jitter','50 ms','100m ms', '200 ms'})
xtickangle(30)
ylim([45 85])     

% figure(); hold on; 
% for z=9:12
%     scatter((z-8)*ones(length(patient_list),1),score_mean(:,z),'MarkerEdgeColor',[0.2 0.2 0.2],'jitter','on')
% end
% errorbar(1:4,mean(score_mean(:,9:12)),std(score_mean(:,9:12),0,1)/sqrt(length(patient_list)),'ks','Linestyle','none','LineWidth',2)
% plot(1:4,50*ones(1,4),'k-.')
% xlim([0.5 4+0.5])
% set(gca,'xtick',1:4,'xticklabel',{'No jitter','50 ms','100m ms', '200 ms'})
% xtickangle(30)
% 



figure(); hold on ;
%bar(mean(score_mean),'EdgeColor','k','FaceColor',[0.9 0.9 0.9])
for z=1:index_final
    scatter(z*ones(length(patient_list),1),score_mean(:,z),'MarkerEdgeColor',[0.2 0.2 0.2],'jitter','on')
end
errorbar(1:index_final,mean(score_mean),std(score_mean,0,1)/sqrt(length(patient_list)),'ks','Linestyle','none','LineWidth',2)
plot(1:index_final,50*ones(1,index_final),'k-.')
xlim([0.5 index_final+0.5])
set(gca,'xtick',1:index_final,'xticklabel',{'Individual C contra','Permutation C contra','Fp/F/C/P contra','Fp/F/C/P ipsi','Permutation Fp/F/C/P contra','Individual All','Individual All 50 trials','Permutation All'})
xtickangle(30)


% Supplementary: individual electrodes =all the same
clear
condition='Temoin'
patient_list = {'Temoin_Laurent_Bis','Temoin_Jean','Temoin_Isabelle','Temoin_Francois','Temoin_6','Temoin_7','Temoin_8','Temoin_9','Temoin_10','Temoin_Sabrina'};
side_movement= {'droite','droite','droite','droite','droite','gauche','droite','droite','droite','droite'};
minimal_number_per_patient_permut=[37,28,30,23,33,30,18,32,29,33] ;
index_final = 6;

for k=1:length(patient_list)
    
    rootpath=strcat('/Volumes/DBS/EEG_Analysis/DATA_EEG/TEMOINS/',patient_list{k},'/Decodage');
    
    cd(rootpath)
    
    for z=1:index_final
        
        if z==1
            saving_name = 'decoding_2mvts_permutation_trials_Fp_contra_decimate_window2_20trials';
        elseif z==2
            saving_name = 'decoding_2mvts_permutation_trials_F_contra_decimate_window2_20trials';
        elseif z==3
            saving_name = 'decoding_2mvts_permutation_trials_Fz_contra_decimate_window2_20trials';
        elseif z==4
            saving_name = 'decoding_2mvts_permutation_trials_C_contra_decimate_window2_20trials';
        elseif z==5
            saving_name = 'decoding_2mvts_permutation_trials_Cz_contra_decimate_window2_20trials';
        elseif z==6
            saving_name = 'decoding_2mvts_permutation_trials_P_contra_decimate_window2_20trials';
        end
        
        load(strcat(rootpath,'/',saving_name,'.mat'))
        variable = score_algo_LDA;
        score_mean(k,z) = mean(variable(:));
        score_std(k,z) = std(variable(:));
        
    end
end

figure(); hold on ;
%bar(mean(score_mean),'EdgeColor','k','FaceColor',[0.9 0.9 0.9])
for z=1:index_final
    scatter(z*ones(length(patient_list),1),score_mean(:,z),'MarkerEdgeColor',[0.2 0.2 0.2],'jitter','on')
end
errorbar(1:index_final,mean(score_mean),std(score_mean,0,1)/sqrt(length(patient_list)),'ks','Linestyle','none','LineWidth',2)
plot(1:index_final,50*ones(1,index_final),'k-.')
xlim([0.5 index_final+0.5])
set(gca,'xtick',1:index_final,'xticklabel',{'Fp','F','Fz','C','Cz','P'})
xtickangle(30)

for z=1:index_final
        [r,p]=signrank(score_mean(:,z)-50)
        pval(z)=r;
end




%% Pour deux mouvements: 
clear

patient_list = {'DBS_1','DBS_2','DBS_3','DBS_4','DBS_5','DBS_6','DBS_7'};
side_movement={'﻿gauche','droite','gauche','droite','droite','droite','droite'}
minimal_number_individual=[167,179,96,105,198,50,75] % per patient for 2 mvts
minimal_number_permutation=[27,29,16,17,33,8,12] % per patient for 2 mvts
age_DBS = [54,46,57,71,45,63,68];
side_DBS=[1,0,1,0,0,0,0]; 
pharmaco=[157,750,300,478,255,856,576];

patient_list_temoins = {'Temoin_Laurent_Bis','Temoin_Jean','Temoin_Isabelle','Temoin_Francois','Temoin_6','Temoin_7','Temoin_8','Temoin_9','Temoin_10','Temoin_Sabrina'};
side_movement= {'droite','droite','droite','droite','droite','gauche','droite','droite','droite','droite'};
minimal_individual_per_patient=[224,220,190,140,215,213,179,226,184,203] %2mvts
minimal_permutation_per_patient=[37,36,31,23,35,35,29,37,30,33] %2mvts

for k=1:length(patient_list)
    
    for z=1:8
    
    if z==1
        saving_name='decoding_2mvts_individual_trials_all_decimate_window2_perpatient';
    elseif z==2
        saving_name='decoding_2mvts_permutation_trials_C_contra_decimate_window2_perpatient';
    elseif z==3
        saving_name='decoding_2mvts_permutation_trials_all_decimate_window2_perpatient';

    elseif z==4
        saving_name='decoding_2mvts_individual_trials_Fp_F_C_P_contra_decimate_window2_perpatient';
    elseif z==5
        saving_name='decoding_2mvts_individual_trials_Fp_F_C_P_contra_decimate_window2_50trials';
    elseif z==6
        saving_name='decoding_2mvts_permutation_trials_Fp_F_C_P_contra_decimate_window2_perpatient';
    
    elseif z==7 
        saving_name='decoding_2mvts_individual_trials_EMG_electrodes_decimate_window_mvt_perpatient';
    elseif z==8
        saving_name='decoding_2mvts_permutation_trials_Fp_F_C_P_contra_decimate_window2_perpatient';
%saving_name='decoding_2mvts_permutation_trials_EMG_electrodes_decimate_window_mvt_perpatient';
    end
    
    %rootpath=strcat('/Volumes/Charlotte/EEG_Analysis/DATA_EEG/DBS/',patient_list{k},'/DBS_Off_Dopa_Off/Decodage');
    rootpath = strcat('/Volumes/DBS/DBS/EEG_Analysis/DATA_EEG/DBS/',patient_list{k},'/DBS_Off_Dopa_Off/Decodage');
    cd(rootpath)
    try
        load(strcat(rootpath,'/',saving_name,'.mat'))
        score_algo_LDA_off = score_algo_LDA(:);
        score_algo_LR_off = score_algo_LR(:);
        score_algo_NN_off = score_algo_NN(:);
        
        Off_Off_score_algo_LDA_mean(k,z) = mean(score_algo_LDA(:));
        Off_Off_score_algo_LDA_std(k,z) = std(score_algo_LDA(:)); %./sqrt(625);
        Off_Off_score_algo_NN_mean(k,z) = mean(score_algo_NN(:));
        Off_Off_score_algo_NN_std(k,z) = std(score_algo_NN(:)) ;%./sqrt(625);
        Off_Off_score_algo_LR_mean(k,z) = mean(score_algo_LR(:));
        Off_Off_score_algo_LR_std(k,z) = std(score_algo_LR(:)) ;%./sqrt(625);
    catch
        Off_Off_score_algo_LDA_mean(k,z)=NaN;
        Off_Off_score_algo_LDA_std(k,z)=NaN;
        Off_Off_score_algo_NN_mean(k,z) = NaN;
        Off_Off_score_algo_NN_std(k,z) = NaN;
        Off_Off_score_algo_LR_mean(k,z) = NaN;
        Off_Off_score_algo_LR_std(k,z)=NaN;
    end
       

    rootpath = strcat('/Volumes/DBS/DBS/EEG_Analysis/DATA_EEG/DBS/',patient_list{k},'/DBS_On_Dopa_On/Decodage');
    %rootpath=strcat('/Volumes/Charlotte/EEG_Analysis/DATA_EEG/DBS/',patient_list{k},'/DBS_On_Dopa_On/Decodage');
    cd(rootpath)
    load(strcat(rootpath,'/',saving_name,'.mat'))
    On_On_score_algo_LDA_mean(k,z) = mean(score_algo_LDA(:));
    On_On_score_algo_LDA_std(k,z) = std(score_algo_LDA(:)); %./sqrt(625);
    On_On_score_algo_NN_mean(k,z) = mean(score_algo_NN(:));
    On_On_score_algo_NN_std(k,z) = std(score_algo_NN(:)) ;%./sqrt(625);
    On_On_score_algo_LR_mean(k,z) = mean(score_algo_LR(:));
    On_On_score_algo_LR_std(k,z) = std(score_algo_LR(:)) ;%./sqrt(625);
            
    score_algo_LDA_on = score_algo_LDA(:);
    score_algo_LR_on = score_algo_LR(:);
    score_algo_NN_on = score_algo_NN(:);

    if z==6
        disp(k)
        disp(mean(score_algo_LDA_on)-mean(score_algo_LDA_off))
        [~,p]=ttest(score_algo_LDA_on,score_algo_LDA_off)
        disp(mean(score_algo_LR_on)-mean(score_algo_LR_off))
        [~,p]=ttest(score_algo_LR_on,score_algo_LR_off)
        disp(mean(score_algo_NN_on)-mean(score_algo_NN_off))
        [~,p]=ttest(score_algo_NN_on,score_algo_NN_off)
    end
    
      
    end
end

patient_list_temoins = {'Temoin_Laurent_Bis','Temoin_Jean','Temoin_Isabelle','Temoin_Francois','Temoin_6','Temoin_7','Temoin_8','Temoin_9','Temoin_10','Temoin_Sabrina'};
side_movement= {'droite','droite','droite','droite','droite','gauche','droite','droite','droite','droite'};
age_temoins=[55,56,45,72,51,69,65,56,62,45];
side_temoins=[0,0,0,0,0,1,0,0,0,0]; 
sex_temoins =[0,0,1,0,0,0,0,0,0,1];

for k=1:length(patient_list_temoins)
    for z=1:11
    if z==1
        saving_name='decoding_2mvts_individual_trials_all_decimate_window2_perpatient';
    elseif z==2
        saving_name='decoding_2mvts_permutation_trials_C_contra_decimate_window2_perpatient';
    elseif z==3
        saving_name='decoding_2mvts_permutation_trials_all_decimate_window2_perpatient';
    elseif z==4
        saving_name='decoding_2mvts_individual_trials_Fp_F_C_P_contra_decimate_window2_perpatient';
    elseif z==5
        saving_name='decoding_2mvts_individual_trials_Fp_F_C_P_contra_decimate_window2_50trials';
    elseif z==6
        saving_name='decoding_2mvts_permutation_trials_Fp_F_C_P_contra_decimate_window2_perpatient'; % window: 304-1024
    elseif z==7 
        saving_name='decoding_2mvts_individual_trials_EMG_decimate_window_no_jitter_perpatient'; % window: 504-1224
        
%     elseif z==8  
%             saving_name = 'decoding_2mvts_permutation_trials_Fp_F_C_P_contra_decimate_window2_perpatient';
%         elseif z==9
%             saving_name = 'decoding_2mvts_permutation_trials_Fp_F_C_P_contra_decimate_window2_jittered_50ms_perpatient';
%         elseif z==10
%             saving_name = 'decoding_2mvts_permutation_trials_Fp_F_C_P_contra_decimate_window2_jittered_100ms_perpatient';
%         elseif z==11
%             saving_name = 'decoding_2mvts_permutation_trials_Fp_F_C_P_contra_decimate_window2_jittered_200ms_perpatient';    
    
    elseif z==8
        saving_name='decoding_2mvts_permutation_trials_EMG_decimate_window_no_jitter_perpatient';
    elseif z==9
        saving_name='decoding_2mvts_permutation_trials_EMG_decimate_window_jittered_50ms_perpatient';
    elseif z==10
        saving_name='decoding_2mvts_permutation_trials_EMG_decimate_window_jittered_100ms_perpatient';
    elseif z==11
        saving_name='decoding_2mvts_permutation_trials_EMG_decimate_window_jittered_200ms_perpatient';
      end

    
    rootpath = strcat('/Volumes/DBS/DBS/EEG_Analysis/DATA_EEG/TEMOINS/',patient_list_temoins{k},'/Decodage');
    %rootpath=strcat('/Volumes/Charlotte/EEG_Analysis/DATA_EEG/TEMOINS/',patient_list_temoins{k},'/Decodage');
    cd(rootpath)
        load(strcat(rootpath,'/',saving_name,'.mat'))
    temoins_score_algo_LDA_mean(k,z) = mean(score_algo_LDA(:));
    temoins_score_algo_LDA_std(k,z) = std(score_algo_LDA(:)); %./sqrt(625);
    temoins_score_algo_NN_mean(k,z) = mean(score_algo_NN(:));
    temoins_score_algo_NN_std(k,z) = std(score_algo_NN(:)) ;%./sqrt(625);
    temoins_score_algo_LR_mean(k,z) = mean(score_algo_LR(:));
    temoins_score_algo_LR_std(k,z) = std(score_algo_LR(:)) ;%./sqrt(625); 
    
    
    
    end    
end



idx=2
table_stats = [nanmean(Off_Off_score_algo_LDA_mean(:,idx)) nanmean(On_On_score_algo_LDA_mean(:,idx)) nanmean(temoins_score_algo_LDA_mean(:,idx)) ;  nanmean(Off_Off_score_algo_NN_mean(:,idx)) nanmean(On_On_score_algo_NN_mean(:,idx)) nanmean(temoins_score_algo_NN_mean(:,idx)) ; nanmean(Off_Off_score_algo_LR_mean(:,idx)) nanmean(On_On_score_algo_LR_mean(:,idx)) nanmean(temoins_score_algo_LR_mean(:,idx))] 
table_stats = [nanstd(Off_Off_score_algo_LDA_mean(:,idx))/sqrt(length(patient_list)) nanstd(On_On_score_algo_LDA_mean(:,idx))/sqrt(length(patient_list)) nanstd(temoins_score_algo_LDA_mean(:,idx))/sqrt(length(patient_list_temoins)) ;  nanstd(Off_Off_score_algo_NN_mean(:,idx))/sqrt(length(patient_list)) nanstd(On_On_score_algo_NN_mean(:,idx))/sqrt(length(patient_list)) nanstd(temoins_score_algo_NN_mean(:,idx))/sqrt(length(patient_list_temoins)) ; nanstd(Off_Off_score_algo_LR_mean(:,idx))/sqrt(length(patient_list)) nanstd(On_On_score_algo_LR_mean(:,idx))/sqrt(length(patient_list)) nanstd(temoins_score_algo_LR_mean(:,idx))/sqrt(length(patient_list_temoins))]

figure()
cchoice={[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
    
for z=1:8
    subplot(2,4,z); idx=z; hold on; 
    for k=1:length(patient_list)
plot(1:2,[Off_Off_score_algo_LDA_mean(k,idx) On_On_score_algo_LDA_mean(k,idx)],'o-','color',cchoice{k})
    end
scatter(3*ones(10,1),temoins_score_algo_LDA_mean(:,idx),'o','MarkerEdgeColor',[0.8 0.8 0.8],"jitter","on")
errorbar(1:2,[nanmean(Off_Off_score_algo_LDA_mean(:,idx)) nanmean(On_On_score_algo_LDA_mean(:,idx))],[nanstd(Off_Off_score_algo_LDA_mean(:,idx))/sqrt(length(patient_list)) nanstd(On_On_score_algo_LDA_mean(:,idx))/sqrt(length(patient_list))],'ks-','LineWidth',2)
errorbar(3,nanmean(temoins_score_algo_LDA_mean(:,idx)),[nanstd(temoins_score_algo_LDA_mean(:,idx))/sqrt(length(patient_list_temoins))],'ks','LineWidth',2)
for k=1:length(patient_list)
plot(5:6,[Off_Off_score_algo_NN_mean(k,idx) On_On_score_algo_NN_mean(k,idx)],'o-','color',cchoice{k})
end
scatter(7*ones(10,1),temoins_score_algo_NN_mean(:,idx),'o','MarkerEdgeColor',[0.8 0.8 0.8],"jitter","on")
errorbar(5:6,[nanmean(Off_Off_score_algo_NN_mean(:,idx)) nanmean(On_On_score_algo_NN_mean(:,idx))],[nanstd(Off_Off_score_algo_NN_mean(:,idx))/sqrt(length(patient_list)) nanstd(On_On_score_algo_NN_mean(:,idx))/sqrt(length(patient_list))],'k-','LineWidth',2)
errorbar(7,nanmean(temoins_score_algo_NN_mean(:,idx)),[nanstd(temoins_score_algo_NN_mean(:,idx))/sqrt(length(patient_list_temoins))],'ks','LineWidth',2)

for k=1:length(patient_list)
plot(8:9,[Off_Off_score_algo_LR_mean(k,idx) On_On_score_algo_LR_mean(k,idx)],'o-','color',cchoice{k})
end
scatter(10*ones(10,1),temoins_score_algo_LR_mean(:,idx),'o','MarkerEdgeColor',[0.8 0.8 0.8],"jitter","on")
errorbar(8:9,[nanmean(Off_Off_score_algo_LR_mean(:,idx)) nanmean(On_On_score_algo_LR_mean(:,idx))],[nanstd(Off_Off_score_algo_LR_mean(:,idx))/sqrt(length(patient_list)) nanstd(On_On_score_algo_LR_mean(:,idx))/sqrt(length(patient_list))],'k-','LineWidth',2)
errorbar(10,nanmean(temoins_score_algo_LR_mean(:,idx)),[nanstd(temoins_score_algo_LR_mean(:,idx))/sqrt(length(patient_list_temoins))],'ks','LineWidth',2)

xlim([0 11])
plot([0 11],[50 50],'k-')

idx=2
[r,p]=signrank(Off_Off_score_algo_LDA_mean(:,idx),On_On_score_algo_LDA_mean(:,idx));
pval(z,1)=r;
[r,p]=signrank(Off_Off_score_algo_NN_mean(:,idx),On_On_score_algo_NN_mean(:,idx));
pval(z,2)=r;
[r,p]=signrank(Off_Off_score_algo_LR_mean(:,idx),On_On_score_algo_LR_mean(:,idx));
pval(z,3)=r;

idx=2
[r]=ranksum(On_On_score_algo_LDA_mean(:,idx),temoins_score_algo_LDA_mean(:,idx))
[r]=ranksum(On_On_score_algo_NN_mean(:,idx),temoins_score_algo_NN_mean(:,idx))
[r]=ranksum(On_On_score_algo_LR_mean(:,idx),temoins_score_algo_LR_mean(:,idx))


if z==1
    title('Individual All')
elseif z==2
    title('Permutation Contra C')
elseif z==3
    title('Permutation All')
elseif z==4 
    title('Individual Fp/F/C/P Contra')
elseif z==5
    title('Individual Fp/F/C/P Contra - 50 trials')
elseif z==6
    title('Permutation Fp/F/C/P Contra')
end

end


% figure EMG for permutation trials 
figure();
cchoice={[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
idx=8; 
subplot(1,3,1); hold on; 
for k=1:length(patient_list)
    plot(1:2,[Off_Off_score_algo_LDA_mean(k,idx) On_On_score_algo_LDA_mean(k,idx)],'o-','color',cchoice{k})
end
for k=1:10
plot(3:6,[temoins_score_algo_LDA_mean(k,idx) temoins_score_algo_LDA_mean(k,idx+1) temoins_score_algo_LDA_mean(k,idx+2) temoins_score_algo_LDA_mean(k,idx+3)],'o-','color',[0.8 0.8 0.8])
end 
errorbar(1:2,[nanmean(Off_Off_score_algo_LDA_mean(:,idx)) nanmean(On_On_score_algo_LDA_mean(:,idx))],[nanstd(Off_Off_score_algo_LDA_mean(:,idx))/sqrt(length(patient_list)) nanstd(On_On_score_algo_LDA_mean(:,idx))/sqrt(length(patient_list))],'ks-','LineWidth',2)
errorbar(3:6,[nanmean(temoins_score_algo_LDA_mean(:,idx)) nanmean(temoins_score_algo_LDA_mean(:,idx+1)) nanmean(temoins_score_algo_LDA_mean(:,idx+2)) nanmean(temoins_score_algo_LDA_mean(:,idx+3))],[nanstd(temoins_score_algo_LDA_mean(:,idx))/sqrt(length(patient_list_temoins)) nanstd(temoins_score_algo_LDA_mean(:,idx+1))/sqrt(length(patient_list_temoins)) nanstd(temoins_score_algo_LDA_mean(:,idx+2))/sqrt(length(patient_list_temoins)) nanstd(temoins_score_algo_LDA_mean(:,idx+3))/sqrt(length(patient_list_temoins))],'k-','LineWidth',2)
set(gca,'xtick',1:6,'xticklabel',{'Off Off','On On','No jitter','50 ms','100m ms', '200 ms'})
xtickangle(30)
xlim([0 7])
ylim([40 100])
%ylim([30 90])
plot([0 7],[50 50],'k-')
title('LDA')
subplot(1,3,2); hold on ; 
for k=1:length(patient_list)
    plot(1:2,[Off_Off_score_algo_NN_mean(k,idx) On_On_score_algo_NN_mean(k,idx)],'o-','color',cchoice{k})
end
for k=1:10
plot(3:6,[temoins_score_algo_NN_mean(k,idx) temoins_score_algo_NN_mean(k,idx+1) temoins_score_algo_NN_mean(k,idx+2) temoins_score_algo_NN_mean(k,idx+3)],'o-','color',[0.8 0.8 0.8])
end
errorbar(1:2,[nanmean(Off_Off_score_algo_NN_mean(:,idx)) nanmean(On_On_score_algo_NN_mean(:,idx))],[nanstd(Off_Off_score_algo_NN_mean(:,idx))/sqrt(length(patient_list)) nanstd(On_On_score_algo_NN_mean(:,idx))/sqrt(length(patient_list))],'k-','LineWidth',2)
errorbar(3:6,[nanmean(temoins_score_algo_NN_mean(:,idx)) nanmean(temoins_score_algo_NN_mean(:,idx+1)) nanmean(temoins_score_algo_NN_mean(:,idx+2)) nanmean(temoins_score_algo_NN_mean(:,idx+3))],[nanstd(temoins_score_algo_NN_mean(:,idx))/sqrt(length(patient_list_temoins)) nanstd(temoins_score_algo_NN_mean(:,idx+1))/sqrt(length(patient_list_temoins)) nanstd(temoins_score_algo_NN_mean(:,idx+2))/sqrt(length(patient_list_temoins)) nanstd(temoins_score_algo_NN_mean(:,idx+3))/sqrt(length(patient_list_temoins))],'k-','LineWidth',2)
xlim([0 7])
ylim([40 100])
%ylim([30 90])
plot([0 7],[50 50],'k-')
title('NN')
xtickangle(30)
set(gca,'xtick',1:6,'xticklabel',{'Off Off','On On','No jitter','50 ms','100m ms', '200 ms'})
subplot(1,3,3); hold on; 
for k=1:length(patient_list)
    plot(1:2,[Off_Off_score_algo_LR_mean(k,idx) On_On_score_algo_LR_mean(k,idx)],'o-','color',cchoice{k})
end
for k=1:10
plot(3:6,[temoins_score_algo_LR_mean(k,idx) temoins_score_algo_LR_mean(k,idx+1) temoins_score_algo_LR_mean(k,idx+2) temoins_score_algo_LR_mean(k,idx+3)],'o-','color',[0.8 0.8 0.8])
end
errorbar(1:2,[nanmean(Off_Off_score_algo_LR_mean(:,idx)) nanmean(On_On_score_algo_LR_mean(:,idx))],[nanstd(Off_Off_score_algo_LR_mean(:,idx))/sqrt(length(patient_list)) nanstd(On_On_score_algo_LR_mean(:,idx))/sqrt(length(patient_list))],'k-','LineWidth',2)
errorbar(3:6,[nanmean(temoins_score_algo_LR_mean(:,idx)) nanmean(temoins_score_algo_LR_mean(:,idx+1)) nanmean(temoins_score_algo_LR_mean(:,idx+2)) nanmean(temoins_score_algo_LR_mean(:,idx+3))],[nanstd(temoins_score_algo_LR_mean(:,idx))/sqrt(length(patient_list_temoins)) nanstd(temoins_score_algo_LR_mean(:,idx+1))/sqrt(length(patient_list_temoins)) nanstd(temoins_score_algo_LR_mean(:,idx+2))/sqrt(length(patient_list_temoins)) nanstd(temoins_score_algo_LR_mean(:,idx+3))/sqrt(length(patient_list_temoins))],'k-','LineWidth',2)
xlim([0 7])
ylim([40 100])
%ylim([30 90])
title('LR')
set(gca,'xtick',1:6,'xticklabel',{'Off Off','On On','No jitter','50 ms','100m ms', '200 ms'})
xtickangle(30)
plot([0 7],[50 50],'k-')
sgtitle('Fp/F/C/P contra electrodes Permutation trials')
    



figure(); 
subplot(2,2,1); hold on; 
idx1=8; idx2=6;
scatter(Off_Off_score_algo_LDA_mean(:,idx1),Off_Off_score_algo_LDA_mean(:,idx2),'r')
scatter(On_On_score_algo_LDA_mean(:,idx1),On_On_score_algo_LDA_mean(:,idx2),'b')
scatter(temoins_score_algo_LDA_mean(:,idx1),temoins_score_algo_LDA_mean(:,idx2),'k')
ylabel('Decoding EEG electrodes')
xlabel('Decoding EMG electrodes')
legend('Off Off','On On','Temoins')
title('LDA, EMG permutation')
[r,p]=corrcoef(temoins_score_algo_LDA_mean(:,idx1),temoins_score_algo_LDA_mean(:,idx2))  %p=0.78; r=0.10
[r,p]=corrcoef(On_On_score_algo_LDA_mean(:,idx1),On_On_score_algo_LDA_mean(:,idx2))  %p=0.63
[r,p]=corrcoef(Off_Off_score_algo_LDA_mean(:,idx1),Off_Off_score_algo_LDA_mean(:,idx2)) %p=0.0755; r=0.70

subplot(2,2,2); hold on; 
idx1=7; idx2=6;
scatter(Off_Off_score_algo_LDA_mean(:,idx1),Off_Off_score_algo_LDA_mean(:,idx2),'r')
scatter(On_On_score_algo_LDA_mean(:,idx1),On_On_score_algo_LDA_mean(:,idx2),'b')
scatter(temoins_score_algo_LDA_mean(:,idx1),temoins_score_algo_LDA_mean(:,idx2),'k')
ylabel('Decoding EEG electrodes')
xlabel('Decoding EMG electrodes')
legend('Off Off','On On','Temoins')
title('LDA, EMG individual')
[r,p]=corrcoef(temoins_score_algo_LDA_mean(:,idx1),temoins_score_algo_LDA_mean(:,idx2)) ; %p=0.27; r=0.38
[r,p]=corrcoef(On_On_score_algo_LDA_mean(:,idx1),On_On_score_algo_LDA_mean(:,idx2)) ; %p=0.735
[r,p]=corrcoef(Off_Off_score_algo_LDA_mean(:,idx1),Off_Off_score_algo_LDA_mean(:,idx2)) ; %p=0.1440; r=0.61

subplot(2,2,3);hold on ;
idx1=8; idx2=6;
scatter(Off_Off_score_algo_LR_mean(:,idx1),Off_Off_score_algo_LR_mean(:,idx2),'r')
scatter(On_On_score_algo_LR_mean(:,idx1),On_On_score_algo_LR_mean(:,idx2),'b')
scatter(temoins_score_algo_LR_mean(:,idx1),temoins_score_algo_LR_mean(:,idx2),'k')
ylabel('Decoding EEG electrodes')
xlabel('Decoding EMG electrodes')
legend('Off Off','On On','Temoins')
title('LR, EMG permutation')

subplot(2,2,4); hold on; 
idx1=8; idx2=6;
scatter(Off_Off_score_algo_NN_mean(:,idx1),Off_Off_score_algo_NN_mean(:,idx2),'r')
scatter(On_On_score_algo_NN_mean(:,idx1),On_On_score_algo_NN_mean(:,idx2),'b')
scatter(temoins_score_algo_NN_mean(:,idx1),temoins_score_algo_NN_mean(:,idx2),'k')
ylabel('Decoding EEG electrodes')
xlabel('Decoding EMG electrodes')
legend('Off Off','On On','Temoins')
title('NN, EMG permutation')





% dependency on number of trials used % for 2 MVTS
temoins_individual = [224,220,190,140,215,213,179,226,184,203];
temoins_permut = [37,36,31,23,35,35,29,37,30,33];
figure(); 
subplot(1,2,1);hold on; 
scatter(temoins_individual,temoins_score_algo_LR_mean(:,1))
scatter(temoins_individual,temoins_score_algo_LR_mean(:,3))
subplot(1,2,2); hold on; 
scatter(temoins_permut,temoins_score_algo_LR_mean(:,2))
scatter(temoins_permut,temoins_score_algo_LR_mean(:,4))
xlabel('Number of training trials')
ylabel('Accuracy LR')

% DBS
minimal_number_individual=[167,179,96,105,198,50,75]; % per patient for 2 mvts
minimal_number_permutation=[27,29,16,17,33,8,12]; % per patient for 2 mvts
figure(); 
subplot(1,2,1);hold on; 
scatter(minimal_number_individual,Off_Off_score_algo_LR_mean(:,1))
scatter(minimal_number_individual,On_On_score_algo_LR_mean(:,1))
ylabel('Accuracy LR')
xlabel('Number of training trials')

figure(); 
subplot(1,2,1); hold on; 
scatter(minimal_number_permutation,Off_Off_score_algo_LR_mean(:,6))
scatter(minimal_number_permutation,On_On_score_algo_LR_mean(:,6))
subplot(1,2,2); hold on; 
scatter(minimal_number_permutation,Off_Off_score_algo_LR_mean(:,3))
scatter(minimal_number_permutation,On_On_score_algo_LR_mean(:,3))
xlabel('Number of training trials')
ylabel('Accuracy LR')
[r,p]=corrcoef(minimal_number_permutation,On_On_score_algo_LR_mean(:,6))

figure(); 
subplot(1,2,1); hold on; 
scatter(minimal_number_permutation,Off_Off_score_algo_LDA_mean(:,6))
scatter(minimal_number_permutation,On_On_score_algo_LDA_mean(:,6))
ylim([30 90])
ylabel('Accuracy LDA')
subplot(1,2,2); hold on; 
scatter(minimal_number_permutation,Off_Off_score_algo_LDA_mean(:,3))
scatter(minimal_number_permutation,On_On_score_algo_LDA_mean(:,3))
xlabel('Number of training trials')
ylabel('Accuracy LDA')
ylim([30 90])
[r,p]=corrcoef(minimal_number_permutation,On_On_score_algo_LDA_mean(:,6))



%% Pour trois mouvements: 
clear ; % ici c'est DBS 5 qui descend à chque fois

patient_list = {'DBS_2','DBS_3','DBS_4','DBS_5','DBS_6','DBS_7'};
side_movement={'droite','gauche','droite','droite','droite','droite'}
minimal_number_individual=[179,62,105,198,46,75]; % per patient for 3 mvts
minimal_number_permutation=[29,10,17,33,7,12] ;

patient_list_temoins = {'Temoin_Laurent_Bis','Temoin_Jean','Temoin_Isabelle','Temoin_Francois','Temoin_6','Temoin_7','Temoin_8','Temoin_9','Temoin_10','Temoin_Sabrina'};
side_movement= {'droite','droite','droite','droite','droite','gauche','droite','droite','droite','droite'};
minimal_individual_per_patient=[224,169,180,140,199,180,112,194,174,198]; %3 mvts
minimal_permutation_per_patient=[37,28,30,23,33,30,18,32,29,33]; %3 mvts


for k=1:length(patient_list)
    
    for z=1:4
        
        if z==1
            saving_name='decoding_3mvts_individual_trials_all_decimate_window2_perpatient';
        elseif z==2
            saving_name='decoding_3mvts_permutation_trials_all_decimate_window2_perpatient';
            
        elseif z==3
            saving_name='decoding_3mvts_individual_trials_Fp_F_C_P_contra_decimate_window2_perpatient';
        elseif z==4
            saving_name='decoding_3mvts_permutation_trials_Fp_F_C_P_contra_decimate_window2_perpatient';
        end
        
        rootpath=strcat('/Volumes/Charlotte/EEG_Analysis/DATA_EEG/DBS/',patient_list{k},'/DBS_Off_Dopa_Off/Decodage');
        cd(rootpath)
        try
            load(strcat(rootpath,'/',saving_name,'.mat'))
            Off_Off_score_algo_LDA_mean(k,z) = mean(score_algo_LDA(:));
            Off_Off_score_algo_LDA_std(k,z) = std(score_algo_LDA(:)); %./sqrt(625);
            Off_Off_score_algo_NN_mean(k,z) = mean(score_algo_NN(:));
            Off_Off_score_algo_NN_std(k,z) = std(score_algo_NN(:)) ;%./sqrt(625);
            Off_Off_score_algo_LR_mean(k,z) = mean(score_algo_LR(:));
            Off_Off_score_algo_LR_std(k,z) = std(score_algo_LR(:)) ;%./sqrt(625);
        catch
            Off_Off_score_algo_LDA_mean(k,z)=NaN;
            Off_Off_score_algo_LDA_std(k,z)=NaN;
            Off_Off_score_algo_NN_mean(k,z) = NaN;
            Off_Off_score_algo_NN_std(k,z) = NaN;
            Off_Off_score_algo_LR_mean(k,z) = NaN;
            Off_Off_score_algo_LR_std(k,z)=NaN;
        end
        
        score_algo_LDA_off = score_algo_LDA(:);
        score_algo_LR_off = score_algo_LR(:);
        score_algo_NN_off = score_algo_NN(:);
        
        
        rootpath=strcat('/Volumes/Charlotte/EEG_Analysis/DATA_EEG/DBS/',patient_list{k},'/DBS_On_Dopa_On/Decodage');
        cd(rootpath)
        load(strcat(rootpath,'/',saving_name,'.mat'))
        On_On_score_algo_LDA_mean(k,z) = mean(score_algo_LDA(:));
        On_On_score_algo_LDA_std(k,z) = std(score_algo_LDA(:)); %./sqrt(625);
        On_On_score_algo_NN_mean(k,z) = mean(score_algo_NN(:));
        On_On_score_algo_NN_std(k,z) = std(score_algo_NN(:)) ;%./sqrt(625);
        On_On_score_algo_LR_mean(k,z) = mean(score_algo_LR(:));
        On_On_score_algo_LR_std(k,z) = std(score_algo_LR(:)) ;%./sqrt(625);
        score_algo_LDA_on = score_algo_LDA(:);
        score_algo_LR_on = score_algo_LR(:);
        score_algo_NN_on = score_algo_NN(:);
        
        
        if z==4
            disp(k)
            disp(mean(score_algo_LDA_on)-mean(score_algo_LDA_off))
            [~,p]=ttest2(score_algo_LDA_on,score_algo_LDA_off)
            disp(mean(score_algo_LR_on)-mean(score_algo_LR_off))
            [~,p]=ttest2(score_algo_LR_on,score_algo_LR_off)
            disp(mean(score_algo_NN_on)-mean(score_algo_NN_off))
            [~,p]=ttest2(score_algo_NN_on,score_algo_NN_off)
        end
        
        
        
        
    end
end

patient_list_temoins = {'Temoin_Laurent_Bis','Temoin_Jean','Temoin_Isabelle','Temoin_Francois','Temoin_6','Temoin_7','Temoin_8','Temoin_9','Temoin_10','Temoin_Sabrina'};
side_movement= {'droite','droite','droite','droite','droite','gauche','droite','droite','droite','droite'};

for k=1:length(patient_list_temoins)
    for z=1:4
    if z==1
        saving_name='decoding_3mvts_individual_trials_all_decimate_window2_perpatient';
    elseif z==2
        saving_name='decoding_3mvts_permutation_trials_all_decimate_window2_perpatient';
    elseif z==3
        saving_name='decoding_3mvts_individual_trials_Fp_F_C_P_contra_decimate_window2_perpatient';
    elseif z==4
        saving_name='decoding_3mvts_permutation_trials_Fp_F_C_P_contra_decimate_window2_perpatient';
    end
    
    rootpath=strcat('/Volumes/Charlotte/EEG_Analysis/DATA_EEG/TEMOINS/',patient_list_temoins{k},'/Decodage');
    cd(rootpath)
        load(strcat(rootpath,'/',saving_name,'.mat'))
    temoins_score_algo_LDA_mean(k,z) = mean(score_algo_LDA(:));
    temoins_score_algo_LDA_std(k,z) = std(score_algo_LDA(:)); %./sqrt(625);
    temoins_score_algo_NN_mean(k,z) = mean(score_algo_NN(:));
    temoins_score_algo_NN_std(k,z) = std(score_algo_NN(:)) ;%./sqrt(625);
    temoins_score_algo_LR_mean(k,z) = mean(score_algo_LR(:));
    temoins_score_algo_LR_std(k,z) = std(score_algo_LR(:)) ;%./sqrt(625);  
    end    
end


figure()
cchoice={[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
    
for z=1:4
    subplot(2,2,z); idx=z; hold on; 
    for k=1:length(patient_list)
plot(1:2,[Off_Off_score_algo_LDA_mean(k,idx) On_On_score_algo_LDA_mean(k,idx)],'o-','color',cchoice{k})
    end
scatter(3*ones(10,1),temoins_score_algo_LDA_mean(:,idx),'o','MarkerEdgeColor',[0.8 0.8 0.8],"jitter","on")
errorbar(1:2,[nanmean(Off_Off_score_algo_LDA_mean(:,idx)) nanmean(On_On_score_algo_LDA_mean(:,idx))],[nanstd(Off_Off_score_algo_LDA_mean(:,idx))/sqrt(length(patient_list)) nanstd(On_On_score_algo_LDA_mean(:,idx))/sqrt(length(patient_list))],'ks-','LineWidth',2)
errorbar(3,nanmean(temoins_score_algo_LDA_mean(:,idx)),[nanstd(temoins_score_algo_LDA_mean(:,idx))/sqrt(length(patient_list_temoins))],'ks','LineWidth',2)

for k=1:length(patient_list)
plot(5:6,[Off_Off_score_algo_NN_mean(k,idx) On_On_score_algo_NN_mean(k,idx)],'o-','color',cchoice{k})
end
scatter(7*ones(10,1),temoins_score_algo_NN_mean(:,idx),'o','MarkerEdgeColor',[0.8 0.8 0.8],"jitter","on")
errorbar(5:6,[nanmean(Off_Off_score_algo_NN_mean(:,idx)) nanmean(On_On_score_algo_NN_mean(:,idx))],[nanstd(Off_Off_score_algo_NN_mean(:,idx))/sqrt(length(patient_list)) nanstd(On_On_score_algo_NN_mean(:,idx))/sqrt(length(patient_list))],'k-','LineWidth',2)
errorbar(7,nanmean(temoins_score_algo_NN_mean(:,idx)),[nanstd(temoins_score_algo_NN_mean(:,idx))/sqrt(length(patient_list_temoins))],'ks','LineWidth',2)

for k=1:length(patient_list)
plot(8:9,[Off_Off_score_algo_LR_mean(k,idx) On_On_score_algo_LR_mean(k,idx)],'o-','color',cchoice{k})
end
scatter(10*ones(10,1),temoins_score_algo_LR_mean(:,idx),'o','MarkerEdgeColor',[0.8 0.8 0.8],"jitter","on")
errorbar(8:9,[nanmean(Off_Off_score_algo_LR_mean(:,idx)) nanmean(On_On_score_algo_LR_mean(:,idx))],[nanstd(Off_Off_score_algo_LR_mean(:,idx))/sqrt(length(patient_list)) nanstd(On_On_score_algo_LR_mean(:,idx))/sqrt(length(patient_list))],'k-','LineWidth',2)
errorbar(10,nanmean(temoins_score_algo_LR_mean(:,idx)),[nanstd(temoins_score_algo_LR_mean(:,idx))/sqrt(length(patient_list_temoins))],'ks','LineWidth',2)

xlim([0 11])
plot([0 11],[33 33],'k-')

[r,p]=ttest2(Off_Off_score_algo_LDA_mean(:,idx),On_On_score_algo_LDA_mean(:,idx));
pval(z,1)=p;
[r,p]=ttest2(Off_Off_score_algo_NN_mean(:,idx),On_On_score_algo_NN_mean(:,idx));
pval(z,2)=p;
[r,p]=ttest2(Off_Off_score_algo_LR_mean(:,idx),On_On_score_algo_LR_mean(:,idx));
pval(z,3)=p;

if z==1
    title('Individual All')
elseif z==2
    title('Permutation All')
elseif z==3 
    title('Individual Fp/F/C/P Contra')
elseif z==4
    title('Permutation Fp/F/C/P Contra')
end

end


idx=2
[r]=ranksum(Off_Off_score_algo_LDA_mean(:,idx),temoins_score_algo_LDA_mean(:,idx))
[r]=ranksum(Off_Off_score_algo_NN_mean(:,idx),temoins_score_algo_NN_mean(:,idx))
[r]=ranksum(Off_Off_score_algo_LR_mean(:,idx),temoins_score_algo_LR_mean(:,idx))

idx=2
[r,p]=signrank(Off_Off_score_algo_LDA_mean(:,idx),On_On_score_algo_LDA_mean(:,idx));
pval(z,1)=r;
[r,p]=signrank(Off_Off_score_algo_NN_mean(:,idx),On_On_score_algo_NN_mean(:,idx));
pval(z,2)=r;
[r,p]=signrank(Off_Off_score_algo_LR_mean(:,idx),On_On_score_algo_LR_mean(:,idx));
pval(z,3)=r;




%% When considering the full fou conditions: 


patient_list = {'DBS_1','DBS_2','DBS_3','DBS_4','DBS_5','DBS_7'};
side_movement={'gauche','droite','gauche','droite','droite','droite'}

minimal_number_individual=[167,179,62,105,198,75]; %per patient for 2 mvts
minimal_number_permutation=[27,29,10,17,33,12]; %per patient for 2 mvts

for k=1:length(patient_list)
    
    for z=1:4
        
        if z==1
            saving_name='decoding_2mvts_individual_trials_all_decimate_window2_perpatient';
        elseif z==2
            saving_name='decoding_2mvts_permutation_trials_all_decimate_window2_perpatient';
            
        elseif z==3
            saving_name='decoding_2mvts_individual_trials_Fp_F_C_P_contra_decimate_window2_perpatient';
        elseif z==4
            saving_name='decoding_2mvts_permutation_trials_Fp_F_C_P_contra_decimate_window2_perpatient';
        end
        
        rootpath=strcat('/Volumes/Charlotte/EEG_Analysis/DATA_EEG/DBS/',patient_list{k},'/DBS_Off_Dopa_Off/Decodage');
        cd(rootpath)
        try
            load(strcat(rootpath,'/',saving_name,'.mat'))
            Off_Off_score_algo_LDA_mean(k,z) = mean(score_algo_LDA(:));
            Off_Off_score_algo_LDA_std(k,z) = std(score_algo_LDA(:)); %./sqrt(625);
            Off_Off_score_algo_NN_mean(k,z) = mean(score_algo_NN(:));
            Off_Off_score_algo_NN_std(k,z) = std(score_algo_NN(:)) ;%./sqrt(625);
            Off_Off_score_algo_LR_mean(k,z) = mean(score_algo_LR(:));
            Off_Off_score_algo_LR_std(k,z) = std(score_algo_LR(:)) ;%./sqrt(625);
        catch
            Off_Off_score_algo_LDA_mean(k,z)=NaN;
            Off_Off_score_algo_LDA_std(k,z)=NaN;
            Off_Off_score_algo_NN_mean(k,z) = NaN;
            Off_Off_score_algo_NN_std(k,z) = NaN;
            Off_Off_score_algo_LR_mean(k,z) = NaN;
            Off_Off_score_algo_LR_std(k,z)=NaN;
        end
        
        rootpath=strcat('/Volumes/Charlotte/EEG_Analysis/DATA_EEG/DBS/',patient_list{k},'/DBS_Off_Dopa_On/Decodage');
        cd(rootpath)
        try
            load(strcat(rootpath,'/',saving_name,'.mat'))
            Off_On_score_algo_LDA_mean(k,z) = mean(score_algo_LDA(:));
            Off_On_score_algo_LDA_std(k,z) = std(score_algo_LDA(:)); %./sqrt(625);
            Off_On_score_algo_NN_mean(k,z) = mean(score_algo_NN(:));
            Off_On_score_algo_NN_std(k,z) = std(score_algo_NN(:)) ;%./sqrt(625);
            Off_On_score_algo_LR_mean(k,z) = mean(score_algo_LR(:));
            Off_On_score_algo_LR_std(k,z) = std(score_algo_LR(:)) ;%./sqrt(625);
            
        catch
            Off_On_score_algo_LDA_mean(k,z)=NaN;
            Off_On_score_algo_LDA_std(k,z)=NaN;
            Off_On_score_algo_NN_mean(k,z) = NaN;
            Off_On_score_algo_NN_std(k,z) = NaN;
            Off_On_score_algo_LR_mean(k,z) = NaN;
            Off_On_score_algo_LR_std(k,z)=NaN;
            
        end
        
        rootpath=strcat('/Volumes/Charlotte/EEG_Analysis/DATA_EEG/DBS/',patient_list{k},'/DBS_On_Dopa_Off/Decodage');
        cd(rootpath)
        load(strcat(rootpath,'/',saving_name,'.mat'))
        On_Off_score_algo_LDA_mean(k,z) = mean(score_algo_LDA(:));
        On_Off_score_algo_LDA_std(k,z) = std(score_algo_LDA(:)); %./sqrt(625);
        On_Off_score_algo_NN_mean(k,z) = mean(score_algo_NN(:));
        On_Off_score_algo_NN_std(k,z) = std(score_algo_NN(:)) ;%./sqrt(625);
        On_Off_score_algo_LR_mean(k,z) = mean(score_algo_LR(:));
        On_Off_score_algo_LR_std(k,z) = std(score_algo_LR(:)) ;%./sqrt(625);
        
        
        
        rootpath=strcat('/Volumes/Charlotte/EEG_Analysis/DATA_EEG/DBS/',patient_list{k},'/DBS_On_Dopa_On/Decodage');
        cd(rootpath)
        load(strcat(rootpath,'/',saving_name,'.mat'))
        On_On_score_algo_LDA_mean(k,z) = mean(score_algo_LDA(:));
        On_On_score_algo_LDA_std(k,z) = std(score_algo_LDA(:)); %./sqrt(625);
        On_On_score_algo_NN_mean(k,z) = mean(score_algo_NN(:));
        On_On_score_algo_NN_std(k,z) = std(score_algo_NN(:)) ;%./sqrt(625);
        On_On_score_algo_LR_mean(k,z) = mean(score_algo_LR(:));
        On_On_score_algo_LR_std(k,z) = std(score_algo_LR(:)) ;%./sqrt(625);
        
    end
end

patient_list_temoins = {'Temoin_Laurent_Bis','Temoin_Jean','Temoin_Isabelle','Temoin_Francois','Temoin_6','Temoin_7','Temoin_8','Temoin_9','Temoin_10','Temoin_Sabrina'};
side_movement= {'droite','droite','droite','droite','droite','gauche','droite','droite','droite','droite'};

for k=1:length(patient_list_temoins)
    for z=1:4
    if z==1
        saving_name='decoding_2mvts_individual_trials_all_decimate_window2_perpatient';
    elseif z==2
        saving_name='decoding_2mvts_permutation_trials_all_decimate_window2_perpatient';
    elseif z==3
        saving_name='decoding_2mvts_individual_trials_Fp_F_C_P_contra_decimate_window2_perpatient';
    elseif z==4
        saving_name='decoding_2mvts_permutation_trials_Fp_F_C_P_contra_decimate_window2_perpatient';
    end
    
    rootpath=strcat('/Volumes/Charlotte/EEG_Analysis/DATA_EEG/TEMOINS/',patient_list_temoins{k},'/Decodage');
    cd(rootpath)
        load(strcat(rootpath,'/',saving_name,'.mat'))
    temoins_score_algo_LDA_mean(k,z) = mean(score_algo_LDA(:));
    temoins_score_algo_LDA_std(k,z) = std(score_algo_LDA(:)); %./sqrt(625);
    temoins_score_algo_NN_mean(k,z) = mean(score_algo_NN(:));
    temoins_score_algo_NN_std(k,z) = std(score_algo_NN(:)) ;%./sqrt(625);
    temoins_score_algo_LR_mean(k,z) = mean(score_algo_LR(:));
    temoins_score_algo_LR_std(k,z) = std(score_algo_LR(:)) ;%./sqrt(625);  
    end    
end


figure()
cchoice={[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
    
for z=1:4
    subplot(2,2,z); idx=z; hold on; 
    for k=1:length(patient_list)
plot(1:4,[Off_Off_score_algo_LDA_mean(k,idx) On_Off_score_algo_LDA_mean(k,idx) Off_On_score_algo_LDA_mean(k,idx)  On_On_score_algo_LDA_mean(k,idx)],'o-','color',cchoice{k})
    end
scatter(5*ones(10,1),temoins_score_algo_LDA_mean(:,idx),'o','MarkerEdgeColor',[0.8 0.8 0.8],"jitter","on")
errorbar(1:4,[nanmean(Off_Off_score_algo_LDA_mean(:,idx)) nanmean(On_Off_score_algo_LDA_mean(:,idx)) nanmean(Off_On_score_algo_LDA_mean(:,idx)) nanmean(On_On_score_algo_LDA_mean(:,idx))],[nanstd(Off_Off_score_algo_LDA_mean(:,idx))/sqrt(length(patient_list)) nanstd(On_Off_score_algo_LDA_mean(:,idx))/sqrt(length(patient_list)) nanstd(Off_On_score_algo_LDA_mean(:,idx))/sqrt(length(patient_list)) nanstd(On_On_score_algo_LDA_mean(:,idx))/sqrt(length(patient_list))],'ks-','LineWidth',2)
errorbar(5,nanmean(temoins_score_algo_LDA_mean(:,idx)),[nanstd(temoins_score_algo_LDA_mean(:,idx))/sqrt(length(patient_list_temoins))],'ks','LineWidth',2)

for k=1:length(patient_list)
    plot(7:10,[Off_Off_score_algo_NN_mean(k,idx) On_Off_score_algo_NN_mean(k,idx) Off_On_score_algo_NN_mean(k,idx)  On_On_score_algo_NN_mean(k,idx)],'o-','color',cchoice{k})
end
scatter(11*ones(10,1),temoins_score_algo_NN_mean(:,idx),'o','MarkerEdgeColor',[0.8 0.8 0.8],"jitter","on")
errorbar(7:10,[nanmean(Off_Off_score_algo_NN_mean(:,idx)) nanmean(On_Off_score_algo_NN_mean(:,idx)) nanmean(Off_On_score_algo_NN_mean(:,idx)) nanmean(On_On_score_algo_NN_mean(:,idx))],[nanstd(Off_Off_score_algo_NN_mean(:,idx))/sqrt(length(patient_list)) nanstd(On_Off_score_algo_NN_mean(:,idx))/sqrt(length(patient_list)) nanstd(Off_On_score_algo_NN_mean(:,idx))/sqrt(length(patient_list)) nanstd(On_On_score_algo_NN_mean(:,idx))/sqrt(length(patient_list))],'ks-','LineWidth',2)
errorbar(11,nanmean(temoins_score_algo_NN_mean(:,idx)),[nanstd(temoins_score_algo_NN_mean(:,idx))/sqrt(length(patient_list_temoins))],'ks','LineWidth',2)

for k=1:length(patient_list)
plot(13:16,[Off_Off_score_algo_LR_mean(k,idx) On_Off_score_algo_LR_mean(k,idx) Off_On_score_algo_LR_mean(k,idx)  On_On_score_algo_LR_mean(k,idx)],'o-','color',cchoice{k})
end
scatter(17*ones(10,1),temoins_score_algo_LR_mean(:,idx),'o','MarkerEdgeColor',[0.8 0.8 0.8],"jitter","on")
errorbar(13:16,[nanmean(Off_Off_score_algo_LR_mean(:,idx)) nanmean(On_Off_score_algo_LR_mean(:,idx)) nanmean(Off_On_score_algo_LR_mean(:,idx)) nanmean(On_On_score_algo_LR_mean(:,idx))],[nanstd(Off_Off_score_algo_LR_mean(:,idx))/sqrt(length(patient_list)) nanstd(On_Off_score_algo_LR_mean(:,idx))/sqrt(length(patient_list)) nanstd(Off_On_score_algo_LR_mean(:,idx))/sqrt(length(patient_list)) nanstd(On_On_score_algo_LR_mean(:,idx))/sqrt(length(patient_list))],'ks-','LineWidth',2)
errorbar(17,nanmean(temoins_score_algo_LR_mean(:,idx)),[nanstd(temoins_score_algo_LR_mean(:,idx))/sqrt(length(patient_list_temoins))],'ks','LineWidth',2)

xlim([0 18])
plot([0 18],[50 50],'k-')

[r,p]=ttest2(Off_Off_score_algo_LDA_mean(:,idx),On_On_score_algo_LDA_mean(:,idx));
pval(z,1)=p;
[r,p]=ttest2(Off_Off_score_algo_NN_mean(:,idx),On_On_score_algo_NN_mean(:,idx));
pval(z,2)=p;
[r,p]=ttest2(Off_Off_score_algo_LR_mean(:,idx),On_On_score_algo_LR_mean(:,idx));
pval(z,3)=p;

if z==1
    title('Individual All')
elseif z==2
    title('Permutation All')
elseif z==3 
    title('Individual Fp/F/C/P Contra')
elseif z==4
    title('Permutation Fp/F/C/P Contra')
end

end







figure()
cchoice={[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
    
for z=1:4
    subplot(2,2,z); idx=z; hold on; 
    for k=1:length(patient_list)
plot(1:4,[Off_Off_score_algo_LDA_mean(k,idx)/Off_Off_score_algo_LDA_mean(k,idx) On_Off_score_algo_LDA_mean(k,idx)/Off_Off_score_algo_LDA_mean(k,idx) Off_On_score_algo_LDA_mean(k,idx)/Off_Off_score_algo_LDA_mean(k,idx)  On_On_score_algo_LDA_mean(k,idx)/Off_Off_score_algo_LDA_mean(k,idx)],'o-','color',cchoice{k})
    end
errorbar(1:4,[nanmean(Off_Off_score_algo_LDA_mean(:,idx)./Off_Off_score_algo_LDA_mean(:,idx)) nanmean(On_Off_score_algo_LDA_mean(:,idx)./Off_Off_score_algo_LDA_mean(:,idx)) nanmean(Off_On_score_algo_LDA_mean(:,idx)./Off_Off_score_algo_LDA_mean(:,idx)) nanmean(On_On_score_algo_LDA_mean(:,idx)./Off_Off_score_algo_LDA_mean(:,idx))],[nanstd(Off_Off_score_algo_LDA_mean(:,idx)./Off_Off_score_algo_LDA_mean(:,idx))/sqrt(length(patient_list)) nanstd(On_Off_score_algo_LDA_mean(:,idx)./Off_Off_score_algo_LDA_mean(:,idx))/sqrt(length(patient_list)) nanstd(Off_On_score_algo_LDA_mean(:,idx)./Off_Off_score_algo_LDA_mean(:,idx))/sqrt(length(patient_list)) nanstd(On_On_score_algo_LDA_mean(:,idx)./Off_Off_score_algo_LDA_mean(:,idx))/sqrt(length(patient_list))],'ks-','LineWidth',2)

for k=1:length(patient_list)
    plot(6:9,[Off_Off_score_algo_NN_mean(k,idx)/Off_Off_score_algo_LDA_mean(k,idx) On_Off_score_algo_NN_mean(k,idx)/Off_Off_score_algo_LDA_mean(k,idx) Off_On_score_algo_NN_mean(k,idx)/Off_Off_score_algo_LDA_mean(k,idx)  On_On_score_algo_NN_mean(k,idx)/Off_Off_score_algo_LDA_mean(k,idx)],'o-','color',cchoice{k})
end
errorbar(6:9,[nanmean(Off_Off_score_algo_NN_mean(:,idx)./Off_Off_score_algo_LDA_mean(:,idx)) nanmean(On_Off_score_algo_NN_mean(:,idx)./Off_Off_score_algo_LDA_mean(:,idx)) nanmean(Off_On_score_algo_NN_mean(:,idx)./Off_Off_score_algo_LDA_mean(:,idx)) nanmean(On_On_score_algo_NN_mean(:,idx)./Off_Off_score_algo_LDA_mean(:,idx))],[nanstd(Off_Off_score_algo_NN_mean(:,idx)./Off_Off_score_algo_LDA_mean(:,idx))/sqrt(length(patient_list)) nanstd(On_Off_score_algo_NN_mean(:,idx)./Off_Off_score_algo_LDA_mean(:,idx))/sqrt(length(patient_list)) nanstd(Off_On_score_algo_NN_mean(:,idx)./Off_Off_score_algo_LDA_mean(:,idx))/sqrt(length(patient_list)) nanstd(On_On_score_algo_NN_mean(:,idx)./Off_Off_score_algo_LDA_mean(:,idx))/sqrt(length(patient_list))],'ks-','LineWidth',2)

for k=1:length(patient_list)
plot(11:14,[Off_Off_score_algo_LR_mean(k,idx)/Off_Off_score_algo_LDA_mean(k,idx) On_Off_score_algo_LR_mean(k,idx)/Off_Off_score_algo_LDA_mean(k,idx) Off_On_score_algo_LR_mean(k,idx)/Off_Off_score_algo_LDA_mean(k,idx)  On_On_score_algo_LR_mean(k,idx)/Off_Off_score_algo_LDA_mean(k,idx)],'o-','color',cchoice{k})
end
errorbar(11:14,[nanmean(Off_Off_score_algo_LR_mean(:,idx)./Off_Off_score_algo_LDA_mean(:,idx)) nanmean(On_Off_score_algo_LR_mean(:,idx)./Off_Off_score_algo_LDA_mean(:,idx)) nanmean(Off_On_score_algo_LR_mean(:,idx)./Off_Off_score_algo_LDA_mean(:,idx)) nanmean(On_On_score_algo_LR_mean(:,idx)./Off_Off_score_algo_LDA_mean(:,idx))],[nanstd(Off_Off_score_algo_LR_mean(:,idx)./Off_Off_score_algo_LDA_mean(:,idx))/sqrt(length(patient_list)) nanstd(On_Off_score_algo_LR_mean(:,idx)./Off_Off_score_algo_LDA_mean(:,idx))/sqrt(length(patient_list)) nanstd(Off_On_score_algo_LR_mean(:,idx)./Off_Off_score_algo_LDA_mean(:,idx))/sqrt(length(patient_list)) nanstd(On_On_score_algo_LR_mean(:,idx)./Off_Off_score_algo_LDA_mean(:,idx))/sqrt(length(patient_list))],'ks-','LineWidth',2)

xlim([0 15])
plot([0 15],[1 1],'k-')

if z==1
    title('Individual All')
elseif z==2
    title('Permutation All')
elseif z==3 
    title('Individual Fp/F/C/P Contra')
elseif z==4
    title('Permutation Fp/F/C/P Contra')
end

end





%% Park 

% retrait de Park 7
patient_list = {'Park_1','Park_2','Park_3','Park_4','Park_5','Park_6','Park_8','Park_9','Park_10'};
side_movement={'gauche','gauche','gauche','droite','gauche','gauche','gauche','droite','gauche'};
minimal_number_individual=[170,195,198,183,220,192,200,171,179]; %for 2mvts
minimal_number_permutation=[28,32,33,30,36,32,33,28,29]; % for 2mvts
LDopa_Park=[460,210,1177,1849,1052,252,750,1106,157]; 
Score_UPDRS_moyen_Off=[0.66,1,3.33,1.33,3.66,1.33,2,3.33,1.33]; 
Score_UPDRS_moyen_On=[0.66,1,2.33,1,2,2,1.66,1.66,1]; 

for k=1:length(patient_list)
    
    for z=1:4
        
        if z==1
            saving_name='decoding_2mvts_individual_trials_all_decimate_window2_perpatient';
        elseif z==2
            saving_name='decoding_2mvts_permutation_trials_all_decimate_window2_perpatient';
            
        elseif z==3
            saving_name='decoding_2mvts_individual_trials_Fp_F_C_P_contra_decimate_window2_perpatient';
        elseif z==4
            saving_name='decoding_2mvts_permutation_trials_Fp_F_C_P_contra_decimate_window2_perpatient';
        end
 
        %rootpath=strcat('/Volumes/Charlotte/EEG_Analysis/DATA_EEG/PARK/',patient_list{k},'/Off/Decodage');
        rootpath=strcat('/Volumes/DBS_COLLAB/EEG_Analysis/DATA_EEG/PARK/',patient_list{k},'/Off/Decodage');
        cd(rootpath)
        load(strcat(rootpath,'/',saving_name,'.mat'))
        Off_score_algo_LDA_mean(k,z) = mean(score_algo_LDA(:));
        Off_score_algo_LDA_std(k,z) = std(score_algo_LDA(:)); %./sqrt(625);
        Off_score_algo_NN_mean(k,z) = mean(score_algo_NN(:));
        Off_score_algo_NN_std(k,z) = std(score_algo_NN(:)) ;%./sqrt(625);
        Off_score_algo_LR_mean(k,z) = mean(score_algo_LR(:));
        Off_score_algo_LR_std(k,z) = std(score_algo_LR(:)) ;%./sqrt(625);
        
        rootpath=strcat('/Volumes/DBS_COLLAB/EEG_Analysis/DATA_EEG/PARK/',patient_list{k},'/On/Decodage');
        cd(rootpath)
        load(strcat(rootpath,'/',saving_name,'.mat'))
        On_score_algo_LDA_mean(k,z) = mean(score_algo_LDA(:));
        On_score_algo_LDA_std(k,z) = std(score_algo_LDA(:)); %./sqrt(625);
        On_score_algo_NN_mean(k,z) = mean(score_algo_NN(:));
        On_score_algo_NN_std(k,z) = std(score_algo_NN(:)) ;%./sqrt(625);
        On_score_algo_LR_mean(k,z) = mean(score_algo_LR(:));
        On_score_algo_LR_std(k,z) = std(score_algo_LR(:)) ;%./sqrt(625);
        
    end
end

patient_list_temoins = {'Temoin_Laurent_Bis','Temoin_Jean','Temoin_Isabelle','Temoin_Francois','Temoin_6','Temoin_7','Temoin_8','Temoin_9','Temoin_10','Temoin_Sabrina'};
side_movement= {'droite','droite','droite','droite','droite','gauche','droite','droite','droite','droite'};

for k=1:length(patient_list_temoins)
    for z=1:4
    if z==1
        saving_name='decoding_2mvts_individual_trials_all_decimate_window2_perpatient';
    elseif z==2
        saving_name='decoding_2mvts_permutation_trials_all_decimate_window2_perpatient';
    elseif z==3
        saving_name='decoding_2mvts_individual_trials_Fp_F_C_P_contra_decimate_window2_perpatient';
    elseif z==4
        saving_name='decoding_2mvts_permutation_trials_Fp_F_C_P_contra_decimate_window2_perpatient';
    end
    
    rootpath=strcat('/Volumes/Charlotte/EEG_Analysis/DATA_EEG/TEMOINS/',patient_list_temoins{k},'/Decodage');
    cd(rootpath)
        load(strcat(rootpath,'/',saving_name,'.mat'))
    temoins_score_algo_LDA_mean(k,z) = mean(score_algo_LDA(:));
    temoins_score_algo_LDA_std(k,z) = std(score_algo_LDA(:)); %./sqrt(625);
    temoins_score_algo_NN_mean(k,z) = mean(score_algo_NN(:));
    temoins_score_algo_NN_std(k,z) = std(score_algo_NN(:)) ;%./sqrt(625);
    temoins_score_algo_LR_mean(k,z) = mean(score_algo_LR(:));
    temoins_score_algo_LR_std(k,z) = std(score_algo_LR(:)) ;%./sqrt(625);  
    end    
end

figure(); 
for k=1:4
    subplot(2,2,k); idx=k; hold on; 

plot(1:2,[Off_score_algo_LDA_mean(:,idx) On_score_algo_LDA_mean(:,idx)],'o-')
errorbar(1:2,[nanmean(Off_score_algo_LDA_mean(:,idx)) nanmean(On_score_algo_LDA_mean(:,idx))],[nanstd(Off_score_algo_LDA_mean(:,idx)) nanstd(On_score_algo_LDA_mean(:,idx))],'k-','LineWidth',2)
scatter(3*ones(10,1),temoins_score_algo_LDA_mean(:,idx),'o','MarkerEdgeColor',[0.8 0.8 0.8],"jitter","on")
errorbar(3,nanmean(temoins_score_algo_LDA_mean(:,idx)),[nanstd(temoins_score_algo_LDA_mean(:,idx))/sqrt(length(patient_list_temoins))],'ks','LineWidth',2)

plot(5:6,[Off_score_algo_NN_mean(:,idx) On_score_algo_NN_mean(:,idx)],'o-')
errorbar(5:6,[nanmean(Off_score_algo_NN_mean(:,idx)) nanmean(On_score_algo_NN_mean(:,idx))],[nanstd(Off_score_algo_NN_mean(:,idx)) nanstd(On_score_algo_NN_mean(:,idx))],'k-','LineWidth',2)
scatter(7*ones(10,1),temoins_score_algo_NN_mean(:,idx),'o','MarkerEdgeColor',[0.8 0.8 0.8],"jitter","on")
errorbar(7,nanmean(temoins_score_algo_NN_mean(:,idx)),[nanstd(temoins_score_algo_NN_mean(:,idx))/sqrt(length(patient_list_temoins))],'ks','LineWidth',2)

plot(9:10,[Off_score_algo_LR_mean(:,idx) On_score_algo_LR_mean(:,idx)],'o-')
errorbar(9:10,[nanmean(Off_score_algo_LR_mean(:,idx)) nanmean(On_score_algo_LR_mean(:,idx))],[nanstd(Off_score_algo_LR_mean(:,idx)) nanstd(On_score_algo_LR_mean(:,idx))],'k-','LineWidth',2)
scatter(11*ones(10,1),temoins_score_algo_LR_mean(:,idx),'o','MarkerEdgeColor',[0.8 0.8 0.8],"jitter","on")
errorbar(11,nanmean(temoins_score_algo_LR_mean(:,idx)),[nanstd(temoins_score_algo_LR_mean(:,idx))/sqrt(length(patient_list_temoins))],'ks','LineWidth',2)

set(gca,'Xtick',[2,5,8],'Xticklabel',{'LDA','NN','MLR'})
ylabel('Decoding accuracy')
plot([0 12],[50 50],'k.-')

if z==1
    title('Individual All')
elseif z==2
    title('Permutation All')
elseif z==3 
    title('Individual Fp/F/C/P Contra')
elseif z==4
    title('Permutation Fp/F/C/P Contra')
end

end


figure();
subplot(1,3,1); hold on; % quasi significatif pour LDA, pour NN(:,1) ou LR(:,1)/LR(:,2) p=0.026
scatter(LDopa_Park,Off_score_algo_LR_mean(:,2))
[r,p]=corrcoef(LDopa_Park,Off_score_algo_LR_mean(:,2)) % significatif pour les deux
mdl=fitlm(LDopa_Park,Off_score_algo_LR_mean(:,2))
plot(mdl)
ylabel('Decoder accuracy (MLR')
xlabel('Equivalence LDOPA treatment')
subplot(1,3,2); hold on; 
[r,p]=corrcoef(LDopa_Park,Off_score_algo_NN_mean(:,2)) % significatif pour les deux
scatter(LDopa_Park,Off_score_algo_NN_mean(:,2))
mdl=fitlm(LDopa_Park,Off_score_algo_NN_mean(:,2))
plot(mdl)
ylabel('Decoder accuracy (NN')
xlabel('Equivalence LDOPA treatment')
subplot(1,3,3); hold on; 
scatter(LDopa_Park,Off_score_algo_LDA_mean(:,2))
ylabel('Decoder accuracy (NN')
xlabel('Equivalence LDOPA treatment')
mdl=fitlm(LDopa_Park,Off_score_algo_LDA_mean(:,2))
plot(mdl)
[r,p]=corrcoef(LDopa_Park,Off_score_algo_LDA_mean(:,2)) % significatif pour les deux


figure();
subplot(1,3,1); hold on; % quasi significatif pour LDA, pour NN(:,1) ou LR(:,1)/LR(:,2) p=0.026
scatter(Score_UPDRS_moyen_Off,Off_score_algo_LR_mean(:,2))
scatter(Score_UPDRS_moyen_Off,Off_score_algo_LR_mean(:,4))
[r,p]=corrcoef(Score_UPDRS_moyen_Off,Off_score_algo_LR_mean(:,2));
ylabel('Decoder accuracy (MLR')
xlabel('Score UPDRS moyen Off')
subplot(1,3,2); hold on; 
scatter(Score_UPDRS_moyen_Off,Off_score_algo_NN_mean(:,2))
scatter(Score_UPDRS_moyen_Off,Off_score_algo_NN_mean(:,4))
[r,p]=corrcoef(Score_UPDRS_moyen_Off,Off_score_algo_NN_mean(:,2));
ylabel('Decoder accuracy (NN')
xlabel('Score UPDRS moyen Off')
subplot(1,3,3); hold on; 
scatter(Score_UPDRS_moyen_Off,Off_score_algo_LDA_mean(:,2))
scatter(Score_UPDRS_moyen_Off,Off_score_algo_LDA_mean(:,4))
[r,p]=corrcoef(LDopa_Park,Off_score_algo_LDA_mean(:,2));
ylabel('Decoder accuracy (NN')
xlabel('Score UPDRS moyen Off')


figure();
subplot(1,3,1); hold on; % quasi significatif pour LDA, pour NN(:,1) ou LR(:,1)/LR(:,2) p=0.026
scatter(Score_UPDRS_moyen_Off-Score_UPDRS_moyen_On,On_score_algo_LR_mean(:,2)-Off_score_algo_LR_mean(:,2))
scatter(Score_UPDRS_moyen_Off-Score_UPDRS_moyen_On,On_score_algo_LR_mean(:,4)-Off_score_algo_LR_mean(:,4))
[r,p]=corrcoef(Score_UPDRS_moyen_Off-Score_UPDRS_moyen_On,Off_score_algo_LR_mean(:,2))
ylabel('Decoder accuracy (MLR')
xlabel(' Delta Score UPDRS moyen Off/On')
subplot(1,3,2); hold on; 
scatter(Score_UPDRS_moyen_Off-Score_UPDRS_moyen_On,On_score_algo_NN_mean(:,1)-Off_score_algo_NN_mean(:,1))
scatter(Score_UPDRS_moyen_Off-Score_UPDRS_moyen_On,On_score_algo_NN_mean(:,2)-Off_score_algo_NN_mean(:,2))
[r,p]=corrcoef(Score_UPDRS_moyen_Off,Off_score_algo_NN_mean(:,2));
ylabel('Decoder accuracy (NN')
xlabel(' Delta Score UPDRS moyen Off/On')
subplot(1,3,3); hold on; 
scatter(Score_UPDRS_moyen_Off-Score_UPDRS_moyen_On,On_score_algo_LDA_mean(:,1)-Off_score_algo_LDA_mean(:,1))
scatter(Score_UPDRS_moyen_Off-Score_UPDRS_moyen_On,On_score_algo_LDA_mean(:,2)-Off_score_algo_LDA_mean(:,2))
[r,p]=corrcoef(LDopa_Park,Off_score_algo_LDA_mean(:,2));
ylabel('Decoder accuracy (NN')
xlabel(' Delta Score UPDRS moyen Off/On')
