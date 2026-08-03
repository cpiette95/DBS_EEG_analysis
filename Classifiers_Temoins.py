#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 27 00:31:10 2021

@author: charlotte.piette
"""
import tensorflow.keras as keras

from tensorflow.keras.wrappers.scikit_learn import KerasClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_score

from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neighbors import NearestCentroid
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn import svm
from sklearn import metrics


import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.io
from scipy import signal
import scipy.signal as ss


## TEMOINS
patient_list = ['Temoin_Laurent_Bis','Temoin_Jean','Temoin_Isabelle','Temoin_Francois','Temoin_6','Temoin_7','Temoin_8','Temoin_9','Temoin_10','Temoin_Sabrina']
side_movement= ['droite','droite','droite','droite','droite','gauche','droite','droite','droite','droite']
#minimal_individual_per_patient=[224,169,180,140,199,180,112,194,174,198]
#minimal_permutation_per_patient=[37,28,30,23,33,30,18,32,29,33] 
minimal_individual_per_patient=[224,220,190,140,215,213,179,226,184,203] # 2mvts
minimal_permutation_per_patient=[37,36,31,23,35,35,29,37,30,33] # 2mvts


condition_name='Temoin'


average_cortical = 1; 
number_permutations = 25;     
label={'Fp1','Fp2','F3','F4','Fz','T3','T4','C3','C4','Cz','P3','P4'};

number_seeds = 5; 
number_splits = 5;
max_i = 300; 
seed_value=[57,67,45,12,89];



for Z in range(len(patient_list)): 
    
    patient_number = patient_list[Z] 
    
    if side_movement[Z]=='gauche':
        electrode_name_1 =4; 
        #electrode_name_2 = 6; 
        #electrode_name_3 = 11; 
        #electrodes=[4,6,11,12];
        electrodes=[4,6,11,14];
        
        #electrodes=[3,4,5,6,7,10,11,12,13,14]; 
        
    else: 
        electrode_name_1=3;
        #electrode_name_2 =5; 
        #electrode_name_3 =10;
        #electrodes=[3,5,10,12]; classique
        electrodes=[3,5,10,13]; # ontra
        
        #electrodes=[3,4,5,6,7,10,11,12,13,14]; #
    
    
    
    if average_cortical == 0: 
        
            os.chdir('/Volumes/DBS/DBS/EEG_Analysis/DATA_EEG/TEMOINS/'+patient_number) #+'/'+condition)  
            extension_name='EMG_based_full_cortical_ERP_filtered'
            name_data='full_cortical_ERP_filtered_artefact_free'
            saving_name = 'decoding_2mvts_individual_trials_C_contra_decimate18_window2_perpatient_training_scores'
        
        
            mat_1 = scipy.io.loadmat(condition_name+'_index_'+extension_name+'.mat')  
            EEG_simple_1 = mat_1[name_data];
            
            mat_2 = scipy.io.loadmat(condition_name+'_poing_'+extension_name+'.mat')  
            EEG_simple_2 = mat_2[name_data];
            
            #mat_3 = scipy.io.loadmat(condition_name+'_marionnette_'+extension_name+'.mat')   
            #EEG_simple_3 = mat_3[name_data];
            
            minimal_number=minimal_individual_per_patient[Z];
         
             
            decimate_factor = 18; temporal_window = np.arange(304,1024);
            
            EEG_simple_1_bis= np.zeros((minimal_number,len(electrodes),int(len(temporal_window)/decimate_factor)))
            for k in range(minimal_number):
                u=0;
                for n in electrodes:
                    EEG_simple_1_bis[k,u,:] = ss.decimate(EEG_simple_1[0,n][k,temporal_window],decimate_factor);
                    u=u+1; 
            EEG_simple_1_bis=EEG_simple_1_bis.reshape((minimal_number,len(electrodes)*EEG_simple_1_bis.shape[2]))    
            
            EEG_simple_2_bis= np.zeros((minimal_number,len(electrodes),int(len(temporal_window)/decimate_factor)))
            for k in range(minimal_number):
                u=0;
                for n in electrodes:
                    EEG_simple_2_bis[k,u,:] = ss.decimate(EEG_simple_2[0,n][k,temporal_window],decimate_factor);
                    u=u+1; 
            EEG_simple_2_bis=EEG_simple_2_bis.reshape((minimal_number,len(electrodes)*EEG_simple_2_bis.shape[2]))    
            
            #EEG_simple_3_bis= np.zeros((minimal_number,len(electrodes),int(len(temporal_window)/decimate_factor)))
            #for k in range(minimal_number):
            #    u=0;
            #    for n in electrodes:
            #        EEG_simple_3_bis[k,u,:] = ss.decimate(EEG_simple_3[0,n][k,temporal_window],decimate_factor);
            #        u=u+1; 
            #EEG_simple_3_bis=EEG_simple_3_bis.reshape((minimal_number,len(electrodes)*EEG_simple_3_bis.shape[2]))                 
                
            Full_data = np.concatenate((EEG_simple_1_bis,EEG_simple_2_bis)); #,EEG_simple_3_bis));
            target = np.concatenate((np.ones((minimal_number,1)),2*np.ones((minimal_number,1)) )); #,3*np.ones((minimal_number,1)) ));

        
        
            train_score_algo_LR = np.zeros((number_seeds*number_splits,1))
            train_score_algo_LDA = np.zeros((number_seeds*number_splits,1))
            train_score_algo_NN = np.zeros((number_seeds*number_splits,1))
            
            score_algo_LR = np.zeros((number_seeds*number_splits,1))
            score_algo_LDA = np.zeros((number_seeds*number_splits,1))
            score_algo_NN = np.zeros((number_seeds*number_splits,1))
            
            
            f=0
            
            number_trials = minimal_number ;
                                            
            for n in range(number_seeds):
            
                seed = seed_value[n]
                
                print(seed)
                np.random.seed(seed)
                
            
                 
                kfold = StratifiedKFold(n_splits=number_splits, shuffle=True, random_state=seed)
                                   
                for train_index, test_index in kfold.split(Full_data,target):
                 
            
                    x_train, x_test = Full_data[train_index,:],Full_data[test_index,:]
                    y_train,y_test = target[train_index],target[test_index]
                           
                    mul_lr = LogisticRegression(multi_class='multinomial', solver='newton-cg',max_iter=max_i)
                    mul_lr.fit(x_train, y_train.ravel())
                    train_score_LR[f] = mul_lr.score(x_train, y_train) * 100  
                    score_algo_LR[f] = mul_lr.score(x_test, y_test)*100
                    print(mul_lr.score(x_test,y_test)*100)
                    #predictions = mul_lr.predict(x_test)
                    #cm_LR = metrics.confusion_matrix(y_test,predictions)        
                    
                    try: 
                        lda = LinearDiscriminantAnalysis(solver='svd')
                        lda.fit(x_train,y_train.ravel())
                        train_score_algo_LDA[f]=lda.score(x_train,y_train)*100
                        score_algo_LDA[f]=lda.score(x_test,y_test)*100
                        print(lda.score(x_test,y_test)*100)
                    except: 
                        pass
                    
                    
                    clf = NearestCentroid(metric='euclidean',shrink_threshold=None)
                    clf.fit(x_train,y_train.ravel())
                    train_score_algo_NN[f]=clf.score(x_train,y_train)*100
                    score_algo_NN[f]=clf.score(x_test,y_test)*100
                    print(clf.score(x_test,y_test)*100)
                   
            
                    f=f+1; 
             
                    
         
            os.chdir('/Volumes/DBS/DBS/EEG_Analysis/DATA_EEG/TEMOINS/'+patient_number+'/Decodage') 
            scipy.io.savemat(saving_name+'.mat',{'score_algo_LDA':score_algo_LDA,'score_algo_LR':score_algo_LR,'score_algo_NN':score_algo_NN})
            scipy.io.savemat(saving_name+'.mat',{'score_algo_LDA':train_score_algo_LDA,'score_algo_LR':train_score_algo_LR,'score_algo_NN':train_score_algo_NN})
    
    
    
    
    
    else: 
        

        os.chdir('/Volumes/DBS/DBS/EEG_Analysis/DATA_EEG/TEMOINS/'+patient_number+'/Average_Permutations')  
        extension_name='average_cortical_EMG_based_ERP_filtered_50Hz'
        name_data='average_cortical_ERP_filtered'
        saving_name = 'decoding_2mvts_permutation_trials_Fp_F_C_P_contra_decimate72_window2_perpatient' 
        saving_name_bis = 'decoding_2mvts_permutation_trials_Fp_F_C_P_contra_decimate72_window2_perpatient_training_scores' 
    
    
        train_score_algo_LR = np.zeros((number_seeds*number_splits,number_permutations))
        train_score_algo_LDA = np.zeros((number_seeds*number_splits,number_permutations))
        train_score_algo_NN = np.zeros((number_seeds*number_splits,number_permutations))

        score_algo_LR = np.zeros((number_seeds*number_splits,number_permutations))
        score_algo_LDA = np.zeros((number_seeds*number_splits,number_permutations))
        score_algo_NN = np.zeros((number_seeds*number_splits,number_permutations))
    
        for v in range(number_permutations):
        
                 
            mat_1 = scipy.io.loadmat(condition_name+'_index_'+extension_name+'_'+str(v+1)+'.mat')  
            EEG_simple_1 = mat_1[name_data];
            
            mat_2 = scipy.io.loadmat(condition_name+'_poing_'+extension_name+'_'+str(v+1)+'.mat')  
            EEG_simple_2 = mat_2[name_data];
            
            #mat_3 = scipy.io.loadmat(condition_name+'_marionnette_'+extension_name+'_'+str(v+1)+'.mat')   
            #EEG_simple_3 = mat_3[name_data];
            
            minimal_number=minimal_permutation_per_patient[Z]; 
                
                 
            decimate_factor = 72; temporal_window = np.arange(304,1024);
            
            EEG_simple_1_bis= np.zeros((minimal_number,len(electrodes),int(len(temporal_window)/decimate_factor)))
            for k in range(minimal_number):
                u=0;
                for n in electrodes:
                    EEG_simple_1_bis[k,u,:] = ss.decimate(EEG_simple_1[0,n][k,temporal_window],decimate_factor);
                    u=u+1; 
            EEG_simple_1_bis=EEG_simple_1_bis.reshape((minimal_number,len(electrodes)*EEG_simple_1_bis.shape[2]))    
            
            EEG_simple_2_bis= np.zeros((minimal_number,len(electrodes),int(len(temporal_window)/decimate_factor)))
            for k in range(minimal_number):
                u=0;
                for n in electrodes:
                    EEG_simple_2_bis[k,u,:] = ss.decimate(EEG_simple_2[0,n][k,temporal_window],decimate_factor);
                    u=u+1; 
            EEG_simple_2_bis=EEG_simple_2_bis.reshape((minimal_number,len(electrodes)*EEG_simple_2_bis.shape[2]))    
            
            # EEG_simple_3_bis= np.zeros((minimal_number,len(electrodes),int(len(temporal_window)/decimate_factor)))
            # for k in range(minimal_number):
            #     u=0;
            #     for n in electrodes:
            #         EEG_simple_3_bis[k,u,:] = ss.decimate(EEG_simple_3[0,n][k,temporal_window],decimate_factor);
            #         u=u+1; 
            # EEG_simple_3_bis=EEG_simple_3_bis.reshape((minimal_number,len(electrodes)*EEG_simple_3_bis.shape[2]))                 
                
            Full_data = np.concatenate((EEG_simple_1_bis,EEG_simple_2_bis)); #,EEG_simple_3_bis));
            target = np.concatenate((np.ones((minimal_number,1)),2*np.ones((minimal_number,1)) )); #,3*np.ones((minimal_number,1)) ));


                
                 
            number_trials = minimal_number ;
            f=0;
                                                 
            for p in range(number_seeds):
                
                    seed = seed_value[p]    
                    print(seed)
                    np.random.seed(seed)
                
                     
                    kfold = StratifiedKFold(n_splits=number_splits, shuffle=True, random_state=seed)
                                       
                    for train_index, test_index in kfold.split(Full_data,target):
                     
                
                        x_train, x_test = Full_data[train_index,:],Full_data[test_index,:]
                        y_train,y_test = target[train_index],target[test_index]
                               
                        mul_lr = LogisticRegression(multi_class='multinomial', solver='newton-cg',max_iter=max_i)
                        mul_lr.fit(x_train, y_train.ravel())
                        train_score_algo_LR[f,v] = mul_lr.score(x_train, y_train)*100
                        score_algo_LR[f,v] = mul_lr.score(x_test, y_test)*100
                        print(mul_lr.score(x_test,y_test)*100)
                        #predictions = mul_lr.predict(x_test)
                        #cm_LR = metrics.confusion_matrix(y_test,predictions)        
                        
                        lda = LinearDiscriminantAnalysis(solver='svd')
                        lda.fit(x_train,y_train.ravel())
                        train_score_algo_LDA[f,v]=lda.score(x_train,y_train)*100
                        score_algo_LDA[f,v]=lda.score(x_test,y_test)*100
                        print(lda.score(x_test,y_test)*100)
                        
                        clf = NearestCentroid(metric='euclidean',shrink_threshold=None)
                        clf.fit(x_train,y_train.ravel())
                        train_score_algo_NN[f,v]=clf.score(x_train,y_train)*100
                        score_algo_NN[f,v]=clf.score(x_test,y_test)*100
                        print(clf.score(x_test,y_test)*100)
                       
            
                        f=f+1; 
                
         
      
        os.chdir('/Volumes/DBS/DBS/EEG_Analysis/DATA_EEG/TEMOINS/'+patient_number+'/Decodage')  
        scipy.io.savemat(saving_name+'.mat',{'score_algo_LDA':score_algo_LDA,'score_algo_LR':score_algo_LR,'score_algo_NN':score_algo_NN})
        scipy.io.savemat(saving_name_bis+'.mat',{'score_algo_LDA':train_score_algo_LDA,'score_algo_LR':train_score_algo_LR,'score_algo_NN':train_score_algo_NN})
    