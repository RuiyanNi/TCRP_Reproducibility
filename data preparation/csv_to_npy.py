from cProfile import label
from turtle import done
from matplotlib.pyplot import axis
import numpy as np
import pandas as pd
import pickle
import glob
import shutil
from os import listdir
from os.path import isfile, join
import numpy as np
import random
import os

data_path = r'D:\UofT\Yr1 Sem2\MBP AI\data\step1'
feature_path = r'D:\UofT\Yr1 Sem2\MBP AI\data\step1\feature_by_tissue'
# label_path = r'D:\UofT\Yr1 Sem2\MBP AI\data\step1\label_by_tissue'
drug = 'Trametinib' # Cetuximab, Erlotinib, Paclitaxel, Tamoxifen, Trametinib
pdx_file = pd.read_csv(r"D:\UofT\Yr1 Sem2\MBP AI\data\step2\new_data\pdxCetuximab_rmNA.csv")
# pdx_file = pdx_file.drop(['Unnamed: 0'], axis=1)
list_pdx = pdx_file.columns

feature_file = [f for f in listdir(feature_path) if isfile(join(feature_path, f))]
# label_file = [f for f in listdir(label_path) if isfile(join(label_path, f))]

def save_feature_description(drug):
    '''convert feature, label, and description - for pretrained data'''
    for files in feature_file:
        '''save feature'''
        if drug+'_feature' in files:
            df = pd.read_csv(os.path.join(feature_path, files))
            # df1 = df.T
            # df = df.drop(['Unnamed: 0'], axis=1)
            list_df = df.columns
            list_df_pdx = [value for value in list_df if value in list_pdx]
            df = df[list_df_pdx]
            for column in df.columns:
                if 'mutation' in column:
                    df[column] = np.nan_to_num(df[column])
                if 'expression' in column:
                    df = df[df[column].notna()]

            cell_line = df['Unnamed: 0'].tolist()
            df = df.drop(['Unnamed: 0'], axis=1)
            df.index, df.columns = [range(df.index.size), range(df.columns.size)]
            save_name = os.path.join(data_path, drug)
            save_name = os.path.join(save_name, files[:-4]+'.npy')
            np.save(save_name, df)
            print(files, str(df.shape))
            
            '''save label'''
            read_label_path = files[:-11] + 'label.csv'
            read_label_path2 = os.path.join(label_path, read_label_path)
            label = pd.read_csv(read_label_path2)
            df_new = label[label['Unnamed: 0'].isin(cell_line)]
            label_save = np.array(df_new['AUC'])
            save_name = os.path.join(data_path, drug)
            save_name = os.path.join(save_name, read_label_path[:-4]+'.npy')
            np.save(save_name, label_save)
            print(read_label_path, str(label_save.shape))

            '''save description'''
            save_name = os.path.join(data_path, drug)
            save_name = os.path.join(save_name, files[:-4]+'_description.npy')
            np.save(save_name, list_df_pdx)
            print(files)
    
# save_feature_description(drug=drug)

def save_feature_description2(drug):
    '''convert feature and description - for PDX'''
    for files in feature_file:
        '''save feature'''
        if drug+'_feature' in files:
            df = pd.read_csv(os.path.join(feature_path, files))
            # df = df.drop(['Unnamed: 0'], axis=1)
            '''save label'''
            label = np.array(df['AUC'])
            save_name = os.path.join(data_path, drug)
            save_name = os.path.join(save_name, files[:-11]+'label.npy')
            np.save(save_name, label)

            list_df = df.columns
            list_df_pdx = [value for value in list_df if value in list_pdx]
            df = df[list_df_pdx]
            '''save feature name'''
            feature_name = df.columns
            feature_name = feature_name[1:]
            feature_name = np.array(feature_name)
            save_name = os.path.join(data_path, drug)
            save_name = os.path.join(save_name, files[:-11]+ 'feature_description.npy')
            np.save(save_name, feature_name)
            '''save feature'''
            df = df.drop(['Unnamed: 0'], axis=1)
            df.index, df.columns = [range(df.index.size), range(df.columns.size)]
            save_name = os.path.join(data_path, drug)
            save_name = os.path.join(save_name, files[:-4]+'.npy')
            np.save(save_name, df)
            print(files, str(df.shape))

# save_feature_description2(drug)

def save_label(drug):
    '''convert label'''
    for files in label_file:
        if drug in files:
            df = pd.read_csv(os.path.join(label_path, files))
            df = df.drop(['Unnamed: 0'], axis=1)
            df.columns = range(df.columns.size)
            save_name = os.path.join(data_path, drug)
            save_name = os.path.join(save_name, files[:-4]+'.npy')
            np.save(save_name, df)
            print(files)

# save_label(drug=drug)

def save_feature_description_pdx(drug, path = r'D:\UofT\Yr1 Sem2\MBP AI\data\step2\new_data'):
    '''convert feature and description - for PDX'''
    file_path = path + '\pdx' + drug + '_rmNA.csv'
    df = pd.read_csv(file_path)
    # df = df.drop(['Unnamed: 0'], axis=1)
    list_df = df.columns
    list_pretrain = np.load(r"D:\UofT\Yr1 Sem2\MBP AI\data\pretrained_data\Cetuximab\aero_dig_tract_Cetuximab_feature_description.npy", allow_pickle=True).tolist()
    list_pretrain_pdx = [value for value in list_df if value in list_pretrain]
    df = df[list_pretrain_pdx]
    '''save feature name'''
    feature_name = df.columns
    feature_name = feature_name[1:]
    feature_name = np.array(feature_name)
    save_name = file_path[:-9] + '_feature_description.npy'
    np.save(save_name, feature_name)
    '''save label'''
    label = np.array(df['Unnamed: 0'])
    save_name = file_path[:-9] + '_label.npy'
    np.save(save_name, label)
    '''save feature'''
    df = df.drop(['Unnamed: 0'], axis=1)
    df.index, df.columns = [range(df.index.size), range(df.columns.size)]
    save_name = file_path[:-9] + '_feature.npy'
    np.save(save_name, df)
    print(df.shape)

save_feature_description_pdx(drug)

def print_nan_row_name(drug):
    for files in feature_file:
        if drug+'_feature' in files:
            df = pd.read_csv(os.path.join(feature_path, files))
            # df1 = df.T
            # df = df.drop(['Unnamed: 0'], axis=1)
            list_df = df.columns
            list_df_pdx = [value for value in list_df if value in list_pdx]
            df = df[list_df_pdx]
            for column in df.columns:
                if 'expression' in column:
                    nan_df = df[df[column].isna()]
                    nan_row_name = str(nan_df['Unnamed: 0'])
            print(nan_row_name)

# print_nan_row_name(drug=drug)