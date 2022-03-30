import time
import argparse
import numpy as np
import random
import torch
import torch.nn.functional as F
import torch.optim as optim
import os
import glob
from torch.autograd import Variable
import sys
import torch.nn as nn
import pickle
import copy
from data_loading import *
from utils import *
from score import *
from inner_loop import InnerLoop
from mlp import mlp
from meta_learner_cv import *
import pandas as pd

# Training settings
# This model uses expresion + somatic mutations as features
# It applies the cross validation framework proposed by Siamese Network
parser = argparse.ArgumentParser()

parser.add_argument('--feature_dic', type=str, default='/cluster/projects/radiomics/Gyn_Autosegmentation/TCRP_data/Trametinib/', help='Feature folder')
parser.add_argument('--PDX_dic', type=str, default='/cluster/projects/radiomics/Gyn_Autosegmentation/TCRP_data/PDX/', help='Feature folder')
# parser.add_argument('--feature_dic', type=str, default='/cluster/home/nir/TCRP/data/Sorafenib/', help='Feature folder')
# parser.add_argument('--PDX_dic', type=str, default='/cluster/home/nir/TCRP/data/PDTC/', help='Feature folder')
parser.add_argument('--model_dic', type=str, default='/cluster/home/nir/TCRP/models/', help='Feature folder')
parser.add_argument('--drug', type=str, default='Trametinib', help='Treated drug')  # Cetuximab, Erlotinib, Paclitaxel, Tamoxifen, Trametinib
parser.add_argument('--seed', type=int, default=1, help='Random seed.')
parser.add_argument('--K', type=int, default=1, help='Perform K shot learning')
parser.add_argument('--meta_batch_size', type=int, default=10, help='Meta-learning batch size, i.e. how many different tasks we need sample')
parser.add_argument('--inner_batch_size', type=int, default=10, help='Batch size for each individual learning job')
parser.add_argument('--num_updates', type=int, default=20, help='Number of training epochs')
parser.add_argument('--num_inner_updates', type=int, default=1, help='Initial learning rate')
parser.add_argument('--num_out_updates', type=int, default=20, help='Final learning rate')
parser.add_argument('--num_trials', type=int, default=20, help='Number of trials for unseen tissue') # repeat 20 times = 20
parser.add_argument('--hidden', type=int, default=20, help='Number of hidden units of NN for single task') # ori = 60
parser.add_argument('--tissue_list', type=str, default='/cluster/projects/radiomics/Gyn_Autosegmentation/TCRP_data/Trametinib_tissue_map.pkl', help='Cell line list for different tissues, used for defining meta-tasks in the meta-learning phase')
parser.add_argument('--meta_lr', type=float, default=0.001, help='Learning rate for meta-learning update') # ori = 0.001
parser.add_argument('--inner_lr', type=float, default=0.001, help='Learning rate for ') # ori = 0.001
parser.add_argument('--tissue_num', type=int, default=12, help='Tissue number evolved in the inner update')
parser.add_argument('--layer', type=int, default=1, help='Number of layers of NN for single task') #ori = 1

args = parser.parse_args()

feature_dic = args.feature_dic
pdx_dic = args.PDX_dic
drug = args.drug
K = args.K
num_trials = args.num_trials
meta_batch_size = args.meta_batch_size
inner_batch_size = args.inner_batch_size
num_updates = args.num_updates
num_inner_updates = args.num_inner_updates

random.seed(args.seed)
np.random.seed(args.seed)
torch.manual_seed(args.seed)

tissue_list = args.tissue_list
print tissue_list
# Load tissue cell line mapping
with open(tissue_list, 'rb') as f:
	tissue_map = pickle.load(f)

model_dic = args.model_dic + 'MODELS/' + drug + '/'

if not os.path.exists(model_dic):
	mkdir_cmd = 'mkdir -p ' + model_dic
	os.system(mkdir_cmd) 
# Load data
#cv_feature_list, cv_label_list, meta_tissue_index_list, test_feature_list, test_label_list, test_tissue_list  = load_data_cell_line(tissue_map, drug, K)
train_feature, train_label, tissue_index_list = load_data(tissue_map, drug, K, path = feature_dic)
# train_feature, train_label, tissue_index_list, drug_test_feature, drug_test_label, _ = load_data_cell_line( tissue_map, drug, tissue, K, path=data_dic )

# PDTC_feature, PDTC_label = load_data_PDTC( args.drug, path = feature_dic )
PDX_feature, PDX_label = load_data_PDX(drug=drug, path=pdx_dic)
#best_param = (2,args.hidden,0.001,0.001,12)

layer, hidden, meta_lr, inner_lr, tissue_num = args.layer, args.hidden, args.meta_lr, args.inner_lr, args.tissue_num

meta_dataset = dataset(train_feature, train_label)
test_dataset = dataset(PDX_feature, PDX_label)

best_train_loss_test_corr_list, best_train_corr_test_corr_list, best_train_corr_test_scorr_list, best_train_scorr_test_scorr_list = [], [], [], []
# test_pred_pear = pd.DataFrame()
# test_pred_spear = pd.DataFrame()
PDX_pred = pd.DataFrame()
corr_save = pd.DataFrame()

best_zero_pear_corr_list = []
best_zero_spearman_corr_list = []
zero_save = pd.DataFrame()


for i in range(num_trials):
	meta_learner = MetaLearner( meta_dataset, test_dataset, K, meta_lr, inner_lr, layer, hidden, tissue_num, meta_batch_size, inner_batch_size, num_updates, num_inner_updates, tissue_index_list, num_trials )

	best_train_loss_test_corr, best_train_corr_test_corr, best_train_corr_test_scorr, best_train_scorr_test_scorr, best_model, best_test_pred_pear, best_test_pred_spear, test_true_label, out_pred, best_zero_pear_corr, best_zero_spearman_corr = meta_learner.train()
	best_train_loss_test_corr_list.append(best_train_loss_test_corr)
	best_train_corr_test_corr_list.append(best_train_corr_test_corr)
	best_train_corr_test_scorr_list.append(best_train_corr_test_scorr)
	best_train_scorr_test_scorr_list.append(best_train_scorr_test_scorr)

	best_zero_pear_corr_list.append(best_zero_pear_corr)
	best_zero_spearman_corr_list.append(best_zero_spearman_corr)
	# save best test prediction in each trial
	# test_pred_pear[i] = np.squeeze(best_test_pred_pear)
	# test_pred_spear[i] = np.squeeze(best_test_pred_spear)
	# test_pred_pear['true label'] = np.squeeze(test_true_label)
	# test_pred_spear['true label'] = np.squeeze(test_true_label)
	PDX_pred[i] = out_pred

    # Please uncomment this line to save your pre-train models
	# torch.save(best_model, model_dic + '/model_'+str(K)+'_trail_' + str(i))

a = np.asarray(best_train_loss_test_corr_list).mean()
b = np.asarray(best_train_corr_test_corr_list).mean() # pearson corr
c = np.asarray(best_train_corr_test_scorr_list).mean()
d = np.asarray(best_train_scorr_test_scorr_list).mean() #spearman corr
save_path_pred = '/cluster/home/nir/TCRP/prediction/' + 'PDX_pred_' + drug + '_k' + str(K) + '.csv'
# save_path_spear = '/cluster/home/nir/TCRP/prediction/' + 'spearman_' + drug + '_k' + str(K) + '.csv'
PDX_pred.to_csv(path_or_buf=save_path_pred, index=False)
# test_pred_spear.to_csv(path_or_buf=save_path_spear, index=False)
corr_save['Pearson'] = best_train_corr_test_corr_list
corr_save['Spearman'] = best_train_scorr_test_scorr_list
save_path_corr = '/cluster/home/nir/TCRP/prediction/' + 'PDX_corr_' + drug + '_k' + str(K) + '.csv'
# corr_save.to_csv(path_or_buf=save_path_corr, index=False)

zero_save['Pearson'] = best_zero_pear_corr_list
zero_save['Spearman'] = best_zero_spearman_corr_list
save_path_zero = '/cluster/home/nir/TCRP/prediction/' + 'PDX_zero_corr_' + drug + '.csv'
zero_save.to_csv(path_or_buf=save_path_zero, index=False)

print 'PDX best_train_loss_test_corr:', float('%.3f'%a), 'best_train_corr_test_corr', float('%.3f'%b), 'best_train_corr_test_scorr', float('%.3f'%c), 'best_train_scorr_test_scorr', float('%.3f'%d)
# print 'best_train_corr_test_corr_list - pearson corr'
# print(best_train_corr_test_corr_list)
# print 'best_train_scorr_test_scorr_list - spearman corr'
# print(best_train_scorr_test_scorr_list)
# # print 'test true label'
# print(test_true_label)

# use the last trail of model for prediction