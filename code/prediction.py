import time
import argparse
import numpy as np
import random
from scipy.fftpack import next_fast_len
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
from mlp import mlp	


parser = argparse.ArgumentParser()

parser.add_argument('--feature_dic', type=str, default='/cluster/projects/radiomics/Gyn_Autosegmentation/TCRP_data/Cetuximab/', help='Feature folder')
parser.add_argument('--PDX_dic', type=str, default='/cluster/projects/radiomics/Gyn_Autosegmentation/TCRP_data/PDX/', help='Feature folder')
parser.add_argument('--model_dic', type=str, default='/cluster/home/nir/TCRP/models/', help='Feature folder')
parser.add_argument('--drug', type=str, default='Cetuximab', help='Treated drug')
parser.add_argument('--seed', type=int, default=19, help='Random seed.')
parser.add_argument('--K', type=int, default=10, help='Perform K shot learning')
parser.add_argument('--meta_batch_size', type=int, default=32, help='Meta-learning batch size, i.e. how many different tasks we need sample')
parser.add_argument('--inner_batch_size', type=int, default=10, help='Batch size for each individual learning job')
parser.add_argument('--num_updates', type=int, default=20, help='Number of training epochs')
parser.add_argument('--num_inner_updates', type=int, default=1, help='Initial learning rate')
parser.add_argument('--num_out_updates', type=int, default=20, help='Final learning rate')
parser.add_argument('--num_trials', type=int, default=20, help='Number of trials for unseen tissue') # repeat 20 times = 20
parser.add_argument('--hidden', type=int, default=10, help='Number of hidden units of NN for single task') # ori = 60
parser.add_argument('--tissue_list', type=str, default='/cluster/projects/radiomics/Gyn_Autosegmentation/TCRP_data/Cetuximab_tissue_map.pkl', help='Cell line list for different tissues, used for defining meta-tasks in the meta-learning phase')
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

# PDTC_feature, PDTC_label = load_data_PDTC( args.drug, path = feature_dic )
PDX_feature, PDX_label = load_data_PDX(drug=drug, path=pdx_dic)
#best_param = (2,args.hidden,0.001,0.001,12)

layer, hidden, meta_lr, inner_lr, tissue_num = args.layer, args.hidden, args.meta_lr, args.inner_lr, args.tissue_num

meta_dataset = dataset(train_feature, train_label)
test_dataset = dataset(PDX_feature, PDX_label)

def get_unseen_data_loader(feature, label, K, batch_size=1):

	index_list = np.random.permutation(feature.shape[0])

	train_index_list = index_list[0:K]
	test_index_list = index_list[K:]

	train_feature = torch.FloatTensor( feature[train_index_list,:] )
	train_label = torch.FloatTensor( label[train_index_list,] )

	test_feature = torch.FloatTensor( feature[test_index_list,:] )
	test_label = torch.FloatTensor( label[test_index_list,] )

	train_dataset = du.TensorDataset( train_feature, train_label )
	test_dataset = du.TensorDataset( test_feature, test_label )

	train_loader = du.DataLoader(train_dataset, batch_size=batch_size)
	train_data_list = []
	for batch_feature, batch_label in train_loader:
		train_data_list.append((batch_feature.cuda(), batch_label.cuda()))
	
	test_loader = du.DataLoader(test_dataset, batch_size=batch_size)
	test_data_list = []
	for batch_feature, batch_label in test_loader:
		test_data_list.append((batch_feature.cuda(), batch_label.cuda()))

	return train_data_list, test_data_list

unseen_train_data, unseen_test_data = get_unseen_data_loader(test_dataset.feature, test_dataset.label, K)

net_path = '/cluster/home/nir/TCRP/models/MODELS/Cetuximab/model_10_trail_19'
net = mlp(503, 1, 10)
net = torch.load(net_path)

print(net)
net.eval()

for i, (in_, target) in enumerate(unseen_test_data):

	input_var = Variable(in_)
	target_var = Variable(target)

	# Second output is hidden
	out, hidden = net(input_var)
	print(str(out), str(hidden), str(target_var))
