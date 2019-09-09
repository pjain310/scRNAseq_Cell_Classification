import os
import numpy as np
import pandas as pd
import time as tm
from sklearn.svm import LinearSVC
from sklearn.preprocessing import RobustScaler
import rpy2.robjects as robjects


def run_SVM(DataPath, LabelsPath, CV_RDataPath, OutputDir, GeneOrderPath = "", NumGenes = 0):
    '''
    run baseline classifier: SVM
    Wrapper script to run an SVM classifier with a linear kernel on a benchmark dataset with 5-fold cross validation,
    outputs lists of true and predicted cell labels as csv files, as well as computation time.

    Parameters
    ----------
    DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes
    as row names and gene names as column names.
    LabelsPath : Cell population annotations file path (.csv).
    CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
    OutputDir : Output directory defining the path of the exported file.
    GeneOrderPath : Gene order file path (.csv) obtained from feature selection,
    defining the genes order for each cross validation fold, default is NULL.
    NumGenes : Number of genes used in case of feature selection (integer), default is 0.
    '''

    # read the Rdata file
    robjects.r['load'](CV_RDataPath)

    nfolds = np.array(robjects.r['n_folds'], dtype = 'int')
    tokeep = np.array(robjects.r['Cells_to_Keep'], dtype = 'bool')
    col = np.array(robjects.r['col_Index'], dtype = 'int')
    col = col - 1
    test_ind = np.array(robjects.r['Test_Idx'])
    train_ind = np.array(robjects.r['Train_Idx'])

    # read the data
    data = pd.read_csv(DataPath,index_col=0,sep=',')
    labels = pd.read_csv(LabelsPath, header=0,index_col=None, sep=',', usecols = col)

    labels = labels.iloc[tokeep]
    data = data.iloc[tokeep]

    # read the feature file
    if (NumGenes > 0):
        features = pd.read_csv(GeneOrderPath,header=0,index_col=None, sep=',')

    # folder with results
    os.chdir(OutputDir)

    # normalize data
    data = np.log1p(data)

    #Use standardscaler to fit data
    transformer = RobustScaler(quantile_range=(10,90)).fit(data)
    data = pd.DataFrame(transformer.transform(data))

    Classifier = LinearSVC()

    tr_time=[]
    ts_time=[]
    truelab = []
    pred = []

    for i in range(np.squeeze(nfolds)):
        test_ind_i = np.array(test_ind[i], dtype = 'int') - 1
        train_ind_i = np.array(train_ind[i], dtype = 'int') - 1

        train=data.iloc[train_ind_i]
        test=data.iloc[test_ind_i]
        y_train=labels.iloc[train_ind_i]
        y_test=labels.iloc[test_ind_i]

        if (NumGenes > 0):
            feat_to_use = features.iloc[0:NumGenes,i]
            train = train.iloc[:,feat_to_use]
            test = test.iloc[:,feat_to_use]

        start=tm.time()
        Classifier.fit(train, y_train.values.ravel())
        tr_time.append(tm.time()-start)

        start=tm.time()
        predicted = Classifier.predict(test)
        print(predicted)
        ts_time.append(tm.time()-start)

        truelab.extend(y_test.values)
        pred.extend(predicted)

    truelab = pd.DataFrame(truelab)
    pred = pd.DataFrame(pred)

    tr_time = pd.DataFrame(tr_time)
    ts_time = pd.DataFrame(ts_time)

    if (NumGenes == 0):
        truelab.to_csv("SVM_True_Labels.csv", index = False)
        pred.to_csv("SVM_Pred_Labels.csv", index = False)
        tr_time.to_csv("SVM_Training_Time.csv", index = False)
        ts_time.to_csv("SVM_Testing_Time.csv", index = False)
    else:
        truelab.to_csv("SVM_" + str(NumGenes) + "_True_Labels.csv", index = False)
        pred.to_csv("SVM_" + str(NumGenes) + "_Pred_Labels.csv", index = False)
        tr_time.to_csv("SVM_" + str(NumGenes) + "_Training_Time.csv", index = False)
        ts_time.to_csv("SVM_" + str(NumGenes) + "_Testing_Time.csv", index = False)

LabelsPath = "~/Desktop/scRNA_Cell_Typing/scRNAseq_Benchmark_datasets/Zheng/Labels.csv"
DataPath = "~/Desktop/scRNA_Cell_Typing/scRNAseq_Benchmark_datasets/Zheng/Filtered_68K_PBMC_data.csv"
CV_RDataPath = "~/Desktop/scRNA_Cell_Typing/CV_folds.RData"
GeneOrderPath = "~/Desktop/scRNA_Cell_Typing/results/og_SVM_results/rank_genes_dropouts.csv"
OutputDir = "results/top5000_SVM_scaled_results"

run_SVM(DataPath, LabelsPath, CV_RDataPath, OutputDir, GeneOrderPath, 5000)
