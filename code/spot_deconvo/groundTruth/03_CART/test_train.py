import os
os.chdir('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')

import pandas as pd
import numpy as np

from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import classification_report
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn import tree, svm

import pyhere
from pathlib import Path
import pickle

df_path = pyhere.here("processed-data", "spot_deconvo", "groundTruth", "02_samui_manual_annotation", "annotations_SCP_processed.csv")
df_path_out = pyhere.here('processed-data', 'spot_deconvo', 'groundTruth', '03_CART', 'annotation_dataset.pkl')

num_cell_types = 5
test_proportion = 0.2 # for training/test split
random_seed = 0

################################################################################
#   Analysis
################################################################################

#-------------------------------------------------------------------------------
#   Preprocess and gather fluorescence data + original manual cell-type labels
#-------------------------------------------------------------------------------

df = pd.read_csv(df_path,header=0)

#   Define the inputs (features we want the model to access) and outputs to the
#   model
x = df.loc[:, ['NeuN', 'TMEM119', 'GFAP', 'OLIG2', 'area', 'label_sample']]
y = df['label']

#   Split data into training and test sets (80%: 20%), evenly stratified across
#   classes and sample ID
x_train, x_test, y_train, y_test = train_test_split(
    x, y, test_size = test_proportion, random_state = random_seed,
    stratify = x['label_sample']
)


#   Remove the column we only used to properly stratify the data
x_train.drop('label_sample', axis = 1, inplace = True)
x_test.drop('label_sample', axis = 1, inplace = True)


#-------------------------------------------------------------------------------
#   Combine the two different sets of manual annotation into one dataset and
#   write to disk
#-------------------------------------------------------------------------------

perc_train = round(100 * x_train.shape[0] / (x_train.shape[0] + x_test.shape[0]), 2)
print(f'Using {x_train.shape[0]} training and {x_test.shape[0]} test examples ({perc_train}% training).')

#   Write the dataset to disk for later use
with open(df_path_out, 'wb') as f:
    pickle.dump((x_train, x_test, y_train, y_test), f)

#-------------------------------------------------------------------------------
#   Try a decision tree classifier
#-------------------------------------------------------------------------------
tuned_parameters = [
    {
        'criterion': ['gini', 'entropy'],
        'max_depth': [2, 3, 4, 5],
        'min_samples_leaf': [1, 5, 10, 15, 20, 25],
        'ccp_alpha': [0, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3]
    }
]

grid = GridSearchCV(
    tree.DecisionTreeClassifier(
        random_state = random_seed, splitter = 'best', class_weight = None
    ),
    tuned_parameters, cv=5, scoring = 'accuracy'
)
grid.fit(x_train, y_train)

#   Compute training and test accuracy on the best model
acc_train = round(100 * grid.best_estimator_.score(x_train, y_train), 1)
acc_test = round(100 * grid.best_estimator_.score(x_test, y_test), 1)
print(f'CART training accuracy: {acc_train}%.')
print(f'CART test accuracy: {acc_test}%.')
print(f'Best params: {grid.best_params_}')

#   Print a more thorough report about training and test scores, making sure
#   performance is good across all classes (since we optimized for accuracy)
labels_train = grid.best_estimator_.predict(x_train)
labels_test = grid.best_estimator_.predict(x_test)
print('Training report:\n', classification_report(y_train, labels_train))
print('Test report:\n', classification_report(y_test, labels_test))

#-------------------------------------------------------------------------------
#   Try logistic regression
#-------------------------------------------------------------------------------

tuned_parameters = [
    {
        'logisticregression__C': np.logspace(-2, 4, 20)
    }
]

pipe = make_pipeline(
    StandardScaler(),
    LogisticRegression(random_state = random_seed, max_iter=10000)
)

grid = GridSearchCV(pipe, tuned_parameters, cv=5, scoring = 'accuracy')
grid.fit(x_train, y_train)

acc_train = round(100 * grid.best_estimator_.score(x_train, y_train), 1)
acc_test = round(100 * grid.best_estimator_.score(x_test, y_test), 1)
print(f'Logistic regression training accuracy: {acc_train}%.')
print(f'Logistic regression test accuracy: {acc_test}%.')
print(f'Best params: {grid.best_params_}')

#   Print a more thorough report about training and test scores, making sure
#   performance is good across all classes (since we optimized for accuracy)
labels_train = grid.best_estimator_.predict(x_train)
labels_test = grid.best_estimator_.predict(x_test)
print('Training report:\n', classification_report(y_train, labels_train))
print('Test report:\n', classification_report(y_test, labels_test))

#-------------------------------------------------------------------------------
#   Try SVM (linear and non-linear kernels)
#-------------------------------------------------------------------------------

tuned_parameters = [
    {
        'svc__kernel': ['linear', 'rbf'],
        'svc__gamma': np.logspace(-5, -1, 10),
        'svc__C': np.logspace(0, 4, 5)
    }
]

pipe = make_pipeline(
    StandardScaler(),
    svm.SVC(random_state = random_seed)
)

grid = GridSearchCV(pipe, tuned_parameters, cv=5, scoring = 'accuracy')
grid.fit(x_train, y_train)

acc_train = round(100 * grid.best_estimator_.score(x_train, y_train), 1)
acc_test = round(100 * grid.best_estimator_.score(x_test, y_test), 1)
print(f'SVM training accuracy: {acc_train}%.')
print(f'SVM test accuracy: {acc_test}%.')
print(f'Best params: {grid.best_params_}')

#   Print a more thorough report about training and test scores, making sure
#   performance is good across all classes (since we optimized for accuracy)
labels_train = grid.best_estimator_.predict(x_train)
labels_test = grid.best_estimator_.predict(x_test)
print('Training report:\n', classification_report(y_train, labels_train))
print('Test report:\n', classification_report(y_test, labels_test))
