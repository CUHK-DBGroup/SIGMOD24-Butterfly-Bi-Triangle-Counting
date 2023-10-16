from sklearn import datasets
from sklearn.model_selection import train_test_split
from sklearn.metrics import precision_score, recall_score, f1_score
import numpy as np
from sklearn import svm
from xgboost import XGBClassifier
in_file = open("xx.txt","r") # input 4-graphlet feature from xx.txt
import time

graph_data = []
graph_target = []
for i in in_file.readlines():
    nm,three_path,bfc,left_clique,right_clique = i.split()
    three_path = float(three_path)
    bfc = float(bfc)
    left_clique = float(left_clique)
    right_clique = float(right_clique)
    sum_graphlets = three_path+bfc+left_clique+right_clique

    graph_data.append([three_path/sum_graphlets,bfc/sum_graphlets,left_clique/sum_graphlets,right_clique/sum_graphlets])

    if "authorship" in nm: # strategy
        graph_target.append(1)
    else:
        graph_target.append(0)


feature_train, feature_test, target_train, target_test = train_test_split(graph_data, graph_target, test_size=0.2, random_state=11)
graph_data = np.array(graph_data)
graph_target = np.array(graph_target)
time_start=time.perf_counter()

svm_classifier = svm.SVC(C=1.0, kernel='linear', decision_function_shape='ovo', gamma=0.01)

svm_classifier.fit(feature_train, target_train)
time_end=time.perf_counter()
print('SVM time cost',time_end-time_start,'s')
print("Train dataset:", svm_classifier.score(feature_train, target_train))
print("Test dataset:", svm_classifier.score(feature_test, target_test))
print("f1_score", f1_score(svm_classifier.predict(feature_test), target_test))
time_start=time.perf_counter()

xgboostModel = XGBClassifier(n_estimators=100, learning_rate= 0.3)

xgboostModel.fit(feature_train, target_train)
time_end=time.perf_counter()
print('XGboost time cost',time_end-time_start,'s')
print('Train dataset: ',xgboostModel.score(feature_train,target_train))
print('Test dataset:',xgboostModel.score(feature_test,target_test))
print("f1_score", f1_score(xgboostModel.predict(feature_test), target_test))