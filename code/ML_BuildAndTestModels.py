#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 16:32:45 2019

@author: martinm
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from xgboost import XGBRegressor, XGBClassifier
from sklearn import metrics
import pickle
import shap
from sklearn.metrics import auc, precision_recall_curve, average_precision_score

sns.set_context('talk')

fsdesforml=pd.read_pickle('./dataframes/ml/fsdesforml.pkl')
fsdesmlx=pd.read_pickle('./dataframes/ml/fsdesmlx.pkl')
fsdestestx=pd.read_pickle('./dataframes/ml/fsdestestx.pkl')
fsdesmly=pd.read_pickle('./dataframes/ml/fsdesmly.pkl')
fsdestesty=pd.read_pickle('./dataframes/ml/fsdestesty.pkl')

hivdbvarforml=pd.read_pickle('./dataframes/ml/hivdbvarforml.pkl')
hivmlx=pd.read_pickle('./dataframes/ml/hivmlx.pkl')
hivtestx=pd.read_pickle('./dataframes/ml/hivtestx.pkl')
hivmly=pd.read_pickle('./dataframes/ml/hivmly.pkl')
hivtesty=pd.read_pickle('./dataframes/ml/hivtesty.pkl')

#####
# Combine features #
#####

aa=pd.read_pickle('./ml/features/aa_allfeaturesenc_fsdes.pkl')
sec=pd.read_pickle('./ml/features/sec_features_fsdesforml.pkl')
tai=pd.read_pickle('./ml/features/tai_features_fsdes.pkl')
slipenc=pd.read_pickle('./ml/features/slipenc_fsdesforml.pkl')

taicols=[str(x)+'_frame0' for x in np.arange(12,26)]+[str(x)+'_framem1' for x in np.arange(12,26)]+[str(x)+'_diffframem1to0' for x in np.arange(12,26)]
aacols=['charged'+str(x) for x in np.arange(10, 14)]+['polar'+str(x) for x in np.arange(10, 14)]+['unpolar'+str(x) for x in np.arange(10, 14)]

seccols=np.arange(42,162)
dgcols=[x for x in fsdesforml.columns if ('dg' in x)&('up' not in x)]
dgcols.remove('dg120')
slipcols=['A','C','G','N','T','Asecond','Csecond','Gsecond','Nsecond','Tsecond','Alast','Clast','Glast','Nlast','Tlast']


feat=tai.join(aa)
feat=feat.join(fsdesforml[[x for x in fsdesforml.columns if 'dg' in x]])
feat=feat.join(sec)
feat=feat.join(slipenc)
feat.to_pickle('./ml/features/all_features_fsdes.pkl')


taihiv=pd.read_pickle('./ml/features/tai_features_hiv.pkl')
aahiv=pd.read_pickle('./ml/features/aa_allfeaturesenc_hiv.pkl')
sechiv=pd.read_pickle('./ml/features/sec_features_hivdbvarforml.pkl')

aacolsforhiv=[x for x in aacols if x in aahiv.columns]

feathiv=hivdbvarforml[[x for x in hivdbvarforml.columns if 'dg' in x]].join(sechiv[seccols])
feathiv=feathiv.join(aahiv[aacolsforhiv])
feathiv=feathiv.join(taihiv[taicols])
feathiv.to_pickle('./ml/features/all_features_hiv.pkl')

####
feat=pd.read_pickle('./ml/features/all_features_fsdes.pkl')
feathiv=pd.read_pickle('./ml/features/all_features_hiv.pkl')
####

wtvalues=pd.read_pickle('./dataframes/wtvalues.pkl')

thrshld=1.296

eventlist=list(wtvalues[wtvalues.gfpm>thrshld].index.values)                     
plist=list(wtvalues[wtvalues.gfpm>thrshld].index.values)                     
nofslist=list(wtvalues[(wtvalues.gfpm<thrshld)&(wtvalues.wavp<thrshld)].index.values)   

#%%
#%% functions

from sklearn.model_selection import KFold, StratifiedKFold

from sklearn.model_selection import RandomizedSearchCV

###################

def ronhivtest_gfp(XY_train, XY_test, cols={'all':feat.columns,'tai':taicols,'aa':aacols,'sec':seccols,'dg':dgcols, 'slip':slipcols}):
    rtest=pd.DataFrame()
    for event in r2scoresoptpeaks1min20.index:
        for col in r2scoresoptpeaks1min20.columns:
            try:
                if (event=='eventlist'):
                    clf=XGBRegressor(**eval(r2scoresoptcondpeaks1min20.at[event, col]))
                    clf.fit(feat.loc[XY_train[[x in eventlist for x in XY_train.fsevent]].index, cols[col]], 
                            XY_train[[x in eventlist for x in XY_train.fsevent]].gfpm_wt)
                    ypred=clf.predict(feathiv.loc[XY_test.index, cols[col]]) 
                    plt.figure(figsize=(3,3))
                    plt.scatter(XY_test.gfpm_wt, ypred, alpha=0.3)
                    plt.title(event + ' ' + col)
                    rtest.loc[event,col]=pearsonr(XY_test.gfpm_wt, ypred)[0]
                elif (event=='full'):
                    clf=XGBRegressor(**eval(r2scoresoptcondpeaks1min20.at[event, col]))
                    clf.fit(feat.loc[XY_train.index, cols[col]], 
                            XY_train.gfpm_wt)
                    ypred=clf.predict(feathiv.loc[XY_test.index, cols[col]]) 
                    plt.figure(figsize=(3,3))
                    plt.scatter(XY_test.gfpm_wt, ypred, alpha=0.3)
                    plt.title(event + ' ' + col)
                    rtest.loc[event,col]=pearsonr(XY_test.gfpm_wt, ypred)[0]
                else:
                    clf=XGBRegressor(**eval(r2scoresoptcondpeaks1min20.at[event, col]))
                    clf.fit(feat.loc[XY_train[(XY_train.fsevent==' '.join(event.split('_')))].index, cols[col]], 
                                              XY_train[(XY_train.fsevent==' '.join(event.split('_')))].gfpm_wt)
                    ypred=clf.predict(feathiv.loc[XY_test.index, cols[col]]) 
                    plt.figure(figsize=(3,3))
                    plt.scatter(XY_test.gfpm_wt, ypred, alpha=0.3)
                    plt.title(event + ' ' + col)
                    rtest.loc[event,col]=pearsonr(XY_test.gfpm_wt, ypred)[0]
            except:
                pass
    return rtest



def rontest_crosspred_allfeatures(XY_train, XY_test, featureset):
    cols={'all':taicols+aacols+list(seccols)+dgcols+list(slipcols),'tai':taicols,'aa':aacols,'sec':seccols,'dg':dgcols, 'slip':slipcols}
    rtest=pd.DataFrame()
    for event in eventlist:
        for event2 in eventlist:
            try:
                col=cols[featureset]
                clf=XGBRegressor(**eval(r2scoresoptcondpeaks1min20.at['_'.join(event.split(' ')), featureset]))
                clf.fit(feat.loc[XY_train[(XY_train.fsevent==' '.join(event.split('_')))].index, col], 
                                          XY_train[(XY_train.fsevent==event)].gfpm_wt)
                ypred=clf.predict(feat.loc[XY_test[(XY_test.fsevent==event2)].index, col]) 
                plt.figure(figsize=(3,3))
                plt.scatter(XY_test[(XY_test.fsevent==event2)].gfpm_wt, ypred, alpha=0.3)
                plt.title('train: ' + event + ', test: ' + event2)
                rtest.loc[event,event2]=pearsonr(XY_test[(XY_test.fsevent==event2)].gfpm_wt, ypred)[0]
            except:
                pass
    return rtest


def rontest_class_crosspred_allfeatures(XY_train, XY_test, featureset):
    cols={'all':taicols+aacols+list(seccols)+dgcols+list(slipcols),'tai':taicols,'aa':aacols,'sec':seccols,'dg':dgcols, 'slip':slipcols}
    rtest=pd.DataFrame()
    for event in eventlist:
        for event2 in eventlist:
            try:
                col=cols[featureset]
                clf=XGBClassifier(**eval(r2scoresoptcondclass.at['_'.join(event.split(' ')), featureset]))
                clf.fit(feat.loc[XY_train[(XY_train.fsevent==' '.join(event.split('_')))].index, col], 
                                          XY_train[(XY_train.fsevent==event)].shifting)
                Xs_test=feat.loc[XY_test[(XY_test.fsevent==event2)].index, col]
                Ys_test=XY_test[(XY_test.fsevent==event2)].shifting
                ypreds=clf.predict_proba(Xs_test) 
                ypred=[x[1] for x in ypreds]
                a, b, _ = metrics.roc_curve(Ys_test, ypred)
                
                f=plt.figure(figsize=(3,3))
                plt.plot(a,b)
                plt.plot([0,1],[0,1], '--',color='gray', linewidth=2)
                plt.xlabel('false positive rate')
                plt.ylabel('true positive rate')
                plt.title('ROC AUC={:.2g}'.format(metrics.roc_auc_score(Ys_test, ypred)), fontsize=14)
                f.savefig('./ml/gbrclass/roc_curve_ontest_allfeatures_eventlist.png',
                          dpi = 300, format='png', bbox_inches='tight', frameon=True)
                rtest.loc[event,event2]=metrics.roc_auc_score(XY_test[(XY_test.fsevent==event2)].shifting, ypred)
            except:
                pass
    return rtest


def rontestperc_crosspred_allfeatures(XY_train, XY_test, featureset):
    cols={'all':taicols+aacols+list(seccols)+dgcols+list(slipcols),'tai':taicols,'aa':aacols,'sec':seccols,'dg':dgcols, 'slip':slipcols}
    rtest=pd.DataFrame()
    for event in eventlist:
        for event2 in eventlist:
            try:
                col=cols[featureset]
                clf=XGBRegressor(**eval(r2scoresoptcondpeaks1min20.at['_'.join(event.split(' ')), featureset]))
                clf.fit(feat.loc[XY_train[(XY_train.fsevent==' '.join(event.split('_')))].index, col], 
                                          XY_train[(XY_train.fsevent==event)].percgfpm_wt)
                ypred=clf.predict(feat.loc[XY_test[(XY_test.fsevent==event2)].index, col]) 
                plt.figure(figsize=(3,3))
                plt.scatter(XY_test[(XY_test.fsevent==event2)].percgfpm_wt, ypred, alpha=0.3)
                plt.title('train: ' + event + ', test: ' + event2)
                rtest.loc[event,event2]=pearsonr(XY_test[(XY_test.fsevent==event2)].percgfpm_wt, ypred)[0]
            except:
                pass
    return rtest

def randsearchhyperparams(Xs, Ys, outputfilename=None, folds=5, param_comb=1000, scrng='r2',
        params={'max_depth':[2,3,4,5,6,7],
        'learning_rate': np.arange(0.02, 0.32, 0.04), 
        'n_estimators':range(20,520,30)}):
    
    '''
    performs Grid Search for hyperparameters (either the default ones or those provided
    as a dictionary in the input) on the data (Xs, Ys) and returns a dataframe of the results
    and prints the best score and the best parameters. It also writes the results to a csv file 
    if an outputfilename is provided
    '''

    skf=KFold(n_splits=folds, shuffle=True, random_state=0)
    xgb=XGBRegressor()
    grid=RandomizedSearchCV(xgb, params, n_iter=int(param_comb), scoring=scrng, 
                                     cv=skf.split(Xs, Ys), verbose=1)
    
    grid.fit(Xs, Ys)
    print grid.best_score_
    print grid.best_params_ 
    griddf=pd.DataFrame(grid.cv_results_)
    
    if outputfilename!=None:
        f=plt.figure(figsize=(5,4))
        sns.heatmap(griddf.pivot_table(index='param_learning_rate', columns='param_max_depth', values='mean_test_score'))
        plt.ylabel('learning_rate')
        plt.xlabel('max_depth')
        f.savefig('/net/mraid08/export/genie/Runs/Martin/prf_mpra/code/ml/xgboptimization/'+outputfilename+'_maxdepth_vs_learningrate.png',
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
        f=plt.figure(figsize=(5,4))
        sns.heatmap(griddf.pivot_table(index='param_n_estimators', columns='param_max_depth', values='mean_test_score'))
        plt.ylabel('n_estimators')
        plt.xlabel('max_depth')
        f.savefig('/net/mraid08/export/genie/Runs/Martin/prf_mpra/code/ml/xgboptimization/'+outputfilename+'_maxdepth_vs_nestimators.png',
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
        f=plt.figure(figsize=(5,4))
        sns.heatmap(griddf.pivot_table(index='param_learning_rate', columns='param_n_estimators', values='mean_test_score'))
        plt.ylabel('learning_rate')
        plt.xlabel('n_estimators')
        f.savefig('/net/mraid08/export/genie/Runs/Martin/prf_mpra/code/ml/xgboptimization/'+outputfilename+'_learningrate_vs_nestimators.png',
              dpi = 300, format='png', bbox_inches='tight', frameon=True)


    if outputfilename!=None:
        griddf.to_csv('/net/mraid08/export/genie/Runs/Martin/prf_mpra/code/ml/xgboptimization/'+outputfilename+'.csv')
    return griddf, grid.best_score_, grid.best_params_ 



def randsearchhyperparams_class(Xs, Ys, outputfilename=None, folds=5, param_comb=1000, scrng='roc_auc',
        params={'max_depth':[2,3,4,5,6,7],
        'learning_rate': np.arange(0.02, 0.32, 0.04), 
        'n_estimators':range(20,520,30)}):
    
    '''
    performs Grid Search for hyperparameters (either the default ones or those provided
    as a dictionary in the input) on the data (Xs, Ys) and returns a dataframe of the results
    and prints the best score and the best parameters. It also writes the results to a csv file 
    if an outputfilename is provided
    '''

    skf=StratifiedKFold(n_splits=folds, shuffle=True, random_state=0)
    xgb=XGBClassifier()
    grid=RandomizedSearchCV(xgb, params, n_iter=int(param_comb), scoring=scrng, 
                                     cv=skf.split(Xs, Ys), verbose=1)
    
    grid.fit(Xs, Ys)
    print grid.best_score_
    print grid.best_params_ 
    griddf=pd.DataFrame(grid.cv_results_)
    
    if outputfilename!=None:
        f=plt.figure(figsize=(5,4))
        sns.heatmap(griddf.pivot_table(index='param_learning_rate', columns='param_max_depth', values='mean_test_score'))
        plt.ylabel('learning_rate')
        plt.xlabel('max_depth')
        f.savefig('/net/mraid08/export/genie/Runs/Martin/prf_mpra/code/ml/xgboptimization/'+outputfilename+'_maxdepth_vs_learningrate.png',
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
        f=plt.figure(figsize=(5,4))
        sns.heatmap(griddf.pivot_table(index='param_n_estimators', columns='param_max_depth', values='mean_test_score'))
        plt.ylabel('n_estimators')
        plt.xlabel('max_depth')
        f.savefig('/net/mraid08/export/genie/Runs/Martin/prf_mpra/code/ml/xgboptimization/'+outputfilename+'_maxdepth_vs_nestimators.png',
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
        f=plt.figure(figsize=(5,4))
        sns.heatmap(griddf.pivot_table(index='param_learning_rate', columns='param_n_estimators', values='mean_test_score'))
        plt.ylabel('learning_rate')
        plt.xlabel('n_estimators')
        f.savefig('/net/mraid08/export/genie/Runs/Martin/prf_mpra/code/ml/xgboptimization/'+outputfilename+'_learningrate_vs_nestimators.png',
              dpi = 300, format='png', bbox_inches='tight', frameon=True)


    if outputfilename!=None:
        griddf.to_csv('/net/mraid08/export/genie/Runs/Martin/prf_mpra/code/ml/xgboptimization/'+outputfilename+'.csv')
    return griddf, grid.best_score_, grid.best_params_ 


#%%

### hyperparameter optimization (should be parallelized)

## only 1 peak and min 20 reads
cols={'all':taicols+aacols+list(seccols)+dgcols+list(slipcols),'tai':taicols,'aa':aacols,'sec':seccols,'dg':dgcols, 'slip':slipcols}

r2scoresoptpeaks1min20=pd.DataFrame()
r2scoresoptcondpeaks1min20=pd.DataFrame()

for col in cols.keys():
    for event in eventlist:
        _, r2scoresoptpeaks1min20.loc[event, col], conds=\
            randsearchhyperparams(feat.loc[fsdesmlx[(fsdesmlx.fsevent==event)&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)].index, cols[col]], 
            fsdesmlx[(fsdesmlx.fsevent==event)&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)].gfpm_wt, 
            param_comb=100)
        r2scoresoptcondpeaks1min20.loc[event, col]=str(conds)

for col in cols.keys():
    _, r2scoresoptpeaks1min20.loc['eventlist', col], conds=\
            randsearchhyperparams(feat.loc[fsdesmlx[[x in eventlist for x in fsdesmlx.fsevent]&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)].index, cols[col]], 
            fsdesmlx[[x in eventlist for x in fsdesmlx.fsevent]&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)].gfpm_wt, 
            param_comb=100)
    r2scoresoptcondpeaks1min20.loc['eventlist', col]=str(conds)


for col in cols.keys():
    _, r2scoresoptpeaks1min20.loc['full', col], conds=\
            randsearchhyperparams(feat.loc[fsdesmlx[(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)].index, cols[col]], 
            fsdesmlx[(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)].gfpm_wt, 
            param_comb=100)
    r2scoresoptcondpeaks1min20.loc['full', col]=str(conds)

r2scoresoptpeaks1min20.to_csv('./forpredictor/r2scoresoptpeaks1min20.csv')
r2scoresoptcondpeaks1min20.to_csv('./forpredictor/r2scoresoptcondpeaks1min20.csv')

####
r2scoresoptpeaks1min20=pd.read_csv('./forpredictor/r2scoresoptpeaks1min20.csv', index_col='Unnamed: 0')
r2scoresoptcondpeaks1min20=pd.read_csv('./forpredictor/r2scoresoptcondpeaks1min20.csv', index_col='Unnamed: 0')

# Predict on test set and save models

cols={'all':taicols+aacols+list(seccols)+list(dgcols)+list(slipcols),'tai':taicols,'aa':aacols,'sec':seccols,'dg':dgcols, 'slip':slipcols}
rontestsetpeak1min20=pd.DataFrame()
for event in r2scoresoptpeaks1min20.index:
    for col in r2scoresoptpeaks1min20.columns:
        if (event=='eventlist'):
            clf=XGBRegressor(**eval(r2scoresoptcondpeaks1min20.at[event, col]))
            clf.fit(feat.loc[fsdesmlx[[x in eventlist for x in fsdesmlx.fsevent]&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]], 
                    fsdesmlx[[x in eventlist for x in fsdesmlx.fsevent]&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].gfpm_wt)
            pickle.dump(clf, open('./forpredictor/models/eventlist_'+col+'_m1.mdl', 'wb'))
            ypred=clf.predict(feat.loc[fsdestestx[[x in eventlist for x in fsdestestx.fsevent]&(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].index, cols[col]]) 
            rontestsetpeak1min20.loc[event,col]=pearsonr(fsdestestx[[x in eventlist for x in fsdestestx.fsevent]&(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].gfpm_wt, ypred)[0]
        elif (event=='full'):
            clf=XGBRegressor(**eval(r2scoresoptcondpeaks1min20.at[event, col]))
            clf.fit(feat.loc[fsdesmlx[(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]], 
                    fsdesmlx[(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].gfpm_wt)
            pickle.dump(clf, open('./forpredictor/models/full_'+col+'_m1.mdl', 'wb'))
            ypred=clf.predict(feat.loc[fsdestestx[(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].index, cols[col]]) 
            rontestsetpeak1min20.loc[event,col]=pearsonr(fsdestestx[(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].gfpm_wt, ypred)[0]
        else:
            if len(fsdesmlx[(fsdesmlx.fsevent==' '.join(event.split('_')))&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)])>100:
                clf=XGBRegressor(**eval(r2scoresoptcondpeaks1min20.at[event, col]))
                clf.fit(feat.loc[fsdesmlx[(fsdesmlx.fsevent==' '.join(event.split('_')))&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]], 
                                          fsdesmlx[(fsdesmlx.fsevent==' '.join(event.split('_')))&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].gfpm_wt)
                pickle.dump(clf, open('./forpredictor/models/'+event+'_'+col+'_m1.mdl', 'wb'))
                ypred=clf.predict(feat.loc[fsdestestx[(fsdestestx.fsevent==' '.join(event.split('_')))&(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].index, cols[col]]) 
                rontestsetpeak1min20.loc[event,col]=pearsonr(fsdestestx[(fsdestestx.fsevent==' '.join(event.split('_')))&(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].gfpm_wt, ypred)[0]

rontestsetpeak1min20.to_csv('./ml/rontestsetpeak1min20_fsdes_reducedtai_aa.csv')

rontestsetpeak1min20=pd.read_csv('./ml/rontestsetpeak1min20_fsdes_reducedtai_aa.csv').set_index('Unnamed: 0')

f=plt.figure(figsize=(6,4))
sns.heatmap(rontestsetpeak1min20[['slip','dg','sec','aa','tai','all']],
    annot=rontestsetpeak1min20[['slip','dg','sec','aa','tai','all']], 
    annot_kws={'fontsize':12})
f.savefig('./figures/ml/heatmap_r_ontestsetpeak1min20heatmaps_gfpmwt.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

cg=sns.clustermap(data=rontestsetpeak1min20[['slip','dg','sec','aa','tai','all']].replace(to_replace=np.nan, value=0), figsize=(6,5), annot=True, 
                  cbar_kws={'ticks':[-1,-0.5,0,0.5,1]},
               annot_kws={'fontsize':12}, fmt='.2f', col_cluster=False)
plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
cg.savefig('./figures/ml/clustermap_r_ontestsetpeak1min20heatmaps_gfpmwt.png',
              dpi = 300, format='png', bbox_inches='tight', frameon=True)

# Scatter plot for all cases

cols={'all':taicols+aacols+list(seccols)+dgcols+list(slipcols),'tai':taicols,'aa':aacols,'sec':seccols,'dg':dgcols, 'slip':slipcols}
for event in r2scoresoptpeaks1min20.index:
    for col in r2scoresoptpeaks1min20.columns:
        if (event=='eventlist'):
            clf=XGBRegressor(**eval(r2scoresoptcondpeaks1min20.at[event, col]))
            clf.fit(feat.loc[fsdesmlx[[x in eventlist for x in fsdesmlx.fsevent]&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]], 
                    fsdesmlx[[x in eventlist for x in fsdesmlx.fsevent]&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].gfpm_wt)
            ypred=clf.predict(feat.loc[fsdestestx[[x in eventlist for x in fsdestestx.fsevent]&(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].index, cols[col]]) 
            f=plt.figure(figsize=(3,3))
            maxval=np.ceil(np.max(list(fsdestestx[[x in eventlist for x in fsdestestx.fsevent]&(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].gfpm_wt.values) + ypred))
            plt.scatter(fsdestestx[[x in eventlist for x in fsdestestx.fsevent]&(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].gfpm_wt, ypred, alpha=0.3)
            plt.title('Pearson r = {:.2f}'.format(pearsonr(fsdestestx[[x in eventlist for x in fsdestestx.fsevent]&(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].gfpm_wt, ypred)[0])
            +'\np = {:.2g}'.format(pearsonr(fsdestestx[[x in eventlist for x in fsdestestx.fsevent]&(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].gfpm_wt, ypred)[1]), fontsize=14)
            plt.xlim(0, 10)
            plt.ylim(0, 10)
            plt.plot([0, maxval],[0, maxval], '--', linewidth=2, color='grey')
            f.savefig('./figures/ml/prediction_sametraining_on_testset_'+'_'.join(event.split(' '))+'_'+col+'.png',\
                      dpi = 300, format='png', bbox_inches='tight', frameon=True)
        elif (event=='full'):
            clf=XGBRegressor(**eval(r2scoresoptcondpeaks1min20.at[event, col]))
            clf.fit(feat.loc[fsdesmlx[(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]], 
                    fsdesmlx[(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].gfpm_wt)
            ypred=clf.predict(feat.loc[fsdestestx[(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].index, cols[col]]) 
            f=plt.figure(figsize=(3,3))
            maxval=np.ceil(np.max(list(fsdestestx[(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].gfpm_wt.values) + ypred))
            plt.scatter(fsdestestx[(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].gfpm_wt, ypred, alpha=0.3)
            plt.title('Pearson r = {:.2f}'.format(pearsonr(fsdestestx[(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].gfpm_wt, ypred)[0])
            +'\np = {:.2g}'.format(pearsonr(fsdestestx[(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].gfpm_wt, ypred)[1]), fontsize=14)
            plt.xlim(0, 10)
            plt.ylim(0, 10)
            plt.plot([0, maxval],[0, maxval], '--', linewidth=2, color='grey')
            f.savefig('./figures/ml/prediction_sametraining_on_testset_'+'_'.join(event.split(' '))+'_'+col+'.png',\
                      dpi = 300, format='png', bbox_inches='tight', frameon=True)
        else:
            if len(fsdesmlx[(fsdesmlx.fsevent==' '.join(event.split('_')))&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)])>100:
                clf=XGBRegressor(**eval(r2scoresoptcondpeaks1min20.at[event, col]))
                clf.fit(feat.loc[fsdesmlx[(fsdesmlx.fsevent==' '.join(event.split('_')))&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]], 
                                          fsdesmlx[(fsdesmlx.fsevent==' '.join(event.split('_')))&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].gfpm_wt)
                ypred=clf.predict(feat.loc[fsdestestx[(fsdestestx.fsevent==' '.join(event.split('_')))&(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].index, cols[col]]) 
                f=plt.figure(figsize=(3,3))
                maxval=np.ceil(np.max(list(fsdestestx[(fsdestestx.fsevent==' '.join(event.split('_')))&(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].gfpm_wt.values) + ypred))
                plt.scatter(fsdestestx[(fsdestestx.fsevent==' '.join(event.split('_')))&(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].gfpm_wt, ypred, alpha=0.3)
                plt.title('Pearson r = {:.2f}'.format(pearsonr(fsdestestx[(fsdestestx.fsevent==' '.join(event.split('_')))&(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].gfpm_wt, ypred)[0])
                +'\np = {:.2g}'.format(pearsonr(fsdestestx[(fsdestestx.fsevent==' '.join(event.split('_')))&(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].gfpm_wt, ypred)[1]), fontsize=14)
                plt.xlim(0, maxval)
                plt.ylim(0, maxval)
                plt.plot([0, maxval],[0, maxval], '--', linewidth=2, color='grey')
                f.savefig('./figures/ml/prediction_sametraining_on_testset_'+'_'.join(event.split(' '))+'_'+col+'.png',\
                          dpi = 300, format='png', bbox_inches='tight', frameon=True)



# plot shap values
cols={'all':taicols+aacols+list(seccols)+dgcols+list(slipcols),'tai':taicols,'aa':aacols,'sec':seccols,'dg':dgcols, 'slip':slipcols}

for event in r2scoresoptpeaks1min20.index:
    for col in r2scoresoptpeaks1min20.columns:
        if (event=='eventlist'):
            clf=XGBRegressor(**eval(r2scoresoptcondpeaks1min20.at[event, col]))
            clf.fit(feat.loc[fsdesmlx[[x in eventlist for x in fsdesmlx.fsevent]&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]], 
                    fsdesmlx[[x in eventlist for x in fsdesmlx.fsevent]&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].gfpm_wt)
            explainer = shap.TreeExplainer(clf)
            shap_values=explainer.shap_values(feat.loc[fsdesmlx[[x in eventlist for x in fsdesmlx.fsevent]&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]])
            f=plt.figure()
            shap.summary_plot(shap_values, feat.loc[fsdesmlx[[x in eventlist for x in fsdesmlx.fsevent]&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]], feature_names=cols[col])
            f.savefig('./figures/ml/shap_summaryplot_feat_peaks1min20reads_'+col+'_eventlist.png',
                      dpi = 300, format='png', bbox_inches='tight', frameon=True)
        elif (event=='full'):
            clf=XGBRegressor(**eval(r2scoresoptcondpeaks1min20.at[event, col]))
            clf.fit(feat.loc[fsdesmlx[(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]], 
                    fsdesmlx[(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].gfpm_wt)
            explainer = shap.TreeExplainer(clf)
            shap_values=explainer.shap_values(feat.loc[fsdesmlx[(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]])
            f=plt.figure()
            shap.summary_plot(shap_values, feat.loc[fsdesmlx[(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]], feature_names=cols[col])
            f.savefig('./figures/ml/shap_summaryplot_feat_peaks1min20reads_'+col+'_full.png',
                      dpi = 300, format='png', bbox_inches='tight', frameon=True)
        else:
            if len(fsdesmlx[(fsdesmlx.fsevent==' '.join(event.split('_')))&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)])>100:
                clf=XGBRegressor(**eval(r2scoresoptcondpeaks1min20.at[event, col]))
                clf.fit(feat.loc[fsdesmlx[(fsdesmlx.fsevent==' '.join(event.split('_')))&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]], 
                                          fsdesmlx[(fsdesmlx.fsevent==' '.join(event.split('_')))&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].gfpm_wt)
                explainer = shap.TreeExplainer(clf)
                shap_values=explainer.shap_values(feat.loc[fsdesmlx[(fsdesmlx.fsevent==' '.join(event.split('_')))&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]])
                f=plt.figure()
                shap.summary_plot(shap_values, feat.loc[fsdesmlx[(fsdesmlx.fsevent==' '.join(event.split('_')))&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]], feature_names=cols[col])
                f.savefig('./figures/ml/shap_summaryplot_feat_peaks1min20reads_'+col+'_'+'_'.join(event.split(' '))+'.png',
                          dpi = 300, format='png', bbox_inches='tight', frameon=True)



#%%

#%% hivdbvar

cols={'all':[x for x in feathiv.columns if x not in list(slipcols)],'all_no_dg':[x for x in feathiv.columns if x not in list(slipcols)+dgcols],
             'all_no_dg_nosec':[x for x in feathiv.columns if x not in list(slipcols)+dgcols+list(seccols)],'tai':taicols,
             'aa':aacols,'sec':seccols,'dg':dgcols}


r2scoresoptpeaks1min20hiv=pd.Series()
r2scoresoptcondpeaks1min20hiv=pd.Series()

for col in cols.keys():
    _, r2scoresoptpeaks1min20hiv.loc[col], conds=\
            randsearchhyperparams(feathiv.loc[hivmlx[(hivmlx.peaksm_wt==1)&(hivmlx.numberreadsm_wt>=20)].index, cols[col]], 
                       hivmlx[(hivmlx.peaksm_wt==1)&(hivmlx.numberreadsm_wt>=20)].gfpm_wt, param_comb=100)
    r2scoresoptcondpeaks1min20hiv.loc[col]=str(conds)


r2scoresoptpeaks1min20hiv.to_csv('./forpredictor/r2scoresoptpeaks1min20hiv.csv')
r2scoresoptcondpeaks1min20hiv.to_csv('./forpredictor/r2scoresoptcondpeaks1min20hiv.csv')

r2scoresoptpeaks1min20hiv=pd.read_csv('./forpredictor/r2scoresoptpeaks1min20hiv.csv', header=None, index_col=0)
r2scoresoptcondpeaks1min20hiv=pd.read_csv('./forpredictor/r2scoresoptcondpeaks1min20hiv.csv', header=None, index_col=0)


rtest=pd.Series()
for col in r2scoresoptpeaks1min20hiv.index:
    try:
        clf=XGBRegressor(**eval(r2scoresoptcondpeaks1min20hiv.loc[col].iloc[0]))
        clf.fit(feathiv.loc[hivmlx[(hivmlx.gfpm_wt>1.3)&(hivmlx.gfpm_wt<5)&(hivmlx.peaksm_wt==1)&(hivmlx.numberreadsm_wt>=20)].index, cols[col]], 
                hivmlx[(hivmlx.gfpm_wt>1.3)&(hivmlx.gfpm_wt<5)&(hivmlx.peaksm_wt==1)&(hivmlx.numberreadsm_wt>=20)].gfpm_wt)
        pickle.dump(clf, open('./forpredictor/models/hiv_variants_'+col+'.mdl', 'wb'))
        ypred=clf.predict(feathiv.loc[hivtestx[(hivtestx.gfpm_wt>1.3)&(hivtestx.gfpm_wt<5)&(hivtestx.peaksm_wt==1)&(hivtestx.numberreadsm_wt>=20)].index, cols[col]]) 
        f=plt.figure(figsize=(3,3))
        plt.scatter(hivtestx[(hivtestx.gfpm_wt>1.3)&(hivtestx.gfpm_wt<5)&(hivtestx.peaksm_wt==1)&(hivtestx.numberreadsm_wt>=20)].gfpm_wt, ypred, alpha=0.3)
        plt.title('Pearson r = {:.2f}'.format(pearsonr(hivtestx[(hivtestx.gfpm_wt>1.3)&(hivtestx.gfpm_wt<5)&(hivtestx.peaksm_wt==1)&(hivtestx.numberreadsm_wt>=20)].gfpm_wt, ypred)[0])
        +'\np = {:.2g}'.format(pearsonr(hivtestx[(hivtestx.gfpm_wt>1.3)&(hivtestx.gfpm_wt<5)&(hivtestx.peaksm_wt==1)&(hivtestx.numberreadsm_wt>=20)].gfpm_wt, ypred)[1]), fontsize=14)
        plt.xlim(1.5,4.5)
        plt.ylim(1.5,4.5)
        plt.plot([1.5,4.5],[1.5,4.5], '--', linewidth=2, color='grey')
        f.savefig('./figures/ml/prediction_training_hivdbvartraining_test_hivdbvartest_'+col+'.png',\
                  dpi = 300, format='png', bbox_inches='tight', frameon=True)
        rtest.loc[col]=pearsonr(hivtestx[(hivtestx.gfpm_wt>1.3)&(hivtestx.gfpm_wt<5)&(hivtestx.peaksm_wt==1)&(hivtestx.numberreadsm_wt>=20)].gfpm_wt, ypred)[0]
    except:
        pass
    

####################


#%%
# predict on each other's test set

# on all
rontestsetpeak1min20gfp_crosspred=rontest_crosspred_allfeatures(fsdesmlx[(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)],
         fsdestestx[(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)], 'all')
rontestsetpeak1min20gfp_crosspred.to_csv('./ml/rontestsetpeak1min20gfp_crosspred.csv')

f=plt.figure(figsize=(7,3))
sns.heatmap(rontestsetpeak1min20gfp_crosspred.drop(['PRRSV - nsp2TF'
                                                        ]).drop(['PRRSV - nsp2TF'], axis=1),
    annot=rontestsetpeak1min20gfp_crosspred.drop(['PRRSV - nsp2TF'
                                                        ]).drop(['PRRSV - nsp2TF'], axis=1), 
    annot_kws={'fontsize':12})
f.savefig('./figures/ml/crosspred_all_heatmaps_gfp.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)



# on slippery site
rontestsetpeak1min20gfp_crosspredslip=rontest_crosspred_allfeatures(fsdesmlx[(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)],
         fsdestestx[(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)], 'slip')
rontestsetpeak1min20gfp_crosspredslip.to_csv('./ml/rontestsetpeak1min20gfp_crosspred_slip.csv')

f=plt.figure(figsize=(7,3))
sns.heatmap(rontestsetpeak1min20gfp_crosspredslip.drop(['PRRSV - nsp2TF', 'west nile virus'
                                                        ]).drop(['PRRSV - nsp2TF', 'west nile virus'], axis=1),
    annot=rontestsetpeak1min20gfp_crosspredslip.drop(['PRRSV - nsp2TF', 'west nile virus'
                                                        ]).drop(['PRRSV - nsp2TF', 'west nile virus'], axis=1), 
    annot_kws={'fontsize':12})
f.savefig('./figures/ml/crosspred_slip_heatmaps_gfp.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

    

# on dg
rontestsetpeak1min20gfp_crosspreddg=rontest_crosspred_allfeatures(fsdesmlx[(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)],
         fsdestestx[(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)], 'dg')
rontestsetpeak1min20gfp_crosspreddg.to_csv('./ml/rontestsetpeak1min20gfp_crosspred_dg.csv')

f=plt.figure(figsize=(7,3))
sns.heatmap(rontestsetpeak1min20gfp_crosspreddg.drop(['PRRSV - nsp2TF'
                                                        ]).drop(['PRRSV - nsp2TF'], axis=1),
    annot=rontestsetpeak1min20gfp_crosspreddg.drop(['PRRSV - nsp2TF'
                                                        ]).drop(['PRRSV - nsp2TF'], axis=1), 
    annot_kws={'fontsize':12})
f.savefig('./figures/ml/crosspred_dg_heatmaps_gfp.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


#%%
### Classification

cols={'all':taicols+aacols+list(seccols)+list(dgcols)+list(slipcols),'tai':taicols,'aa':aacols,'sec':seccols,'dg':dgcols, 'slip':slipcols}



r2scoresoptclass=pd.DataFrame()
r2scoresoptcondclass=pd.DataFrame()

for col in cols.keys():
    for event in eventlist:
        _, r2scoresoptclass.loc[event, col], conds=\
            randsearchhyperparams_class(feat.loc[fsdesmlx[(fsdesmlx.fsevent==event)&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)].index, cols[col]], 
            fsdesmlx[(fsdesmlx.fsevent==event)&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)].shifting, 
            param_comb=500)
        r2scoresoptcondclass.loc[event, col]=str(conds)

for col in cols.keys():
    _, r2scoresoptclass.loc['eventlist', col], conds=\
            randsearchhyperparams_class(feat.loc[fsdesmlx[[x in eventlist for x in fsdesmlx.fsevent]&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)].index, cols[col]], 
            fsdesmlx[[x in eventlist for x in fsdesmlx.fsevent]&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)].shifting, 
            param_comb=500)
    r2scoresoptcondclass.loc['eventlist', col]=str(conds)


for col in cols.keys():
    _, r2scoresoptclass.loc['full', col], conds=\
            randsearchhyperparams_class(feat.loc[fsdesmlx[(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)].index, cols[col]], 
            fsdesmlx[(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)].shifting, 
            param_comb=500)
    r2scoresoptcondclass.loc['full', col]=str(conds)


r2scoresoptcondclass.to_csv('./ml/xgboptimization/r2scoresoptcondclass_peaks1min20.csv')

r2scoresoptcondclass=pd.read_csv('./forpredictor/r2scoresoptcondclass_peaks1min20.csv', index_col='Unnamed: 0')

rontestset=pd.DataFrame()
for event in r2scoresoptcondclass.index:
    for col in r2scoresoptcondclass.columns:
        if (event=='full'):
            clf=XGBClassifier(**eval(r2scoresoptcondclass.at[event, col]))
            clf.fit(feat.loc[fsdesmlx[(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]], 
                                      fsdesmlx[(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].shifting)
            pickle.dump(clf, open('./forpredictor/models/classifier_full_'+col+'_m1.mdl', 'wb'))
            ypreds=clf.predict_proba(feat.loc[fsdestestx[(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].index, cols[col]]) 
            ypred=[x[1] for x in ypreds]
            rontestset.loc[event,col]=metrics.roc_auc_score(fsdestestx[(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].shifting, ypred)
        elif (event=='eventlist'):
            clf=XGBClassifier(**eval(r2scoresoptcondclass.at[event, col]))
            clf.fit(feat.loc[fsdesmlx[[x in eventlist for x in fsdesmlx.fsevent]&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]], 
                    fsdesmlx[[x in eventlist for x in fsdesmlx.fsevent]&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].shifting)
            pickle.dump(clf, open('./forpredictor/models/classifier_eventlist_'+col+'_m1.mdl', 'wb'))
            ypreds=clf.predict_proba(feat.loc[fsdestestx[[x in eventlist for x in fsdestestx.fsevent]&(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].index, cols[col]]) 
            ypred=[x[1] for x in ypreds]
            rontestset.loc[event,col]=metrics.roc_auc_score(fsdestestx[[x in eventlist for x in fsdestestx.fsevent]&(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].shifting, ypred)
        else:
            if len(fsdesmlx[(fsdesmlx.fsevent==' '.join(event.split('_')))&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)])>100:
                clf=XGBClassifier(**eval(r2scoresoptcondclass.at[event, col]))
                clf.fit(feat.loc[fsdesmlx[(fsdesmlx.fsevent==' '.join(event.split('_')))&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]], 
                                          fsdesmlx[(fsdesmlx.fsevent==' '.join(event.split('_')))&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].shifting)
                pickle.dump(clf, open('./forpredictor/models/classifier_'+event+'_'+col+'_m1.mdl', 'wb'))
                ypreds=clf.predict_proba(feat.loc[fsdestestx[(fsdestestx.fsevent==' '.join(event.split('_')))&(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].index, cols[col]]) 
                ypred=[x[1] for x in ypreds]
                rontestset.loc[event,col]=metrics.roc_auc_score(fsdestestx[(fsdestestx.fsevent==' '.join(event.split('_')))&(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].shifting, ypred)


f=plt.figure(figsize=(6,4))
sns.heatmap(rontestset[['slip','dg','sec','aa','tai','all']],
    annot=rontestset[['slip','dg','sec','aa','tai','all']], 
    annot_kws={'fontsize':12})
f.savefig('./figures/ml/heatmap_rocauc_ontestsetpeak1min20heatmaps_shifting.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


cg=sns.clustermap(data=rontestset[['slip','dg','sec','aa','tai','all']], figsize=(6,5), annot=True, 
                  cbar_kws={'ticks':[-1,-0.5,0,0.5,0.6,0.7,0.8,1]},
               annot_kws={'fontsize':12}, fmt='.2f', col_cluster=False)
plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
cg.savefig('./figures/ml/clustermap_rocauc_ontestsetpeak1min20heatmaps_shifting.png',
              dpi = 300, format='png', bbox_inches='tight', frameon=True)

# Plot ROC curves

col='all'
for event in r2scoresoptcondclass.index:
    if (event=='full'):
        cols={'all':taicols+aacols+list(seccols)+list(dgcols)+list(slipcols),'tai':taicols,'aa':aacols,'sec':seccols,'dg':dgcols, 'slip':slipcols}
        clf=XGBClassifier(**eval(r2scoresoptcondclass.at[event, col]))
        clf.fit(feat.loc[fsdesmlx[(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]], 
                                  fsdesmlx[(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].shifting)
        Xs_test=feat.loc[fsdestestx[(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].index, cols[col]]
        Ys_test=fsdestestx[(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].shifting
        ypreds=clf.predict_proba(Xs_test) 
        ypred=[x[1] for x in ypreds]
        a, b, _ = metrics.roc_curve(Ys_test, ypred)
        
        f=plt.figure(figsize=(3,3))
        plt.plot(a,b)
        plt.plot([0,1],[0,1], '--',color='gray', linewidth=2)
        plt.xlabel('false positive rate')
        plt.ylabel('true positive rate')
        plt.title('ROC AUC={:.2g}'.format(metrics.roc_auc_score(Ys_test, ypred)), fontsize=14)
        f.savefig('./figures/ml/roc_curve_ontest_allfeatures_full.png',
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    elif (event=='eventlist'):
        cols={'all':taicols+aacols+list(seccols)+list(dgcols)+list(slipcols),'tai':taicols,'aa':aacols,'sec':seccols,'dg':dgcols, 'slip':slipcols}
        clf=XGBClassifier(**eval(r2scoresoptcondclass.at[event, col]))
        clf.fit(feat.loc[fsdesmlx[[x in eventlist for x in fsdesmlx.fsevent]&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]], 
                fsdesmlx[[x in eventlist for x in fsdesmlx.fsevent]&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].shifting)
        Xs_test=feat.loc[fsdestestx[[x in eventlist for x in fsdestestx.fsevent]&(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].index, cols[col]]
        Ys_test=fsdestestx[[x in eventlist for x in fsdestestx.fsevent]&(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].shifting
        ypreds=clf.predict_proba(Xs_test) 
        ypred=[x[1] for x in ypreds]
        a, b, _ = metrics.roc_curve(Ys_test, ypred)
        
        f=plt.figure(figsize=(3,3))
        plt.plot(a,b)
        plt.plot([0,1],[0,1], '--',color='gray', linewidth=2)
        plt.xlabel('false positive rate')
        plt.ylabel('true positive rate')
        plt.title('ROC AUC={:.2g}'.format(metrics.roc_auc_score(Ys_test, ypred)), fontsize=14)
        f.savefig('./figures/ml/roc_curve_ontest_allfeatures_eventlist.png',
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    else:
        cols={'all':taicols+aacols+list(seccols)+list(dgcols)+list(slipcols),'tai':taicols,'aa':aacols,'sec':seccols,'dg':dgcols, 'slip':slipcols}
        clf=XGBClassifier(**eval(r2scoresoptcondclass.at[event, col]))
        clf.fit(feat.loc[fsdesmlx[(fsdesmlx.fsevent==' '.join(event.split('_')))&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]], 
                                  fsdesmlx[(fsdesmlx.fsevent==' '.join(event.split('_')))&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].shifting)
        Xs_test=feat.loc[fsdestestx[(fsdestestx.fsevent==' '.join(event.split('_')))&(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].index, cols[col]]
        Ys_test=fsdestestx[(fsdestestx.fsevent==' '.join(event.split('_')))&(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].shifting
        ypreds=clf.predict_proba(Xs_test) 
        ypred=[x[1] for x in ypreds]
        a, b, _ = metrics.roc_curve(Ys_test, ypred)
        
        f=plt.figure(figsize=(3,3))
        plt.plot(a,b)
        plt.plot([0,1],[0,1], '--',color='gray', linewidth=2)
        plt.xlabel('false positive rate')
        plt.ylabel('true positive rate')
        plt.title('ROC AUC={:.2g}'.format(metrics.roc_auc_score(Ys_test, ypred)), fontsize=14)
        f.savefig('./figures/ml/roc_curve_ontest_allfeatures_'+event+'.png',
              dpi = 300, format='png', bbox_inches='tight', frameon=True)


cols={'all':taicols+aacols+list(seccols)+list(dgcols)+list(slipcols),'tai':taicols,'aa':aacols,'sec':seccols,'dg':dgcols, 'slip':slipcols}
col='all'
event='eventlist'
clf=XGBClassifier(**eval(r2scoresoptcondclass.at[event, col]))
clf.fit(feat.loc[fsdesmlx[(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]], 
                          fsdesmlx[(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].shifting)
Xs_test=feat.loc[fsdestestx[(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].index, cols[col]]
Ys_test=fsdestestx[(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].shifting

prec, rec, thresh=precision_recall_curve(Ys_test, [x[1] for x in clf.predict_proba(Xs_test)])

f=plt.figure(figsize=(3,2))
plt.plot(rec,prec)
plt.xlabel('recall')
plt.ylabel('precision')
plt.title('auPR='+'{:.2g}'.format(auc(rec,prec)), fontsize=12)
f.savefig('./figures/ml/prec_rec_curve_ontest_allfeatures_eventlist.png',
      dpi = 300, format='png', bbox_inches='tight', frameon=True)

average_precision_score(Ys_test, [x[1] for x in clf.predict_proba(Xs_test)])
# 0.7995829679773897
auc(rec,prec)
# 0.7995829679773897

col='all'
event='full'
clf=XGBClassifier(**eval(r2scoresoptcondclass.at[event, col]))
clf.fit(feat.loc[fsdesmlx[(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]], 
                          fsdesmlx[(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].shifting)
Xs_test=feat.loc[fsdestestx[(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].index, cols[col]]
Ys_test=fsdestestx[(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)].shifting

prec, rec, thresh=precision_recall_curve(Ys_test, [x[1] for x in clf.predict_proba(Xs_test)])

f=plt.figure(figsize=(3,2))
plt.plot(rec,prec)
plt.xlabel('recall')
plt.ylabel('precision')
plt.title('auPR='+'{:.2g}'.format(auc(rec,prec)), fontsize=12)
f.savefig('./figures/ml/prec_rec_curve_ontest_allfeatures_full.png',
      dpi = 300, format='png', bbox_inches='tight', frameon=True)

average_precision_score(Ys_test, [x[1] for x in clf.predict_proba(Xs_test)])
# 0.80349298365391086
auc(rec,prec)
# 0.80349298365391086


# plot shap values
cols={'all':taicols+aacols+list(seccols)+list(dgcols)+list(slipcols),'tai':taicols,'aa':aacols,'sec':seccols,'dg':dgcols, 'slip':slipcols}

for event in r2scoresoptcondclass.index:
    for col in r2scoresoptcondclass.columns:
        if (event=='eventlist'):
            clf=XGBClassifier(**eval(r2scoresoptcondclass.at[event, col]))
            clf.fit(feat.loc[fsdesmlx[[x in eventlist for x in fsdesmlx.fsevent]&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]], 
                    fsdesmlx[[x in eventlist for x in fsdesmlx.fsevent]&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].shifting)
            explainer = shap.TreeExplainer(clf)
            shap_values=explainer.shap_values(feat.loc[fsdesmlx[[x in eventlist for x in fsdesmlx.fsevent]&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]])
            f=plt.figure()
            shap.summary_plot(shap_values, feat.loc[fsdesmlx[[x in eventlist for x in fsdesmlx.fsevent]&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]], feature_names=cols[col])
            f.savefig('./figures/ml/shap_summaryplot_feat_class_peaks1min20reads_'+col+'_eventlist.png',
                      dpi = 300, format='png', bbox_inches='tight', frameon=True)
        elif (event=='full'):
            clf=XGBClassifier(**eval(r2scoresoptcondclass.at[event, col]))
            clf.fit(feat.loc[fsdesmlx[(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]], 
                                      fsdesmlx[(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].shifting)
            explainer = shap.TreeExplainer(clf)
            shap_values=explainer.shap_values(feat.loc[fsdesmlx[(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]])
            f=plt.figure()
            shap.summary_plot(shap_values, feat.loc[fsdesmlx[(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]], feature_names=cols[col])
            f.savefig('./figures/ml/shap_summaryplot_feat_class_peaks1min20reads_'+col+'_full.png',
                      dpi = 300, format='png', bbox_inches='tight', frameon=True)
        else:
            clf=XGBClassifier(**eval(r2scoresoptcondclass.at[event, col]))
            clf.fit(feat.loc[fsdesmlx[(fsdesmlx.fsevent==' '.join(event.split('_')))&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]], 
                                      fsdesmlx[(fsdesmlx.fsevent==' '.join(event.split('_')))&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].shifting)
            explainer = shap.TreeExplainer(clf)
            shap_values=explainer.shap_values(feat.loc[fsdesmlx[(fsdesmlx.fsevent==' '.join(event.split('_')))&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]])
            f=plt.figure()
            shap.summary_plot(shap_values, feat.loc[fsdesmlx[(fsdesmlx.fsevent==' '.join(event.split('_')))&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]], feature_names=cols[col])
            f.savefig('./figures/ml/shap_summaryplot_feat_class_peaks1min20reads_'+col+'_'+'_'.join(event.split(' '))+'.png',
                      dpi = 300, format='png', bbox_inches='tight', frameon=True)




event='full'
col='dg'
if (event=='eventlist'):
    clf=XGBClassifier(**eval(r2scoresoptcondclass.at[event, col]))
    clf.fit(feat.loc[fsdesmlx[[x in eventlist for x in fsdesmlx.fsevent]&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]], 
            fsdesmlx[[x in eventlist for x in fsdesmlx.fsevent]&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].shifting)
    explainer = shap.TreeExplainer(clf)
    shap_values=explainer.shap_values(feat.loc[fsdesmlx[[x in eventlist for x in fsdesmlx.fsevent]&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]])
    Xs=feat.loc[fsdesmlx[[x in eventlist for x in fsdesmlx.fsevent]&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]]
elif (event=='full'):
    clf=XGBClassifier(**eval(r2scoresoptcondclass.at[event, col]))
    clf.fit(feat.loc[fsdesmlx[(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]], 
                              fsdesmlx[(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].shifting)
    explainer = shap.TreeExplainer(clf)
    shap_values=explainer.shap_values(feat.loc[fsdesmlx[(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]])
    Xs=feat.loc[fsdesmlx[(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]]
else:
    clf=XGBClassifier(**eval(r2scoresoptcondclass.at[event, col]))
    clf.fit(feat.loc[fsdesmlx[(fsdesmlx.fsevent==' '.join(event.split('_')))&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]], 
                              fsdesmlx[(fsdesmlx.fsevent==' '.join(event.split('_')))&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].shifting)
    explainer = shap.TreeExplainer(clf)
    shap_values=explainer.shap_values(feat.loc[fsdesmlx[(fsdesmlx.fsevent==' '.join(event.split('_')))&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]])
    Xs=feat.loc[fsdesmlx[(fsdesmlx.fsevent==' '.join(event.split('_')))&(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)].index, cols[col]]


for j in range(10,110,10):
    f=plt.figure(figsize=(4,3))
    shap.dependence_plot("dg30",shap_values,Xs, show=False, interaction_index="dg"+str(j))
    f.savefig('./figures/ml/shap_dependenceplot_'+col+'_'+'_'.join(event.split(' '))+'dg30_and_dg'+str(j)+'.png',
              dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(4,3))
shap.dependence_plot("dg60",shap_values,Xs, show=False, interaction_index="dg_downstream")
f.savefig('./figures/ml/shap_dependenceplot_class_'+col+'_'+'_'.join(event.split(' '))+'dg60_and_dgdownstream.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(4,3))
shap.dependence_plot("dg_downstream",shap_values,Xs, show=False, interaction_index="dg60")
f.savefig('./figures/ml/shap_dependenceplot_class_'+col+'_'+'_'.join(event.split(' '))+'dgdownstream_and_dg60.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(4,3))
plt.scatter(fsdesforml.dg_downstream, fsdesforml.rnaperdna, alpha=0.1)
plt.xlabel('MFE (120 nt downstream region)')
plt.ylabel('log2(RNA/DNA reads)')
plt.title('Pearson r='+'{:.2g}'.format(pearsonr(fsdesforml[['dg_downstream','rnaperdna']].replace(to_replace=-np.inf, value=np.nan).dropna().dg_downstream, 
                   fsdesforml[['dg_downstream','rnaperdna']].replace(to_replace=-np.inf, value=np.nan).dropna().rnaperdna)[0])+
    ', p='+'{:.2g}'.format(pearsonr(fsdesforml[['dg_downstream','rnaperdna']].replace(to_replace=-np.inf, value=np.nan).dropna().dg_downstream, 
                   fsdesforml[['dg_downstream','rnaperdna']].replace(to_replace=-np.inf, value=np.nan).dropna().rnaperdna)[1]), fontsize=12)
f.savefig('./figures/ml/correlation_dgdownstream_rnaperdna.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


# predict on each other's test set

# on all
rontestsetclass_crosspred=rontest_class_crosspred_allfeatures(fsdesmlx[(fsdesmlx.peaksm_wt==1)&(fsdesmlx.numberreadsm_wt>=20)&(fsdesmlx.laststopinframem1<12)],
         fsdestestx[(fsdestestx.peaksm_wt==1)&(fsdestestx.numberreadsm_wt>=20)&(fsdestestx.laststopinframem1<12)], 'all')
rontestsetclass_crosspred.to_csv('./ml/rontestsetclass_crosspred.csv')


f=plt.figure(figsize=(7,3))
sns.heatmap(rontestsetclass_crosspred.drop(['PRRSV - nsp2TF'
                                                        ]).drop(['PRRSV - nsp2TF'], axis=1),
    annot=rontestsetclass_crosspred.drop(['PRRSV - nsp2TF'
                                                        ]).drop(['PRRSV - nsp2TF'], axis=1), 
    annot_kws={'fontsize':12})
f.savefig('./figures/ml/rontestsetclass_crosspred_heatmaps_shifting.png',
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

