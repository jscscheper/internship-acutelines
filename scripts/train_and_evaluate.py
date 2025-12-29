"""
This module contains functions that were used in the project
regarding the discovery of mitochondria-related genes. With
this module, we tried through feature selection efforts by
using RFECV and our own functions to identify these
biomarkers. In addition, there are train, cross-validation,
and various visualization functions here to accomplish the
development of a prediction model with these biomarkers.
Features include 210 different mitochondria-related genes
and 3 class variables, namely a sepsis severity variable
('sepsis_severity'), an endotype ('cluster'), and a 
combination set between High + Intermediate vs. Low 
phenotypes ('HighInt_Low').

Usage of functions is highly encouraged in the accompanied 
Juptyer notebook ('Feature_Selection.ipynb').
"""

__author__ = "Jamie Scheper"
__version__ = "v1.0"
__date__ = "25-01-2024"
__email__ = "j.s.c.scheper@st.hanze.nl"

from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import scikitplot as skplt

from scipy.interpolate import interp1d
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFECV, SelectKBest, mutual_info_classif
from sklearn.metrics import (auc, classification_report, confusion_matrix,
                             f1_score, precision_score, recall_score, 
                             roc_auc_score, roc_curve)
from sklearn.model_selection import (LeaveOneOut, StratifiedKFold,
                                     cross_val_predict, cross_validate)
from sklearn.preprocessing import LabelBinarizer, LabelEncoder, label_binarize


def train(data, model, *, cls='cluster', scoring=None, plot=False):
    """
    Trains a specified model using provided training data and makes predictions on test data. 
    This function allows customization of scoring metrics and offers an option to plot a ROC curve 
    for multi-class classification problems. For binary classification, it returns false and true 
    positive rates if the plot parameter is enabled.

    Parameters:
        data (a tuple): A tuple containing training and test datasets 
                       (format: (x_train, x_test, y_train, y_test))
        model (a classifier): Model expected to be trained and used for classification.
        cls (str, optional): The name of the class variable. Default is 'cluster'.
        scoring (None, optional): List of scoring metrics to be used.
                               Default is None, then uses: recall_score, precision_score, f1_score
        plot (bool, optional): whether to plot ROC curves (standard is False).
    
    Returns:
        A tuple containing:
        - predictions (list): a list with predictions regarding test dataset.
        - y_test (numpy array): a numpy array with the real values of class variable.
        - scores (dict): dictionary of scores with all the given scoring metrics.
        - auc_res (None or dict): 
                        None for multi-class problems; 
                        dict of false positive rate, true positive rate, auc if binary class
    """
    # If scoring is not given, provide scoring metrics
    if scoring is None:
        scoring = [recall_score, precision_score, f1_score]

    # Unpack data
    x_train, x_test, y_train, y_test = data

    # If data consists of multiple class columns, get the right one
    if not isinstance(y_train, pd.Series):
        y_train, y_test = y_train[cls], y_test[cls]

    # How many classes?
    n_classes = len(np.unique(y_train))

    model.fit(x_train, y_train)
    prediction = model.predict(x_test)

    scores = {}
    # calculate the metrics, each in their own way
    for metric in scoring:
        if callable(metric):
            metric_name = metric.__name__
            if metric_name in ['recall_score', 'precision_score', 'f1_score']:
                scores[metric_name] = [metric(y_train, model.predict(x_train), average='weighted'), 
                metric(y_test, prediction, average='weighted')]
            # get log_loss via predict_proba
            elif metric_name in ['log_loss']:
                scores[metric_name] = [metric(y_train, model.predict_proba(x_train)),
                                       metric(y_test, model.predict_proba(x_test))]
            elif metric_name == 'roc_auc_score' and n_classes > 2:
                #  For multi-class, use one-vs-rest strategy
                scores[metric_name] = roc_auc_score(y_test, model.predict_proba(x_test), 
                multi_class='ovr', average='weighted')
                
    auc_res = None
    if plot:
        # binarize (even for a binary class label)
        label_binarizer = LabelBinarizer().fit(y_train)
        y_onehot_test = label_binarizer.transform(y_test) # one-vs-rest approach
        auc_res = calculate_roc(model, x_test, y_onehot_test, n_classes, label_binarizer.classes_)
        
    return prediction, y_test, scores, auc_res


def calculate_roc(model, x_test, y_test, n_classes, names):

    """
    Calculates and plots the ROC curve of a given model based on multi-class
    problems. Otherwise, return the false positive rate, true positive rate, 
    and area under the curve score.

    Parameters:
        model (a classifier): model to be used in calculating and plotting ROC
        x_test (a pd dataframe): counts used to predict (test set)
        y_test (a numpy array): true labels of class variable (test set)
        n_classes (int): number of classes
        names (list of str): a list of string with the nonbinarized names of class values

    Returns (only for binary classes):
        - fpr - false positive rate
        - tpr - true positive rate
        - auc_score - area under the curve score
    """
    # prediction based on probabilities
    prediction = model.predict_proba(x_test)

    # only for multi-class problems
    if n_classes > 2:
        _, ax = plt.subplots()
        colors = sns.color_palette(n_colors=n_classes)

        mean_fpr = np.linspace(0, 1, 100) # evenly spaced points to calculate mean tpr
        tprs = []
        n = len(y_test)
        weights = [np.sum(y_test[:, i])/n for i in range(n_classes)] # class weights

        # go per class in a one-vs-rest approach
        for class_number in range(n_classes):
            fpr, tpr, _ = roc_curve(y_test[:, class_number], prediction[:, class_number])
            auc_score = auc(fpr, tpr)
            ax.plot(fpr, tpr, color=colors[class_number], 
            label=f'ROC curve for {names[class_number]} (AUC = {round(auc_score, 2)})')
            
            intp_tpr = np.interp(mean_fpr, fpr, tpr)
            intp_tpr[0] = 0.0 # begin at 0
            tprs.append(intp_tpr*weights[class_number]) # weigh based on class proportion

        weighted_tpr = np.sum(tprs, axis=0) # weighed average tpr
        weighted_auc = auc(mean_fpr, weighted_tpr) # average auc score
        ax.plot(mean_fpr, weighted_tpr, color='black', linestyle='dashed',
            label=f'ROC curve weighted mean (AUC = {round(weighted_auc, 2)})')
        plt.plot([0, 1], [0, 1], color='red', linestyle='dashed') # diagonal line
        plt.xlabel("False positive rate")
        plt.ylabel("True positive rate")
        plt.title(f"ROC curve for {type(model).__name__}")
        plt.legend()

        return fpr, tpr, auc_score
    else:
        fpr, tpr, _ = roc_curve(y_test, prediction[:, 1])
        auc_score = auc(fpr, tpr)
    
        return fpr, tpr, auc_score


def visualize_features(data, *, method="count", feature_name="cluster", 
                       hue=None, y=None, stat="percent"):
    """
    Visualize features by either count, bar, or box pot given a dataset.
    Parameters:
        data (pd dataframe) - clinical metadata dataset
        method (str, optional) - which method to plot with
            (can either be 'count' (standard), 'bar', 'boxplot').
        feature_name (str, optional) - the feature to be on x axis (cluster is standard).
        hue (str, optional) - variable to be used as the color encoding (standard is None).
        y (str, optional) - the feature to be on y axis. 
            Only availble for method='bar' (standard is None).
        stat (str, optional) - what kind of stat to use. 
            Only available for method='count' (standard is percent)

    Returns: -
    """
    try:
        if method == "count":
            sns.countplot(data=data, x=feature_name, hue=hue, stat=stat)

            # Add labels
            plt.title(f"Count of feature {feature_name} {f'on {hue}' if hue is not None else ''}")
            plt.xlabel(f'{feature_name}')
            plt.ylabel('Percentage')
        elif method == "bar":
            sns.barplot(x=feature_name, y=y, hue=hue, data=data)
            # Add labels
            plt.title(f'Bar plot of feature {feature_name} based on {hue}')

        elif method == "boxplot":
            sns.boxplot(x=feature_name, y=y, data=data, hue=hue)

            plt.title(f'Boxplot of {feature_name}')

    except ValueError as e:
        print(f"Error: {e} (use either method='count', method='boxplot', or method='bar')")


def visualize_results(data, *, confusion_mat=True, summary_table=True, plot=True):
    """
    A function to visualize test and cross-validation results. Function recognizes
    either and summarizes the results of all evaluated models at once. 
    Can make a confusion matrix and summary plot depicting scoring metrics.
    Plotting these is optional.

    Parameters:
        data (dictionary of lists) - dataframe with evaluated models' scoring metrics
        confusion_mat (optional, bool) - whether to calculate a confusion matrix
        summary_table (optional, bool) - whether to construct a table with scoring metrics
        plot (optional, bool) - whether to plot confusion matrix and/or summary table

    Returns: -
    """
    # Confusion matrix
    if confusion_mat:
        col_names = ["Classifier", "Prediction", "Real"]
        # flatten the df; only keep predicted and real valuations
        flattened_df = pd.DataFrame([[name, params[0].tolist(), params[1].tolist()] 
                                   for name, params in data.items()], 
                                   # explore means inner list to rows in pd dataframe
                                  columns=col_names).explode(["Prediction", "Real"])
        if plot:
            plot_confusion_matrix(flattened_df)
    
    if summary_table:
        table = defaultdict(list) # creates key entry instead of KeyError; produces empty list
        
        # Recalibrate our scoring table
        for name, params in data.items():
            scores = params[2]
            table['Classifier'].append(name)
            for metric, values in scores.items():
                # remove redundancies
                if metric not in ["fit_time", "score_time", "roc_auc_score"] \
                    and not np.isnan(values).any():
                    metric_name = metric.replace("score", "").replace("test", "mean").\
                        replace("_", " ").capitalize()
                    mean_score = np.mean(values) if len(values) > 2 else values[1]
                    table[metric_name].append(mean_score)

        if plot and len(table) > 1:
            plot_summary_table(table)


def plot_confusion_matrix(table):
    """
    Plots a confusion matrix from a filtered table. Also prints a 
    classification report, indicating several scoring metrics.

    Parameters:
        table (pd dataframe) - a pandas dataframe with:
            - classifier name
            - predicted values
            - real values
    
    Returns: -
    """
    n_classifiers = table['Classifier'].unique()
    fig, ax = plt.subplots(2, 3)
    fig.suptitle("Confusion matrices for all used classifiers")

    for index, classifier in enumerate(n_classifiers):
        # filter df
        labels = np.unique(table["Real"])
        # checks whether classifier is the right one, then drop from original table
        filtered_df = table[table['Classifier'] == classifier].drop('Classifier', axis=1)
        cm = confusion_matrix(filtered_df["Real"], filtered_df["Prediction"])

        # In addition, print a classification report
        print(f"The resulting classification report for {classifier} \n", 
              classification_report(filtered_df["Real"], filtered_df["Prediction"],))

        # plot the confusion matrix without color bar
        ax.flat[index].set_title(f"{classifier}")
        # ax.flat[index] on the right subplot
        sns.heatmap(cm, ax=ax.flat[index], annot=True, cbar=False,
                    xticklabels=labels, yticklabels=labels)


def plot_summary_table(table):
    """
    A small function to plot a table from a pandas dataframe.

    Parameters:
        table (pd dataframe) - a pandas dataframe with:
            - classifier name
            - scoring metrics
    
    Returns: -
    """
    df_table = pd.DataFrame(table).set_index("Classifier").round(decimals=2)

    _, ax = plt.subplots()
    ax.axis('off') # no axis lines (because table)

    pd.plotting.table(ax, df_table, loc='center', cellLoc='center')


def cross_val(data, model, *, cls='cluster', cv=StratifiedKFold(n_splits=5), 
              scoring='accuracy', plot=False):
    """
    A function to cross-validate with very specific settings.
    Only cross-validates on training data.
    The function allows for a plot that visualizes 
    the ROC curve for all class values in a one-vs-rest approach.

    Parameters:
        data (a tuple): A tuple containing training and test datasets 
                       (format: (x_train, x_test, y_train, y_test))
        model (a classifier): Model expected to be trained and used for classification.
        cls (str, optional): The name of the class variable. Default is 'cluster'.
        cv (a cross-validation method): cross-validation method to use. 
            Standard is StratifiedKFold(n_splits=5).
        scoring (str, optional): List of scoring metrics to be used.
                               Default is 'accuracy'
        plot (bool, optional): whether to plot ROC curves in a one-vs-rest approach

    Returns:
        A tuple containing:
        - predictions (list): a list with predictions regarding test dataset.
        - y_test (numpy array): a numpy array with the real values of class variable.
        - scores (dict): dictionary of scores with all the given scoring metrics.
        - auc_res (None): None since the function visualize_results expects four inputs.
    """
    if plot:
        cv_plot_per_class(data, model, cls=cls, cv=cv)
    
    # we only need training data
    x_train, _, y_train, _ = data
    # Throw away other than needed class
    if not isinstance(y_train, pd.Series):
        y_train = y_train[cls]

    scores = cross_validate(model, x_train, y_train, cv=cv, scoring=scoring)
    prediction = cross_val_predict(model, x_train, y_train, cv=cv)
    
    # One since the function visualize_results expects four inputs!
    return prediction, y_train, scores, None
    

def cv_plot_per_class(data, model, *, cls='cluster', cv=StratifiedKFold(n_splits=5)):
    """
    Function plots an ROC per class valuation in regarding the class variable. The ROC
    curves are displayed in the same plot, with the confidence interval. ROC curves are
    calculated in a one-vs-rest approach and only for the training data! Training data
    is split into training and test (cross-validation).

    Parameters:
        data data (a tuple): A tuple containing training and test datasets 
                       (format: (x_train, x_test, y_train, y_test))
        model (a classifier): Model expected to be trained and used for classification.
        cls (str, optional): The name of the class variable. Default is 'cluster'.
        cv (a cross-validation method): cross-validation method to use. 
            Standard is StratifiedKFold(n_splits=5).
    
    Returns: -
    """
    x_train, _, y_train, _ = data
    # Throw away other than needed class
    if not isinstance(y_train, pd.Series):
        y_train = y_train[cls]
    
    # Split data
    x, y = x_train, y_train

    # Binarize for multi-class
    unique_classes = np.unique(y)
    y_bin = label_binarize(y_train, classes=unique_classes)
    n_classes = len(unique_classes)

    _, ax = plt.subplots()
    # diagonal line
    plt.plot([0, 1], [0, 1], color='black', linestyle='dashed')

    # prepare a linear space for plotting and calculating mean tpr
    mean_fpr = np.linspace(0, 1, 100)
    
    # loop per class
    for n in range(n_classes):
        tprs, aucs = [], []
        
        # split into training and test
        for training, test in cv.split(x, y):
            x_train, x_test = x.iloc[training], x.iloc[test]
            # more than two classes, use [:, n]
            model_fit = model.fit(x_train, y_bin[training][:, n] if n_classes > 2 else y_bin[training])
            prediction = model_fit.predict_proba(x_test) # get probabilities for roc
            fpr, tpr, _ = roc_curve(y_bin[test][:, n] if n_classes >2 else y_bin[test], prediction[:, 1])

            # calc the tpr, fpr
            interp = interp1d(fpr, tpr)
            result = interp(mean_fpr)
            result[0] = 0 # start line at 0
    
            tprs.append(result)
            aucs.append(auc(fpr, tpr))

        # Some stats! (mean and standard deviations)
        mean_tpr, sd_tpr = np.mean(tprs, axis=0), np.std(tprs, axis=0)
        mean_auc, sd_auc = auc(mean_fpr, mean_tpr), np.std(aucs, axis=0)

        ax.step(mean_fpr, mean_tpr, label=f"""Class {unique_classes[n]} 
            (AUC = {mean_auc:.2f} +/- {sd_auc:.2f})""")
        # fill in the confidence interval
        ax.fill_between(mean_fpr, mean_tpr-sd_tpr, mean_tpr+sd_tpr, color="grey", alpha=0.15)
        
    ax.legend(loc="lower right")
    plt.title("ROC Curves")
    plt.xlabel("False positive rate")
    plt.ylabel("True positive rate")
    plt.show()
    

def train_and_cross(data, model, *, cls='cluster', cv=LeaveOneOut(), scoring=['accuracy']):
    """
    Train and cross-validate the training data via an independent route to
    validate findings (package: skplt).

    Parameters:
        data data (a tuple): A tuple containing training and test datasets 
                       (format: (x_train, x_test, y_train, y_test))
        model (a classifier): Model expected to be trained and used for classification.
        cls (str, optional): The name of the class variable. Default is 'cluster'.
        cv (a cross-validation method, optional): cross-validation method to use. 
            Standard is LeaveOneOut().
        scoring (str, optional): scoring metric to use. Standard is accuracy.

    Returns: -
    """
    x_train, _, y_train, _ = data
    # Throw away other than needed class
    if not isinstance(y_train, pd.Series):
        y_train = y_train[cls]
    
    # Split data
    x, y = x_train, y_train
    
    # Plot for each scoring metric
    for score in scoring:
        skplt.estimators.plot_learning_curve(
        model, x, y, cv=cv, shuffle=False, scoring=score,
        n_jobs=-1, title_fontsize="large",
        text_fontsize="large", title=f"{model} Learning Curve based on {score}")


def select_features_cv(data, model, class_label, *, 
                       cv=StratifiedKFold(n_splits=5), scoring='accuracy',
                       n_features=None, njobs=-1):
    """
    Select the best feature sets based on a scale of parameters with RFECV.

    Returns scoring on all possibilities - selection of best-fitted feature set follows 
    after training, testing and cross-validating said feature sets.

    Parameters:
        data data (a tuple): A tuple containing training and test datasets 
                       (format: (x_train, x_test, y_train, y_test))
        model: machine-learning model to use
        class_label: feature to predict
        [Optional]
        cv: cross-validation method to use (default: StratifiedKFold(n_split=5)
        scoring: scoring metric(s) to use (default: accuracy)
        n_features: how many features does a dataframe need
        njobs: how many cpus to use (-1 = all available)

    Returns a dataframe with given scores for all possible feature combinations
    """
    # Prepare dataset
    x_train, _, y_train, _ = data
    # Throw away other than needed class
    y_train = y_train[class_label]
    
    # Split data
    x, y = x_train, y_train

    # specifcy the RFECV
    rfecv = RFECV(
        model,
        cv=cv,
        scoring=scoring,
        n_jobs=-1
    )

    rfecv.fit(x, y)
        
    print(f"Optimal number of features {rfecv.n_features_}")

    mean_scores = rfecv.cv_results_["mean_test_score"]
    std_scores = rfecv.cv_results_["std_test_score"]
        
    _, ax = plt.subplots()
    ax.errorbar(range(len(mean_scores)), mean_scores, yerr=std_scores, alpha=0.6)
    plt.ylabel(f"{scoring}")
    plt.xlabel("Feature")
    plt.title(f"""Recursive feature selection with cross validation on {type(model).__name__}""")
    plt.show()

    reduced_dataframe = rfecv.transform(x)

    return reduced_dataframe, rfecv.support_, rfecv.ranking_

def assess_feature_importances(data, models, *, selector=RandomForestClassifier(), cls='cluster', 
                                                plot=False, batch_size=None, stop=0):
    """
    A function that assesses feature importance. Based on the selector, which can be either:
        - random forest: uses feature_importance_ attribute (Gini impurity)
        - logistic regression with L1 or L2: coef_ (coefficients)
        - model independent (when parameter selector is None): mutual information!
    We perform feature selection based on those scores via batches. Batch sizes can be specified
    by user and where to stop adding features. Different sized feature sets are trained with and
    cross-validated.

    Parameters:
        data (a tuple): A tuple containing training and test datasets 
                       (format: (x_train, x_test, y_train, y_test))
        models (list of classifiers) - classifiers to train the reduced features on
        selector (classifier, optional) - which model to use as feature selector 
                (either RandomForestClassifier (standard), LogisticRegression, or None)
        cls (str, optional) - The name of the class variable. Default is 'cluster'.
        plot (bool, optional) - whether to plot results
        batch_size (int, optional) - how big batches need to be (positive int)
        stop (int, optional) - at what size to stop adding features 
            (0 means all features and is standard).

    Returns:
        if batch_size is not None and above 0: 
            names of features (ranked),
            test scores and cross-validated results
        else:
            names of features (ranked)  
    """
    # Prepare dataset
    x_train, x_test, y_train, y_test = data

    # Throw away other than needed class
    y_train_tmp = y_train[cls]

    # binarize if selector is ElasticNet
    if type(selector).__name__ ==  'ElasticNet':
        label_encoder = LabelEncoder()
        y_train_tmp = label_encoder.fit_transform(y_train_tmp)

    # Split data
    x, y = x_train, y_train_tmp
    
    if selector is not None:
        selector.fit(x, y)
    
    if type(selector).__name__ == "RandomForestClassifier":
        # Gini impurity
        important_features = selector.feature_importances_
        df = pd.DataFrame({"features": x.columns, "importance": important_features})
        df = df.sort_values('importance', ascending=False, kind='stable') 
        
    elif type(selector).__name__ in ['LogisticRegression', 'ElasticNet']:
        coefficients = selector.coef_.flatten() # create a 1D array
        df = pd.DataFrame({"features": x.columns, "importance": coefficients, "coeff_abs": np.abs(coefficients)})
        # remove coefficients equal to 0
        df = df[df['coeff_abs'] != 0].sort_values('coeff_abs', ascending=False, kind='stable')
    
    if plot:
        # plot result of ranking features based in importance
        names = plot_top_features(df, selector, stop)

    # train and cross-validate with different sized feature dataframes
    if batch_size is not None and not batch_size <= 0:
        if selector is None:
            # perform mutual information if selector is None
            selector = SelectKBest(mutual_info_classif, k=stop).fit(x, y)
            feature_scores = selector.scores_
            names = x.columns

            df = pd.DataFrame({'features': names, 'Score': feature_scores})
            df = df.sort_values(by='Score', ascending=False)
        
        # execute this for every model at once
        result = [execute_model(x, x_test, y_train, y_test, df, model, stop, batch_size, cls) for name, model in models.items()]

        return names, result

    return names


def execute_model(x_train, x_test, y_train, y_test, df, model, stop, batch_size, cls):
    """
    Execute a model on different sized feature dataframes, as specified in the assess_feature_importances
    function. Loops through the gene/features by using a batch size and stops adding at a specified point
    via parameter stop.

    Parameters:
        x_train - gene features (train set)
        x_test - class labels (train set)
        y_train - gene features (test set)
        y_test - class labels (test set)
        df - the ranked features based on gini impurity, coefficients or mutual information
        model - which classifier to train and cross-validate on
        stop - at what size do we stop
        batch_size - how big are batch sizes (e.g., 5 = add 5 featurs per loop)
        cls - class variable

    Returns:
        res_train - dict of test results per batch 
        res_cross - dict of cross-validated results per batch 
    """
    res_train, res_cross = defaultdict(list), defaultdict(list)

    for index in range(0, stop, batch_size):
        # select genes based on batch size from the ranked dataframe
        selected_genes = df['features'].iloc[0:index+batch_size].tolist()
        # reduce the counts datasets for training and testing
        tmp_df_train = x_train[[column for column in x_train.columns if column in selected_genes]]
        tmp_df_test = x_test[[column for column in x_test.columns if column in selected_genes]]
        
        new_data = tmp_df_train, tmp_df_test, y_train, y_test

        # Add to a defaultdict and only retain scores (through '[2]')
        res_train[f"train_{index + batch_size}"] = train(new_data, model, cls=cls, scoring=[recall_score])[2]
        res_cross[f"cross_{index + batch_size}"] = cross_val(new_data, model, cls=cls)[2]
    
    return res_train, res_cross


def plot_top_features(dataset, selector, top_n=50):
    """
    Function that visualizes the ranking established based on a selector.

    Parameters:
        dataset - a ranked dataframe from the function assess_feature_importances
        selector - used selector that made the ranking
        top_n - how many features do we include in the visualization

    Returns: a list with the top_n features
    """
    selector_name = type(selector).__name__
    top_features = dataset.head(top_n) # get top_n features
    top_features.plot.barh(x='features', y='importance', legend=None)

    plt.title(f"Top {top_n} features based on {selector_name}")
    plt.axvline(x=0, color="grey")
    plt.xlabel(f"{'Coefficient value' if selector_name != 'RandomForestClassifier' else 'Gini importance'}")
    plt.ylabel("Features")
    plt.yticks(fontsize=8)
    plt.savefig("ridge.png")
    plt.show()

    return top_features


def plot_table(data, *, batch_size=5):
    """
    Plot a table based on the results from the "assesses_feature_importance" function.
    Plots a different table, depending if given training or cross-validation data.
    Has the ability to plot a line plot of training and test scores.

    Parameters:
        data (tuple) - a tuple consisting of two dictonarnies with training and 
                       cross-validation scores
        batch_size (int, optional) - an int with the used batch size

    Returns: -
    """

    # Transpose table
    data = pd.DataFrame(data).T
    
    for column, values in data.items():
        if column not in ['fit_time', 'score_time']:
            # pd.Series to pd.DataFrame and round values
            tmp_df = pd.DataFrame(values.tolist()).round(2)
            
            # remove prefix from training and cross labels
            tmp_df.index = [label.removeprefix('train_').removeprefix('cross_') for label in values.index]

            # cross-validation results
            if len(tmp_df.columns) > 2:
                fig, ax = plt.subplots(1, figsize=(15, 8))
                # get the same amount of columns as there have been cross-validation folds
                tmp_df.columns = [f'Fold {int(col) + 1}' for col in tmp_df.columns]
                tmp_df_ = tmp_df.T
                # calculate the mean score of all folds
                tmp_df['Mean score'] = [round(np.mean(value), 2) for value in values]
                ax.set_xlabel('Gene set size')
                # rotate x-axis labels to 45 degrees
                ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation=45, ha='right')
                ax.set_ylabel(f'{column}')
                pd.plotting.table(ax, tmp_df, loc='center', cellLoc='center')
                ax.axis('off') # no outline
            # training results
            else:
                fig, ax = plt.subplots(1, 2, figsize=(15, 8))
                tmp_df.columns = ["Train", "Test"]
                ax[0].plot(tmp_df.index, tmp_df["Train"], marker='o')
                ax[0].plot(tmp_df.index, tmp_df["Test"], marker='x')
                ax[0].set_xlabel('Gene set size')
                # rotate x-axis labels to 45 degrees
                ax[0].set_xticks(ax[0].get_xticks(), ax[0].get_xticklabels(), rotation=45, ha='right')
                ax[0].set_ylabel(f'{column}')
                ax[0].legend(["Train", "Test"])
    
                ax[1].axis('off') # no outline
                fig.suptitle(f"Scores of metric {column}")
                pd.plotting.table(ax[1], tmp_df, loc='center', cellLoc='center')
            
            plt.show()


def check_name(model_name):
    """
    A simple helper function used in the function roc_in_one to change the name
    of OneVsRestClassifier to SVC.

    Parameters:
        model_name (str) - a given model's name
    
    Returns: the model_name, adjusted when model name is OneVsRestClassifier
    """
    if model_name == 'OneVsRestClassifier':
        return 'SVC'

    return model_name


def roc_in_one(data):
    """
    Plots multiple ROC curves, extracted from the train function into one plot.

    Parameters:
        data (pd dataframe) - a dict with lists regarding predicting information of testing models

    Returns: -
    """

    for model, model_data in data.items():
        # retrieve the data (scoring metrics)
        points = model_data[3]
        # use the check_name function to adjust the name of SVC
        plt.plot(points[0], points[1], label=f"{check_name(model)} (AUC = {points[2]:.2f})")

    plt.title("ROC curves for all optimized classifiers")
    plt.xlabel("False positive rate")
    plt.ylabel("True positive rate")
    plt.plot([0, 1], [0, 1], color='black', linestyle='dashed')
    plt.legend()
    plt.grid(True)
    plt.show()
