__version__ = "0.4.0"
import os
import sys
import xgboost as xgb
import numpy as np
import pandas as pd
import argparse
from sklearn import preprocessing
from sklearn.model_selection import (
    StratifiedShuffleSplit,
    StratifiedKFold,
    ShuffleSplit,
    KFold,
    GridSearchCV,
)
import time
from sklearn.metrics import (
    roc_curve,
    auc,
    average_precision_score,
    f1_score,
    roc_auc_score,
    r2_score,
    mean_squared_error,
    median_absolute_error,
    explained_variance_score,
    accuracy_score,
    precision_score,
    recall_score,
)
import matplotlib.pyplot as plt
from math import sqrt
import json
from tqdm import tqdm
import warnings
import seaborn as sns

warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)

# Import SHAP for model interpretability
try:
    import shap
    SHAP_AVAILABLE = True
except ImportError:
    SHAP_AVAILABLE = False
    print("Warning: SHAP library not found. SHAP values will not be computed.")
    print("To install: pip install shap")


class BlankLinesHelpFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        return super()._split_lines(text, width) + [""]


def load_json_file(file_path):
    """Safely load a JSON file with proper error handling."""
    try:
        with open(file_path, "r") as f:
            return json.load(f)
    except FileNotFoundError:
        print(f"Error: File not found: {file_path}")
        return None
    except json.JSONDecodeError:
        print(f"Error: Invalid JSON format in file: {file_path}")
        return None
    except Exception as e:
        print(f"Error loading file {file_path}: {str(e)}")
        return None


def create_default_param_grid(is_tuning=False):
    """Create default parameter grid if JSON file not found or invalid."""
    if is_tuning:
        # Larger grid for tuning runs
        return {
            "max_depth": [3, 5, 7, 9, 11],
            "learning_rate": [0.01, 0.05, 0.1, 0.2, 0.3],
            "n_estimators": [50, 100, 200, 300],
            "subsample": [0.6, 0.8, 1.0],
            "colsample_bytree": [0.6, 0.8, 1.0]
        }
    else:
        # Smaller grid for regular runs
        return {
            "max_depth": [3, 5, 7],
            "learning_rate": [0.05, 0.1, 0.2],
            "n_estimators": [100, 200],
            "subsample": [0.8, 1.0],
            "colsample_bytree": [0.8, 1.0]
        }


def compute_shap_values(model, X_test, feat_names):
    """Compute and visualize SHAP values with enhanced styling for model interpretation."""
    if not SHAP_AVAILABLE:
        return None
    
    try:
        # Apply seaborn styling
        sns.set_style("whitegrid")
        sns.set_context("talk")
        
        # Create SHAP explainer
        explainer = shap.Explainer(model)
        
        # Calculate SHAP values
        shap_values = explainer(X_test)
        
        # Summary plot with enhanced styling
        plt.figure(figsize=(12, 10))
        shap.summary_plot(shap_values, X_test, feature_names=feat_names, show=False, plot_size=(12, 10))
        plt.title("SHAP Feature Importance", fontsize=16, pad=15)
        plt.tight_layout()
        
        # Waterfall plot with enhanced styling
        plt.figure(figsize=(12, 8))
        shap.waterfall_plot(shap_values[0], show=False, max_display=15)
        plt.title("SHAP Waterfall Plot (First Instance)", fontsize=16, pad=15)
        plt.tight_layout()
        
        # Bar plot with enhanced styling
        plt.figure(figsize=(12, 8))
        shap.plots.bar(shap_values, show=False, max_display=15)
        plt.title("Mean SHAP Value Impact", fontsize=16, pad=15)
        plt.tight_layout()
        
        # Create dependence plots for top features
        # Get mean absolute shap values for ranking features
        mean_abs_shap = np.abs(shap_values.values).mean(0)
        feature_importance = pd.DataFrame(list(zip(feat_names, mean_abs_shap)), 
                                         columns=['feature','importance'])
        feature_importance = feature_importance.sort_values('importance', ascending=False)
        
        # Plot dependence plots for top 3 features
        for i in range(min(3, len(feature_importance))):
            plt.figure(figsize=(10, 6))
            feature_idx = list(feat_names).index(feature_importance['feature'].iloc[i])
            shap.dependence_plot(feature_idx, shap_values.values, X_test, 
                               feature_names=feat_names, show=False)
            plt.title(f"SHAP Dependence Plot - {feature_importance['feature'].iloc[i]}", 
                     fontsize=16, pad=15)
            plt.tight_layout()
            
        return shap_values
        
    except Exception as e:
        print(f"Error computing SHAP values: {str(e)}")
        return None


def main(args):
    descriptor_message = (
        "Use the --help flag to get more information about XGBeast.\n\nYou are using XGBeast "
        + __version__
        + "!\n\n\n"
    )
    #  COMMAND LINE ARGUMENTS
    parser = argparse.ArgumentParser(
        prog="XGBeast",
        formatter_class=BlankLinesHelpFormatter,
        epilog=descriptor_message,
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s {version}".format(version=__version__),
    )
    parser.add_argument(
        "-name",
        "--output_name",
        help="name to be included in the model output folder name (optional)\n",
        default="",
        type=str,
    )
    parser.add_argument(
        "-path",
        "--path_to_input_data",
        help='path to the folder containing the required "input_data" directory (mandatory)',
        type=str,
    )
    parser.add_argument(
        "-x",
        "--mode",
        help='Choose classification ("class") or regression ("reg")',
        choices=["reg", "class"],
    )
    parser.add_argument(
        "-permute",
        "--permute_y",
        help='Run the model in "PERMUTED mode". This is used a sanity check / statistical control. Will randomly shuffle the y / target at each iteration.',
        dest="permute_mode",
        action="store_true",
    )
    parser.add_argument(
        "-t",
        "--threads",
        help="number of threads to use (default = all available threads / -1)",
        type=int,
        default=-1,
    )
    parser.add_argument(
        "-n",
        "--stability_runs",
        help="number of shuffles / stability runs (default = 100)",
        type=int,
        default=100,
    )
    parser.add_argument(
        "-test",
        "--test_fraction",
        help="fraction of dataset to use as test (default = 0.2)",
        type=float,
        default=0.2,
    )
    parser.add_argument(
        "-cv",
        "--cv_folds",
        help="number of Cross-Validation folds (default = 5)",
        type=int,
        default=5,
    )
    parser.add_argument(
        "-stop",
        "--early_stopping_rounds",
        help="value of early stopping rounds (default = 50)",
        type=int,
        default=50,
    )
    parser.add_argument(
        "-rand_seed",
        "--random_seed_value",
        help="Random seed value. Must be an integer. (default = 512)",
        type=int,
        default=512,
    )
    parser.add_argument(
        "-param",
        "--parameter_grid_file",
        help="Path to file containing the XGBoost parameter grid dictionary"
        '(default = dictionary found in the "parameter_grid.json" file)',
        type=str,
        default="param_grid.json",
    )
    parser.add_argument(
        "-scale",
        "--scale_data",
        help="Scale input data (predictors).\n"
        'Options: "none" (no scaling; default); "scale" (zero-mean, unit-variance scaling); \n\n"'
        'robust-scale" (scales so that median becomes 0 and IQR becomes 1; may be prefferable for data with outliers); '
        '"range" (range scale, all features scaled to range from 0 to 1; may be prefferable for sparse data).',
        type=str,
        default="none",
    )
    parser.add_argument(
        "-no-rv",
        "--no-rand-vars",
        dest="add_random_variables",
        help="2 random variables are (by default) added as predictors, as a sanity check and to serve as a benchmark for the relative importance of other features. "
        'These variables are randomized at each model iteration. The "-no-rv" prevents these random variables from being added.',
        action="store_false",
    )
    parser.add_argument(
        "-no-dt",
        "--do_not_add_date_time",
        help='Date and Time are (by default)s added as a suffix to the output folder name. The "-no-dt" flag allows this to be switched off.',
        dest="add_date_time",
        action="store_false",
    )
    parser.add_argument(
        "-tune",
        "--tune-parameters",
        dest="tune_parameters",
        help="This will run the model on the input data using a low number of shuffles (n=5) but on a large parameter grid. "
        'Will save best parameters of each iteration to a file: "best_parameters.txt".'
        'The large grid is stored in "grid_tuning.json".',
        action="store_true",
    )
    parser.add_argument(
        "-s",
        "--silent",
        help="Limits verbosity. Will hide additional metrics and info at each iteration.",
        dest="verbosity",
        action="store_false",
    )
    parser.add_argument(
        "-shap",
        "--compute-shap",
        dest="compute_shap",
        help="Compute SHAP values for model interpretability. Requires the shap library.",
        action="store_true",
    )
    parser.add_argument(
        "-no-conda",
        "--no-conda-env",
        dest="save_conda_env",
        help="Disable saving conda environment information. Useful if conda is not installed.",
        action="store_false",
    )

    parser.set_defaults(add_random_variables=True)
    parser.set_defaults(add_date_time=True)
    parser.set_defaults(tune_parameters=False)
    parser.set_defaults(verbosity=True)
    parser.set_defaults(permute_mode=False)
    parser.set_defaults(compute_shap=False)
    parser.set_defaults(save_conda_env=True)
    args = parser.parse_args(args[1:])

    if args.path_to_input_data is None:
        print(
            '\n\nYou must specify a <<< PATH >>> to the folder with the input data!\n\nThis folder most countain a subfolder "input_data" containing: \n    "X_data.txt" containing the features / predictors\n    "y_reg.txt" file for regression models OR "y_binary.txt for classification models (containing the predicted outcome) \n    "feat_ids.txt" file, containing the feature names.\n\nAll files must be in tab-delimited format.\n\n\nUse the "-path" flag!\n\n'
        )
        print(descriptor_message)
        return

    if args.mode is None:
        print(
            '\n\nYou must specify a <<< MODE >>> for the model!\n\nThis can be either CLASSIFICATION ("class") or REGRESSION ("reg")\n\nUse the "-x" flag follower by one of the 2 options!\n\n'
        )
        print(descriptor_message)
        return

    # classification OR regression
    if args.mode == "class":
        model_type = "XGB_class"
    else:
        model_type = "XGB_reg"

    name = args.output_name

    if args.add_date_time:
        # append date / time to output folder name (default, can be switched off with -dt / --datetime False)
        from datetime import datetime

        now = datetime.now()
        dt_string = now.strftime("%Y_%m_%d__%H-%M-%S")
        name = name + "_" + dt_string

    nthreads = args.threads

    # Set environment variables to control threading behavior
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"

    # random seed
    rand_seed = args.random_seed_value

    # model configuration
    N_stab = args.stability_runs
    test_size = args.test_fraction
    train_size = 1 - test_size
    CV_fold = args.cv_folds
    early_stopping_rounds = args.early_stopping_rounds
    permute_mode = args.permute_mode
    scale_data = args.scale_data
    tune_parameters = args.tune_parameters
    compute_shap = args.compute_shap and SHAP_AVAILABLE

    if permute_mode:
        name = name + "_PERMUTED"

    # Build paths correctly without changing the working directory
    base_path = args.path_to_input_data
    input_data_path = os.path.join(base_path, "input_data")
    
    # if it's a tuning run, do a few iterations with a large parameter grid
    # print best_parameters to both screen and file
    if tune_parameters:
        name = name + "_TUNING"
        N_stab = 5
        
        # Try loading the grid_tuning.json from script dir or data dir
        script_dir = os.path.dirname(os.path.abspath(__file__))
        param_file_options = [
            os.path.join(script_dir, "grid_tuning.json"),
            os.path.join(base_path, "grid_tuning.json")
        ]
        
        param_grid = None
        for file_path in param_file_options:
            param_grid = load_json_file(file_path)
            if param_grid:
                print(f"Loaded parameter grid from: {file_path}")
                break
                
        if not param_grid:
            print("Warning: Could not find grid_tuning.json, using default tuning grid.")
            param_grid = create_default_param_grid(is_tuning=True)
            
        best_parameters_list = []
    else:
        # Regular run - load parameter grid
        if os.path.isabs(args.parameter_grid_file):
            param_file = args.parameter_grid_file
        else:
            param_file = os.path.join(base_path, args.parameter_grid_file)
            
        param_grid = load_json_file(param_file)
        if not param_grid:
            print(f"Warning: Could not load parameter grid from {param_file}, using default grid.")
            param_grid = create_default_param_grid(is_tuning=False)
    
    # OUTPUT FOLDER NAME
    output_folder = os.path.join(base_path, "output_" + model_type + "_" + name)

    if args.verbosity:
        verbosity = 1
    else:
        verbosity = 0

    np.random.seed(rand_seed)
    start_time = time.time()

    # make output folder if it does not already exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder, exist_ok=True)
        print("\nDirectory " + output_folder + " Created!\n")
    else:
        print("\nDirectory " + output_folder + " already exists!\n")
    print("\nData loading! \n")

    # INPUT DATA
    # load X
    X_filepath = os.path.join(input_data_path, "X_data.txt")
    try:
        X = np.loadtxt(X_filepath, dtype=float, delimiter="\t")
    except Exception as e:
        print(f"Error loading X data: {str(e)}")
        return
        
    # load y
    if model_type == "XGB_class":
        y_filepath = os.path.join(input_data_path, "y_binary.txt")
        try:
            y_int = np.loadtxt(y_filepath, dtype=int, delimiter="\t")
            y = np.asarray(y_int, dtype=int)
        except Exception as e:
            print(f"Error loading y data: {str(e)}")
            return
    else:
        y_filepath = os.path.join(input_data_path, "y_reg.txt")
        try:
            y_float = np.loadtxt(y_filepath, dtype=float, delimiter="\t")
            y = np.asarray(y_float, dtype=float)
        except Exception as e:
            print(f"Error loading y data: {str(e)}")
            return

    # load feature names
    feat_names_filepath = os.path.join(input_data_path, "feat_ids.txt")
    try:
        feat_names_df = pd.read_csv(
            feat_names_filepath, sep="\t", dtype=str, header=None, usecols=[0]
        )
    except Exception as e:
        print(f"Error loading feature names: {str(e)}")
        return

    list_feat_importances, list_feat_names = [], []
    feat_names = np.asarray(feat_names_df.values, dtype=str)

    # add random variables names to list of features
    if args.add_random_variables:
        feat_names = np.append(feat_names, "random_variable1")
        feat_names = np.append(feat_names, "random_variable2")
    else:
        feat_names = feat_names.flatten()

    # scale data (optional, default is "no scaling')
    if scale_data != "none":
        if scale_data == "scale":
            X = preprocessing.scale(X)
        elif scale_data == "range":
            min_max_scaler = preprocessing.MinMaxScaler()
            X = min_max_scaler.fit_transform(X)
        elif scale_data == "robust-scale":
            X = preprocessing.robust_scale(X)
    # this prints out X data to file after optional scaling
    # pd.DataFrame(X).to_csv(os.path.join(output_folder, "X_data_as_used_in_model.csv"))

    # metrics
    if model_type == "XGB_class":
        (
            tprs, aucs, roc_auc_scores, f1s, avg_precision_scores, 
            balanced_accuracy_scores, mean_fp, accuracies, 
            precisions, recalls
        ) = ([] for _ in range(10))
        r = np.linspace(0, 1, 100)
        # stratified shuffle split for classification
        ShufSpl = StratifiedShuffleSplit(
            N_stab, test_size=test_size, train_size=train_size, random_state=rand_seed
        )
    else:
        r2s, expl_vars, rmses, maes = [], [], [], []
        ShufSpl = ShuffleSplit(
            N_stab, test_size=test_size, train_size=train_size, random_state=rand_seed
        )

    # initiate iterator
    i = 0

    # style for plots
    plt.style.use("ggplot")
    params = {
        "legend.fontsize": "medium",
        "figure.figsize": (10, 10),
        "axes.labelsize": "medium",
        "axes.titlesize": "medium",
        "xtick.labelsize": "medium",
        "ytick.labelsize": "medium",
    }
    plt.rcParams.update(params)
    plt.rcParams["font.family"] = "sans-serif"

    # Save the best model from all iterations for SHAP computation
    best_overall_model = None
    best_overall_score = -np.inf
    best_X_test = None

    # start iterations...
    for samples, test in tqdm(list(ShufSpl.split(X, y)), desc="Model Iterations"):
        i += 1
        print(
            "\nRunning iteration ---> "
            + " " * 3
            + str(i)
            + " out of "
            + str(N_stab)
            + " " * 3
            + "*" * 20
        )
        if verbosity:
            print(
                "\nTest set ratio = "
                + str(test_size)
                + "."
                + "\nNo. of CV = "
                + str(CV_fold)
                + ".\n"
            )

        if scale_data != "none":
            if scale_data == "scale":
                X = preprocessing.scale(X)
                print("\nZero-mean unit-variance scaling of features...")
            elif scale_data == "range":
                min_max_scaler = preprocessing.MinMaxScaler()
                print("\nRange scaling of features...")
                X = min_max_scaler.fit_transform(X)
            elif scale_data == "robust-scale":
                print("\nRobust scaling of features...")
                X = preprocessing.robust_scale(X)
            else:
                print("\nNo scaling...")

        # permute y if model is in permuted-mode
        if permute_mode:
            np.random.shuffle(y)
            print(
                "\nThis model is running in PERMUTED Mode!\nY is shuffled at each iteration.\n"
            )
        if verbosity:
            print(
                "\nDataset has", X.shape[0], "examples and", X.shape[1], "features.\n"
            )

        # split to train / test sets
        X_train, y_train = X[samples], y[samples]
        X_test, y_test = X[test], y[test]

        if tune_parameters:
            print("\n *** This is a parameter tuning run ...\n")
            print("Optimizing paramaters based on the grid:\n")
            print(param_grid)

        # add random variables
        if args.add_random_variables:
            # add random var1
            # random variables are randomized every iteration
            np.random.seed(i)
            random_column1_train = np.random.random(size=len(y_train)).reshape(
                len(y_train), 1
            )
            X_train = np.append(X_train, random_column1_train, axis=1)
            random_column1_test = np.random.random(size=len(y_test)).reshape(
                len(y_test), 1
            )
            X_test = np.append(X_test, random_column1_test, axis=1)
            # add random var2
            np.random.seed(i)
            random_column2_train = np.random.random(size=len(y_train)).reshape(
                len(y_train), 1
            )
            X_train = np.append(X_train, random_column2_train, axis=1)
            random_column2_test = np.random.random(size=len(y_test)).reshape(
                len(y_test), 1
            )
            X_test = np.append(X_test, random_column2_test, axis=1)
            if verbosity:
                print(
                    "\nData size after appending 2 random variables:",
                    X_train.shape[1],
                    "\n ",
                )
                print("Random variables are randomized at each test/train shuffle!")
            else:
                print("\n2 random variables were added!")

        if verbosity:
            print("\nTrain size = ", len(samples), "examples.\n")
            print("Test size = ", len(test), "examples\n")

        if model_type == "XGB_class":
            fit_params = {
                "verbose": False,
                "eval_set": [(X_train, y_train), (X_test, y_test)],
            }
            xgbclf = xgb.XGBClassifier(
                n_jobs=1,
                objective="binary:logistic",
                tree_method="auto",
                eval_metric="auc",
                early_stopping_rounds=early_stopping_rounds,
            )
            skf = StratifiedKFold(n_splits=CV_fold)
            scoring = "roc_auc"
        else:
            fit_params = {
                "verbose": False,
                "eval_set": [(X_train, y_train), (X_test, y_test)],
            }
            xgbclf = xgb.XGBRegressor(
                n_jobs=1,
                objective="reg:squarederror",
                tree_method="auto",
                eval_metric="rmse",
                early_stopping_rounds=early_stopping_rounds,
            )
            skf = KFold(n_splits=CV_fold)
            scoring = "r2"

        random_search = GridSearchCV(
            xgbclf,
            param_grid,
            n_jobs=nthreads,
            cv=skf,
            scoring=scoring,
            refit=True,
            verbose=verbosity,
        )
        
        # fit model
        try:
            random_search.fit(X_train, y_train, **fit_params)
        except Exception as e:
            print(f"Error fitting model: {str(e)}")
            continue

        # get best model on training set
        best_model = random_search.best_estimator_

        # test model
        y_pred = best_model.predict(X_test)

        if verbosity | tune_parameters:
            print("\nBest parameters were:\n")
            print(random_search.best_params_)
            print("\n")

        if tune_parameters:
            best_parameters_list.append(random_search.best_params_)

        # get feature importance
        feature_importance = best_model.feature_importances_
        list_feat_names.append(feat_names)
        list_feat_importances.append(feature_importance)

        if model_type == "XGB_class":
            fpr, tpr, _ = roc_curve(y_test, y_pred, pos_label=1)
            mean_fpr = np.linspace(0, 1, 100)
            auc_roc1 = auc(fpr, tpr)
            aucs.append(auc_roc1)
            tprs.append(np.interp(mean_fpr, fpr, tpr))
            tprs[-1][0] = 0.0
            roc_auc__score = roc_auc_score(y_test, y_pred)
            print("\n ROC-AUC-score: %.3f" % roc_auc__score)
            roc_auc_scores.append(roc_auc__score)

            # Additional classification metrics
            acc = accuracy_score(y_test, y_pred)
            accuracies.append(acc)
            print(" Accuracy:       %.3f" % acc)
            
            prec = precision_score(y_test, y_pred)
            precisions.append(prec)
            print(" Precision:      %.3f" % prec)
            
            rec = recall_score(y_test, y_pred)
            recalls.append(rec)
            print(" Recall:         %.3f" % rec)

            f1 = f1_score(y_test, y_pred)
            if verbosity:
                print(" F1-score:       %.3f" % f1)
            f1s.append(f1)

            avg_precision_score_val = average_precision_score(y_test, y_pred)
            if verbosity:
                print(" Avg precision:  %.3f" % avg_precision_score_val)
            avg_precision_scores.append(avg_precision_score_val)

            # Update best model for SHAP if this is the best so far
            if roc_auc__score > best_overall_score:
                best_overall_score = roc_auc__score
                best_overall_model = best_model
                best_X_test = X_test.copy()

            end_time = time.time()
            run_time = end_time - start_time

            print("\n\nrun_time: %.2f seconds" % run_time)
            print("\n" + "*" * 40 + "\n")
        else:
            r2 = r2_score(y_test, y_pred)
            print("\n R2-score:                %.3f" % r2)
            r2s.append(r2)

            expl_var = explained_variance_score(y_test, y_pred)
            print(" Explained Variance:      %.3f" % expl_var)
            expl_vars.append(expl_var)

            rmse = sqrt(mean_squared_error(y_test, y_pred))
            if verbosity:
                print(" Root Mean Squared Error: %.3f" % rmse)
            rmses.append(rmse)

            mae = median_absolute_error(y_test, y_pred)
            if verbosity:
                print(" Median Absolute Error:   %.3f" % mae)
            maes.append(mae)

            # Update best model for SHAP if this is the best so far
            if r2 > best_overall_score:
                best_overall_score = r2
                best_overall_model = best_model
                best_X_test = X_test.copy()

            end_time = time.time()
            run_time = end_time - start_time
            if verbosity:
                print("\n\nRuntime: %.2f seconds" % run_time)
            print("\n")
            print("\n" + "*" * 40 + "\n")

    # feature importance
    feature_importance = np.mean(np.array(list_feat_importances), axis=0)
    feature_importance = np.round(
        100.0 * (feature_importance / feature_importance.max()), 2
    )
    d = {"FeatName": feat_names, "RelFeatImp": feature_importance}
    df = pd.DataFrame(data=d)
    df = df.sort_values(by=["RelFeatImp"], ascending=False)
    pd.DataFrame.to_csv(
        df, os.path.join(output_folder, "feature_importance.txt"), sep="\t", index=False
    )

    # Compute SHAP values if requested
    if compute_shap and best_overall_model is not None and best_X_test is not None:
        print("\nComputing SHAP values for model interpretation...")
        shap_values = compute_shap_values(best_overall_model, best_X_test, feat_names)
        
        if shap_values is not None:
            # Save SHAP plots
            plt.figure(1)  # Summary plot
            plt.savefig(os.path.join(output_folder, "shap_summary_plot.pdf"), bbox_inches="tight", dpi=300)
            plt.savefig(os.path.join(output_folder, "shap_summary_plot.png"), bbox_inches="tight", dpi=300)
            plt.close()
            
            plt.figure(2)  # Waterfall plot
            plt.savefig(os.path.join(output_folder, "shap_waterfall_plot.pdf"), bbox_inches="tight", dpi=300)
            plt.savefig(os.path.join(output_folder, "shap_waterfall_plot.png"), bbox_inches="tight", dpi=300)
            plt.close()
            
            plt.figure(3)  # Bar plot
            plt.savefig(os.path.join(output_folder, "shap_bar_plot.pdf"), bbox_inches="tight", dpi=300)
            plt.savefig(os.path.join(output_folder, "shap_bar_plot.png"), bbox_inches="tight", dpi=300)
            plt.close()
            
            # Save dependence plots if they were created
            for i in range(min(3, len(feature_importance))):
                plt.figure(4+i)
                plt.savefig(os.path.join(output_folder, f"shap_dependence_plot_{i+1}.pdf"), bbox_inches="tight", dpi=300)
                plt.savefig(os.path.join(output_folder, f"shap_dependence_plot_{i+1}.png"), bbox_inches="tight", dpi=300)
                plt.close()
            
            # Save SHAP values as CSV
            shap_df = pd.DataFrame(shap_values.values, columns=feat_names)
            shap_df.to_csv(os.path.join(output_folder, "shap_values.csv"), index=False)
            print("SHAP values and plots saved to output folder.")

    if model_type == "XGB_class":
        # AUC plot
        def create_enhanced_auc_plot(mean_fpr, tprs, aucs, output_folder):
            """Create an enhanced AUC plot with better styling using seaborn."""
            # Set up the figure with seaborn styling
            plt.figure(figsize=(10, 8))
            sns.set_style("whitegrid")
            sns.set_context("talk")
            
            # Calculate metrics
            mean_tpr = np.mean(tprs, axis=0)
            mean_tpr[-1] = 1.0
            mean_auc = auc(mean_fpr, mean_tpr)
            std_auc = np.std(aucs)
            std_tpr = np.std(tprs, axis=0)
            tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
            tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
            
            # Random chance line
            plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='#CC0000', 
                     label='Random chance', alpha=0.8)
            
            # Mean ROC curve
            plt.plot(mean_fpr, mean_tpr, color='#0066CC', 
                     label=f'Mean ROC (AUC = {mean_auc:.2f} ± {std_auc:.2f})',
                     lw=3)
            
            # Standard deviation area
            plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='#0066CC', 
                            alpha=0.2, label='± 1 std. dev.')
            
            # Formatting
            plt.xlim([-0.05, 1.05])
            plt.ylim([-0.05, 1.05])
            plt.xlabel('False Positive Rate', fontsize=14, labelpad=10)
            plt.ylabel('True Positive Rate', fontsize=14, labelpad=10)
            plt.title('Receiver Operating Characteristic (ROC)', fontsize=16, pad=15)
            
            # Enhanced legend with shadow and border
            legend = plt.legend(loc='lower right', fontsize=12, frameon=True, 
                               fancybox=True, framealpha=0.9, facecolor='white', 
                               edgecolor='#CCCCCC')
            
            # Add gridlines
            plt.grid(True, linestyle='--', alpha=0.3)
            
            # Remove top and right spines
            sns.despine(offset=10)
            
            plt.tight_layout()
            plt.savefig(os.path.join(output_folder, "ROC_curve.pdf"), bbox_inches="tight", dpi=300)
            plt.close()

        create_enhanced_auc_plot(mean_fpr, tprs, aucs, output_folder)

        # classification metrics
        median_roc_auc_score = round(np.median(roc_auc_scores), 3)
        median_f1 = round(np.median(f1s), 3)
        median_avg_precision_score = round(np.median(avg_precision_scores), 3)
        median_accuracy = round(np.median(accuracies), 3)
        median_precision = round(np.median(precisions), 3)
        median_recall = round(np.median(recalls), 3)
        
        print("\nDone!\n")
        print(f"***  Median ROC-AUC score was {median_roc_auc_score}  ***")
        print(f"***  Median Accuracy was {median_accuracy}  ***")
        print(f"***  Median Precision was {median_precision}  ***")
        print(f"***  Median Recall was {median_recall}  ***")
        print(f"***  Median F1-score was {median_f1}  ***")
        print(f"***  Median Average Precision score was {median_avg_precision_score}  ***")

        # save results per shuffle
        rez = pd.DataFrame(
            np.column_stack([roc_auc_scores, accuracies, precisions, recalls, f1s, avg_precision_scores]),
            columns=["ROC_AUC_scores", "Accuracy", "Precision", "Recall", "F1-score", "Average_precision_score"],
        )
        pd.DataFrame.to_csv(
            rez,
            os.path.join(output_folder, "model_results_per_iteration.txt"),
            sep="\t",
            index=False,
        )

        # save median metrics over all shuffles
        rr = pd.DataFrame(
            np.column_stack([median_roc_auc_score, median_accuracy, median_precision, 
                            median_recall, median_f1, median_avg_precision_score]),
            columns=["Median AUC", "Median Accuracy", "Median Precision", 
                    "Median Recall", "Median F1-score", "Median Avg. Precision Score"],
        )
        pd.DataFrame.to_csv(
            rr,
            os.path.join(output_folder, "aggregated_metrics_classification.txt"),
            sep="\t",
            index=False,
        )
    else:
        # regression metrics
        median_r2 = round(np.median(r2s), 4)
        print(f"***  Median R2 was {median_r2}  ***")
        median_expl_var = round(np.median(expl_vars), 4)
        print(f"***  Median Explained Variance was {median_expl_var}  ***")
        median_rmse = round(np.median(rmses), 4)
        print(f"***  Median RMSE was {median_rmse}  ***")
        median_mae = round(np.median(maes), 4)
        print(f"***  Median Median Absolute Error was {median_mae}  ***")

        # save results per shuffle
        rez = pd.DataFrame(
            np.column_stack([r2s, expl_vars, rmses, maes]),
            columns=[
                "R2",
                "Explained Variance",
                "RMSE",
                "MAE",
            ],
        )
        pd.DataFrame.to_csv(
            rez,
            os.path.join(output_folder, "model_results_per_iteration.txt"),
            sep="\t",
            index=False,
        )
        # save mean Explained Variance over all shuffles
        rr = pd.DataFrame(
            np.column_stack([median_r2, median_expl_var, median_rmse, median_mae]),
            columns=[
                "Median R2",
                "Median Explained Variance",
                "Median RMSE",
                "Median MAE",
            ],
        )
        pd.DataFrame.to_csv(
            rr,
            os.path.join(output_folder, "aggregated_metrics_regression.txt"),
            sep="\t",
            index=False,
        )

    # plot feature importance with enhanced styling
    def create_enhanced_feature_importance_plot(feature_importance, feat_names, output_folder):
        """Create a visually appealing feature importance plot using seaborn styling"""
        # Create figure and axis explicitly
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Change from whitegrid to white to remove gray background
        sns.set_style("white")
        sns.set_context("talk")
        
        sorted_idx = np.argsort(feature_importance)
        top_sorted_idx = sorted_idx[-30:] if len(sorted_idx) > 30 else sorted_idx
        pos = np.arange(top_sorted_idx.shape[0]) + 0.5
        
        # Create gradient colors based on importance
        importance_norm = (feature_importance[top_sorted_idx] - min(feature_importance[top_sorted_idx])) / \
                          (max(feature_importance[top_sorted_idx]) - min(feature_importance[top_sorted_idx]))
        colors = plt.cm.viridis(importance_norm)
        
        # Create horizontal bar plot - use ax explicitly
        bars = ax.barh(pos, feature_importance[top_sorted_idx], align="center", 
                    color=colors, alpha=0.8, edgecolor='none')
        
        # Add light grid lines but keep white background
        ax.grid(axis='x', linestyle='--', alpha=0.3, color='#CCCCCC')
        ax.set_facecolor('white')
        fig.patch.set_facecolor('white')
        
        # Add gradient color bar with ax parameter
        sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, 
                                  norm=plt.Normalize(vmin=min(feature_importance[top_sorted_idx]), 
                                                    vmax=max(feature_importance[top_sorted_idx])))
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, orientation='vertical', pad=0.01)
        cbar.set_label('Relative Importance Score', rotation=270, labelpad=20, fontsize=12)
        
        # Formatting - use ax explicitly
        ax.set_yticks(pos)
        ax.set_yticklabels(feat_names[top_sorted_idx], fontsize=11)
        ax.set_xlabel("Relative Importance", fontsize=14, labelpad=10)
        title = "Top 30 Most Important Features" if len(sorted_idx) > 30 else "Feature Importance Ranking"
        ax.set_title(title, fontsize=16, pad=15)
        
        # Remove spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(True)
        ax.spines['bottom'].set_visible(True)
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_folder, "feature_importance.pdf"), bbox_inches="tight", dpi=300)
        plt.savefig(os.path.join(output_folder, "feature_importance.png"), bbox_inches="tight", dpi=300)
        plt.close(fig)

    # Call the enhanced feature importance plotting function
    create_enhanced_feature_importance_plot(feature_importance, feat_names, output_folder)

    # save all parameters to file
    params = {**param_grid}
    params["random variables added"] = args.add_random_variables
    params["early_stopping_rounds"] = early_stopping_rounds
    params["N_stab"] = N_stab
    params["test_size"] = test_size
    params["CV_fold"] = CV_fold
    params["rand_seed"] = rand_seed
    params["model_type"] = model_type
    params["name"] = name

    with open(os.path.join(output_folder, "all_model_parameters.txt"), "w") as f:
        f.write(str(params))

    # save conda environment to file that can be used to re-create the environment
    if args.save_conda_env:
        try:
            import conda.cli.python_api as Conda
            
            with open(os.path.join(output_folder, "conda_environment.txt"), "w") as f:
                (stdout_str, stderr_str, return_code_int) = Conda.run_command(
                    Conda.Commands.LIST,
                    ["--explicit"],
                    stdout=f,
                    stderr=sys.stderr,
                )
            print("Conda environment saved successfully.")
        except ImportError:
            print("Warning: conda not found, skipping conda environment export.")
        except Exception as e:
            print(f"Error saving conda environment: {str(e)}")

    # save system configuration to file
    try:
        with open(os.path.join(output_folder, "system_info.txt"), "w") as f:
            import platform
            import psutil

            f.write("=" * 40 + "System Information" + "=" * 40)
            uname = platform.uname()
            f.write("\n\nSystem: " + uname.system)
            f.write("\nNode Name: " + uname.node)
            f.write("\nRelease: " + uname.release)
            f.write("\nVersion: " + uname.version)
            f.write("\nProcessor: " + uname.processor + "\n\n")
            # let's print CPU information
            f.write("=" * 40 + "CPU Info" + "=" * 40)
            # number of cores
            f.write("\n\nPhysical cores: " + str(psutil.cpu_count(logical=False)))
            f.write("\n\nTotal cores: " + str(psutil.cpu_count(logical=True)))
    except Exception as e:
        print(f"Error saving system information: {str(e)}")

    # if tuning run, print best parameters from each iteration to file
    if tune_parameters:
        try:
            with open(
                os.path.join(output_folder + "_best_parameters.txt"), "w", encoding="utf-8"
            ) as best_parameters_file:
                for dic in best_parameters_list:
                    json.dump(dic, best_parameters_file)
                    best_parameters_file.write("\n")
        except Exception as e:
            print(f"Error saving best parameters: {str(e)}")

    end_time = time.time()
    run_time = end_time - start_time
    if run_time < 60:
        print("\n\nTotal Runtime: %.2f seconds" % run_time)
    elif run_time < 3600:
        run_minutes = run_time / 60
        print("\n\nTotal Runtime: %.2f minutes" % run_minutes)
    else:
        run_hours = run_time / 3600
        print("\n\nTotal Runtime: %.2f hours" % run_hours)
    print("\n")
    print("\n" + "*" * 40 + "\n")


if __name__ == "__main__":
    main(sys.argv)