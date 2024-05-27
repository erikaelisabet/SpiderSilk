import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestRegressor
from sklearn.pipeline import Pipeline

# path = "C:/Users/46705/Documents/SpiderSilk/data/post_filtering/arakawa_ara.csv"
# df = pd.read_csv(path, sep=",")

# path_arakawa = "C:/Users/46705/Documents/SpiderSilk/data/post_filtering/arakawa_ara.csv"
# arakawa = pd.read_csv(path_arakawa, sep= ",")
path_MaG = "C:/Users/46705/Documents/SpiderSilk/data/post_filtering/MaG_ara.csv"
path_GC = "C:/Users/46705/Documents/SpiderSilk/data/post_filtering/count_mag_ara.csv"
Gene_c = pd.read_csv(path_GC, sep= ",")
mag =pd.read_csv(path_MaG, sep= ",")

df_1 = mag.iloc[:, :7]
slice_df = mag.iloc[4:, 7:-29]
slice_numeric = slice_df.apply(pd.to_numeric, errors='coerce')
boolean_df = slice_numeric.applymap(lambda x: 0 if x < 2 else 1)
df_3 = mag.iloc[:, -29:]
boolean= pd.concat([df_1, boolean_df, df_3], axis=1)

df= mag
# props = ["strain_at_break", "toughness", "tensile_strength", "young's_modulus"]
props = ["tensile_strength", "young's_modulus"]
for prop in props: 
    features_df = df.iloc[4:, 7:-29]
    properties_df = df.iloc[4:, :][[prop]]
    X = np.array(features_df)
    y = np.array(properties_df)


    pipeline = Pipeline([
        ("model", RandomForestRegressor(random_state=8))
    ])


    param_grid = {
        # "model__n_estimators": [10, 20, 30],
        # "model__max_depth": [5, 10]
        "model__n_estimators": [200, 300, 400, 500],
        "model__max_depth": [5, 10, 15, 20, None] 
    }

    search_RF = GridSearchCV(pipeline, param_grid, cv=5, scoring="neg_mean_squared_error", verbose=3)
    search_RF.fit(X, y)
    results = search_RF.cv_results_

    mean_test_scores = -results['mean_test_score']
    params = results['params']

    # max_depth_values = [5, 10]
    # n_estimators_values = [10, 20, 30]

    max_depth_values = [5, 10, 15, 20, None]
    n_estimators_values = [200, 300, 400, 500]

    plt.figure(figsize=(10, 6))

    for n_estimators in n_estimators_values:
        scores = []
        for max_depth in max_depth_values:
            for mean_test_score, param in zip(mean_test_scores, params):
                if param["model__n_estimators"] == n_estimators and param["model__max_depth"] == max_depth:
                    scores.append(mean_test_score)
                    break
            else:
                scores.append(None) 
        
        plt.plot(max_depth_values, scores, marker='o', linestyle='-', label=f'n_estimators={n_estimators}')

    lowest_mse_index = np.argmin(np.abs(mean_test_scores))
    best_params = search_RF.best_params_

    plt.scatter([params[lowest_mse_index]["model__max_depth"]], [mean_test_scores[lowest_mse_index]], color='r', label='Lowest MSE', zorder=5)
    plt.xlabel('Max Depth')
    plt.ylabel('Mean Test Score (Negative MSE)')
    plt.title(f'GridSearchCV Results {prop}')
    plt.legend()
    plt.tight_layout()
    plt.show()

    plt.savefig(f'C:/Users/46705/Documents/SpiderSilk/output/RF/para_opt_{prop}.png')
    # param_grid = {
    #     "model__n_estimators": [200, 300, 400, 500],
    #     "model__max_depth": [5, 10, 15, 20, None]   
    # }


    # search_RF = GridSearchCV(pipeline, param_grid, cv=5, scoring="neg_mean_squared_error", verbose=3)
    # search_RF.fit(X, y)
    # results = search_RF.cv_results_
    # mean_test_scores = results['mean_test_score']
    # params = results['params']
    # param_values = [str(param) for param in params]


    # lowest_mse_index = np.argmin(np.abs(mean_test_scores))
    # best_params = search_RF.best_params_

    # plt.figure(figsize=(10, 6))
    # plt.plot(range(len(param_values)), mean_test_scores, marker='o', linestyle='-', color='b')
    # plt.scatter(lowest_mse_index, mean_test_scores[lowest_mse_index], color='r', label='Lowest MSE', zorder=5)
    # plt.xlabel('Parameter Combination')
    # plt.ylabel('Mean Test Score (Negative MSE)')
    # plt.title(f'GridSearchCV Results {prop}')
    # plt.xticks(range(len(param_values)), param_values, rotation=45)
    # plt.tight_layout()
    # plt.legend()

    # plt.savefig(f'C:/Users/46705/Documents/SpiderSilk/output/RF/para_opt_{prop}.png')

    # plt.show()

    # print("Best Parameters:", best_params)

