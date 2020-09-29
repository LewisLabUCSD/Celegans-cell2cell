import scipy

import numpy as np
import pandas as pd
import seaborn as sns

def get_correlation(physical_distance_matrix, cci_distance_matrix, corr_type='spearman'):
    cci_ = cci_distance_matrix.copy()
    np.fill_diagonal(cci_.values, 0.0)
    if (cci_.values.transpose() == cci_.values).all():
        phy_distance_vector = scipy.spatial.distance.squareform(physical_distance_matrix)
        cci_distance_vector = scipy.spatial.distance.squareform(cci_)
    else:
        phy_distance_vector = physical_distance_matrix.values.flatten()
        cci_distance_vector = cci_.values.flatten()
    if corr_type == 'spearman':
        return scipy.stats.spearmanr(phy_distance_vector, cci_distance_vector)
    elif corr_type == 'pearson':
        return scipy.stats.pearsonr(phy_distance_vector, cci_distance_vector)
    else:
        raise NotImplementedError("Correlation {} is not implemented.".format(corr_type))

def get_correlation_plot(physical_distance_matrix, cci_distance_matrix, transform_type=None, p=None):
    cci_ = cci_distance_matrix.copy()
    np.fill_diagonal(cci_.values, 0.0)
    if (cci_.values.transpose() == cci_.values).all():
        phy_distance_vector = scipy.spatial.distance.squareform(physical_distance_matrix)
        cci_distance_vector = scipy.spatial.distance.squareform(cci_)
    else:
        phy_distance_vector = physical_distance_matrix.values.flatten()
        cci_distance_vector = cci_.values.flatten()
    data = np.array(list(zip(phy_distance_vector, cci_distance_vector)))
    df = pd.DataFrame(data, columns=['X', 'Y'])
    x = 'X'
    y = 'Y'
    if transform_type is not None:
        x = 'Transformed'
        if transform_type == 'log':
            df[x] = df['X'].apply(lambda x: np.log(x))
        elif transform_type == 'reciprocal_squared':
            if p is None:
                p = 2
            df[x] = df[['X', 'Y']].apply(lambda row: 1.0 if (row['X'] < 1) else 1.0 / (row['X'] ** p), axis=1)
        elif transform_type == 'power':
            if p is None:
                p = 2
            df[x] = df[['X', 'Y']].apply(lambda row: row['X'] ** p, axis=1)
        elif transform_type == 'exp_decay':
            if p is None:
                p = 200
            df[x] = df[['X', 'Y']].apply(lambda row: np.exp((-1*(row['X']**2))/(2.*(p**2))), axis=1)
        else:
            raise ValueError('transform_type not implemented')
    df = df.replace([np.inf, -np.inf], np.nan).dropna()

    fig = sns.regplot(data=df, x=x, y=y, color='r', fit_reg=True, scatter_kws={'s': 8, 'marker': 'o', 'color': 'b'})
    corr = {'pearson': scipy.stats.pearsonr(df[x].values, df[y].values)[0],
            'spearman': scipy.stats.spearmanr(df[x].values, df[y].values)[0]
            }
    return fig, df, corr
