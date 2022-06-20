import numpy as np
import pandas as pd
import COPASI


def dropConstantCol(df):
    return df.loc[:, (df != df.iloc[0]).any()]


def classifyPathway(df, e=0.1, E1='ET_PNCA', E2='ET_NAMPT'):
    df.loc[df[E1] <= e, 'pathway'] = E2[3:]
    df.loc[df[E2] <= e, 'pathway'] = E1[3:]
    return df


def appendDataframe(df1, df2, E1='PNCA', E2='NAMPT'):
    df1.loc[:, 'pathway'] = E1
    df2.loc[:, 'pathway'] = E2
    return df1.append(df2)


def selectPathway(pathway: str):
    if pathway == 'PNCA':
        inactive = 'NAMPT'
    elif pathway == 'NAMPT':
        inactive = 'PNCA'
    else:
        raise ValueError('pathway expected PNCA or NAMPT')
    return inactive


def setValue(dataModel: COPASI.CDataModel, parameterName: str, parameterValue: float):
    model = dataModel.getModel()
    mv = model.getModelValue(parameterName)
    if mv is None:
        raise ValueError(f'{parameterName} not found')
    else:
        mv.setInitialValue(parameterValue)
        assert mv.getInitialValue(
        ) == parameterValue, f'parameter value could not be set for {parameterName}'


def arrheniusEquation(T, A, Ea):
    R = 8.3145
    A = A * np.ones(len(T))
    Ea = Ea * np.ones(len(T))
    return A * np.exp(-Ea/(R * T))


def Ea():
    Ea = {'old': {'NaMN': 249.4,
                  'NAD': 90.2,
                  'NAR': 149.3,
                  'NMN': 93.9,
                  'NR': 106.1},
          'new': {'NaMN': 142.370957,
                  'NAD': 89.749386,
                  'NAR': 138.489012,
                  'NMN': 103.536009,
                  'NR': 72.358094}}
    return Ea


def A():
    A = {'old': {'NaMN': 4.978e31,
                 'NAD': 2.531e9,
                 'NAR': 4.36e17,
                 'NMN': 1.361e10,
                 'NR': 2.463e12},
         'new': {'NaMN': 2.085151e13,
                 'NAD': 2.197221e6,
                 'NAR': 1.224697e13,
                 'NMN': 3.639423e8,
                 'NR': 2.680513e4}}
    return A


def prepET(filepath_optET, objective):
    enames = {'ET_NRK1': 'NadR', 'ET_NT5': 'SurE', 'ET_NMNAT1': 'NadD', 'ET_NAPRT': 'PncB', 'ET_NADS': 'NadE',
              'ET_NAMPT': 'Nampt', 'ET_PNCA': 'PncA', 'ET_PNP':'PNP', 'ET_SIRT':'NCE'}
    df = pd.read_csv(filepath_optET, sep='\t', index_col=0)
    if not 'objective' in df.columns:
        df['objective'] = objective
    df = df.rename(columns=enames)
    try:
        df.pathway = df.pathway.replace({'NAMPT': 'Nampt', 'PNCA': 'PncA'})
    except AttributeError:
        df.loc[(df.PncA < 1), 'pathway'] = 'Nampt'
        df.loc[(df.Nampt < 1), 'pathway'] = 'PncA'
    df1 = df.replace([np.inf, -np.inf], np.nan).dropna(axis=0).reset_index()
    return df1
