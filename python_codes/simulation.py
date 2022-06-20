from tabnanny import verbose
import pycoexp.tasks
import pandas as pd
import re
import os
import numpy as np
import shutil
import scipy.optimize as optimize
import utilities as u
import plot as pl
task = pycoexp.tasks.tasks()


def setEa_andA(filepath_CPSmodel, Ea, A, filepath_updated):
    A['NaMN'] = A['NAMN']
    Ea['NaMN'] = Ea['NAMN']
    dataModel = task.init_dataModel(filepath_CPSmodel)
    for compound in ['NAD', 'NaMN', 'NAR', 'NMN', 'NR']:
        try:
            u.setValue(dataModel=dataModel, parameterName='Ea_' +
                       compound, parameterValue=Ea[compound])
            u.setValue(dataModel=dataModel, parameterName='A_' +
                       compound, parameterValue=A[compound])
        except ValueError:
            pass
    return dataModel.saveModel(filepath_updated, True)


def initializeET(filepath_CPSmodel, initialConc: float, newpars: dict, filepath_updated, prefix='ET'):
    dataModel = task.init_dataModel(filepath_CPSmodel=filepath_CPSmodel)
    enzymes = ['NAPRT', 'NMNAT1', 'NADS', 'NRK1',
               'NT5', 'PNP', 'SIRT', 'PNCA', 'NAMPT']
    epars = [f'{prefix}_{e}' for e in enzymes]
    if isinstance(initialConc, (int, float)):
        evalues = {par: initialConc for par in epars}
    elif initialConc is None:
        model = dataModel.getModel()
        evalues = {par: model.getModelValue(
            par).getInitialValue() for par in epars}
    else:
        raise ValueError(
            f'initialConc expects data of type int or float but got type {type(initialConc)}')
    if newpars is None:
        pass
    else:
        evalues.update(newpars)
    for par in evalues:
        u.setValue(dataModel=dataModel, parameterName=par,
                   parameterValue=evalues[par])
    return dataModel.saveModel(filepath_updated, True)


def fitArrheniusEquation(thermolysisRate, **kwargs):
    for compound in set(thermolysisRate.Compound):
        xdata = thermolysisRate.loc[thermolysisRate.Compound ==
                                    compound]['Temperature (K)']
        ydata = thermolysisRate.loc[thermolysisRate.Compound ==
                                    compound]['rate (1/s)']
        if compound in kwargs.get('weighted_compounds', ['NMN', 'NAMN', 'ATP']):
            if 'data_indices' in kwargs:
                data_indices = kwargs['data_indices']
            else:
                raise ValueError('missing argument: data_indices')
            sigma = np.ones(len(xdata))
            sigma[data_indices] = kwargs.get('weight', 0.01)
            popt, pcov = optimize.curve_fit(
                f=u.arrheniusEquation, xdata=xdata, ydata=ydata, maxfev=100000, sigma=sigma)
        else:
            popt, pcov = optimize.curve_fit(
                f=u.arrheniusEquation, xdata=xdata, ydata=ydata, maxfev=100000)
        thermolysisRate.loc[thermolysisRate.Compound ==
                            compound, 'prefactor_A'] = popt[0]/1000.
        thermolysisRate.loc[thermolysisRate.Compound ==
                            compound, 'activation_energy (KJ/mol)'] = popt[1]/1000.
    return thermolysisRate


def scan(parameterName, filepath_CPSmodel, **kwargs):
    '''
    scan a range of temperature to calculate steady state
    :param filepath_CPSmodel: file path to the COPASI model
    :type filepath_CPSmodel: str
    :param pathway: 'PNCA' or 'NAMPT'
    :type pathway: str
    :param initialConc: initial enzyme concentration
    :type initialConc: float or int
    :param kwargs: T0, lower bound for the temperature range; T1, upper bound for the temperature range; steps,
    number of steps between the lower and upper bound
    :type kwargs: float, float, int
    :return: steady state concentrations and fluxes
    :rtype: tuple of pandas dataframes
    '''
    # create a folder to save CPS models
    foldername = kwargs.get('foldername', 'tmp/')
    os.makedirs(foldername + f'{parameterName}_scan', exist_ok=True)

    p0 = kwargs.get('p0', 30.)
    p1 = kwargs.get('p1', 90.)
    steps = kwargs.get('steps', 5)
    prange = np.linspace(p0, p1, int(steps))
    conc, flux = [], []
    for i in range(int(steps)):
        #print(f'{parameterName}: {prange[i]}')
        filepath_cpsmodel = os.path.join(
            foldername + f'{parameterName}_scan', f'{parameterName}_{prange[i]}.cps')
        initialize = kwargs.get('initialize', True)
        if initialize == True:
            # initialize enzyme concentrations in the model
            initializeET(filepath_CPSmodel=filepath_CPSmodel, newpars=kwargs.get('newpars', None),
                         initialConc=kwargs.get('initialConc', None),
                         filepath_updated=filepath_cpsmodel)
        else:
            pass
        newEa_andA = kwargs.get('newEa_andA', False)
        if newEa_andA == True:
            setEa_andA(filepath_CPSmodel=filepath_cpsmodel, Ea=kwargs['Ea'], A=kwargs['A'],
                       filepath_updated=filepath_cpsmodel)
        else:
            pass
        dataModel = task.init_dataModel(filepath_CPSmodel=filepath_cpsmodel)
        u.setValue(dataModel=dataModel, parameterName=parameterName,
                   parameterValue=prange[i])
        dataModel.saveModel(filepath_cpsmodel, True)
        concentrations, fluxes = task.steadystate(
            filepath_CPSmodel=filepath_cpsmodel)
        conc.append(concentrations)
        flux.append(fluxes)

    conc = pd.DataFrame(conc)
    conc[parameterName] = prange
    flux = pd.DataFrame(flux)
    flux[parameterName] = prange
    deleteModels = kwargs.get('deleteModels', False)
    if deleteModels == True:
        shutil.rmtree(foldername + f'{parameterName}_scan')
    return conc.set_index(parameterName), flux.set_index(parameterName)


def scan_and_optimize(parameterName, filepath_CPSmodel, **kwargs):
    # create a folder to save CPS models and results
    foldername = kwargs.get('foldername', 'tmp_opt/')
    os.makedirs(foldername + f'{parameterName}_scan', exist_ok=True)
    os.makedirs(foldername + f'{parameterName}_scan/results', exist_ok=True)
    T0 = kwargs.get('T0', 20.)
    T1 = kwargs.get('T1', 90.)
    steps = kwargs.get('steps', 5)
    n = kwargs.get('n', 5)
    Trange = np.linspace(T0, T1, int(steps))
    for i in range(int(steps)):
        filepath_cpsmodel = os.path.join(
            foldername + f'{parameterName}_scan', f'{parameterName}_{Trange[i]}.cps')
        initialize = kwargs.get('initialize', True)
        if initialize == True:
            # initialize enzyme concentrations in the model
            initializeET(filepath_CPSmodel=filepath_CPSmodel, newpars=kwargs.get('newpars', None),
                         initialConc=kwargs.get('initialConc', None),
                         filepath_updated=filepath_cpsmodel)
        else:
            pass
        newEa_andA = kwargs.get('newEa_andA', False)
        if newEa_andA == True:
            setEa_andA(filepath_CPSmodel=filepath_cpsmodel, Ea=kwargs['Ea'], A=kwargs['A'],
                       filepath_updated=filepath_cpsmodel)
        else:
            pass
        dataModel = task.init_dataModel(filepath_CPSmodel=filepath_cpsmodel)
        u.setValue(dataModel=dataModel, parameterName=parameterName,
                   parameterValue=Trange[i])
        dataModel.saveModel(filepath_cpsmodel, True)
        Sol, E = [], []
        for j in range(n):
            value, result = task.optimization(
                filepath_CPSmodel=filepath_cpsmodel)
            Sol.append(value)
            E.append(result)
        df = pd.DataFrame(np.array(E), columns=readOptReport(
            foldername+f'{parameterName}_scan/test-results-opt.txt'))
        df['Sol'] = Sol
        df.to_csv(
            foldername + f'{parameterName}_scan/results/{parameterName}_'+str(Trange[i])+'.csv', sep='\t')
    deleteModels = kwargs.get('deleteModels', False)
    if deleteModels == True:
        shutil.rmtree(foldername + f'{parameterName}_scan')
    return readOptResults(foldername + f'{parameterName}_scan/results/')


def getnewCPS(filepath_CPSmodel, Optdf, foldername='tmp_opt/optimized'):
    # create a folder to save CPS models
    os.makedirs(foldername, exist_ok=True)

    for i in Optdf.index:
        # create dataModel object
        dataModel = task.init_dataModel(filepath_CPSmodel=filepath_CPSmodel)
        for j in range(len(Optdf.T)):
            par = Optdf.columns[j]
            value = Optdf.T[i][Optdf.columns[j]]
            u.setValue(dataModel=dataModel, parameterName=par,
                       parameterValue=value)
        filepath_cpsmodel = os.path.join(
            foldername, f'optimized_{i}.cps')
        dataModel.saveModel(filepath_cpsmodel, True)


def getsteadystates(filepath_CPSmodel, Optdf, **kwargs):
    Optdf = Optdf.rename(columns={'Temperature': 'temperature'})
    foldername = 'tmp_opt/optimized'
    # create a folder to save CPS models
    os.makedirs(foldername, exist_ok=True)

    conc, flux = [], []
    temp = []
    for i in range(len(Optdf)):
        # create dataModel object
        dataModel = task.init_dataModel(filepath_CPSmodel=filepath_CPSmodel)
        for j in range(len(Optdf.T)):
            par = Optdf.columns[j]
            value = Optdf.T[i][Optdf.columns[j]]
            u.setValue(dataModel=dataModel, parameterName=par,
                       parameterValue=value)
        filepath_cpsmodel = os.path.join(
            foldername, f'optimized_{Optdf.T.columns[i]}.cps')
        dataModel.saveModel(filepath_cpsmodel, True)
    #     newEa_andA = kwargs.get('newEa_andA', False)
    #     if newEa_andA == True:
    #         setEa_andA(filepath_CPSmodel=filepath_cpsmodel, Ea=kwargs['Ea'], A=kwargs['A'],
    #                    filepath_updated=filepath_cpsmodel)
    #     else:
    #         pass
    #     try:
    #         concentrations, fluxes = task.steadystate(
    #             filepath_CPSmodel=filepath_cpsmodel)
    #         conc.append(concentrations)
    #         flux.append(fluxes)
    #         temp.append(Optdf.iloc[i]['temperature'])
    #     except AssertionError:
    #         print(f'Steady-state not found for {i}')
    #         pass

    # conc = pd.DataFrame(conc)
    # conc['temperature'] = temp
    # flux = pd.DataFrame(flux)
    # flux['temperature'] = temp
    # deleteModels = kwargs.get('deleteModels', False)
    # if deleteModels == True:
    #     shutil.rmtree(foldername)
    # return conc.set_index('temperature'), flux.set_index('temperature')


def readOptReport(filepathOptReport):
    file = open(filepathOptReport, 'r')
    lines = file.readlines()
    optItems = []
    for i in range(len(lines)):
        if 'List of Optimization Items:' in lines[i]:
            for j in lines[i + 1:]:
                if not 'List of Constraint Items:' in j:
                    optItems.append(j)
                else:
                    break
        else:
            pass
    return [item.split('<=')[1].strip(' ')[len('Values['):-len('].InitialValue')] for item in optItems[:-1]]


def readOptResults(foldername):
    df = pd.DataFrame()
    for file in os.listdir(foldername):
        if file.endswith('.csv') and file.startswith('temperature'):
            df_ = pd.read_csv(foldername+file, delimiter='\t')
            df_['Temperature'] = float(file[:-4].split('_')[1])
            df = df.append(df_)
    df = df.drop(df.columns[0], axis=1).set_index('Temperature')
    return df


def calculateATPconsFlux(F):
    F['ATP consumption'] = F['NAPRT'] + 2 * F['NADS'] + F['NMNAT1-NaMN'] + F['NAMPT'] + F['NMNAT1-NMN'] + \
        F['NRK1-NMN'] + F['NRK1-NaMN']
    F['NAD production'] = F['NADS'] + F['NMNAT1-NMN']
    F['ATPconsNADprod'] = F['ATP consumption']/F['NAD production']
    return F


def atpconsNADflux(F):
    F['ATP consumption(NADA)'] = F.NAPRT + 2 * F.NADS + F['NMNAT1-NaMN'] + F.NAMPT + F['NMNAT1-NMN'] + F['NRK1-NMN'] + \
        F['NRK1-NaMN']
    F['NAD production(NADA)'] = F.NADS + F['NMNAT1-NMN']
    F['ATP consumption(NAMPT)'] = F['NAPRT[NAMPT]'] + 2 * F['NADS[NAMPT]'] + F['NMNAT1-NaMN[0]'] + F['NAMPT[NAMPT]'] + \
        F['NMNAT1-NMN[0]'] + F['NRK1-NMN[0]'] + F['NRK1-NaMN[0]']
    F['NAD production(NAMPT)'] = F['NADS[NAMPT]'] + F['NMNAT1-NMN[0]']

    F['PNCA'] = F['ATP consumption(NADA)'] / F['NAD production(NADA)']
    F['NAMPT'] = F['ATP consumption(NAMPT)'] / F['NAD production(NAMPT)']
    return F


def NADAvsNAMPT(df, units=' (mM)'):
    dfC_ = df.loc[:, (df != df.iloc[0]).any()]
    nada, nampt = '{NADA}', '{NAMPT}'
    dfC_NADA = dfC_.filter(like=nada)
    dfC_NAMPT = dfC_.filter(like=nampt)
    dfC_NADA = dfC_NADA.rename(
        columns={i: i.split(nada)[0]+units for i in dfC_NADA.columns})
    dfC_NAMPT = dfC_NAMPT.rename(
        columns={i: i.split(nampt)[0]+units for i in dfC_NAMPT.columns})
    dfC_NADA['pathway'] = 'NADA'
    dfC_NAMPT['pathway'] = 'NAMPT'
    dfC11 = dfC_NADA.append(dfC_NAMPT)
    return dfC11


def Kole_exp(exp, genes):
    fdf = exp.loc[exp['Locus Tag'].isin(genes.Gene_name)]
    Kole_genes = genes.set_index('Gene_name').join(fdf.set_index('Locus Tag'))
    Kole = Kole_genes.reset_index().drop(['Description', 'Feature ID', 'UniProtID', 'Gene_name', 'Gene Annotation',
                                          'COGs', 'COG Categories', 'Pfams', 'TIGRfams', 'PC1', 'PC2',
                                          'Length of vector'], axis=1)
    Kole = Kole.set_index('Name')
    Kole.columns = [int(i) for i in [re.findall('\d+', i)[0]
                                     for i in Kole.columns]]
    Kole = Kole.T.sort_index().T
    return Kole


def performMCA(Optdf, filepath_CPSmodel, system_variable, foldername='./', deleteModels=False):
    # create a folder to save CPS models
    os.makedirs(foldername+'mca', exist_ok=True)

    CC = pd.DataFrame()

    for i in range(len(Optdf)):
        # create dataModel object
        dataModel = task.init_dataModel(filepath_CPSmodel=filepath_CPSmodel)
        model = dataModel.getModel()
        for j in range(len(Optdf.T)):
            mv = model.getModelValue(Optdf.columns[j])
            assert mv != None, "Parameter {} not found".format(
                Optdf.columns[j])
            mv.setInitialValue(Optdf.T[i][Optdf.columns[j]])
            assert mv.getInitialValue() == Optdf.T[i][Optdf.columns[j]]

        #assert count == len(df), "Successfully integrated {} of {} genes".format(count, len(df))
        filepath_cpsmodel = os.path.join(
            foldername, 'mca/optimized_' + str(Optdf.T.columns[i]) + '.cps')
        dataModel.saveModel(filepath_cpsmodel, True)

        model = dataModel.getModel()
        temp = model.getModelValue('temperature')

        cc = task.mca(filepath_CPSmodel=filepath_cpsmodel,
                      system_variable=system_variable, verbose=True)
        cc['temperature'] = round(temp.getInitialValue(), 2)
        CC = CC.append(cc)

    if deleteModels == True:
        shutil.rmtree(foldername)

    return CC


def get_steadystate_CF(df, model='../models/new/NAD_biosynthesis_PN_mATPNAD_NAD.cps'):
    CP, FP = getsteadystates(filepath_CPSmodel=model,
                             Optdf=df[df.pathway == 'PNCA'].drop(['Sol', 'objective', 'pathway'], axis=1))
    CP['pathway'] = 'PNCA'
    FP['pathway'] = 'PNCA'
    CN, FN = getsteadystates(filepath_CPSmodel=model,
                             Optdf=df[df.pathway == 'NAMPT'].drop(['Sol', 'objective', 'pathway'], axis=1))
    CN['pathway'] = 'NAMPT'
    FN['pathway'] = 'NAMPT'
    C = pd.concat((CP, CN))
    F = pd.concat((FP, FN))
    F = calculateATPconsFlux(F)
    return C, F
