from fileinput import filename
import numpy as np
import seaborn as sns
from itertools import combinations
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import utilities as u
import simulation as sim
import matplotlib.ticker as mticker


def plot_corr(dfP, dfN, figsize=(15, 6), fsize=20, filename=False, sharey=True):
    corr = dfP.corr()
    mask = np.zeros_like(corr)
    mask[np.triu_indices_from(mask)] = True
    fig, ax = plt.subplots(1, 2, figsize=figsize, sharey=sharey)
    ax[0].set_title('Pathway: PNCA', fontsize=fsize, loc='left')
    sns.heatmap(data=corr, cmap='seismic', center=0.0,
                annot=False, mask=mask, ax=ax[0], cbar=True),
    # cbar_kws={'label': 'Pearsons correlation\ncoefficient (r)'});

    corr = dfN.corr()
    mask = np.zeros_like(corr)
    mask[np.triu_indices_from(mask)] = True
    ax[1].set_title('Pathway: NAMPT', fontsize=fsize, loc='left')
    sns.heatmap(data=corr, cmap='seismic', center=0.0, annot=False, mask=mask,
                ax=ax[1], cbar_kws={'label': 'Pearsons correlation\ncoefficient (r)'})
    plt.tight_layout()
    if filename:
        fig.savefig(filename, dpi=300)
    return plt.show()


def subplot_lineplot(df, layout: tuple, figsize: tuple, xlabel='Temperature ($\degree$C)', color='blue',
                     palette=['#FF0000', '#1E90FF'], **kwargs):
    fig, axn = plt.subplots(layout[0], layout[1], figsize=figsize, sharex=kwargs.get(
        'sharex', False), sharey=kwargs.get('sharey', False))
    for i, ax in enumerate(axn.flat):
        if i < len(df.columns):
            if 'hue' in kwargs:
                lp = sns.lineplot(data=df, y=df[df.columns[i]], x=df.index, ax=ax, palette=palette,
                                  hue=kwargs['hue'], err_style=kwargs.get('err_style', "bars"))
                lp.get_legend().remove()
            else:
                lp = sns.lineplot(
                    data=df, y=df.columns[i], x=df.index, ax=ax, color=color)
            if kwargs.get('sharey', False) == True:
                lp.set(xlabel=xlabel)
                ax.set_title(df.columns[i])
            else:
                lp.set(xlabel=xlabel, ylabel=df.columns[i])
        else:
            break
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right', framealpha=1.0)
    plt.tight_layout()
    if 'filename' in kwargs:
        fig.savefig(kwargs['filename'], dpi=300)
    return plt.show()


def arrheniusEquation(T, A, Ea):
    R = 8.3145
    A = A * np.ones(len(T))
    Ea = Ea * np.ones(len(T))
    return A * np.exp(-Ea/(R * T))


def plot_thermolysis(drates, save=False, figsize=(10, 8), filename='../images/thermolysis_rates.png'):
    legend_map = {'NR': 'NR ($E_a$ = 72.41 KJ/mol)',
                  'NAR': 'NAR ($E_a$ = 138.6 KJ/mol)',
                  'NMN': 'NMN ($E_a$ = 103.62 KJ/mol)',
                  'NaMN': 'NaMN ($E_a$ = 142.49 KJ/mol)',
                  'NAD': 'NAD ($E_a$ = 123.12 KJ/mol)',
                  'NAAD': 'NAAD'}
    fig, ax = plt.subplots(figsize=figsize)
    drates.loc[:, 'Compound'] = [i.upper() for i in drates.Compound]
    data = drates[~drates.Compound.isin(['NAM', 'NIA', 'ATP'])]
    sns.scatterplot(data=data, x='Temperature (K)', y='Rate (1/s)', hue=data['Compound'],
                    hue_order=['NR', 'NMN', 'NAD', 'NAR', 'NAMN', 'NAAD'], ax=ax)
    plt.legend(ncol=2, title='Compound')
    for c in ['NR', 'NMN', 'NAD', 'NAR', 'NAMN', 'NAAD']:
        ax.plot(np.linspace(320, 365, 100),
                arrheniusEquation(np.linspace(320, 365, 100), A=data[data.Compound == c].prefactor_A.mean()*1e3,
                                  Ea=data[data.Compound == c]['activation_energy (KJ/mol)'].mean()*1e3))
    ax.ticklabel_format(style='sci', axis='y',
                        scilimits=(-3, 6), useMathText=True)

    plt.tight_layout()
    if save == True:
        plt.savefig(filename, dpi=300)
    return plt.show()


def subplots_concentration(df, layout: tuple, figsize: tuple, xlabel='Temperature ($\degree$C)', color='blue',
                           palette=['#FF0000', '#1E90FF'], **kwargs):
    mets = [i+' (mM)' for i in ['NAD', 'NA', 'NAMN',
                                'NR', 'NAM', 'NMN', 'NAAD', 'NAR']]
    fig, axn = plt.subplots(layout[0], layout[1], figsize=figsize, sharex=kwargs.get(
        'sharex', False), sharey=kwargs.get('sharey', False))
    for i, ax in enumerate(axn.flat):
        if i < len(df.columns):
            if 'hue' in kwargs:
                lp = sns.lineplot(data=df, y=mets[i], x=df.temperature, ax=ax, palette=palette,
                                  hue=kwargs['hue'], err_style=kwargs.get('err_style', "bars"))
                lp.lines[0].set_linestyle('dashed')
                lp.get_legend().remove()
            else:
                lp = sns.lineplot(
                    data=df, y=mets[i], x=df.temperature, ax=ax, color=color)
            if kwargs.get('sharey', False) == True:
                lp.set(xlabel=xlabel)
                ax.set_title(df.columns[i])
            else:
                lp.set(xlabel=xlabel, ylabel=mets[i])
        else:
            break
        ax.ticklabel_format(style='scientific', axis='y',
                            scilimits=(-1, 1), useMathText=True)
        if mets[i] in [j+' (mM)' for j in ['NAD', ]]:
            ax.set_ylim([-5e-4, 5e-1])
        elif mets[i] in [j+' (mM)' for j in ['NMN', ]]:
            ax.set_ylim([-5e-4, 5e-2])

    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right', framealpha=1.0)
    plt.tight_layout()
    if 'filename' in kwargs:
        fig.savefig(kwargs['filename'], dpi=300)
    return plt.show()


def subplot_fluxes(F, save=False, filename='', **kwargs):
    F = sim.calculateATPconsFlux(F).reset_index()
    F.pathway = F.pathway.replace({'NAMPT': 'Nampt', 'PNCA': 'PncA'})

    y = ['ATP consumption', 'NAD production', 'ATPconsNADprod']
    y_label = ['$J_{ATP}$', '$J_{NAD}$', '$J_{ATP}/J_{NAD}$']
    y_lim = [kwargs.get('ylim1', [-5e-6, 2e-4]), kwargs.get('ylim2', [-1e-6, 5e-5]),
             kwargs.get('ylim3', [-2e-1, 6.5])]

    fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    for i, ax in enumerate(ax.flat):
        sns.lineplot(data=F, x=F.temperature, y=y[i], hue='pathway', hue_order=['PncA', 'Nampt'],
                     palette=['#FF0000', '#1E90FF'], ax=ax, err_style='bars', dashes=True)
        ax.set_ylabel(y_label[i], fontsize=20)
        ax.set_xlabel('Temperature ($\degree$C)', fontsize=20)
        ax.ticklabel_format(style='scientific', axis='y',
                            scilimits=(-3, 6), useMathText=True)
        ax.set_ylim(y_lim[i])
        ax.legend(loc='upper left')
    plt.tight_layout()
    if save == True:
        fig.savefig(filename, dpi=300)
    return plt.show()


def subplots(df_, scaling, objective, title_fsize=20, **kwargs):
    g = sns.jointplot(data=df_, x=df_[
                      "Nampt"]*scaling, y=df_["PncA"]*scaling, hue="objective", space=0.5, height=5)
    g.set_axis_labels('Nampt (nM)', 'PncA (nM)')
    axins1 = inset_axes(g.ax_joint, width="120%", height="110%",
                        bbox_to_anchor=(1.5, 0, 1, 1),
                        bbox_transform=g.ax_joint.transAxes, loc=3, borderpad=0)
    axins2 = inset_axes(g.ax_joint, width="120%", height="110%",
                        bbox_to_anchor=(3.0, 0, 1, 1),
                        bbox_transform=g.ax_joint.transAxes, loc=3, borderpad=0)
    cols = ['PncA', 'Nampt', 'NadE', 'PncB', 'NadD']
    fdf = df_.loc[(df_.objective == objective) & (df_.pathway == 'PncA')]
    fdf = fdf.set_index('Temperature')[cols]*scaling
    lp = sns.lineplot(data=fdf, dashes=False, markers=True, ax=axins1,
                      palette=['#FF0000', '#1E90FF', '#FFA500', '#2E8B57', '#000080'], ci=98)
    lp.legend(ncol=2)
    axins1.set_ylabel('Concentration (nM)')
    axins1.set_title('Pathway: PncA', fontsize=title_fsize)
    axins1.set_ylim(kwargs.get('ylim1', [-2, 1e2]))
    fdf = df_.loc[(df_.objective == objective)
                  & (df_.pathway == 'Nampt')]
    fdf = fdf.set_index('Temperature')[cols]*scaling
    lp = sns.lineplot(data=fdf, dashes=False, markers=True, ax=axins2,
                      palette=['#FF0000', '#1E90FF', '#FFA500', '#2E8B57', '#000080'], ci=98)
    lp.legend(ncol=2)
    axins2.set_ylabel('Concentration (nM)')
    axins2.set_title('Pathway: Nampt', fontsize=title_fsize)
    axins2.set_ylim(kwargs.get('ylim2', [-2, 1e2]))
    return plt.show()


def controlcoefficients(CC, save=False, filename=''):
    """plots concentration control coefficients as subplots

    Args:
        CC (pandas.DataFrame): control coefficients

    Returns:
        _type_: shows plot
    """
    fig, axn = plt.subplots(3, 4, figsize=(15, 9))
    for i, ax in enumerate(axn.flat):
        sp = sns.lineplot(
            data=CC, x=CC['temperature'], y=CC[CC.columns[1+i]], hue='metabolite', ax=ax)
        ax.ticklabel_format(style='scientific', axis='y',
                            scilimits=(-3, 6), useMathText=True)
        sp.get_legend().remove()
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, ncol=7, loc='upper right', framealpha=1.0)
    plt.tight_layout()
    if save == True:
        fig.savefig(filename, dpi=300)
    return plt.show()
