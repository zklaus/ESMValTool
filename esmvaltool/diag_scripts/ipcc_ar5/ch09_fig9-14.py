# -*- coding: utf-8 -*-

"""Diagnostic script to plot figure 9.42a of IPCC AR5 chapter 9.

Description
-----------
Calculate and plot the following quantities with regards to sea
surface temperature: zonal mean error, equatorial mean error,
equatorial mean.  The errors are calculated agains the reference given
in the namelist.  Equatorial here means between 5 degrees north and 5
degrees south.  This has been modelled after IPCC AR5 WG1 Ch. 9,
Fig. 9.14.

Author
------
Klaus Zimmermann (SMHI, Sweden)

Project
-------
CRESCENDO
"""

import os

import iris
import matplotlib
from matplotlib.ticker import MultipleLocator
import numpy as np

# ESMValTool python packages
from auxiliary import info, warning, error  # noqa: F401
from smhi import regrid_esmpy
from smhi.pipeline import (build_var_constraint,
                           load,
                           LoadingStep, ProcessingStep)
from smhi.tools import (
    calc_error,
    ensure_dir_exists,
    get_plot_options,
    LATITUDE_FORMATTER,
    LONGITUDE_FORMATTER,
    multi_model_merge,
    prepare_project_info,
    tag_output_file,
)

matplotlib.use("Agg")
matplotlib.rcParams.update({'font.size': 18})
import iris.plot as iplt  # noqa: E402 # has to be done after use('Agg')
import matplotlib.pyplot as plt  # noqa: E402 # has to be done after use('Agg')


def clim(cube, reference=None, guess_bounds=True):
    cc = cube.collapsed('time', iris.analysis.MEAN)
    cc.convert_units('degC')
    if guess_bounds:
        for coord_name in ['latitude', 'longitude']:
            coord = cc.coord(coord_name)
            if coord.bounds is None:
                try:
                    coord.guess_bounds()
                except:  # noqa: E722
                    error('Guessing bounds for {} failed.'
                          ''.format(cc.attributes['model']))
                    raise
    if reference is not None:
        cc = regrid_esmpy(cc, reference)
    return cc


def zonal_mean(cube, **kwargs):
    zm = cube.collapsed('longitude', iris.analysis.MEAN)
    zm.long_name = 'Zonal mean SST'
    return zm


def equatorial(cube, **kwargs):
    constraint = iris.Constraint(latitude=lambda cell: -5. <= cell.point <= 5.)
    equ = cube.extract(constraint)
    equ = equ.collapsed('latitude', iris.analysis.MEAN)
    equ.long_name = 'Equatorial SST'
    lon = equ.coord('longitude').points
    equ.data.mask[np.logical_and(98. <= lon, lon <= 121.)] = True
    return equ


def setup_chains(E):
    data_dir = os.path.join(E.get_work_dir(), E.get_diag_script_name())
    ensure_dir_exists(data_dir)
    loader = LoadingStep(constraint=build_var_constraint('tos'))
    extractor = ProcessingStep(clim,
                               data_dir=data_dir,
                               source=loader)
    zm = ProcessingStep(zonal_mean,
                        data_dir=data_dir,
                        source=extractor)
    zm_error = ProcessingStep(calc_error,
                              prefix='zonal-mean-error',
                              data_dir=data_dir,
                              source=zm)
    equ = ProcessingStep(equatorial,
                         data_dir=data_dir,
                         source=extractor)
    equ_error = ProcessingStep(calc_error,
                               prefix='equatorial-error',
                               data_dir=data_dir,
                               source=equ)
    chains = {
        'tos': [zm_error, equ, equ_error],
    }
    return chains


def prepare_data(E, modelconfig):
    chains = setup_chains(E)
    models = load(E, chains)
    return models


def setup_figure():
    fig = plt.figure(figsize=(18, 15))
    ax = np.array(
        [[fig.add_axes([0.10, 0.56, 0.30, 0.35]),
          fig.add_axes([0.50, 0.56, 0.30, 0.35])],
         [fig.add_axes([0.10, 0.10, 0.30, 0.35]),
          fig.add_axes([0.50, 0.10, 0.30, 0.35])]]
    )
    return fig, ax


def plot_zonal_mean_errors_ensemble(ax, zonal_mean_errors, ref_line_style):
    ax.set_title('(a) Zonal mean SST error CMIP5')
    ax.yaxis.set_label_text(u'SST error (°C)')
    ax.yaxis.set_minor_locator(MultipleLocator(.5))
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.xaxis.set_major_locator(MultipleLocator(30))
    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.xaxis.set_major_formatter(LATITUDE_FORMATTER)
    ax.set_ylim(-5., 5.)
    ax.set_xlim(-90., 90.)
    ax.tick_params(which='both', direction='in')
    ax.xaxis.set_label_text(u'Latitude')
    ls = []
    labels = []
    cl = multi_model_merge(zonal_mean_errors)
    for e in zonal_mean_errors:
        ls.append(iplt.plot(e, axes=ax)[0])
        labels.append(e.attributes['model'])
    ensemble_mean = cl.collapsed('model', iris.analysis.MEAN)
    m = iplt.plot(ensemble_mean, axes=ax, **ref_line_style)[0]
    ls = [m] + ls
    labels = ['CMIP5 mean'] + labels
    return (ls, labels)


def plot_equatorial_errors(ax, equatorial_errors, ref_line_style):
    ax.set_title('(b) Equatorial SST error CMIP5')
    ax.yaxis.set_label_text(u'SST error (°C)')
    ax.yaxis.set_minor_locator(MultipleLocator(.5))
    ax.xaxis.set_minor_locator(MultipleLocator(30))
    ax.xaxis.set_major_locator(MultipleLocator(60))
    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
    ax.set_ylim(-5., 5.)
    ax.set_xlim(25., 360.)
    ax.tick_params(which='both', direction='in')
    ax.xaxis.set_label_text(u'Longitude')
    for e in equatorial_errors:
        iplt.plot(e, label=e.attributes['model'], axes=ax)
    cl = multi_model_merge(equatorial_errors)
    ensemble_mean = cl.collapsed('model', iris.analysis.MEAN)
    iplt.plot(ensemble_mean, label='CMIP5 mean', axes=ax, **ref_line_style)


def plot_zonal_mean_errors_comparison(ax, zonal_mean_errors, ref_line_style):
    ax.set_title('(c) Zonal mean SST error CMIP5')
    ax.yaxis.set_label_text(u'SST error (°C)')
    ax.yaxis.set_minor_locator(MultipleLocator(.5))
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.xaxis.set_major_locator(MultipleLocator(30))
    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.xaxis.set_major_formatter(LATITUDE_FORMATTER)
    ax.set_ylim(-5., 5.)
    ax.set_xlim(-90., 90.)
    ax.tick_params(which='both', direction='in')
    ax.xaxis.set_label_text(u'Latitude')
    lat = zonal_mean_errors[0].coord('latitude').points
    data = np.ma.vstack([m.data for m in zonal_mean_errors])
    std = data.std(axis=0)
    avg = data.mean(axis=0)
    ax.fill_between(lat, avg-std, avg+std, alpha=.5)
    ax.plot(lat, avg, **ref_line_style)


def plot_equatorials(ax, equatorials, ref_line_style):
    ax.set_title('(d) Equatorial SST CMIP5')
    ax.yaxis.set_label_text(u'SST (°C)')
    ax.yaxis.set_minor_locator(MultipleLocator(.5))
    ax.xaxis.set_minor_locator(MultipleLocator(30))
    ax.xaxis.set_major_locator(MultipleLocator(60))
    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
    ax.set_ylim(22., 31.)
    ax.set_xlim(25., 360.)
    ax.tick_params(which='both', direction='in')
    ax.xaxis.set_label_text(u'Longitude')
    reference = equatorials[0]
    lon = reference.coord('longitude').points
    data = np.ma.vstack([m.data for m in equatorials[1:]])
    std = data.std(axis=0)
    avg = data.mean(axis=0)
    ax.fill_between(lon, avg-std, avg+std, alpha=.5)
    ax.plot(lon, avg, **ref_line_style)
    ls = ax.plot(lon, reference.data, 'k', **ref_line_style)
    return (ls, ['HadISST'])


def draw_legend(fig, ls, labels):
    return fig.legend(ls, labels,
                      loc='upper left',
                      fontsize=16.,
                      bbox_to_anchor=(.81, .92))


def produce_plots(E, modelconfig, data):
    zonal_mean_errors = data['zonal-mean-error_cube'].values[~data['isref']]
    equatorials = data['equatorial_cube'].values[~data['isref']]
    equatorial_errors = data['equatorial-error_cube'].values[~data['isref']]
    cfg, output_path, ref_line_style = get_plot_options(modelconfig, E,
                                                        'flato13_fig9-14')
    plot_dir = os.path.dirname(output_path)
    ensure_dir_exists(plot_dir, 'plot')
    fig, ax = setup_figure()
    ls, labels = plot_zonal_mean_errors_ensemble(ax[0, 0],
                                                 zonal_mean_errors,
                                                 ref_line_style)
    plot_equatorial_errors(ax[0, 1], equatorial_errors, ref_line_style)
    plot_zonal_mean_errors_comparison(ax[1, 0],
                                      zonal_mean_errors, ref_line_style)
    ref_ls, ref_labels = plot_equatorials(ax[1, 1],
                                          equatorials, ref_line_style)
    ls = ref_ls + ls
    labels = ref_labels + labels
    draw_legend(fig, ls, labels)
    fig.savefig(output_path)
    tags = [
        "DM_global", "DM_trop",
        "PT_zonal", "PT_pro",
        "ST_diff", "ST_mean", "ST_stddev", "ST_clim",
    ]
    tag_output_file(output_path, E, "A_zimm_kl",
                    "SST error as both zonal and equatorial mean. "
                    "Similar to Flato et al. 2013, fig. 9.14.",
                    tags)


def main(project_info):
    """
    Arguments
        project_info : Dictionary containing project information

    Description
        This is the main routine of the diagnostic.
    """

    E, modelconfig = prepare_project_info(project_info,
                                          authors=["A_zimm_kl"],
                                          diagnostics=["D_flato13ipcc"],
                                          projects=["P_crescendo"])
    data = prepare_data(E, modelconfig)
    if E.get_write_plots():
        produce_plots(E, modelconfig, data)
