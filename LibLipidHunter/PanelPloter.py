# -*- coding: utf-8 -*-
# Copyright 2015-2017 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

from __future__ import division

import os

import matplotlib

matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4'] = 'PySide'

from matplotlib import pyplot as plt
import matplotlib.patches as patches
import matplotlib as mpl
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from matplotlib._png import read_png
import pandas as pd


def plot_spectra(mz_se, xic_dct, ident_info_dct, spec_info_dct, specific_check_dct, isotope_checker_dct, isotope_score,
                 save_img_as=None, ms1_precision=50e-6):
    ms2_pr_mz = mz_se['mz']
    ms1_obs = mz_se['MS1_obs_mz']
    lib_mz = mz_se['Lib_mz']
    abbr_id = mz_se['Abbreviation']
    func_id = mz_se['function']
    ms1_pr_ppm = mz_se['ppm']
    # _usr_formula = mz_se['Formula']
    # _usr_ms2_function = mz_se['function']
    # _usr_ms2_scan_id = mz_se['scan_id']
    # _usr_rt = mz_se['rt']
    # _usr_abbr_bulk = mz_se['Abbreviation']
    # _usr_pl_class = mz_se['Class']

    ms1_delta = lib_mz * ms1_precision

    ms1_pr_i = spec_info_dct['ms1_i']
    ms1_pr_mz = spec_info_dct['ms1_mz']
    ms1_rt = spec_info_dct['ms1_rt']
    ms2_rt = spec_info_dct['ms2_rt']
    ms1_df = spec_info_dct['ms1_df']
    ms2_df = spec_info_dct['ms2_df']

    ms_zoom_query_str = ' %.2f < mz < %.2f' % (ms1_obs - 2.1, ms1_obs + 2.1)
    ms_zoom_df = ms1_df.query(ms_zoom_query_str)

    print ('Start looking for MS2 PR m/z %f @ MS1 best PR m/z %f with lib m/z %f'
           % (ms2_pr_mz, ms1_obs, lib_mz))

    xic_df = xic_dct[ms1_obs]

    # Generate A4 image in landscape
    fig, pic_array = plt.subplots(nrows=3, ncols=2, figsize=(11.692, 8.267), sharex=False,
                                  sharey=False)
    # Make better spacing between subplots
    plt.tight_layout()
    xic_pic = pic_array[0, 0]
    msms_pic = pic_array[0, 1]
    ms_pic = pic_array[1, 0]
    msms_low_pic = pic_array[1, 1]
    ms_zoom_pic = pic_array[2, 0]
    msms_high_pic = pic_array[2, 1]

    xic_pic.tick_params(axis='both', which='major', labelsize=10)
    msms_pic.tick_params(axis='both', which='major', labelsize=10)
    ms_pic.tick_params(axis='both', which='major', labelsize=10)
    msms_low_pic.tick_params(axis='both', which='major', labelsize=10)
    ms_zoom_pic.tick_params(axis='both', which='major', labelsize=10)
    msms_high_pic.tick_params(axis='both', which='major', labelsize=10)

    # ms spectrum start
    ms_pic.stem(ms1_df['mz'].tolist(), ms1_df['i'].tolist(), 'grey', markerfmt=' ')

    dash_i = [ms1_df['i'].max()]
    markerline, stemlines, baseline = ms_pic.stem([ms1_pr_mz], dash_i, markerfmt=' ')
    plt.setp(stemlines, color='cyan', linewidth=5, alpha=0.3)
    markerline, stemlines, baseline = ms_pic.stem([ms1_pr_mz], [ms1_pr_i], markerfmt='D')
    plt.setp(stemlines, color='magenta')
    plt.setp(markerline, markerfacecolor='magenta', markersize=4, markeredgewidth=0)

    ms_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    ms_pic.set_xlabel("m/z", fontsize=10, labelpad=-1)
    ms_pic.set_ylabel("Intensity", fontsize=10)
    ms_pic.set_ylim([0, max(ms1_df['i'].tolist()) * 1.3])

    # add annotation
    _ms_pkl_top_df = ms1_df.sort_values(by='i', ascending=False).head(10)
    _ms_pkl_top_peak_list = zip(_ms_pkl_top_df['mz'].tolist(), _ms_pkl_top_df['i'].tolist())
    for _ms_pkl_top_peak in _ms_pkl_top_peak_list:
        _ms_pkl_top_peak_str = '%.4f' % _ms_pkl_top_peak[0]
        _ms_pkl_top_peak_y = _ms_pkl_top_peak[1]
        ms_pic.text(_ms_pkl_top_peak[0], _ms_pkl_top_peak_y, _ms_pkl_top_peak_str, fontsize=6)

    m0_theo_box = patches.Rectangle((lib_mz - ms1_delta, 0), 2 * ms1_delta, ms1_pr_i,
                                    facecolor='cyan', edgecolor="none")

    ms_zoom_pic.add_patch(m0_theo_box)

    # isotope region | if any peak in M-1.0034

    m_pre_theo_box = patches.Rectangle((lib_mz - 1.0034 - ms1_delta, 0), 2 * ms1_delta, ms1_pr_i,
                                       facecolor=(1.0, 0.0, 0.0, 0.4), edgecolor="none")
    ms_zoom_pic.add_patch(m_pre_theo_box)

    ms_zoom_pic.set_xlim([ms1_pr_mz - 1.5, ms1_pr_mz + 2.1])
    ms_zoom_pic.set_ylim([0, max(ms_zoom_df['i'].tolist()) * 1.3])
    ms_zoom_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0), fontsize=10)
    ms_zoom_pic.ticklabel_format(axis='x', useOffset=False, fontsize=10)
    ms_zoom_pic.set_xlabel('m/z', fontsize=10, labelpad=-1)
    ms_zoom_pic.set_ylabel('Intensity', fontsize=10)

    ms_zoom_pic.stem(ms_zoom_df['mz'].tolist(), ms_zoom_df['i'].tolist(), 'grey', markerfmt=' ')
    markerline, stemlines, baseline = ms_zoom_pic.stem([ms1_pr_mz], [ms1_pr_i],
                                                       'magenta', markerfmt='D'
                                                       )
    plt.setp(markerline, markerfacecolor='magenta', markeredgecolor='none', markeredgewidth=0,
             markersize=6, alpha=0.8)
    ms_zoom_pic.text(ms1_pr_mz, ms1_pr_i, '%.4f' % float(ms1_pr_mz),
                     color='magenta', fontsize=6
                     )
    markerline, stemlines, baseline = ms_zoom_pic.stem([lib_mz], [ms1_pr_i], '--', markerfmt='o')
    plt.setp(markerline, markerfacecolor='orange', markersize=6, markeredgewidth=0, alpha=0.9)
    plt.setp(stemlines, color='orange', alpha=0.8)
    ms_zoom_pic.text(lib_mz - 0.7, ms1_pr_i, 'Calc m/z: %.4f' % lib_mz, color='orange', fontsize=6)

    # isotope region | highlight the 1st isotope
    m1_dct = isotope_checker_dct[1]
    m1_theo_mz = m1_dct['theo_mz']
    # m1_theo_i = m1_dct['theo_i']
    m1_obs_mz = m1_dct['obs_mz']
    m1_obs_i = m1_dct['obs_i']
    m1_theo_r = m1_dct['theo_ratio']
    # m1_obs_r = m1_dct['obs_ratio']

    # theo range box
    m1_theo_box = patches.Rectangle((m1_theo_mz - ms1_delta, 0), 2 * ms1_delta, m1_theo_r * ms1_pr_i,
                                    facecolor=(0.1, 1.0, 1.0, 0.6), edgecolor="none")
    ms_zoom_pic.add_patch(m1_theo_box)

    markerline, stemlines, baseline = ms_zoom_pic.stem([m1_theo_mz], [m1_theo_r * ms1_pr_i], '--',
                                                       markerfmt='o')
    plt.setp(stemlines, color='orange', alpha=0.8)
    plt.setp(markerline, markerfacecolor='orange', markersize=6, markeredgewidth=0, alpha=0.9)
    ms_zoom_pic.text(m1_theo_mz - 0.93, m1_theo_r * ms1_pr_i, 'Calc 1st isotope: %.4f' % m1_theo_mz,
                     color='orange', fontsize=6)
    ms_zoom_pic.text(m1_obs_mz, m1_obs_i, '%.4f' % m1_obs_mz, color='magenta', fontsize=6)

    # isotope region | highlight the 2nd isotope
    m2_dct = isotope_checker_dct[2]
    m2_theo_mz = m2_dct['theo_mz']
    # m2_theo_i = m2_dct['theo_i']
    m2_obs_mz = m2_dct['obs_mz']
    m2_obs_i = m2_dct['obs_i']
    m2_theo_r = m2_dct['theo_ratio']
    # m2_obs_r = m2_dct['obs_ratio']
    m2_theo_box = patches.Rectangle((m2_theo_mz - ms1_delta, 0), 2 * ms1_delta, m2_theo_r * ms1_pr_i,
                                    facecolor=(0.2, 1.0, 1.0, 0.6), edgecolor="none")
    ms_zoom_pic.add_patch(m2_theo_box)
    markerline, stemlines, baseline = ms_zoom_pic.stem([m2_theo_mz], [m2_theo_r * ms1_pr_i], '--',
                                                       markerfmt='o')
    plt.setp(stemlines, color='orange', alpha=0.8)
    plt.setp(markerline, markerfacecolor='orange', markersize=6, markeredgewidth=0, alpha=0.9)
    ms_zoom_pic.text(m2_theo_mz - 0.93, m2_theo_r * ms1_pr_i, 'Calc 2nd isotope: %.4f' % m2_theo_mz,
                     color='orange', fontsize=6)
    ms_zoom_pic.text(m2_obs_mz, m2_obs_i, '%.4f' % m2_obs_mz, color='magenta', fontsize=6)
    ms_zoom_pic.text(m2_theo_mz, max(ms_zoom_df['i'].tolist()) * 1.2, 'Isotope score = % f' % isotope_score,
                     verticalalignment='top', horizontalalignment='right',
                     color='magenta', fontsize=8)

    # XIC spectrum start
    xic_rt_lst = xic_df['rt'].tolist()
    xic_i_lst = xic_df['i'].tolist()
    xic_pic.plot(xic_rt_lst, xic_i_lst, alpha=0.7, color='grey')
    xic_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    markerline, stemlines, baseline = xic_pic.stem([ms1_rt], [max(xic_i_lst)], markerfmt=' ')
    plt.setp(stemlines, color='magenta', linewidth=3, alpha=0.3)
    markerline, stemlines, baseline = xic_pic.stem([ms2_rt], [max(xic_i_lst)], '--', markerfmt=' ')
    plt.setp(stemlines, color='blue', linewidth=2, alpha=0.3)
    xic_pic.text(ms1_rt - 0.3, max(xic_i_lst) * 0.98, 'MS', fontsize=8, color='magenta')
    xic_pic.text(ms2_rt, max(xic_i_lst) * 0.98, 'MS/MS', fontsize=8, color='blue')
    xic_pic.set_xlabel("RT (min)", fontsize=10, labelpad=-1)
    xic_pic.set_ylabel("Intensity", fontsize=10)

    # prepare DataFrame for msms zoomed plot
    # plot color markers for zoomed MS2 first. Then overlay with zoomed spectra. Plot full ms2 in the last step.
    _msms_low_df = ms2_df.query('mz <= 350')
    _msms_high_df = ms2_df.query('mz > 350')

    _ident_table_df = ident_info_dct['SCORE_INFO']
    _fa_table_df = ident_info_dct['FA_INFO']
    _lyso_table_df = ident_info_dct['LYSO_INFO']

    if _ident_table_df.shape[0] > 0:
        _ident_table_df = _ident_table_df[['Proposed_structures', 'Score']]
        ident_col_labels = _ident_table_df.columns.values.tolist()
        ident_row_labels = _ident_table_df.index.tolist()
        ident_table_vals = map(list, _ident_table_df.values)
        # ident_col_width_lst = [0.03 * len(str(x)) for x in ident_col_labels]
        ident_col_width_lst = [0.3, 0.1]
        ident_table = ms_pic.table(cellText=ident_table_vals, rowLabels=ident_row_labels,
                                   colWidths=ident_col_width_lst,
                                   colLabels=ident_col_labels, loc='upper right')
        ident_table.set_fontsize(6)

    if _fa_table_df.shape[0] > 0:
        fa_row_color_lst = []
        for _i_fa, _fa_se in _fa_table_df.iterrows():
            # color of stemlines is tuple of R, G, B, alpha from 0 to 1
            _rgb_color = ((20 * _i_fa) / 255, (255 - 5 * _i_fa) / 255,
                          (255 - 5 * _i_fa) / 255, 0.6 - 0.02 * _i_fa
                          )
            markerline, stemlines, baseline = msms_low_pic.stem([_fa_se['mz']], [_fa_se['i']], markerfmt='D')
            plt.setp(stemlines, color=_rgb_color, linewidth=3)
            plt.setp(markerline, markerfacecolor=_rgb_color, markersize=5, markeredgewidth=0)
            markerline, stemlines, baseline = msms_pic.stem([_fa_se['mz']], [_fa_se['i']], markerfmt='D')
            plt.setp(stemlines, color=_rgb_color, linewidth=3)
            plt.setp(markerline, markerfacecolor=_rgb_color, markersize=5, markeredgewidth=0)
            _msms_low_peak_str = '%.4f' % _fa_se['mz']
            _msms_low_peak_y = float(_fa_se['i'])
            # can set rotation = 90 or 0
            msms_low_pic.text(_fa_se['mz'], _msms_low_peak_y, _msms_low_peak_str, fontsize=6)

            fa_row_color_lst.append(_rgb_color)
        fa_col_labels = _fa_table_df.columns.values.tolist()
        fa_row_labels = _fa_table_df.index.tolist()
        fa_table_vals = map(list, _fa_table_df.values)
        # fa_col_width_lst = [0.025 * len(str(x)) for x in fa_col_labels]
        fa_col_width_lst = [0.3, 0.1, 0.1, 0.1]
        fa_table = msms_pic.table(cellText=fa_table_vals, rowLabels=fa_row_labels,
                                  colWidths=fa_col_width_lst, rowColours=fa_row_color_lst,
                                  colLabels=fa_col_labels, loc='upper center')
        fa_table.set_fontsize(6)

    if _lyso_table_df.shape[0] > 0:
        lyso_row_color_lst = []
        for _i_lyso, _lyso_se in _lyso_table_df.iterrows():
            # color of stemlines is tuple of R, G, B, alpha from 0 to 1
            _rgb_color = ((20 * _i_lyso) / 255, (255 - 5 * _i_lyso) / 255,
                          (255 - 5 * _i_lyso) / 255, 0.6 - 0.02 * _i_lyso
                          )
            markerline, stemlines, baseline = msms_high_pic.stem([_lyso_se['mz']], [_lyso_se['i']],
                                                                 markerfmt="D"
                                                                 )
            plt.setp(stemlines, color=_rgb_color, linewidth=3)
            plt.setp(markerline, markerfacecolor=_rgb_color, markersize=5, markeredgewidth=0)
            markerline, stemlines, baseline = msms_pic.stem([_lyso_se['mz']], [_lyso_se['i']], markerfmt="D")
            plt.setp(stemlines, color=_rgb_color, linewidth=3)
            plt.setp(markerline, markerfacecolor=_rgb_color, markersize=5, markeredgewidth=0)
            _msms_high_peak_str = '%.4f' % _lyso_se['mz']
            _msms_high_peak_y = float(_lyso_se['i'])
            msms_high_pic.text(_lyso_se['mz'], _msms_high_peak_y, _msms_high_peak_str, fontsize=6)

            lyso_row_color_lst.append(_rgb_color)

        lyso_col_labels = _lyso_table_df.columns.values.tolist()
        lyso_row_labels = _lyso_table_df.index.tolist()
        lyso_table_vals = map(list, _lyso_table_df.values)
        # lyso_col_width_lst = [0.01 * len(str(x)) for x in lyso_col_labels]
        lyso_col_width_lst = [0.3, 0.1, 0.1, 0.1]

        lyso_table = msms_high_pic.table(cellText=lyso_table_vals, rowLabels=lyso_row_labels,
                                         colWidths=lyso_col_width_lst, rowColours=lyso_row_color_lst,
                                         colLabels=lyso_col_labels, loc='upper center')
        lyso_table.set_fontsize(6)

    # msms spectrum zoomed < 350 start
    if _msms_low_df.shape[0] > 0:
        msms_low_pic.stem(_msms_low_df['mz'].tolist(),
                          _msms_low_df['i'].tolist(),
                          'black', lw=2, markerfmt=' ')
        msms_low_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        msms_low_pic.set_xlabel("m/z", fontsize=10, labelpad=-1)
        msms_low_pic.set_ylabel("Intensity", fontsize=10)
        msms_low_pic.set_xlim([min(_msms_low_df['mz'].tolist()) - 1, 350])
        msms_low_pic.set_ylim([0, max(_msms_low_df['i'].tolist()) * 1.3])
    else:
        pass

    # msms spectrum zoomed > 350 start
    if _msms_high_df.shape[0] > 0:
        msms_high_pic.stem(_msms_high_df['mz'].tolist(),
                           _msms_high_df['i'].tolist(),
                           'black', lw=4, markerfmt=' ')
        msms_high_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        msms_high_pic.set_xlabel("m/z", fontsize=10, labelpad=-1)
        msms_high_pic.set_ylabel("Intensity", fontsize=10)
        msms_high_pic.set_xlim([350, ms2_pr_mz + 20])
        msms_high_pic.set_ylim([0, max(_msms_high_df['i'].tolist()) * 1.3])

        # add annotations
        _top_msms_high_df = _msms_high_df.sort_values(by='i', ascending=False)
        _top_msms_high_df = _top_msms_high_df.head(10)
        _msms_high_peak_list = zip(_top_msms_high_df['mz'].tolist(),
                                   _top_msms_high_df['i'].tolist())
        for _msms_high_peak in _msms_high_peak_list:
            _msms_high_peak_str = '%.4f' % _msms_high_peak[0]
            _msms_high_peak_y = _msms_high_peak[1]
            msms_high_pic.text(_msms_high_peak[0], _msms_high_peak_y, _msms_high_peak_str, fontsize=6)
    else:
        pass

    # add specific ion info

    if 'OTHER_FRAG' in specific_check_dct.keys():
        other_frag_df = specific_check_dct['OTHER_FRAG']
        for _idx, _frag_se in other_frag_df.iterrows():
            _frag_mz = _frag_se['mz']
            _frag_i = _frag_se['i']
            _frag_class = _frag_se['LABEL']
            _frag_i_x = min(_frag_i * 5, 0.3 * max(_msms_low_df['i'].tolist()))
            _frag_i = sorted([max(_msms_low_df['i'].tolist()) * 1.1, _frag_i * 1.1, _frag_i_x])[1]
            markerline, stemlines, baseline = msms_low_pic.stem([_frag_mz], [_frag_i], markerfmt=' ')
            plt.setp(stemlines, color='red', linewidth=3, alpha=0.4)
            markerline, stemlines, baseline = msms_pic.stem([_frag_mz], [_frag_i], markerfmt=' ')
            plt.setp(stemlines, color='red', linewidth=3, alpha=0.4)
            msms_low_pic.text(_frag_mz, _frag_i, _frag_class, fontsize=8, color='red')

    if 'OTHER_NL' in specific_check_dct.keys():
        other_nl_df = specific_check_dct['OTHER_NL']
        for _idx, _nl_se in other_nl_df.iterrows():
            _nl_mz = _nl_se['mz']
            _nl_i = _nl_se['i']
            _nl_class = _nl_se['LABEL']
            _nl_i_x = min(_nl_i * 5, 0.3 * max(_msms_high_df['i'].tolist()))
            _nl_i = sorted([max(_msms_high_df['i'].tolist()) * 1.1, _nl_i * 1.1, _nl_i_x])[1]
            markerline, stemlines, baseline = msms_high_pic.stem([_nl_mz], [_nl_i], markerfmt=' ')
            plt.setp(stemlines, color='red', linewidth=3, alpha=0.4)
            markerline, stemlines, baseline = msms_pic.stem([_nl_mz], [_nl_i], markerfmt=' ')
            plt.setp(stemlines, color='red', linewidth=3, alpha=0.4)
            msms_high_pic.text(_nl_mz, _nl_i, _nl_class, fontsize=8, color='red')

    if 'TARGET_FRAG' in specific_check_dct.keys():
        target_frag_df = specific_check_dct['TARGET_FRAG']
        for _idx, _frag_se in target_frag_df.iterrows():
            _frag_mz = _frag_se['mz']
            _frag_i = _frag_se['i']
            _frag_class = _frag_se['LABEL']
            _frag_i_x = min(_frag_i * 5, 0.5 * max(_msms_low_df['i'].tolist()))
            _frag_i = sorted([max(_msms_low_df['i'].tolist()) * 1.1, _frag_i * 1.1, _frag_i_x])[1]
            markerline, stemlines, baseline = msms_low_pic.stem([_frag_mz], [_frag_i], markerfmt=' ')
            plt.setp(stemlines, color='green', linewidth=3, alpha=0.4)
            markerline, stemlines, baseline = msms_pic.stem([_frag_mz], [_frag_i], markerfmt=' ')
            plt.setp(stemlines, color='green', linewidth=3, alpha=0.4)
            msms_low_pic.text(_frag_mz, _frag_i, _frag_class, fontsize=8, color='green')

    if 'TARGET_NL' in specific_check_dct.keys():
        target_nl_df = specific_check_dct['TARGET_NL']
        for _idx, _nl_se in target_nl_df.iterrows():
            _nl_mz = _nl_se['mz']
            _nl_i = _nl_se['i']
            _nl_class = _nl_se['LABEL']
            _nl_i_x = min(_nl_i * 5, 0.3 * max(_msms_low_df['i'].tolist()))
            _nl_i = sorted([max(_msms_high_df['i'].tolist()) * 1.1, _nl_i * 1.1, _nl_i_x])[1]
            markerline, stemlines, baseline = msms_high_pic.stem([_nl_mz], [_nl_i], markerfmt=' ')
            plt.setp(stemlines, color='green', linewidth=3, alpha=0.4)
            markerline, stemlines, baseline = msms_pic.stem([_nl_mz], [_nl_i], markerfmt=' ')
            plt.setp(stemlines, color='green', linewidth=3, alpha=0.4)
            msms_high_pic.text(_nl_mz, _nl_i, _nl_class, fontsize=8, color='green')

    # msms spectrum start
    msms_pic.stem(ms2_df['mz'].tolist(), ms2_df['i'].tolist(), 'black', markerfmt=' ')
    msms_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    msms_pic.set_xlabel("m/z", fontsize=10, labelpad=-1)
    msms_pic.set_ylabel("Intensity", fontsize=10)
    msms_pic.set_xlim([min(ms2_df['mz'].tolist()) - 1, ms2_pr_mz + 20])
    msms_pic.set_ylim([0, max(ms2_df['i'].tolist()) * 1.3])

    # set title
    xic_title_str = 'XIC of m/z %.4f | %s @ m/z %.4f ppm=%.2f' % (ms1_pr_mz, abbr_id, lib_mz, ms1_pr_ppm)
    ms_title_str = 'MS @ %.3f min' % ms1_rt
    ms_zoom_title_str = 'MS zoomed'
    msms_title_str = ('MS/MS of m/z %.4f | DDA Top %d @ %.3f min' % (ms2_pr_mz, func_id, ms2_rt))
    msms_low_str = 'MS/MS zoomed below m/z 350'
    msms_high_str = 'MS/MS zoomed above m/z 350'

    xic_pic.set_title(xic_title_str, color='b', fontsize=10, y=0.98)
    ms_pic.set_title(ms_title_str, color='b', fontsize=10, y=0.98)
    ms_zoom_pic.set_title(ms_zoom_title_str, color='b', fontsize=10, y=0.98)
    msms_pic.set_title(msms_title_str, color='b', fontsize=10, y=0.98)
    msms_low_pic.set_title(msms_low_str, color='b', fontsize=10, y=0.98)
    msms_high_pic.set_title(msms_high_str, color='b', fontsize=10, y=0.98)

    print ('>>> >>> >>> try to plot >>> >>> >>>')

    plt.savefig(save_img_as, dpi=300)
    print ('=====> Image saved as: %s' % save_img_as)
    plt.close()
    isotope_checker = 0
    return isotope_checker, isotope_score
