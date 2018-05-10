# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2017  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
# LipidHunter is Dual-licensed
#     For academic and non-commercial use: `GPLv2 License` Please read more information by the following link:
#         [The GNU General Public License version 2] (https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
#     For commercial use:
#         please contact the SysMedOs_team by email.
# Please cite our publication in an appropriate form.
# Ni, Zhixu, Georgia Angelidou, Mike Lange, Ralf Hoffmann, and Maria Fedorova.
# "LipidHunter identifies phospholipids by high-throughput processing of LC-MS and shotgun lipidomics datasets."
# Analytical Chemistry (2017).
# DOI: 10.1021/acs.analchem.7b01126
#
# For more info please contact:
#     SysMedOs_team: oxlpp@bbz.uni-leipzig.de
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de
#     Developer Georgia Angelidou georgia.angelidou@uni-leipzig.de

from __future__ import division
from __future__ import print_function

import io
import pandas as pd
from PIL import Image
import time
import matplotlib

matplotlib.use('agg')
from matplotlib import pyplot as plt
import matplotlib.patches as patches
# import matplotlib as mpl
# from matplotlib.offsetbox import AnnotationBbox, OffsetImage
# from matplotlib._png import read_png

from concurrent.futures import ThreadPoolExecutor


def plot_spectra(abbr, mz_se, xic_dct, ident_info_dct, spec_info_dct, isotope_score_info_dct, specific_dct,
                 formula_charged, charge, core_count, save_img_as=None, img_type='png', dpi=300, vendor='waters',
                 ms1_precision=50e-6):
    ms2_pr_mz = mz_se['MS2_PR_mz']
    ms1_obs = mz_se['MS1_obs_mz']
    ms1_xic_mz = mz_se['MS1_XIC_mz']
    lib_mz = mz_se['Lib_mz']
    func_id = mz_se['DDA_rank']
    ms1_pr_ppm = mz_se['ppm']

    obs_info_df = ident_info_dct['INFO']
    obs_fa_df = ident_info_dct['OBS_FA']
    obs_lyso_df = ident_info_dct['OBS_LYSO']
    obs_ident_df = ident_info_dct['IDENT']

    isotope_score = isotope_score_info_dct['isotope_score']
    isotope_checker_dct = isotope_score_info_dct['isotope_checker_dct']
    m2_score = isotope_score_info_dct['m2_score']
    m2_checker_dct = isotope_score_info_dct['m2_checker_dct']
    deconv_lst = isotope_score_info_dct['deconv_lst']

    print(core_count, '>>> Start to plot %s -> MS2 PR m/z %.4f @ MS1 best PR m/z %.4f with lib m/z %.4f'
          % (abbr, ms2_pr_mz, ms1_obs, lib_mz))

    if len(deconv_lst) == 4:
        pass
    else:
        deconv_lst = [0, 0, 0, 0]

    ms1_delta = lib_mz * ms1_precision

    ms1_pr_i = spec_info_dct['ms1_i']
    ms1_pr_mz = spec_info_dct['ms1_mz']
    ms1_rt = spec_info_dct['ms1_rt']
    ms2_rt = spec_info_dct['ms2_rt']
    ms1_df = spec_info_dct['ms1_df']
    ms2_df = spec_info_dct['ms2_df']

    ms1_df = ms1_df.sort_values(by='mz', ascending='True')
    dash_i = [ms1_df['i'].max()]

    ms_zoom_query_str = ' %.2f < mz < %.2f' % (ms1_obs - 1.5, ms1_obs + 3.55)
    ms_zoom_df = ms1_df.query(ms_zoom_query_str)
    try:
        ms_zoom_bp_i = max(ms_zoom_df['i'].values.tolist())
    except ValueError:
        ms_zoom_bp_i = 0

    xic_df = xic_dct[ms1_xic_mz]

    xic_rt_lst = xic_df['rt'].values.tolist()
    xic_i_lst = xic_df['i'].values.tolist()

    # if ms_zoom_bp_i > 0 and len(xic_rt_lst) > 0 and len(xic_i_lst) > 0:

    # cut lower peaks to accelerate plotting time
    m1_dct = isotope_checker_dct[1]
    m1_theo_mz = m1_dct['theo_mz']
    m1_theo_i = m1_dct['theo_i']
    m1_obs_mz = m1_dct['obs_mz']
    m1_obs_i = m1_dct['obs_i']

    if ms1_df['i'].max() >= 10000 and ms1_df.shape[0] >= 500:
        ms1_min = ms1_df['i'].min()
        ms1_max = ms1_df['i'].max()
        ms1_top1000_i = sorted(ms1_df['i'].values.tolist(), reverse=True)[499]
        ms1_plot_th = min(m1_obs_i, 3 * ms1_min, ms1_max * 0.01, 1000, ms1_top1000_i)
        ms1_plot_th = max(ms1_plot_th, ms1_top1000_i)
        # print(core_count, m1_obs_i, 3 * ms1_min, ms1_max * 0.01, 1000, ms1_top1000_i)
        ms1_df = ms1_df.query('i >= %f' % ms1_plot_th)
        print(core_count, 'Plot full MS1 with abs intensity filter > %f' % ms1_plot_th)
    if ms2_df['i'].max() >= 1000 and ms2_df.shape[0] >= 500:
        ms2_min = ms2_df['i'].min()
        ms2_max = ms2_df['i'].max()

        ms2_top1000_i = sorted(ms2_df['i'].values.tolist(), reverse=True)[499]
        ms2_min_lst = [3 * ms2_min, ms2_max * 0.01, 10, ms2_top1000_i]
        ms2_plot_th = max(min(ms2_min_lst), ms2_top1000_i)

        # print(core_count, ms2_min_lst)
        ms2_plot_th -= 1
        if ms2_plot_th > 0:
            ms2_df = ms2_df.query('i >= %f' % ms2_plot_th)
            print(core_count, 'Plot full MS/MS with abs intensity filter > %f' % ms2_plot_th)

    _msms_low_df = ms2_df.query('mz <= 400')
    _msms_high_df = ms2_df.query('mz > 400')
    _msms_high_df = _msms_high_df.query('mz < %.4f' % (ms2_pr_mz + 1))

    # Generate A4 image in landscape
    fig, pic_array = plt.subplots(nrows=3, ncols=2, figsize=(11.692, 8.267), sharex='none', sharey='none')
    # Make better spacing between subplots
    plt.tight_layout()

    _msms_max = ms2_df['i'].max()

    # remove identified FA
    ident_peak_lst = obs_ident_df['obs_abbr'].values.tolist()
    ident_peak_lst.sort()
    ident_peak_lst = set(ident_peak_lst)
    frag_idx_lst = []
    mg_idx_lst = []
    dg_Na_idx_lst = []
    nl_idx_lst = []
    obs_mg_df = obs_fa_df.loc[obs_fa_df['TYPE'] == 'MG']
    obs_mg_df.reset_index(drop=True, inplace=True)
    # Need this because we need to keep their original position in the DF to remove them later
    obs_mg_idx = obs_fa_df.index[obs_fa_df['TYPE'] == 'MG'].tolist()
    obs_dg_na_df = obs_lyso_df.loc[obs_lyso_df['TYPE'] == 'NL_Na']
    obs_dg_na_df.reset_index(drop=True, inplace=True)
    obs_dg_na_idx = obs_lyso_df.index[obs_lyso_df['TYPE'] == 'NL_Na'].tolist()

    if isinstance(obs_fa_df, pd.DataFrame):
        for _idx, _ion_se in obs_fa_df.iterrows():
            # if _ion_se['obs_abbr'] in ident_peak_lst:
            if _ion_se['obs_abbr'] in ident_peak_lst and _ion_se['TYPE'] == 'FA':
                frag_idx_lst.append(_idx)
    if isinstance(obs_mg_df, pd.DataFrame):
        for _idx, _ion_se in obs_mg_df.iterrows():
            if _ion_se['obs_abbr'] in ident_peak_lst:
                mg_idx_lst.append(_idx)
    if isinstance(obs_lyso_df, pd.DataFrame):
        for _idx, _ion_se in obs_lyso_df.iterrows():
            # if _ion_se['obs_abbr'] in ident_peak_lst:
            if _ion_se['obs_abbr'] in ident_peak_lst and _ion_se['TYPE'] == 'NL':
                nl_idx_lst.append(_idx)
    if isinstance(obs_dg_na_df, pd.DataFrame):
        for _idx, _ion_se in obs_dg_na_df.iterrows():
            if _ion_se['obs_abbr'] in ident_peak_lst and _ion_se['TYPE'] == 'NL_Na':
                dg_Na_idx_lst.append(_idx)

    # from the below we remove all of the MG for TG because there is another table for them
    plt_obs_fa_df = pd.DataFrame(obs_fa_df.drop(frag_idx_lst + obs_mg_idx))
    plt_obs_mg_df = pd.DataFrame(obs_mg_df.drop(mg_idx_lst))
    plt_obs_lyso_df = pd.DataFrame(obs_lyso_df.drop(nl_idx_lst + obs_dg_na_idx))
    plt_obs_dg_na_df = pd.DataFrame(obs_dg_na_df.drop(dg_Na_idx_lst))

    # add specific ion info
    txt_props = {'ha': 'left', 'va': 'bottom'}
    obs_ident_df['i_r'] = obs_ident_df['i'] * 1.025
    if not plt_obs_fa_df.empty:
        plt_obs_fa_df['i_r'] = obs_fa_df['i'] * 1.025
    else:
        plt_obs_fa_df = False
    if not plt_obs_mg_df.empty:
        plt_obs_mg_df['i_r'] = obs_mg_df['i'] * 1.025
    else:
        plt_obs_mg_df = False
    if not plt_obs_lyso_df.empty:
        plt_obs_lyso_df['i_r'] = obs_lyso_df['i'] * 1.075
    else:
        plt_obs_lyso_df = False
    if not plt_obs_dg_na_df.empty:
        plt_obs_dg_na_df['i_r'] = obs_dg_na_df['i'] * 1.075
    else:
        plt_obs_dg_na_df = False

    if 'OTHER_FRAG' in list(specific_dct.keys()):
        other_frag_df = specific_dct['OTHER_FRAG']
        other_frag_df['i_r'] = other_frag_df['i'] * 1.4
    else:
        other_frag_df = False
    if 'OTHER_NL' in list(specific_dct.keys()):
        other_nl_df = specific_dct['OTHER_NL']
        other_nl_df['i_r'] = other_nl_df['i'] * 1.2
    else:
        other_nl_df = False
    if 'TARGET_FRAG' in list(specific_dct.keys()):
        target_frag_df = specific_dct['TARGET_FRAG']
        target_frag_df['i_r'] = target_frag_df['i'] * 1.4
    else:
        target_frag_df = False
    if 'TARGET_NL' in list(specific_dct.keys()):
        target_nl_df = specific_dct['TARGET_NL']
        target_nl_df['i_r'] = target_nl_df['i'] * 1.2
    else:
        target_nl_df = False

    def plot_xic():
        # _t_img_0 = time.time()
        # print(core_count, 'start to plot XIC ...')
        xic_pic = pic_array[0, 0]
        xic_pic.tick_params(axis='both', which='major', labelsize=10)
        xic_rt_min = min(xic_rt_lst)
        xic_rt_max = max(xic_rt_lst)
        xic_rt_label_shift = (xic_rt_max - xic_rt_min) * 0.04
        xic_pic.plot(xic_rt_lst, xic_i_lst, alpha=0.7, color='grey')
        xic_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        _marker_line, _stem_lines, _base_line = xic_pic.stem([ms1_rt], [max(xic_i_lst)], markerfmt=' ')
        plt.setp(_stem_lines, color=(0.0, 0.5, 0.9, 0.7), linewidth=3, alpha=0.3)
        _marker_line, _stem_lines, _base_line = xic_pic.stem([ms2_rt], [max(xic_i_lst)], '--', markerfmt=' ')
        plt.setp(_stem_lines, color='blue', linewidth=2, alpha=0.3)
        xic_pic.text(ms1_rt - xic_rt_label_shift, max(xic_i_lst) * 0.98, 'MS', fontsize=8, color=(0.0, 0.5, 1.0, 1.0))
        xic_pic.text(ms2_rt, max(xic_i_lst) * 0.98, 'MS/MS', fontsize=8, color='blue')
        xic_pic.set_xlabel("Scan time (min)", fontsize=10, labelpad=-1)
        xic_pic.set_ylabel("Intensity", fontsize=10)
        xic_pic.set_xlim([xic_rt_min, xic_rt_max])
        xic_pic.set_ylim([0, max(xic_i_lst) * 1.1])
        xic_title_str = 'XIC of m/z %.4f | @ m/z %.4f ppm=%.2f' % (ms1_pr_mz, lib_mz, ms1_pr_ppm)
        xic_pic.set_title(xic_title_str, color='b', fontsize=8, y=0.98)
        # print(core_count, 'plot XIC in ', time.time() - _t_img_0)

    def plot_ms():
        # _t_img_0 = time.time()
        # print(core_count, 'start to plot MS ...')
        ms_pic = pic_array[1, 0]
        ms_pic.tick_params(axis='both', which='major', labelsize=10)

        if ms1_df.shape[0] > 800:
            ms_pic.plot(ms1_df['mz'].values.tolist(), ms1_df['i'].values.tolist(), 'grey', lw=1)
        else:
            ms_pic.stem(ms1_df['mz'].values.tolist(), ms1_df['i'].values.tolist(), 'grey', markerfmt=' ')

        _marker_line, _stem_lines, _base_line = ms_pic.stem([ms1_pr_mz], dash_i, markerfmt=' ')
        plt.setp(_stem_lines, color='#00ccff', linewidth=5, alpha=0.15)
        _marker_line, _stem_lines, _base_line = ms_pic.stem([ms1_pr_mz], [ms1_pr_i], markerfmt='D')
        plt.setp(_stem_lines, color=(0.0, 0.5, 0.9, 1.0))
        plt.setp(_marker_line, markerfacecolor=(0.0, 0.5, 0.9, 1.0), markersize=4, markeredgewidth=0)

        ms_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        ms_pic.set_xlabel("m/z", fontsize=10, labelpad=-1)
        ms_pic.set_ylabel("Intensity", fontsize=10)
        ms_pic.set_ylim([0, max(ms1_df['i'].values.tolist()) * 1.3])

        # add annotation
        _ms_pkl_top_df = ms1_df.sort_values(by='i', ascending=False).head(10)
        _ms_pkl_top_peak_list = list(zip(_ms_pkl_top_df['mz'].values.tolist(), _ms_pkl_top_df['i'].values.tolist()))
        for _ms_pkl_top_peak in _ms_pkl_top_peak_list:
            _ms_pkl_top_peak_str = '%.4f' % _ms_pkl_top_peak[0]
            _ms_pkl_top_peak_y = _ms_pkl_top_peak[1]
            ms_pic.text(_ms_pkl_top_peak[0], _ms_pkl_top_peak_y, _ms_pkl_top_peak_str, fontsize=6)

        ms_title_str = 'MS @ %.3f min | %s' % (ms1_rt, abbr)
        ms_pic.set_title(ms_title_str, color='b', fontsize=8, y=0.98)
        # print(core_count, 'plot MS in ', time.time() - _t_img_0)

    def plot_ms_zoom():
        # _t_img_0 = time.time()
        # print(core_count, 'start to plot MS zoom ...')
        ms_zoom_pic = pic_array[2, 0]
        ms_zoom_pic.tick_params(axis='both', which='major', labelsize=10)
        # isotope region | if any peak in M-1.0034
        m_pre_theo_box = patches.Rectangle((lib_mz - 1.0034 - ms1_delta, 0), 2 * ms1_delta, ms1_pr_i,
                                           facecolor=(1.0, 0.0, 0.0, 0.6), edgecolor='none')
        ms_zoom_pic.add_patch(m_pre_theo_box)
        ms_zoom_offset_i = ms_zoom_bp_i * 0.1

        ms_zoom_pic.set_xlim([ms1_pr_mz - 1.5, ms1_pr_mz + 3.55])
        ms_zoom_pic.set_ylim([0, ms_zoom_bp_i * 1.45])
        ms_zoom_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0), fontsize=10)
        ms_zoom_pic.ticklabel_format(axis='x', useOffset=False, fontsize=10)
        ms_zoom_pic.set_xlabel('m/z', fontsize=10, labelpad=-1)
        ms_zoom_pic.set_ylabel('Intensity', fontsize=10)

        plt_ms_zoom_df = ms_zoom_df.sort_values(by='mz', ascending='True')
        if vendor == 'waters':
            ms_zoom_pic.plot(plt_ms_zoom_df['mz'].values.tolist(), plt_ms_zoom_df['i'].values.tolist(),
                             'grey', lw=1, zorder=1)
        elif vendor == 'thermo':
            ms_zoom_pic.stem(plt_ms_zoom_df['mz'].values.tolist(), plt_ms_zoom_df['i'].values.tolist(),
                             'grey', markerfmt=' ')  # zorder=1
        else:
            ms_zoom_pic.plot(plt_ms_zoom_df['mz'].values.tolist(), plt_ms_zoom_df['i'].values.tolist(),
                             'grey', lw=1, zorder=1)

        _marker_l, _stem_l, _base_l = ms_zoom_pic.stem([ms1_pr_mz], [ms1_pr_i],
                                                       color=(0.0, 0.5, 0.9, 1.0), markerfmt='D')  # zorder=20
        plt.setp(_marker_l, markerfacecolor=(0.0, 0.5, 1.0, 1.0), markeredgecolor='none', markeredgewidth=0,
                 markersize=6, alpha=0.8)
        ms_zoom_pic.text(ms1_pr_mz + 0.06, ms1_pr_i, '%.4f' % float(ms1_pr_mz),
                         color=(0.0, 0.5, 1.0, 1.0), fontsize=6)
        _marker_l, _stem_l, _base_l = ms_zoom_pic.stem([lib_mz], [ms1_pr_i], '--', markerfmt='o')  # zorder=21
        plt.setp(_marker_l, markerfacecolor='#ff6600', markersize=6, markeredgewidth=0)
        plt.setp(_stem_l, color='#ff6600')
        ms_zoom_pic.text(lib_mz - 0.15, ms1_pr_i + ms_zoom_offset_i, '[M+0]', color='#ff6600', fontsize=6)
        ms_zoom_pic.text(lib_mz - 0.71, ms1_pr_i, 'Calc: %.4f' % lib_mz, color='#ff6600', fontsize=6)

        # isotope region | highlight the 1st isotope
        # theo range box
        m1_theo_base_box = patches.Rectangle((m1_theo_mz - ms1_delta, 0), 2 * ms1_delta, deconv_lst[1],
                                             facecolor=(1.0, 0.0, 0.0, 0.6), edgecolor='none', zorder=1)
        ms_zoom_pic.add_patch(m1_theo_base_box)
        m1_theo_box = patches.Rectangle((m1_theo_mz - ms1_delta, deconv_lst[1]), 2 * ms1_delta,
                                        m1_theo_i - deconv_lst[1],
                                        facecolor=(0, 0.8, 1.0, 0.6), edgecolor='none', zorder=1)
        ms_zoom_pic.add_patch(m1_theo_box)

        _marker_l, _stem_l, _base_l = ms_zoom_pic.stem([m1_theo_mz], [m1_theo_i], '--', markerfmt='o') # zorder=22
        plt.setp(_stem_l, color='#ff6600')
        plt.setp(_marker_l, markerfacecolor='#ff6600', markersize=6, markeredgewidth=0)
        ms_zoom_pic.text(m1_theo_mz - 0.15, m1_theo_i + ms_zoom_offset_i, '[M+1]', color='#ff6600', fontsize=6)
        ms_zoom_pic.text(m1_theo_mz - 0.71, m1_theo_i, 'Calc: %.4f' % m1_theo_mz,
                         color='#ff6600', fontsize=6)
        ms_zoom_pic.text(m1_obs_mz + 0.04, m1_obs_i, '%.4f' % m1_obs_mz, color=(0.0, 0.5, 0.9, 1.0), fontsize=6)

        opt_box_lst = []

        # isotope region | highlight the 2nd isotope
        if 2 in list(isotope_checker_dct.keys()):
            m2_dct = isotope_checker_dct[2]
            m2_theo_mz = m2_dct['theo_mz']
            m2_theo_i = m2_dct['theo_i']
            m2_obs_mz = m2_dct['obs_mz']
            m2_obs_i = m2_dct['obs_i']
            # m2_theo_r = m2_dct['theo_ratio']
            # # m2_obs_r = m2_dct['obs_ratio']
            # m2_theo_box = patches.Rectangle((m2_theo_mz - ms1_delta, 0), 2 * ms1_delta, m2_theo_i,
            #                                 facecolor=(0.2, 1.0, 1.0, 0.3), edgecolor='none', zorder=9)
            m2_theo_box = patches.Rectangle((m2_theo_mz - ms1_delta, 0), 2 * ms1_delta, m2_theo_i,
                                            facecolor=(0, 0.8, 1.0, 0.6), edgecolor='none', zorder=1)
            ms_zoom_pic.add_patch(m2_theo_box)
            opt_box_lst.append(ms_zoom_pic)
            _marker_l, _stem_l, _base_l = ms_zoom_pic.stem([m2_theo_mz], [m2_theo_i], '--', markerfmt='o')  # zorder=23
            plt.setp(_stem_l, color='#ff6600', alpha=0.8)
            plt.setp(_marker_l, markerfacecolor='#ff6600', markersize=6, markeredgewidth=0, alpha=0.9)
            ms_zoom_pic.text(m2_theo_mz - 0.15, m2_theo_i + ms_zoom_offset_i, '[M+2]', color='#ff6600', fontsize=6)
            plt.setp(_marker_l, markerfacecolor='#ff6600', markersize=6, markeredgewidth=0, alpha=0.9)
            ms_zoom_pic.text(m2_theo_mz - 0.71, m2_theo_i, 'Calc: %.4f' % m2_theo_mz,
                             color='#ff6600', fontsize=6)
            ms_zoom_pic.text(m2_obs_mz + 0.04, m2_obs_i, '%.4f' % m2_obs_mz, color=(0.0, 0.5, 0.9, 1.0), fontsize=6)

        if len(list(m2_checker_dct.keys())) > 0:
            for _mh2 in list(m2_checker_dct.keys()):
                mh2_dct = m2_checker_dct[_mh2]
                mh2_theo_mz = mh2_dct['theo_mz']
                mh2_theo_i = mh2_dct['theo_i']
                mh2_obs_mz = mh2_dct['obs_mz']
                mh2_obs_i = mh2_dct['obs_i']
                decon_idx = _mh2 + 2
                # mh2_theo_r = mh2_dct['theo_ratio']
                # mh2_obs_r = mh2_dct['obs_ratio']
                mh2_theo_base_box = patches.Rectangle((mh2_theo_mz - ms1_delta, 0), 2 * ms1_delta,
                                                      deconv_lst[decon_idx],
                                                      facecolor=(0, 0.8, 1.0, 0.6), edgecolor='none', zorder=1)
                ms_zoom_pic.add_patch(mh2_theo_base_box)
                # opt_box_lst.append(mh2_theo_base_box)
                mh2_theo_box = patches.Rectangle((mh2_theo_mz - ms1_delta, deconv_lst[decon_idx]),
                                                 2 * ms1_delta, mh2_theo_i - deconv_lst[decon_idx],
                                                 facecolor=(1.0, 0.0, 0.0, 0.6), edgecolor='none', zorder=1)
                ms_zoom_pic.add_patch(mh2_theo_box)
                _marker_l, _stem_l, _base_l = ms_zoom_pic.stem([mh2_theo_mz], [mh2_theo_i], '--',
                                                               markerfmt='o')  # zorder=24
                plt.setp(_stem_l, color='red', alpha=0.8)
                plt.setp(_marker_l, markerfacecolor='red', markersize=6, markeredgewidth=0, alpha=0.9)
                if _mh2 == 0:
                    _mh2_name = ''
                else:
                    _mh2_name = '+%i' % _mh2
                ms_zoom_pic.text(mh2_theo_mz - 0.2 - 0.05 * _mh2, mh2_theo_i + ms_zoom_offset_i,
                                 '[M+2H%s]' % _mh2_name, color='red', fontsize=6)
                ms_zoom_pic.text(mh2_theo_mz - 0.71, mh2_theo_i, 'Calc: %.4f' % mh2_theo_mz,
                                 color='red', fontsize=6)
                ms_zoom_pic.text(mh2_obs_mz + 0.04, mh2_obs_i, '%.4f' % mh2_obs_mz, color='red', fontsize=6)

            # plot the M+2H isotope score

            ms_zoom_pic.text(m1_theo_mz + 2.5, max(ms_zoom_bp_i - ms_zoom_offset_i, ms1_pr_i * 0.8),
                             '[M+2H] Isotope score = %.1f' % m2_score,
                             verticalalignment='top', horizontalalignment='right',
                             color='red', fontsize=7)

        # plot the isotope score
        ms_zoom_pic.text(m1_theo_mz + 0.2, ms_zoom_bp_i + 3 * ms_zoom_offset_i,
                         'Isotope score = %.1f' % isotope_score,
                         verticalalignment='top', horizontalalignment='right',
                         color=(0.0, 0.5, 0.9, 1.0), fontsize=10)

        m0_theo_base_box = patches.Rectangle((lib_mz - ms1_delta, 0), 2 * ms1_delta, deconv_lst[0],
                                             facecolor=(1.0, 0.0, 0.0, 0.6), edgecolor='none', zorder=1)
        ms_zoom_pic.add_patch(m0_theo_base_box)
        m0_theo_box = patches.Rectangle((lib_mz - ms1_delta, deconv_lst[0]), 2 * ms1_delta, ms1_pr_i - deconv_lst[0],
                                        facecolor=(0, 0.8, 1.0, 0.6), edgecolor='none', zorder=1)
        ms_zoom_pic.add_patch(m0_theo_box)

        ms_zoom_title_str = 'Theoretical isotopic distribution for %s %s' % (formula_charged, charge)
        ms_zoom_pic.set_title(ms_zoom_title_str, color='b', fontsize=8, y=0.98)

        # print(core_count, 'plot MS zoom in ', time.time() - _t_img_0)

    def plot_msms():
        # _t_img_0 = time.time()
        # print(core_count, 'start to plot FULL MS/MS ...')
        msms_pic = pic_array[0, 1]

        # plot identification table
        ident_col_labels = ('Proposed_structure', 'Score')
        _ident_table_df = pd.DataFrame(data={'Proposed_structure': obs_info_df['DISCRETE_ABBR'].values.tolist(),
                                             'Score': obs_info_df['RANK_SCORE'].values.tolist()})
        _ident_table_df.sort_values(by=['Score'], ascending=False, inplace=True)
        _ident_table_df.drop_duplicates(subset=['Proposed_structure', 'Score'], keep='first', inplace=True)
        ident_table_vals = list(map(list, _ident_table_df.values))
        ident_col_width_lst = [0.6, 0.15]
        msms_pic.tick_params(axis='both', which='major', labelsize=10)
        ident_table = msms_pic.table(cellText=ident_table_vals, colWidths=ident_col_width_lst,
                                     colLabels=ident_col_labels, loc='upper center', cellLoc='center')
        ident_table.set_fontsize(8)

        # plot MS/MS
        msms_pic.stem(ms2_df['mz'].values.tolist(), ms2_df['i'].values.tolist(),
                      'black', markerfmt=' ', basefmt='k-')  # zorder=10
        msms_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        msms_pic.set_xlabel("m/z", fontsize=10, labelpad=-1)
        msms_pic.set_ylabel("Intensity", fontsize=10)
        if min(ms2_df['mz'].values.tolist()) > 400:
            msms_pic.set_xlim([min(ms2_df['mz'].values.tolist()) - 100, ms2_pr_mz + 20])
        elif min(ms2_df['mz'].values.tolist()) - 10 > 0:
            msms_pic.set_xlim([min(ms2_df['mz'].values.tolist()) - 10, ms2_pr_mz + 20])
        else:
            msms_pic.set_xlim([min(ms2_df['mz'].values.tolist()) - 1, ms2_pr_mz + 20])
        msms_pic.set_ylim([0, _msms_max * 1.5])

        if obs_ident_df is not False:
            marker_l, stem_l, base_l = msms_pic.stem(obs_ident_df['mz'], obs_ident_df['i'], markerfmt=' ')
            plt.setp(stem_l, color=(0, 0.7, 1.0, 0.7), linewidth=3)
        else:
            pass
        if other_frag_df is not False:
            marker_l, stem_l, base_l = msms_pic.stem(other_frag_df['mz'], other_frag_df['i'], markerfmt=' ')
            plt.setp(stem_l, color=(0.8, 0.0, 0.0, 0.7), linewidth=3, alpha=0.4)
        else:
            pass
        if other_nl_df is not False:
            marker_l, stem_l, base_l = msms_pic.stem(other_nl_df['mz'], other_nl_df['i'], markerfmt=' ')
            plt.setp(stem_l, color=(0.8, 0.0, 0.0, 0.7), linewidth=3, alpha=0.4)
        else:
            pass
        if plt_obs_fa_df is not False:
            marker_l, stem_l, base_l = msms_pic.stem(plt_obs_fa_df['mz'], plt_obs_fa_df['i'], markerfmt=' ')
            plt.setp(stem_l, color=(0.8, 0.0, 0.0, 0.7), linewidth=3, alpha=0.4)
        else:
            pass
        if plt_obs_lyso_df is not False:
            marker_l, stem_l, base_l = msms_pic.stem(plt_obs_lyso_df['mz'], plt_obs_lyso_df['i'], markerfmt=' ')
            plt.setp(stem_l, color=(0.8, 0.0, 0.0, 0.7), linewidth=3, alpha=0.4)
        else:
            pass
        if target_frag_df is not False:
            marker_l, stem_l, base_l = msms_pic.stem(target_frag_df['mz'], target_frag_df['i'], markerfmt=' ')
            plt.setp(stem_l, color=(0.0, 0.5, 0.9, 0.7), linewidth=3, alpha=0.7)
        else:
            pass
        if target_nl_df is not False:
            marker_l, stem_l, base_l = msms_pic.stem(target_nl_df['mz'], target_nl_df['i'], markerfmt=' ')
            plt.setp(stem_l, color=(0.0, 0.5, 0.9, 0.7), linewidth=3, alpha=0.7)
        else:
            pass

        msms_title_str = ('MS/MS for m/z %.4f | DDA rank %d @ %.3f min' % (ms2_pr_mz, func_id, ms2_rt))
        msms_pic.set_title(msms_title_str, color='b', fontsize=8, y=0.98)

        # print(core_count, 'plot FULL MSMS in ', time.time() - _t_img_0)

    def plot_msms_low():
        # _t_img_0 = time.time()
        # print(core_count, 'start to plot MS/MS <= m/z 400 ...')
        msms_low_pic = pic_array[1, 1]
        msms_low_pic.tick_params(axis='both', which='major', labelsize=10)

        # plot fa frag identification table
        if not obs_fa_df.empty:
            _fa_col_labels = ('#', 'identification', 'm/z', 'ppm', 'i (%)')
            # Create a second table with the MG fragments for TG
            if obs_fa_df.loc[obs_fa_df['TYPE'] == 'MG']['obs_abbr'].shape[0] > 0:
                _fa_table_df2 = pd.DataFrame(data={
                    '#': obs_fa_df.loc[obs_fa_df['TYPE'] == 'MG']['obs_rank'].astype(int).values.tolist(),
                    'identification': obs_fa_df.loc[obs_fa_df['TYPE'] == 'MG']['obs_abbr'].values.tolist(),
                    'm/z': obs_fa_df.loc[obs_fa_df['TYPE'] == 'MG']['obs_mz'].values.tolist(),
                    'ppm': obs_fa_df.loc[obs_fa_df['TYPE'] == 'MG']['obs_ppm'].values.tolist(),
                    'i (%)': obs_fa_df.loc[obs_fa_df['TYPE'] == 'MG']['obs_i_r'].values.tolist()})
                _fa_table_df2 = _fa_table_df2.reindex(columns=_fa_col_labels)
                _fa_table_vals2 = list(map(list, _fa_table_df2.values))
                _fa_col_width_lst2 = [0.03, 0.2, 0.10, 0.06, 0.06]
                _fa_table2 = msms_low_pic.table(cellText=_fa_table_vals2, colWidths=_fa_col_width_lst2,
                                                colLabels=_fa_col_labels, loc='upper right', cellLoc='center')
                _fa_table2.set_fontsize(5)
                for _frag_idx in mg_idx_lst:
                    for _r in [0, 1, 2, 3, 4]:
                        _cell = _fa_table2.get_celld()[(_frag_idx + 1, _r)]
                        _cell.set_color((0, 0.7, 1.0, 0.4))
            else:
                _fa_table2 = False
            # Create the table with the Free FA fragments for TG and PL
            if obs_fa_df.loc[obs_fa_df['TYPE'] == 'FA']['obs_abbr'].values.tolist():
                _fa_table_df = pd.DataFrame(
                    data={'#': obs_fa_df.loc[obs_fa_df['TYPE'] == 'FA']['obs_rank'].astype(int).values.tolist(),
                          'identification': obs_fa_df.loc[obs_fa_df['TYPE'] == 'FA']['obs_abbr'].values.tolist(),
                          'm/z': obs_fa_df.loc[obs_fa_df['TYPE'] == 'FA']['obs_mz'].values.tolist(),
                          'ppm': obs_fa_df.loc[obs_fa_df['TYPE'] == 'FA']['obs_ppm'].values.tolist(),
                          'i (%)': obs_fa_df.loc[obs_fa_df['TYPE'] == 'FA']['obs_i_r'].values.tolist()})
                _fa_table_df = _fa_table_df.reindex(columns=_fa_col_labels)
                _fa_table_vals = list(map(list, _fa_table_df.values))
                _fa_col_width_lst = [0.03, 0.2, 0.10, 0.06, 0.06]
                _fa_table = msms_low_pic.table(cellText=_fa_table_vals, colWidths=_fa_col_width_lst,
                                               colLabels=_fa_col_labels, loc='upper left', cellLoc='center')

                _fa_table.set_fontsize(5)
                for _frag_idx in frag_idx_lst:
                    for _r in [0, 1, 2, 3, 4]:
                        _cell = _fa_table.get_celld()[(_frag_idx + 1, _r)]
                        _cell.set_color((0, 0.7, 1.0, 0.4))
            else:
                _fa_table = False

        # msms spectrum zoomed < 400 start
        if not _msms_low_df.empty:
            msms_low_pic.stem(_msms_low_df['mz'].values.tolist(),
                              _msms_low_df['i'].values.tolist(),
                              'black', markerfmt=' ')
            msms_low_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            msms_low_pic.set_xlabel("m/z", fontsize=10, labelpad=-1)
            msms_low_pic.set_ylabel("Intensity", fontsize=10)
            msms_low_pic.set_xlim([min(_msms_low_df['mz'].values.tolist()) - 1, 400])
            msms_low_pic.set_ylim([0, max(_msms_low_df['i'].values.tolist()) * 2])
        else:
            pass

        if isinstance(obs_ident_df, pd.DataFrame):
            low_obs_ident_df = obs_ident_df[obs_ident_df['mz'] <= 400]
            low_obs_ident_df.is_copy = False
            if not low_obs_ident_df.empty:
                marker_l, stem_l, base_l = msms_low_pic.stem(low_obs_ident_df['mz'],
                                                             low_obs_ident_df['i_r'], markerfmt=' ')
                plt.setp(stem_l, color=(0, 0.7, 1.0, 0.4), linewidth=3)
                for _i_idx, _ident_se in low_obs_ident_df.iterrows():
                    _ident_mz = _ident_se['mz']
                    _ident_i_r = _ident_se['i_r']
                    msms_low_pic.text(_ident_mz, _ident_i_r, _ident_se['obs_label'], txt_props,
                                      fontsize=8, color=(0, 0.6, 1.0, 1.0), rotation=60, weight='bold')
            else:
                pass
        else:
            pass

        if isinstance(plt_obs_fa_df, pd.DataFrame):
            marker_l, stem_l, base_l = msms_low_pic.stem(plt_obs_fa_df['mz'], plt_obs_fa_df['i'], markerfmt=' ')
            plt.setp(stem_l, color=(0.8, 0.0, 0.0, 0.3), linewidth=3, alpha=0.4)
            for _i_f_idx, _frag_se in plt_obs_fa_df.iterrows():
                _frag_mz = _frag_se['mz']
                _frag_i_r = _frag_se['i_r']
                msms_low_pic.text(_frag_mz, _frag_i_r, _frag_se['obs_label'], txt_props, fontsize=6,
                                  color=(0.8, 0.0, 0.0, 1), rotation=60)
        else:
            pass

        if isinstance(plt_obs_mg_df, pd.DataFrame):
            marker_l, stem_l, base_l = msms_low_pic.stem(plt_obs_mg_df['mz'], plt_obs_mg_df['i'], markerfmt=' ')
            plt.setp(stem_l, color=(0.8, 0.0, 0.0, 0.3), linewidth=3, alpha=0.4)
            for _i_f_idx, _frag_se in plt_obs_mg_df.iterrows():
                _frag_mz = _frag_se['mz']
                _frag_i_r = _frag_se['i_r']
                msms_low_pic.text(_frag_mz, _frag_i_r, _frag_se['obs_label'], txt_props, fontsize=6,
                                  color=(0.8, 0.0, 0.0, 1), rotation=60)
        else:
            pass
        # add specific ion info
        if isinstance(other_frag_df, pd.DataFrame):
            marker_l, stem_l, base_l = msms_low_pic.stem(other_frag_df['mz'], other_frag_df['i'], markerfmt=' ')
            plt.setp(stem_l, color=(0.8, 0.0, 0.0, 0.5), linewidth=3)
            for _o_f_idx, _frag_se in other_frag_df.iterrows():
                _frag_mz = _frag_se['mz']
                _frag_i_r = _frag_se['i_r']
                _frag_class = _frag_se['LABEL']
                msms_low_pic.text(_frag_mz, _frag_i_r, _frag_class, fontsize=7, color=(0.8, 0.0, 0.0, 1))
        else:
            pass

        if isinstance(target_frag_df, pd.DataFrame):
            marker_l, stem_l, base_l = msms_low_pic.stem(target_frag_df['mz'], target_frag_df['i'], markerfmt=' ')
            plt.setp(stem_l, color=(0.0, 0.5, 0.9, 0.6), linewidth=3)
            for _t_f_idx, _frag_se in target_frag_df.iterrows():
                _frag_mz = _frag_se['mz']
                _frag_i_r = _frag_se['i_r']
                _frag_class = _frag_se['LABEL']
                msms_low_pic.text(_frag_mz, _frag_i_r, _frag_class, fontsize=8, color=(0.0, 0.5, 0.9, 1.0),
                                  weight='bold')
        else:
            pass

        # msms_low_pic.set_ylim([0, _msms_max * 1.5])
        msms_low_str = 'MS/MS zoomed below m/z 400'
        msms_low_pic.set_title(msms_low_str, color='b', fontsize=8, y=0.98)

        # print(core_count, 'plot MSMS <= 400 in ', time.time() - _t_img_0)

    def plot_msms_high():
        # _t_img_0 = time.time()
        # print(core_count, 'start to plot MS/MS > m/z 400 ...')
        msms_high_pic = pic_array[2, 1]
        msms_high_pic.tick_params(axis='both', which='major', labelsize=10)

        # plot fa nl identification table
        # For TG in [M+Na]+ there will be 2 different plots in this section
        if not obs_lyso_df.empty:
            _lyso_col_labels = ('#', 'identification', 'm/z', 'ppm', 'i (%)')
            if obs_lyso_df.loc[obs_lyso_df['TYPE'] == 'NL_Na']['obs_abbr'].shape[0] > 0:
                _lyso_table_df2 = pd.DataFrame(
                    data={'#': obs_lyso_df.loc[obs_lyso_df['TYPE'] == 'NL_Na']['obs_rank'].astype(int).values.tolist(),
                          'identification': obs_lyso_df.loc[obs_lyso_df['TYPE'] == 'NL_Na']['obs_abbr'].values.tolist(),
                          'm/z': obs_lyso_df.loc[obs_lyso_df['TYPE'] == 'NL_Na']['obs_mz'].values.tolist(),
                          'ppm': obs_lyso_df.loc[obs_lyso_df['TYPE'] == 'NL_Na']['obs_ppm'].values.tolist(),
                          'i (%)': obs_lyso_df.loc[obs_lyso_df['TYPE'] == 'NL_Na']['obs_i_r'].values.tolist()})
                _lyso_table_df2 = _lyso_table_df2.reindex(columns=_lyso_col_labels)
                _lyso_table_vals2 = list(map(list, _lyso_table_df2.values))
                _lyso_col_width_lst2 = [0.03, 0.2, 0.10, 0.06, 0.06]
                _lyso_table2 = msms_high_pic.table(cellText=_lyso_table_vals2, colWidths=_lyso_col_width_lst2,
                                                   colLabels=_lyso_col_labels, loc='upper right', cellLoc='center',
                                                   bbox=[0.5, 1 - (0.067 * (len(_lyso_table_vals2) + 1)), 0.45,
                                                         0.067 * (len(_lyso_table_vals2) + 1)])

                # set back ground of the table
                _lyso_table2.set_fontsize(5)
                for _nl_idx in dg_Na_idx_lst:
                    for _r in [0, 1, 2, 3, 4]:
                        if _lyso_table2 is False:
                            pass
                        else:
                            _cell = _lyso_table2.get_celld()[(_nl_idx + 1, _r)]
                            _cell.set_color((0, 0.7, 1.0, 0.4))
            else:
                _lyso_table2 = False
            if obs_lyso_df.loc[obs_lyso_df['TYPE'] == 'NL']['obs_abbr'].shape[0] > 0:
                _lyso_table_df = pd.DataFrame(
                    data={'#': obs_lyso_df.loc[obs_lyso_df['TYPE'] == 'NL']['obs_rank'].astype(int).values.tolist(),
                          'identification': obs_lyso_df.loc[obs_lyso_df['TYPE'] == 'NL']['obs_abbr'].values.tolist(),
                          'm/z': obs_lyso_df.loc[obs_lyso_df['TYPE'] == 'NL']['obs_mz'].values.tolist(),
                          'ppm': obs_lyso_df.loc[obs_lyso_df['TYPE'] == 'NL']['obs_ppm'].values.tolist(),
                          'i (%)': obs_lyso_df.loc[obs_lyso_df['TYPE'] == 'NL']['obs_i_r'].values.tolist()})
                _lyso_table_df = _lyso_table_df.reindex(columns=_lyso_col_labels)
                _lyso_table_vals = list(map(list, _lyso_table_df.values))
                _lyso_col_width_lst = [0.03, 0.2, 0.10, 0.06, 0.06]
                if charge in ['[M+H]+', '[M+NH4]+'] and abbr[0:2] == 'TG':
                    _lyso_table = msms_high_pic.table(cellText=_lyso_table_vals, colWidths=_lyso_col_width_lst,
                                                      colLabels=_lyso_col_labels, loc='upper center', cellLoc='center',
                                                      bbox=[0.4, 1 - (0.067 * (len(_lyso_table_vals) + 1)), 0.5, 0.067*(len(_lyso_table_vals)+1)])
                elif charge in ['[M+Na]+'] and abbr[0:2] == 'TG':
                    _lyso_table = msms_high_pic.table(cellText=_lyso_table_vals, colWidths=_lyso_col_width_lst,
                                                      colLabels=_lyso_col_labels, loc='upper left', cellLoc='center',
                                                      bbox=[0, 1 - (0.067 * (len(_lyso_table_vals) + 1)), 0.45,
                                                                    0.067 * (len(_lyso_table_vals) + 1)])
                else:
                    _lyso_table = msms_high_pic.table(cellText=_lyso_table_vals, colWidths=_lyso_col_width_lst,
                                                      colLabels=_lyso_col_labels, loc='upper center', cellLoc='center')

                # set back ground of the table
                _lyso_table.set_fontsize(5)
                for _nl_idx in nl_idx_lst:
                    for _r in [0, 1, 2, 3, 4]:
                        if _lyso_table is False:
                            pass
                        else:
                            _cell = _lyso_table.get_celld()[(_nl_idx + 1, _r)]
                            _cell.set_color((0, 0.7, 1.0, 0.4))
            else:
                _lyso_table = False

        # msms spectrum zoomed > 400 start
        if not _msms_high_df.empty:
            msms_high_max = _msms_high_df['i'].max()
            msms_high_pic.stem(_msms_high_df['mz'].values.tolist(), _msms_high_df['i'].values.tolist(),
                               'black', markerfmt=' ')
            msms_high_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            msms_high_pic.set_xlabel("m/z", fontsize=10, labelpad=-1)
            msms_high_pic.set_ylabel("Intensity", fontsize=10)
            msms_high_pic.set_xlim([400, ms2_pr_mz + 20])
            msms_high_pic.set_ylim([0, msms_high_max * 1.3])
        else:
            msms_high_max = ms2_df['i'].max()
            msms_high_pic.set_ylim([0, msms_high_max * 1.25])

        if isinstance(obs_ident_df, pd.DataFrame):
            high_obs_ident_df = obs_ident_df[obs_ident_df['mz'] > 400]
            high_obs_ident_df.is_copy = False
            high_obs_ident_df['i_r'] = high_obs_ident_df['i'] * 1.025
            if not high_obs_ident_df.empty:
                marker_l, stem_l, base_l = msms_high_pic.stem(high_obs_ident_df['mz'],
                                                              high_obs_ident_df['i_r'], markerfmt=' ')
                plt.setp(stem_l, color=(0, 0.7, 1.0, 0.4), linewidth=3)
                for _i_idx, _ident_se in high_obs_ident_df.iterrows():
                    _ident_mz = _ident_se['mz']
                    _ident_i_r = _ident_se['i_r']
                    msms_high_pic.text(_ident_mz, _ident_i_r, _ident_se['obs_label'], txt_props,
                                       fontsize=8, color=(0, 0.6, 1.0, 1.0), rotation=60, weight='bold')
            else:
                pass
        else:
            pass

        if isinstance(plt_obs_lyso_df, pd.DataFrame):
            marker_l, stem_l, base_l = msms_high_pic.stem(plt_obs_lyso_df['mz'], plt_obs_lyso_df['i'], markerfmt=' ')
            plt.setp(stem_l, color=(0.8, 0.0, 0.0, 0.3), linewidth=3, alpha=0.4)
            for _i_nl_idx, _nl_se in plt_obs_lyso_df.iterrows():
                _nl_mz = _nl_se['mz']
                _nl_i_r = _nl_se['i_r']
                msms_high_pic.text(_nl_mz, _nl_i_r, _nl_se['obs_label'], txt_props, fontsize=6,
                                   color=(0.8, 0.0, 0.0, 1), rotation=60)
        else:
            pass

        if isinstance(plt_obs_dg_na_df, pd.DataFrame):
            marker_l, stem_l, base_l = msms_high_pic.stem(plt_obs_dg_na_df['mz'], plt_obs_dg_na_df['i'], markerfmt=' ')
            plt.setp(stem_l, color=(0.8, 0.0, 0.0, 0.3), linewidth=3, alpha=0.4)
            for _i_nl_idx, _nl_se in plt_obs_dg_na_df.iterrows():
                _nl_mz = _nl_se['mz']
                _nl_i_r = _nl_se['i_r']
                msms_high_pic.text(_nl_mz, _nl_i_r, _nl_se['obs_label'], txt_props, fontsize=6,
                                   color=(0.8, 0.0, 0.0, 1), rotation=60)
        else:
            pass

        if isinstance(other_nl_df, pd.DataFrame):
            marker_l, stem_l, base_l = msms_high_pic.stem(other_nl_df['mz'], other_nl_df['i'], markerfmt=' ')
            plt.setp(stem_l, color=(0.8, 0.0, 0.0, 0.5), linewidth=3)
            for _o_nl_idx, _nl_se in other_nl_df.iterrows():
                _nl_mz = _nl_se['mz']
                _nl_i = _nl_se['i']
                _nl_class = _nl_se['LABEL']
                _nl_i_r = _nl_i * 1.2
                msms_high_pic.text(_nl_mz, _nl_i_r, _nl_class, fontsize=7, color=(0.8, 0.0, 0.0, 1))
        else:
            pass

        if isinstance(target_nl_df, pd.DataFrame):
            marker_l, stem_l, base_l = msms_high_pic.stem(target_nl_df['mz'], target_nl_df['i_r'], markerfmt=' ')
            plt.setp(stem_l, color=(0.0, 0.5, 0.9, 0.6), linewidth=3)
            for _t_nl_idx, _nl_se in target_nl_df.iterrows():
                _nl_mz = _nl_se['mz']
                _nl_i_r = _nl_se['i_r']
                _nl_class = _nl_se['LABEL']
                # _nl_i_r = _nl_i * 1.2
                msms_high_pic.text(_nl_mz, _nl_i_r, _nl_class, fontsize=8, color=(0.0, 0.5, 0.9, 1.0), weight='bold')
        else:
            pass

        msms_high_str = 'MS/MS zoomed above m/z 400'
        msms_high_pic.set_title(msms_high_str, color='b', fontsize=8, y=0.98)
        # print(core_count, 'plot MSMS > 400 ', time.time() - _t_img_0)

    try:
        tasks = [plot_msms(), plot_msms_low(), plot_msms_high(), plot_xic(), plot_ms(), plot_ms_zoom()]
        with ThreadPoolExecutor(max_workers=2) as executor:
            for _task in tasks:
                executor.submit(_task)

        plt.savefig(save_img_as, type=img_type, dpi=dpi)
        print(core_count, '=====> Image saved as: %s' % save_img_as)
        plt.close()
    except Exception as e:
        print(e)
        print('Use single thread and try again ...')
        plot_msms()
        plot_msms_low()
        plot_msms_high()
        plot_xic()
        plot_ms()
        plot_ms_zoom()
        plt.savefig(save_img_as, type=img_type, dpi=dpi)
        print(core_count, '=====> Image saved as: %s' % save_img_as)
        plt.close()


def gen_plot(param_dct_lst, core_count, img_type='png', dpi=300, vendor='waters', ms1_precision=50e-6):
    core_count = 'Core #{core}'.format(core=core_count)

    if isinstance(param_dct_lst, list):
        img_counter = 1
        tot_img_count = len(param_dct_lst)
        for param_dct in param_dct_lst:
            abbr = param_dct['abbr']
            mz_se = param_dct['mz_se']
            xic_dct = param_dct['xic_dct']
            ident_info_dct = param_dct['ident_info_dct']
            spec_info_dct = param_dct['spec_info_dct']
            isotope_score_info_dct = param_dct['isotope_score_info_dct']
            specific_dct = param_dct['specific_dct']
            formula_charged = param_dct['formula_charged']
            charge = param_dct['charge']
            save_img_as = param_dct['save_img_as']
            print('%s ==> image: %i / %i' % (core_count, img_counter, tot_img_count))
            try:
                plot_spectra(abbr, mz_se, xic_dct, ident_info_dct, spec_info_dct, isotope_score_info_dct, specific_dct,
                             formula_charged, charge, core_count, save_img_as=save_img_as, img_type=img_type,
                             dpi=dpi, vendor=vendor, ms1_precision=ms1_precision)
            except Exception as e:
                print(e)

            img_counter += 1

    else:
        if isinstance(param_dct_lst, dict):
            param_dct = param_dct_lst
            abbr = param_dct['abbr']
            mz_se = param_dct['mz_se']
            xic_dct = param_dct['xic_dct']
            ident_info_dct = param_dct['ident_info_dct']
            spec_info_dct = param_dct['spec_info_dct']
            isotope_score_info_dct = param_dct['isotope_score_info_dct']
            specific_dct = param_dct['specific_dct']
            formula_charged = param_dct['formula_charged']
            charge = param_dct['charge']
            save_img_as = param_dct['save_img_as']
            try:
                plot_spectra(abbr, mz_se, xic_dct, ident_info_dct, spec_info_dct, isotope_score_info_dct, specific_dct,
                             formula_charged, charge, core_count, save_img_as=save_img_as, img_type=img_type,
                             dpi=dpi, vendor=vendor, ms1_precision=ms1_precision)
            except Exception as e:
                print(e)
