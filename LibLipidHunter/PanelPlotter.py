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

    print(core_count, '[STATUS] >>> Start to plot %s -> MS2 PR m/z %.4f @ MS1 best PR m/z %.4f with lib m/z %.4f'
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
    if not ms_zoom_df.empty:
        if 'i' in ms_zoom_df.columns.tolist():
            ms_zoom_bp_i = max(ms_zoom_df['i'].values.tolist())
        else:
            if not ms1_df.empty:
                if 'i' in ms1_df.columns.tolist():
                    ms_zoom_bp_i = ms1_df['i'].max()
                else:
                    ms_zoom_bp_i = 0
            else:
                ms_zoom_bp_i = 0
    else:
        if not ms1_df.empty:
            if 'i' in ms1_df.columns.tolist():
                ms_zoom_bp_i = ms1_df['i'].max()
            else:
                ms_zoom_bp_i = 0
        else:
            ms_zoom_bp_i = 0

    xic_df = xic_dct[ms1_xic_mz]

    xic_rt_lst = xic_df['rt'].values.tolist()
    xic_i_lst = xic_df['i'].values.tolist()

    # if ms_zoom_bp_i > 0 and len(xic_rt_lst) > 0 and len(xic_i_lst) > 0:

    # cut lower peaks to accelerate plotting time
    try:
        m1_dct = isotope_checker_dct[1]
        m1_theo_mz = m1_dct['theo_mz']
        m1_theo_i = m1_dct['theo_i']
        m1_obs_mz = m1_dct['obs_mz']
        m1_obs_i = m1_dct['obs_i']
    except KeyError:
        m1_theo_mz = ms1_obs
        m1_theo_i = ms1_df['i'].max()
        m1_obs_mz = ms1_obs
        m1_obs_i = ms1_df['i'].max()
        ms1_pr_mz = ms1_obs
        ms_zoom_bp_i = ms1_df['i'].max()

    if ms1_df['i'].max() >= 10000 and ms1_df.shape[0] >= 500:
        ms1_min = ms1_df['i'].min()
        ms1_max = ms1_df['i'].max()
        ms1_top1000_i = sorted(ms1_df['i'].values.tolist(), reverse=True)[499]
        ms1_plot_th = min(m1_obs_i, 3 * ms1_min, ms1_max * 0.01, 1000, ms1_top1000_i)
        ms1_plot_th = max(ms1_plot_th, ms1_top1000_i)
        # print(core_count, m1_obs_i, 3 * ms1_min, ms1_max * 0.01, 1000, ms1_top1000_i)
        ms1_df = ms1_df.query('i >= %f' % ms1_plot_th)
        print(core_count, '[INFO] --> Plot full MS1 with abs intensity filter > %f' % ms1_plot_th)
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
            print(core_count, '[INFO] --> Plot full MS/MS with abs intensity filter > %f' % ms2_plot_th)

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
    # obs_ident_df['i_r'] = obs_ident_df['i'] * 1.025
    if not plt_obs_fa_df.empty:
        pass
        # plt_obs_fa_df['i_r'] = obs_fa_df['i'] * 1.025
    else:
        plt_obs_fa_df = False
    if not plt_obs_mg_df.empty:
        pass
        # plt_obs_mg_df['i_r'] = obs_mg_df['i'] * 1.025
    else:
        plt_obs_mg_df = False
    if not plt_obs_lyso_df.empty:
        pass
        # plt_obs_lyso_df['i_r'] = obs_lyso_df['i'] * 1.075
    else:
        plt_obs_lyso_df = False
    if not plt_obs_dg_na_df.empty:
        pass
        # plt_obs_dg_na_df['i_r'] = obs_dg_na_df['i'] * 1.075
    else:
        plt_obs_dg_na_df = False

    if 'OTHER_FRAG' in list(specific_dct.keys()):
        other_frag_df = specific_dct['OTHER_FRAG']
        # other_frag_df['i_r'] = other_frag_df['i'] * 1.4
    else:
        other_frag_df = False
    if 'OTHER_NL' in list(specific_dct.keys()):
        other_nl_df = specific_dct['OTHER_NL']
        # other_nl_df['i_r'] = other_nl_df['i'] * 1.2
    else:
        other_nl_df = False
    if 'TARGET_FRAG' in list(specific_dct.keys()):
        target_frag_df = specific_dct['TARGET_FRAG']
        # target_frag_df['i_r'] = target_frag_df['i'] * 1.4
    else:
        target_frag_df = False
    if 'TARGET_NL' in list(specific_dct.keys()):
        target_nl_df = specific_dct['TARGET_NL']
        # target_nl_df['i_r'] = target_nl_df['i'] * 1.2
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
        plt.setp(_stem_lines, color=(0.0, 0.4, 1.0, 0.3), linewidth=3)
        _marker_line, _stem_lines, _base_line = xic_pic.stem([ms2_rt], [max(xic_i_lst)], '--', markerfmt=' ')
        plt.setp(_stem_lines, color=(0.0, 0.65, 1.0, 0.95), linewidth=2)
        xic_pic.text(ms1_rt - xic_rt_label_shift, max(xic_i_lst) * 0.98, 'MS',
                     fontsize=8, color=(0.0, 0.4, 1.0, 1.0), weight='bold')
        xic_pic.text(ms2_rt, max(xic_i_lst) * 0.98, 'MS/MS',
                     fontsize=8, color=(0.0, 0.65, 1.0, 1.0), weight='bold')
        xic_pic.set_xlabel("Scan time (min)", fontsize=7, labelpad=-1)
        xic_pic.set_ylabel("Intensity", fontsize=7)
        xic_pic.set_xlim([xic_rt_min, xic_rt_max])
        xic_pic.set_ylim([0, max(xic_i_lst) * 1.1])
        xic_title_str = 'XIC of m/z %.4f @ %s m/z %.4f ppm=%.2f' % (ms1_pr_mz, abbr, lib_mz, ms1_pr_ppm)
        xic_pic.set_title(xic_title_str, color=(0.0, 0.4, 1.0, 1.0), fontsize=8, y=0.98)
        # print(core_count, 'plot XIC in ', time.time() - _t_img_0)

    def plot_ms():
        # _t_img_0 = time.time()
        # print(core_count, 'start to plot MS ...')
        ms_pic = pic_array[1, 0]
        ms_pic.tick_params(axis='both', which='major', labelsize=10)

        if ms1_df.shape[0] > 700:
            ms_pic.plot(ms1_df['mz'].values.tolist(), ms1_df['i'].values.tolist(), 'grey', lw=0.6)
        else:
            marker_l, stem_l, base_l = ms_pic.stem(ms1_df['mz'].values.tolist(), ms1_df['i'].values.tolist(),
                                                   markerfmt=' ')
            plt.setp(stem_l, color='grey', lw=0.6)
            plt.setp(base_l, visible=False)
        _marker_l, _stem_l, _base_l = ms_pic.stem([ms1_pr_mz], dash_i, markerfmt=' ')
        plt.setp(_stem_l, color=(0.0, 0.4, 1.0, 0.2), linewidth=4)
        plt.setp(_base_l, visible=False)
        _marker_l, _stem_l, _base_l = ms_pic.stem([ms1_pr_mz], [ms1_pr_i], markerfmt='D')
        # plt.setp(_stem_lines, color=(0.4, 1.0, 0.8, 1.0))
        plt.setp(_marker_l, markerfacecolor=(0.3, 0.9, 1.0, 0.8), markersize=5, markeredgewidth=0)
        plt.setp(_stem_l, visible=False)
        plt.setp(_base_l, visible=False)
        ms_pic.text(ms1_pr_mz, ms1_pr_i, '%.4f' % ms1_pr_mz,
                    fontsize=7, color=(0.0, 0.4, 1.0, 1.0))
        ms_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        ms_pic.set_xlabel("m/z", fontsize=7, labelpad=-1)
        ms_pic.set_ylabel("Intensity", fontsize=7)
        ms_pic.set_ylim([0, max(ms1_df['i'].values.tolist()) * 1.3])

        # add annotation
        _ms_pkl_top_df = ms1_df.sort_values(by='i', ascending=False).head(10)
        _ms_pkl_top_peak_list = list(zip(_ms_pkl_top_df['mz'].values.tolist(), _ms_pkl_top_df['i'].values.tolist()))
        for _ms_pkl_top_peak in _ms_pkl_top_peak_list:
            if _ms_pkl_top_peak[0] < ms1_pr_mz - 5 or _ms_pkl_top_peak[0] > ms1_pr_mz + 5:  # avoid overlay with pr_mz
                _ms_pkl_top_peak_str = '%.4f' % _ms_pkl_top_peak[0]
                _ms_pkl_top_peak_y = _ms_pkl_top_peak[1]
                ms_pic.text(_ms_pkl_top_peak[0], _ms_pkl_top_peak_y, _ms_pkl_top_peak_str, fontsize=6)

        ms_title_str = 'MS @ %.3f min ' % ms1_rt
        ms_pic.set_title(ms_title_str, color=(0.0, 0.4, 1.0, 1.0), fontsize=8, y=0.98)
        # print(core_count, 'plot MS in ', time.time() - _t_img_0)

    def plot_ms_zoom():
        # _t_img_0 = time.time()
        # print(core_count, 'start to plot MS zoom ...')
        ms_zoom_pic = pic_array[2, 0]
        ms_zoom_pic.tick_params(axis='both', which='major', labelsize=10)

        ms_zoom_pic.set_xlim([ms1_pr_mz - 1.5, ms1_pr_mz + 3.55])
        ms_zoom_pic.set_ylim([0, ms_zoom_bp_i * 1.45])
        ms_zoom_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        ms_zoom_pic.ticklabel_format(axis='x', useOffset=False)
        ms_zoom_pic.set_xlabel('m/z', fontsize=7, labelpad=-1)
        ms_zoom_pic.set_ylabel('Intensity', fontsize=7)

        plt_ms_zoom_df = ms_zoom_df.sort_values(by='mz', ascending='True')
        theo_i_bar_color = (0, 0.4, 1.0, 0.5)
        theo_i2_bar_color = (0, 0.4, 1.0, 0.3)
        if vendor == 'waters':
            ms_zoom_pic.plot(plt_ms_zoom_df['mz'].values.tolist(), plt_ms_zoom_df['i'].values.tolist(),
                             'grey', lw=1, zorder=1)
        elif vendor == 'thermo':
            marker_l, stem_l, base_l = ms_zoom_pic.stem(plt_ms_zoom_df['mz'].values.tolist(),
                                                        plt_ms_zoom_df['i'].values.tolist(), markerfmt=' ')
            plt.setp(stem_l, color='grey', lw=0.75, alpha=0.6)
            plt.setp(stem_l, zorder=1)
            plt.setp(base_l, visible=False)
            theo_i_bar_color = (0.2, 0.8, 1.0, 0.8)
            theo_i2_bar_color = (0.2, 0.8, 1.0, 0.5)
        else:
            ms_zoom_pic.plot(plt_ms_zoom_df['mz'].values.tolist(), plt_ms_zoom_df['i'].values.tolist(),
                             'grey', lw=1, zorder=1)

        # isotope region | if any peak in M-1.0034

        if deconv_lst[0] > 0:
            pseudo_m0_theo_df = plt_ms_zoom_df.query('%.4f <= mz <= %4f'
                                                     % ((lib_mz - 1.0034 - ms1_delta - ms1_delta),
                                                        (lib_mz - 1.0034 - ms1_delta + ms1_delta)))
            pseudo_m0_i = pseudo_m0_theo_df['i'].max()
            pseudo_m0_theo_i_box = patches.Rectangle((lib_mz - 1.0034 - ms1_delta, 0),
                                                     2 * ms1_delta, pseudo_m0_i,
                                                     facecolor=(1.0, 0.0, 0.0, 0.5), edgecolor='none', zorder=9)
            ms_zoom_pic.add_patch(pseudo_m0_theo_i_box)
            if pseudo_m0_i > 0.05 * ms1_pr_i:
                ms_zoom_pic.text(lib_mz - 1.0034 - ms1_delta, pseudo_m0_i,
                                 '%.4f' % (lib_mz - 1.0034 - ms1_delta),
                                 color='red', fontsize=7, alpha=0.8)

        # if ms1_precision < 19e-6:
        #     m_pre_theo_box = patches.Rectangle((lib_mz - 1.0034 - ms1_delta, 0), 2 * ms1_delta, ms1_pr_i,
        #                                        facecolor=(1.0, 0.0, 0.0, 0.5), edgecolor='none', zorder=1)
        # else:
        #     m_pre_theo_box = patches.Rectangle((lib_mz - 1.0034 - ms1_delta, 0), 2 * ms1_delta, ms1_pr_i,
        #                                        facecolor=(1.0, 0.0, 0.0, 0.5), edgecolor=(1.0, 0.0, 0.0, 0.5),
        #                                        zorder=1, fill=False, linewidth=0.6, linestyle='dotted')
        # ms_zoom_pic.add_patch(m_pre_theo_box)

        ms_zoom_offset_i = ms_zoom_bp_i * 0.15
        m0_theo_base_box = patches.Rectangle((lib_mz - ms1_delta, 0), 2 * ms1_delta, deconv_lst[0],
                                             facecolor=(1.0, 0.0, 0.0, 0.5), edgecolor='none', zorder=1)
        ms_zoom_pic.add_patch(m0_theo_base_box)
        m0_theo_box = patches.Rectangle((lib_mz - ms1_delta, deconv_lst[0]), 2 * ms1_delta, ms1_pr_i - deconv_lst[0],
                                        facecolor=theo_i_bar_color, edgecolor='none', zorder=10)
        ms_zoom_pic.add_patch(m0_theo_box)

        _marker_l, _stem_l, _base_l = ms_zoom_pic.stem([ms1_pr_mz], [ms1_pr_i], markerfmt='D')  # zorder=20
        plt.setp(_stem_l, visible=False)
        plt.setp(_base_l, visible=False)
        plt.setp(_marker_l, markerfacecolor=(0.2, 0.8, 1.0, 0.8), markeredgecolor='none', markeredgewidth=0,
                 markersize=6, lw=0, zorder=21)

        ms_zoom_pic.text(ms1_pr_mz + 0.06, ms1_pr_i, '%.4f' % ms1_pr_mz,
                         color=(0.0, 0.4, 1.0, 1.0), fontsize=7)

        ms_zoom_pic.text(lib_mz - 0.06, ms1_pr_i + ms_zoom_offset_i, '[M]', color=(0.0, 0.4, 1.0, 1.0), fontsize=7)
        ms_zoom_pic.text(lib_mz - 0.78, ms1_pr_i, 'Calc: %.4f' % lib_mz, color=(0.2, 0.8, 1.0, 1.0), fontsize=7)

        # isotope region | highlight the 1st isotope
        # theo range box
        m1_theo_base_box = patches.Rectangle((m1_theo_mz - ms1_delta, 0), 2 * ms1_delta, deconv_lst[1],
                                             facecolor=(1.0, 0.0, 0.0, 0.5), edgecolor='none', zorder=11)
        ms_zoom_pic.add_patch(m1_theo_base_box)
        m1_theo_box = patches.Rectangle((m1_theo_mz - ms1_delta, deconv_lst[1]), 2 * ms1_delta,
                                        m1_theo_i - deconv_lst[1],
                                        facecolor=theo_i_bar_color, edgecolor='none', zorder=12)
        ms_zoom_pic.add_patch(m1_theo_box)

        _marker_l, _stem_l, _base_l = ms_zoom_pic.stem([m1_theo_mz], [m1_theo_i], '--', markerfmt='o')
        plt.setp(_marker_l, markerfacecolor=(0.2, 0.8, 1.0, 0.8), markersize=6, markeredgewidth=0, zorder=22)
        plt.setp(_stem_l, visible=False)
        plt.setp(_base_l, visible=False)
        ms_zoom_pic.text(m1_theo_mz - 0.15, m1_theo_i + ms_zoom_offset_i, '[M+1]',
                         color=(0.0, 0.4, 1.0, 1.0), fontsize=7)
        ms_zoom_pic.text(m1_theo_mz - 0.78, m1_theo_i, 'Calc: %.4f' % m1_theo_mz,
                         color=(0.2, 0.8, 1.0, 1.0), fontsize=7)
        ms_zoom_pic.text(m1_obs_mz + 0.04, m1_obs_i, '%.4f' % m1_obs_mz, color=(0.0, 0.4, 1.0, 1.0), fontsize=7)

        opt_box_lst = []

        # isotope region | highlight the 2nd isotope
        if 2 in list(isotope_checker_dct.keys()):
            m2_dct = isotope_checker_dct[2]
            m2_theo_mz = m2_dct['theo_mz']
            m2_theo_i = m2_dct['theo_i']

            if 'obs_mz' in list(m2_dct.keys()) and 'obs_i' in list(m2_dct.keys()):
                m2_obs_mz = m2_dct['obs_mz']
                m2_obs_i = m2_dct['obs_i']
                ms_zoom_pic.text(m2_obs_mz + 0.04, m2_obs_i, '%.4f' % m2_obs_mz,
                                 color=(0.0, 0.4, 1.0, 1.0), fontsize=7)

                m2_theo_box = patches.Rectangle((m2_theo_mz - ms1_delta, 0), 2 * ms1_delta, m2_theo_i,
                                                facecolor=theo_i_bar_color, edgecolor='none', zorder=13)
                ms_zoom_pic.add_patch(m2_theo_box)
                opt_box_lst.append(ms_zoom_pic)
                _marker_l, _stem_l, _base_l = ms_zoom_pic.stem([m2_theo_mz], [m2_theo_i], '--', markerfmt='o')
                plt.setp(_marker_l, markerfacecolor=(0.2, 0.8, 1.0, 0.8), markersize=6, markeredgewidth=0, zorder=23)
                plt.setp(_stem_l, visible=False)
                plt.setp(_base_l, visible=False)
                ms_zoom_pic.text(m2_theo_mz - 0.15, m2_theo_i + ms_zoom_offset_i, '[M+2]',
                                 color=(0.0, 0.4, 1.0, 1.0), fontsize=7)

            else:
                if len(list(m2_checker_dct.keys())) > 0:
                    pass
                else:
                    ms_zoom_pic.text(m2_theo_mz + 0.025, m2_theo_i, '[M+2] predicted',
                                     color=(0.0, 0.4, 1.0, 1.0), fontsize=5)
                _marker_l, _stem_l, _base_l = ms_zoom_pic.stem([m2_theo_mz], [m2_theo_i], '--', markerfmt='o')
                plt.setp(_marker_l, markerfacecolor=(0.2, 0.8, 1.0, 0.8), markersize=6, markeredgewidth=0, zorder=23)
                plt.setp(_stem_l, color=(0.0, 0.65, 1.0, 0.8), linewidth=1.75)
                plt.setp(_base_l, visible=False)

            ms_zoom_pic.text(m2_theo_mz - 0.78, m2_theo_i, 'Calc: %.4f' % m2_theo_mz,
                             color=(0.2, 0.8, 1.0, 1.0), fontsize=7)

        if len(list(m2_checker_dct.keys())) > 0:
            for _mh2 in list(m2_checker_dct.keys()):
                mh2_dct = m2_checker_dct[_mh2]
                mh2_theo_mz = mh2_dct['theo_mz']
                mh2_theo_i = mh2_dct['theo_i']
                mh2_obs_mz = mh2_dct['obs_mz']
                mh2_obs_i = mh2_dct['obs_i']
                decon_idx = _mh2 + 2

                mh2_theo_base_box = patches.Rectangle((mh2_theo_mz - ms1_delta, 0), 2 * ms1_delta,
                                                      deconv_lst[decon_idx],
                                                      facecolor=theo_i2_bar_color, edgecolor='none', zorder=14)
                ms_zoom_pic.add_patch(mh2_theo_base_box)
                ms_zoom_pic.text(mh2_theo_mz + 0.025, deconv_lst[decon_idx], '[M+%i] predicted' % decon_idx,
                                 color=(0.0, 0.4, 1.0, 1.0), fontsize=5)

                mh2_theo_box = patches.Rectangle((mh2_theo_mz - ms1_delta, deconv_lst[decon_idx]),
                                                 2 * ms1_delta, mh2_theo_i - deconv_lst[decon_idx],
                                                 facecolor=(1.0, 0.0, 0.0, 0.5), edgecolor='none', zorder=15)
                ms_zoom_pic.add_patch(mh2_theo_box)
                _marker_l, _stem_l, _base_l = ms_zoom_pic.stem([mh2_theo_mz], [mh2_theo_i], '--', markerfmt='o')
                plt.setp(_stem_l, color='red', alpha=0.6)
                plt.setp(_marker_l, markerfacecolor='red', markersize=5, markeredgewidth=0, alpha=0.6, zorder=24)
                plt.setp(_stem_l, visible=False)
                if _mh2 == 0:
                    _mh2_name = ''
                else:
                    _mh2_name = '+%i' % _mh2
                ms_zoom_pic.text(mh2_theo_mz - 0.2 - 0.05 * _mh2, mh2_theo_i + ms_zoom_offset_i,
                                 '[M+2H%s]' % _mh2_name, color='red', fontsize=6)
                ms_zoom_pic.text(mh2_theo_mz - 0.78, mh2_theo_i, 'Calc: %.4f' % mh2_theo_mz,
                                 color='red', fontsize=7)
                ms_zoom_pic.text(mh2_obs_mz + 0.04, mh2_obs_i, '%.4f' % mh2_obs_mz, color='red', fontsize=7)

            # plot the M+2H isotope score

            ms_zoom_pic.text(m1_theo_mz + 2.5, max(ms_zoom_bp_i - ms_zoom_offset_i, ms1_pr_i * 0.8),
                             '[M+2H] Isotope score = %.1f' % m2_score,
                             verticalalignment='top', horizontalalignment='right',
                             color='red', fontsize=7)

        # plot the isotope score
        ms_zoom_pic.text(m1_theo_mz + 0.2, ms_zoom_bp_i + 2.5 * ms_zoom_offset_i,
                         'Isotope score = %.1f' % isotope_score,
                         verticalalignment='top', horizontalalignment='right',
                         color=(0.0, 0.4, 1.0, 1.0), fontsize=10)

        ms_zoom_title_str = 'Isotopic distribution: %s  Charge: %s  Formula: %s' % (abbr, charge, formula_charged)
        ms_zoom_pic.set_title(ms_zoom_title_str, color=(0.0, 0.4, 1.0, 1.0), fontsize=8, y=0.98)

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
        marker_l, stem_l, base_l = msms_pic.stem(ms2_df['mz'].values.tolist(), ms2_df['i'].values.tolist(),
                                                 markerfmt=' ', basefmt='k-')  # zorder=10
        plt.setp(stem_l, color='grey', linewidth=0.6)
        plt.setp(base_l, visible=False)
        msms_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        msms_pic.set_xlabel("m/z", fontsize=7, labelpad=-1)
        msms_pic.set_ylabel("Intensity", fontsize=7)
        if min(ms2_df['mz'].values.tolist()) > 400:
            msms_pic.set_xlim([min(ms2_df['mz'].values.tolist()) - 100, ms2_pr_mz + 20])
        elif min(ms2_df['mz'].values.tolist()) - 10 > 0:
            msms_pic.set_xlim([min(ms2_df['mz'].values.tolist()) - 10, ms2_pr_mz + 20])
        else:
            msms_pic.set_xlim([min(ms2_df['mz'].values.tolist()) - 1, ms2_pr_mz + 20])
        msms_pic.set_ylim([0, _msms_max * 1.5])

        if obs_ident_df is not False:
            marker_l, stem_l, base_l = msms_pic.stem(obs_ident_df['mz'], obs_ident_df['i'], markerfmt=' ')
            plt.setp(stem_l, color=(0, 0.7, 1.0, 0.6), linewidth=1.2)
            plt.setp(base_l, visible=False)
        else:
            pass
        if other_frag_df is not False:
            marker_l, stem_l, base_l = msms_pic.stem(other_frag_df['mz'], other_frag_df['i'], markerfmt=' ')
            plt.setp(stem_l, color=(0.8, 0.0, 0.0, 0.4), linewidth=1.2)
            plt.setp(base_l, visible=False)
            for _o_f_idx, _frag_se in other_frag_df.iterrows():
                _frag_mz = _frag_se['mz']
                _frag_i_r = _frag_se['i']
                _frag_class = _frag_se['LABEL']
                msms_pic.text(_frag_mz, _frag_i_r, _frag_class, txt_props,
                              fontsize=7, color=(0.8, 0.0, 0.0, 1), rotation=60)
        else:
            pass
        if other_nl_df is not False:
            marker_l, stem_l, base_l = msms_pic.stem(other_nl_df['mz'], other_nl_df['i'], markerfmt=' ')
            plt.setp(stem_l, color=(0.8, 0.0, 0.0, 0.4), linewidth=1.2)
            plt.setp(base_l, visible=False)
            for _o_nl_idx, _nl_se in other_nl_df.iterrows():
                _nl_mz = _nl_se['mz']
                _nl_i_r = _nl_se['i']
                _nl_class = _nl_se['LABEL']
                # _nl_i_r = _nl_i_r * 1.2
                msms_pic.text(_nl_mz, _nl_i_r, _nl_class, txt_props,
                              fontsize=7, color=(0.8, 0.0, 0.0, 1), rotation=60)
        else:
            pass
        if plt_obs_fa_df is not False:
            marker_l, stem_l, base_l = msms_pic.stem(plt_obs_fa_df['mz'], plt_obs_fa_df['i'], markerfmt=' ')
            plt.setp(stem_l, color=(0.8, 0.0, 0.0, 0.4), linewidth=1.2)
            plt.setp(base_l, visible=False)
        else:
            pass
        if plt_obs_lyso_df is not False:
            marker_l, stem_l, base_l = msms_pic.stem(plt_obs_lyso_df['mz'], plt_obs_lyso_df['i'], markerfmt=' ')
            plt.setp(stem_l, color=(0.8, 0.0, 0.0, 0.4), linewidth=1.2)
            plt.setp(base_l, visible=False)
        else:
            pass
        if target_frag_df is not False:
            marker_l, stem_l, base_l = msms_pic.stem(target_frag_df['mz'], target_frag_df['i'], markerfmt=' ')
            plt.setp(stem_l, color=(0.0, 0.45, 1.0, 0.6), linewidth=1.2, alpha=0.7)
            plt.setp(base_l, visible=False)
            for _t_f_idx, _frag_se in target_frag_df.iterrows():
                _frag_mz = _frag_se['mz']
                _frag_i_r = _frag_se['i']
                _frag_class = _frag_se['LABEL']
                msms_pic.text(_frag_mz, _frag_i_r, _frag_class, txt_props,
                              fontsize=8, color=(0.0, 0.5, 1.0, 1.0), weight='bold', rotation=60)
        else:
            pass
        if target_nl_df is not False:
            marker_l, stem_l, base_l = msms_pic.stem(target_nl_df['mz'], target_nl_df['i'], markerfmt=' ')
            plt.setp(stem_l, color=(0.0, 0.45, 1.0, 0.6), linewidth=1.2)
            plt.setp(base_l, visible=False)
            for _t_nl_idx, _nl_se in target_nl_df.iterrows():
                _nl_mz = _nl_se['mz']
                _nl_i_r = _nl_se['i']
                _nl_class = _nl_se['LABEL']
                # _nl_i_r = _nl_i * 1.2
                msms_pic.text(_nl_mz, _nl_i_r, _nl_class, txt_props,
                              fontsize=8, color=(0.0, 0.5, 1.0, 1.0), weight='bold', rotation=60)
        else:
            pass

        msms_title_str = ('MS/MS for m/z %.4f | DDA rank %d @ %.3f min' % (ms2_pr_mz, func_id, ms2_rt))
        msms_pic.set_title(msms_title_str, color=(0.0, 0.65, 1.0, 1.0), fontsize=8, y=0.98)

        # print(core_count, 'plot FULL MSMS in ', time.time() - _t_img_0)

    def plot_msms_low():
        # _t_img_0 = time.time()
        # print(core_count, 'start to plot MS/MS <= m/z 400 ...')
        msms_low_pic = pic_array[1, 1]
        msms_low_pic.tick_params(axis='both', which='major', labelsize=10)
        msms_low_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        msms_low_pic.set_xlabel("m/z", fontsize=7, labelpad=-1)
        msms_low_pic.set_ylabel("Intensity", fontsize=7)

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
            marker_l, stem_l, base_l = msms_low_pic.stem(_msms_low_df['mz'].values.tolist(),
                                                         _msms_low_df['i'].values.tolist(), markerfmt=' ')
            plt.setp(stem_l, color='grey', linewidth=0.6)
            plt.setp(base_l, visible=False)
            msms_low_pic.set_xlim([min(_msms_low_df['mz'].values.tolist()) - 1, 400])
            msms_low_pic.set_ylim([0, max(_msms_low_df['i'].values.tolist()) * 2])

        else:
            msms_low_pic.set_xlim([100, 400])
            msms_low_pic.set_ylim([0, ms2_df['i'].max()])

        if isinstance(obs_ident_df, pd.DataFrame):
            low_obs_ident_df = obs_ident_df[obs_ident_df['mz'] <= 400]
            if not low_obs_ident_df.empty:
                low_obs_ident_df.is_copy = False
                marker_l, stem_l, base_l = msms_low_pic.stem(low_obs_ident_df['mz'],
                                                             low_obs_ident_df['i'], markerfmt=' ')
                plt.setp(stem_l, color=(0, 0.7, 1.0, 0.6), linewidth=1.2)
                plt.setp(base_l, visible=False)
                for _i_idx, _ident_se in low_obs_ident_df.iterrows():
                    _ident_mz = _ident_se['mz']
                    _ident_i_r = _ident_se['i']
                    msms_low_pic.text(_ident_mz, _ident_i_r, _ident_se['obs_label'], txt_props,
                                      fontsize=8, color=(0, 0.6, 1.0, 1.0), weight='bold', rotation=60)
            else:
                pass
        else:
            pass

        if isinstance(plt_obs_fa_df, pd.DataFrame):
            low_plt_obs_fa_df = plt_obs_fa_df[plt_obs_fa_df['mz'] <= 400]
            if not low_plt_obs_fa_df.empty:
                low_plt_obs_fa_df.is_copy = False
                marker_l, stem_l, base_l = msms_low_pic.stem(low_plt_obs_fa_df['mz'], low_plt_obs_fa_df['i'],
                                                             markerfmt=' ')
                plt.setp(stem_l, color=(0.8, 0.0, 0.0, 0.4), linewidth=1.2)
                plt.setp(base_l, visible=False)
                for _i_f_idx, _frag_se in low_plt_obs_fa_df.iterrows():
                    _frag_mz = _frag_se['mz']
                    _frag_i_r = _frag_se['i']
                    msms_low_pic.text(_frag_mz, _frag_i_r, _frag_se['obs_label'], txt_props, fontsize=6,
                                      color=(0.8, 0.0, 0.0, 1), rotation=60)
        else:
            pass

        if isinstance(plt_obs_mg_df, pd.DataFrame):
            low_plt_obs_mg_df = plt_obs_mg_df[plt_obs_mg_df['mz'] <= 400]
            if not low_plt_obs_mg_df.empty:
                low_plt_obs_mg_df.is_copy = False
                marker_l, stem_l, base_l = msms_low_pic.stem(low_plt_obs_mg_df['mz'], low_plt_obs_mg_df['i'],
                                                             markerfmt=' ')
                plt.setp(stem_l, color=(0.8, 0.0, 0.0, 0.4), linewidth=1.2, alpha=0.4)
                plt.setp(base_l, visible=False)
                for _i_f_idx, _frag_se in low_plt_obs_mg_df.iterrows():
                    _frag_mz = _frag_se['mz']
                    _frag_i_r = _frag_se['i']
                    msms_low_pic.text(_frag_mz, _frag_i_r, _frag_se['obs_label'], txt_props,
                                      fontsize=6, color=(0.8, 0.0, 0.0, 1), rotation=60)
        else:
            pass
        # add specific ion info
        if isinstance(other_frag_df, pd.DataFrame):
            low_other_frag_df = other_frag_df[other_frag_df['mz'] <= 400]
            if not low_other_frag_df.empty:
                low_other_frag_df.is_copy = False
                marker_l, stem_l, base_l = msms_low_pic.stem(low_other_frag_df['mz'], low_other_frag_df['i'],
                                                             markerfmt=' ')
                plt.setp(stem_l, color=(0.8, 0.0, 0.0, 0.6), linewidth=1.2)
                plt.setp(base_l, visible=False)
                for _o_f_idx, _frag_se in low_other_frag_df.iterrows():
                    _frag_mz = _frag_se['mz']
                    _frag_i_r = _frag_se['i']
                    _frag_class = _frag_se['LABEL']
                    msms_low_pic.text(_frag_mz, _frag_i_r, _frag_class, txt_props,
                                      fontsize=7, color=(0.8, 0.0, 0.0, 1), rotation=60)
        else:
            pass

        if isinstance(target_frag_df, pd.DataFrame):
            marker_l, stem_l, base_l = msms_low_pic.stem(target_frag_df['mz'], target_frag_df['i'], markerfmt=' ')
            plt.setp(stem_l, color=(0.0, 0.45, 1.0, 0.6), linewidth=1.2)
            plt.setp(base_l, visible=False)
            for _t_f_idx, _frag_se in target_frag_df.iterrows():
                _frag_mz = _frag_se['mz']
                _frag_i_r = _frag_se['i']
                _frag_class = _frag_se['LABEL']
                msms_low_pic.text(_frag_mz, _frag_i_r, _frag_class, txt_props,
                                  fontsize=8, color=(0.0, 0.5, 1.0, 1.0), weight='bold', rotation=60)
        else:
            pass

        # Check missing part from higher m/z
        if isinstance(plt_obs_lyso_df, pd.DataFrame):
            low_plt_obs_lyso_df = plt_obs_lyso_df[plt_obs_lyso_df['mz'] <= 400]
            if not low_plt_obs_lyso_df.empty:
                low_plt_obs_lyso_df.is_copy = False
                marker_l, stem_l, base_l = msms_low_pic.stem(low_plt_obs_lyso_df['mz'], low_plt_obs_lyso_df['i'],
                                                             markerfmt=' ')
                plt.setp(stem_l, color=(0.8, 0.0, 0.0, 0.4), linewidth=1.2, alpha=0.4)
                plt.setp(base_l, visible=False)
                for _i_nl_idx, _nl_se in low_plt_obs_lyso_df.iterrows():
                    _nl_mz = _nl_se['mz']
                    _nl_i_r = _nl_se['i']
                    msms_low_pic.text(_nl_mz, _nl_i_r, _nl_se['obs_label'], txt_props, fontsize=6,
                                      color=(0.8, 0.0, 0.0, 1), rotation=60)
        else:
            pass

        if isinstance(plt_obs_dg_na_df, pd.DataFrame):
            low_plt_obs_dg_na_df = plt_obs_dg_na_df[plt_obs_dg_na_df['mz'] <= 400]
            if not low_plt_obs_dg_na_df.empty:
                low_plt_obs_dg_na_df.is_copy = False
                marker_l, stem_l, base_l = msms_low_pic.stem(low_plt_obs_dg_na_df['mz'], low_plt_obs_dg_na_df['i'],
                                                             markerfmt=' ')
                plt.setp(stem_l, color=(0.8, 0.0, 0.0, 0.4), linewidth=1.2, alpha=0.4)
                plt.setp(base_l, visible=False)
                for _i_nl_idx, _nl_se in low_plt_obs_dg_na_df.iterrows():
                    _nl_mz = _nl_se['mz']
                    _nl_i_r = _nl_se['i']
                    msms_low_pic.text(_nl_mz, _nl_i_r, _nl_se['obs_label'], txt_props, fontsize=6,
                                      color=(0.8, 0.0, 0.0, 1), rotation=60)
        else:
            pass

        if isinstance(other_nl_df, pd.DataFrame):
            low_other_nl_df = other_nl_df[other_nl_df['mz'] <= 400]
            if not low_other_nl_df.empty:
                low_other_nl_df.is_copy = False
                marker_l, stem_l, base_l = msms_low_pic.stem(low_other_nl_df['mz'], low_other_nl_df['i'],
                                                             markerfmt=' ')
                plt.setp(stem_l, color=(0.8, 0.0, 0.0, 0.6), linewidth=1.2)
                plt.setp(base_l, visible=False)
                for _o_nl_idx, _nl_se in low_other_nl_df.iterrows():
                    _nl_mz = _nl_se['mz']
                    _nl_i_r = _nl_se['i']
                    _nl_class = _nl_se['LABEL']
                    # _nl_i_r = _nl_i_r * 1.2
                    msms_low_pic.text(_nl_mz, _nl_i_r, _nl_class, txt_props,
                                      fontsize=7, color=(0.8, 0.0, 0.0, 1), rotation=60)
        else:
            pass

        # msms_low_pic.set_ylim([0, _msms_max * 1.5])
        msms_low_str = 'MS/MS zoomed below m/z 400'
        msms_low_pic.set_title(msms_low_str, color=(0.0, 0.65, 1.0, 1.0), fontsize=8, y=0.98)

        # print(core_count, 'plot MSMS <= 400 in ', time.time() - _t_img_0)

    def plot_msms_high():
        # _t_img_0 = time.time()
        # print(core_count, 'start to plot MS/MS > m/z 400 ...')
        msms_high_pic = pic_array[2, 1]
        msms_high_pic.tick_params(axis='both', which='major', labelsize=10)
        msms_high_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        msms_high_pic.set_xlabel("m/z", fontsize=7, labelpad=-1)
        msms_high_pic.set_ylabel("Intensity", fontsize=7)
        msms_high_pic.set_xlim([400, ms2_pr_mz + 20])

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
                                                      bbox=[0.4, 1 - (0.067 * (len(_lyso_table_vals) + 1)), 0.5,
                                                            0.067 * (len(_lyso_table_vals) + 1)])
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
            marker_l, stem_l, base_l = msms_high_pic.stem(_msms_high_df['mz'].values.tolist(),
                                                          _msms_high_df['i'].values.tolist(), markerfmt=' ')
            plt.setp(stem_l, color='grey', linewidth=0.6)
            plt.setp(base_l, visible=False)
            msms_high_pic.set_ylim([0, msms_high_max * 1.3])
        else:
            msms_high_max = ms2_df['i'].max()
            msms_high_pic.set_ylim([0, msms_high_max * 1.25])

        if isinstance(obs_ident_df, pd.DataFrame):
            high_obs_ident_df = obs_ident_df[obs_ident_df['mz'] > 400]
            if not high_obs_ident_df.empty:
                high_obs_ident_df.is_copy = False
                # high_obs_ident_df['i'] = high_obs_ident_df['i'] * 1.025
                marker_l, stem_l, base_l = msms_high_pic.stem(high_obs_ident_df['mz'],
                                                              high_obs_ident_df['i'], markerfmt=' ')
                plt.setp(stem_l, color=(0, 0.7, 1.0, 0.6), linewidth=1.2)
                plt.setp(base_l, visible=False)
                for _i_idx, _ident_se in high_obs_ident_df.iterrows():
                    _ident_mz = _ident_se['mz']
                    _ident_i_r = _ident_se['i']
                    msms_high_pic.text(_ident_mz, _ident_i_r, _ident_se['obs_label'], txt_props,
                                       fontsize=8, color=(0, 0.6, 1.0, 1.0), weight='bold', rotation=60)
            else:
                pass
        else:
            pass

        if isinstance(plt_obs_lyso_df, pd.DataFrame):
            high_plt_obs_lyso_df = plt_obs_lyso_df[plt_obs_lyso_df['mz'] > 400]
            if not high_plt_obs_lyso_df.empty:
                high_plt_obs_lyso_df.is_copy = False
                marker_l, stem_l, base_l = msms_high_pic.stem(high_plt_obs_lyso_df['mz'], high_plt_obs_lyso_df['i'],
                                                              markerfmt=' ')
                plt.setp(stem_l, color=(0.8, 0.0, 0.0, 0.4), linewidth=1.2)
                plt.setp(base_l, visible=False)
                for _i_nl_idx, _nl_se in high_plt_obs_lyso_df.iterrows():
                    _nl_mz = _nl_se['mz']
                    _nl_i_r = _nl_se['i']
                    msms_high_pic.text(_nl_mz, _nl_i_r, _nl_se['obs_label'], txt_props,
                                       fontsize=6, color=(0.8, 0.0, 0.0, 1), rotation=60)
        else:
            pass

        if isinstance(plt_obs_dg_na_df, pd.DataFrame):
            high_plt_obs_dg_na_df = plt_obs_dg_na_df[plt_obs_dg_na_df['mz'] > 400]
            if not high_plt_obs_dg_na_df.empty:
                high_plt_obs_dg_na_df.is_copy = False
                marker_l, stem_l, base_l = msms_high_pic.stem(high_plt_obs_dg_na_df['mz'], high_plt_obs_dg_na_df['i'],
                                                              markerfmt=' ')
                plt.setp(stem_l, color=(0.8, 0.0, 0.0, 0.4), linewidth=1.2)
                plt.setp(base_l, visible=False)
                for _i_nl_idx, _nl_se in high_plt_obs_dg_na_df.iterrows():
                    _nl_mz = _nl_se['mz']
                    _nl_i_r = _nl_se['i']
                    msms_high_pic.text(_nl_mz, _nl_i_r, _nl_se['obs_label'], txt_props,
                                       fontsize=6, color=(0.8, 0.0, 0.0, 1), rotation=60)
        else:
            pass

        if isinstance(other_nl_df, pd.DataFrame):
            high_other_nl_df = other_nl_df[other_nl_df['mz'] > 400]
            if not high_other_nl_df.empty:
                high_other_nl_df.is_copy = False
                marker_l, stem_l, base_l = msms_high_pic.stem(high_other_nl_df['mz'], high_other_nl_df['i'],
                                                              markerfmt=' ')
                plt.setp(stem_l, color=(0.8, 0.0, 0.0, 0.6), linewidth=1.2)
                plt.setp(base_l, visible=False)
                for _o_nl_idx, _nl_se in high_other_nl_df.iterrows():
                    _nl_mz = _nl_se['mz']
                    _nl_i_r = _nl_se['i']
                    _nl_class = _nl_se['LABEL']
                    # _nl_i_r = _nl_i_r * 1.2
                    msms_high_pic.text(_nl_mz, _nl_i_r, _nl_class, txt_props,
                                       fontsize=7, color=(0.8, 0.0, 0.0, 1), rotation=60)
        else:
            pass
        # add specific ion info
        if isinstance(target_nl_df, pd.DataFrame):
            marker_l, stem_l, base_l = msms_high_pic.stem(target_nl_df['mz'], target_nl_df['i'], markerfmt=' ')
            plt.setp(stem_l, color=(0.0, 0.45, 1.0, 0.6), linewidth=1.2)
            plt.setp(base_l, visible=False)
            for _t_nl_idx, _nl_se in target_nl_df.iterrows():
                _nl_mz = _nl_se['mz']
                _nl_i_r = _nl_se['i']
                _nl_class = _nl_se['LABEL']
                # _nl_i_r = _nl_i * 1.2
                msms_high_pic.text(_nl_mz, _nl_i_r, _nl_class, txt_props,
                                   fontsize=8, color=(0.0, 0.5, 1.0, 1.0), weight='bold', rotation=60)
        else:
            pass

        # Check missing part from lower m/z

        if isinstance(plt_obs_fa_df, pd.DataFrame):
            high_plt_obs_fa_df = plt_obs_fa_df[plt_obs_fa_df['mz'] > 400]
            if not high_plt_obs_fa_df.empty:
                high_plt_obs_fa_df.is_copy = False
                marker_l, stem_l, base_l = msms_high_pic.stem(high_plt_obs_fa_df['mz'], high_plt_obs_fa_df['i'],
                                                              markerfmt=' ')
                plt.setp(stem_l, color=(0.8, 0.0, 0.0, 0.4), linewidth=1.2, alpha=0.4)
                plt.setp(base_l, visible=False)
                for _i_f_idx, _frag_se in high_plt_obs_fa_df.iterrows():
                    _frag_mz = _frag_se['mz']
                    _frag_i_r = _frag_se['i']
                    msms_high_pic.text(_frag_mz, _frag_i_r, _frag_se['obs_label'], txt_props,
                                       fontsize=6, color=(0.8, 0.0, 0.0, 1), rotation=60)
        else:
            pass

        if isinstance(plt_obs_mg_df, pd.DataFrame):
            high_plt_obs_mg_df = plt_obs_mg_df[plt_obs_mg_df['mz'] > 400]
            if not high_plt_obs_mg_df.empty:
                high_plt_obs_mg_df.is_copy = False
                marker_l, stem_l, base_l = msms_high_pic.stem(high_plt_obs_mg_df['mz'], high_plt_obs_mg_df['i'],
                                                              markerfmt=' ')
                plt.setp(stem_l, color=(0.8, 0.0, 0.0, 0.4), linewidth=1.2, alpha=0.4)
                plt.setp(base_l, visible=False)
                for _i_f_idx, _frag_se in high_plt_obs_mg_df.iterrows():
                    _frag_mz = _frag_se['mz']
                    _frag_i_r = _frag_se['i']
                    msms_high_pic.text(_frag_mz, _frag_i_r, _frag_se['obs_label'], txt_props,
                                       fontsize=6, color=(0.8, 0.0, 0.0, 1), rotation=60)
        else:
            pass

        if isinstance(other_frag_df, pd.DataFrame):
            high_other_frag_df = other_frag_df[other_frag_df['mz'] > 400]
            if not high_other_frag_df.empty:
                high_other_frag_df.is_copy = False
                marker_l, stem_l, base_l = msms_high_pic.stem(high_other_frag_df['mz'], high_other_frag_df['i'],
                                                              markerfmt=' ')
                plt.setp(stem_l, color=(0.8, 0.0, 0.0, 0.6), linewidth=1.2)
                plt.setp(base_l, visible=False)
                for _o_f_idx, _frag_se in high_other_frag_df.iterrows():
                    _frag_mz = _frag_se['mz']
                    _frag_i_r = _frag_se['i']
                    _frag_class = _frag_se['LABEL']
                    msms_high_pic.text(_frag_mz, _frag_i_r, _frag_class, txt_props,
                                       fontsize=7, color=(0.8, 0.0, 0.0, 1), rotation=60)
        else:
            pass

        msms_high_str = 'MS/MS zoomed above m/z 400'
        msms_high_pic.set_title(msms_high_str, color=(0.0, 0.65, 1.0, 1.0), fontsize=8, y=0.98)
        # print(core_count, 'plot MSMS > 400 ', time.time() - _t_img_0)

    # all individual sub plot func finished
    # start to generate images

    try:
        tasks = [plot_msms(), plot_msms_low(), plot_msms_high(), plot_xic(), plot_ms(), plot_ms_zoom()]
        with ThreadPoolExecutor(max_workers=2) as executor:
            for _task in tasks:
                executor.submit(_task)

        plt.savefig(save_img_as, type=img_type, dpi=dpi)
        print(core_count, '[OUTPUT] ==> Image saved as: %s' % save_img_as)
        plt.close()
    except Exception as e:
        print('[INFO] --> Use single thread and try again ...', e)
        plot_msms()
        plot_msms_low()
        plot_msms_high()
        plot_xic()
        plot_ms()
        plot_ms_zoom()
        plt.savefig(save_img_as, type=img_type, dpi=dpi)
        print(core_count, '[INFO] --> Image saved as: %s' % save_img_as)
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
            print('%s [STATUS] >>> image: %i / %i' % (core_count, img_counter, tot_img_count))

            try:
                plot_spectra(abbr, mz_se, xic_dct, ident_info_dct, spec_info_dct, isotope_score_info_dct, specific_dct,
                             formula_charged, charge, core_count, save_img_as=save_img_as, img_type=img_type,
                             dpi=dpi, vendor=vendor, ms1_precision=ms1_precision)
            except Exception as e:
                print(core_count, '[EXCEPTION] !!! gen_plot failed to save images from data list ...', e)

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
                print(core_count, '[EXCEPTION] !!! gen_plot failed to save image ...', e)
