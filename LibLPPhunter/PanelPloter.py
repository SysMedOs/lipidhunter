# -*- coding: utf-8 -*-
# Copyright 2015-2016 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

from __future__ import division
import os

from matplotlib import pyplot as plt
import matplotlib as mpl
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from matplotlib._png import read_png
import pandas as pd

from StructureScore import AssignStructure


def plot_spectra(mz_se, xic_dct, ident_info_dct, ms1_rt, ms2_rt, ms1_df, ms2_df,
                 target_frag_df, target_nl_df, other_frag_df, other_nl_df, save_img_as=None):
    pr_mz = mz_se['mz']
    ms1_obs = mz_se['MS1_obs_mz']
    # _usr_rt = mz_se['rt']
    # _usr_abbr_bulk = mz_se['Abbreviation']
    # _usr_pl_class = mz_se['Class']
    # _usr_formula = mz_se['Formula']
    # _usr_ms2_function = mz_se['function']
    # _usr_ms2_scan_id = mz_se['scan_id']

    print ('Start looking for m/z', pr_mz)

    abbr_id = mz_se['Abbreviation']
    _mzlib_id = mz_se['Lib_mz']
    _func_id = mz_se['function']
    # _scan_id = mz_se['scan_id']

    _ms1_pr_df = ms1_df[ms1_df['mz'] == ms1_obs]
    # _ms1_pr_i = 0
    _ppm = 0.0
    # check if ms1_obs is present in MS survey scan with abs i >= 1000
    if _ms1_pr_df.shape[0] > 0:
        _ms1_pr_i = _ms1_pr_df['i'].tolist()[0]
        _ms1_pr_mz = _ms1_pr_df['mz'].tolist()[0]
        _mz_zoom_query_str = ' %.2f < mz < %.2f' % (ms1_obs - 2.1, ms1_obs + 2.1)
        _ms_zoom_df = ms1_df.query(_mz_zoom_query_str)
        _ms_zoom_df = _ms_zoom_df.sort_values(by='i', ascending=False)

        _mz_isotope_query_str = ' %.2f < mz < %.2f' % (ms1_obs - 1.5, ms1_obs - 0.1)
        _isotope_df = ms1_df.query(_mz_isotope_query_str)
        _isotope_df = _isotope_df.sort_values(by='i', ascending=False)

        print('_ms1_pr_df shape', _ms1_pr_df.shape)
        print(_ms1_pr_df)

        # check if i of ms_obs is higher than ms_obs - 1.5
        if _isotope_df.shape[0] == 0 or _ms1_pr_df.get_value(_ms1_pr_df.index[0],
                                                             'i') >= _isotope_df.get_value(_isotope_df.index[0], 'i'):
            print('>>> Isotope checker === ==> PASSED >>> >>> >>>')

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
            ms_pic.stem(ms1_df['mz'].tolist(), ms1_df['i'].tolist(), 'black', markerfmt=" ")
            # markerline, stemlines, baseline =
            # plt.setp(stemlines, color='grey')
            _dash_i = [ms1_df['i'].max()]
            # print(_dash_i)
            # print(_ms1_pr_df['mz'])

            markerline, stemlines, baseline = ms_pic.stem(_ms1_pr_df['mz'], _dash_i, ':', markerfmt=" ")
            plt.setp(stemlines, color='cyan', alpha=0.5)
            markerline, stemlines, baseline = ms_pic.stem(_ms1_pr_df['mz'], _ms1_pr_df['i'], markerfmt=" ")
            plt.setp(stemlines, color='cyan')

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

            markerline, stemlines, baseline = ms_zoom_pic.stem(_ms1_pr_df['mz'], _ms1_pr_df['i'], markerfmt=" ")
            # markerline, stemlines, baseline
            plt.setp(stemlines, color='cyan', linewidth=12, alpha=0.3)

            ms_zoom_pic.stem(_ms_zoom_df['mz'].tolist(), _ms_zoom_df['i'].tolist(), 'black', markerfmt=" ")
            # # markerline, stemlines, baseline
            # plt.setp(stemlines, color='grey')

            ms_zoom_pic.set_xlim([pr_mz - 2.1, pr_mz + 2.1])
            ms_zoom_pic.set_ylim([0, max(_ms_zoom_df['i'].tolist()) * 1.3])
            ms_zoom_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0), fontsize=10)
            ms_zoom_pic.ticklabel_format(axis='x', useOffset=False, fontsize=10)
            ms_zoom_pic.set_xlabel("m/z", fontsize=10, labelpad=-1)
            ms_zoom_pic.set_ylabel("Intensity", fontsize=10)

            # add annotation
            # _ms_zoom_peak_list = zip(_ms_zoom_df['mz'].tolist(), _ms_zoom_df['i'].tolist())
            # print('_ms_zoom_df')
            # print(_ms_zoom_df)

            _ms_zoom_top_df = _ms_zoom_df.sort_values(by='i', ascending=False).head(10)
            _ms_zoom_peak_list = zip(_ms_zoom_top_df['mz'].tolist(),
                                     _ms_zoom_top_df['i'].tolist())

            for _ms_zoom_peak in _ms_zoom_peak_list:
                if _ms_zoom_peak[0] != _ms1_pr_mz:
                    _ms_zoom_peak_str = '%.4f' % _ms_zoom_peak[0]
                    _ms_zoom_peak_y = _ms_zoom_peak[1]
                    ms_zoom_pic.text(_ms_zoom_peak[0], _ms_zoom_peak_y, _ms_zoom_peak_str,
                                     fontsize=6)

            markerline, stemlines, baseline = ms_zoom_pic.stem(_ms1_pr_df['mz'], _ms1_pr_df['i'],
                                                               'magenta', markerfmt="D"
                                                               )
            plt.setp(markerline, markerfacecolor='magenta', markersize=5, markeredgewidth=0)
            ms_zoom_pic.text(_ms1_pr_df['mz'], _ms1_pr_df['i'], '%.4f' % float(_ms1_pr_df['mz']),
                             color='magenta', fontsize=6
                             )

            # XIC spectrum start
            xic_rt_lst = xic_df['rt'].tolist()
            xic_i_lst = xic_df['i'].tolist()
            xic_pic.plot(xic_rt_lst, xic_i_lst, alpha=0.3)
            xic_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            xic_pic.stem([ms1_rt], [max(xic_i_lst)], ':', markerfmt=" ")
            xic_pic.stem([ms2_rt], [max(xic_i_lst)], ':', markerfmt=" ")
            xic_pic.text(ms1_rt - 0.3, max(xic_i_lst) * 0.98, 'MS', fontsize=10)
            xic_pic.text(ms2_rt, max(xic_i_lst) * 0.98, 'MS/MS', fontsize=10)
            xic_pic.set_xlabel("RT (min)", fontsize=10, labelpad=-1)
            xic_pic.set_ylabel("Intensity", fontsize=10)

            # msms spectrum start
            msms_pic.stem(ms2_df['mz'].tolist(), ms2_df['i'].tolist(), 'black', markerfmt=" ")
            msms_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            msms_pic.set_xlabel("m/z", fontsize=10, labelpad=-1)
            msms_pic.set_ylabel("Intensity", fontsize=10)
            msms_pic.set_xlim([min(ms2_df['mz'].tolist()) - 1, pr_mz + 20])
            msms_pic.set_ylim([0, max(ms2_df['i'].tolist()) * 1.3])

            # # add annotations
            # _topms2_df = ms2_df.sort_values(by='i', ascending=False)
            # _topms2_df = _topms2_df.head(20)
            # _msms_peak_list = zip(_topms2_df['mz'].tolist(), _topms2_df['i'].tolist())
            # for _msms_peak in _msms_peak_list:
            #     _msms_peak_str = '%.4f' % _msms_peak[0]
            #     _msms_peak_y = (max(ms2_df['i'].tolist()) * 1.3) * 0.175 + _msms_peak[1]
            #     msms_pic.text(_msms_peak[0], _msms_peak_y, _msms_peak_str,
            #                   rotation=90, fontsize=6)

            # prepare DataFrame for msms zoomed plot
            _msms_low_df = ms2_df.query('mz <= 350')
            _msms_high_df = ms2_df.query('mz > 350')

            _ident_table_df = ident_info_dct['SCORE_INFO']
            _fa_table_df = ident_info_dct['FA_INFO']
            _lyso_table_df = ident_info_dct['LYSO_INFO']

            if _ident_table_df.shape[0] > 0:
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
                    markerline, stemlines, baseline = msms_low_pic.stem([_fa_se['mz']], [_fa_se['i']],
                                                                        markerfmt="D"
                                                                        )
                    # color of stemlines is tuple of R, G, B, alpha from 0 to 1
                    _rgb_color = ((20*_i_fa)/255, (255-5*_i_fa)/255, (255-5*_i_fa)/255, 0.6-0.02*_i_fa)
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
                    markerline, stemlines, baseline = msms_high_pic.stem([_lyso_se['mz']], [_lyso_se['i']],
                                                                         markerfmt="D"
                                                                         )
                    # color of stemlines is tuple of R, G, B, alpha from 0 to 1
                    _rgb_color = ((20*_i_lyso)/255, (255-5*_i_lyso)/255, (255-5*_i_lyso)/255, 0.6-0.02*_i_lyso)
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

            # add specific ion info
            if other_frag_df.shape[0] > 0:
                for _idx, _frag_se in other_frag_df.iterrows():
                    _frag_mz = _frag_se['mz']
                    _frag_i = _frag_se['i']
                    _frag_class = _frag_se['CLASS']
                    _frag_i_x = min(_frag_i * 5, 0.3 * max(_msms_low_df['i'].tolist()))
                    _frag_i = sorted([max(_msms_low_df['i'].tolist()) * 1.1, _frag_i * 1.1, _frag_i_x])[1]
                    markerline, stemlines, baseline = msms_low_pic.stem([_frag_mz], [_frag_i], markerfmt=" ")
                    plt.setp(stemlines, color='orange', linewidth=3, alpha=0.4)
                    msms_low_pic.text(_frag_mz, _frag_i, _frag_class, fontsize=10, color='orange')

            if other_nl_df.shape[0] > 0:
                for _idx, _nl_se in other_nl_df.iterrows():
                    _nl_mz = _nl_se['mz']
                    _nl_i = _nl_se['i']
                    _nl_class = _nl_se['CLASS']
                    _nl_i_x = min(_nl_i * 5, 0.3 * max(_msms_high_df['i'].tolist()))
                    _nl_i = sorted([max(_msms_high_df['i'].tolist()) * 1.1, _nl_i * 1.1, _nl_i_x])[1]
                    markerline, stemlines, baseline = msms_high_pic.stem([_nl_mz], [_nl_i], markerfmt=" ")
                    plt.setp(stemlines, color='orange', linewidth=3, alpha=0.4)
                    msms_high_pic.text(_nl_mz, _nl_i, _nl_class, fontsize=10, color='orange')

            if target_frag_df.shape[0] > 0:
                for _idx, _frag_se in target_frag_df.iterrows():
                    _frag_mz = _frag_se['mz']
                    _frag_i = _frag_se['i']
                    _frag_class = _frag_se['CLASS']
                    _frag_i_x = min(_frag_i * 5, 0.5 * max(_msms_low_df['i'].tolist()))
                    _frag_i = sorted([max(_msms_low_df['i'].tolist()) * 1.1, _frag_i * 1.1, _frag_i_x])[1]
                    markerline, stemlines, baseline = msms_low_pic.stem([_frag_mz], [_frag_i], markerfmt=" ")
                    plt.setp(stemlines, color='green', linewidth=3, alpha=0.4)
                    msms_low_pic.text(_frag_mz, _frag_i, _frag_class, fontsize=10, color='green')

            if target_nl_df.shape[0] > 0:
                for _idx, _nl_se in target_nl_df.iterrows():
                    _nl_mz = _nl_se['mz']
                    _nl_i = _nl_se['i']
                    _nl_class = _nl_se['CLASS']
                    _nl_i_x = min(_nl_i * 5, 0.3 * max(_msms_low_df['i'].tolist()))
                    _nl_i = sorted([max(_msms_high_df['i'].tolist()) * 1.1, _nl_i * 1.1, _nl_i_x])[1]
                    markerline, stemlines, baseline = msms_high_pic.stem([_nl_mz], [_nl_i], markerfmt=" ")
                    # markerline, stemlines, baseline
                    plt.setp(stemlines, color='green', linewidth=3, alpha=0.4)
                    msms_high_pic.text(_nl_mz, _nl_i, _nl_class, fontsize=10, color='green')

            # msms spectrum zoomed > 350 start
            if _msms_high_df.shape[0] > 0:
                msms_high_pic.stem(_msms_high_df['mz'].tolist(),
                                   _msms_high_df['i'].tolist(),
                                   'black', lw=4, markerfmt=" ")
                msms_high_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
                msms_high_pic.set_xlabel("m/z", fontsize=10, labelpad=-1)
                msms_high_pic.set_ylabel("Intensity", fontsize=10)
                msms_high_pic.set_xlim([350, pr_mz + 20])
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

            # msms spectrum zoomed < 350 start
            if _msms_low_df.shape[0] > 0:
                msms_low_pic.stem(_msms_low_df['mz'].tolist(),
                                  _msms_low_df['i'].tolist(),
                                  'black', lw=2, markerfmt=" ")
                msms_low_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
                msms_low_pic.set_xlabel("m/z", fontsize=10, labelpad=-1)
                msms_low_pic.set_ylabel("Intensity", fontsize=10)
                msms_low_pic.set_xlim([min(_msms_low_df['mz'].tolist()) - 1, 350])
                msms_low_pic.set_ylim([0, max(_msms_low_df['i'].tolist()) * 1.3])
                # PC HG form DOI: 10.1021/acs.jproteome.5b00169
                # for _usr_mz in [168.0358, 224.0693, 242.0798]:
                #     _query_str = '%f - %f <= mz <= %f + %f' % (_usr_mz, 0.25, _usr_mz, 0.25)
                #     _snms2_df = ms2_df.query(_query_str)
                #     _snms2_df = _snms2_df.sort(columns='i', ascending=False)
                #     # print _snms2_df.head(2)
                #     try:
                #         msms_low_pic.stem([_snms2_df['mz'].tolist()[0]],
                #                           [_snms2_df['i'].tolist()[0]], 'red', lw=4, markerfmt=" ")
                #         _hg_checker += 1
                #     except:
                #         print('no sn fit')

                # # add annotations
                # _top_msms_low_df = _msms_low_df.sort_values(by='i', ascending=False)
                # _top_msms_low_df = _top_msms_low_df.head(20)
                # _msms_low_peak_list = zip(_top_msms_low_df['mz'].tolist(),
                #                           _top_msms_low_df['i'].tolist())
                # for _msms_low_peak in _msms_low_peak_list:
                #     _msms_low_peak_str = '%.4f' % _msms_low_peak[0]
                #     _msms_low_peak_y = ((max(_msms_low_df['i'].tolist()) * 1.3)
                #                         * 0.175 + _msms_low_peak[1]
                #                         )
                #     msms_low_pic.text(_msms_low_peak[0], _msms_low_peak_y, _msms_low_peak_str,
                #                       rotation=90, fontsize=6)
            else:
                pass

            # set title
            _ppm = 1e6 * (ms1_obs - _mzlib_id) / _mzlib_id
            xic_title_str = 'XIC of m/z %.4f | %s @ m/z %.4f ppm=%.2f' % (ms1_obs, abbr_id, _mzlib_id, _ppm)
            ms_title_str = 'MS @ %.3f min' % ms1_rt
            ms_zoom_title_str = 'MS zoomed'
            msms_title_str = ('MS/MS of m/z %.4f | DDA Top %d @ %.3f min' % (pr_mz, _func_id, ms2_rt))
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
            return _ms1_pr_i, _ppm, isotope_checker

        else:
            print ('!!!!!! Isotope checker !=!=!= FAILED !!!!!!')
            plt.close()
            print ('!!!!!! Not identified !!!!!! >>> >>> >>>')
            isotope_checker = 1
            return _ms1_pr_i, _ppm, isotope_checker
