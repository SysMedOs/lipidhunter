# -*- coding: utf-8 -*-
# Copyright 2014-16 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release.
# For more info please contact zhixu.ni@uni-leipzig.de

import os

from matplotlib import pyplot as plt
import matplotlib as mpl
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from matplotlib._png import read_png
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

from FAindicator import FAindicator
from FAindicator import Lyso_indicator
from StructureScore import AssignStructure


class Spectra_Ploter(object):
    def plot_all(self, mz1get_lst, mz2get_lst, sdf_df, xic_dct,
                 ms_spectra_dct, msms_spectra_dct, usr_fa_csv, Pl_type, path=None):
        """

        :param mz2get_lst: A list of m/z values to be searched.
        :type mz2get_lst: list[float]
        :param sdf_df: The dataframe contains the extracted info from ``sdf`` files
        :type sdf_df: pandas.DataFrame
        :param xic_dct:
        :type xic_dct: dict
        :param ms_spectra_dct:
        :type ms_spectra_dct: dict
        :param msms_spectra_dct:
        :type msms_spectra_dct: dict
        :param path: The path to save generated images.
        :type path: str
        :return:
        """

        if path is None:
            _path = os.getcwd()
        else:
            if os.path.exists(path):
                _path = os.path.abspath(path)
            else:
                _path = os.getcwd()

        # print _path

        _fa_checker = 0
        _hg_checker = 0

        fa_indicator = FAindicator(usr_fa_csv)
        lyso_indicator = Lyso_indicator(usr_fa_csv)

        for mz2get in mz2get_lst:
            print ('Start looking for m/z', mz2get)
            _idx_ms2 = mz2get_lst.index(mz2get)
            _ms = mz1get_lst[_idx_ms2]
            if _ms in xic_dct.keys() and xic_dct[_ms].shape[0] > 0:
                xic_df = xic_dct[_ms]

                # print ('xic_df')
                # print (xic_df)
                _ms_spectra_dct = ms_spectra_dct[_ms]
                _msms_spectra_dct = msms_spectra_dct[mz2get]
                # print ('_msms_spectra_dct', _msms_spectra_dct)

                _query_str = 'mz == %.6f' % mz2get
                _sdf_df = sdf_df.query(_query_str)
                print ('_sdf_df.shape', _sdf_df.shape)
                if _sdf_df.shape[0] == 0:
                    # _query_str = '%.6f <= mz <=%.6f' % (mz2get*0.999999, mz2get*1.000001)
                    # set ppm
                    _query_str = '%.6f <= mz <= %.6f' % (mz2get - 0.9, mz2get + 0.9)
                    _sdf_df = sdf_df.query(_query_str)
                    if _sdf_df.shape[0] > 0:
                        # print ('MS2 PR in 1 ppm')
                        print ('MS2 PR in +/- 0.9 m/z')
                        _sdf_df = sdf_df.query(_query_str)
                        print ('New_sdf_df.shape', _sdf_df.shape)
                    else:
                        print ('MS2 PR NOT found!!!!')
                else:
                    print ('MS2 PR just fit')
                # print ('_sdf_df', _sdf_df)

                for index, row in _sdf_df.iterrows():
                    _hmdb_id = row['Abbreviation']
                    _mzlib_id = row['Lib_mz']
                    _func_id = row['function']
                    _scan_id = row['scan_id']

                    # try:

                    xic_mz_lst = xic_df['rt'].tolist()
                    xic_i_lst = xic_df['i'].tolist()

                    ms_rt_lst = _ms_spectra_dct.keys()
                    ms_rt_lst.sort()
                    print(ms_rt_lst)
                    # print('_msms_spectra_dct', _msms_spectra_dct)
                    for _msms in _msms_spectra_dct.keys():
                        print ('_msms', _msms)
                        _pr_mz = _msms[0]
                        _pr_rt = _msms[1]
                        # print (_pr_mz, _pr_rt, ms_rt_lst)
                        _msms_df = _msms_spectra_dct[_msms]

                        _ms_rt_lst = []
                        for _ms_rt in ms_rt_lst:
                            if 0 < _pr_rt - _ms_rt < 0.05:  # MAX DDA top 12 x 0.25s = 3s = 0.05 min
                                _ms_rt_lst.append(_ms_rt)
                        _ms_rt_lst.sort()
                        print('_ms_rt_lst', _ms_rt_lst)
                        if len(_ms_rt_lst) > 0:
                            _ms_rt_pr = _ms_rt_lst[-1]

                            # print _ms_rt_lst
                            print('MS @ %.2f for DDA @ %.2f' % (_ms_rt_pr, _pr_rt))

                            # skip the image if MS and MS/MS time do not match
                            # if _pr_rt - _ms_rt_pr > 0.2:
                            #     print 'MS and MS/MS time NOT match.'
                            #     pass
                            #
                            # if len(_msms_df['mz'].tolist()) < 1 or len(_msms_df['i'].tolist())<1:
                            #     print 'No peaks in MS/MS'

                            if len(_msms_df['mz'].tolist()) > 1 or len(_msms_df['i'].tolist()) > 1:
                                print ('try to plot---------------------------------------------------------->')
                                _auto_ident_chker = 0
                                # _ms_pkl_lst = _ms_spectra_dct[_ms_rt_pr][1]
                                # _ms_pkl_df = pd.DataFrame(data=_ms_pkl_lst, columns=['mz', 'i'])

                                _ms_pkl_df = _ms_spectra_dct[_ms_rt_pr][1]
                                print(_ms_pkl_df.head())
                                _mz_zoom_query_str = ' %.2f < mz < %.2f' % (mz2get - 2.1, mz2get + 2.1)
                                _ms_zoom_df = _ms_pkl_df.query(_mz_zoom_query_str)
                                _ms_zoom_df = _ms_zoom_df.sort_values(by='i', ascending=False)

                                _max_ms_iso = _ms_zoom_df['mz'].tolist()[0]
                                if mz2get - 0.75 < _max_ms_iso < mz2get + 0.75:
                                    pass
                                else:
                                    print ('Isotopes--->>>--->>>--->>>pass>>>')
                                    break

                                # Generate A4 image in landscape
                                fig, pic_array = plt.subplots(nrows=3, ncols=2, figsize=(11.692, 8.267), sharex=False,
                                                              sharey=False)
                                # Make better spacing between subplots
                                plt.tight_layout()
                                # label_size = 8
                                # mpl.rcParams['xtick.labelsize'] = label_size
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

                                ms_pic.stem(_ms_pkl_df['mz'].tolist(), _ms_pkl_df['i'].tolist(),
                                            'blue', lw=4, markerfmt=" ")
                                _ms_pr_lst = _ms_spectra_dct[_ms_rt_pr][0]
                                _dash_i = [max(_ms_pkl_df['i'].tolist())] * len(_ms_pr_lst['i'])
                                ms_pic.stem(_ms_pr_lst['mz'], _dash_i, ':', 'yellow', markerfmt=" ")
                                ms_pic.stem(_ms_pr_lst['mz'], _ms_pr_lst['i'], 'red', lw=4, markerfmt=" ")

                                ms_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
                                ms_pic.set_xlabel("m/z", fontsize=10, labelpad=-1)
                                ms_pic.set_ylabel("Intensity", fontsize=10)
                                ms_pic.set_ylim([0, max(_ms_pkl_df['i'].tolist()) * 1.3])

                                # add annotation
                                _ms_pkl_top_df = _ms_pkl_df.sort_values(by='i', ascending=False).head(10)
                                _ms_pkl_top_peak_list = zip(_ms_pkl_top_df['mz'].tolist(), _ms_pkl_top_df['i'].tolist())
                                for _ms_pkl_top_peak in _ms_pkl_top_peak_list:
                                    _ms_pkl_top_peak_str = '%.4f' % _ms_pkl_top_peak[0]
                                    _ms_pkl_top_peak_y = ((max(_ms_pkl_top_df['i'].tolist()) * 1.3) * 0.175 +
                                                          _ms_pkl_top_peak[1])
                                    ms_pic.text(_ms_pkl_top_peak[0], _ms_pkl_top_peak_y, _ms_pkl_top_peak_str,
                                                rotation=90, fontsize=6)

                                ms_zoom_pic.stem(_ms_zoom_df['mz'].tolist(), _ms_zoom_df['i'].tolist(),
                                                 'black', lw=4, markerfmt=" ")
                                _ms_pr_lst = _ms_spectra_dct[_ms_rt_pr][0]
                                ms_zoom_pic.stem(_ms_pr_lst['mz'], _ms_pr_lst['i'], 'red', lw=4, markerfmt=" ")

                                ms_zoom_pic.set_xlim([mz2get - 2.1, mz2get + 2.1])
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
                                    _ms_zoom_peak_str = '%.4f' % _ms_zoom_peak[0]
                                    _ms_zoom_peak_y = (max(_ms_zoom_df['i'].tolist()) * 1.3) * 0.005 + _ms_zoom_peak[1]
                                    ms_zoom_pic.text(_ms_zoom_peak[0], _ms_zoom_peak_y, _ms_zoom_peak_str,
                                                     fontsize=6)

                                # XIC spectrum start
                                xic_pic.plot(xic_mz_lst, xic_i_lst, alpha=0.3)
                                xic_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
                                xic_pic.stem([_ms_rt_pr], [max(xic_i_lst)], ':', markerfmt=" ")
                                xic_pic.stem([_pr_rt], [max(xic_i_lst)], ':', markerfmt=" ")
                                xic_pic.text(_ms_rt_pr - 0.3, max(xic_i_lst) * 0.98, 'MS', fontsize=10)
                                xic_pic.text(_pr_rt, max(xic_i_lst) * 0.98, 'MS/MS', fontsize=10)
                                xic_pic.set_xlabel("RT (min)", fontsize=10, labelpad=-1)
                                xic_pic.set_ylabel("Intensity", fontsize=10)

                                # msms spectrum start
                                msms_pic.stem(_msms_df['mz'].tolist(), _msms_df['i'].tolist(), 'blue', lw=2,
                                              markerfmt=" ")
                                msms_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
                                msms_pic.set_xlabel("m/z", fontsize=10, labelpad=-1)
                                msms_pic.set_ylabel("Intensity", fontsize=10)
                                msms_pic.set_xlim([min(_msms_df['mz'].tolist()) - 1, _pr_mz + 20])
                                msms_pic.set_ylim([0, max(_msms_df['i'].tolist()) * 1.3])

                                # add annotations
                                _top_msms_df = _msms_df.sort_values(by='i', ascending=False)
                                _top_msms_df = _top_msms_df.head(20)
                                _msms_peak_list = zip(_top_msms_df['mz'].tolist(), _top_msms_df['i'].tolist())
                                for _msms_peak in _msms_peak_list:
                                    _msms_peak_str = '%.4f' % _msms_peak[0]
                                    _msms_peak_y = (max(_msms_df['i'].tolist()) * 1.3) * 0.175 + _msms_peak[1]
                                    msms_pic.text(_msms_peak[0], _msms_peak_y, _msms_peak_str,
                                                  rotation=90, fontsize=6)
                                # _pr_mol = Chem.MolFromSmiles(_pr_smi)
                                # AllChem.Compute2DCoords(_pr_mol)
                                #
                                # # arr_hand = read_png('/path/to/this/image.png')
                                # _pr_img = Draw.MolToImage(_pr_mol, size=(450, 300))
                                # imagebox = OffsetImage(_pr_img)
                                # xy = [max(_msms_df['mz'].tolist()) * 0.5, max(_msms_df['i'].tolist())]
                                # _pr_ab = AnnotationBbox(imagebox, xy, xycoords='data', xybox=(90., -60.),
                                #                         boxcoords="offset points", frameon=True)
                                # msms_pic.add_artist(_pr_ab)

                                # for _usr_mz in [_sn1_mz, _sn2_mz]:
                                #     _msms_top_df = _msms_df.sort(columns='i', ascending=False)
                                #     _msms_top_df = _msms_top_df.head(10)
                                #     _query_str = '%f - %f <= mz <= %f + %f' % (_usr_ms2_pr_mz, 0.25, _usr_ms2_pr_mz, 0.25)
                                #     _sn_msms_df = _msms_top_df.query(_query_str)
                                #     _sn_msms_df = _sn_msms_df.sort(columns='i', ascending=False)
                                #     # print _sn_msms_df.head(2)
                                #     try:
                                #         msms_pic.stem([_sn_msms_df['mz'].tolist()[0]], [_sn_msms_df['i'].tolist()[0]],
                                #                       'red', lw=4, markerfmt=" ")
                                #         _fa_checker += 1
                                #     except:
                                #         print('no sn fit')

                                # msms spectrum zoomed below 355 start
                                _msms_low_df = _msms_df.query('mz < 350')
                                if len(_msms_low_df['mz'].tolist()) > 0:
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
                                    #     _query_str = '%f - %f <= mz <= %f + %f' % (_usr_ms2_pr_mz, 0.25, _usr_ms2_pr_mz, 0.25)
                                    #     _sn_msms_df = _msms_df.query(_query_str)
                                    #     _sn_msms_df = _sn_msms_df.sort(columns='i', ascending=False)
                                    #     # print _sn_msms_df.head(2)
                                    #     try:
                                    #         msms_low_pic.stem([_sn_msms_df['mz'].tolist()[0]],
                                    #                           [_sn_msms_df['i'].tolist()[0]], 'red', lw=4, markerfmt=" ")
                                    #         _hg_checker += 1
                                    #     except:
                                    #         print('no sn fit')

                                    # add annotations
                                    _top_msms_low_df = _msms_low_df.sort_values(by='i', ascending=False)
                                    _top_msms_low_df = _top_msms_low_df.head(20)
                                    _msms_low_peak_list = zip(_top_msms_low_df['mz'].tolist(),
                                                              _top_msms_low_df['i'].tolist())
                                    for _msms_low_peak in _msms_low_peak_list:
                                        _msms_low_peak_str = '%.4f' % _msms_low_peak[0]
                                        _msms_low_peak_y = (max(_msms_low_df['i'].tolist()) * 1.3) * 0.175 + \
                                                           _msms_low_peak[1]
                                        msms_low_pic.text(_msms_low_peak[0], _msms_low_peak_y, _msms_low_peak_str,
                                                          rotation=90, fontsize=6)
                                else:
                                    pass

                                # get identification
                                (fa_df, fa_info_dct) = fa_indicator.indicate(_msms_df)

                                if fa_info_dct is not None:
                                    # print fa_info
                                    fa_count = len(fa_info_dct['mz'])
                                    fa_count_lst = range(fa_count)
                                    _fa_name_lst = fa_info_dct['name']
                                    _fa_mz_lst = fa_info_dct['mz']
                                    _fa_i_lst = fa_info_dct['i']
                                    _fa_ppm_lst = fa_info_dct['ppm']
                                    _fa_delta_lst = fa_info_dct['D']

                                    fa_info_df = pd.DataFrame(fa_info_dct, columns=['name', 'mz', 'i', 'ppm'])
                                    fa_info_df = fa_info_df.sort_values(by='i', ascending=False)
                                    fa_info_df = fa_info_df.round({'mz': 4})
                                    # fa_info_df.index = range(1, len(fa_info['mz']) + 1)
                                    table_info_df = fa_info_df

                                    # print table_info_df.head(5)
                                    old_idx_lst = table_info_df.index.tolist()
                                    table_info_df.index = range(1, len(old_idx_lst) + 1)

                                    col_labels = table_info_df.columns.tolist()
                                    # col_labels = mz_info_df.head().tolist()
                                    row_labels = table_info_df.index.tolist()
                                    table_vals = map(list, table_info_df.values)
                                    try:
                                        # the rectangle is where I want to place the table
                                        fa_table = msms_pic.table(cellText=table_vals, rowLabels=row_labels,
                                                                  colWidths=[.10] * len(col_labels),
                                                                  colLabels=col_labels, loc='upper center')
                                        fa_table.set_fontsize(5)
                                        _auto_ident_chker += 1
                                    except IndexError:
                                        pass
                                        # table_props = the_table.properties()
                                        # table_cells = table_props['child_artists']
                                        # for cell in table_cells:
                                        #     cell.set_height(0.12)

                                # msms spectrum zoomed above 350 start
                                _msms_high_df = _msms_df.query('mz > 350')
                                if len(_msms_high_df['mz'].tolist()) > 0:
                                    msms_high_pic.stem(_msms_high_df['mz'].tolist(),
                                                       _msms_high_df['i'].tolist(),
                                                       'black', lw=4, markerfmt=" ")
                                    msms_high_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
                                    msms_high_pic.set_xlabel("m/z", fontsize=10, labelpad=-1)
                                    msms_high_pic.set_ylabel("Intensity", fontsize=10)
                                    msms_high_pic.set_xlim([350, _pr_mz + 20])
                                    msms_high_pic.set_ylim([0, max(_msms_high_df['i'].tolist()) * 1.3])

                                    # add annotations
                                    _top_msms_high_df = _msms_high_df.sort_values(by='i', ascending=False)
                                    _top_msms_high_df = _top_msms_high_df.head(20)
                                    _msms_high_peak_list = zip(_top_msms_high_df['mz'].tolist(),
                                                               _top_msms_high_df['i'].tolist())
                                    for _msms_high_peak in _msms_high_peak_list:
                                        _msms_high_peak_str = '%.4f' % _msms_high_peak[0]
                                        _msms_high_peak_y = ((max(_msms_high_df['i'].tolist()) * 1.3) * 0.175 +
                                                             _msms_high_peak[1])
                                        msms_high_pic.text(_msms_high_peak[0], _msms_high_peak_y, _msms_high_peak_str,
                                                           rotation=90, fontsize=6)
                                else:
                                    pass

                                # get Lyso identification
                                (lyso_df, lyso_info_dct) = lyso_indicator.indicate(_msms_df, _ms, PLtype=Pl_type)

                                if lyso_info_dct is not None:
                                    # print lyso_info
                                    lyso_count = len(lyso_info_dct['mz'])
                                    lyso_count_lst = range(lyso_count)
                                    _lyso_name_lst = lyso_info_dct['name']
                                    _lyso_mz_lst = lyso_info_dct['mz']
                                    _lyso_i_lst = lyso_info_dct['i']
                                    _lyso_ppm_lst = lyso_info_dct['ppm']
                                    _lyso_delta_lst = lyso_info_dct['D']

                                    lyso_info_df = pd.DataFrame(lyso_info_dct, columns=['name', 'mz', 'i', 'ppm'])
                                    lyso_info_df = lyso_info_df.sort_values(by='i', ascending=False)
                                    lyso_info_df = lyso_info_df.round({'mz': 4})

                                    # print table_info_df.head(5)
                                    old_idx_lst = lyso_info_df.index.tolist()
                                    lyso_info_df.index = range(1, len(old_idx_lst) + 1)

                                    col_labels = lyso_info_df.columns.tolist()
                                    # col_labels = mz_info_df.head().tolist()
                                    row_labels = lyso_info_df.index.tolist()
                                    table_vals = map(list, lyso_info_df.values)
                                    # plot lyso table
                                    try:
                                        # the rectangle is where I want to place the table
                                        lyso_table = msms_high_pic.table(cellText=table_vals, rowLabels=row_labels,
                                                                         colWidths=[.22, .1, .1, .1],
                                                                         colLabels=col_labels, loc='upper center')
                                        lyso_table.set_fontsize(7)
                                        _auto_ident_chker += 1
                                    except IndexError:
                                        pass

                                    # get assignment
                                    ident_struct = AssignStructure()
                                    _match_fa_df = pd.DataFrame()
                                    usr_std_fa_df = pd.read_csv(usr_fa_csv)
                                    usr_std_fa_df['C'].astype(int)
                                    usr_std_fa_df['DB'].astype(int)
                                    for _ident_idx, _ident_row in usr_std_fa_df.iterrows():
                                        pre_ident_fa_df = usr_std_fa_df
                                        pre_ident_fa_df['abs'] = 0
                                        tmp_fa = pre_ident_fa_df.ix[_ident_idx]['FA']
                                        if tmp_fa in fa_info_dct['fa']:
                                            _tmp_fa_idx = fa_info_dct['fa'].index(tmp_fa)
                                            _tmp_i = fa_info_dct['abs'][_tmp_fa_idx]
                                            pre_ident_fa_df.set_value(_ident_idx, 'abs', _tmp_i)
                                            _match_fa_df = _match_fa_df.append(pre_ident_fa_df.ix[_ident_idx])
                                    _match_lyso_df = pd.DataFrame()
                                    for _ident_idx, _ident_row in usr_std_fa_df.iterrows():
                                        pre_lyso_fa_df = usr_std_fa_df
                                        pre_lyso_fa_df['type'] = ''
                                        pre_lyso_fa_df['abs'] = 0
                                        tmp_lyso = pre_lyso_fa_df.ix[_ident_idx]['FA']
                                        lyso_info_zip_lst = zip(lyso_info_dct['type'], lyso_info_dct['fa'])
                                        if ('Lyso-H2O', tmp_lyso) in lyso_info_zip_lst:
                                            # _tmp_lyso_idx = lyso_info_dct['fa'].index(tmp_lyso)
                                            # _tmp_type = lyso_info_dct['type'][_tmp_lyso_idx]
                                            _tmp_lyso_idx = lyso_info_zip_lst.index(('Lyso-H2O', tmp_lyso))
                                            _tmp_i = lyso_info_dct['abs'][_tmp_lyso_idx]
                                            pre_lyso_fa_df.set_value(_ident_idx, 'type', 'Lyso-H2O')
                                            pre_lyso_fa_df.set_value(_ident_idx, 'abs', _tmp_i)
                                            _match_lyso_df = _match_lyso_df.append(pre_lyso_fa_df.ix[_ident_idx])
                                        if ('Lyso', tmp_lyso) in lyso_info_zip_lst:
                                            _tmp_lyso_idx = lyso_info_zip_lst.index(('Lyso', tmp_lyso))
                                            _tmp_i = lyso_info_dct['abs'][_tmp_lyso_idx]
                                            pre_lyso_fa_df.set_value(_ident_idx, 'type', 'Lyso')
                                            pre_lyso_fa_df.set_value(_ident_idx, 'abs', _tmp_i)
                                            _match_lyso_df = _match_lyso_df.append(pre_lyso_fa_df.ix[_ident_idx])

                                    print ('_match_fa_df', _match_fa_df.shape, '_match_lyso_df', _match_lyso_df.shape)

                                    _ident_df = ident_struct.check(_hmdb_id, _match_fa_df,
                                                                   _match_lyso_df, usr_std_fa_df)
                                    print ('_ident_df')
                                    print (_ident_df)

                                    if _ident_df.shape[0] > 0:
                                        # print fa_info

                                        _ident_table_df = _ident_df.loc[:, ['Abbr', 'Score']]
                                        _ident_table_df = _ident_table_df.sort_values(by=['Score', 'Abbr'],
                                                                                      ascending=[False, True])
                                        # print table_info_df.head(5)
                                        old_idx_lst = _ident_table_df.index.tolist()
                                        _ident_table_df.index = range(1, len(old_idx_lst) + 1)

                                        col_labels = _ident_table_df.columns.tolist()
                                        # col_labels = mz_info_df.head().tolist()
                                        row_labels = _ident_table_df.index.tolist()
                                        table_vals = map(list, _ident_table_df.values)

                                        try:
                                            # the rectangle is where I want to place the table
                                            fa_table = ms_pic.table(cellText=table_vals, rowLabels=row_labels,
                                                                    colWidths=[.2] * len(col_labels),
                                                                    colLabels=col_labels, loc='upper left')
                                            fa_table.set_fontsize(6)
                                            _auto_ident_chker += 1
                                        except:
                                            pass

                                    # if _auto_ident_chker > 0:
                                    # set title
                                    xic_title_str = 'XIC of m/z %.4f Abbr.: %s @ %.4f' % (_ms, _hmdb_id, _mzlib_id)
                                    ms_title_str = 'MS @ %.3f min [top 1000]' % _ms_rt_pr
                                    ms_zoom_title_str = 'MS zoomed'
                                    msms_title_str = ('MS/MS of m/z %.4f @ fuc %d - %.3f min [top 500]' %
                                                      (_pr_mz, _func_id, _pr_rt))
                                    msms_low_str = 'MS/MS zoomed below 350'
                                    msms_high_str = 'MS/MS zoomed above 350'

                                    xic_pic.set_title(xic_title_str, color='b', fontsize=10, y=0.98)
                                    ms_pic.set_title(ms_title_str, color='b', fontsize=10, y=0.98)
                                    ms_zoom_pic.set_title(ms_zoom_title_str, color='b', fontsize=10, y=0.98)
                                    msms_pic.set_title(msms_title_str, color='b', fontsize=10, y=0.98)
                                    msms_low_pic.set_title(msms_low_str, color='b', fontsize=10, y=0.98)
                                    msms_high_pic.set_title(msms_high_str, color='b', fontsize=10, y=0.98)

                                    # if _fa_checker > 0 and _hg_checker > 0:
                                    _n_hmdb_id = _hmdb_id
                                    try:
                                        _n_hmdb_id = _n_hmdb_id.replace('(', '[')
                                        _n_hmdb_id = _n_hmdb_id.replace(')', ']')
                                        _n_hmdb_id = _n_hmdb_id.replace(':', '-')
                                        _n_hmdb_id = _n_hmdb_id.replace('\\', '_')
                                        _n_hmdb_id = _n_hmdb_id.replace('/', '_')
                                    except:
                                        print ('Nothing to replace in ID %s' % _hmdb_id)
                                    print ('_n_hmdb_id', _n_hmdb_id)
                                    image_name_str = _path + r'/mz%.4f_%.3fmin_%s.png' % (_pr_mz, _pr_rt, _n_hmdb_id)

                                    plt.savefig(image_name_str, dpi=300)
                                    print (image_name_str, '===> Saved!')
                                    plt.close()
                                    _fa_checker = 0
                                    _hg_checker = 0
                                else:
                                    pass
                            else:
                                plt.close()
                                print ('Not identified==========>>>>>>>>>')

                                # except KeyError:
                                #     print ('error')
                                #     pass

                        else:
                            break
