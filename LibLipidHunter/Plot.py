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


class Spectra_Ploter(object):

    def plot_all(self, mz2get_lst, sdf_df, xic_dct, ms_spectra_dct, msms_spectra_dct, path=None):
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

        for mz2get in mz2get_lst:

            xic_df = xic_dct[mz2get]
            _ms_spectra_dct = ms_spectra_dct[mz2get]
            _msms_spectra_dct = msms_spectra_dct[mz2get]

            _query_str = 'pr_mz ==' + str(mz2get)
            _sdf_df = sdf_df.query(_query_str)
            # print _sdf_df

            for index, row in _sdf_df.iterrows():
                _hmdb_id = row['hmdb_id']
                _sn1_mz = row['sn1_mz']
                _sn2_mz = row['sn2_mz']
                _pr_smi = row['pr_smi']

                try:

                    xic_mz_lst = xic_df['rt'].tolist()
                    xic_i_lst = xic_df['i'].tolist()

                    ms_rt_lst = _ms_spectra_dct.keys()
                    ms_rt_lst.sort()
                    # print ms_rt_lst
                    for _msms in _msms_spectra_dct.keys():
                        _pr_mz = _msms[0]
                        _pr_rt = _msms[1]
                        _msms_df = _msms_spectra_dct[_msms]

                        _ms_rt_lst = []
                        for _ms_rt in ms_rt_lst:
                            if _ms_rt < _pr_rt:
                                _ms_rt_lst.append(_ms_rt)
                        _ms_rt_lst.sort()
                        _ms_rt_pr = _ms_rt_lst[-1]

                        # print _ms_rt_lst
                        print 'MS @ %.2f for DDA @ %.2f' % (_ms_rt_pr, _pr_rt)

                        # skip the image if MS and MS/MS time do not match
                        if _pr_rt - _ms_rt_pr > 0.2:
                            print 'MS and MS/MS time NOT match.'
                            pass

                        if len(_msms_df['mz'].tolist()) < 1 or len(_msms_df['i'].tolist())<1:
                            print 'No peaks in MS/MS'

                        else:
                            # Generate A4 image in landscape
                            fig, pic_array = plt.subplots(nrows=3, ncols=2, figsize=(11.692, 8.267), sharex=False, sharey=False)
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
                            _ms_pkl_lst = _ms_spectra_dct[_ms_rt_pr][1]
                            _ms_pkl_df = pd.DataFrame(data=_ms_pkl_lst, columns=['mz', 'i'])
                            ms_pic.stem(_ms_pkl_df['mz'].tolist(), _ms_pkl_df['i'].tolist(),
                                        'blue', lw=4, markerfmt=" ")
                            _ms_pr_lst = _ms_spectra_dct[_ms_rt_pr][0]
                            _dash_i = [max(_ms_pkl_df['i'].tolist())] * len(_ms_pr_lst[0])
                            ms_pic.stem(_ms_pr_lst[0], _dash_i, ':', 'yellow', markerfmt=" ")
                            ms_pic.stem(_ms_pr_lst[0], _ms_pr_lst[1], 'red', lw=4, markerfmt=" ")

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

                            _mz_zoom_query_str = ' %.2f < mz < %.2f' % (mz2get - 1.1, mz2get + 1.1)
                            _ms_zoom_df = _ms_pkl_df.query(_mz_zoom_query_str)
                            ms_zoom_pic.stem(_ms_zoom_df['mz'].tolist(), _ms_zoom_df['i'].tolist(),
                                             'black', lw=4, markerfmt=" ")
                            _ms_pr_lst = _ms_spectra_dct[_ms_rt_pr][0]
                            ms_zoom_pic.stem(_ms_pr_lst[0], _ms_pr_lst[1], 'red', lw=4, markerfmt=" ")

                            ms_zoom_pic.set_xlim([mz2get - 1.1, mz2get + 1.1])
                            ms_zoom_pic.set_ylim([0, max(_ms_zoom_df['i'].tolist()) * 1.3])
                            ms_zoom_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0), fontsize=10)
                            ms_zoom_pic.ticklabel_format(axis='x', useOffset=False, fontsize=10)
                            ms_zoom_pic.set_xlabel("m/z", fontsize=10, labelpad=-1)
                            ms_zoom_pic.set_ylabel("Intensity", fontsize=10)

                            # add annotation
                            _ms_zoom_peak_list = zip(_ms_zoom_df['mz'].tolist(), _ms_zoom_df['i'].tolist())
                            for _ms_zoom_peak in _ms_zoom_peak_list:
                                _ms_zoom_peak_str = '%.4f' % _ms_zoom_peak[0]
                                _ms_zoom_peak_y = (max(_ms_zoom_df['i'].tolist()) * 1.3) * 0.175 + _ms_zoom_peak[1]
                                ms_zoom_pic.text(_ms_zoom_peak[0], _ms_zoom_peak_y, _ms_zoom_peak_str,
                                                 rotation=90, fontsize=6)

                            # XIC spectrum start
                            xic_pic.plot(xic_mz_lst, xic_i_lst, alpha=0.3)
                            xic_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
                            xic_pic.stem([_ms_rt_pr], [max(xic_i_lst)], ':', markerfmt=" ")
                            xic_pic.stem([_pr_rt], [max(xic_i_lst)], ':', markerfmt=" ")
                            xic_pic.text(_ms_rt_pr - 0.3, max(xic_i_lst)*0.98, 'MS', fontsize=10)
                            xic_pic.text(_pr_rt, max(xic_i_lst)*0.98, 'MS/MS', fontsize=10)
                            xic_pic.set_xlabel("RT (min)", fontsize=10, labelpad=-1)
                            xic_pic.set_ylabel("Intensity", fontsize=10)

                            # msms spectrum start
                            msms_pic.stem(_msms_df['mz'].tolist(), _msms_df['i'].tolist(), 'blue', lw=2, markerfmt=" ")
                            msms_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
                            msms_pic.set_xlabel("m/z", fontsize=10, labelpad=-1)
                            msms_pic.set_ylabel("Intensity", fontsize=10)
                            msms_pic.set_xlim([min(_msms_df['mz'].tolist()) - 1, _pr_mz + 20])
                            msms_pic.set_ylim([0, max(_msms_df['i'].tolist()) * 1.3])

                            # add annotations
                            _msms_peak_list = zip(_msms_df['mz'].tolist(), _msms_df['i'].tolist())
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

                            for _mz in [_sn1_mz, _sn2_mz]:
                                _msms_top_df = _msms_df.sort(columns='i', ascending=False)
                                _msms_top_df = _msms_top_df.head(10)
                                _query_str = '%f - %f <= mz <= %f + %f' % (_mz, 0.25, _mz, 0.25)
                                _sn_msms_df = _msms_top_df.query(_query_str)
                                _sn_msms_df = _sn_msms_df.sort(columns='i', ascending=False)
                                # print _sn_msms_df.head(2)
                                try:
                                    msms_pic.stem([_sn_msms_df['mz'].tolist()[0]], [_sn_msms_df['i'].tolist()[0]],
                                                  'red', lw=4, markerfmt=" ")
                                    _fa_checker += 1
                                except:
                                    print('no sn fit')

                            # msms spectrum zoomed below 255 start
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
                                for _mz in [168.0358, 224.0693, 242.0798]:
                                    _query_str = '%f - %f <= mz <= %f + %f' % (_mz, 0.25, _mz, 0.25)
                                    _sn_msms_df = _msms_df.query(_query_str)
                                    _sn_msms_df = _sn_msms_df.sort(columns='i', ascending=False)
                                    # print _sn_msms_df.head(2)
                                    try:
                                        msms_low_pic.stem([_sn_msms_df['mz'].tolist()[0]],
                                                          [_sn_msms_df['i'].tolist()[0]], 'red', lw=4, markerfmt=" ")
                                        _hg_checker += 1
                                    except:
                                        print('no sn fit')

                                # add annotations
                                _msms_low_peak_list = zip(_msms_low_df['mz'].tolist(), _msms_low_df['i'].tolist())
                                for _msms_low_peak in _msms_low_peak_list:
                                    _msms_low_peak_str = '%.4f' % _msms_low_peak[0]
                                    _msms_low_peak_y = (max(_msms_low_df['i'].tolist()) * 1.3) * 0.175 + _msms_low_peak[1]
                                    msms_low_pic.text(_msms_low_peak[0], _msms_low_peak_y, _msms_low_peak_str,
                                                      rotation=90, fontsize=6)
                            else:
                                pass

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
                                _msms_high_peak_list = zip(_msms_high_df['mz'].tolist(), _msms_high_df['i'].tolist())
                                for _msms_high_peak in _msms_high_peak_list:
                                    _msms_high_peak_str = '%.4f' % _msms_high_peak[0]
                                    _msms_high_peak_y = ((max(_msms_high_df['i'].tolist()) * 1.3) * 0.175 +
                                                         _msms_high_peak[1])
                                    msms_high_pic.text(_msms_high_peak[0], _msms_high_peak_y, _msms_high_peak_str,
                                                       rotation=90, fontsize=6)
                            else:
                                pass

                            # set title
                            xic_title_str = 'XIC of m/z %.3f HMDB_ID: %s' % (mz2get, _hmdb_id)
                            ms_title_str = 'MS @ %.2f min [top 200]' % _ms_rt_pr
                            ms_zoom_title_str = 'MS zoomed'
                            msms_title_str = 'MS/MS of m/z %.2f @ %.2f min [top 200]' % (_pr_mz, _pr_rt)
                            msms_low_str = 'MS/MS zoomed below 350'
                            msms_high_str = 'MS/MS zoomed above 350'

                            xic_pic.set_title(xic_title_str, color='b', fontsize=10, y=0.98)
                            ms_pic.set_title(ms_title_str, color='b', fontsize=10, y=0.98)
                            ms_zoom_pic.set_title(ms_zoom_title_str, color='b', fontsize=10, y=0.98)
                            msms_pic.set_title(msms_title_str, color='b', fontsize=10, y=0.98)
                            msms_low_pic.set_title(msms_low_str, color='b', fontsize=10, y=0.98)
                            msms_high_pic.set_title(msms_high_str, color='b', fontsize=10, y=0.98)

                            if _fa_checker > 0 and _hg_checker > 0:
                                image_name_str = _path + r'/mz%.4f_%.2fmin-%s.png' % (_pr_mz, _pr_rt, _hmdb_id)

                                plt.savefig(image_name_str, dpi=300)
                                print (image_name_str, '===> Saved!')
                                plt.close()
                                _fa_checker = 0
                                _hg_checker = 0
                            else:
                                pass
                except KeyError:
                    pass
