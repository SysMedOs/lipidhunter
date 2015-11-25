# -*- coding: utf-8 -*-
# Copyright 2014-16 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig          
# The software is currently  under development and is not ready to be released. 
# A suitable license will be chosen before the official release.               
# For more info please contact zhixu.ni@uni-leipzig.de

from matplotlib import pyplot as plt
import pandas as pd


class Spectra_Ploter(object):

    def plot_all(self, mz2get, xic_df, ms_spectra_dct, msms_spectra_dct):

        xic_mz_lst = xic_df['rt'].tolist()
        xic_i_lst = xic_df['i'].tolist()

        # print xic_mz_lst
        # print xic_i_lst

        ms_rt_lst = ms_spectra_dct.keys()
        ms_rt_lst.sort()
        print ms_rt_lst
        for _msms in msms_spectra_dct.keys():
            _pr_mz = _msms[0]
            _pr_rt = _msms[1]
            _msms_df = msms_spectra_dct[_msms]

            fig, pic_array = plt.subplots(nrows=3, ncols=2, figsize=(16, 12), sharex=False, sharey=False)
            xic_pic = pic_array[0, 0]
            msms_pic = pic_array[0, 1]
            ms_pic = pic_array[1, 0]
            msms_low_pic = pic_array[1, 1]
            ms_zoom_pic = pic_array[2, 0]
            msms_high_pic = pic_array[2, 1]

            # ms spectrum start
            _ms_rt_lst = []
            for _ms_rt in ms_rt_lst:
                if _ms_rt < _pr_rt:
                    _ms_rt_lst.append(_ms_rt)
            _ms_rt_lst.sort()
            _ms_rt_pr = _ms_rt_lst[-1]
            print _ms_rt_lst
            print 'MS scan for this DDA @', _ms_rt_pr

            _ms_pkl_lst = ms_spectra_dct[_ms_rt_pr][1]
            _ms_pkl_df = pd.DataFrame(data=_ms_pkl_lst, columns=['mz', 'i'])
            ms_pic.stem(_ms_pkl_df['mz'].tolist(), _ms_pkl_df['i'].tolist(), 'blue', lw=4, markerfmt=" ")
            _ms_pr_lst = ms_spectra_dct[_ms_rt_pr][0]
            _dash_i = [max(_ms_pkl_df['i'].tolist())] * len(_ms_pr_lst[0])
            ms_pic.stem(_ms_pr_lst[0], _dash_i, ':', 'yellow', markerfmt=" ")
            ms_pic.stem(_ms_pr_lst[0], _ms_pr_lst[1], 'red', lw=4, markerfmt=" ")

            ms_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            ms_pic.set_xlabel("m/z")
            ms_pic.set_ylabel("Intensity")

            _mz_zoom_query_str = ' %.2f < mz < %.2f' % (mz2get - 0.75, mz2get + 2)
            _ms_zoom_df = _ms_pkl_df.query(_mz_zoom_query_str)
            # ms_zoom_pic.stem(_ms_zoom_df['mz'].tolist(), _ms_zoom_df['i'].tolist(), 'black', lw=4, markerfmt=" ")
            _ms_pr_lst = ms_spectra_dct[_ms_rt_pr][0]
            ms_zoom_pic.stem(_ms_pr_lst[0], _ms_pr_lst[1], 'red', lw=4, markerfmt=" ")

            ms_zoom_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            ms_zoom_pic.ticklabel_format(axis='x', useOffset=False)
            ms_zoom_pic.set_xlabel("m/z")
            ms_zoom_pic.set_ylabel("Intensity")
            ms_zoom_pic.set_xlim([mz2get - 0.75, mz2get + 2])

            # XIC spectrum start
            xic_pic.plot(xic_mz_lst, xic_i_lst, alpha=0.3)
            xic_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            xic_pic.stem([_ms_rt_pr], [max(xic_i_lst)], ':', markerfmt=" ")
            xic_pic.stem([_pr_rt], [max(xic_i_lst)], '-.', markerfmt=" ")
            xic_pic.text(_ms_rt_pr - 0.3, max(xic_i_lst), 'MS')
            xic_pic.text(_pr_rt, max(xic_i_lst), 'MS/MS')
            xic_pic.set_xlabel("RT (min)")
            xic_pic.set_ylabel("Intensity")

            # msms spectrum start
            msms_pic.stem(_msms_df['mz'].tolist(), _msms_df['i'].tolist(), 'blue', lw=4, markerfmt=" ")
            msms_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            msms_pic.set_xlabel("m/z")
            msms_pic.set_ylabel("Intensity")

            # msms spectrum zoomed below 255 start
            _msms_low_df = _msms_df.query('mz < 255')
            if len(_msms_low_df['mz'].tolist()) > 0:
                msms_low_pic.stem(_msms_low_df['mz'].tolist(),
                                  _msms_low_df['i'].tolist(),
                                  'black', lw=4, markerfmt=" ")
                msms_low_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
                msms_low_pic.set_xlabel("m/z")
                msms_low_pic.set_ylabel("Intensity")
            else:
                pass

            # msms spectrum zoomed above 350 start
            _msms_high_df = _msms_df.query('mz > 350')
            if len(_msms_high_df['mz'].tolist()) > 0:
                msms_high_pic.stem(_msms_high_df['mz'].tolist(),
                                   _msms_high_df['i'].tolist(),
                                   'black', lw=4, markerfmt=" ")
                msms_high_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
                msms_high_pic.set_xlabel("m/z")
                msms_high_pic.set_ylabel("Intensity")
            else:
                pass

            # set title

            xic_title_str = 'XIC of m/z %.3f' % mz2get
            ms_title_str = 'MS @ %.2f min [top 200]' % _ms_rt_pr
            msms_title_str = 'MS/MS of m/z %.2f @ %.2f min [top 200]' % (_pr_mz, _pr_rt)
            msms_low_str = 'MS/MS zoomed below 255'
            # msms_high_str = 'MS/MS zoomed above 350'

            xic_pic.set_title(xic_title_str, color='b')
            ms_pic.set_title(ms_title_str, color='b')
            msms_pic.set_title(msms_title_str, color='b')
            msms_low_pic.set_title(msms_low_str, color='b')
            # msms_high_pic.set_title(msms_high_str, color='b')

            image_name_str = 'mz%.4f_%.2fmin.png' % (_pr_mz, _pr_rt)

            plt.savefig(image_name_str, dpi=300)
            print image_name_str, '===> Saved!'
            plt.close()


