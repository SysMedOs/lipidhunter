# -*- coding: utf-8 -*-
# Copyright 2014-16 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig          
# The software is currently  under development and is not ready to be released. 
# A suitable license will be chosen before the official release.               
# For more info please contact zhixu.ni@uni-leipzig.de

from matplotlib import pyplot as plt


class Spectra_Ploter(object):

    def plot_dual(self, xic_df, msms_spectra_dct):

        xic_mz_lst = xic_df['rt'].tolist()
        xic_i_lst = xic_df['i'].tolist()

        # print xic_mz_lst
        # print xic_i_lst

        for _msms in msms_spectra_dct.keys():
            _pr_mz = _msms[0]
            _pr_rt = _msms[1]
            _msms_df = msms_spectra_dct[_msms]

            fig, (xic_pic, msms_pic, msms_zoom_pic) = plt.subplots(nrows=3, figsize=(12, 18))

            # XIC spectrum start
            xic_pic.plot(xic_mz_lst, xic_i_lst, alpha=0.3)
            xic_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            xic_pic.stem([_pr_rt], [max(xic_i_lst)], '-.')
            xic_pic.set_xlabel("RT (min)")
            xic_pic.set_ylabel("Intensity")

            # msms spectrum start
            msms_pic.stem(_msms_df['mz'].tolist(), _msms_df['i'].tolist(), 'blue', lw=4, markerfmt=" ")
            msms_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            msms_pic.set_xlabel("m/z")
            msms_pic.set_ylabel("Intensity")

            # msms spectrum start
            _msms_zoom_df = _msms_df.query('mz < 255')
            if len(_msms_zoom_df['mz'].tolist()) > 0:
                msms_zoom_pic.stem(_msms_zoom_df['mz'].tolist(), _msms_zoom_df['i'].tolist(), 'black', lw=4, markerfmt=" ")
                msms_zoom_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
                msms_zoom_pic.set_xlabel("m/z")
                msms_zoom_pic.set_ylabel("Intensity")
            else:
                pass

            xic_title_str = 'XIC of m/z %.4f with MS/MS @ %.2f min' % (_pr_mz, _pr_rt)
            image_name_str = 'mz%.4f_%.2fmin.png' % (_pr_mz, _pr_rt)

            xic_pic.set_title(xic_title_str, color='b')

            plt.savefig(image_name_str, dpi=300)
            print image_name_str, '===> Saved!'
            plt.close()


