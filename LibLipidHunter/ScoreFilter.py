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
#
# For more info please contact:
#     SysMedOs_team: oxlpp@bbz.uni-leipzig.de
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de
#     Developer Georgia Angelidou georgia.angelidou@uni-leipzig.de
#

import pandas as pd


def check_peaks(score_df, fa_ident_df, lyso_ident_df, lyso_w_ident_df, mg_w_ident_df, score_filter=20):

    #mg_w_ident_df = pd.DataFrame()
    if fa_ident_df.shape[0] or lyso_ident_df.shape[0] or lyso_w_ident_df.shape[0] or mg_w_ident_df.shape[0]:
        # score_df = score_df[['Lipid_species', 'Score']]
        # score_df = score_df.rename({'Lipid_species': 'Proposed structures'})
        score_df = score_df.query('Score >= %.2f' % score_filter)
        if score_df.shape[0]:
            score_df = score_df.sort_values(by='Score', ascending=False)
            score_df = score_df.reset_index(drop=True)
            score_df.index += 1


            # format fa info DataFrame
            fa_ident_df = fa_ident_df[['Proposed_structures', 'mz', 'i', 'ppm']].reset_index(drop=True)
            # fa_ident_df = fa_ident_df.rename({'Lipid_species': 'Identified species'})
            fa_ident_df = fa_ident_df.round({'mz': 4, 'ppm': 2})
            _fa_i_lst = []
            for _idx, _fa_se in fa_ident_df.iterrows():
                _fa_i_lst.append('%.2e' % float(_fa_se['i']))
            fa_ident_df.loc[:, 'i'] = _fa_i_lst
            fa_ident_df.index += 1


            # merge Lyso and Lyso - H2O
            lyso_ident_df = lyso_ident_df.append(lyso_w_ident_df)
            lyso_ident_df = lyso_ident_df.append(mg_w_ident_df)
            if lyso_ident_df.shape[0] > 0:
                lyso_ident_df = lyso_ident_df.sort_values(by='i', ascending=False)
                lyso_ident_df = lyso_ident_df[['Proposed_structures', 'mz', 'i', 'ppm']].reset_index(drop=True)
                # lyso_ident_df = lyso_ident_df.rename({'Lipid_species': 'Identified species'})
                lyso_ident_df = lyso_ident_df.round({'mz': 4, 'ppm': 2})
                _lyso_i_lst = []
                for _idx, _lyso_se in lyso_ident_df.iterrows():
                    _lyso_i_lst.append('%.2e' % float(_lyso_se['i']))
                lyso_ident_df.loc[:, 'i'] = _lyso_i_lst
                lyso_ident_df.index += 1

                lyso_ident_df = lyso_ident_df.drop_duplicates()
            else:
                lyso_ident_df = pd.DataFrame()

            usr_ident_info_dct = {'SCORE_INFO': score_df, 'FA_INFO': fa_ident_df, 'LYSO_INFO': lyso_ident_df}
        else:

            usr_ident_info_dct = {'SCORE_INFO': pd.DataFrame(), 'FA_INFO': pd.DataFrame(), 'LYSO_INFO': pd.DataFrame()}
    else:
        usr_ident_info_dct = {'SCORE_INFO': pd.DataFrame(), 'FA_INFO': pd.DataFrame(), 'LYSO_INFO': pd.DataFrame()}

    return usr_ident_info_dct
