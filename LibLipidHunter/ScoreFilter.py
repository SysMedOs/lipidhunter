# -*- coding: utf-8 -*-
# Copyright 2015-2017 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

import pandas as pd


def check_peaks(score_df, fa_ident_df, lyso_ident_df, lyso_w_ident_df, score_filter=20):

    if fa_ident_df.shape[0] or lyso_ident_df.shape[0]:
        # score_df = score_df[['Lipid_species', 'Score']]
        # score_df = score_df.rename({'Lipid_species': 'Proposed structures'})
        score_df = score_df.query('Score >= %.2f' % score_filter)
        score_df = score_df.sort_values(by='Score', ascending=False)
        score_df = score_df.reset_index(drop=True)
        score_df.index += 1
        print(score_df)

        # format fa info DataFrame
        if fa_ident_df.shape[0] > 0:
            fa_ident_df = fa_ident_df[['Proposed_structures', 'mz', 'i', 'ppm']].reset_index(drop=True)
            # fa_ident_df = fa_ident_df.rename({'Lipid_species': 'Identified species'})
            fa_ident_df = fa_ident_df.round({'mz': 4, 'ppm': 2})
            _fa_i_lst = []
            for _idx, _fa_se in fa_ident_df.iterrows():
                _fa_i_lst.append('%.2e' % float(_fa_se['i']))
            fa_ident_df.loc[:, 'i'] = _fa_i_lst
            fa_ident_df.index += 1
            print(fa_ident_df)

        # merge Lyso and Lyso - H2O
        lyso_ident_df = lyso_ident_df.append(lyso_w_ident_df)
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
            print(lyso_ident_df)
        else:
            lyso_ident_df = pd.DataFrame()

        usr_ident_info_dct = {'SCORE_INFO': score_df, 'FA_INFO': fa_ident_df, 'LYSO_INFO': lyso_ident_df}
        print usr_ident_info_dct
    else:
        usr_ident_info_dct= {'SCORE_INFO':pd.DataFrame(), 'FA_INFO': '', 'LYSO_INFO': ''}
    return usr_ident_info_dct
