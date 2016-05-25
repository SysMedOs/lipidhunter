# -*- coding: utf-8 -*-
# Copyright 2014-16 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release.
# For more info please contact zhixu.ni@uni-leipzig.de
import pandas as pd


class FAindicator(object):

    def __init__(self, fa_csv):

        self.fa_df = pd.read_csv(fa_csv, index_col=0)
        self.fa_lst = self.fa_df.index.tolist()

    def indicate(self, spec_df, charge=None, top=None, ppm=None):

        if charge in ['PL[+]_mz', 'PL[Na+]_mz', 'PL[-]_mz', 'PL[FA-]_mz']:
            pass
        else:
            charge = 'PL[-]_mz'

        if charge == 'PL[-]_mz':
            fa_mz_typ_lst = ['[M-H]-', '[M-W-H]-']

        if charge == 'PL[FA-]_mz':
            fa_mz_typ_lst = ['[M-H]-', '[M+FA-H]-', '[M-W-H]-', 'M+FA-W-H]-']

        else:
            fa_mz_typ_lst = ['[M-H]-', '[M-W-H]-']

        if top > 0:
            pass
        else:
            top = 50

        if ppm > 0:
            ppm_neg = 0 - ppm
        else:
            ppm = 50.0
            ppm_neg = -50.0

        match_lst = []
        match_mz_lst = []
        match_i_lst = []
        match_name_lst = []
        match_delta_lst = []
        match_ppm_lst = []
        spec_lmw_df = spec_df
        # spec_lmw_df = spec_lmw_df[spec_lmw_df['mz'] < 430.0]
        # spec_lmw_df = spec_lmw_df[spec_lmw_df['mz'] > 220.0]
        spec_lmw_df = spec_lmw_df[(spec_lmw_df['mz'] > 220.0) & (spec_lmw_df['mz'] < 430.0)]
        spec_lmw_df = spec_lmw_df.sort_values(by='i', ascending=False)
        spec_lmw_df = spec_lmw_df.tail(top)
        for idx in spec_lmw_df.index.tolist():
            tmp_mz = spec_lmw_df.loc[idx, 'mz']

            for fa in self.fa_lst:
                for mz_typ in fa_mz_typ_lst:
                    # print 'mz_target_mz =', mz_target_mz
                    mz_target_mz = self.fa_df.loc[fa, mz_typ]
                    if mz_target_mz - 0.5 <= tmp_mz <= mz_target_mz + 0.5:
                        tmp_ppm = ((tmp_mz-mz_target_mz)/mz_target_mz) * 1000000
                        if ppm_neg <= tmp_ppm <= ppm:
                            match_lst.append(idx)
                            tmp_i = spec_lmw_df.loc[idx, 'i']
                            match_i_lst.append('%.2e' % tmp_i)
                            match_mz_lst.append(tmp_mz)
                            fa_label = fa + mz_typ
                            match_name_lst.append(fa_label)
                            tmp_delta = tmp_mz-mz_target_mz
                            tmp_delta_txt = '%.2f Da' % tmp_delta
                            match_delta_lst.append(tmp_delta_txt)
                            tmp_ppm_txt = '%.2f ppm' % tmp_ppm
                            match_ppm_lst.append(tmp_ppm_txt)
                        else:
                            pass
                    else:
                        pass
        print('match_mz_lst', match_mz_lst)

        if len(match_lst) > 0:
            m_df = spec_df.ix[match_lst]

        else:
            m_df = pd.DataFrame()

        m_info_dct = {'name': match_name_lst, 'mz': match_mz_lst, 'i': match_i_lst,
                      'D': match_delta_lst, 'ppm': match_ppm_lst}

        # print m_df

        return m_df, m_info_dct

class Lyso_indicator(object):

    def __init__(self, fa_csv):

        self.fa_df = pd.read_csv(fa_csv, index_col=0)
        self.fa_lst = self.fa_df.index.tolist()

    def indicate(self, spec_df, prmz, top=None, ppm=None):

        lyso_typ_lst = ['mass', 'NL-H2O']

        if top > 0:
            pass
        else:
            top = 50

        if ppm > 0:
            ppm_neg = 0 - ppm
        else:
            ppm = 50.0
            ppm_neg = -50.0

        match_lst = []
        match_mz_lst = []
        match_i_lst = []
        match_name_lst = []
        match_delta_lst = []
        match_ppm_lst = []
        spec_lmw_df = spec_df
        spec_lmw_df = spec_lmw_df[spec_lmw_df['mz'] > 400.0]
        spec_lmw_df = spec_lmw_df.sort_values(by='i', ascending=False)
        spec_lmw_df = spec_lmw_df.tail(top)
        for idx in spec_lmw_df.index.tolist():
            tmp_mz = spec_lmw_df.loc[idx, 'mz']

            for fa in self.fa_lst:
                for mz_typ in lyso_typ_lst:
                    # print 'mz_target_mz =', mz_target_mz
                    mz_target_mz = prmz - self.fa_df.loc[fa, mz_typ]
                    if mz_target_mz - 0.5 <= tmp_mz <= mz_target_mz + 0.5:
                        tmp_ppm = ((tmp_mz-mz_target_mz)/mz_target_mz) * 1000000
                        if ppm_neg <= tmp_ppm <= ppm:
                            match_lst.append(idx)
                            tmp_i = spec_lmw_df.loc[idx, 'i']
                            match_i_lst.append('%.2e' % tmp_i)
                            match_mz_lst.append(tmp_mz)
                            if mz_typ == 'mass':
                                fa_label = '[Lyso%s-H2O-H]-' % fa
                                match_name_lst.append(fa_label)
                                tmp_delta = tmp_mz - mz_target_mz
                                tmp_delta_txt = '%.2f Da' % tmp_delta
                                match_delta_lst.append(tmp_delta_txt)
                                tmp_ppm_txt = '%.2f ppm' % tmp_ppm
                                match_ppm_lst.append(tmp_ppm_txt)
                            if mz_typ == 'NL-H2O':
                                fa_label = '[Lyso%s-H]-' % fa
                                match_name_lst.append(fa_label)
                                tmp_delta = tmp_mz - mz_target_mz
                                tmp_delta_txt = '%.2f Da' % tmp_delta
                                match_delta_lst.append(tmp_delta_txt)
                                tmp_ppm_txt = '%.2f ppm' % tmp_ppm
                                match_ppm_lst.append(tmp_ppm_txt)
                            else:
                                pass

                        else:
                            pass
                    else:
                        pass
        print('match_mz_lst', match_mz_lst)

        if len(match_lst) > 0:
            m_df = spec_df.ix[match_lst]

        else:
            m_df = pd.DataFrame()

        m_info_dct = {'name': match_name_lst, 'mz': match_mz_lst, 'i': match_i_lst,
                      'D': match_delta_lst, 'ppm': match_ppm_lst}

        # print m_df

        return m_df, m_info_dct
