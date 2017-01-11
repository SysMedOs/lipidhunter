# -*- coding: utf-8 -*-
# Copyright 2015-2016 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

import re

import pandas as pd


class ScoreGenerator:
    def __init__(self, fa_def_df, weight_df, key_frag_df, lipid_type):
        self.fa_def_df = fa_def_df
        self.weight_df = weight_df
        self.target_frag_df = key_frag_df.query('CLASS == "%s" and TYPE == "FRAG"' % lipid_type)
        self.target_nl_df = key_frag_df.query('CLASS == "%s" and TYPE == "NL"' % lipid_type)
        self.other_frag_df = key_frag_df.query('CLASS != "%s" and TYPE == "FRAG"' % lipid_type)
        self.other_nl_df = key_frag_df.query('CLASS != "%s" and TYPE == "NL"' % lipid_type)
        self.lipid_type = lipid_type

    @staticmethod
    def get_pr_mz(charge_type, mz_lib):

        pr_mz = 0.0

        if charge_type in ['[M-H]-', '[M+HCOO]-', '[M+FA-H]-', '[M+CH3COO]-', '[M+OAc-H]-', '[M+AcOH-H]-']:
            charge_mode = 'NEG'
            if charge_type == '[M-H]-':
                pr_mz = mz_lib
            elif charge_type in ['[M+HCOO]-', '[M+FA-H]-']:
                pr_mz = mz_lib - 46.005480  # - HCOOH
            elif charge_type == '[M+CH3COO]-':
                pr_mz = mz_lib - 60.021130  # - CH3COOH

        elif charge_type in ['[M+H]+', '[M+Na]+', '[M+NH4]+', '[M+K]+']:
            charge_mode = 'POS'
            if charge_type == '[M+H]+':
                pr_mz = mz_lib
            elif charge_type == '[M+Na]+':
                pr_mz = mz_lib - 22.989770 + 1.007825  # - Na + H
            elif charge_type == '[M+NH4]+':
                pr_mz = mz_lib - 17.026549  # - NH3
            elif charge_type == '[M+K]+':
                pr_mz = mz_lib - 38.963708 + 1.007825  # - K + H
        else:
            charge_mode = 'NEG'
            pr_mz = mz_lib

        return pr_mz, charge_mode

    @staticmethod
    def decode_abbr(abbr):

        pl_checker = re.compile(r'(P[ACEGSI])([(])(.*)([)])')
        pip_checker = re.compile(r'(PIP)([(])(.*)([)])')
        tg_checker = re.compile(r'(TG)([(])(.*)([)])')
        fa_checker = re.compile(r'(\d{1,2})([:])(\d)')
        fa_o_checker = re.compile(r'(O-)(\d{1,2})([:])(\d)')
        fa_p_checker = re.compile(r'(P-)(\d{1,2})([:])(\d)')

        # Check PL Type
        _pl_typ = ''
        bulk_fa_typ = ''
        bulk_fa_linker = ''
        bulk_fa_c = 0
        bulk_fa_db = 0
        lyso_fa_linker_dct = {'sn1': '', 'sn2': ''}

        if pl_checker.match(abbr):
            print ('PL')
            pl_re_chk = pl_checker.match(abbr)
            pl_typ_lst = pl_re_chk.groups()
            _pl_typ = pl_typ_lst[0]
            bulk_fa_typ = pl_typ_lst[2]
        if pip_checker.match(abbr):
            print ('PIP')
            pip_re_chk = pip_checker.match(abbr)
            pip_typ_lst = pip_re_chk.groups()
            _pl_typ = pip_typ_lst[0]
            bulk_fa_typ = pip_typ_lst[2]
        if tg_checker.match(abbr):
            print ('TG')
            pip_re_chk = pip_checker.match(abbr)
            pip_typ_lst = pip_re_chk.groups()
            _pl_typ = pip_typ_lst[0]
            bulk_fa_typ = pip_typ_lst[2]
        if fa_checker.match(abbr):
            print ('FA')
            _pl_typ = 'FA'
            bulk_fa_typ = abbr
        if fa_o_checker.match(abbr):
            print ('FA')
            _pl_typ = 'FA'
            bulk_fa_typ = abbr
        if fa_p_checker.match(abbr):
            print ('FA')
            _pl_typ = 'FA'
            bulk_fa_typ = abbr

        print(bulk_fa_typ)

        if fa_checker.match(bulk_fa_typ):
            bulk_fa_linker = 'A-A-'
            lyso_fa_linker_dct = {'A': ''}
            fa_chk = fa_checker.match(bulk_fa_typ)
            bulk_fa_lst = fa_chk.groups()
            bulk_fa_c = bulk_fa_lst[0]
            bulk_fa_db = bulk_fa_lst[2]
        elif fa_o_checker.match(bulk_fa_typ):
            bulk_fa_linker = 'O-A-'
            lyso_fa_linker_dct = {'O': '', 'A': 'O-'}  # link of the other sn after NL of this sn
            fa_chk = fa_o_checker.match(bulk_fa_typ)
            bulk_fa_lst = fa_chk.groups()
            bulk_fa_c = bulk_fa_lst[1]
            bulk_fa_db = bulk_fa_lst[3]
        elif fa_p_checker.match(bulk_fa_typ):
            bulk_fa_linker = 'P-A-'
            lyso_fa_linker_dct = {'P': '', 'A': 'P-'}  # link of the other sn after NL of this sn
            fa_chk = fa_p_checker.match(bulk_fa_typ)
            bulk_fa_lst = fa_chk.groups()
            bulk_fa_c = bulk_fa_lst[1]
            bulk_fa_db = bulk_fa_lst[3]

        bulk_fa_c = int(bulk_fa_c)
        bulk_fa_db = int(bulk_fa_db)

        lipid_info_dct = {'TYPE': _pl_typ, 'LINK': bulk_fa_linker, 'C': bulk_fa_c, 'DB': bulk_fa_db,
                          'LYSO_LINK': lyso_fa_linker_dct}

        return lipid_info_dct

    def get_fa_search(self, abbr, charge_type, mz_lib, ms2_df, ms2_precision=500e-6, ms2_threshold=100):

        fa_ident_df = pd.DataFrame()
        lyso_ident_df = pd.DataFrame()
        lyso_w_ident_df = pd.DataFrame()

        lipid_info_dct = self.decode_abbr(abbr)
        pl_typ = lipid_info_dct['TYPE']
        bulk_fa_c = lipid_info_dct['C']
        bulk_fa_db = lipid_info_dct['DB']
        bulk_fa_linker = lipid_info_dct['LINK']
        lyso_fa_linker_dct = lipid_info_dct['LYSO_LINK']

        calc_pr_mz, charge_mode = self.get_pr_mz(charge_type, mz_lib)

        if abbr[:2] in ['PA', 'PC', 'PE', 'PG', 'PI', 'PS', 'SM']:
            lipid_type = 'PL'
        elif abbr[:2] in ['TA', 'TG', 'DA', 'DG', 'MA', 'MG']:
            lipid_type = 'GL'
        else:
            lipid_type = 'PL'

        if lipid_type == 'PL' and charge_mode == 'NEG':
            fa_chk_df = self.fa_def_df[['FA', 'Link', 'C', 'DB', 'mass', '[M-H]-', 'NL-H2O']]
            fa_chk_df = fa_chk_df.rename(columns={'[M-H]-': 'sn', 'mass': 'NL'})
            fa_chk_df['M-sn'] = calc_pr_mz - fa_chk_df['NL-H2O']
            fa_chk_df['M-(sn-H2O)'] = calc_pr_mz - fa_chk_df['NL']
            fa_chk_df['Lipid_species'] = ''

            for _i, _fa_se in fa_chk_df.iterrows():

                _fa_abbr = _fa_se['FA']
                _fa_link = _fa_se['Link']
                _fa_c = _fa_se['C']
                _fa_db = _fa_se['DB']

                for _frag_type in ['sn', 'M-sn', 'M-(sn-H2O)']:
                    _frag_mz = _fa_se[_frag_type]
                    _frag_mz_low = _frag_mz - _frag_mz * ms2_precision
                    _frag_mz_high = _frag_mz + _frag_mz * ms2_precision
                    _frag_mz_query_code = '%f <= mz <= %f' % (_frag_mz_low, _frag_mz_high)

                    _frag_df = ms2_df.query(_frag_mz_query_code)

                    if _frag_df.shape[0] > 0:
                        _frag_df.loc[:, 'ppm'] = 1e6 * (_frag_df['mz'] - _frag_mz) / _frag_mz
                        _frag_df.loc[:, 'ppm_abs'] = _frag_df['ppm'].abs()
                        _frag_df.loc[:, 'FA'] = _fa_abbr

                        if _frag_df.shape[0] > 1:
                            _frag_i_df = _frag_df.sort_values(by='i', ascending=False).head(1)
                            _frag_ppm_df = _frag_df.sort_values(by='ppm_abs').head(1)
                            _frag_df = _frag_i_df
                            if _frag_ppm_df['i'].tolist() == _frag_i_df['i'].tolist():
                                pass
                            else:
                                _frag_df = _frag_df.append(_frag_ppm_df)

                        if _frag_type == 'sn':
                            _frag_df.loc[:, 'Lipid_species'] = '[FA%s-H]-' % _fa_abbr
                            fa_ident_df = fa_ident_df.append(_frag_df)
                        elif _frag_type == 'M-sn':
                            if _fa_link in lyso_fa_linker_dct.keys():
                                if bulk_fa_db - _fa_db >=0:
                                    _fa_lyso_link = lyso_fa_linker_dct[_fa_link]
                                    _frag_df.loc[:, 'Lipid_species'] = '[Lyso%s(%s%i:%i)-H2O-H]-' % (pl_typ,
                                                                                                     _fa_lyso_link,
                                                                                                     bulk_fa_c - _fa_c,
                                                                                                     bulk_fa_db - _fa_db
                                                                                                     )
                                    lyso_ident_df = lyso_ident_df.append(_frag_df)
                        elif _frag_type == 'M-(sn-H2O)':
                            if _fa_link in lyso_fa_linker_dct.keys():
                                if bulk_fa_db - _fa_db >= 0:
                                    _fa_lyso_link = lyso_fa_linker_dct[_fa_link]
                                    _frag_df.loc[:, 'Lipid_species'] = '[Lyso%s(%s%i:%i)-H]-' % (pl_typ, _fa_lyso_link,
                                                                                                 bulk_fa_c - _fa_c,
                                                                                                 bulk_fa_db - _fa_db
                                                                                                 )

                                    lyso_w_ident_df = lyso_w_ident_df.append(_frag_df)

        # format the output DataFrame
        if fa_ident_df.shape[0] > 0:
            fa_ident_df = fa_ident_df.query('i > %f' % ms2_threshold)
            fa_ident_df = fa_ident_df[['Lipid_species', 'FA', 'mz', 'i', 'ppm']].reset_index(drop=True)
            fa_ident_df = fa_ident_df.sort_values(by='i', ascending=False).head(10)

        if lyso_ident_df.shape[0] > 0:
            lyso_ident_df = lyso_ident_df.query('i > %f' % ms2_threshold)
            lyso_ident_df = lyso_ident_df[['Lipid_species', 'FA', 'mz', 'i', 'ppm']].reset_index(drop=True)
            lyso_ident_df = lyso_ident_df.sort_values(by='i', ascending=False).head(5)

        if lyso_w_ident_df.shape[0] > 0:
            lyso_w_ident_df = lyso_w_ident_df.query('i > %f' % ms2_threshold)
            lyso_w_ident_df = lyso_w_ident_df[['Lipid_species', 'FA', 'mz', 'i', 'ppm']].reset_index(drop=True)
            lyso_w_ident_df = lyso_w_ident_df.sort_values(by='i', ascending=False).head(5)

        return fa_ident_df, lyso_ident_df, lyso_w_ident_df

    def get_structure(self, abbr):

        lipid_abbr_lst = []
        lipid_sn1_lst = []
        lipid_sn2_lst = []
        db_sn1_lst = []
        db_sn2_lst = []

        lipid_info_dct = self.decode_abbr(abbr)
        pl_typ = lipid_info_dct['TYPE']
        bulk_fa_c = lipid_info_dct['C']
        bulk_fa_db = lipid_info_dct['DB']
        bulk_fa_linker = lipid_info_dct['LINK']

        for _i, _fa_se in self.fa_def_df.iterrows():
            # FA, Link, C, DB
            _fa_abbr = _fa_se['FA']
            _fa_link = _fa_se['Link']
            _fa_c = _fa_se['C']
            _fa_db = _fa_se['DB']

            if _fa_db <= bulk_fa_db and _fa_c <= bulk_fa_c:
                if _fa_link == bulk_fa_linker[0:1]:
                    _rest_fa_link = bulk_fa_linker[2]
                    _rest_fa_c = bulk_fa_c - _fa_c
                    _rest_fa_db = bulk_fa_db - _fa_db

                    _rest_fa_df = self.fa_def_df.query('Link == "%s" and C == %i and DB == %i'
                                                       % (_rest_fa_link, _rest_fa_c, _rest_fa_db)
                                                       )

                    if _rest_fa_df.shape[0] == 1:
                        # if _fa_link == 'A-':
                        #     _fa_link = ''
                        # else:
                        #     pass
                        _rest_fa_abbr = _rest_fa_df['FA'].tolist()[0]
                        lipid_abbr = '%s(%s_%s)' % (pl_typ, _fa_abbr, _rest_fa_abbr)

                        lipid_abbr_lst.append(lipid_abbr)
                        lipid_sn1_lst.append(_fa_abbr)
                        lipid_sn2_lst.append(_rest_fa_abbr)
                        db_sn1_lst.append(_fa_db)
                        db_sn2_lst.append(_rest_fa_db)

        lipid_abbr_df = pd.DataFrame(data={'Lipid_species': lipid_abbr_lst, 'sn1_abbr': lipid_sn1_lst,
                                           'sn2_abbr': lipid_sn2_lst, 'sn1_DB': db_sn1_lst, 'sn2_DB': db_sn2_lst})

        lipid_abbr_df = lipid_abbr_df.query('sn1_DB <=sn2_DB')
        lipid_abbr_df = lipid_abbr_df[['Lipid_species', 'sn1_abbr', 'sn2_abbr']]

        return lipid_abbr_df

    def get_match(self, abbr, charge_type, mz_lib, ms2_df, ms2_precision=500e-6, ms2_threshold=100):

        fa_ident_df, lyso_ident_df, lyso_w_ident_df = self.get_fa_search(abbr, charge_type, mz_lib, ms2_df,
                                                                         ms2_precision=ms2_precision,
                                                                         ms2_threshold=ms2_threshold
                                                                         )

        lipid_abbr_df = self.get_structure(abbr)

        weight_type_lst = ['sn1', 'sn2', 'M-sn1', 'M-sn2', 'M-(sn1-H2O)', 'M-(sn2-H2O)']
        weight_dct = {}
        for _type in weight_type_lst:
            lipid_abbr_df[_type] = 0
        # lipid_abbr_df['Score'] = 0
        if fa_ident_df.shape[0] > 0:
            fa_ident_lst = fa_ident_df['FA'].tolist()

            try:
                lyso_ident_lst = lyso_ident_df['FA'].tolist()
            except KeyError:
                lyso_ident_lst = []
            try:
                lyso_w_ident_lst = lyso_w_ident_df['FA'].tolist()
            except KeyError:
                lyso_w_ident_lst = []

            self.weight_df['mz'] = 0.0
            for _i, _weight_se in self.weight_df.iterrows():
                _type = _weight_se['Type']
                _weight = _weight_se['Weight']
                # _pr_h = _weight_se['[M-H]-']
                # _pr_formate = _weight_se['[M+HCOO]-']
                # _fa = _weight_se['FA']
                # _water = _weight_se['H2O']
                weight_dct[_type] = _weight

            for _i_abbr, _abbr_se in lipid_abbr_df.iterrows():
                # _pl_abbr = _abbr_se['Lipid_abbr']
                _sn1_abbr = _abbr_se['sn1_abbr']
                _sn2_abbr = _abbr_se['sn2_abbr']

                if _sn1_abbr in fa_ident_lst:
                    lipid_abbr_df.set_value(_i_abbr, 'sn1',
                                            weight_dct['sn1'] * (10 - fa_ident_lst.index(_sn1_abbr)) / 10)
                if _sn2_abbr in fa_ident_lst:
                    lipid_abbr_df.set_value(_i_abbr, 'sn2',
                                            weight_dct['sn2'] * (10 - fa_ident_lst.index(_sn2_abbr)) / 10)
                if _sn1_abbr in lyso_ident_lst:
                    lipid_abbr_df.set_value(_i_abbr, 'M-sn1',
                                            weight_dct['M-sn1'] * (10 - lyso_ident_lst.index(_sn1_abbr)) / 10)
                if _sn2_abbr in lyso_ident_lst:
                    lipid_abbr_df.set_value(_i_abbr, 'M-sn2',
                                            weight_dct['M-sn2'] * (10 - lyso_ident_lst.index(_sn2_abbr)) / 10)
                if _sn1_abbr in lyso_w_ident_lst:
                    lipid_abbr_df.set_value(_i_abbr, 'M-(sn1-H2O)',
                                            weight_dct['M-(sn1-H2O)'] * (10 - lyso_w_ident_lst.index(_sn1_abbr)) / 10)
                if _sn2_abbr in lyso_w_ident_lst:
                    lipid_abbr_df.set_value(_i_abbr, 'M-(sn2-H2O)',
                                            weight_dct['M-(sn2-H2O)'] * (10 - lyso_w_ident_lst.index(_sn2_abbr)) / 10)

                lipid_abbr_df['Score'] = lipid_abbr_df[weight_type_lst].sum(axis=1, numeric_only=True)
        else:
            print('!!!!!! NO FA identified =====>--> Skip >>> >>>')
        return lipid_abbr_df, fa_ident_df, lyso_ident_df, lyso_w_ident_df

    def get_specific_peaks(self, mz_lib, ms2_df, ms2_precision=50e-6, ms2_threshold=10):

        _target_frag_df = pd.DataFrame()
        _target_nl_df = pd.DataFrame()
        _other_frag_df = pd.DataFrame()
        _other_nl_df = pd.DataFrame()

        for _i, _frag_se in self.target_frag_df.iterrows():

            _frag_mz = _frag_se['EXACTMASS']
            _frag_class = _frag_se['CLASS']

            _frag_mz_low = _frag_mz - _frag_mz * ms2_precision
            _frag_mz_high = _frag_mz + _frag_mz * ms2_precision
            _frag_mz_query_code = '%f <= mz <= %f and i > %f' % (_frag_mz_low, _frag_mz_high, ms2_threshold)

            _frag_df = ms2_df.query(_frag_mz_query_code)

            if _frag_df.shape[0] > 0:
                _frag_df = _frag_df.sort_values(by='i', ascending=False)
                _frag_df.loc[:, 'CLASS'] = _frag_class
                _target_frag_df = _target_frag_df.append(_frag_df.head(1))

        for _i, _frag_se in self.other_frag_df.iterrows():

            _frag_mz = _frag_se['EXACTMASS']
            _frag_class = _frag_se['CLASS']

            _frag_mz_low = _frag_mz - _frag_mz * ms2_precision
            _frag_mz_high = _frag_mz + _frag_mz * ms2_precision
            _frag_mz_query_code = '%f <= mz <= %f and i > %f' % (_frag_mz_low, _frag_mz_high, ms2_threshold)

            _frag_df = ms2_df.query(_frag_mz_query_code)

            if _frag_df.shape[0] > 0:
                _frag_df = _frag_df.sort_values(by='i', ascending=False)
                _frag_df.loc[:, 'CLASS'] = _frag_class
                _other_frag_df = _other_frag_df.append(_frag_df.head(1))

        for _i, _nl_se in self.target_nl_df.iterrows():

            _nl_mz = _nl_se['EXACTMASS']
            _nl_class = _nl_se['CLASS']

            _nl_mz_low = mz_lib - _nl_mz - _nl_mz * ms2_precision
            _nl_mz_high = mz_lib - _nl_mz + _nl_mz * ms2_precision
            _nl_mz_query_code = '%f <= mz <= %f and i > %f' % (_nl_mz_low, _nl_mz_high, ms2_threshold)

            _nl_df = ms2_df.query(_nl_mz_query_code)

            if _nl_df.shape[0] > 0:
                _nl_df = _nl_df.sort_values(by='i', ascending=False)
                _nl_df.loc[:, 'CLASS'] = _nl_class
                _target_nl_df = _target_nl_df.append(_nl_df.head(1))

        for _i, _nl_se in self.other_nl_df.iterrows():

            _nl_mz = _nl_se['EXACTMASS']
            _nl_class = _nl_se['CLASS']

            _nl_mz_low = mz_lib - _nl_mz - _nl_mz * ms2_precision
            _nl_mz_high = mz_lib - _nl_mz + _nl_mz * ms2_precision
            _nl_mz_query_code = '%f <= mz <= %f and i > %f' % (_nl_mz_low, _nl_mz_high, ms2_threshold)

            _nl_df = ms2_df.query(_nl_mz_query_code)

            if _nl_df.shape[0] > 0:
                _nl_df = _nl_df.sort_values(by='i', ascending=False)
                _nl_df.loc[:, 'CLASS'] = _nl_class
                _other_nl_df = _other_nl_df.append(_nl_df.head(1))

        return _target_frag_df, _target_nl_df, _other_frag_df, _other_nl_df


if __name__ == '__main__':
    fa_list_csv = r'D:\LPPhunter\FA_list.csv'
    score_cfg = r'D:\LPPhunter\Score_cfg.xlsx'

    usr_fa_def_df = pd.read_csv(fa_list_csv)
    usr_weight_df = pd.read_excel(score_cfg)
    usr_fa_def_df['C'] = usr_fa_def_df['C'].astype(int)
    usr_fa_def_df['DB'] = usr_fa_def_df['DB'].astype(int)

    # for usr_abbr in ['PC(O-36:3)', 'PC(P-38:4)', 'PC(40:5)']:
    #     struc_df = get_structure(usr_fa_def_df, usr_abbr)
    #     print(struc_df)

    usr_lp = 'PC(40:5)'
    score_calc = ScoreGenerator(usr_fa_def_df, usr_weight_df)
    struc_df = score_calc.get_structure(usr_lp)
    print(struc_df)

    usr_fa_ident_df = usr_fa_def_df.head(8)
    usr_lyso_ident_df = usr_fa_def_df.head(4)
    usr_lyso_w_ident_df = usr_fa_def_df.head(2)

    fa_abbr_lst = ['PC(40:5)', '20:5', 'O-18:1', 'P-16:0']
    for _fa in fa_abbr_lst:
        print(_fa)
        fa_info_dct = score_calc.decode_abbr(_fa)
        print(fa_info_dct)
