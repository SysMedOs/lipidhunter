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

from __future__ import division

import re
from LibLipidHunter.AbbrElemCalc import BulkAbbrFormula
import pandas as pd


class ScoreGenerator:
    def __init__(self, fa_def_df, weight_df, key_frag_df, lipid_type, ion_charge='[M-H]-'):
        self.fa_def_df = fa_def_df
        self.weight_df = weight_df
        self.target_frag_df = key_frag_df.query(r'CLASS == "%s" and TYPE == "FRAG" and PR_CHARGE == "%s"'
                                                % (lipid_type, ion_charge))
        self.target_nl_df = key_frag_df.query(r'CLASS == "%s" and TYPE == "NL" and PR_CHARGE == "%s"'
                                              % (lipid_type, ion_charge))

        if ion_charge in ['[M-H]-', '[M+HCOO]-', '[M+CH3COO]-', '[M+FA-H]-', '[M+OAc]-']:
            charge_mode = 'NEG'
        elif ion_charge in ['[M+H]+', '[M+NH4]+']:
            charge_mode = 'POS'
        else:
            charge_mode = 'NEG'

        self.other_frag_df = key_frag_df.query('CLASS != "%s" and TYPE == "FRAG" and CHARGE_MODE == "%s"'
                                               % (lipid_type, charge_mode))
        self.other_nl_df = key_frag_df.query('CLASS != "%s" and TYPE == "NL" and CHARGE_MODE == "%s"'
                                             % (lipid_type, charge_mode))
        self.lipid_type = lipid_type

    @staticmethod
    def get_pr_mz(charge_type, mz_lib):

        pr_mz = 0.0

        if charge_type in ['[M-H]-', '[M+HCOO]-', '[M+FA-H]-', '[M+CH3COO]-', '[M+OAc]-', '[M+AcOH-H]-']:
            charge_mode = 'NEG'
            if charge_type == '[M-H]-':
                pr_mz = mz_lib
            elif charge_type in ['[M+HCOO]-', '[M+FA-H]-']:
                pr_mz = mz_lib - 46.005480  # - HCOOH
            elif charge_type in ['[M+CH3COO]-', '[M+OAc]-']:
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
        fa_checker = re.compile(r'(\d{1,2})([:])(\d{1,2})')
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
            print('PL')
            pl_re_chk = pl_checker.match(abbr)
            pl_typ_lst = pl_re_chk.groups()
            _pl_typ = pl_typ_lst[0]
            bulk_fa_typ = pl_typ_lst[2]
        if pip_checker.match(abbr):
            print('PIP')
            pip_re_chk = pip_checker.match(abbr)
            pip_typ_lst = pip_re_chk.groups()
            _pl_typ = pip_typ_lst[0]
            bulk_fa_typ = pip_typ_lst[2]
        if tg_checker.match(abbr):
            print('TG')
            tg_re_chk = tg_checker.match(abbr)
            tg_typ_lst = tg_re_chk.groups()
            _pl_typ = tg_typ_lst[0]
            bulk_fa_typ = tg_typ_lst[2]
        if fa_checker.match(abbr):
            print('FA')
            _pl_typ = 'FA'
            bulk_fa_typ = abbr
        if fa_o_checker.match(abbr):
            print('FA')
            _pl_typ = 'FA'
            bulk_fa_typ = abbr
        if fa_p_checker.match(abbr):
            print('FA')
            _pl_typ = 'FA'
            bulk_fa_typ = abbr

        print(bulk_fa_typ)

        if fa_checker.match(bulk_fa_typ):
            if _pl_typ == "TG":
                bulk_fa_linker = 'A-A-A-'
            else:
                bulk_fa_linker = 'A-A-'
            lyso_fa_linker_dct = {'A': ''}
            fa_chk = fa_checker.match(bulk_fa_typ)
            bulk_fa_lst = fa_chk.groups()
            bulk_fa_c = bulk_fa_lst[0]
            bulk_fa_db = bulk_fa_lst[2]
        elif fa_o_checker.match(bulk_fa_typ):
            if _pl_typ == "TG":
                bulk_fa_linker = 'O-A-A-'
            else:
                bulk_fa_linker = 'O-A-'
            lyso_fa_linker_dct = {'O': '', 'A': 'O-'}  # link of the other sn after NL of this sn
            fa_chk = fa_o_checker.match(bulk_fa_typ)
            bulk_fa_lst = fa_chk.groups()
            bulk_fa_c = bulk_fa_lst[1]
            bulk_fa_db = bulk_fa_lst[3]
        elif fa_p_checker.match(bulk_fa_typ):
            if _pl_typ == "TG":
                bulk_fa_linker = 'P-A-A-'
            else:
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

    def get_fa_search(self, abbr, charge_type, mz_lib, ms2_df, ms2_precision=500e-6,
                      ms2_threshold=100, ms2_infopeak_threshold=0.02):

        fa_ident_df = pd.DataFrame()
        lyso_ident_df = pd.DataFrame()
        lyso_w_ident_df = pd.DataFrame()
        mg_w_ident_df = pd.DataFrame()

        lipid_info_dct = self.decode_abbr(abbr)
        pl_typ = lipid_info_dct['TYPE']
        bulk_fa_c = lipid_info_dct['C']
        bulk_fa_db = lipid_info_dct['DB']
        # bulk_fa_linker = lipid_info_dct['LINK']
        lyso_fa_linker_dct = lipid_info_dct['LYSO_LINK']

        # use the max threshold from abs & relative intensity settings
        ms2_basepeak_i = ms2_df['i'].max()
        ms2_info_i = ms2_basepeak_i * ms2_infopeak_threshold
        ms2_threshold = max(ms2_threshold, ms2_info_i)

        calc_pr_mz, charge_mode = self.get_pr_mz(charge_type, mz_lib)
        calc_pr_mz = mz_lib

        if abbr[:2] in ['PA', 'PC', 'PE', 'PG', 'PI', 'PS', 'SM']:
            lipid_type = 'PL'
        elif abbr[:2] in ['TA', 'TG', 'DA', 'DG', 'MA', 'MG']:
            lipid_type = 'GL'
        else:
            lipid_type = 'PL'

        if lipid_type == 'PL' and charge_mode == 'NEG':
            print ('negative)')
            fa_chk_df = self.fa_def_df[['FA', 'Link', 'C', 'DB', 'mass', '[M-H]-', 'NL-H2O']]
            fa_chk_df = fa_chk_df.rename(columns={'[M-H]-': 'sn', 'mass': 'NL'})

            if abbr[:2] == 'PC' and charge_type == '[M+HCOO]-':
                fa_chk_df['[M-H]-sn'] = calc_pr_mz - fa_chk_df['NL-H2O'] - 60.021130  # - CH3COOH for PC
                fa_chk_df['[M-H]-sn-H2O'] = calc_pr_mz - fa_chk_df['NL'] - 60.021130  # - CH3COOH for PC
                fa_chk_df['Proposed_structures'] = ''
                lyso_hg_mod = '-CH3'
            elif abbr[:2] == 'PC' and charge_type == '[M+OAc]-':
                fa_chk_df['[M-H]-sn'] = calc_pr_mz - fa_chk_df['NL-H2O'] - 74.036780  # - CH3COOCH3 for PC
                fa_chk_df['[M-H]-sn-H2O'] = calc_pr_mz - fa_chk_df['NL'] - 74.036780  # - CH3COOCH3 for PC
                fa_chk_df['Proposed_structures'] = ''
                lyso_hg_mod = '-CH3'

            elif abbr[:2] == 'PS' and charge_type == '[M-H]-':
                fa_chk_df['[M-H]-sn'] = calc_pr_mz - fa_chk_df['NL-H2O'] - 87.032029  # - C3H5NO2 for PS
                fa_chk_df['[M-H]-sn-H2O'] = calc_pr_mz - fa_chk_df['NL'] - 87.032029  # - C3H5NO2 for PS
                fa_chk_df['Proposed_structures'] = ''
                lyso_hg_mod = '-87(Ser)'

            else:
                # Loss of FA-18, -OH remains on Glycerol back bone
                fa_chk_df['[M-H]-sn'] = calc_pr_mz - fa_chk_df['NL-H2O']
                # Loss of FA as full acid, -OH remains on FA NL
                fa_chk_df['[M-H]-sn-H2O'] = calc_pr_mz - fa_chk_df['NL']
                fa_chk_df['Proposed_structures'] = ''
                lyso_hg_mod = ''

            fa_abbr_lst = fa_chk_df['FA'].tolist()

            for _i, _fa_se in fa_chk_df.iterrows():

                _fa_abbr = _fa_se['FA']
                _fa_link = _fa_se['Link']
                _fa_c = _fa_se['C']
                _fa_db = _fa_se['DB']

                for _frag_type in ['sn', '[M-H]-sn', '[M-H]-sn-H2O']:
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
                            _frag_df = _frag_i_df.copy()
                            if _frag_ppm_df['i'].tolist() == _frag_i_df['i'].tolist():
                                pass
                            else:
                                _frag_df = _frag_i_df.append(_frag_ppm_df)

                        if _frag_type == 'sn':
                            if _fa_link != 'O' and _fa_link != 'P':
                                _frag_df.loc[:, 'Proposed_structures'] = 'FA %s [M-H]-' % _fa_abbr
                                fa_ident_df = fa_ident_df.append(_frag_df)
                        elif _frag_type == '[M-H]-sn':
                            if _fa_link in lyso_fa_linker_dct.keys():
                                if bulk_fa_db - _fa_db >= 0:
                                    _fa_lyso_link = lyso_fa_linker_dct[_fa_link]
                                    _fa_lyso_str = '%s%i:%i' % (_fa_lyso_link, bulk_fa_c - _fa_c, bulk_fa_db - _fa_db)
                                    # keep theoretical common Lyso PL only
                                    if _fa_lyso_str in fa_abbr_lst:
                                        _frag_df.loc[:, 'Proposed_structures'] = ('L%s %s [M%s-H]-'
                                                                                  % (pl_typ, _fa_lyso_str, lyso_hg_mod)
                                                                                  )
                                        lyso_ident_df = lyso_ident_df.append(_frag_df)
                        elif _frag_type == '[M-H]-sn-H2O':
                            if _fa_link in lyso_fa_linker_dct.keys():
                                if bulk_fa_db - _fa_db >= 0:
                                    _fa_lyso_link = lyso_fa_linker_dct[_fa_link]
                                    _fa_lyso_str = '%s%i:%i' % (_fa_lyso_link, bulk_fa_c - _fa_c, bulk_fa_db - _fa_db)
                                    # keep theoretical common Lyso PL only
                                    if _fa_lyso_str in fa_abbr_lst:
                                        _frag_df.loc[:, 'Proposed_structures'] = ('L%s %s [M%s-H2O-H]-'
                                                                                  % (pl_typ, _fa_lyso_str, lyso_hg_mod)
                                                                                  )
                                        lyso_w_ident_df = lyso_w_ident_df.append(_frag_df)

        elif lipid_type == 'GL' and charge_mode == 'POS':

            print ('positive')
            if charge_type == '[M+NH4]+':
                fa_chk_df = self.fa_def_df[['FA', 'Link', 'C', 'DB', 'mass', '[M+H]+', 'NL-H2O']]
                fa_chk_df = fa_chk_df.rename(columns={'[M+H]+': 'sn', 'mass': 'NL'})
                fa_chk_df['[M+H]-sn'] = calc_pr_mz - fa_chk_df['NL-H2O'] - 17.026549  # - NH3 adduct
                fa_chk_df['[M+H]-sn-H2O'] = calc_pr_mz - fa_chk_df['NL'] - 17.026549  # - NH3 adduct
                fa_chk_df['Proposed_structures'] = ''
                lyso_hg_mod = '-NH3'
            elif charge_type == '[M+Na]+':
                fa_chk_df = self.fa_def_df[['FA', 'Link', 'C', 'DB', 'mass', '[M+H]+', 'NL-H2O', 'NL+Na']]
                fa_chk_df = fa_chk_df.rename(columns={'[M+H]+': 'sn', 'mass': 'NL'})
                fa_chk_df['[M+Na]-RCOOH'] = calc_pr_mz - fa_chk_df['NL']
                fa_chk_df['[M+H]-RCOONa'] = calc_pr_mz - fa_chk_df['NL+Na']
                fa_chk_df['Proposed_structures'] = ''
                lyso_hg_mod = ''
                print ('This is the calculate precursor')
                print calc_pr_mz
            else:
                fa_chk_df = self.fa_def_df[['FA', 'Link', 'C', 'DB', 'mass', '[M+H]+', 'NL-H2O']]
                fa_chk_df = fa_chk_df.rename(columns={'[M+H]+': 'sn', 'mass': 'NL'})
                fa_chk_df['[M+H]-sn'] = calc_pr_mz - fa_chk_df['NL-H2O']
                fa_chk_df['[M+H]-sn-H2O'] = calc_pr_mz - fa_chk_df['NL']
                fa_chk_df['Proposed_structures'] = ''
                lyso_hg_mod = ''

            fa_abbr_lst = fa_chk_df['FA'].tolist()
            print ('This is the calculate precursor')
            print (calc_pr_mz)
            for _i, _fa_se in fa_chk_df.iterrows():

                _fa_abbr = _fa_se['FA']
                _fa_link = _fa_se['Link']
                _fa_c = _fa_se['C']
                _fa_db = _fa_se['DB']
                if charge_type in ['[M+H]+', '[M+NH4]+']:
                    for _frag_type in ['sn', '[M+H]-sn', '[M+H]-sn-H2O', '[M+H]-2sn-H2O']:
                        if _frag_type is not '[M+H]-2sn-H2O':
                            _frag_mz = _fa_se[_frag_type]
                            _frag_mz_low = _frag_mz - _frag_mz * ms2_precision
                            _frag_mz_high = _frag_mz + _frag_mz * ms2_precision
                            _frag_mz_query_code = '%f <= mz <= %f' % (_frag_mz_low, _frag_mz_high)
                            _frag_df = ms2_df.query(_frag_mz_query_code)

                            if _frag_df.shape[0] > 0 :

                                _frag_df.loc[:, 'ppm'] = 1e6 * (_frag_df['mz'] - _frag_mz) / _frag_mz
                                _frag_df.loc[:, 'ppm_abs'] = _frag_df['ppm'].abs()
                                _frag_df.loc[:, 'FA'] = _fa_abbr

                                if _frag_df.shape[0] > 1:

                                    _frag_i_df = _frag_df.sort_values(by='i', ascending=False).head(1)
                                    _frag_ppm_df = _frag_df.sort_values(by='ppm_abs').head(1)
                                    _frag_df = _frag_i_df.copy()
                                    if _frag_ppm_df['i'].tolist() == _frag_i_df['i'].tolist():
                                        pass
                                    else:
                                        _frag_df = _frag_i_df.append(_frag_ppm_df)
                                if _frag_type == 'sn':
                                    if _fa_link != 'O' and _fa_link != 'P':
                                        _frag_df.loc[:, 'Proposed_structures'] = 'FA %s [M+H]+' % _fa_abbr
                                        fa_ident_df = fa_ident_df.append(_frag_df)
                                elif _frag_type == '[M+H]-sn':
                                    if _fa_link in lyso_fa_linker_dct.keys():
                                        if bulk_fa_db - _fa_db >= 0:
                                            _fa_lyso_link = lyso_fa_linker_dct[_fa_link]
                                            _fa_lyso_str = '%s%i:%i' % (_fa_lyso_link, bulk_fa_c - _fa_c, bulk_fa_db - _fa_db)
                                            # keep theoretical common Lyso PL only
                                            #if _fa_lyso_str in fa_abbr_lst:
                                            _frag_df.loc[:, 'Proposed_structures'] = ('DG %s [M%s+H]+'
                                                                                          % (_fa_lyso_str, lyso_hg_mod)
                                                                                          )
                                            lyso_ident_df = lyso_ident_df.append(_frag_df)
                                elif _frag_type == '[M+H]-sn-H2O':

                                    if _fa_link in lyso_fa_linker_dct.keys():
                                        if bulk_fa_db - _fa_db >= 0:
                                            _fa_lyso_link = lyso_fa_linker_dct[_fa_link]
                                            _fa_lyso_str = '%s%i:%i' % (_fa_lyso_link, bulk_fa_c - _fa_c, bulk_fa_db - _fa_db)
                                            # keep theoretical common Lyso PL only
                                            #if _fa_lyso_str in fa_abbr_lst:

                                            _frag_df.loc[:, 'Proposed_structures'] = ('DG %s [M%s-H2O+H]+'
                                                                                          % (_fa_lyso_str, lyso_hg_mod)
                                                                                          )
                                            lyso_w_ident_df = lyso_w_ident_df.append(_frag_df)
                        elif _frag_type == '[M+H]-2sn-H2O':
                            for _i2, _fa_se2 in fa_chk_df.iterrows():
                                _fa_abbr2 = _fa_se2['FA']
                                _fa_link2 = _fa_se2['Link']
                                _fa_c2 = _fa_se2['C']
                                _fa_db2 = _fa_se2['DB']
                                print ('Is it here')
                                if not _fa_link == _fa_link2 == 'A':
                                    print ('Welll')
                                    _frag_mz2 = _fa_se['[M+H]-sn'] + _fa_se2['NL']
                                    print ('Maybe not')
                                    print ms2_precision

                                    _frag_mz_low2= _frag_mz2 - _frag_mz2*ms2_precision
                                    _frag_mz_high2 = _frag_mz2 + _frag_mz2*ms2_precision
                                    _frag_mz_query_code2='%f <= mz <= %f' % (_frag_mz_low2, _frag_mz_high2)

                                    _frag_df = ms2_df.query(_frag_mz_query_code2)

                                    if _frag_df.shape[0] > 0:
                                        _frag_df.loc[:, 'ppm'] = 1e6*(_frag_df['mz']-_frag_mz2)/_frag_mz2
                                        _frag_df.loc[:, 'ppm_abs'] = _frag_df['ppm'].abs()
                                        _frag_df.loc[:, 'FA'] = _fa_abbr
                                        _frag_df.loc[:, 'FA2'] = _fa_abbr2

                                        if _frag_df.shape[0] > 1:
                                            # print (_frag_df)
                                            _frag_i_df = _frag_df.sort_values(by='i', ascending=False).head(1)
                                            # print (_frag_i_df)
                                            _frag_ppm_df = _frag_df.sort_values(by='ppm_abs').head(1)
                                            # print (_frag_ppm_df)
                                            _frag_df = _frag_i_df.copy()
                                            # if _frag_ppm_df['i'].tolist() == _frag_i_df['i'].tolist():
                                            #     pass
                                            # else:
                                            #     print ('4')
                                            #     _frag_df = _frag_i_df.append(_frag_ppm_df)
                                            if _fa_link in lyso_fa_linker_dct.keys() and _fa_link2 in lyso_fa_linker_dct.keys():
                                                if bulk_fa_db - _fa_db - _fa_db2 >= 0:
                                                    if _fa_link == _fa_link2  and  _fa_link == 'A':
                                                        _fa_mg_link = 'P-'
                                                    else:
                                                        _fa_mg_link = ''
                                                    _fa_mg_str = '%s%i:%i' % (_fa_mg_link, bulk_fa_c - _fa_c - _fa_c2, bulk_fa_db - _fa_db - _fa_db2)
                                                    if _fa_mg_str in fa_abbr_lst:
                                                        # print ('So I guess somewhere around')
                                                        _frag_df.loc[:, 'Proposed_structures'] = ('MG%s %s [M%s+H]+' % (pl_typ, _fa_mg_str, lyso_hg_mod))
                                                        mg_w_ident_df = mg_w_ident_df.append(_frag_df)
                                                        print mg_w_ident_df
                                                        exit()
                elif charge_type in ['[M+Na]+']:
                    for _frag_type in ['sn', '[M+H]-RCOONa', '[M+Na]-RCOOH', '[M+H]-RCOOH-RCOONa']:
                        if _frag_type is not '[M+H]-RCOOH-RCOONa':
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
                                    _frag_df = _frag_i_df.copy()
                                    if _frag_ppm_df['i'].tolist() == _frag_i_df['i'].tolist():
                                        pass
                                    else:
                                        _frag_df = _frag_i_df.append(_frag_ppm_df)
                                if _frag_type == 'sn':
                                    if _fa_link != 'O' and _fa_link != 'P':
                                        _frag_df.loc[:, 'Proposed_structures'] = 'FA %s [M+H]+' % _fa_abbr
                                        fa_ident_df = fa_ident_df.append(_frag_df)
                                elif _frag_type == '[M+H]-RCOONa':
                                    if _fa_link in lyso_fa_linker_dct.keys():
                                        if bulk_fa_db - _fa_db >= 0:
                                            _fa_lyso_link = lyso_fa_linker_dct[_fa_link]
                                            _fa_lyso_str = '%s%i:%i' % (
                                            _fa_lyso_link, bulk_fa_c - _fa_c, bulk_fa_db - _fa_db)
                                            # keep theoretical common Lyso PL only
                                            # if _fa_lyso_str in fa_abbr_lst:
                                            _frag_df.loc[:, 'Proposed_structures'] = ('DG %s [M%s+H-RCOONa]+'
                                                                                      % (_fa_lyso_str, lyso_hg_mod)
                                                                                      )
                                            lyso_ident_df = lyso_ident_df.append(_frag_df)
                                elif _frag_type == '[M+Na]-RCOOH':

                                    if _fa_link in lyso_fa_linker_dct.keys():
                                        if bulk_fa_db - _fa_db >= 0:
                                            _fa_lyso_link = lyso_fa_linker_dct[_fa_link]
                                            _fa_lyso_str = '%s%i:%i' % (
                                            _fa_lyso_link, bulk_fa_c - _fa_c, bulk_fa_db - _fa_db)
                                            # keep theoretical common Lyso PL only
                                            # if _fa_lyso_str in fa_abbr_lst:

                                            _frag_df.loc[:, 'Proposed_structures'] = ('DG %s [M%s+Na-RCOOH]+'
                                                                                      % (_fa_lyso_str, lyso_hg_mod)
                                                                                      )
                                            lyso_w_ident_df = lyso_w_ident_df.append(_frag_df)
                        elif _frag_type == '[M+H]-RCOOH-RCOONa':
                            for _i2, _fa_se2 in fa_chk_df.iterrows():
                                _fa_abbr2 = _fa_se2['FA']
                                _fa_link2 = _fa_se2['Link']
                                _fa_c2 = _fa_se2['C']
                                _fa_db2 = _fa_se2['DB']
                                print ('Is it here')
                                if not _fa_link == _fa_link2 == 'A':
                                    print ('Welll')
                                    _frag_mz2 = _fa_se['[M+H]-RCOONa'] + _fa_se2['NL']
                                    print ('Maybe not')
                                    print ms2_precision

                                    _frag_mz_low2 = _frag_mz2 - _frag_mz2 * ms2_precision
                                    _frag_mz_high2 = _frag_mz2 + _frag_mz2 * ms2_precision
                                    _frag_mz_query_code2 = '%f <= mz <= %f' % (_frag_mz_low2, _frag_mz_high2)

                                    _frag_df = ms2_df.query(_frag_mz_query_code2)

                                    if _frag_df.shape[0] > 0:
                                        _frag_df.loc[:, 'ppm'] = 1e6 * (_frag_df['mz'] - _frag_mz2) / _frag_mz2
                                        _frag_df.loc[:, 'ppm_abs'] = _frag_df['ppm'].abs()
                                        _frag_df.loc[:, 'FA'] = _fa_abbr
                                        _frag_df.loc[:, 'FA2'] = _fa_abbr2

                                        if _frag_df.shape[0] > 1:
                                            # print (_frag_df)
                                            _frag_i_df = _frag_df.sort_values(by='i', ascending=False).head(1)
                                            # print (_frag_i_df)
                                            _frag_ppm_df = _frag_df.sort_values(by='ppm_abs').head(1)
                                            # print (_frag_ppm_df)
                                            _frag_df = _frag_i_df.copy()
                                            # if _frag_ppm_df['i'].tolist() == _frag_i_df['i'].tolist():
                                            #     pass
                                            # else:
                                            #     print ('4')
                                            #     _frag_df = _frag_i_df.append(_frag_ppm_df)
                                            if _fa_link in lyso_fa_linker_dct.keys() and _fa_link2 in lyso_fa_linker_dct.keys():
                                                if bulk_fa_db - _fa_db - _fa_db2 >= 0:
                                                    if _fa_link == _fa_link2 and _fa_link == 'A':
                                                        _fa_mg_link = 'P-'
                                                    else:
                                                        _fa_mg_link = ''
                                                    _fa_mg_str = '%s%i:%i' % (_fa_mg_link, bulk_fa_c - _fa_c - _fa_c2,
                                                                              bulk_fa_db - _fa_db - _fa_db2)
                                                    if _fa_mg_str in fa_abbr_lst:
                                                        # print ('So I guess somewhere around')
                                                        _frag_df.loc[:, 'Proposed_structures'] = (
                                                        'MG%s %s [M%s+H]+' % (pl_typ, _fa_mg_str, lyso_hg_mod))
                                                        mg_w_ident_df = mg_w_ident_df.append(_frag_df)
                                                        print mg_w_ident_df
                                                        exit()
        # format the output DataFrame
        if fa_ident_df.shape[0] > 0:
            fa_ident_df = fa_ident_df.query('i > %f' % ms2_threshold)
            # print ('FA identification')
            # print fa_ident_df
            # exit()
            proposed_str_lst = {}
            if lipid_type == 'GL':
                print ('oooooooooo please')
                # exit()
                fa_ident_df['Flag'] = 1
                fa_ident_df = fa_ident_df[['Proposed_structures', 'FA', 'mz', 'i',
                                           'ppm', 'ppm_abs', 'Flag']].reset_index(drop=True)
            else:
                fa_ident_df['Flag'] = 1
                fa_ident_df = fa_ident_df[['Proposed_structures', 'FA', 'mz', 'i',
                                           'ppm', 'ppm_abs', 'Flag']].reset_index(drop=True)
            fa_ident_df = fa_ident_df.sort_values(by=['i', 'ppm_abs'], ascending=[False, True])
            fa_ident_df = fa_ident_df.drop_duplicates(['FA'], keep='first')
            fa_ident_df = fa_ident_df.sort_values(by='i', ascending=False).head(10)
        # print ('Finish with the FA')
        if lyso_ident_df.shape[0] > 0:
            lyso_found_dct = {}
            lyso_ident_df = lyso_ident_df.query('i > %f' % ms2_threshold)
            # print ('Lyso Identification')
            # print lyso_ident_df
            # exit()
            if lipid_type == 'GL':
                # print ('Is there any more of this ')
                # exit()
                lyso_ident_df['Flag'] = 1
                lyso_ident_df = lyso_ident_df.loc[lyso_ident_df['Flag'] == 1][['Proposed_structures', 'FA', 'mz', 'i',
                                                                               'ppm', 'ppm_abs',
                                                                               'Flag']].reset_index(drop=True)
            else:
                lyso_ident_df['Flag'] = 1
                lyso_ident_df = lyso_ident_df.loc[lyso_ident_df['Flag'] == 1][['Proposed_structures', 'FA', 'mz', 'i',
                                                                               'ppm', 'ppm_abs',
                                                                               'Flag']].reset_index(drop=True)
            lyso_ident_df = lyso_ident_df.sort_values(by=['i', 'ppm_abs'], ascending=[False, True])
            lyso_ident_df = lyso_ident_df.drop_duplicates(['FA'], keep='first')
            lyso_ident_df = lyso_ident_df.sort_values(by='i', ascending=False).head(5)

        # print ('Finish with the Lyso')
        if lyso_w_ident_df.shape[0] > 0:
            lyso_w_dct = {}
            lyso_w_ident_df = lyso_w_ident_df.query('i > %f' % ms2_threshold).reset_index(drop=True)
            # print ('Lyso Water')
            # print lyso_w_ident_df
            # exit()
            if lipid_type == 'GL':
                # print ('realllyyyyyy now')
                # exit()
                lyso_w_ident_df['Flag'] = 1
                lyso_w_ident_df = lyso_w_ident_df[['Proposed_structures', 'FA', 'mz', 'i',
                                                   'ppm', 'ppm_abs', 'Flag']].reset_index(drop=True)
            else:
                lyso_w_ident_df['Flag'] = 1
                lyso_w_ident_df = lyso_w_ident_df[['Proposed_structures', 'FA', 'mz', 'i',
                                                   'ppm', 'ppm_abs', 'Flag']].reset_index(drop=True)
            lyso_w_ident_df = lyso_w_ident_df.sort_values(by=['i', 'ppm_abs'], ascending=[False, True])
            lyso_w_ident_df = lyso_w_ident_df.drop_duplicates(['FA'], keep='first')
            lyso_w_ident_df = lyso_w_ident_df.sort_values(by='i', ascending=False).head(5)

        # print ('Finish with the Lyso W')
        if mg_w_ident_df.shape[0] > 0:
            mg_w_ident_df = mg_w_ident_df.query('i > %f' % ms2_threshold)
            # print ('MG identification')
            # print mg_w_ident_df
            # exit()
            if lipid_type == 'GL':
                mg_w_ident_df['Flag'] = 1
                mg_w_ident_df = mg_w_ident_df.loc[mg_w_ident_df['Flag'] == 1][['Proposed_structures', 'FA', 'FA2', 'mz', 'i', 'ppm', 'ppm_abs', 'Flag']].reset_index(drop = True)
            else:
                mg_w_ident_df['Flag'] = 1
                mg_w_ident_df =mg_w_ident_df.loc[mg_w_ident_df['Flag'] == 1][['Proposed_structures', 'FA', 'mz', 'i', 'ppm', 'ppm_abs', 'Flag']].reset_index(drop = True)

            mg_w_ident_df = mg_w_ident_df.sort_values(by=['i', 'ppm_abs'], ascending=[False, True])
            mg_w_ident_df = mg_w_ident_df.drop_duplicates(['FA'], keep='first')
            mg_w_ident_df = mg_w_ident_df.sort_values(by='i', ascending=False).head(5)

        print ('Values return')
        return fa_ident_df, lyso_ident_df, lyso_w_ident_df, mg_w_ident_df

    def get_structure(self, abbr):

        lipid_abbr_lst = []
        lipid_sn1_lst = []
        lipid_sn2_lst = []
        lipid_sn3_lst = []
        db_sn1_lst = []
        db_sn2_lst = []
        db_sn3_lst = []
        #abbr='TG(46:3)'
        #abbr = 'TG(O-58:8)'
        print(abbr)
        lipid_info_dct = self.decode_abbr(abbr)
        pl_typ = lipid_info_dct['TYPE']
        bulk_fa_c = lipid_info_dct['C']
        bulk_fa_db = lipid_info_dct['DB']
        bulk_fa_linker = lipid_info_dct['LINK']
        #print ('kouuuukou')
        #print lipid_info_dct
        #print ('hehehhehehehe')
        #print self.fa_def_df
        #exit()

        if abbr[:2] in ['PE', 'PA', 'PC', 'PI', 'PS', 'PG']:
            for _i, _fa_se in self.fa_def_df.iterrows():
                # FA, Link, C, DB
                _fa_abbr = _fa_se['FA']
                _fa_link = _fa_se['Link']
                _fa_c = _fa_se['C']
                _fa_db = _fa_se['DB']

                if _fa_db <= bulk_fa_db and _fa_c <= bulk_fa_c:
                    if _fa_link == bulk_fa_linker[0:1]:
                        print bulk_fa_linker
                        # print bulk_fa_linker[0:1]
                        # exit()
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
            lipid_abbr_df = pd.DataFrame(data={'Proposed_structures': lipid_abbr_lst, 'sn1_abbr': lipid_sn1_lst,
                                               'sn2_abbr': lipid_sn2_lst, 'sn1_DB': db_sn1_lst, 'sn2_DB': db_sn2_lst})

            lipid_abbr_df = lipid_abbr_df.query('sn1_DB <=sn2_DB')
            lipid_abbr_df = lipid_abbr_df[['Proposed_structures', 'sn1_abbr', 'sn2_abbr']]
            # print lipid_abbr_df
            # exit()

            return lipid_abbr_df

        elif abbr[:2] in ['TG'] and bulk_fa_linker not in ['O', 'P']:
            allsnList=[]
            print bulk_fa_linker
            print bulk_fa_linker[0:1]
            print bulk_fa_linker[2]
            print bulk_fa_linker[4]
            for _i, _fa_se in self.fa_def_df.iterrows():
                # FA, Link, C, DB
                _fa_abbr = _fa_se['FA']
                _fa_link = _fa_se['Link']
                _fa_c = _fa_se['C']
                _fa_db = _fa_se['DB']
                for _i2, _fa_se2 in self.fa_def_df.iterrows():
                    if _i2 >= _i:
                        _fa_abbr2 = _fa_se2['FA']
                        _fa_link2 = _fa_se2['Link']
                        _fa_c2 = _fa_se2['C']
                        _fa_db2 = _fa_se2['DB']
                        _fa_db_sum = _fa_db + _fa_db2
                        _fa_c_sum = _fa_c + _fa_c2
                        if _fa_db_sum <= bulk_fa_db and _fa_c_sum <= bulk_fa_c:
                            if _fa_link == bulk_fa_linker[0:1] and _fa_link2 == bulk_fa_linker[2]:
                                _rest_fa_link = bulk_fa_linker[4]

                                _rest_fa_c = bulk_fa_c - _fa_c_sum
                                _rest_fa_db = bulk_fa_db - _fa_db_sum

                                _rest_fa_df = self.fa_def_df.query('Link == "%s" and C == %i and DB == %i'
                                                                   % (_rest_fa_link, _rest_fa_c, _rest_fa_db)
                                                                   )

                                if _rest_fa_df.shape[0] == 1:
                                    # if _fa_link == 'A-':
                                    #     _fa_link = ''
                                    # else:
                                    #     pass

                                    _rest_fa_abbr = _rest_fa_df['FA'].tolist()[0]
                                    lipid_abbr = '%s(%s_%s_%s)' % (pl_typ, _fa_abbr, _fa_abbr2, _rest_fa_abbr)
                                    newEntry=[_fa_abbr, _fa_abbr2, _rest_fa_abbr]
                                    if sorted(newEntry) not in allsnList:
                                        allsnList.append(sorted(newEntry))
                                        lipid_abbr_lst.append(lipid_abbr)
                                        lipid_sn1_lst.append(_fa_abbr)
                                        lipid_sn2_lst.append(_fa_abbr2)
                                        lipid_sn3_lst.append(_rest_fa_abbr)
                                        db_sn1_lst.append(_fa_db)
                                        db_sn2_lst.append(_fa_db2)
                                        db_sn3_lst.append(_rest_fa_db)
                                    else:
                                        pass
            lipid_abbr_df = pd.DataFrame(data={'Proposed_structures': lipid_abbr_lst, 'sn1_abbr': lipid_sn1_lst,
                                               'sn2_abbr': lipid_sn2_lst, 'sn3_abbr': lipid_sn3_lst, 'sn1_DB': db_sn1_lst, 'sn2_DB': db_sn2_lst, 'sn3_DB': db_sn3_lst})

            #lipid_abbr_df = lipid_abbr_df.query('sn1_DB <=sn2_DB')
            lipid_abbr_df = lipid_abbr_df[['Proposed_structures', 'sn1_abbr', 'sn2_abbr', 'sn3_abbr']]

            # exit()
            return lipid_abbr_df
            # print "Daisy"
            # exit()

        else:

            lipid_sn3_lst = []
            db_sn3_lst = []
            # _fa_compination_3 = []

            lipid_abbr_df = pd.DataFrame(data={'Proposed_structures': lipid_abbr_lst, 'sn1_abbr': lipid_sn1_lst,
                                               'sn2_abbr': lipid_sn2_lst, 'sn3_abbr': lipid_sn3_lst,
                                               'sn1_DB': db_sn1_lst,
                                               'sn2_DB': db_sn2_lst, 'sn3_DB': db_sn3_lst})

            lipid_abbr_df = lipid_abbr_df[['Proposed_structures', 'sn1_abbr', 'sn2_abbr', 'sn3_abbr']]

            return lipid_abbr_df

    def get_match(self, abbr, charge_type, mz_lib, ms2_df, ms2_precision=500e-6,
                  ms2_threshold=100, ms2_infopeak_threshold=0.02, rank_mode=True):

        match_reporter = 0
        ms2_max_i = ms2_df['i'].max()

        fa_ident_df, lyso_ident_df, lyso_w_ident_df, mg_w_ident_df = self.get_fa_search(abbr, charge_type, mz_lib, ms2_df,
                                                                         ms2_precision=ms2_precision,
                                                                         ms2_threshold=ms2_threshold,
                                                                         ms2_infopeak_threshold=ms2_infopeak_threshold
                                                                         )
        print ('So, did it go')
        lipid_abbr_df = self.get_structure(abbr)
        sodiumFlag=0
        if abbr[:2] in ['TG']:
            # print "Kipors"
            # exit()
            if charge_type == '[M+Na]+':
                sodiumFlag=1
                weight_type_lst = ['sn1', 'sn2', 'sn3',  '[M+H]-sn1', '[M+H]-sn2', '[M+H]-sn3',
                               '[M+Na]-sn1', '[M+Na]-sn2', '[M+Na]-sn3', '[M+H]-(sn1+sn2)-H2O', '[M+H]-(sn1+sn3)-H2O', '[M+H]-(sn2+sn3)-H2O']
            else:
                weight_type_lst = ['sn1', 'sn2', 'sn3',  '[M+H]-sn1', '[M+H]-sn2', '[M+H]-sn3',
                               '[M+H]-sn1-H2O', '[M+H]-sn2-H2O', '[M+H]-sn3-H2O', '[M+H]-(sn1+sn2)-H2O', '[M+H]-(sn1+sn3)-H2O', '[M+H]-(sn2+sn3)-H2O']
        else:
            weight_type_lst = ['sn1', 'sn2', '[M-H]-sn1', '[M-H]-sn2',
                               '[M-H]-sn1-H2O', '[M-H]-sn2-H2O']
        weight_dct = {}
        for _type in weight_type_lst:
            lipid_abbr_df[_type] = 0
        if fa_ident_df.shape[0] > 0 or lyso_ident_df.shape[0] > 0 or lyso_w_ident_df.shape[0] > 0 or mg_w_ident_df.shape[0] > 0:
            combine_all_lst = pd.DataFrame()
            try:
                fa_ident_df['Type'] = 'FA'
                fa_ident_lst = fa_ident_df.loc[fa_ident_df['Flag'] == 1]['FA'].tolist()
                fa_i_lst = fa_ident_df.loc[fa_ident_df['Flag'] == 1]['i'].tolist()
                combine_all_lst = combine_all_lst.append(fa_ident_df.loc[fa_ident_df['Flag'] == 1])
            except KeyError:
                fa_ident_lst = []
                fa_i_lst = []

            try:
                lyso_ident_df['Type'] = 'Lyso'
                lyso_ident_lst = lyso_ident_df.loc[lyso_ident_df['Flag'] == 1]['FA'].tolist()
                lyso_i_lst = lyso_ident_df.loc[lyso_ident_df['Flag'] == 1]['i'].tolist()
                combine_all_lst = combine_all_lst.append(lyso_ident_df)
            except KeyError:
                lyso_ident_lst = []
                lyso_i_lst = []
            # print lyso_ident_df
            try:
                lyso_w_ident_df['Type'] = 'LysoW'
                lyso_w_ident_lst = lyso_w_ident_df.loc[lyso_w_ident_df['Flag'] == 1]['FA'].tolist()
                lyso_w_i_lst = lyso_w_ident_df.loc[lyso_w_ident_df['Flag'] == 1]['i'].tolist()
                combine_all_lst = combine_all_lst.append(
                    lyso_w_ident_df[['Proposed_structures', 'FA', 'mz', 'i', 'ppm', 'ppm_abs', 'Flag', 'Type']])
            except KeyError:
                lyso_w_ident_lst = []
                lyso_w_i_lst = []
            print mg_w_ident_df
            if mg_w_ident_df.shape[0] > 0 :
                print mg_w_ident_df
            mg_w_ident_lst2 = []
            ##########################################
            #
            # Check maybe this is were the problem is and thats why there is no output for the MG
            #
            ################################################
            try:
                mg_w_ident_df['Type'] = 'MG'

                mg_w_ident_lst = mg_w_ident_df.loc[mg_w_ident_df['Flag'] == 1][['FA', 'FA2']]

                count=0
                for _i, _g in mg_w_ident_df.iterrows():
                    mg_w_ident_lst2.append([_g['FA'], _g['FA2']])
                print mg_w_ident_df
                print lipid_abbr_df
                exit()
                # if lipid_abbr_df.shape[0] > 0:
                #     print mg_w_ident_lst2
                #     print ('Lets see')
                #     print mg_w_ident_lst
                #     print lipid_abbr_df

                #mg_w_ident_lst2 = mg_w_ident_df.loc[mg_w_ident_df['Flag'] == 1]['FA2'].tolist()
                mg_w_i_lst = mg_w_ident_df.loc[mg_w_ident_df['Flag'] == 1]['i'].tolist()
            except KeyError:
                mg_w_ident_lst = ()
                mg_w_ident_lst2 = ()
                mg_w_i_lst = ()
            # print lyso_w_ident_df
            self.weight_df['mz'] = 0.0
            for _i, _weight_se in self.weight_df.iterrows():
                _type = _weight_se['Type']
                _weight = _weight_se['Weight']
                weight_dct[_type] = _weight
            if abbr[:2] in ['TG']:
                #print "Miria"
                # exit()
                # print lipid_abbr_df
                for _i_abbr, _abbr_se in lipid_abbr_df.iterrows():
                    # _pl_abbr = _abbr_se['Lipid_abbr']
                    _sn1_abbr = _abbr_se['sn1_abbr']
                    _sn2_abbr = _abbr_se['sn2_abbr']
                    _sn3_abbr = _abbr_se['sn3_abbr']

                    if _sn1_abbr in fa_ident_lst:
                        print fa_ident_df
                        _rank_sn1 = fa_ident_lst.index(_sn1_abbr)
                        r_sn1_i = 100 * fa_i_lst[_rank_sn1] / ms2_max_i
                        lipid_abbr_df.set_value(_i_abbr, 'i_sn1', r_sn1_i)
                        if rank_mode is True:
                            lipid_abbr_df.set_value(_i_abbr, 'sn1', weight_dct['sn1'] * (10 - _rank_sn1) / 10)
                        else:
                            lipid_abbr_df.set_value(_i_abbr, 'sn1', weight_dct['sn1'] * r_sn1_i * 0.01)

                    if _sn2_abbr in fa_ident_lst:
                        _rank_sn2 = fa_ident_lst.index(_sn2_abbr)
                        r_sn2_i = 100 * fa_i_lst[_rank_sn2] / ms2_max_i
                        lipid_abbr_df.set_value(_i_abbr, 'i_sn2', r_sn2_i)
                        if rank_mode is True:
                            lipid_abbr_df.set_value(_i_abbr, 'sn2', weight_dct['sn2'] * (10 - _rank_sn2) / 10)
                        else:
                            lipid_abbr_df.set_value(_i_abbr, 'sn2', weight_dct['sn2'] * r_sn2_i * 0.01)

                    if _sn3_abbr in fa_ident_lst:
                        _rank_sn3 = fa_ident_lst.index(_sn3_abbr)
                        r_sn3_i = 100 * fa_i_lst[_rank_sn3] / ms2_max_i
                        lipid_abbr_df.set_value(_i_abbr, 'i_sn3', r_sn3_i)
                        if rank_mode is True:
                            lipid_abbr_df.set_value(_i_abbr, 'sn3', weight_dct['sn3'] * (10 - _rank_sn3) / 10)
                        else:
                            lipid_abbr_df.set_value(_i_abbr, 'sn3', weight_dct['sn3'] * r_sn3_i * 0.01)
                    if _sn1_abbr in lyso_ident_lst:
                        _rank_l_sn1 = lyso_ident_lst.index(_sn1_abbr)
                        r_lyso1_i = 100 * lyso_i_lst[_rank_l_sn1] / ms2_max_i
                        lipid_abbr_df.set_value(_i_abbr, 'i_[M+H]-sn1', r_lyso1_i)
                        if rank_mode is True:
                            lipid_abbr_df.set_value(_i_abbr, '[M+H]-sn1',
                                                    weight_dct['[M+H]-sn1'] * (10 - _rank_l_sn1) / 10)
                        else:
                            lipid_abbr_df.set_value(_i_abbr, '[M+H]-sn1', weight_dct['[M+H]-sn1'] * r_lyso1_i * 0.01)
                    if _sn2_abbr in lyso_ident_lst:
                        _rank_l_sn2 = lyso_ident_lst.index(_sn2_abbr)
                        r_lyso2_i = 100 * lyso_i_lst[_rank_l_sn2] / ms2_max_i
                        lipid_abbr_df.set_value(_i_abbr, 'i_[M+H]-sn2', r_lyso2_i)
                        if rank_mode is True:
                            lipid_abbr_df.set_value(_i_abbr, '[M+H]-sn2',
                                                    weight_dct['[M+H]-sn2'] * (10 - _rank_l_sn2) / 10)
                        else:
                            lipid_abbr_df.set_value(_i_abbr, '[M+H]-sn2', weight_dct['[M+H]-sn2'] * r_lyso2_i * 0.01)
                    if _sn3_abbr in lyso_ident_lst:
                        _rank_l_sn3 = lyso_ident_lst.index(_sn3_abbr)
                        r_lyso3_i = 100 * lyso_i_lst[_rank_l_sn3] / ms2_max_i
                        lipid_abbr_df.set_value(_i_abbr, 'i_[M+H]-sn3', r_lyso3_i)
                        if rank_mode is True:
                            lipid_abbr_df.set_value(_i_abbr, '[M+H]-sn3',
                                                    weight_dct['[M+H]-sn3'] * (10 - _rank_l_sn3) / 10)
                        else:
                            lipid_abbr_df.set_value(_i_abbr, '[M+H]-sn3', weight_dct['[M+H]-sn3'] * r_lyso3_i * 0.01)
                    if _sn1_abbr in lyso_w_ident_lst:
                        _rank_lw_sn1 = lyso_w_ident_lst.index(_sn1_abbr)
                        r_lyso_w1_i = 100 * lyso_w_i_lst[_rank_lw_sn1] / ms2_max_i
                        lipid_abbr_df.set_value(_i_abbr, 'i_[M+H]-sn1-H2O', r_lyso_w1_i)
                        #############################################
                        #
                        # Can go in a function to avoid the double check
                        #
                        ##############################################
                        if rank_mode is True:
                            if sodiumFlag == 1:
                                lipid_abbr_df.set_value(_i_abbr, '[M+Na]-sn1',
                                                        weight_dct['[M+Na]-sn1'] * (10 - _rank_lw_sn1) / 10)
                            else:
                                lipid_abbr_df.set_value(_i_abbr, '[M+H]-sn1-H2O',
                                                    weight_dct['[M+H]-sn1-H2O'] * (10 - _rank_lw_sn1) / 10)
                        else:
                            if sodiumFlag == 1:
                                lipid_abbr_df.set_value(_i_abbr, '[M+Na]-sn1',
                                                        weight_dct['[M+Na]-sn1'] * r_lyso_w1_i * 0.01)
                            else:
                                lipid_abbr_df.set_value(_i_abbr, '[M+H]-sn1-H2O',
                                                    weight_dct['[M+H]-sn1-H2O'] * r_lyso_w1_i * 0.01)
                    if _sn2_abbr in lyso_w_ident_lst:
                        _rank_lw_sn2 = lyso_w_ident_lst.index(_sn2_abbr)
                        r_lyso_w2_i = 100 * lyso_w_i_lst[_rank_lw_sn2] / ms2_max_i
                        lipid_abbr_df.set_value(_i_abbr, 'i_[M+H]-sn2-H2O', r_lyso_w2_i)
                        if rank_mode is True:
                            if sodiumFlag == 1:
                                lipid_abbr_df.set_value(_i_abbr, '[M+Na]-sn2',
                                                        weight_dct['[M+Na]-sn2'] * (10 - _rank_lw_sn2) / 10)
                            else:
                                lipid_abbr_df.set_value(_i_abbr, '[M+H]-sn2-H2O',
                                                    weight_dct['[M+H]-sn2-H2O'] * (10 - _rank_lw_sn2) / 10)
                        else:
                            if sodiumFlag ==1:
                                lipid_abbr_df.set_value(_i_abbr, '[M+Na]-sn2',
                                                        weight_dct['[M+Na]-sn2'] * r_lyso_w2_i * 0.01)
                            else:
                                lipid_abbr_df.set_value(_i_abbr, '[M+H]-sn2-H2O',
                                                    weight_dct['[M+H]-sn2-H2O'] * r_lyso_w2_i * 0.01)
                    if _sn3_abbr in lyso_w_ident_lst:
                        _rank_lw_sn3 = lyso_w_ident_lst.index(_sn3_abbr)
                        r_lyso_w3_i = 100 * lyso_w_i_lst[_rank_lw_sn3] / ms2_max_i
                        lipid_abbr_df.set_value(_i_abbr, 'i_[M+H]-sn3-H2O', r_lyso_w3_i)
                        if rank_mode is True:
                            if sodiumFlag == 1:
                                lipid_abbr_df.set_value(_i_abbr, '[M+Na]-sn3',
                                                        weight_dct['[M+Na]-sn3'] * (10 - _rank_lw_sn3) / 10)
                            else:
                                lipid_abbr_df.set_value(_i_abbr, '[M+H]-sn3-H2O',
                                                    weight_dct['[M+H]-sn3-H2O'] * (10 - _rank_lw_sn3) / 10)
                        else:
                            if sodiumFlag == 1:
                                lipid_abbr_df.set_value(_i_abbr, '[M+Na]-sn3',
                                                        weight_dct['[M+Na]-sn3'] * r_lyso_w3_i * 0.01)
                            else:
                                lipid_abbr_df.set_value(_i_abbr, '[M+H]-sn3-H2O',
                                                    weight_dct['[M+H]-sn3-H2O'] * r_lyso_w3_i * 0.01)

                    geo= [_sn1_abbr, _sn2_abbr]
                    if geo in mg_w_ident_lst2:
                        print ('sn1,sn2')
                        exit()
                        _rank_mgw_sn1_sn2 = mg_w_ident_lst2.index([_sn1_abbr, _sn2_abbr])
                        r_mg_w_i = 100 * mg_w_i_lst[_rank_mgw_sn1_sn2] / ms2_max_i
                        lipid_abbr_df.set_value(_i_abbr, 'i_[M+H]-(sn1+sn2)-H2O', r_mg_w_i)
                        if rank_mode is True:
                            lipid_abbr_df.set_value(_i_abbr, '[M+H]-(sn1+sn2)-H2O',
                                                    weight_dct['[M+H]-(sn1+sn2)-H2O'] * (10 - _rank_mgw_sn1_sn2) / 10)
                        else:
                            lipid_abbr_df.set_value(_i_abbr, '[M+H]-(sn1+sn2)-H2O',
                                                    weight_dct['[M+H]-(sn1+sn2)-H2O'] * r_mg_w_i * 0.01)

                    geo2 = [_sn2_abbr, _sn1_abbr]
                    if geo2 in mg_w_ident_lst2:
                        print ('sn2, sn1')
                        exit()
                        _rank_mgw_sn1_sn2 = mg_w_ident_lst2.index([_sn2_abbr, _sn1_abbr])
                        r_mg_w_i = 100 * mg_w_i_lst[_rank_mgw_sn1_sn2] / ms2_max_i
                        lipid_abbr_df.set_value(_i_abbr, 'i_[M+H]-(sn1+sn2)-H2O', r_mg_w_i)
                        if rank_mode is True:
                            lipid_abbr_df.set_value(_i_abbr, '[M+H]-(sn1+sn2)-H2O',
                                                    weight_dct['[M+H]-(sn1+sn2)-H2O'] * (
                                                    10 - _rank_mgw_sn1_sn2) / 10)
                        else:
                            lipid_abbr_df.set_value(_i_abbr, '[M+H]-(sn1+sn2)-H2O',
                                                    weight_dct['[M+H]-(sn1+sn2)-H2O'] * r_mg_w_i * 0.01)

                    geo3= [_sn2_abbr, _sn3_abbr]
                    if geo3 in mg_w_ident_lst2:
                        print ('sn2,sn3')
                        exit()
                        _rank_mgw_sn2_sn3 = mg_w_ident_lst2.index([_sn2_abbr, _sn3_abbr])
                        r_mg_w_i = 100 * mg_w_i_lst[_rank_mgw_sn2_sn3] / ms2_max_i
                        lipid_abbr_df.set_value(_i_abbr, 'i_[M+H]-(sn2+sn3)-H2O', r_mg_w_i)
                        if rank_mode is True:
                            lipid_abbr_df.set_value(_i_abbr, '[M+H]-(sn2+sn3)-H2O',
                                                    weight_dct['[M+H]-(sn2+sn3)-H2O'] * (
                                                        10 - _rank_mgw_sn2_sn3) / 10)
                        else:
                            lipid_abbr_df.set_value(_i_abbr, '[M+H]-(sn2+sn3)-H2O',
                                                    weight_dct['[M+H]-(sn2+sn3)-H2O'] * r_mg_w_i * 0.01)

                    geo4= [_sn3_abbr, _sn2_abbr]
                    if geo4 in mg_w_ident_lst2:
                        print ('sn3,sn2')
                        exit()
                        _rank_mgw_sn2_sn3 = mg_w_ident_lst2.index([_sn3_abbr, _sn2_abbr])
                        r_mg_w_i = 100 * mg_w_i_lst[_rank_mgw_sn2_sn3] / ms2_max_i
                        lipid_abbr_df.set_value(_i_abbr, 'i_[M+H]-(sn2+sn3)-H2O', r_mg_w_i)
                        if rank_mode is True:
                            lipid_abbr_df.set_value(_i_abbr, '[M+H]-(sn2+sn3)-H2O',
                                                    weight_dct['[M+H]-(sn2+sn3)-H2O'] * (
                                                        10 - _rank_mgw_sn2_sn3) / 10)
                        else:
                            lipid_abbr_df.set_value(_i_abbr, '[M+H]-(sn2+sn3)-H2O',
                                                    weight_dct['[M+H]-(sn2+sn3)-H2O'] * r_mg_w_i * 0.01)

                    geo5 =[_sn1_abbr, _sn3_abbr]
                    if geo5 in mg_w_ident_lst2:
                        print ('sn1, sn3')
                        exit()
                        _rank_mgw_sn1_sn3 = mg_w_ident_lst2.index([_sn1_abbr, _sn3_abbr])
                        r_mg_w_i = 100 * mg_w_i_lst[_rank_mgw_sn1_sn3] / ms2_max_i
                        lipid_abbr_df.set_value(_i_abbr, 'i_[M+H]-(sn1+sn3)-H2O', r_mg_w_i)
                        if rank_mode is True:
                            lipid_abbr_df.set_value(_i_abbr, '[M+H]-(sn1+sn3)-H2O',
                                                    weight_dct['[M+H]-(sn1+sn3)-H2O'] * (
                                                    10 - _rank_mgw_sn1_sn3) / 10)
                        else:
                            lipid_abbr_df.set_value(_i_abbr, '[M+H]-(sn1+sn3)-H2O',
                                                    weight_dct['[M+H]-(sn1+sn3)-H2O'] * r_mg_w_i * 0.01)

                    geo6= [_sn3_abbr, _sn1_abbr]
                    if geo6 in mg_w_ident_lst2:
                        print ('sn3, sn1')
                        _rank_mgw_sn1_sn3 = mg_w_ident_lst2.index([_sn3_abbr, _sn1_abbr])
                        r_mg_w_i = 100 * mg_w_i_lst[_rank_mgw_sn1_sn3] / ms2_max_i
                        lipid_abbr_df.set_value(_i_abbr, 'i_[M+H]-(sn1+sn3)-H2O', r_mg_w_i)
                        if rank_mode is True:
                            lipid_abbr_df.set_value(_i_abbr, '[M+H]-(sn1+sn3)-H2O',
                                                    weight_dct['[M+H]-(sn1+sn3)-H2O'] * (
                                                        10 - _rank_mgw_sn1_sn3) / 10)
                        else:
                            lipid_abbr_df.set_value(_i_abbr, '[M+H]-(sn1+sn3)-H2O',
                                                    weight_dct['[M+H]-(sn1+sn3)-H2O'] * r_mg_w_i * 0.01)

                                    # print lipid_abbr_df
                print ('Hellooooooooo')

            else:

                for _i_abbr, _abbr_se in lipid_abbr_df.iterrows():
                    # _pl_abbr = _abbr_se['Lipid_abbr']
                    _sn1_abbr = _abbr_se['sn1_abbr']
                    _sn2_abbr = _abbr_se['sn2_abbr']

                    if _sn1_abbr in fa_ident_lst:
                        _rank_sn1 = fa_ident_lst.index(_sn1_abbr)
                        r_sn1_i = 100 * fa_i_lst[_rank_sn1] / ms2_max_i
                        lipid_abbr_df.set_value(_i_abbr, 'i_sn1', r_sn1_i)
                        if rank_mode is True:
                            lipid_abbr_df.set_value(_i_abbr, 'sn1', weight_dct['sn1'] * (10 - _rank_sn1) / 10)
                        else:
                            lipid_abbr_df.set_value(_i_abbr, 'sn1', weight_dct['sn1'] * r_sn1_i * 0.01)
                    if _sn2_abbr in fa_ident_lst:
                        _rank_sn2 = fa_ident_lst.index(_sn2_abbr)
                        r_sn2_i = 100 * fa_i_lst[_rank_sn2] / ms2_max_i
                        lipid_abbr_df.set_value(_i_abbr, 'i_sn2', r_sn2_i)
                        if rank_mode is True:
                            lipid_abbr_df.set_value(_i_abbr, 'sn2', weight_dct['sn2'] * (10 - _rank_sn2) / 10)
                        else:
                            lipid_abbr_df.set_value(_i_abbr, 'sn2', weight_dct['sn2'] * r_sn2_i * 0.01)
                    if _sn1_abbr in lyso_ident_lst:
                        _rank_l_sn1 = lyso_ident_lst.index(_sn1_abbr)
                        r_lyso1_i = 100 * lyso_i_lst[_rank_l_sn1] / ms2_max_i
                        lipid_abbr_df.set_value(_i_abbr, 'i_[M-H]-sn1', r_lyso1_i)
                        if rank_mode is True:
                            lipid_abbr_df.set_value(_i_abbr, '[M-H]-sn1',
                                                    weight_dct['[M-H]-sn1'] * (10 - _rank_l_sn1) / 10)
                        else:
                            lipid_abbr_df.set_value(_i_abbr, '[M-H]-sn1', weight_dct['[M-H]-sn1'] * r_lyso1_i * 0.01)
                    if _sn2_abbr in lyso_ident_lst:
                        _rank_l_sn2 = lyso_ident_lst.index(_sn2_abbr)
                        r_lyso2_i = 100 * lyso_i_lst[_rank_l_sn2] / ms2_max_i
                        lipid_abbr_df.set_value(_i_abbr, 'i_[M-H]-sn2', r_lyso2_i)
                        if rank_mode is True:
                            lipid_abbr_df.set_value(_i_abbr, '[M-H]-sn2',
                                                    weight_dct['[M-H]-sn2'] * (10 - _rank_l_sn2) / 10)
                        else:
                            lipid_abbr_df.set_value(_i_abbr, '[M-H]-sn2', weight_dct['[M-H]-sn2'] * r_lyso2_i * 0.01)
                    if _sn1_abbr in lyso_w_ident_lst:
                        _rank_lw_sn1 = lyso_w_ident_lst.index(_sn1_abbr)
                        r_lyso_w1_i = 100 * lyso_w_i_lst[_rank_lw_sn1] / ms2_max_i
                        lipid_abbr_df.set_value(_i_abbr, 'i_[M-H]-sn1-H2O', r_lyso_w1_i)
                        if rank_mode is True:
                            lipid_abbr_df.set_value(_i_abbr, '[M-H]-sn1-H2O',
                                                    weight_dct['[M-H]-sn1-H2O'] * (10 - _rank_lw_sn1) / 10)
                        else:
                            lipid_abbr_df.set_value(_i_abbr, '[M-H]-sn1-H2O',
                                                    weight_dct['[M-H]-sn1-H2O'] * r_lyso_w1_i * 0.01)
                    if _sn2_abbr in lyso_w_ident_lst:
                        _rank_lw_sn2 = lyso_w_ident_lst.index(_sn2_abbr)
                        r_lyso_w2_i = 100 * lyso_w_i_lst[_rank_lw_sn2] / ms2_max_i
                        lipid_abbr_df.set_value(_i_abbr, 'i_[M-H]-sn2-H2O', r_lyso_w2_i)
                        if rank_mode is True:
                            lipid_abbr_df.set_value(_i_abbr, '[M-H]-sn2-H2O',
                                                    weight_dct['[M-H]-sn2-H2O'] * (10 - _rank_lw_sn2) / 10)
                        else:
                            lipid_abbr_df.set_value(_i_abbr, '[M-H]-sn2-H2O',
                                                    weight_dct['[M-H]-sn2-H2O'] * r_lyso_w2_i * 0.01)

            lipid_abbr_df['Score'] = lipid_abbr_df[weight_type_lst].sum(axis=1, numeric_only=True)
            match_reporter = 1


        else:
            print('!!!!!! NO FA identified =====>--> Skip >>> >>>')
        print ('matxh')
        match_info_dct = {'MATCH_INFO': match_reporter, 'SCORE_INFO': lipid_abbr_df, 'FA_INFO': fa_ident_df,
                          'LYSO_INFO': lyso_ident_df, 'LYSO_W_INFO': lyso_w_ident_df, 'MG_W_INFO' : mg_w_ident_df}
        print (match_info_dct)
        return match_info_dct

    def get_specific_peaks(self, mz_lib, ms2_df, ms2_precision=50e-6, ms2_threshold=10,
                           ms2_hginfo_threshold=0.02, vendor='waters'):

        _target_frag_df = pd.DataFrame()
        _target_nl_df = pd.DataFrame()
        _other_frag_df = pd.DataFrame()
        _other_nl_df = pd.DataFrame()

        ms2_max_i = ms2_df['i'].max()
        ms2_hginfo_abs_i = ms2_max_i * ms2_hginfo_threshold
        ms2_threshold = max(ms2_threshold, ms2_hginfo_abs_i)

        for _i, _frag_se in self.target_frag_df.iterrows():

            _frag_mz = _frag_se['EXACTMASS']
            _frag_class = _frag_se['CLASS']
            _frag_label = _frag_se['LABEL']
            if vendor == 'thermo':
                if mz_lib < 450.0:
                    seg_shift = 0.0
                elif 450.0 <= mz_lib < 600.0:
                    seg_shift = 0.1
                elif 600.0 <= mz_lib < 750.0:
                    seg_shift = 0.2
                elif 750.0 <= mz_lib < 900.0:
                    seg_shift = 0.3
                elif 900.0 <= mz_lib < 1200.0:
                    seg_shift = 0.4
                elif 1200.0 <= mz_lib:
                    seg_shift = 0.5
                else:
                    seg_shift = 0.0
                _delta = _frag_mz * ms2_precision + seg_shift
                _frag_mz_low = _frag_mz - _delta
                _frag_mz_high = _frag_mz + _delta
            else:
                _frag_mz_low = _frag_mz - _frag_mz * ms2_precision
                _frag_mz_high = _frag_mz + _frag_mz * ms2_precision
            _frag_mz_query_code = '%f <= mz <= %f and i > %f' % (_frag_mz_low, _frag_mz_high, ms2_threshold)

            _frag_df = ms2_df.query(_frag_mz_query_code)

            if _frag_df.shape[0] > 0:
                _frag_df = _frag_df.sort_values(by='i', ascending=False)
                _frag_df.loc[:, 'CLASS'] = _frag_class
                _frag_df.loc[:, 'LABEL'] = _frag_label
                _frag_df.loc[:, _frag_label] = 100 * _frag_df['i'] / ms2_max_i
                _target_frag_df = _target_frag_df.append(_frag_df.head(1))

        for _i, _frag_se in self.other_frag_df.iterrows():

            _frag_mz = _frag_se['EXACTMASS']
            _frag_class = _frag_se['CLASS']
            _frag_label = _frag_se['LABEL']

            if vendor == 'thermo':
                if _frag_mz < 450.0:
                    seg_shift = 0.0
                elif 450.0 <= mz_lib < 600.0:
                    seg_shift = 0.1
                elif 600.0 <= mz_lib < 750.0:
                    seg_shift = 0.2
                elif 750.0 <= mz_lib < 900.0:
                    seg_shift = 0.3
                elif 900.0 <= mz_lib < 1200.0:
                    seg_shift = 0.4
                elif 1200.0 <= mz_lib:
                    seg_shift = 0.5
                else:
                    seg_shift = 0.0

                _delta = _frag_mz * ms2_precision + seg_shift
                _frag_mz_low = _frag_mz - _delta
                _frag_mz_high = _frag_mz + _delta

            else:
                _frag_mz_low = _frag_mz - _frag_mz * ms2_precision
                _frag_mz_high = _frag_mz + _frag_mz * ms2_precision
            _frag_mz_query_code = '%f <= mz <= %f and i > %f' % (_frag_mz_low, _frag_mz_high, ms2_threshold)

            _frag_df = ms2_df.query(_frag_mz_query_code)

            if _frag_df.shape[0] > 0:
                _frag_df = _frag_df.sort_values(by='i', ascending=False)
                _frag_df.loc[:, 'CLASS'] = _frag_class
                _frag_df.loc[:, 'LABEL'] = _frag_label
                _frag_df.loc[:, _frag_label] = 100 * _frag_df['i'] / ms2_max_i
                _other_frag_df = _other_frag_df.append(_frag_df.head(1))

        for _i, _nl_se in self.target_nl_df.iterrows():

            _nl_mz = _nl_se['EXACTMASS']
            _nl_class = _nl_se['CLASS']
            _nl_label = _nl_se['LABEL']

            if vendor == 'thermo':
                if mz_lib < 450.0:
                    seg_shift = 0.0
                elif 450.0 <= mz_lib < 600.0:
                    seg_shift = 0.1
                elif 600.0 <= mz_lib < 750.0:
                    seg_shift = 0.2
                elif 750.0 <= mz_lib < 900.0:
                    seg_shift = 0.3
                elif 900.0 <= mz_lib < 1200.0:
                    seg_shift = 0.4
                elif 1200.0 <= mz_lib:
                    seg_shift = 0.5
                else:
                    seg_shift = 0.0

                _delta = (mz_lib - _nl_mz) * ms2_precision * ms2_precision + seg_shift
                _nl_mz_low = mz_lib - _nl_mz - _delta
                _nl_mz_high = mz_lib - _nl_mz + _delta

            else:
                _nl_mz_low = mz_lib - _nl_mz - _nl_mz * ms2_precision
                _nl_mz_high = mz_lib - _nl_mz + _nl_mz * ms2_precision

            _nl_mz_query_code = '%f <= mz <= %f and i > %f' % (_nl_mz_low, _nl_mz_high, ms2_threshold)

            _nl_df = ms2_df.query(_nl_mz_query_code)

            if _nl_df.shape[0] > 0:
                _nl_df = _nl_df.sort_values(by='i', ascending=False)
                _nl_df.loc[:, 'CLASS'] = _nl_class
                _nl_df.loc[:, 'LABEL'] = _nl_label
                _nl_df.loc[:, _nl_label] = 100 * _nl_df['i'] / ms2_max_i
                _target_nl_df = _target_nl_df.append(_nl_df.head(1))

        for _i, _nl_se in self.other_nl_df.iterrows():

            _nl_mz = _nl_se['EXACTMASS']
            _nl_class = _nl_se['CLASS']
            _nl_label = _nl_se['LABEL']

            if vendor == 'thermo':
                if mz_lib < 450.0:
                    seg_shift = 0.0
                elif 450.0 <= mz_lib < 600.0:
                    seg_shift = 0.1
                elif 600.0 <= mz_lib < 750.0:
                    seg_shift = 0.2
                elif 750.0 <= mz_lib < 900.0:
                    seg_shift = 0.3
                elif 900.0 <= mz_lib < 1200.0:
                    seg_shift = 0.4
                elif 1200.0 <= mz_lib:
                    seg_shift = 0.5
                else:
                    seg_shift = 0.0

                _delta = (mz_lib - _nl_mz) * ms2_precision + seg_shift
                _nl_mz_low = mz_lib - _nl_mz - _delta
                _nl_mz_high = mz_lib - _nl_mz + _delta

            else:
                _nl_mz_low = mz_lib - _nl_mz - _nl_mz * ms2_precision
                _nl_mz_high = mz_lib - _nl_mz + _nl_mz * ms2_precision

            _nl_mz_query_code = '%f <= mz <= %f and i > %f' % (_nl_mz_low, _nl_mz_high, ms2_threshold)

            _nl_df = ms2_df.query(_nl_mz_query_code)

            if _nl_df.shape[0] > 0:
                _nl_df = _nl_df.sort_values(by='i', ascending=False)
                _nl_df.loc[:, 'CLASS'] = _nl_class
                _nl_df.loc[:, 'LABEL'] = _nl_label
                _nl_df.loc[:, _nl_label] = 100 * _nl_df['i'] / ms2_max_i
                _other_nl_df = _other_nl_df.append(_nl_df.head(1))

        specific_ion_dct = {}
        if _target_frag_df.shape[0] > 0:
            specific_ion_dct['TARGET_FRAG'] = _target_frag_df
        if _target_nl_df.shape[0] > 0:
            specific_ion_dct['TARGET_NL'] = _target_nl_df
        if _other_frag_df.shape[0] > 0:
            specific_ion_dct['OTHER_FRAG'] = _other_frag_df
        if _other_nl_df.shape[0] > 0:
            specific_ion_dct['OTHER_NL'] = _other_nl_df

        return specific_ion_dct


if __name__ == '__main__':
    pass
