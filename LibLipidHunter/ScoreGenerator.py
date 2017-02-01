# -*- coding: utf-8 -*-
# Copyright 2015-2017 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de
# Changes: 17/1/2017 Georgia 14:40

from __future__ import division

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
        print (abbr)
        lipid_checker = re.compile(r'(\w{2,3})\((\w\-)*((\d{1,2})\:(\d{1,2}))\)')

        # Check PL Type
        _pl_typ = ''
        bulk_fa_typ = ''
        bulk_fa_linker = ''
        bulk_fa_c = 0
        bulk_fa_db = 0
        lyso_fa_linker_dct = {'sn1': '', 'sn2': ''}

######### Georgia Changes 17-1-2017
        pl_re_chk2 = lipid_checker.match(abbr)
        pl_typ_lst2 = pl_re_chk2.groups()
        if pl_typ_lst2[0] in ['PA', 'PE', 'PC', 'PS', 'PG', 'PI']:
            _pl_typ=pl_typ_lst2[0]
            if pl_typ_lst2[1] == 'O-':
                bulk_fa_linker = 'O-A-'
                lyso_fa_linker_dct = {'O':'', 'A':'O-'}
                #bulk_fa_c=pl_typ_lst2[3]
                #bulk_fa_db=pl_typ_lst2[4]
            elif pl_typ_lst2[1] == 'P-':
                bulk_fa_linker = 'P-A-'
                lyso_fa_linker_dct = {'P': '', 'A':'P-'}
                #bulk_fa_c = pl_typ_lst2[3]
                #bulk_fa_db = pl_typ_lst2[4]
            elif pl_typ_lst2[1]:
                bulk_fa_linker = 'A-A-'
                lyso_fa_linker_dct={'A':''}
            else:
                bulk_fa_linker = 'A-A-'
                lyso_fa_linker_dct = {'A': ''}
        else:
            if pl_typ_lst2[1] == 'O-':
                bulk_fa_linker = 'O'
            elif pl_typ_lst2[1] == 'P-':
                bulk_fa_linker = 'P'
            else:
                bulk_fa_linker = 'A'
            _pl_typ= pl_typ_lst2[0]
        bulk_fa_c = pl_typ_lst2[3]
        bulk_fa_db = pl_typ_lst2[4]
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
            fa_chk_df['M-sn'] = calc_pr_mz - fa_chk_df['NL']
            fa_chk_df['M-(sn-H2O)'] = calc_pr_mz - fa_chk_df['NL-H2O']
            fa_chk_df['Proposed_structures'] = ''

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
                            _frag_df.loc[:, 'Proposed_structures'] = 'FA %s [M-H]-' % _fa_abbr
                            fa_ident_df = fa_ident_df.append(_frag_df)
                        elif _frag_type == 'M-sn':
                            if _fa_link in lyso_fa_linker_dct.keys():
                                if bulk_fa_db - _fa_db >= 0:
                                    _fa_lyso_link = lyso_fa_linker_dct[_fa_link]
                                    _frag_df.loc[:, 'Proposed_structures'] = 'Lyso%s(%s%i:%i) [M-H2O-H]-' % (pl_typ,
                                                                                                     _fa_lyso_link,
                                                                                                     bulk_fa_c - _fa_c,
                                                                                                     bulk_fa_db - _fa_db
                                                                                                     )
                                    lyso_ident_df = lyso_ident_df.append(_frag_df)
                        elif _frag_type == 'M-(sn-H2O)':
                            if _fa_link in lyso_fa_linker_dct.keys():
                                if bulk_fa_db - _fa_db >= 0:
                                    _fa_lyso_link = lyso_fa_linker_dct[_fa_link]
                                    _frag_df.loc[:, 'Proposed_structures'] = 'Lyso%s(%s%i:%i) [M-H]-' % (pl_typ, _fa_lyso_link,
                                                                                                 bulk_fa_c - _fa_c,
                                                                                                 bulk_fa_db - _fa_db
                                                                                                 )
                                    lyso_w_ident_df = lyso_w_ident_df.append(_frag_df)

        elif lipid_type == 'GL' and charge_mode == 'POS':
            _fa_compination=[]
            ####################################################
            #   Georgia changes 18-19/1/2017
            #####################################################
            for _i, _fa_se in self.fa_def_df.iterrows():
                fa_chk_df = self.fa_def_df[['FA', 'Link', 'C', 'DB', 'mass', '[M+H]+']]
                fa_chk_df = fa_chk_df.rename(columns={'FA': 'FA2', 'Link': 'Link2', 'C':'C2', 'DB': 'DB2','[M+H]+': 'sn2', 'mass': 'NL2'})
                fa_chk_df['FA'] = _fa_se['FA']
                fa_chk_df['C'] = _fa_se['C']
                fa_chk_df['DB'] = _fa_se['DB']
                fa_chk_df['Link'] = _fa_se['Link']
                fa_chk_df['NL'] = _fa_se['mass']
                fa_chk_df['sn'] = _fa_se['[M+H]+']
                fa_chk_df['M-sn'] = calc_pr_mz - _fa_se['mass']
                fa_chk_df['M-sn-sn2'] = calc_pr_mz - _fa_se['mass'] - fa_chk_df['NL2']+18.010565
                fa_chk_df['Abbr'] = _fa_se['FA'] + '/' + fa_chk_df['FA2']
                fa_chk_df['Proposed_structures'] = ''

                for _frag_type in ['sn', 'M-sn', 'sn2', 'M-sn-sn2']:
                    if _frag_type in ['sn', 'M-sn']:
                        _frag_mz = fa_chk_df.loc[_i, _frag_type]
                        _frag_mz_low = _frag_mz - _frag_mz * ms2_precision
                        _frag_mz_high = _frag_mz + _frag_mz * ms2_precision
                        _frag_mz_query_code = '%f <= mz <= %f' % (_frag_mz_low, _frag_mz_high)
                        _frag_df = ms2_df.query(_frag_mz_query_code)
                        if _frag_df.shape[0] > 0:
                            _frag_df.loc[:, 'ppm'] = 1e6 * (_frag_df['mz'] - _frag_mz) / _frag_mz
                            _frag_df.loc[:, 'ppm_abs'] = _frag_df['ppm'].abs()
                            _frag_df.loc[:, 'FA'] = fa_chk_df.loc[_i, 'FA']
                            if _frag_df.shape[0] > 1:
                                _frag_i_df = _frag_df.sort_values(by='i', ascending=False).head(1)
                                _frag_ppm_df = _frag_df.sort_values(by='ppm_abs').head(1)
                                _frag_df = _frag_i_df
                                if _frag_ppm_df['i'].tolist() == _frag_i_df['i'].tolist():
                                    pass
                                else:
                                    _frag_df = _frag_df.append(_frag_ppm_df)
                            if _frag_type == 'sn':
                                _frag_df.loc[:, 'Proposed_structures'] = 'FA %s [M+H]+' % fa_chk_df.loc[_i,'FA']
                                fa_ident_df = fa_ident_df.append(_frag_df)
                            elif _frag_type == 'M-sn':
                                #####################################
                                # For now is not working
                                ################################
                                #if _fa_se['Link'] in lyso_fa_linker_dct.keys():
                                if bulk_fa_db - fa_chk_df.loc[_i, 'DB'] >= 0:
                                        #_fa_lyso_link = lyso_fa_linker_dct[_fa_se['Link']]
                                        _frag_df.loc[:, 'Proposed_structures'] = 'DG-%i:%i [M+H]+' % ( bulk_fa_c- fa_chk_df.loc[_i,'C'],int(bulk_fa_db - fa_chk_df.loc[_i,'DB']))
                                        _frag_df.loc[:, 'Short_name'] = 'DG-%i:%i' % (bulk_fa_c- fa_chk_df.loc[_i,'C'],int(bulk_fa_db - fa_chk_df.loc[_i,'DB']))
                                        lyso_ident_df = lyso_ident_df.append(_frag_df)
                    elif _frag_type in ['M-sn-sn2']:
                        for _i2, _fa_se2 in fa_chk_df.iterrows():
                            if sorted([_fa_se2['FA'], _fa_se2['FA2']]) not in _fa_compination:
                                _fa_compination.append(sorted([_fa_se2['FA'], _fa_se2['FA2']]))
                                _frag_mz = _fa_se2[_frag_type]
                                _frag_mz_low = _frag_mz - _frag_mz * ms2_precision
                                _frag_mz_high = _frag_mz + _frag_mz * ms2_precision
                                _frag_mz_query_code = '%f <= mz <= %f' % (_frag_mz_low, _frag_mz_high)
                                _frag_df = ms2_df.query(_frag_mz_query_code)
                                if _frag_df.shape[0] > 0:
                                    _frag_df.loc[:, 'ppm']=1e6*(_frag_df['mz']-_frag_mz)/_frag_mz
                                    _frag_df.loc[:, 'ppm_abs']= _frag_df['ppm'].abs()
                                    _frag_df.loc[:, 'FA'] = _fa_se2['Abbr']
                                    _frag_df.loc[:, 'First'] = _fa_se2['FA']
                                    _frag_df.loc[:, 'Second'] = _fa_se2['FA2']
                                    if _frag_df.shape[0] > 1:
                                        _frag_i_df = _frag_df.sort_values(by='i', ascending=False).head(1)
                                        _frag_ppm_df =_frag_df.sort_values(by='ppm_abs').head(1)
                                        _frag_df=_frag_i_df
                                        if _frag_ppm_df['i'].tolist() == _frag_i_df['i'].tolist():
                                            pass
                                        else:
                                            _frag_df = _frag_df.append(_frag_ppm_df)
                                    ###################################################
                                    #   For now this is deactivate
                                    ###################################################
                                    #if _fa_se2['Link'] in lyso_fa_linker_dct.keys():
                                    if bulk_fa_db - _fa_se2['DB'] - _fa_se2['DB2'] >= 0:
                                            #_fa_lyso_link = lyso_fa_linker_dct[_fa_se2['Link']]
                                            #####################################################################################
                                            #   Note: 29/1/2017:
                                            #   It doesn't give the link in the information
                                            #   Fix for later
                                            #################################################################################
                                            _frag_df.loc[:, 'Proposed_structures']='MG-%i:%i[M+H]+' % (bulk_fa_c - _fa_se2['C'] - _fa_se2['C2'], int(bulk_fa_db - _fa_se2['DB'] - _fa_se2['DB2']))
                                            _frag_df.loc[:, 'Short_name']='MG-%i:%i' % (bulk_fa_c- _fa_se2['C'] - _fa_se2['C2'], bulk_fa_db-_fa_se2['DB']-_fa_se2['DB2'])
                                            lyso_w_ident_df = lyso_w_ident_df.append(_frag_df)
        # format the output DataFrame
        if fa_ident_df.shape[0] > 0:
            fa_ident_df = fa_ident_df.query('i > %f' % ms2_threshold)
            proposed_str_lst={}
            if lipid_type == 'GL':
                fa_ident_df['Flag'] = 0
                for _ident, _info in fa_ident_df.iterrows():
                    if fa_ident_df.loc[_ident, 'Proposed_structures'] in proposed_str_lst.keys():
                        if _info['ppm_abs'] < proposed_str_lst[_info['Proposed_structures']][1]:
                            fa_ident_df.loc[_ident, 'Flag']=1
                            fa_ident_df.loc[proposed_str_lst[_info['Proposed_structures']][0], 'Flag']= 0
                            proposed_str_lst[_info['Proposed_structures']] = [_ident, _info['ppm_abs']]
                    else:
                        fa_ident_df.loc[_ident, 'Flag']=1
                        proposed_str_lst[_info['Proposed_structures']]=[_ident, _info['ppm_abs']]
                fa_ident_df = fa_ident_df.loc[fa_ident_df['Flag']==1][['Proposed_structures', 'FA', 'mz', 'i', 'ppm', 'ppm_abs', 'Flag']].reset_index(drop=True)
            else:
                fa_ident_df['Flag'] = 1
                fa_ident_df = fa_ident_df[
                    ['Proposed_structures', 'FA', 'mz', 'i', 'ppm', 'ppm_abs', 'Flag']].reset_index(drop=True)
            fa_ident_df = fa_ident_df.sort_values(by='i', ascending=False).head(10)

        if lyso_ident_df.shape[0] > 0:
            lyso_found_dct={}
            lyso_ident_df = lyso_ident_df.query('i > %f' % ms2_threshold)
            if lipid_type == 'GL':
                lyso_ident_df['Flag']=0
                for _ident, _inf_g in lyso_ident_df.iterrows():
                    if lyso_ident_df.loc[_ident, 'Short_name'] in lyso_found_dct.keys():
                        if _inf_g['ppm_abs'] < lyso_found_dct[_inf_g['Short_name']][1]:
                            lyso_ident_df.loc[_ident, 'Flag'] = 1
                            lyso_ident_df.loc[lyso_found_dct[_inf_g['Short_name']][0], 'Flag']=0
                            lyso_found_dct[_inf_g['Short_name']] = [_ident,_inf_g['ppm_abs']]
                    else:
                        lyso_ident_df.loc[_ident, 'Flag']=1
                        lyso_found_dct[_inf_g['Short_name']] = [_ident, _inf_g['ppm_abs']]
                lyso_ident_df = lyso_ident_df.loc[lyso_ident_df['Flag'] == 1][['Proposed_structures', 'FA', 'mz', 'i', 'ppm', 'ppm_abs', 'Flag']].reset_index(drop=True)
            else:
                lyso_ident_df['Flag'] = 1
                lyso_ident_df = lyso_ident_df.loc[lyso_ident_df['Flag'] == 1][
                    ['Proposed_structures', 'FA', 'mz', 'i', 'ppm', 'ppm_abs', 'Flag']].reset_index(drop=True)
            lyso_ident_df = lyso_ident_df.sort_values(by='i', ascending=False).head(5)

        if lyso_w_ident_df.shape[0] > 0:
            lyso_w_dct={}
            lyso_w_ident_df = lyso_w_ident_df.query('i > %f' % ms2_threshold).reset_index(drop=True)
            if lipid_type == 'GL':
                lyso_w_ident_df['Flag']=0
                for _ident, _inf in lyso_w_ident_df.iterrows():
                    if _inf['Short_name'] in lyso_w_dct.keys():
                        if _inf['ppm_abs'] < lyso_w_dct[_inf['Short_name']][1]:
                            lyso_w_ident_df.loc[_ident, 'Flag']=1
                            lyso_w_ident_df.loc[lyso_w_dct[_inf['Short_name']][0], 'Flag'] = 0
                            lyso_w_dct[_inf['Short_name']]=[_ident, _inf['ppm_abs']]
                    else:
                        lyso_w_ident_df.loc[_ident, 'Flag']=1
                        lyso_w_dct[_inf['Short_name']]=[_ident, _inf['ppm_abs']]
                lyso_w_ident_df = lyso_w_ident_df.loc[lyso_w_ident_df['Flag']==1][['Proposed_structures', 'FA', 'mz', 'i', 'ppm', 'ppm_abs', 'Flag', 'First', 'Second']].reset_index(drop=True)

            else:
                lyso_w_ident_df['Flag'] = 1
                lyso_w_ident_df = lyso_w_ident_df[
                    ['Proposed_structures', 'FA', 'mz', 'i', 'ppm', 'ppm_abs', 'Flag']].reset_index(
                    drop=True)
            lyso_w_ident_df = lyso_w_ident_df.sort_values(by='i', ascending=False).head(5)


        return fa_ident_df, lyso_ident_df, lyso_w_ident_df

    def get_structure(self, abbr):

        lipid_abbr_lst = []
        lipid_sn1_lst = []
        lipid_sn2_lst = []
        db_sn1_lst = []
        db_sn2_lst = []
        #abbr='TG(46:6)'
        print abbr

        lipid_info_dct = self.decode_abbr(abbr)
        pl_typ = lipid_info_dct['TYPE']
        bulk_fa_c = lipid_info_dct['C']
        bulk_fa_db = lipid_info_dct['DB']
        bulk_fa_linker = lipid_info_dct['LINK']

        if abbr[:2] in ['PE', 'PA', 'PC', 'PI', 'PS', 'PG']:
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
            lipid_abbr_df = pd.DataFrame(data={'Proposed_structures': lipid_abbr_lst, 'sn1_abbr': lipid_sn1_lst,
                                               'sn2_abbr': lipid_sn2_lst, 'sn1_DB': db_sn1_lst, 'sn2_DB': db_sn2_lst})

            lipid_abbr_df = lipid_abbr_df.query('sn1_DB <=sn2_DB')
            lipid_abbr_df = lipid_abbr_df[['Proposed_structures', 'sn1_abbr', 'sn2_abbr']]
            return lipid_abbr_df
        elif abbr[:2] in ['TG'] and bulk_fa_linker not in ['O', 'P']:
            #print bulk_fa_linker
            lipid_sn3_lst = []
            db_sn3_lst = []
            _fa_compination_3=[]
            for _i, _fa_se in self.fa_def_df.iterrows():
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

                        if int(_fa_db + _fa_db2) <= bulk_fa_db and int(_fa_c + _fa_c2) <= bulk_fa_c:
                            #This is for a later use when the bulk_fa_linker works
                            #if _fa_link == bulk_fa_linker[0:1]:
                            #if bulk_fa_linker == 'O' or bulk_fa_linker == 'P':
                             #   pass
                            #else:
                                #_rest_fa_link=bulk_fa_linker[2]
                                _rest_fa_c=bulk_fa_c - _fa_c - _fa_c2
                                _rest_fa_db=bulk_fa_db - _fa_db - _fa_db2
                                #_rest_fa_df=self.fa_def_df.query('Link == "%s" and C == %i and DB == %i' %(_rest_fa_link, _rest_fa_c, _rest_fa_db))
                                _rest_fa_df=self.fa_def_df.query('C == %i and DB == %i' %(_rest_fa_c, _rest_fa_db))


                                if _rest_fa_df.shape[0] == 1:
                                    _rest_fa_abbr = _rest_fa_df['FA'].tolist()[0]
                                    if sorted([_fa_abbr, _fa_abbr2, _rest_fa_abbr]) not in _fa_compination_3:
                                        _fa_compination_3.append(sorted([_fa_abbr, _fa_abbr2, _rest_fa_abbr]))
                                        lipid_abbr = '%s (%i:%i_%i:%i_%i:%i)' % (abbr[:2], _fa_se['C'], _fa_se['DB'], _fa_se2['C'], _fa_se2['DB'], _rest_fa_c, _rest_fa_db)
                                        lipid_abbr_lst.append(lipid_abbr)
                                        lipid_sn1_lst.append(_fa_abbr)
                                        lipid_sn2_lst.append(_fa_abbr2)
                                        lipid_sn3_lst.append(_rest_fa_abbr)
                                        db_sn1_lst.append(_fa_db)
                                        db_sn2_lst.append(_fa_db2)
                                        db_sn3_lst.append(_rest_fa_db)
            lipid_abbr_df = pd.DataFrame(data={'Proposed_structures': lipid_abbr_lst, 'sn1_abbr': lipid_sn1_lst,
                                                   'sn2_abbr': lipid_sn2_lst, 'sn3_abbr': lipid_sn3_lst,
                                                   'sn1_DB': db_sn1_lst,
                                                   'sn2_DB': db_sn2_lst, 'sn3_DB': db_sn3_lst})
            lipid_abbr_df = lipid_abbr_df[['Proposed_structures', 'sn1_abbr', 'sn2_abbr', 'sn3_abbr']]
            return lipid_abbr_df
        else:
                lipid_sn3_lst = []
                db_sn3_lst = []
                #_fa_compination_3 = []

                lipid_abbr_df = pd.DataFrame(data= {'Proposed_structures': lipid_abbr_lst, 'sn1_abbr': lipid_sn1_lst,
                                                'sn2_abbr':lipid_sn2_lst, 'sn3_abbr': lipid_sn3_lst, 'sn1_DB': db_sn1_lst,
                                                'sn2_DB': db_sn2_lst, 'sn3_DB': db_sn3_lst})
                lipid_abbr_df=lipid_abbr_df[['Proposed_structures', 'sn1_abbr', 'sn2_abbr', 'sn3_abbr']]

                return lipid_abbr_df

    def get_match(self, abbr, charge_type, mz_lib, ms2_df, ms2_precision=500e-6, ms2_threshold=100):

        match_reporter = 0
        ms2_max_i = ms2_df['i'].max()
        #abbr = 'TG(46:6)'
        fa_ident_df, lyso_ident_df, lyso_w_ident_df = self.get_fa_search(abbr, charge_type, mz_lib, ms2_df,
                                                                         ms2_precision=ms2_precision,
                                                                         ms2_threshold=ms2_threshold
                                                                         )
        lipid_abbr_df = self.get_structure(abbr)

        if abbr[:2] in ['TG']:
            weight_type_lst = ['sn1', 'sn2', 'sn3', 'M-sn1', 'M-sn2', 'M-sn3', 'M-(sn1+sn2)', 'M-(sn2+sn3)', 'M-(sn1+sn3)']
        else:
            weight_type_lst = ['sn1', 'sn2', 'M-sn1', 'M-sn2', 'M-(sn1-H2O)', 'M-(sn2-H2O)']
        weight_dct = {}
        for _type in weight_type_lst:
            lipid_abbr_df[_type] = 0
        if fa_ident_df.shape[0] > 0 or lyso_ident_df.shape[0] > 0  or lyso_w_ident_df.shape[0] > 0:
            combine_all_lst=pd.DataFrame()
            try:
                fa_ident_df['Type']='FA'
                fa_ident_lst = fa_ident_df.loc[fa_ident_df['Flag']== 1]['FA'].tolist()
                fa_i_lst = fa_ident_df.loc[fa_ident_df['Flag']== 1]['i'].tolist()
                combine_all_lst=combine_all_lst.append(fa_ident_df.loc[fa_ident_df['Flag']==1])
            except KeyError:
                fa_ident_lst = []
                fa_i_lst = []

            try:
                lyso_ident_df['Type']='Lyso'
                lyso_ident_lst = lyso_ident_df.loc[lyso_ident_df['Flag'] == 1]['FA'].tolist()
                lyso_i_lst = lyso_ident_df.loc[lyso_ident_df['Flag'] == 1]['i'].tolist()
                combine_all_lst=combine_all_lst.append(lyso_ident_df)
            except KeyError:
                lyso_ident_lst = []
                lyso_i_lst = []

            try:
                lyso_w_ident_df['Type']='LysoW'
                lyso_w_ident_lst = lyso_w_ident_df.loc[lyso_w_ident_df['Flag'] == 1]['FA'].tolist()
                lyso_w_i_lst = lyso_w_ident_df.loc[lyso_w_ident_df['Flag'] == 1]['i'].tolist()
                combine_all_lst = combine_all_lst.append(lyso_w_ident_df[['Proposed_structures', 'FA', 'mz', 'i', 'ppm', 'ppm_abs', 'Flag', 'Type']])
            except KeyError:
                lyso_w_ident_lst = []
                lyso_w_i_lst = []
            #print combine_all_lst
            self.weight_df['mz'] = 0.0
            for _i, _weight_se in self.weight_df.iterrows():
                _type = _weight_se['Type']
                _weight = _weight_se['Weight']
                weight_dct[_type] = _weight
            if abbr[:2] in ['TG']:
                for _i_abbr, _abbr_se in lipid_abbr_df.iterrows():
                    _sn1_abbr = _abbr_se['sn1_abbr']
                    _sn2_abbr = _abbr_se['sn2_abbr']
                    _sn3_abbr = _abbr_se['sn3_abbr']
                    sn2 =_sn1_abbr+'/'+_sn2_abbr
                    sn1 = _sn1_abbr+'/'+_sn3_abbr
                    sn3= _sn2_abbr+'/'+_sn3_abbr
                    sn4 = _sn2_abbr+'/'+_sn1_abbr
                    sn5 = _sn3_abbr+'/'+_sn1_abbr
                    sn6 = _sn3_abbr+'/'+_sn2_abbr
                    _ident_group=combine_all_lst.groupby(['mz'])
                    for name, group in _ident_group:
                        if group.shape[0]:
                            _ident_group_FA = group.sort_values(by='ppm_abs', ascending=True).head(1).reset_index(drop=True)
                            if _ident_group_FA.loc[0, 'FA'] in [_sn2_abbr]:
                                if _ident_group_FA.loc[0, 'FA'] in fa_ident_lst:
                                    _rank_sn1 = fa_ident_lst.index(_ident_group_FA.loc[0, 'FA'])
                                    lipid_abbr_df.set_value(_i_abbr, 'i_sn1',
                                                            100 * fa_i_lst[_rank_sn1] / ms2_max_i)
                                    lipid_abbr_df.set_value(_i_abbr, 'sn1',
                                                            weight_dct['sn1'] * (10 - _rank_sn1) / 10)
                                if _ident_group_FA.loc[0, 'FA'] in lyso_ident_lst:
                                    _rank_l_sn1 = lyso_ident_lst.index(_ident_group_FA.loc[0, 'FA'])
                                    lipid_abbr_df.set_value(_i_abbr, 'i_M-sn1',
                                                            100 * lyso_i_lst[_rank_l_sn1] / ms2_max_i)
                                    lipid_abbr_df.set_value(_i_abbr, 'M-sn1',
                                                            weight_dct['M-sn1'] * (10 - _rank_l_sn1) / 10)
                                if _ident_group_FA.loc[0, 'FA'] in [_sn1_abbr]:
                                        if _ident_group_FA.loc[0, 'FA'] in fa_ident_lst:
                                            _rank_sn1 = fa_ident_lst.index(_ident_group_FA.loc[0, 'FA'])
                                            lipid_abbr_df.set_value(_i_abbr, 'i_sn2',
                                                                        100 * fa_i_lst[_rank_sn1] / ms2_max_i)
                                            lipid_abbr_df.set_value(_i_abbr, 'sn2',
                                                                        weight_dct['sn2'] * (10 - _rank_sn1) / 10)
                                        if _ident_group_FA.loc[0, 'FA'] in lyso_ident_lst:
                                            _rank_l_sn1 = lyso_ident_lst.index(_ident_group_FA.loc[0, 'FA'])
                                            lipid_abbr_df.set_value(_i_abbr, 'i_M-sn2',
                                                                        100 * lyso_i_lst[_rank_l_sn1] / ms2_max_i)
                                            lipid_abbr_df.set_value(_i_abbr, 'M-sn2',
                                                                        weight_dct['M-sn2'] * (10 - _rank_l_sn1) / 10)
                                if _ident_group_FA.loc[0, 'FA'] in [_sn3_abbr]:
                                        if _ident_group_FA.loc[0, 'FA'] in fa_ident_lst:
                                            _rank_sn1 = fa_ident_lst.index(_ident_group_FA.loc[0, 'FA'])
                                            lipid_abbr_df.set_value(_i_abbr, 'i_sn3',
                                                                        100 * fa_i_lst[_rank_sn1] / ms2_max_i)
                                            lipid_abbr_df.set_value(_i_abbr, 'sn3',
                                                                        weight_dct['sn3'] * (10 - _rank_sn1) / 10)
                                        if _ident_group_FA.loc[0, 'FA'] in lyso_ident_lst:
                                            _rank_l_sn1 = lyso_ident_lst.index(_ident_group_FA.loc[0, 'FA'])
                                            lipid_abbr_df.set_value(_i_abbr, 'i_M-sn3',
                                                                        100 * lyso_i_lst[_rank_l_sn1] / ms2_max_i)
                                            lipid_abbr_df.set_value(_i_abbr, 'M-sn3',
                                                                        weight_dct['M-sn3'] * (10 - _rank_l_sn1) / 10)

                                #_ident_group_FA=group[group['FA'].isin([_sn2_abbr, _sn3_abbr, _sn1_abbr])]
                                #_ident_group_FA=_ident_group_FA.append(group[group['FA'].isin([sn1, sn2, sn3, sn4, sn5, sn6])])
                                #if _ident_group_FA.shape[0]:
                                 #   _ident_group_FA2 = _ident_group_FA[_ident_group_FA['FA'].isin([sn1, sn2, sn3, sn4, sn5, sn6])]
                                 #   _ident_group_FA2=_ident_group_FA2.append(_ident_group_FA[_ident_group_FA['Type'] =='FA'])
                                 #   _ident_group_FA2= _ident_group_FA2.append(_ident_group_FA[_ident_group_FA['Type']== 'Lyso'])
                                 #   _ident_group_FA2=_ident_group_FA2.sort_values(by='ppm_abs', ascending=True).head(1).reset_index(drop=True)

                                if _ident_group_FA.loc[0, 'FA'] in lyso_w_ident_lst:
                                        _fa_info= re.split('/', _ident_group_FA.loc[0, 'FA'])
                                        #print _fa_info
                                        #geo=_ident_group_FA2.loc[0, 'FA']
                                        #print geo[:4]
                                        if (_fa_info[0] == _sn1_abbr and _fa_info[1] == _sn2_abbr) or (_fa_info[0] == _sn2_abbr and _fa_info== _sn1_abbr):
                                            _rank_lw_sn1 = lyso_w_ident_lst.index(_ident_group_FA.loc[0,'FA'])
                                            lipid_abbr_df.set_value(_i_abbr, 'i_M-(sn1+sn2)', 2*100*lyso_w_i_lst[_rank_lw_sn1]/ms2_max_i)
                                            lipid_abbr_df.set_value(_i_abbr, 'M-(sn1+sn2)', weight_dct['M-(sn1+sn2)']*2*(10 - _rank_l_sn1)/10)
                                        elif (_fa_info[0] == _sn1_abbr and _fa_info[1] == _sn3_abbr) or (_fa_info[0] == _sn3_abbr and _fa_info[1] == _sn1_abbr):
                                            _rank_lw_sn3 = lyso_w_ident_lst.index(_ident_group_FA.loc[0, 'FA'])
                                            lipid_abbr_df.set_value(_i_abbr, 'i_M-(sn1+sn3)',
                                                                    2 * 100 * lyso_w_i_lst[_rank_lw_sn3] / ms2_max_i)
                                            lipid_abbr_df.set_value(_i_abbr, 'M-(sn1+sn3)',
                                                                    weight_dct['M-(sn1+sn3)'] * 2 * (
                                                                    10 - _rank_lw_sn3) / 10)
                                        elif (_fa_info[0] == _sn2_abbr and _fa_info[1] == _sn3_abbr) or (_fa_info[0] == _sn3_abbr and _fa_info[1] == _sn2_abbr):
                                            _rank_lw_sn2 = lyso_w_ident_lst.index(_ident_group_FA.loc[0, 'FA'])
                                            lipid_abbr_df.set_value(_i_abbr, 'i_M-(sn2+sn3)',
                                                                    2 * 100 * lyso_w_i_lst[_rank_lw_sn2] / ms2_max_i)
                                            lipid_abbr_df.set_value(_i_abbr, 'M-(sn2+sn3)',
                                                                    weight_dct['M-(sn2+sn3)'] * 2 * (
                                                                        10 - _rank_lw_sn2) / 10)

                    lipid_abbr_df['Score'] = lipid_abbr_df[weight_type_lst].sum(axis=1, numeric_only=True)
                    match_reporter = 1
            else:
                for _i_abbr, _abbr_se in lipid_abbr_df.iterrows():
                    # _pl_abbr = _abbr_se['Lipid_abbr']
                    _sn1_abbr = _abbr_se['sn1_abbr']
                    _sn2_abbr = _abbr_se['sn2_abbr']

                    if _sn1_abbr in fa_ident_lst:
                        _rank_sn1 = fa_ident_lst.index(_sn1_abbr)
                        lipid_abbr_df.set_value(_i_abbr, 'i_sn1', 100 * fa_i_lst[_rank_sn1] / ms2_max_i)
                        lipid_abbr_df.set_value(_i_abbr, 'sn1',
                                                weight_dct['sn1'] * (10 - _rank_sn1) / 10
                                                )
                    if _sn2_abbr in fa_ident_lst:
                        _rank_sn2 = fa_ident_lst.index(_sn2_abbr)
                        lipid_abbr_df.set_value(_i_abbr, 'i_sn2', 100 * fa_i_lst[_rank_sn2] / ms2_max_i)
                        lipid_abbr_df.set_value(_i_abbr, 'sn2',
                                                weight_dct['sn2'] * (10 - _rank_sn2) / 10
                                                )
                    if _sn1_abbr in lyso_ident_lst:
                        _rank_l_sn1 = lyso_ident_lst.index(_sn1_abbr)
                        lipid_abbr_df.set_value(_i_abbr, 'i_M-sn1', 100 * lyso_i_lst[_rank_l_sn1] / ms2_max_i)
                        lipid_abbr_df.set_value(_i_abbr, 'M-sn1',
                                                weight_dct['M-sn1'] * (10 - _rank_l_sn1) / 10)
                    if _sn2_abbr in lyso_ident_lst:
                        _rank_l_sn2 = lyso_ident_lst.index(_sn2_abbr)
                        lipid_abbr_df.set_value(_i_abbr, 'i_M-sn2', 100 * lyso_i_lst[_rank_l_sn2] / ms2_max_i)
                        lipid_abbr_df.set_value(_i_abbr, 'M-sn2',
                                                weight_dct['M-sn2'] * (10 - _rank_l_sn2) / 10)
                    if _sn1_abbr in lyso_w_ident_lst:
                        _rank_lw_sn1 = lyso_w_ident_lst.index(_sn1_abbr)
                        lipid_abbr_df.set_value(_i_abbr, 'i_M-(sn1-H2O)', 100 * lyso_w_i_lst[_rank_lw_sn1] / ms2_max_i)
                        lipid_abbr_df.set_value(_i_abbr, 'M-(sn1-H2O)',
                                                weight_dct['M-(sn1-H2O)'] * (10 - _rank_lw_sn1) / 10)
                    if _sn2_abbr in lyso_w_ident_lst:
                        _rank_lw_sn2 = lyso_w_ident_lst.index(_sn2_abbr)
                        lipid_abbr_df.set_value(_i_abbr, 'i_M-(sn2-H2O)', 100 * lyso_w_i_lst[_rank_lw_sn2] / ms2_max_i)
                        lipid_abbr_df.set_value(_i_abbr, 'M-(sn2-H2O)',
                                                weight_dct['M-(sn2-H2O)'] * (10 - _rank_lw_sn2) / 10)
                lipid_abbr_df['Score'] = lipid_abbr_df[weight_type_lst].sum(axis=1, numeric_only=True)
                match_reporter = 1
        else:
            print('!!!!!! NO FA identified =====>--> Skip >>> >>>')

        match_info_dct = {'MATCH_INFO': match_reporter, 'SCORE_INFO': lipid_abbr_df, 'FA_INFO': fa_ident_df,
                          'LYSO_INFO': lyso_ident_df, 'LYSO_W_INFO': lyso_w_ident_df}
        match_factor = match_info_dct['MATCH_INFO']

        return match_info_dct

    def get_specific_peaks(self, mz_lib, ms2_df, ms2_precision=50e-6, ms2_threshold=10):

        _target_frag_df = pd.DataFrame()
        _target_nl_df = pd.DataFrame()
        _other_frag_df = pd.DataFrame()
        _other_nl_df = pd.DataFrame()

        ms2_max_i = ms2_df['i'].max()

        for _i, _frag_se in self.target_frag_df.iterrows():

            _frag_mz = _frag_se['EXACTMASS']
            _frag_class = _frag_se['CLASS']
            _frag_label = _frag_se['LABEL']

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
    fa_list_csv = r'C:\Users\UserPc\Desktop\PhD\LipidHunter_beta_exe\lipidhunter\ConfigurationFiles\FA_list.csv'
    score_cfg = r'C:\Users\UserPc\Desktop\PhD\LipidHunter_beta_exe\lipidhunter\ConfigurationFiles\Score_cfg.xlsx'
    lipid_specifc_cfg = r'C:\Users\UserPc\Desktop\PhD\LipidHunter_beta_exe\lipidhunter\ConfigurationFiles\PL_specific_ion_cfg.xlsx'

    usr_fa_def_df = pd.read_csv(fa_list_csv)
    usr_weight_df = pd.read_excel(score_cfg)
    usr_key_frag_df = pd.read_excel(lipid_specifc_cfg)
    usr_fa_def_df['C'] = usr_fa_def_df['C'].astype(int)
    usr_fa_def_df['DB'] = usr_fa_def_df['DB'].astype(int)

    # for usr_abbr in ['PC(O-36:3)', 'PC(P-38:4)', 'PC(40:5)']:
    #     struc_df = get_structure(usr_fa_def_df, usr_abbr)
    #     print(struc_df)

    usr_lp = 'PC(40:5)'
    usr_lipid_type='PC'
    score_calc = ScoreGenerator(usr_fa_def_df, usr_weight_df, usr_key_frag_df, usr_lipid_type)
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
