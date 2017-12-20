# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2017  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     SysMedOs_team: oxlpp@bbz.uni-leipzig.de
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import itertools

import pandas as pd

from AbbrParser import AbbrParser
from AbbrElemCalc import ElemCalc


class LipidComposer:

    def __init__(self):
        pa_hg_elem = {'C': 0, 'H': 3, 'O': 4, 'P': 1, 'N': 0}
        pc_hg_elem = {'C': 5, 'H': 14, 'O': 4, 'P': 1, 'N': 1}
        pe_hg_elem = {'C': 2, 'H': 8, 'O': 4, 'P': 1, 'N': 1}
        pg_hg_elem = {'C': 3, 'H': 9, 'O': 6, 'P': 1, 'N': 0}
        pi_hg_elem = {'C': 6, 'H': 13, 'O': 9, 'P': 1, 'N': 0}
        pip_hg_elem = {'C': 6, 'H': 14, 'O': 12, 'P': 2, 'N': 0}
        ps_hg_elem = {'C': 3, 'H': 8, 'O': 6, 'P': 1, 'N': 1}
        tg_hg_elem = {'C': 0, 'H': 0, 'O': 0, 'P': 0, 'N': 0}

        self.lipid_hg_lst = ['PA', 'PC', 'PE', 'PG', 'PI', 'PS', 'PIP', 'TG']

        self.lipid_hg_elem_dct = {'PA': pa_hg_elem, 'PC': pc_hg_elem, 'PE': pe_hg_elem, 'PG': pg_hg_elem,
                                  'PI': pi_hg_elem, 'PS': ps_hg_elem, 'PIP': pip_hg_elem, 'TG': tg_hg_elem}

        self.glycerol_bone_elem_dct = {'C': 3, 'H': 2}
        self.link_o_elem_dct = {'O': -1, 'H': 2}
        self.link_p_elem_dct = {'O': -1}

    def calc_fa_df(self, lipid_class, usr_fa_df):

        sn_units_lst = []

        header_lst = usr_fa_df.columns.values.tolist()

        if lipid_class in ['PA', 'PC', 'PE', 'PG', 'PI', 'PS', 'DG']:
            if 'PL' in header_lst and 'FattyAcid' in header_lst:
                pl_sn1_df = usr_fa_df.query('PL == "T" and sn1 == "T"')
                pl_sn2_df = usr_fa_df.query('PL == "T" and sn2 == "T"')

                pl_sn1_lst = pl_sn1_df['FattyAcid'].tolist()
                pl_sn2_lst = pl_sn2_df['FattyAcid'].tolist()

                sn_units_lst = [pl_sn1_lst, pl_sn2_lst]

        elif lipid_class in ['TG']:
            if 'TG' in header_lst and 'FattyAcid' in header_lst:
                tg_sn1_df = usr_fa_df.query('TG == "T" and sn1 == "T"')
                tg_sn2_df = usr_fa_df.query('TG == "T" and sn2 == "T"')
                tg_sn3_df = usr_fa_df.query('TG == "T" and sn3 == "T"')

                tg_sn1_lst = tg_sn1_df['FattyAcid'].tolist()
                tg_sn2_lst = tg_sn2_df['FattyAcid'].tolist()
                tg_sn3_lst = tg_sn3_df['FattyAcid'].tolist()

                sn_units_lst = [tg_sn1_lst, tg_sn2_lst, tg_sn3_lst]

        elif lipid_class in ['CL']:
            if 'CL' in header_lst and 'FattyAcid' in header_lst:
                cl_sn1_df = usr_fa_df.query('CL == "T" and sn1 == "T"')
                cl_sn2_df = usr_fa_df.query('CL == "T" and sn2 == "T"')
                cl_sn3_df = usr_fa_df.query('CL == "T" and sn3 == "T"')
                cl_sn4_df = usr_fa_df.query('CL == "T" and sn4 == "T"')

                cl_sn1_lst = cl_sn1_df['FattyAcid'].tolist()
                cl_sn2_lst = cl_sn2_df['FattyAcid'].tolist()
                cl_sn3_lst = cl_sn3_df['FattyAcid'].tolist()
                cl_sn4_lst = cl_sn4_df['FattyAcid'].tolist()

                sn_units_lst = [cl_sn1_lst, cl_sn2_lst, cl_sn3_lst, cl_sn4_lst]

        return sn_units_lst

    def gen_all_comb(self, lipid_class, usr_fa_df):
        sn_units_lst = self.calc_fa_df(lipid_class, usr_fa_df)

        if lipid_class in ['PA', 'PC', 'PE', 'PG', 'PI', 'PS', 'DG'] and len(sn_units_lst) == 2:
            sn_comb_lst = list(itertools.product(sn_units_lst[0], sn_units_lst[1]))
            lipid_template = '{}'
        elif lipid_class == 'TG' and len(sn_units_lst) == 3:
            sn_comb_lst = list(itertools.product(sn_units_lst[0], sn_units_lst[1], sn_units_lst[2]))
        elif lipid_class == 'CL' and len(sn_units_lst) == 4:
            sn_comb_lst = list(itertools.product(sn_units_lst[0], sn_units_lst[1], sn_units_lst[2], sn_units_lst[3]))
        else:
            sn_comb_lst = []

        sn_comb_lite_lst = []
        sn_comb_rm_lst = []

        for _comb in sn_comb_lst:
            _rev_comb = tuple(reversed(list(_comb)))
            if _comb not in sn_comb_lite_lst and _rev_comb not in sn_comb_lite_lst:
                sn_comb_lite_lst.append(_comb)
            else:
                sn_comb_rm_lst.append(_comb)

        lipid_comb_dct = {}

        if lipid_class in ['PA', 'PC', 'PE', 'PG', 'PI', 'PS', 'DG'] and len(sn_comb_lite_lst) > 0:
            for _comb_lite in sn_comb_lite_lst:
                _lipid_abbr = '{pl}({sn1}_{sn2})'.format(pl=lipid_class, sn1=_comb_lite[0].strip('FA'),
                                                         sn2=_comb_lite[1].strip('FA'))
                lipid_comb_dct[_lipid_abbr] = {'CLASS': lipid_class, 'SN1': _comb_lite[0], 'SN2': _comb_lite[1]}
        else:
            pass

        # print(sn_units_lst)
        # print(len(sn_comb_lst))
        # print(sn_comb_lst)
        # print(len(sn_comb_lite_lst))
        # print(sn_comb_lite_lst)
        # print(len(sn_comb_rm_lst))
        # print(sn_comb_rm_lst)
        # print(len(lipid_comb_dct.keys()))
        # print(lipid_comb_dct)

        return lipid_comb_dct

    def compose_lipid(self, param_dct):

        lipid_class = param_dct['lipid_type']
        lipid_charge = param_dct['charge_mode']
        usr_fa_df = pd.read_excel(param_dct['fa_whitelist'])
        usr_fa_df = usr_fa_df.fillna(value='F')
        print(usr_fa_df)
        lipid_comb_dct = self.gen_all_comb(lipid_class, usr_fa_df)

        abbr_parser = AbbrParser()
        elem_calc = ElemCalc()
        for _lipid in lipid_comb_dct.keys():
            _lipid_dct = lipid_comb_dct[_lipid]
            _sn1_abbr = _lipid_dct['SN1']
            _sn2_abbr = _lipid_dct['SN2']
            _sn1_info_dct = abbr_parser.get_fa_info(_sn1_abbr)
            _sn2_info_dct = abbr_parser.get_fa_info(_sn2_abbr)

            for _sn1_k in _sn1_info_dct.keys():
                lipid_comb_dct[_lipid]['SN1_' + _sn1_k] = _sn1_info_dct[_sn1_k]

            for _sn2_k in _sn1_info_dct.keys():
                lipid_comb_dct[_lipid]['SN2_' + _sn2_k] = _sn2_info_dct[_sn2_k]

            lipid_comb_dct[_lipid]['M_DB'] = _sn1_info_dct['DB'] + _sn2_info_dct['DB']
            if _sn1_info_dct['LINK'] in ['FA', 'A']:
                lipid_bulk_str = '{pl}({c}:{db})'.format(pl=lipid_class,
                                                         c=_sn1_info_dct['C'] + _sn2_info_dct['C'],
                                                         db=lipid_comb_dct[_lipid]['M_DB'])
            else:
                lipid_bulk_str = '{pl}({lk}{c}:{db})'.format(pl=lipid_class, lk=_sn1_info_dct['LINK'],
                                                             c=_sn1_info_dct['C'] + _sn2_info_dct['C'],
                                                             db=lipid_comb_dct[_lipid]['M_DB'])

            lipid_comb_dct[_lipid]['Bulk_ABBR'] = lipid_bulk_str

            _lipid_formula, _lipid_elem_dct = elem_calc.get_formula(lipid_bulk_str)
            lipid_comb_dct[_lipid]['FORMULA'] = _lipid_formula
            lipid_comb_dct[_lipid]['EXACTMASS'] = elem_calc.get_exactmass(_lipid_elem_dct)
            for _elem_k in _lipid_elem_dct.keys():
                lipid_comb_dct[_lipid]['M_' + _elem_k] = _lipid_elem_dct[_elem_k]

            # charged
            _chg_lipid_formula, _chg_lipid_elem_dct = elem_calc.get_formula(lipid_bulk_str, charge=lipid_charge)
            lipid_comb_dct[_lipid][lipid_charge + '_FORMULA'] = _chg_lipid_formula
            lipid_comb_dct[_lipid][lipid_charge + '_MZ'] = elem_calc.get_exactmass(_lipid_elem_dct)

        lipid_master_df = pd.DataFrame(lipid_comb_dct).T

        return lipid_master_df


if __name__ == '__main__':
    fa_lst_file = r'../ConfigurationFiles/FA_Whitelist.xlsx'

    usr_param_dct = {'fa_whitelist': fa_lst_file, 'lipid_type': 'PE', 'charge_mode': '[M-H]-'}

    composer = LipidComposer()
    usr_lipid_master_df = composer.compose_lipid(param_dct=usr_param_dct)

    print(usr_lipid_master_df.shape)
    print(usr_lipid_master_df.head())
    print(usr_lipid_master_df.tail())

    master_xlsx = r'../Temp/LipidMaster_Whitelist.xlsx'

    usr_lipid_master_df.to_excel(master_xlsx)
