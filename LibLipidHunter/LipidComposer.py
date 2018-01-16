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
# Ni, Zhixu, Georgia Angelidou, Mike Lange, Ralf Hoffmann, and Maria Fedorova.
# "LipidHunter identifies phospholipids by high-throughput processing of LC-MS and shotgun lipidomics datasets."
# Analytical Chemistry (2017).
# DOI: 10.1021/acs.analchem.7b01126
#
# For more info please contact:
#     SysMedOs_team: oxlpp@bbz.uni-leipzig.de
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de
#     Developer Georgia Angelidou georgia.angelidou@uni-leipzig.de

from __future__ import print_function

import itertools

import pandas as pd

from LibLipidHunter.LipidNomenclature import NameParserFA
from LibLipidHunter.AbbrElemCalc import ElemCalc
from LibLipidHunter.ParallelFunc import ppm_window_para


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
        self.elem_dct = {'H': [1.0078250321, 0.999885],
                         'D': [2.0141017780, 0.0001157],
                         'C': [12.0, 0.9893],
                         'N': [14.0030740052, 0.99632],
                         'O': [15.9949146221, 0.99757],
                         'Na': [22.98976967, 1.0],
                         'P': [30.97376151, 1.0],
                         'S': [31.97207069, 0.9493],
                         'K': [38.9637069, 0.932581]}

    @staticmethod
    def calc_fa_df(lipid_class, fa_df):

        sn_units_lst = []

        header_lst = fa_df.columns.values.tolist()

        if lipid_class in ['PA', 'PC', 'PE', 'PG', 'PI', 'PS', 'DG']:
            if 'PL' in header_lst and 'FattyAcid' in header_lst:
                pl_sn1_df = fa_df.query('PL == "T" and sn1 == "T"')
                pl_sn2_df = fa_df.query('PL == "T" and sn2 == "T"')

                pl_sn1_lst = pl_sn1_df['FattyAcid'].tolist()
                pl_sn2_lst = pl_sn2_df['FattyAcid'].tolist()

                sn_units_lst = [pl_sn1_lst, pl_sn2_lst]

        elif lipid_class in ['TG']:
            if 'TG' in header_lst and 'FattyAcid' in header_lst:
                tg_sn1_df = fa_df.query('TG == "T" and sn1 == "T"')
                tg_sn2_df = fa_df.query('TG == "T" and sn2 == "T"')
                tg_sn3_df = fa_df.query('TG == "T" and sn3 == "T"')

                tg_sn1_lst = tg_sn1_df['FattyAcid'].tolist()
                tg_sn2_lst = tg_sn2_df['FattyAcid'].tolist()
                tg_sn3_lst = tg_sn3_df['FattyAcid'].tolist()

                sn_units_lst = [tg_sn1_lst, tg_sn2_lst, tg_sn3_lst]

        elif lipid_class in ['CL']:
            if 'CL' in header_lst and 'FattyAcid' in header_lst:
                cl_sn1_df = fa_df.query('CL == "T" and sn1 == "T"')
                cl_sn2_df = fa_df.query('CL == "T" and sn2 == "T"')
                cl_sn3_df = fa_df.query('CL == "T" and sn3 == "T"')
                cl_sn4_df = fa_df.query('CL == "T" and sn4 == "T"')

                cl_sn1_lst = cl_sn1_df['FattyAcid'].tolist()
                cl_sn2_lst = cl_sn2_df['FattyAcid'].tolist()
                cl_sn3_lst = cl_sn3_df['FattyAcid'].tolist()
                cl_sn4_lst = cl_sn4_df['FattyAcid'].tolist()

                sn_units_lst = [cl_sn1_lst, cl_sn2_lst, cl_sn3_lst, cl_sn4_lst]

        return sn_units_lst

    def calc_fa_query(self, lipid_type, fa_whitelist, ms2_ppm=100):

        usr_fa_df = pd.read_excel(fa_whitelist)
        usr_fa_df = usr_fa_df.fillna(value='F')

        sn_units_lst = self.calc_fa_df(lipid_type, usr_fa_df)
        fa_abbr_lst = []
        for _s in sn_units_lst:
            fa_abbr_lst.extend(_s)
        fa_abbr_lst = sorted(list(set(fa_abbr_lst)))

        abbr_parser = NameParserFA()
        elem_calc = ElemCalc()
        usr_fa_dct = {}
        for _fa_abbr in fa_abbr_lst:
            _fa_info_dct = abbr_parser.get_fa_info(_fa_abbr)
            _lipid_formula, _lipid_elem_dct = elem_calc.get_formula(_fa_abbr)
            _fa_info_dct['ABBR'] = _fa_abbr
            _fa_info_dct['FORMULA'] = _lipid_formula
            _fa_info_dct['EXACTMASS'] = elem_calc.get_exactmass(_lipid_elem_dct)
            usr_fa_dct[_fa_abbr] = _fa_info_dct

        usr_fa_df = pd.DataFrame(usr_fa_dct).T.copy()
        usr_fa_df.is_copy = False

        for _fa_ion in ['[FA-H]-', '[FA-H2O-H]-', '[FA-H2O+H]+']:
            usr_fa_df['%s_MZ_LOW' % _fa_ion] = ppm_window_para(usr_fa_df['%s_MZ' % _fa_ion].values.tolist(),
                                                               ms2_ppm * -1)
            usr_fa_df['%s_MZ_HIGH' % _fa_ion] = ppm_window_para(usr_fa_df['%s_MZ' % _fa_ion].values.tolist(), ms2_ppm)
            usr_fa_df['%s_Q' % _fa_ion] = (usr_fa_df['%s_MZ_LOW' % _fa_ion].astype(str) + ' <= mz <= '
                                           + usr_fa_df['%s_MZ_HIGH' % _fa_ion].astype(str))

        if lipid_type in ['PA', 'PC', 'PE', 'PG', 'PI', 'PS', 'PIP']:

            lyso_type_dct = {'[L%s-H]-' % lipid_type: 'EXACTMASS', '[L%s-H2O-H]-' % lipid_type: '[FA-H2O]_MZ'}

            lyso_base_elem_dct = self.lipid_hg_elem_dct[lipid_type]
            # create the phospholipid head group
            for _e in self.glycerol_bone_elem_dct.keys():
                lyso_base_elem_dct[_e] += self.glycerol_bone_elem_dct[_e]

            # the element here is with no Hydroxyl on sn1 and sn2, here a [M-H]- is already considered
            lyso_base_mz = elem_calc.get_exactmass(lyso_base_elem_dct) + 1.0078250321 + 15.9949146221

            if lipid_type in ['PC', 'SM']:
                lyso_base_mz -= (12.0 + 2 * 1.0078250321)  # LPC loss one -CH3 from HG (one H already remove above)

            for _lyso_ion in lyso_type_dct.keys():
                if lipid_type in ['PC', 'SM']:

                    if lyso_type_dct[_lyso_ion] == 'EXACTMASS':
                        usr_fa_df['%s_ABBR' % _lyso_ion] = ('[L' + lipid_type + '(' + usr_fa_df['ABBR'].str.strip('FA')
                                                            + ')-CH3]-')
                    elif lyso_type_dct[_lyso_ion] == '[FA-H2O]_MZ':
                        usr_fa_df['%s_ABBR' % _lyso_ion] = ('[L' + lipid_type + '(' + usr_fa_df['ABBR'].str.strip('FA')
                                                            + ')-H2O-CH3]-')
                    else:
                        usr_fa_df['%s_ABBR' % _lyso_ion] = 'ERROR'
                else:
                    if lyso_type_dct[_lyso_ion] == 'EXACTMASS':
                        usr_fa_df['%s_ABBR' % _lyso_ion] = ('[L' + lipid_type + '(' + usr_fa_df['ABBR'].str.strip('FA')
                                                            + ')-H]-')
                    elif lyso_type_dct[_lyso_ion] == '[FA-H2O]_MZ':
                        usr_fa_df['%s_ABBR' % _lyso_ion] = ('[L' + lipid_type + '(' + usr_fa_df['ABBR'].str.strip('FA')
                                                            + ')-H2O-H]-')
                    else:
                        usr_fa_df['%s_ABBR' % _lyso_ion] = 'ERROR'

                usr_fa_df['%s_MZ' % _lyso_ion] = lyso_base_mz + usr_fa_df[lyso_type_dct[_lyso_ion]]
                usr_fa_df['%s_MZ_LOW' % _lyso_ion] = ppm_window_para(usr_fa_df['%s_MZ' % _lyso_ion].values.tolist(),
                                                                     ms2_ppm * -1)
                usr_fa_df['%s_MZ_HIGH' % _lyso_ion] = ppm_window_para(usr_fa_df['%s_MZ' % _lyso_ion].values.tolist(),
                                                                      ms2_ppm)
                usr_fa_df['%s_Q' % _lyso_ion] = (usr_fa_df['%s_MZ_LOW' % _lyso_ion].astype(str) + ' <= mz <= '
                                                 + usr_fa_df['%s_MZ_HIGH' % _lyso_ion].astype(str))
        elif lipid_type in ['TG']:
            # TODO(georgia.angelidou@uni-leipzig.de): create the section for theuniue fragments when there is TG
            mg_type_dct = {'[MG-H2O+H]+': 'EXACTMASS'}
            mg_base_elem_dct = self.lipid_hg_elem_dct[lipid_type]

            for _e in self.glycerol_bone_elem_dct.keys():
                mg_base_elem_dct[_e] += self.glycerol_bone_elem_dct[_e]

            # Calculate the rest of monoglycerol after the neutral loss of the FA in protonated form
            mg_base_elem_dct = elem_calc.get_exactmass(mg_base_elem_dct) + 15.9949146221 + (3*1.0078250321)
            for _mg_ion in mg_type_dct.keys():
                usr_fa_df['%s_ABBR' % _mg_ion] = ('[MG(' + usr_fa_df['ABBR'].str.strip('FA') + ')-H2O+H]+')
                usr_fa_df['%s_MZ' % _mg_ion] = mg_base_elem_dct + usr_fa_df[mg_type_dct[_mg_ion]]
                usr_fa_df['%s_MZ_LOW' % _mg_ion] = ppm_window_para(usr_fa_df['%s_MZ' % _mg_ion].values.tolist(),
                                                                   ms2_ppm * -1)
                usr_fa_df['%s_MZ_HIGH' % _mg_ion] = ppm_window_para(usr_fa_df['%s_MZ' % _mg_ion].values.tolist(),
                                                                    ms2_ppm)
                usr_fa_df['%s_Q' % _mg_ion] = (usr_fa_df['%s_MZ_LOW' % _mg_ion].astype(str) + ' <= mz <= ' + usr_fa_df[
                    '%s_MZ_HIGH' % _mg_ion].astype(str))
        return usr_fa_df

    def gen_all_comb(self, lipid_class, usr_fa_df, position=False):

        sn_units_lst = self.calc_fa_df(lipid_class, usr_fa_df)

        if lipid_class in ['PA', 'PC', 'PE', 'PG', 'PI', 'PS', 'DG'] and len(sn_units_lst) == 2:
            sn_comb_lst = list(itertools.product(sn_units_lst[0], sn_units_lst[1]))
            # lipid_template = '{}'
        elif lipid_class == 'TG' and len(sn_units_lst) == 3:
            sn_comb_lst = list(itertools.product(sn_units_lst[0], sn_units_lst[1], sn_units_lst[2]))
        elif lipid_class == 'CL' and len(sn_units_lst) == 4:
            sn_comb_lst = list(itertools.product(sn_units_lst[0], sn_units_lst[1], sn_units_lst[2], sn_units_lst[3]))
        else:
            sn_comb_lst = []

        sn_comb_lite_lst = []
        sn_comb_rm_lst = []

        if position is False:
            for _comb in sn_comb_lst:
                _rev_comb = tuple(sorted(list(_comb)))
                if _comb not in sn_comb_lite_lst and _rev_comb not in sn_comb_lite_lst:
                    sn_comb_lite_lst.append(_comb)
                else:
                    sn_comb_rm_lst.append(_comb)
                    # sn_comb_rm_lst.append(_rev_comb)
        else:
            sn_comb_lite_lst = sn_comb_lst

        lipid_comb_dct = {}

        if lipid_class in ['PA', 'PC', 'PE', 'PG', 'PI', 'PS', 'DG'] and len(sn_comb_lite_lst) > 0:
            for _comb_lite in sn_comb_lite_lst:
                _lipid_abbr = '{pl}({sn1}_{sn2})'.format(pl=lipid_class, sn1=_comb_lite[0].strip('FA'),
                                                         sn2=_comb_lite[1].strip('FA'))
                lipid_comb_dct[_lipid_abbr] = {'CLASS': lipid_class, 'SN1': _comb_lite[0], 'SN2': _comb_lite[1],
                                               'DISCRETE_ABBR': _lipid_abbr}
        elif lipid_class in ['TG'] and len(sn_comb_lite_lst) > 0:
            for _comb_lite in sn_comb_lite_lst:
                _lipid_abbr = '{pl}({sn1}_{sn2}_{sn3})'.format(pl=lipid_class, sn1=_comb_lite[0].strip('FA'),
                                                               sn2=_comb_lite[1].strip('FA'),
                                                               sn3=_comb_lite[2].strip('FA'))

                lipid_comb_dct[_lipid_abbr] = {'CLASS': lipid_class, 'SN1': _comb_lite[0], 'SN2': _comb_lite[1],
                                               'SN3': _comb_lite[2], 'DISCRETE_ABBR': _lipid_abbr}
        else:
            pass
        # print(sn_units_lst)
        # print(len(sn_comb_lst))
        # print(sn_comb_lst)
        # print(len(sn_comb_lite_lst))
        # print(sn_comb_lite_lst)
        # print(len(sn_comb_rm_lst))
        # print(sorted(sn_comb_rm_lst))
        # print(len(lipid_comb_dct.keys()))
        # print(lipid_comb_dct)

        return lipid_comb_dct

    @staticmethod
    def calc_fragments(lipid_dct):

        # m_formula = lipid_dct['FORMULA']
        m_exactmass = lipid_dct['EXACTMASS']
        m_class = lipid_dct['CLASS']
        h_exactmass = 1.0078250321
        ch3_exactmass = 12.0 + 3 * 1.0078250321
        nl_water = 2 * 1.0078250321 + 15.9949146221
        sn1_abbr = lipid_dct['SN1'].strip('FA')
        sn2_abbr = lipid_dct['SN2'].strip('FA')

        sn1_exactmass = lipid_dct['SN1_EXACTMASS']
        sn2_exactmass = lipid_dct['SN2_EXACTMASS']

        #################################################33
        #
        #   Note:
        #   Why calculated the below fragments when there were all ready calculate before in the above section
        #   In calc_fa_query library    Line 126
        #
        #####################################################
        if m_class in ['PA', 'PE', 'PG', 'PI', 'PS']:
            sn1_abbr = lipid_dct['SN1'].strip('FA')
            sn2_abbr = lipid_dct['SN2'].strip('FA')

            sn1_exactmass = lipid_dct['SN1_EXACTMASS']
            sn2_exactmass = lipid_dct['SN2_EXACTMASS']
            lyso_str = 'L' + m_class

            lipid_dct['[LPL(SN1)-H]-_ABBR'] = '[%s(%s)-H]-' % (lyso_str, sn1_abbr)
            lipid_dct['[LPL(SN2)-H]-_ABBR'] = '[%s(%s)-H]-' % (lyso_str, sn2_abbr)
            lipid_dct['[LPL(SN1)-H2O-H]-_ABBR'] = '[%s(%s)-H2O-H]-' % (lyso_str, sn1_abbr)
            lipid_dct['[LPL(SN2)-H2O-H]-_ABBR'] = '[%s(%s)-H2O-H]-' % (lyso_str, sn2_abbr)

            lipid_dct['[LPL(SN1)-H]-_MZ'] = round(m_exactmass - (sn2_exactmass - nl_water) - h_exactmass, 6)
            lipid_dct['[LPL(SN2)-H]-_MZ'] = round(m_exactmass - (sn1_exactmass - nl_water) - h_exactmass, 6)
            lipid_dct['[LPL(SN1)-H2O-H]-_MZ'] = round(m_exactmass - sn2_exactmass - h_exactmass, 6)
            lipid_dct['[LPL(SN2)-H2O-H]-_MZ'] = round(m_exactmass - sn1_exactmass - h_exactmass, 6)

        elif m_class in ['PC']:
            sn1_abbr = lipid_dct['SN1'].strip('FA')
            sn2_abbr = lipid_dct['SN2'].strip('FA')

            sn1_exactmass = lipid_dct['SN1_EXACTMASS']
            sn2_exactmass = lipid_dct['SN2_EXACTMASS']

            lyso_str = 'L' + m_class
            # The abbr. here is not exactly correct due to the compatibility issues with ranks core calc functions
            lipid_dct['[LPL(SN1)-H]-_ABBR'] = '[%s(%s)-CH3]-' % (lyso_str, sn1_abbr)
            lipid_dct['[LPL(SN2)-H]-_ABBR'] = '[%s(%s)-CH3]-' % (lyso_str, sn2_abbr)
            lipid_dct['[LPL(SN1)-H2O-H]-_ABBR'] = '[%s(%s)-H2O-CH3]-' % (lyso_str, sn1_abbr)
            lipid_dct['[LPL(SN2)-H2O-H]-_ABBR'] = '[%s(%s)-H2O-CH3]-' % (lyso_str, sn2_abbr)

            lipid_dct['[LPL(SN1)-H]-_MZ'] = round(m_exactmass - (sn2_exactmass - nl_water) - ch3_exactmass, 6)
            lipid_dct['[LPL(SN2)-H]-_MZ'] = round(m_exactmass - (sn1_exactmass - nl_water) - ch3_exactmass, 6)
            lipid_dct['[LPL(SN1)-H2O-H]-_MZ'] = round(m_exactmass - sn2_exactmass - ch3_exactmass, 6)
            lipid_dct['[LPL(SN2)-H2O-H]-_MZ'] = round(m_exactmass - sn1_exactmass - ch3_exactmass, 6)

        else:
            ############### Here maybe should get the DG fragments
            # TODO(georgia.angelidou@uni-leipzig.de): create the section for theuniue fragments when there is TG
            #   Missing the fragments for the sodium adduct
            sn3_abbr = lipid_dct['SN3'].strip('FA')
            sn3_exactmass = lipid_dct['SN3_EXACTMASS']

            mg_str = 'MG'
            ######## Deactivated this part. Instead calculate the DG
            # lipid_dct['[MG(SN1)-H2O+H]+_ABBR'] = '[%s(%s)-H2O+H]+' % (mg_str, sn1_abbr)
            # lipid_dct['[MG(SN2)-H2O+H]+_ABBR'] = '[%s(%s)-H2O+H]+' % (mg_str, sn2_abbr)
            # lipid_dct['[MG(SN3)-H2O+H]+_ABBR'] = '[%s(%s)-H2O+H]+' % (mg_str, sn3_abbr)

            #lipid_dct['MG(SN1)-H2O+H]+_MZ'] = round(m_exactmass - ())

            dg_str = 'M'
            ###########################################################3
            #
            #   Ask if it should be FA(16:0) or without the letters FA
            #
            ###############################################################
            lipid_dct['[M-(SN1)+H]+_ABBR'] = '[%s-FA(%s)+H]+' % (dg_str, sn1_abbr)
            lipid_dct['[M-(SN2)+H]+_ABBR'] = '[%s-FA(%s)+H]+' % (dg_str, sn2_abbr)
            lipid_dct['[M-(SN3)+H]+_ABBR'] = '[%s-FA(%s)+H]+' % (dg_str, sn3_abbr)

            lipid_dct['[M-(SN1)+H]+_MZ'] = round(m_exactmass - sn1_exactmass + h_exactmass, 6)
            lipid_dct['[M-(SN2)+H]+_MZ'] = round(m_exactmass - sn2_exactmass + h_exactmass, 6)
            lipid_dct['[M-(SN3)+H]+_MZ'] = round(m_exactmass - sn3_exactmass + h_exactmass, 6)

            lipid_dct['[M-(SN1-H2O)+H]+_ABBR'] = '[%s-(FA(%s)-H2O)+H]+' % (dg_str, sn1_abbr)
            lipid_dct['[M-(SN2-H2O)+H]+_ABBR'] = '[%s-(FA(%s)-H2O)+H]+' % (dg_str, sn2_abbr)
            lipid_dct['[M-(SN3-H2O)+H]+_ABBR'] = '[%s-(FA(%s)-H2O)+H]+' % (dg_str, sn3_abbr)

            lipid_dct['[M-(SN1-H2O)+H]+_MZ'] = round(m_exactmass - (sn1_exactmass - nl_water) + h_exactmass, 6)
            lipid_dct['[M-(SN2-H2O)+H]+_MZ'] = round(m_exactmass - (sn2_exactmass - nl_water) + h_exactmass, 6)
            lipid_dct['[M-(SN3-H2O)+H]+_MZ'] = round(m_exactmass - (sn3_exactmass - nl_water) + h_exactmass, 6)

        return lipid_dct

    def compose_lipid(self, param_dct):

        lipid_class = param_dct['lipid_type']
        lipid_charge = param_dct['charge_mode']
        if param_dct['exact_position'] == 'TRUE':
            position_set = True
        else:
            position_set = False

        usr_fa_df = pd.read_excel(param_dct['fa_whitelist'])
        usr_fa_df = usr_fa_df.fillna(value='F')
        print('=== ==> --> FA white list loaded >>>')
        # print(usr_fa_df)
        lipid_comb_dct = self.gen_all_comb(lipid_class, usr_fa_df, position=position_set)

        lipid_info_dct = {}

        abbr_parser = NameParserFA()
        elem_calc = ElemCalc()
        # TODO(georgia.angelidou@uni-leipzig.de): need to add the lipid class to support more than 2 sn
        for _lipid in lipid_comb_dct.keys():
            _lipid_dct = lipid_comb_dct[_lipid]

            _sn1_abbr = _lipid_dct['SN1']
            _sn2_abbr = _lipid_dct['SN2']
            _sn1_info_dct = abbr_parser.get_fa_info(_sn1_abbr)
            _sn2_info_dct = abbr_parser.get_fa_info(_sn2_abbr)

            for _sn1_k in _sn1_info_dct.keys():
                _lipid_dct['SN1_' + _sn1_k] = _sn1_info_dct[_sn1_k]

            for _sn2_k in _sn1_info_dct.keys():
                _lipid_dct['SN2_' + _sn2_k] = _sn2_info_dct[_sn2_k]

            if lipid_class in ['PA', 'PC', 'PE', 'PG', 'PI', 'PS', 'DG']:
                _lipid_dct['M_DB'] = _sn1_info_dct['DB'] + _sn2_info_dct['DB']
                # TODO(georgia.angelidou@uni-leipzig.de): not important (just keep in mind for future correction)
                # consideration the case if the user choose the sn2 position for the different link types
                # or in that they may mistype something
                if _sn1_info_dct['LINK'] in ['FA', 'A']:
                    lipid_bulk_str = '{pl}({c}:{db})'.format(pl=lipid_class,
                                                             c=_sn1_info_dct['C'] + _sn2_info_dct['C'],
                                                             db=lipid_comb_dct[_lipid]['M_DB'])
                else:
                    lipid_bulk_str = '{pl}({lk}{c}:{db})'.format(pl=lipid_class, lk=_sn1_info_dct['LINK'],
                                                                 c=_sn1_info_dct['C'] + _sn2_info_dct['C'],
                                                                 db=lipid_comb_dct[_lipid]['M_DB'])
            elif lipid_class in ['TG']:
                _sn3_abbr = _lipid_dct['SN3']
                _sn3_info_dct = abbr_parser.get_fa_info(_sn3_abbr)
                for _sn3_k in _sn3_info_dct.keys():
                    _lipid_dct['SN3_' + _sn3_k] = _sn3_info_dct[_sn3_k]

                _lipid_dct['M_DB'] = _sn3_info_dct['DB'] + _sn2_info_dct['DB'] + _sn1_info_dct['DB']
                # Note: For TG in the current default not consider the different lipids with other type of bond
                # If stay like this need to be mention in somewhere for the user
                lipid_bulk_str = '{pl}({c}:{db})'.format(pl=lipid_class,
                                                         c=_sn1_info_dct['C'] + _sn2_info_dct['C'] + _sn3_info_dct['C'],
                                                         db=lipid_comb_dct[_lipid]['M_DB'])

            _lipid_dct['BULK_ABBR'] = lipid_bulk_str

            _lipid_formula, _lipid_elem_dct = elem_calc.get_formula(lipid_bulk_str)

            _lipid_dct['FORMULA'] = _lipid_formula
            _lipid_dct['EXACTMASS'] = elem_calc.get_exactmass(_lipid_elem_dct)
            for _elem_k in _lipid_elem_dct.keys():
                _lipid_dct['M_' + _elem_k] = _lipid_elem_dct[_elem_k]

            # charged
            _chg_lipid_formula, _chg_lipid_elem_dct = elem_calc.get_formula(lipid_bulk_str, charge=lipid_charge)
            _lipid_dct[lipid_charge + '_FORMULA'] = _chg_lipid_formula
            _lipid_dct[lipid_charge + '_MZ'] = elem_calc.get_exactmass(_chg_lipid_elem_dct)

            # fragments

            _lipid_dct = self.calc_fragments(_lipid_dct)

            lipid_info_dct[_lipid] = _lipid_dct
            del _lipid_dct

        lipid_master_df = pd.DataFrame(lipid_comb_dct).T
        lipid_master_df.reset_index(drop=True, inplace=True)

        return lipid_master_df


if __name__ == '__main__':

    fa_lst_file = r'../ConfigurationFiles/FA_Whitelist.xlsx'

    # Note:
    # exact position means to consider the poition from the FA white list that the user give but,
    # in the case that the user define 2 different FA for both positions then:
    # When it is false it will give only one option
    # and when it is TRUE to give both compinations that these 2 FA an make (incase of phospholipids)
    usr_param_dct = {'fa_whitelist': fa_lst_file, 'lipid_type': 'TG', 'charge_mode': '[M+H]+','exact_position': 'FALSE'}

    composer = LipidComposer()
    usr_lipid_master_df = composer.compose_lipid(param_dct=usr_param_dct)

    # print(usr_lipid_master_df.shape)
    # print(usr_lipid_master_df.head())
    # print(usr_lipid_master_df.tail())

    master_xlsx = r'../Temp/LipidMaster_Whitelist_TG.xlsx'
    fa_xlsx = r'../Temp/LipidMaster_FAlist.xlsx'

    calc_fa_df = composer.calc_fa_query(usr_param_dct['lipid_type'], r'../ConfigurationFiles/FA_Whitelist.xlsx',
                                        ms2_ppm=50)

    print(calc_fa_df)
    usr_lipid_master_df.to_excel(master_xlsx)
    calc_fa_df.to_excel(fa_xlsx)
