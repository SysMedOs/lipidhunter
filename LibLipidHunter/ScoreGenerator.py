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
from __future__ import division

import re
from LibLipidHunter.AbbrElemCalc import ElemCalc
from LibLipidHunter.IsotopeHunter import IsotopeHunter
from FAwhiteList import FA_list
import pandas as pd
import itertools


# All the class ProposedStructure can be remove.
# This part know is include in the ScoreGenerator Class


class ProposedStructure:
    def __init__(self, fa_def_df, lipid_type, ion_charge='[M+H]+'):
        self.fa_def_df = fa_def_df
        if ion_charge in ['[M+H]+', '[M+NH4]+', '[M+Na]+']:
            charge_mode = 'POS'
        elif ion_charge in ['[M-H]-', '[M+HCOO]-']:
            charge_mode = 'NEG'
        else:
            charge_mode = 'NEG'
        self.charge_mode = charge_mode


class ScoreGenerator:
    def __init__(self, param_dct, weight_df, key_frag_df, lipid_type, checked_info_df, ion_charge='[M-H]-',
                 ms2_ppm=200):

        # new params
        self.weight_dct = weight_df.to_dict(orient='index')
        self.weight_type_lst = self.weight_dct.keys()
        self.weight_group_lst = list(set(weight_df['Group'].tolist()))
        print('self.weight_group_lst', self.weight_group_lst)

        # end new params

    #
    #         #####################################################################################
    #         #
    #         #   Here should get the part where the program calculates all the possible columns for the FA white list
    #         #
    #         ########################################################################################
    # #         fa_def_df['M+2_i_isomer'] = ''
    # #         for _i, _row in fa_def_df.iterrows():
    # #             elem_dct, elem_comp = FA_list().elemental_composition(_row['Link'], _row['C'], _row['DB'])
    # # #            elem_dct = IsotopeHunter().get_elements(_row['elem'])
    # #             mass_calc = IsotopeHunter().get_mono_mz(elem_dct)
    # #             print mass_calc
    # #             print ('kwioowew')
    # #             fa_def_df.set_value(_i, 'mass', mass_calc)
    # #             print ('wowowje')
    # #             print ('daisy')
    # #             value_dict={}
    # #             print elem_dct
    # #             fa_def_df.set_value(_i, 'elem', elem_comp)
    # #             if ion_charge in ['[M-H]-', '[M+OAc]-', '[M+HCOO]-']:
    # #                 print ('olalala')
    # #                 value_dict = FA_list().negative_column(elem_dct)
    # #                 fa_def_df.set_value(_i, '[M-H]-', value_dict['neg_value'])
    # #                 fa_def_df.set_value(_i, '[M-W-H]-', value_dict['neg_value_w'])
    # #                 fa_def_df.set_value(_i, 'NL-H2O', value_dict['m_water'])
    # #             elif ion_charge in ['[M+H]+', '[M+NH4]+']:
    # #                 print ('doggy')
    # #                 print elem_dct
    # #                 value_dict = FA_list().pos_geo(elem_dct, ion_charge)
    # #                 fa_def_df.set_value(_i, '[RCO]+', value_dict['pos_value'])
    # #                 fa_def_df.set_value(_i, 'NL-H2O', value_dict['m_water'])
    # #                 fa_def_df.set_value(_i, '[RCO+74]+', value_dict['mg_value'])
    # #             elif ion_charge in ['[M+Na]+']:
    # #                 value_dict = FA_list().pos_geo(elem_dct, ion_charge)
    # #                 fa_def_df.set_value(_i, '[RCO]+', value_dict['pos_value'])
    # #                 fa_def_df.set_value(_i, 'NL-H2O', value_dict['m_water'])
    # #                 fa_def_df.set_value(_i, '[RCO+74]+', value_dict['mg_value'])
    # #                 fa_def_df.set_value(_i, 'RCOONa', value_dict['m_Na'])
    # #
    # #             m_pre1_isotope_pattern_df = IsotopeHunter().get_isotope_mz(elem_dct, isotope_number=3)
    # #             print ('So is this the problem')
    # #             fa_def_df.set_value(_i, 'M+2_i_isomer', m_pre1_isotope_pattern_df.iloc[2]['ratio'])
    # #
    # #         self.fa_def_df = fa_def_df
    #         ##############################
    #         #
    #         #   Need to be removed
    #
    #         #
    #         #####################################
    #         print self.fa_def_df
    #
    #         self.weight_df = weight_df
    #         self.target_frag_df = key_frag_df.query(r'CLASS == "%s" and TYPE == "FRAG" and PR_CHARGE == "%s"'
    #                                                 % (lipid_type, ion_charge))
    #         self.target_nl_df = key_frag_df.query(r'CLASS == "%s" and TYPE == "NL" and PR_CHARGE == "%s"'
    #                                               % (lipid_type, ion_charge))
    #
    #         if ion_charge in ['[M-H]-', '[M+HCOO]-', '[M+CH3COO]-', '[M+FA-H]-', '[M+OAc]-']:
    #             charge_mode = 'NEG'
    #         elif ion_charge in ['[M+H]+', '[M+NH4]+', '[M+Na]+']:
    #             charge_mode = 'POS'
    #         else:
    #             charge_mode = 'NEG'
    #         self.charge_mode = charge_mode
    #         self.other_frag_df = key_frag_df.query('CLASS != "%s" and TYPE == "FRAG" and CHARGE_MODE == "%s"'
    #                                                % (lipid_type, charge_mode))
    #         self.other_nl_df = key_frag_df.query('CLASS != "%s" and TYPE == "NL" and CHARGE_MODE == "%s"'
    #                                              % (lipid_type, charge_mode))
    #         self.lipid_type = lipid_type
    #
    #     ############################################################################
    #     #   New way to predict the TG structures base on the white list.
    #     #   For now support only TG but can be change for PL also
    #     #   This new function creates all the possible structures from the beggining and it just find all the possible
    #     #   matches at once and it doesnt need to run for each precursor
    #     #   This way the speed of the program is update
    #     ###########################################################################
    #     def Propose_pre_str(self):
    #         lipid_abbr_lst = []
    #         lipid_sn1_lst = []
    #         lipid_sn2_lst = []
    #         lipid_sn3_lst = []
    #         c_sn1_lst = []
    #         c_sn2_lst = []
    #         c_sn3_lst = []
    #         db_sn1_lst = []
    #         db_sn2_lst = []
    #         db_sn3_lst = []
    #         c_total_lst = []
    #         db_total_lst = []
    #         link_total_lst = []
    #
    #
    #         if self.charge_mode == 'NEG':
    #             print('this is neg')
    #         else:
    #             fa_list = list(self.fa_def_df['FA'])
    #             proposed_str_lst = list(itertools.combinations_with_replacement(fa_list, 3))
    #             for _info in proposed_str_lst:
    #                 _pre_str = 'TG(' + _info[0] + '_' + _info[1] + '_' + _info[2] + ')'
    #                 lipid_abbr_lst.append(_pre_str)
    #                 lipid_sn1_lst.append(_info[0])
    #                 lipid_sn2_lst.append(_info[1])
    #                 lipid_sn3_lst.append(_info[2])
    #                 # It checks the link of the sn1
    #                 if _info[0][0:2] == 'O-':
    #                     _link_sn1 = 'O'
    #                 elif _info[0][0:2] == 'P-':
    #                     _link_sn1 = 'P'
    #                 else:
    #                     _link_sn1 = 'A'
    #                 # It checks the link of the sn2
    #                 if _info[1][0:2] == 'O-':
    #                     _link_sn2 = 'O'
    #                 elif _info[1][0:2] == 'P-':
    #                     _link_sn2 = 'P'
    #                 else:
    #                     _link_sn2 = 'A'
    #                 # It checks the link of the sn3
    #                 if _info[2][0:2] == 'O-':
    #                     _link_sn3 = 'O'
    #                 elif _info[2][0:2] == 'P-':
    #                     _link_sn3 = 'P'
    #                 else:
    #                     _link_sn3 = 'A'
    #                 # It checks the link type of the different FA so it will but first the O and P links first and then all the rest.
    #                 # This is done to create the TG link type and it will avoid to cause any probles in the later steps
    #                 # If the link was A-A-O- it want identify anything
    #                 if _link_sn1 == 'O' or _link_sn1 == 'P':
    #                     _total_link = _link_sn1 + '-'
    #                     if _link_sn2 == 'O' or _link_sn2 == 'P':
    #                         _total_link = _total_link + _link_sn2 + '-' + _link_sn3 + '-'
    #                     else:
    #                         _total_link = _total_link + _link_sn3 + '-' + _link_sn2 + '-'
    #                 elif _link_sn2 == 'O' or _link_sn2 == 'P':
    #                     _total_link = _link_sn2 + '-'
    #                     if _link_sn1 == 'O' or _link_sn1 == 'P':
    #                         _total_link = _total_link + _link_sn1 + '-' + _link_sn3 + '-'
    #                     else:
    #                         _total_link = _total_link + _link_sn3 + '-' + _link_sn1 + '-'
    #                 elif _link_sn3 == 'O' or _link_sn3 == 'P':
    #                     _total_link = _link_sn3 + '-'
    #                     if _link_sn1 == 'O' or _link_sn1 == 'P':
    #                         _total_link = _total_link + _link_sn1 + '-' + _link_sn2 + '-'
    #                     else:
    #                         _total_link = _total_link + _link_sn2 + '-' + _link_sn1 + '-'
    #                 else:
    #                     _total_link = _link_sn1 + '-' + _link_sn2 + '-' + _link_sn3 + '-'
    #
    #                 c_sn1_lst.append(_info[0].replace('O-', '').replace('P-', '')[0:2])
    #                 c_sn2_lst.append(_info[1].replace('O-', '').replace('P-', '')[0:2])
    #                 c_sn3_lst.append(_info[2].replace('O-', '').replace('P-', '')[0:2])
    #                 db_sn1_lst.append(_info[0].replace('O-', '').replace('P-', '')[3:5])
    #                 db_sn2_lst.append(_info[1].replace('O-', '').replace('P-', '')[3:5])
    #                 db_sn3_lst.append(_info[2].replace('O-', '').replace('P-', '')[3:5])
    #                 _sn1 = _info[0].replace('O-', '').replace('P-', '')
    #                 _sn2 = _info[1].replace('O-', '').replace('P-', '')
    #                 _sn3 = _info[2].replace('O-', '').replace('P-', '')
    #                 _sn1_lst = _sn1.split(':')
    #                 _sn2_lst = _sn2.split(':')
    #                 _sn3_lst = _sn3.split(':')
    #                 # _total_c = int(_sn1[0:2]) + int(_sn2[0:2]) + int(_sn3[0:2])
    #                 # _total_db = int(_sn1[3:5]) + int(_sn2[3:5]) + int(_sn3[3:5])
    #                 _total_c = int(_sn1_lst[0]) + int(_sn2_lst[0]) + int(_sn3_lst[0])
    #                 _total_db = int(_sn1_lst[1]) + int(_sn2_lst[1]) + int(_sn3_lst[1])
    #                 c_total_lst.append(_total_c)
    #                 db_total_lst.append(_total_db)
    #                 link_total_lst.append(_total_link)
    #
    #         lipid_abbr_df = pd.DataFrame(data={'Proposed_structures':lipid_abbr_lst, 'sn1_abbr':lipid_sn1_lst,
    #                                            'sn1_C': c_sn1_lst, 'sn1_DB': db_sn1_lst, 'sn2_abbr': lipid_sn2_lst,
    #                                            'sn2_C': c_sn2_lst, 'sn2_DB': db_sn2_lst, 'sn3_abbr': lipid_sn3_lst,
    #                                            'sn3_C': c_sn3_lst, 'sn3_DB': db_sn3_lst, 'total_C': c_total_lst,
    #                                            'total_DB': db_total_lst, 'link': link_total_lst})
    #         self.lipid_abbr_df = lipid_abbr_df
    #
    #     @staticmethod
    #     def get_pr_mz(charge_type, mz_lib):
    #
    #         pr_mz = 0.0
    #
    #         if charge_type in ['[M-H]-', '[M+HCOO]-', '[M+FA-H]-', '[M+CH3COO]-', '[M+OAc]-', '[M+AcOH-H]-']:
    #             charge_mode = 'NEG'
    #             if charge_type == '[M-H]-':
    #                 pr_mz = mz_lib
    #             elif charge_type in ['[M+HCOO]-', '[M+FA-H]-']:
    #                 pr_mz = mz_lib - 46.005480  # - HCOOH
    #             elif charge_type in ['[M+CH3COO]-', '[M+OAc]-']:
    #                 pr_mz = mz_lib - 60.021130  # - CH3COOH
    #
    #         elif charge_type in ['[M+H]+', '[M+Na]+', '[M+NH4]+', '[M+K]+']:
    #             charge_mode = 'POS'
    #             if charge_type == '[M+H]+':
    #                 pr_mz = mz_lib
    #             elif charge_type == '[M+Na]+':
    #                 pr_mz = mz_lib - 22.989770 + 1.007825  # - Na + H
    #             elif charge_type == '[M+NH4]+':
    #                 pr_mz = mz_lib - 17.026549  # - NH3
    #             elif charge_type == '[M+K]+':
    #                 pr_mz = mz_lib - 38.963708 + 1.007825  # - K + H
    #         else:
    #             charge_mode = 'NEG'
    #             pr_mz = mz_lib
    #
    #         return pr_mz, charge_mode
    #
    #     @staticmethod
    #     def decode_abbr(abbr):
    #         ###################################################################
    #         #
    #         #   The below regular expression ccan be compine in one
    #         #
    #         ##################################################################
    #         pl_checker = re.compile(r'(P[ACEGSI])([(])(.*)([)])')
    #         pip_checker = re.compile(r'(PIP)([(])(.*)([)])')
    #         tg_checker = re.compile(r'(TG)([(])(.*)([)])')
    #         fa_checker = re.compile(r'(\d{1,2})([:])(\d{1,2})')
    #         fa_o_checker = re.compile(r'(O-)(\d{1,2})([:])(\d)')
    #         fa_p_checker = re.compile(r'(P-)(\d{1,2})([:])(\d)')
    #
    #         # Check PL Type
    #         _pl_typ = ''
    #         bulk_fa_typ = ''
    #         bulk_fa_linker = ''
    #         bulk_fa_c = 0
    #         bulk_fa_db = 0
    #         lyso_fa_linker_dct = {'sn1': '', 'sn2': ''}
    #
    #         if pl_checker.match(abbr):
    #             print('PL')
    #             pl_re_chk = pl_checker.match(abbr)
    #             pl_typ_lst = pl_re_chk.groups()
    #             _pl_typ = pl_typ_lst[0]
    #             bulk_fa_typ = pl_typ_lst[2]
    #         if pip_checker.match(abbr):
    #             print('PIP')
    #             pip_re_chk = pip_checker.match(abbr)
    #             pip_typ_lst = pip_re_chk.groups()
    #             _pl_typ = pip_typ_lst[0]
    #             bulk_fa_typ = pip_typ_lst[2]
    #         if tg_checker.match(abbr):
    #             print('TG')
    #             tg_re_chk = tg_checker.match(abbr)
    #             tg_typ_lst = tg_re_chk.groups()
    #             _pl_typ = tg_typ_lst[0]
    #             bulk_fa_typ = tg_typ_lst[2]
    #         if fa_checker.match(abbr):
    #             print('FA')
    #             _pl_typ = 'FA'
    #             bulk_fa_typ = abbr
    #         if fa_o_checker.match(abbr):
    #             print('FA')
    #             _pl_typ = 'FA'
    #             bulk_fa_typ = abbr
    #         if fa_p_checker.match(abbr):
    #             print('FA')
    #             _pl_typ = 'FA'
    #             bulk_fa_typ = abbr
    #
    #         if fa_checker.match(bulk_fa_typ):
    #             if _pl_typ == "TG":
    #                 bulk_fa_linker = 'A-A-A-'
    #             else:
    #                 bulk_fa_linker = 'A-A-'
    #             lyso_fa_linker_dct = {'A': ''}
    #             fa_chk = fa_checker.match(bulk_fa_typ)
    #             bulk_fa_lst = fa_chk.groups()
    #             bulk_fa_c = bulk_fa_lst[0]
    #             bulk_fa_db = bulk_fa_lst[2]
    #         elif fa_o_checker.match(bulk_fa_typ):
    #             if _pl_typ == "TG":
    #                 bulk_fa_linker = 'O-A-A-'
    #             else:
    #                 bulk_fa_linker = 'O-A-'
    #             lyso_fa_linker_dct = {'O': '', 'A': 'O-'}  # link of the other sn after NL of this sn
    #             fa_chk = fa_o_checker.match(bulk_fa_typ)
    #             bulk_fa_lst = fa_chk.groups()
    #             bulk_fa_c = bulk_fa_lst[1]
    #             bulk_fa_db = bulk_fa_lst[3]
    #         elif fa_p_checker.match(bulk_fa_typ):
    #             if _pl_typ == "TG":
    #                 bulk_fa_linker = 'P-A-A-'
    #             else:
    #                 bulk_fa_linker = 'P-A-'
    #             lyso_fa_linker_dct = {'P': '', 'A': 'P-'}  # link of the other sn after NL of this sn
    #             fa_chk = fa_p_checker.match(bulk_fa_typ)
    #             bulk_fa_lst = fa_chk.groups()
    #             bulk_fa_c = bulk_fa_lst[1]
    #             bulk_fa_db = bulk_fa_lst[3]
    #
    #         bulk_fa_c = int(bulk_fa_c)
    #         bulk_fa_db = int(bulk_fa_db)
    #
    #         lipid_info_dct = {'TYPE': _pl_typ, 'LINK': bulk_fa_linker, 'C': bulk_fa_c, 'DB': bulk_fa_db,
    #                           'LYSO_LINK': lyso_fa_linker_dct}
    #         return lipid_info_dct
    #
    #     def get_fa_search(self, abbr, charge_type, mz_lib, ms2_df, ms2_precision=500e-6,
    #                       ms2_threshold=100, ms2_infopeak_threshold=0.02):
    #
    #         fa_ident_df = pd.DataFrame()
    #         lyso_ident_df = pd.DataFrame()      # For TG contains information about [TG - FA(RCO) -> [M+H]+ and [M+NH4]+] or [TG - FA(RCOOH)  -> [M+Na]+] => DG structures
    #         lyso_w_ident_df = pd.DataFrame()    # For TG contains information : [TG - FA (RCOOH) -> [M+H]+ and [M+NH4]+] or [TG - FA(RCOONa) -> [M+Na]+] => DG structures
    #         mg_w_ident_df = pd.DataFrame()      # only for TG structures: [TG - FA(RCO) - FA(RCOOH) -> [M+H]+ and [M+NH4]+] or [TG - FA(RCOOH) - FA(RCOONa) -> [M+Na]+]
    #
    #         lipid_info_dct = self.decode_abbr(abbr)
    #         pl_typ = lipid_info_dct['TYPE']
    #         bulk_fa_c = lipid_info_dct['C']
    #         bulk_fa_db = lipid_info_dct['DB']
    #         # bulk_fa_linker = lipid_info_dct['LINK']
    #         lyso_fa_linker_dct = lipid_info_dct['LYSO_LINK']
    #
    #         # _usr_formula_charged, usr_elem_charged_dct = ElemCalc().get_formula(abbr, charge=charge_type)
    #
    #         # use the max threshold from abs & relative intensity settings
    #         ms2_basepeak_i = ms2_df['i'].max()
    #         ms2_info_i = ms2_basepeak_i * ms2_infopeak_threshold
    #         ms2_threshold = max(ms2_threshold, ms2_info_i)
    #
    #         calc_pr_mz, charge_mode = self.get_pr_mz(charge_type, mz_lib)
    #         calc_pr_mz = mz_lib
    #
    #         if abbr[:2] in ['PA', 'PC', 'PE', 'PG', 'PI', 'PS', 'SM']:
    #             lipid_type = 'PL'
    #         elif abbr[:2] in ['TA', 'TG', 'DA', 'DG', 'MA', 'MG']:
    #             lipid_type = 'GL'
    #         else:
    #             lipid_type = 'PL'
    #
    #         if lipid_type == 'PL' and charge_mode == 'NEG':
    #             print ('negative')
    #             fa_chk_df = self.fa_def_df[['FA', 'Link', 'C', 'DB', 'mass', '[M-H]-', 'NL-H2O']]
    #             fa_chk_df = fa_chk_df.rename(columns={'[M-H]-': 'sn', 'mass': 'NL'})
    #
    #             if abbr[:2] == 'PC' and charge_type == '[M+HCOO]-':
    #                 fa_chk_df['[M-H]-sn'] = calc_pr_mz - fa_chk_df['NL-H2O'] - 60.021130  # - CH3COOH for PC
    #                 fa_chk_df['[M-H]-sn-H2O'] = calc_pr_mz - fa_chk_df['NL'] - 60.021130  # - CH3COOH for PC
    #                 fa_chk_df['Proposed_structures'] = ''
    #                 lyso_hg_mod = '-CH3'
    #             elif abbr[:2] == 'PC' and charge_type == '[M+OAc]-':
    #                 fa_chk_df['[M-H]-sn'] = calc_pr_mz - fa_chk_df['NL-H2O'] - 74.036780  # - CH3COOCH3 for PC
    #                 fa_chk_df['[M-H]-sn-H2O'] = calc_pr_mz - fa_chk_df['NL'] - 74.036780  # - CH3COOCH3 for PC
    #                 fa_chk_df['Proposed_structures'] = ''
    #                 lyso_hg_mod = '-CH3'
    #
    #             elif abbr[:2] == 'PS' and charge_type == '[M-H]-':
    #                 fa_chk_df['[M-H]-sn'] = calc_pr_mz - fa_chk_df['NL-H2O'] - 87.032029  # - C3H5NO2 for PS
    #                 fa_chk_df['[M-H]-sn-H2O'] = calc_pr_mz - fa_chk_df['NL'] - 87.032029  # - C3H5NO2 for PS
    #                 fa_chk_df['Proposed_structures'] = ''
    #                 lyso_hg_mod = '-87(Ser)'
    #
    #             else:
    #                 # Loss of FA-18, -OH remains on Glycerol back bone
    #                 fa_chk_df['[M-H]-sn'] = calc_pr_mz - fa_chk_df['NL-H2O']
    #                 # Loss of FA as full acid, -OH remains on FA NL
    #                 fa_chk_df['[M-H]-sn-H2O'] = calc_pr_mz - fa_chk_df['NL']
    #                 fa_chk_df['Proposed_structures'] = ''
    #                 lyso_hg_mod = ''
    #
    #             fa_abbr_lst = fa_chk_df['FA'].tolist()
    #
    #             for _i, _fa_se in fa_chk_df.iterrows():
    #
    #                 _fa_abbr = _fa_se['FA']
    #                 _fa_link = _fa_se['Link']
    #                 _fa_c = _fa_se['C']
    #                 _fa_db = _fa_se['DB']
    #
    #                 for _frag_type in ['sn', '[M-H]-sn', '[M-H]-sn-H2O']:
    #                     _frag_mz = _fa_se[_frag_type]
    #                     _frag_mz_low = _frag_mz - _frag_mz * ms2_precision
    #                     _frag_mz_high = _frag_mz + _frag_mz * ms2_precision
    #                     _frag_mz_query_code = '%f <= mz <= %f' % (_frag_mz_low, _frag_mz_high)
    #
    #                     _frag_df = ms2_df.query(_frag_mz_query_code)
    #
    #                     if _frag_df.shape[0] > 0:
    #                         _frag_df.loc[:, 'ppm'] = 1e6 * (_frag_df['mz'] - _frag_mz) / _frag_mz
    #                         _frag_df.loc[:, 'ppm_abs'] = _frag_df['ppm'].abs()
    #                         _frag_df.loc[:, 'FA'] = _fa_abbr
    #
    #                         if _frag_df.shape[0] > 1:
    #                             _frag_i_df = _frag_df.sort_values(by='i', ascending=False).head(1)
    #                             _frag_ppm_df = _frag_df.sort_values(by='ppm_abs').head(1)
    #                             _frag_df = _frag_i_df.copy()
    #                             if _frag_ppm_df['i'].tolist() == _frag_i_df['i'].tolist():
    #                                 pass
    #                             else:
    #                                 _frag_df = _frag_i_df.append(_frag_ppm_df)
    #
    #                         if _frag_type == 'sn':
    #                             if _fa_link != 'O' and _fa_link != 'P':
    #                                 _frag_df.loc[:, 'Proposed_structures'] = 'FA %s [M-H]-' % _fa_abbr
    #                                 fa_ident_df = fa_ident_df.append(_frag_df)
    #                         elif _frag_type == '[M-H]-sn':
    #                             if _fa_link in lyso_fa_linker_dct.keys():
    #                                 if bulk_fa_db - _fa_db >= 0:
    #                                     _fa_lyso_link = lyso_fa_linker_dct[_fa_link]
    #                                     _fa_lyso_str = '%s%i:%i' % (_fa_lyso_link, bulk_fa_c - _fa_c, bulk_fa_db - _fa_db)
    #                                     # keep theoretical common Lyso PL only
    #                                     if _fa_lyso_str in fa_abbr_lst:
    #                                         _frag_df.loc[:, 'Proposed_structures'] = ('L%s %s [M%s-H]-'
    #                                                                                   % (pl_typ, _fa_lyso_str, lyso_hg_mod)
    #                                                                                   )
    #                                         lyso_ident_df = lyso_ident_df.append(_frag_df)
    #                         elif _frag_type == '[M-H]-sn-H2O':
    #                             if _fa_link in lyso_fa_linker_dct.keys():
    #                                 if bulk_fa_db - _fa_db >= 0:
    #                                     _fa_lyso_link = lyso_fa_linker_dct[_fa_link]
    #                                     _fa_lyso_str = '%s%i:%i' % (_fa_lyso_link, bulk_fa_c - _fa_c, bulk_fa_db - _fa_db)
    #                                     # keep theoretical common Lyso PL only
    #                                     if _fa_lyso_str in fa_abbr_lst:
    #                                         _frag_df.loc[:, 'Proposed_structures'] = ('L%s %s [M%s-H2O-H]-'
    #                                                                                   % (pl_typ, _fa_lyso_str, lyso_hg_mod)
    #                                                                                   )
    #                                         lyso_w_ident_df = lyso_w_ident_df.append(_frag_df)
    #
    #         elif lipid_type == 'GL' and charge_mode == 'POS':
    #             print ('positive')
    #             if charge_type == '[M+NH4]+':
    #                 fa_chk_df = self.fa_def_df[['FA', 'Link', 'C', 'DB', 'mass', '[RCO]+', 'NL-H2O', '[RCO+74]+']]
    #                 fa_chk_df = fa_chk_df.rename(columns={'[RCO]+': 'sn', 'mass': 'NL', '[RCO+74]+': 'snGL'})
    #                 fa_chk_df['[M+H]-sn'] = calc_pr_mz - fa_chk_df['NL-H2O'] - 17.026549  # - NH3 adduct
    #                 fa_chk_df['[M+H]-sn-H2O'] = calc_pr_mz - fa_chk_df['NL'] - 17.026549  # - NH3 adduct
    #                 fa_chk_df['Proposed_structures'] = ''
    #                 lyso_hg_mod = '-NH3'
    #             elif charge_type == '[M+Na]+':
    #                 fa_chk_df = self.fa_def_df[['FA', 'Link', 'C', 'DB', 'mass', '[RCO]+', 'NL-H2O', 'RCOONa','[RCO+74]+']]
    #                 fa_chk_df = fa_chk_df.rename(columns={'[RCO]+': 'sn', 'mass': 'NL', '[RCO+74]+': 'snGL'})
    #                 fa_chk_df['[M+Na]-RCOOH'] = calc_pr_mz - fa_chk_df['NL']
    #                 fa_chk_df['[M+H]-RCOONa'] = calc_pr_mz - fa_chk_df['RCOONa']
    #                 fa_chk_df['Proposed_structures'] = ''
    #                 lyso_hg_mod = ''
    #             else:
    #                 fa_chk_df = self.fa_def_df[['FA', 'Link', 'C', 'DB', 'mass', '[RCO]+', 'NL-H2O', '[RCO+74]+']]
    #                 fa_chk_df = fa_chk_df.rename(columns={'[RCO]+': 'sn', 'mass': 'NL', '[RCO+74]+': 'snGL'})
    #                 fa_chk_df['[M+H]-sn'] = calc_pr_mz - fa_chk_df['NL-H2O']
    #                 fa_chk_df['[M+H]-sn-H2O'] = calc_pr_mz - fa_chk_df['NL']
    #                 fa_chk_df['Proposed_structures'] = ''
    #                 lyso_hg_mod = ''
    #
    #             fa_abbr_lst = fa_chk_df['FA'].tolist()
    #             print ('This is the calculate precursor')
    #             print (calc_pr_mz)
    #             for _i, _fa_se in fa_chk_df.iterrows():
    #
    #                 _fa_abbr = _fa_se['FA']
    #                 _fa_link = _fa_se['Link']
    #                 _fa_c = _fa_se['C']
    #                 _fa_db = _fa_se['DB']
    #                 if charge_type in ['[M+H]+', '[M+NH4]+']:
    #                     for _frag_type in ['sn', '[M+H]-sn', '[M+H]-sn-H2O', 'snGL']:
    #                         if _frag_type is not '[M+H]-2sn-H2O':
    #                             _frag_mz = _fa_se[_frag_type]
    #                             _frag_mz_low = _frag_mz - _frag_mz * ms2_precision
    #                             _frag_mz_high = _frag_mz + _frag_mz * ms2_precision
    #                             _frag_mz_query_code = '%f <= mz <= %f' % (_frag_mz_low, _frag_mz_high)
    #                             _frag_df = ms2_df.query(_frag_mz_query_code)
    #
    #                             if _frag_df.shape[0] > 0 :
    #
    #                                 _frag_df.loc[:, 'ppm'] = 1e6 * (_frag_df['mz'] - _frag_mz) / _frag_mz
    #                                 _frag_df.loc[:, 'ppm_abs'] = _frag_df['ppm'].abs()
    #                                 _frag_df.loc[:, 'FA'] = _fa_abbr
    #
    #                                 if _frag_df.shape[0] > 1:
    #
    #                                     _frag_i_df = _frag_df.sort_values(by='i', ascending=False).head(1)
    #                                     _frag_ppm_df = _frag_df.sort_values(by='ppm_abs').head(1)
    #                                     _frag_df = _frag_i_df.copy()
    #                                     if _frag_ppm_df['i'].tolist() == _frag_i_df['i'].tolist():
    #                                         pass
    #                                     else:
    #                                         _frag_df = _frag_i_df.append(_frag_ppm_df)
    #                                 if _frag_type == 'sn':
    #                                     if _fa_link != 'O' and _fa_link != 'P':
    #                                         _frag_df.loc[:, 'Proposed_structures'] = 'FA %s [RCO]+' % _fa_abbr
    #                                         fa_ident_df = fa_ident_df.append(_frag_df)
    #                                 elif _frag_type == '[M+H]-sn':
    #                                     if _fa_link in lyso_fa_linker_dct.keys():
    #                                         if bulk_fa_db - _fa_db >= 0:
    #                                             _fa_lyso_link = lyso_fa_linker_dct[_fa_link]
    #                                             _fa_lyso_str = '%s%i:%i' % (_fa_lyso_link, bulk_fa_c - _fa_c, bulk_fa_db - _fa_db)
    #                                             # keep theoretical common Lyso PL only
    #                                             #if _fa_lyso_str in fa_abbr_lst:
    #                                             _frag_df.loc[:, 'Proposed_structures'] = ('DG %s [M%s+H]+'
    #                                                                                           % (_fa_lyso_str, lyso_hg_mod)
    #                                                                                           )
    #                                             lyso_ident_df = lyso_ident_df.append(_frag_df)
    #                                 elif _frag_type == '[M+H]-sn-H2O':
    #
    #                                     if _fa_link in lyso_fa_linker_dct.keys():
    #                                         if bulk_fa_db - _fa_db >= 0:
    #                                             _fa_lyso_link = lyso_fa_linker_dct[_fa_link]
    #                                             _fa_lyso_str = '%s%i:%i' % (_fa_lyso_link, bulk_fa_c - _fa_c, bulk_fa_db - _fa_db)
    #                                             # keep theoretical common Lyso PL only
    #                                             #if _fa_lyso_str in fa_abbr_lst:
    #
    #                                             _frag_df.loc[:, 'Proposed_structures'] = ('DG %s [M%s-H2O+H]+'
    #                                                                                           % (_fa_lyso_str, lyso_hg_mod)
    #                                                                                           )
    #                                             lyso_w_ident_df = lyso_w_ident_df.append(_frag_df)
    #                                 elif _frag_type == 'snGL':
    #                                     if _fa_link in lyso_fa_linker_dct.keys():
    #                                         if bulk_fa_db - _fa_db >= 0:
    #                                             if _fa_link is not 'A':
    #                                                 _fa_lyso_link = _fa_link + '-'
    #                                             else:
    #                                                 _fa_lyso_link = ''
    #                                             _fa_lyso_str = '%s%i:%i' % (_fa_lyso_link, _fa_c, _fa_db)
    #                                             _frag_df.loc[:, 'Proposed_structures'] = ('MG %s [M%s+H]+' % (_fa_lyso_str, lyso_hg_mod))
    #                                             mg_w_ident_df = mg_w_ident_df.append(_frag_df)
    #                         elif _frag_type == '[M+H]-2sn-H2O':
    #                             for _i2, _fa_se2 in fa_chk_df.iterrows():
    #                                 _fa_abbr2 = _fa_se2['FA']
    #                                 _fa_link2 = _fa_se2['Link']
    #                                 _fa_c2 = _fa_se2['C']
    #                                 _fa_db2 = _fa_se2['DB']
    #                                 if not _fa_link == _fa_link2 == 'A':
    #                                     _frag_mz2 = _fa_se['[M+H]-sn'] + _fa_se2['NL']
    #
    #                                     _frag_mz_low2= _frag_mz2 - _frag_mz2*ms2_precision
    #                                     _frag_mz_high2 = _frag_mz2 + _frag_mz2*ms2_precision
    #                                     _frag_mz_query_code2='%f <= mz <= %f' % (_frag_mz_low2, _frag_mz_high2)
    #
    #                                     _frag_df = ms2_df.query(_frag_mz_query_code2)
    #
    #                                     if _frag_df.shape[0] > 0:
    #                                         _frag_df.loc[:, 'ppm'] = 1e6*(_frag_df['mz']-_frag_mz2)/_frag_mz2
    #                                         _frag_df.loc[:, 'ppm_abs'] = _frag_df['ppm'].abs()
    #                                         _frag_df.loc[:, 'FA'] = _fa_abbr
    #                                         _frag_df.loc[:, 'FA2'] = _fa_abbr2
    #
    #                                         if _frag_df.shape[0] > 1:
    #                                             # print (_frag_df)
    #                                             _frag_i_df = _frag_df.sort_values(by='i', ascending=False).head(1)
    #                                             # print (_frag_i_df)
    #                                             _frag_ppm_df = _frag_df.sort_values(by='ppm_abs').head(1)
    #                                             # print (_frag_ppm_df)
    #                                             _frag_df = _frag_i_df.copy()
    #                                             # if _frag_ppm_df['i'].tolist() == _frag_i_df['i'].tolist():
    #                                             #     pass
    #                                             # else:
    #                                             #     print ('4')
    #                                             #     _frag_df = _frag_i_df.append(_frag_ppm_df)
    #                                             if _fa_link in lyso_fa_linker_dct.keys() and _fa_link2 in lyso_fa_linker_dct.keys():
    #                                                 if bulk_fa_db - _fa_db - _fa_db2 >= 0:
    #                                                     if _fa_link == _fa_link2  and  _fa_link == 'A':
    #                                                         _fa_mg_link = 'P-'
    #                                                     else:
    #                                                         _fa_mg_link = ''
    #                                                     _fa_mg_str = '%s%i:%i' % (_fa_mg_link, bulk_fa_c - _fa_c - _fa_c2, bulk_fa_db - _fa_db - _fa_db2)
    #                                                     if _fa_mg_str in fa_abbr_lst:
    #                                                         # print ('So I guess somewhere around')
    #                                                         _frag_df.loc[:, 'Proposed_structures'] = ('MG%s %s [M%s+H]+' % (pl_typ, _fa_mg_str, lyso_hg_mod))
    #                                                         mg_w_ident_df = mg_w_ident_df.append(_frag_df)
    #                                                         print(mg_w_ident_df)
    #                 elif charge_type in ['[M+Na]+']:
    #                     for _frag_type in ['sn', '[M+H]-RCOONa', '[M+Na]-RCOOH', 'snGL']:
    #                         if _frag_type is not '[M+H]-RCOOH-RCOONa':
    #                             _frag_mz = _fa_se[_frag_type]
    #                             _frag_mz_low = _frag_mz - _frag_mz * ms2_precision
    #                             _frag_mz_high = _frag_mz + _frag_mz * ms2_precision
    #                             _frag_mz_query_code = '%f <= mz <= %f' % (_frag_mz_low, _frag_mz_high)
    #                             _frag_df = ms2_df.query(_frag_mz_query_code)
    #
    #                             if _frag_df.shape[0] > 0:
    #
    #                                 _frag_df.loc[:, 'ppm'] = 1e6 * (_frag_df['mz'] - _frag_mz) / _frag_mz
    #                                 _frag_df.loc[:, 'ppm_abs'] = _frag_df['ppm'].abs()
    #                                 _frag_df.loc[:, 'FA'] = _fa_abbr
    #
    #                                 if _frag_df.shape[0] > 1:
    #
    #                                     _frag_i_df = _frag_df.sort_values(by='i', ascending=False).head(1)
    #                                     _frag_ppm_df = _frag_df.sort_values(by='ppm_abs').head(1)
    #                                     _frag_df = _frag_i_df.copy()
    #                                     if _frag_ppm_df['i'].tolist() == _frag_i_df['i'].tolist():
    #                                         pass
    #                                     else:
    #                                         _frag_df = _frag_i_df.append(_frag_ppm_df)
    #                                 if _frag_type == 'sn':
    #                                     if _fa_link != 'O' and _fa_link != 'P':
    #                                         _frag_df.loc[:, 'Proposed_structures'] = 'FA %s [RCONL+Na]+' % _fa_abbr
    #                                         fa_ident_df = fa_ident_df.append(_frag_df)
    #                                 elif _frag_type == '[M+H]-RCOONa':
    #                                     if _fa_link in lyso_fa_linker_dct.keys():
    #                                         if bulk_fa_db - _fa_db >= 0:
    #                                             _fa_lyso_link = lyso_fa_linker_dct[_fa_link]
    #                                             _fa_lyso_str = '%s%i:%i' % (
    #                                             _fa_lyso_link, bulk_fa_c - _fa_c, bulk_fa_db - _fa_db)
    #                                             # keep theoretical common Lyso PL only
    #                                             # if _fa_lyso_str in fa_abbr_lst:
    #                                             _frag_df.loc[:, 'Proposed_structures'] = ('DG %s [M%s+H-RCOONa]+'
    #                                                                                       % (_fa_lyso_str, lyso_hg_mod)
    #                                                                                       )
    #                                             lyso_ident_df = lyso_ident_df.append(_frag_df)
    #                                 elif _frag_type == '[M+Na]-RCOOH':
    #
    #                                     if _fa_link in lyso_fa_linker_dct.keys():
    #                                         if bulk_fa_db - _fa_db >= 0:
    #                                             _fa_lyso_link = lyso_fa_linker_dct[_fa_link]
    #                                             _fa_lyso_str = '%s%i:%i' % (
    #                                             _fa_lyso_link, bulk_fa_c - _fa_c, bulk_fa_db - _fa_db)
    #                                             # keep theoretical common Lyso PL only
    #                                             # if _fa_lyso_str in fa_abbr_lst:
    #
    #                                             _frag_df.loc[:, 'Proposed_structures'] = ('DG %s [M%s+Na-RCOOH]+'
    #                                                                                       % (_fa_lyso_str, lyso_hg_mod)
    #                                                                                       )
    #                                             lyso_w_ident_df = lyso_w_ident_df.append(_frag_df)
    #                                 elif _frag_type == 'snGL':
    #                                     if _fa_link in lyso_fa_linker_dct.keys():
    #                                         if bulk_fa_db - _fa_db >= 0:
    #                                             if _fa_link is not 'A':
    #                                                 _fa_lyso_link = _fa_link + '-'
    #                                             else:
    #                                                 _fa_lyso_link = ''
    #                                             _fa_lyso_link = _fa_link + '-'
    #                                             _fa_lyso_str = '%s%i:%i' % (_fa_lyso_link, _fa_c, _fa_db)
    #                                             _frag_df.loc[:, 'Proposed_structures'] = ('MG %s [M%s+H]+' % (_fa_lyso_str, lyso_hg_mod))
    #                                             mg_w_ident_df = mg_w_ident_df.append(_frag_df)
    #
    #                         elif _frag_type == '[M+H]-RCOOH-RCOONa':
    #                             for _i2, _fa_se2 in fa_chk_df.iterrows():
    #                                 _fa_abbr2 = _fa_se2['FA']
    #                                 _fa_link2 = _fa_se2['Link']
    #                                 _fa_c2 = _fa_se2['C']
    #                                 _fa_db2 = _fa_se2['DB']
    #                                 if not _fa_link == _fa_link2 == 'A':
    #                                     _frag_mz2 = _fa_se['[M+H]-RCOONa'] + _fa_se2['NL']
    #
    #                                     _frag_mz_low2 = _frag_mz2 - _frag_mz2 * ms2_precision
    #                                     _frag_mz_high2 = _frag_mz2 + _frag_mz2 * ms2_precision
    #                                     _frag_mz_query_code2 = '%f <= mz <= %f' % (_frag_mz_low2, _frag_mz_high2)
    #
    #                                     _frag_df = ms2_df.query(_frag_mz_query_code2)
    #
    #                                     if _frag_df.shape[0] > 0:
    #                                         _frag_df.loc[:, 'ppm'] = 1e6 * (_frag_df['mz'] - _frag_mz2) / _frag_mz2
    #                                         _frag_df.loc[:, 'ppm_abs'] = _frag_df['ppm'].abs()
    #                                         _frag_df.loc[:, 'FA'] = _fa_abbr
    #                                         _frag_df.loc[:, 'FA2'] = _fa_abbr2
    #
    #                                         if _frag_df.shape[0] > 1:
    #                                             # print (_frag_df)
    #                                             _frag_i_df = _frag_df.sort_values(by='i', ascending=False).head(1)
    #                                             # print (_frag_i_df)
    #                                             _frag_ppm_df = _frag_df.sort_values(by='ppm_abs').head(1)
    #                                             # print (_frag_ppm_df)
    #                                             _frag_df = _frag_i_df.copy()
    #                                             # if _frag_ppm_df['i'].tolist() == _frag_i_df['i'].tolist():
    #                                             #     pass
    #                                             # else:
    #                                             #     print ('4')
    #                                             #     _frag_df = _frag_i_df.append(_frag_ppm_df)
    #                                             if _fa_link in lyso_fa_linker_dct.keys() and _fa_link2 in lyso_fa_linker_dct.keys():
    #                                                 if bulk_fa_db - _fa_db - _fa_db2 >= 0:
    #                                                     if _fa_link == _fa_link2 and _fa_link == 'A':
    #                                                         _fa_mg_link = 'P-'
    #                                                     else:
    #                                                         _fa_mg_link = ''
    #                                                     _fa_mg_str = '%s%i:%i' % (_fa_mg_link, bulk_fa_c - _fa_c - _fa_c2,
    #                                                                               bulk_fa_db - _fa_db - _fa_db2)
    #                                                     if _fa_mg_str in fa_abbr_lst:
    #                                                         # print ('So I guess somewhere around')
    #                                                         _frag_df.loc[:, 'Proposed_structures'] = (
    #                                                         'MG%s %s [M%s+H]+' % (pl_typ, _fa_mg_str, lyso_hg_mod))
    #                                                         mg_w_ident_df = mg_w_ident_df.append(_frag_df)
    #                                                         print(mg_w_ident_df)
    #
    #         # format the output DataFrame
    #         if fa_ident_df.shape[0] > 0:
    #             fa_ident_df = fa_ident_df.query('i > %f' % ms2_threshold)
    #             # print ('FA identification')
    #             # print fa_ident_df
    #             # exit()
    #             proposed_str_lst = {}
    #             if lipid_type == 'GL':
    #                 # exit()
    #                 fa_ident_df['Flag'] = 1
    #                 fa_ident_df['Type']='FA'
    #                 fa_ident_df = fa_ident_df[['Proposed_structures', 'FA', 'mz', 'i',
    #                                            'ppm', 'ppm_abs', 'Flag']].reset_index(drop=True)
    #             else:
    #                 fa_ident_df['Flag'] = 1
    #                 fa_ident_df = fa_ident_df[['Proposed_structures', 'FA', 'mz', 'i',
    #                                            'ppm', 'ppm_abs', 'Flag']].reset_index(drop=True)
    #
    #             fa_ident_df = fa_ident_df.sort_values(by=['i', 'ppm_abs'], ascending=[False, True])
    #             fa_ident_df = fa_ident_df.drop_duplicates(['FA'], keep='first')
    #             fa_ident_df = fa_ident_df.sort_values(by='i', ascending=False).head(10)
    #         if lyso_ident_df.shape[0] > 0:
    #             lyso_found_dct = {}
    #             lyso_ident_df = lyso_ident_df.query('i > %f' % ms2_threshold)
    #             if lipid_type == 'GL':
    #                 lyso_ident_df['Flag'] = 1
    #                 lyso_ident_df['Type'] = 'Lyso'
    #                 lyso_ident_df = lyso_ident_df.loc[lyso_ident_df['Flag'] == 1][['Proposed_structures', 'FA', 'mz', 'i',
    #                                                                                'ppm', 'ppm_abs',
    #                                                                                'Flag']].reset_index(drop=True)
    #
    #             else:
    #                 lyso_ident_df['Flag'] = 1
    #                 lyso_ident_df = lyso_ident_df.loc[lyso_ident_df['Flag'] == 1][['Proposed_structures', 'FA', 'mz', 'i',
    #                                                                                'ppm', 'ppm_abs',
    #                                                                                'Flag']].reset_index(drop=True)
    #             lyso_ident_df = lyso_ident_df.sort_values(by=['i', 'ppm_abs'], ascending=[False, True])
    #             lyso_ident_df = lyso_ident_df.drop_duplicates(['FA'], keep='first')
    #             lyso_ident_df = lyso_ident_df.sort_values(by='i', ascending=False).head(5)
    #
    #         if lyso_w_ident_df.shape[0] > 0:
    #             lyso_w_dct = {}
    #             lyso_w_ident_df = lyso_w_ident_df.query('i > %f' % ms2_threshold).reset_index(drop=True)
    #             if lipid_type == 'GL':
    #                 lyso_w_ident_df['Flag'] = 1
    #                 lyso_w_ident_df['Type'] = 'Lyso_W'
    #                 lyso_w_ident_df = lyso_w_ident_df[['Proposed_structures', 'FA', 'mz', 'i',
    #                                                    'ppm', 'ppm_abs', 'Flag']].reset_index(drop=True)
    #             else:
    #                 lyso_w_ident_df['Flag'] = 1
    #                 lyso_w_ident_df = lyso_w_ident_df[['Proposed_structures', 'FA', 'mz', 'i',
    #                                                    'ppm', 'ppm_abs', 'Flag']].reset_index(drop=True)
    #             lyso_w_ident_df = lyso_w_ident_df.sort_values(by=['i', 'ppm_abs'], ascending=[False, True])
    #             lyso_w_ident_df = lyso_w_ident_df.drop_duplicates(['FA'], keep='first')
    #             lyso_w_ident_df = lyso_w_ident_df.sort_values(by='i', ascending=False).head(5)
    #         if mg_w_ident_df.shape[0] > 0:
    #             mg_w_ident_df = mg_w_ident_df.query('i > %f' % ms2_threshold)
    #             if lipid_type == 'GL':
    #                 mg_w_ident_df['Flag'] = 1
    #                 mg_w_ident_df['Type'] = 'MG'
    #                 mg_w_ident_df = mg_w_ident_df.loc[mg_w_ident_df['Flag'] == 1][['Proposed_structures', 'FA', 'mz', 'i', 'ppm', 'ppm_abs', 'Flag']].reset_index(drop = True)
    #             else:
    #                 mg_w_ident_df['Flag'] = 1
    #                 mg_w_ident_df =mg_w_ident_df.loc[mg_w_ident_df['Flag'] == 1][['Proposed_structures', 'FA', 'mz', 'i', 'ppm', 'ppm_abs', 'Flag']].reset_index(drop = True)
    #             mg_w_ident_df = mg_w_ident_df.sort_values(by=['i', 'ppm_abs'], ascending=[False, True])
    #             mg_w_ident_df = mg_w_ident_df.drop_duplicates(['FA'], keep='first')
    #             mg_w_ident_df = mg_w_ident_df.sort_values(by='i', ascending=False).head(5)
    #
    #         return fa_ident_df, lyso_ident_df, lyso_w_ident_df, mg_w_ident_df
    #
    #     def get_structure(self, abbr):
    #
    #         lipid_abbr_lst = []
    #         lipid_sn1_lst = []
    #         lipid_sn2_lst = []
    #         lipid_sn3_lst = []
    #         db_sn1_lst = []
    #         db_sn2_lst = []
    #         db_sn3_lst = []
    #         # abbr='TG(46:3)'
    #         # abbr = 'TG(O-58:8)'
    #         print(abbr)
    #         lipid_info_dct = self.decode_abbr(abbr)
    #         pl_typ = lipid_info_dct['TYPE']
    #         bulk_fa_c = lipid_info_dct['C']
    #         bulk_fa_db = lipid_info_dct['DB']
    #         bulk_fa_linker = lipid_info_dct['LINK']
    #
    #         #######################################
    #         #   Also can change like the TG prediction
    #         #######################################
    #         if abbr[:2] in ['PE', 'PA', 'PC', 'PI', 'PS', 'PG']:
    #             for _i, _fa_se in self.fa_def_df.iterrows():
    #                 # FA, Link, C, DB
    #                 _fa_abbr = _fa_se['FA']
    #                 _fa_link = _fa_se['Link']
    #                 _fa_c = _fa_se['C']
    #                 _fa_db = _fa_se['DB']
    #
    #                 if _fa_db <= bulk_fa_db and _fa_c <= bulk_fa_c:
    #                     if _fa_link == bulk_fa_linker[0:1]:
    #                         _rest_fa_link = bulk_fa_linker[2]
    #                         _rest_fa_c = bulk_fa_c - _fa_c
    #                         _rest_fa_db = bulk_fa_db - _fa_db
    #
    #                         _rest_fa_df = self.fa_def_df.query('Link == "%s" and C == %i and DB == %i'
    #                                                            % (_rest_fa_link, _rest_fa_c, _rest_fa_db)
    #                                                            )
    #
    #                         if _rest_fa_df.shape[0] == 1:
    #                             _rest_fa_abbr = _rest_fa_df['FA'].tolist()[0]
    #                             lipid_abbr = '%s(%s_%s)' % (pl_typ, _fa_abbr, _rest_fa_abbr)
    #
    #                             lipid_abbr_lst.append(lipid_abbr)
    #                             lipid_sn1_lst.append(_fa_abbr)
    #                             lipid_sn2_lst.append(_rest_fa_abbr)
    #                             db_sn1_lst.append(_fa_db)
    #                             db_sn2_lst.append(_rest_fa_db)
    #             lipid_abbr_df = pd.DataFrame(data={'Proposed_structures': lipid_abbr_lst, 'sn1_abbr': lipid_sn1_lst,
    #                                                'sn2_abbr': lipid_sn2_lst, 'sn1_DB': db_sn1_lst, 'sn2_DB': db_sn2_lst})
    #
    #             lipid_abbr_df = lipid_abbr_df.query('sn1_DB <=sn2_DB')
    #             lipid_abbr_df = lipid_abbr_df[['Proposed_structures', 'sn1_abbr', 'sn2_abbr']]
    #             # print lipid_abbr_df
    #             # exit()
    #
    #             return lipid_abbr_df
    #
    #         elif abbr[:2] in ['TG'] and bulk_fa_linker not in ['O', 'P']:
    #             allsnList=[]
    #             final_info_df = self.lipid_abbr_df.query('link == "%s" and total_C == %i and total_DB == %i' % (bulk_fa_linker, int(bulk_fa_c), int(bulk_fa_db)))
    #             lipid_abbr_df = final_info_df[['Proposed_structures', 'sn1_abbr', 'sn2_abbr', 'sn3_abbr']]
    #
    #             return lipid_abbr_df
    #         else:
    #
    #             lipid_sn3_lst = []
    #             db_sn3_lst = []
    #             # _fa_compination_3 = []
    #
    #             lipid_abbr_df = pd.DataFrame(data={'Proposed_structures': lipid_abbr_lst, 'sn1_abbr': lipid_sn1_lst,
    #                                                'sn2_abbr': lipid_sn2_lst, 'sn3_abbr': lipid_sn3_lst,
    #                                                'sn1_DB': db_sn1_lst,
    #                                                'sn2_DB': db_sn2_lst, 'sn3_DB': db_sn3_lst})
    #
    #             lipid_abbr_df = lipid_abbr_df[['Proposed_structures', 'sn1_abbr', 'sn2_abbr', 'sn3_abbr']]
    #
    #             return lipid_abbr_df
    #
    #     def get_match(self, abbr, charge_type, mz_lib, ms2_df, ms2_precision=500e-6,
    #                   ms2_threshold=100, ms2_infopeak_threshold=0.02, rank_mode=True):
    #
    #         match_reporter = 0
    #         ms2_max_i = ms2_df['i'].max()
    #         formula, formula_dct = ElemCalc().get_formula(abbr, charge=charge_type)
    #
    #         fa_ident_df, lyso_ident_df, lyso_w_ident_df, mg_w_ident_df = self.get_fa_search(abbr, charge_type, mz_lib, ms2_df,
    #                                                                          ms2_precision=ms2_precision,
    #                                                                          ms2_threshold=ms2_threshold,
    #                                                                          ms2_infopeak_threshold=ms2_infopeak_threshold
    #                                                                          )
    #         lipid_abbr_df = self.get_structure(abbr)
    #         sodiumFlag=0
    #         if abbr[:2] in ['TG']:
    #             # print "Kipors"
    #             # exit()
    #             if charge_type == '[M+Na]+':
    #                 sodiumFlag=1
    #                 weight_type_lst = ['sn1', 'sn2', 'sn3',  '[M+H]-sn1', '[M+H]-sn2', '[M+H]-sn3',
    #                                '[M+Na]-sn1', '[M+Na]-sn2', '[M+Na]-sn3', '[M+H]-(sn1+sn2)-H2O', '[M+H]-(sn1+sn3)-H2O', '[M+H]-(sn2+sn3)-H2O']
    #             else:
    #                 weight_type_lst = ['sn1', 'sn2', 'sn3',  '[M+H]-sn1', '[M+H]-sn2', '[M+H]-sn3',
    #                                '[M+H]-sn1-H2O', '[M+H]-sn2-H2O', '[M+H]-sn3-H2O', '[M+H]-(sn1+sn2)-H2O', '[M+H]-(sn1+sn3)-H2O', '[M+H]-(sn2+sn3)-H2O']
    #         else:
    #             weight_type_lst = ['sn1', 'sn2', '[M-H]-sn1', '[M-H]-sn2',
    #                                '[M-H]-sn1-H2O', '[M-H]-sn2-H2O']
    #         weight_dct = {}
    #         for _type in weight_type_lst:
    #             lipid_abbr_df[_type] = 0
    #         if fa_ident_df.shape[0] > 0 or lyso_ident_df.shape[0] > 0 or lyso_w_ident_df.shape[0] > 0 or mg_w_ident_df.shape[0] > 0:
    #             combine_all_lst = pd.DataFrame()
    #             try:
    #                 fa_ident_df['Type'] = 'FA'
    #                 fa_ident_lst = fa_ident_df.loc[fa_ident_df['Flag'] == 1]['FA'].tolist()
    #                 fa_i_lst = fa_ident_df.loc[fa_ident_df['Flag'] == 1]['i'].tolist()
    #                 combine_all_lst = combine_all_lst.append(fa_ident_df.loc[fa_ident_df['Flag'] == 1], ignore_index=True)
    #             except KeyError:
    #                 fa_ident_lst = []
    #                 fa_i_lst = []
    #
    #             try:
    #                 lyso_ident_df['Type'] = 'Lyso'
    #                 lyso_ident_lst = lyso_ident_df.loc[lyso_ident_df['Flag'] == 1]['FA'].tolist()
    #                 lyso_i_lst = lyso_ident_df.loc[lyso_ident_df['Flag'] == 1]['i'].tolist()
    #                 combine_all_lst = combine_all_lst.append(lyso_ident_df.loc[lyso_ident_df['Flag'] == 1], ignore_index=True)
    #             except KeyError:
    #                 lyso_ident_lst = []
    #                 lyso_i_lst = []
    #             # print lyso_ident_df
    #             try:
    #                 lyso_w_ident_df['Type'] = 'LysoW'
    #                 lyso_w_ident_lst = lyso_w_ident_df.loc[lyso_w_ident_df['Flag'] == 1]['FA'].tolist()
    #                 lyso_w_i_lst = lyso_w_ident_df.loc[lyso_w_ident_df['Flag'] == 1]['i'].tolist()
    #                 combine_all_lst = combine_all_lst.append(lyso_w_ident_df.loc[lyso_w_ident_df['Flag'] == 1], ignore_index=True)
    #             except KeyError:
    #                 lyso_w_ident_lst = []
    #                 lyso_w_i_lst = []
    #             try:
    #                 mg_w_ident_df['Type'] = 'MG'
    #                 mg_w_ident_lst = mg_w_ident_df.loc[mg_w_ident_df['Flag'] == 1][['FA']]
    #                 mg_w_i_lst = mg_w_ident_df.loc[mg_w_ident_df['Flag'] == 1]['i'].tolist()
    #                 combine_all_lst = combine_all_lst.append(mg_w_ident_df.loc[mg_w_ident_df['Flag'] == 1], ignore_index=True)
    #             except KeyError:
    #                 mg_w_ident_lst = ()
    #                 mg_w_i_lst = ()
    #
    #             combine_all_lst = combine_all_lst.sort_values(by=['mz', 'i', 'ppm_abs'], ascending=[True, False, True])
    #             combine_all_lst = combine_all_lst.drop_duplicates(['mz'], keep = 'first')
    #             print(combine_all_lst)
    #
    #             ####################################
    #             #
    #             #   Maybe in a statement for the differences in the different modes
    #             #
    #             for _i_comb, _row_comb in combine_all_lst.iterrows():
    #                 _mz_low_combine = _row_comb['mz'] - _row_comb['mz']*ms2_precision
    #                 _mz_high_combine = _row_comb['mz'] + _row_comb['mz']*ms2_precision
    #                 _query_comb = '%f <= mz <= %f' % (_mz_low_combine, _mz_high_combine)
    #                 combine_small = combine_all_lst.query(_query_comb)
    #                 combine_small_index_list = combine_small.index.tolist()
    #                 if combine_small.shape[0] > 1:
    #                     combine_small = combine_small.sort_values(by='i', ascending=False).head(1)
    #                     small_combine_index = combine_small.index.tolist()
    #                     for _i_pos in combine_small_index_list:
    #                         if _i_pos not in small_combine_index:
    #                             combine_all_lst = combine_all_lst.drop([_i_pos])
    #
    #                 _query_comb2 = '%f <= mz <= %f' % (_row_comb['mz'] - 2.5, _row_comb['mz'] + 2.5)
    #                 combine_small2 = combine_all_lst.query(_query_comb2)
    #                 combine_small2 = combine_small2.sort_values(by='mz', ascending=True)
    #                 if combine_small2.shape[0] <3 and combine_small2.shape[0] > 1:
    #                     base_m1_i = combine_small2.iloc[0]['i']
    #                     base_m2_i = combine_small2.iloc[1]['i']
    #                     formula_dg_dct = formula_dct.copy()
    #                     # formula_dg_dct_dha = formula_dct.copy()
    #                     # formula_dg_dct_pa = formula_dct.copy()
    #                     if base_m1_i > base_m2_i:
    #                         if combine_small2.iloc[0]['Type'] == 'LysoW':
    #                             fa_m1_i = combine_small2.iloc[0]['FA']
    #                             fa_elem_mass = self.fa_def_df.loc[self.fa_def_df['FA'] == str(fa_m1_i)]['elem']
    #                             fa_formula_dct = IsotopeHunter().get_elements(fa_elem_mass.iloc[0])
    #                             formula_dg=''
    #                             for k in formula_dct.keys():
    #                                 if formula_dct[k] == 0:
    #                                     pass
    #                                 else:
    #                                     if k in fa_formula_dct.keys():
    #                                         formula_dg_dct[k]= formula_dct[k] - fa_formula_dct[k]
    #                                     else:
    #                                         formula_dg_dct[k] = formula_dct[k]
    #
    #                             for k in ['C', 'H', 'N', 'O', 'Na']:
    #                                 if k in formula_dg_dct.keys():
    #                                     if formula_dg_dct[k] == 0:
    #                                         pass
    #                                     else:
    #                                         formula_dg = formula_dg + k + str(formula_dg_dct[k])
    #
    #                             fa_elem_mass_dha = self.fa_def_df.loc[self.fa_def_df['FA'] == '22:06']['elem']
    #                             fa_formula_dct_dha = IsotopeHunter().get_elements(fa_elem_mass_dha.iloc[0])
    #                             formula_dg_dct_dha={}
    #                             for k in formula_dct.keys():
    #                                 if formula_dct[k] == 0:
    #                                     pass
    #                                 else:
    #                                     if k in fa_formula_dct_dha.keys():
    #                                         formula_dg_dct_dha[k]= formula_dct[k] - fa_formula_dct_dha[k]
    #                                     else:
    #                                         formula_dg_dct_dha[k] = formula_dct[k]
    #                             isotope_dha = IsotopeHunter().get_isotope_mz(formula_dg_dct_dha, isotope_number=2)
    #
    #                             fa_elem_mass_pa = self.fa_def_df.loc[self.fa_def_df['FA'] == '14:00']['elem']
    #                             fa_formula_dct_pa = IsotopeHunter().get_elements(fa_elem_mass_pa.iloc[0])
    #                             formula_dg_dct_pa={}
    #                             for k in formula_dct.keys():
    #                                 if formula_dct[k] == 0:
    #                                     pass
    #                                 else:
    #                                     if k in fa_formula_dct_pa.keys():
    #                                         formula_dg_dct_pa[k] = formula_dct[k] - fa_formula_dct_pa[k]
    #                                     else:
    #                                         formula_dg_dct_pa[k] = formula_dct[k]
    #                             isotope_pa = IsotopeHunter().get_isotope_mz(formula_dg_dct_pa, isotope_number=2)
    #
    #                             pra_m2_i_ratio = ((100 * base_m2_i)/base_m1_i)/100
    #                             if pra_m2_i_ratio <= isotope_pa.iloc[2]['ratio'] and pra_m2_i_ratio >= isotope_dha.iloc[2]['ratio']:
    #                                 isotope_flag = IsotopeHunter().get_isotope_fragments(combine_small2.iloc[1]['mz'],
    #                                                                                   combine_small2.iloc[1]['i'], formula_dg,
    #                                                                                   ms2_df, isotope_number=2,
    #                                                                                   ms1_precision=ms2_precision)
    #                                 if isotope_flag == 1 :
    #                                     combine_small2_index = combine_small2.index.tolist()
    #                                     combine_all_lst = combine_all_lst.drop([combine_small2_index[1]])
    #                             else:
    #                                 print ('Nop an isotope')
    #
    #                         elif combine_small2.iloc[0]['Type'] == 'MG':
    #                             fa_m1_i = combine_small2.iloc[0]['FA']
    #                             fa_elem_mass = self.fa_def_df.loc[self.fa_def_df['FA'] == str(fa_m1_i)]['elem']
    #                             fa_formula_dct = IsotopeHunter().get_elements(fa_elem_mass.iloc[0])
    #                             fa_formula_dct['C'] = int(fa_formula_dct['C']) + int(3)
    #                             fa_formula_dct['H'] = int(fa_formula_dct['H']) + int(5)
    #                             fa_formula_dct['O'] = int(fa_formula_dct['O']) + int(1)
    #
    #                             fa_elem_mass_dha = self.fa_def_df.loc[self.fa_def_df['FA'] == '22:06']['elem']
    #                             fa_formula_dct_dha = IsotopeHunter().get_elements(fa_elem_mass_dha.iloc[0])
    #                             fa_formula_dct_dha['C'] = int(fa_formula_dct_dha['C']) + int(3)
    #                             fa_formula_dct_dha['H'] = int(fa_formula_dct_dha['H']) + int(5)
    #                             fa_formula_dct_dha['O'] = int(fa_formula_dct_dha['O']) + int(1)
    #                             isotope_dha = IsotopeHunter().get_isotope_mz(fa_formula_dct_dha, isotope_number=2)
    #
    #                             fa_elem_mass_pa = self.fa_def_df.loc[self.fa_def_df['FA'] == '14:00']['elem']
    #                             fa_formula_dct_pa = IsotopeHunter().get_elements(fa_elem_mass_pa.iloc[0])
    #                             fa_formula_dct_pa['C'] = int(fa_formula_dct_pa['C']) + int(3)
    #                             fa_formula_dct_pa['H'] = int(fa_formula_dct_pa['H']) + int(5)
    #                             fa_formula_dct_pa['O'] = int(fa_formula_dct_pa['O']) + int(1)
    #                             isotope_pa = IsotopeHunter().get_isotope_mz(fa_formula_dct_pa, isotope_number=2)
    #
    #                             pra_m2_i_ratio = ((100 * base_m2_i) / base_m1_i)/100
    #                             if pra_m2_i_ratio <= isotope_pa.iloc[2]['ratio'] and pra_m2_i_ratio >= isotope_dha.iloc[2]['ratio']:
    #                                 print('Probably is an isotope')
    #                                 print(combine_small2)
    #                         elif combine_small2.iloc[0]['Type'] == 'Lyso':
    #                             fa_m1_i = combine_small2.iloc[0]['FA']
    #                             fa_elem_mass = self.fa_def_df.loc[self.fa_def_df['FA'] == str(fa_m1_i)]['elem']
    #                             fa_formula_dct = IsotopeHunter().get_elements(fa_elem_mass.iloc[0])
    #                             fa_formula_dct['H'] = fa_formula_dct['H'] - 2
    #                             fa_formula_dct['O'] = fa_formula_dct['O'] -1
    #                             formula_dg = ''
    #                             formula_dg_dct = {}
    #                             for k in formula_dct.keys():
    #                                 if formula_dct[k] == 0:
    #                                     pass
    #                                 else:
    #                                     if k in fa_formula_dct.keys():
    #                                         formula_dg_dct[k] = formula_dct[k] - fa_formula_dct[k]
    #                                     else:
    #                                         if k == 'Na':
    #                                             pass
    #                                         else:
    #                                             formula_dg_dct[k] = formula_dct[k]
    #
    #                             for k in ['C', 'H', 'N', 'O', 'Na']:
    #                                 if k in formula_dg_dct.keys():
    #                                     if formula_dg_dct[k] == 0:
    #                                         pass
    #                                     else:
    #                                         formula_dg = formula_dg + k + str(formula_dg_dct[k])
    #
    #                             fa_elem_mass_dha = self.fa_def_df.loc[self.fa_def_df['FA'] == '22:06']['elem']
    #                             fa_formula_dct_dha = IsotopeHunter().get_elements(fa_elem_mass_dha.iloc[0])
    #                             formula_dg_dct_dha = {}
    #                             fa_formula_dct_dha['O'] = int(fa_formula_dct_dha['O']) - 1
    #                             fa_formula_dct_dha['H'] = int(fa_formula_dct_dha['H']) - 2
    #                             for k in formula_dct.keys():
    #                                 if formula_dct[k] == 0:
    #                                     pass
    #                                 else:
    #                                     if k in fa_formula_dct_dha.keys():
    #                                         formula_dg_dct_dha[k] = formula_dct[k] - fa_formula_dct_dha[k]
    #                                     else:
    #                                         if k == 'Na':
    #                                             pass
    #                                         else:
    #                                             formula_dg_dct_dha[k] = formula_dct[k]
    #                             isotope_dha = IsotopeHunter().get_isotope_mz(formula_dg_dct_dha, isotope_number=2)
    #
    #                             fa_elem_mass_pa = self.fa_def_df.loc[self.fa_def_df['FA'] == '14:00']['elem']
    #                             fa_formula_dct_pa = IsotopeHunter().get_elements(fa_elem_mass_pa.iloc[0])
    #                             fa_formula_dct_pa['O'] = int(fa_formula_dct_pa['O']) - 1
    #                             fa_formula_dct_pa['H'] = int(fa_formula_dct_pa['H']) - 2
    #                             formula_dg_dct_pa = {}
    #                             for k in formula_dct.keys():
    #                                 if formula_dct[k] == 0:
    #                                     pass
    #                                 else:
    #                                     if k in fa_formula_dct_pa.keys():
    #                                         formula_dg_dct_pa[k] = formula_dct[k] - fa_formula_dct_pa[k]
    #                                     else:
    #                                         if k == 'Na':
    #                                             pass
    #                                         else:
    #                                             formula_dg_dct_pa[k] = formula_dct[k]
    #                             isotope_pa = IsotopeHunter().get_isotope_mz(formula_dg_dct_pa, isotope_number=2)
    #
    #                             pra_m2_i_ratio = ((100 * base_m2_i) / base_m1_i) / 100
    #                             if pra_m2_i_ratio <= isotope_pa.iloc[2]['ratio'] and pra_m2_i_ratio >= isotope_dha.iloc[2]['ratio']:
    #                                 print('Check in case of an isotope')
    #                                 isotope_flag = IsotopeHunter().get_isotope_fragments(combine_small2.iloc[1]['mz'], combine_small2.iloc[1]['i'], formula_dg,
    #                                                                                ms2_df, isotope_number=2,
    #                                                                                ms1_precision=ms2_precision)
    #                                 if isotope_flag == 1:
    #                                     combine_small2_index = combine_small2.index.tolist()
    #                                     combine_all_lst = combine_all_lst.drop([combine_small2_index[1]])
    #                             else:
    #                                 print('Nop an isotope')
    #                     else:
    #                         print('It is not an isotope')
    #                 else:
    #                     pass
    #
    #
    #             fa_ident_df = combine_all_lst[combine_all_lst['Type']=='FA'].sort_values(by='i', ascending=False)
    #             fa_ident_df = fa_ident_df.reset_index(drop=True)
    #             lyso_ident_df = combine_all_lst[combine_all_lst['Type'] == 'Lyso'].sort_values(by='i', ascending=False)
    #             lyso_ident_df = lyso_ident_df.reset_index(drop=True)
    #             lyso_w_ident_df = combine_all_lst[combine_all_lst['Type'] == 'LysoW'].sort_values(by='i', ascending=False)
    #             lyso_w_ident_df = lyso_w_ident_df.reset_index(drop=True)
    #             mg_w_ident_df = combine_all_lst[combine_all_lst['Type'] == 'MG'].sort_values(by='i', ascending=False)
    #             mg_w_ident_df = mg_w_ident_df.reset_index(drop=True)
    #
    #             # print lyso_w_ident_df
    #             self.weight_df['mz'] = 0.0
    #             for _i, _weight_se in self.weight_df.iterrows():
    #                 _type = _weight_se['Type']
    #                 _weight = _weight_se['Weight']
    #                 weight_dct[_type] = _weight
    #             ######################################################################
    #             #
    #             #   Also the changes and the calculation must be done before this part
    #             #   Because the part below here start the calculation of the score and thats is why is need to be done before in this part
    #             ########################################################################
    #             if abbr[:2] in ['TG']:
    #                 #print "Miria"
    #                 # exit()
    #                 # print lipid_abbr_df
    #                 for _i_abbr, _abbr_se in lipid_abbr_df.iterrows():
    #                     # _pl_abbr = _abbr_se['Lipid_abbr']
    #                     _sn1_abbr = _abbr_se['sn1_abbr']
    #                     _sn2_abbr = _abbr_se['sn2_abbr']
    #                     _sn3_abbr = _abbr_se['sn3_abbr']
    #
    #                     if _sn1_abbr in fa_ident_lst:
    #                         _rank_sn1 = fa_ident_lst.index(_sn1_abbr)
    #                         r_sn1_i = 100 * fa_i_lst[_rank_sn1] / ms2_max_i
    #                         lipid_abbr_df.set_value(_i_abbr, 'i_sn1', r_sn1_i)
    #                         if rank_mode is True:
    #                             lipid_abbr_df.set_value(_i_abbr, 'sn1', weight_dct['sn1'] * (10 - _rank_sn1) / 10)
    #                         else:
    #                             lipid_abbr_df.set_value(_i_abbr, 'sn1', weight_dct['sn1'] * r_sn1_i * 0.01)
    #
    #                     if _sn2_abbr in fa_ident_lst:
    #                         _rank_sn2 = fa_ident_lst.index(_sn2_abbr)
    #                         r_sn2_i = 100 * fa_i_lst[_rank_sn2] / ms2_max_i
    #                         lipid_abbr_df.set_value(_i_abbr, 'i_sn2', r_sn2_i)
    #                         if rank_mode is True:
    #                             lipid_abbr_df.set_value(_i_abbr, 'sn2', weight_dct['sn2'] * (10 - _rank_sn2) / 10)
    #                         else:
    #                             lipid_abbr_df.set_value(_i_abbr, 'sn2', weight_dct['sn2'] * r_sn2_i * 0.01)
    #
    #                     if _sn3_abbr in fa_ident_lst:
    #                         _rank_sn3 = fa_ident_lst.index(_sn3_abbr)
    #                         r_sn3_i = 100 * fa_i_lst[_rank_sn3] / ms2_max_i
    #                         lipid_abbr_df.set_value(_i_abbr, 'i_sn3', r_sn3_i)
    #                         if rank_mode is True:
    #                             lipid_abbr_df.set_value(_i_abbr, 'sn3', weight_dct['sn3'] * (10 - _rank_sn3) / 10)
    #                         else:
    #                             lipid_abbr_df.set_value(_i_abbr, 'sn3', weight_dct['sn3'] * r_sn3_i * 0.01)
    #                     if _sn1_abbr in lyso_ident_lst:
    #                         _rank_l_sn1 = lyso_ident_lst.index(_sn1_abbr)
    #                         r_lyso1_i = 100 * lyso_i_lst[_rank_l_sn1] / ms2_max_i
    #                         lipid_abbr_df.set_value(_i_abbr, 'i_[M+H]-sn1', r_lyso1_i)
    #                         if rank_mode is True:
    #                             lipid_abbr_df.set_value(_i_abbr, '[M+H]-sn1',
    #                                                     weight_dct['[M+H]-sn1'] * (10 - _rank_l_sn1) / 10)
    #                         else:
    #                             lipid_abbr_df.set_value(_i_abbr, '[M+H]-sn1', weight_dct['[M+H]-sn1'] * r_lyso1_i * 0.01)
    #                     if _sn2_abbr in lyso_ident_lst:
    #                         _rank_l_sn2 = lyso_ident_lst.index(_sn2_abbr)
    #                         r_lyso2_i = 100 * lyso_i_lst[_rank_l_sn2] / ms2_max_i
    #                         lipid_abbr_df.set_value(_i_abbr, 'i_[M+H]-sn2', r_lyso2_i)
    #                         if rank_mode is True:
    #                             lipid_abbr_df.set_value(_i_abbr, '[M+H]-sn2',
    #                                                     weight_dct['[M+H]-sn2'] * (10 - _rank_l_sn2) / 10)
    #                         else:
    #                             lipid_abbr_df.set_value(_i_abbr, '[M+H]-sn2', weight_dct['[M+H]-sn2'] * r_lyso2_i * 0.01)
    #                     if _sn3_abbr in lyso_ident_lst:
    #                         _rank_l_sn3 = lyso_ident_lst.index(_sn3_abbr)
    #                         r_lyso3_i = 100 * lyso_i_lst[_rank_l_sn3] / ms2_max_i
    #                         lipid_abbr_df.set_value(_i_abbr, 'i_[M+H]-sn3', r_lyso3_i)
    #                         if rank_mode is True:
    #                             lipid_abbr_df.set_value(_i_abbr, '[M+H]-sn3',
    #                                                     weight_dct['[M+H]-sn3'] * (10 - _rank_l_sn3) / 10)
    #                         else:
    #                             lipid_abbr_df.set_value(_i_abbr, '[M+H]-sn3', weight_dct['[M+H]-sn3'] * r_lyso3_i * 0.01)
    #                     if _sn1_abbr in lyso_w_ident_lst:
    #                         _rank_lw_sn1 = lyso_w_ident_lst.index(_sn1_abbr)
    #                         r_lyso_w1_i = 100 * lyso_w_i_lst[_rank_lw_sn1] / ms2_max_i
    #                         lipid_abbr_df.set_value(_i_abbr, 'i_[M+H]-sn1-H2O', r_lyso_w1_i)
    #                         #############################################
    #                         #
    #                         # Can go in a function to avoid the double check
    #                         #
    #                         ##############################################
    #                         if rank_mode is True:
    #                             if sodiumFlag == 1:
    #                                 lipid_abbr_df.set_value(_i_abbr, '[M+Na]-sn1',
    #                                                         weight_dct['[M+Na]-sn1'] * (10 - _rank_lw_sn1) / 10)
    #                             else:
    #                                 lipid_abbr_df.set_value(_i_abbr, '[M+H]-sn1-H2O',
    #                                                     weight_dct['[M+H]-sn1-H2O'] * (10 - _rank_lw_sn1) / 10)
    #                         else:
    #                             if sodiumFlag == 1:
    #                                 lipid_abbr_df.set_value(_i_abbr, '[M+Na]-sn1',
    #                                                         weight_dct['[M+Na]-sn1'] * r_lyso_w1_i * 0.01)
    #                             else:
    #                                 lipid_abbr_df.set_value(_i_abbr, '[M+H]-sn1-H2O',
    #                                                     weight_dct['[M+H]-sn1-H2O'] * r_lyso_w1_i * 0.01)
    #                     if _sn2_abbr in lyso_w_ident_lst:
    #                         _rank_lw_sn2 = lyso_w_ident_lst.index(_sn2_abbr)
    #                         r_lyso_w2_i = 100 * lyso_w_i_lst[_rank_lw_sn2] / ms2_max_i
    #                         lipid_abbr_df.set_value(_i_abbr, 'i_[M+H]-sn2-H2O', r_lyso_w2_i)
    #                         if rank_mode is True:
    #                             if sodiumFlag == 1:
    #                                 lipid_abbr_df.set_value(_i_abbr, '[M+Na]-sn2',
    #                                                         weight_dct['[M+Na]-sn2'] * (10 - _rank_lw_sn2) / 10)
    #                             else:
    #                                 lipid_abbr_df.set_value(_i_abbr, '[M+H]-sn2-H2O',
    #                                                     weight_dct['[M+H]-sn2-H2O'] * (10 - _rank_lw_sn2) / 10)
    #                         else:
    #                             if sodiumFlag ==1:
    #                                 lipid_abbr_df.set_value(_i_abbr, '[M+Na]-sn2',
    #                                                         weight_dct['[M+Na]-sn2'] * r_lyso_w2_i * 0.01)
    #                             else:
    #                                 lipid_abbr_df.set_value(_i_abbr, '[M+H]-sn2-H2O',
    #                                                     weight_dct['[M+H]-sn2-H2O'] * r_lyso_w2_i * 0.01)
    #                     if _sn3_abbr in lyso_w_ident_lst:
    #                         _rank_lw_sn3 = lyso_w_ident_lst.index(_sn3_abbr)
    #                         r_lyso_w3_i = 100 * lyso_w_i_lst[_rank_lw_sn3] / ms2_max_i
    #                         lipid_abbr_df.set_value(_i_abbr, 'i_[M+H]-sn3-H2O', r_lyso_w3_i)
    #                         if rank_mode is True:
    #                             if sodiumFlag == 1:
    #                                 lipid_abbr_df.set_value(_i_abbr, '[M+Na]-sn3',
    #                                                         weight_dct['[M+Na]-sn3'] * (10 - _rank_lw_sn3) / 10)
    #                             else:
    #                                 lipid_abbr_df.set_value(_i_abbr, '[M+H]-sn3-H2O',
    #                                                     weight_dct['[M+H]-sn3-H2O'] * (10 - _rank_lw_sn3) / 10)
    #                         else:
    #                             if sodiumFlag == 1:
    #                                 lipid_abbr_df.set_value(_i_abbr, '[M+Na]-sn3',
    #                                                         weight_dct['[M+Na]-sn3'] * r_lyso_w3_i * 0.01)
    #                             else:
    #                                 lipid_abbr_df.set_value(_i_abbr, '[M+H]-sn3-H2O',
    #                                                     weight_dct['[M+H]-sn3-H2O'] * r_lyso_w3_i * 0.01)
    #                     if _sn1_abbr in mg_w_ident_lst:
    #                         _rank_mgw_sn1 = mg_w_ident_lst.index(_sn1_abbr)
    #                         r_mg_w_i = 100 * mg_w_i_lst[_rank_mgw_sn1]/ms2_max_i
    #                         lipid_abbr_df.set_value(_i_abbr, 'i_[M+H]-(sn2+sn3)-H2O', r_mg_w_i)
    #                         if rank_mode is True:
    #                             lipid_abbr_df.set_value(_i_abbr, '[M+H]-(sn2+sn3)-H2O', weight_dct['[M+H]-(sn2+sn3)-H2O'] * (10 - _rank_mgw_sn1) / 10)
    #                         else:
    #                             lipid_abbr_df.set_value(_i_abbr, '[M+H]-(sn2+sn3)-H2O', weight_dct['[M+H]-(sn2+sn3)-H2O'] * r_mg_w_i * 0.01)
    #
    #                     if _sn2_abbr in mg_w_ident_df:
    #                         _rank_mgw_sn2 = mg_w_ident_lst.index(_sn2_abbr)
    #                         r_mg_w_i = 100 * mg_w_i_lst[_rank_mgw_sn2]/ms2_max_i
    #                         lipid_abbr_df.set_value(_i_abbr, 'i_[M+H]-(sn1+sn3)-H2O', r_mg_w_i)
    #                         if rank_mode is True:
    #                             lipid_abbr_df.set_value(_i_abbr, '[M+H]-(sn1+sn3)-H2O', weight_dct['[M+H]-(sn1+sn3)-H2O'] * (10 - _rank_mgw_sn2) / 10)
    #                         else:
    #                             lipid_abbr_df.set_value(_i_abbr, '[M+H]-(sn1+sn3)-H2O', weight_dct['[M+H]-(sn1+sn3)-H2O'] * r_mg_w_i * 0.01)
    #
    #                     if _sn3_abbr in mg_w_ident_df:
    #                         _rank_mgw_sn3 = mg_w_ident_lst.index(_sn3_abbr)
    #                         r_mg_w_i = 100 * mg_w_i_lst[_rank_mgw_sn3] / ms2_max_i
    #                         lipid_abbr_df.set_value(_i_abbr, 'i_[M+H]-(sn1+sn2)-H2O', r_mg_w_i)
    #                         if rank_mode is True:
    #                             lipid_abbr_df.set_value(_i_abbr, '[M+H]-(sn1+sn2)-H2O',
    #                                                     weight_dct['[M+H]-(sn1+sn2)-H2O'] * (
    #                                                     10 - _rank_mgw_sn3) / 10)
    #                         else:
    #                             lipid_abbr_df.set_value(_i_abbr, '[M+H]-(sn1+sn2)-H2O',
    #                                                     weight_dct['[M+H]-(sn1+sn2)-H2O'] * r_mg_w_i * 0.01)
    #
    #             else:
    #
    #                 for _i_abbr, _abbr_se in lipid_abbr_df.iterrows():
    #                     # _pl_abbr = _abbr_se['Lipid_abbr']
    #                     _sn1_abbr = _abbr_se['sn1_abbr']
    #                     _sn2_abbr = _abbr_se['sn2_abbr']
    #
    #                     if _sn1_abbr in fa_ident_lst:
    #                         _rank_sn1 = fa_ident_lst.index(_sn1_abbr)
    #                         r_sn1_i = 100 * fa_i_lst[_rank_sn1] / ms2_max_i
    #                         lipid_abbr_df.set_value(_i_abbr, 'i_sn1', r_sn1_i)
    #                         if rank_mode is True:
    #                             lipid_abbr_df.set_value(_i_abbr, 'sn1', weight_dct['sn1'] * (10 - _rank_sn1) / 10)
    #                         else:
    #                             lipid_abbr_df.set_value(_i_abbr, 'sn1', weight_dct['sn1'] * r_sn1_i * 0.01)
    #                     if _sn2_abbr in fa_ident_lst:
    #                         _rank_sn2 = fa_ident_lst.index(_sn2_abbr)
    #                         r_sn2_i = 100 * fa_i_lst[_rank_sn2] / ms2_max_i
    #                         lipid_abbr_df.set_value(_i_abbr, 'i_sn2', r_sn2_i)
    #                         if rank_mode is True:
    #                             lipid_abbr_df.set_value(_i_abbr, 'sn2', weight_dct['sn2'] * (10 - _rank_sn2) / 10)
    #                         else:
    #                             lipid_abbr_df.set_value(_i_abbr, 'sn2', weight_dct['sn2'] * r_sn2_i * 0.01)
    #                     if _sn1_abbr in lyso_ident_lst:
    #                         _rank_l_sn1 = lyso_ident_lst.index(_sn1_abbr)
    #                         r_lyso1_i = 100 * lyso_i_lst[_rank_l_sn1] / ms2_max_i
    #                         lipid_abbr_df.set_value(_i_abbr, 'i_[M-H]-sn1', r_lyso1_i)
    #                         if rank_mode is True:
    #                             lipid_abbr_df.set_value(_i_abbr, '[M-H]-sn1',
    #                                                     weight_dct['[M-H]-sn1'] * (10 - _rank_l_sn1) / 10)
    #                         else:
    #                             lipid_abbr_df.set_value(_i_abbr, '[M-H]-sn1', weight_dct['[M-H]-sn1'] * r_lyso1_i * 0.01)
    #                     if _sn2_abbr in lyso_ident_lst:
    #                         _rank_l_sn2 = lyso_ident_lst.index(_sn2_abbr)
    #                         r_lyso2_i = 100 * lyso_i_lst[_rank_l_sn2] / ms2_max_i
    #                         lipid_abbr_df.set_value(_i_abbr, 'i_[M-H]-sn2', r_lyso2_i)
    #                         if rank_mode is True:
    #                             lipid_abbr_df.set_value(_i_abbr, '[M-H]-sn2',
    #                                                     weight_dct['[M-H]-sn2'] * (10 - _rank_l_sn2) / 10)
    #                         else:
    #                             lipid_abbr_df.set_value(_i_abbr, '[M-H]-sn2', weight_dct['[M-H]-sn2'] * r_lyso2_i * 0.01)
    #                     if _sn1_abbr in lyso_w_ident_lst:
    #                         _rank_lw_sn1 = lyso_w_ident_lst.index(_sn1_abbr)
    #                         r_lyso_w1_i = 100 * lyso_w_i_lst[_rank_lw_sn1] / ms2_max_i
    #                         lipid_abbr_df.set_value(_i_abbr, 'i_[M-H]-sn1-H2O', r_lyso_w1_i)
    #                         if rank_mode is True:
    #                             lipid_abbr_df.set_value(_i_abbr, '[M-H]-sn1-H2O',
    #                                                     weight_dct['[M-H]-sn1-H2O'] * (10 - _rank_lw_sn1) / 10)
    #                         else:
    #                             lipid_abbr_df.set_value(_i_abbr, '[M-H]-sn1-H2O',
    #                                                     weight_dct['[M-H]-sn1-H2O'] * r_lyso_w1_i * 0.01)
    #                     if _sn2_abbr in lyso_w_ident_lst:
    #                         _rank_lw_sn2 = lyso_w_ident_lst.index(_sn2_abbr)
    #                         r_lyso_w2_i = 100 * lyso_w_i_lst[_rank_lw_sn2] / ms2_max_i
    #                         lipid_abbr_df.set_value(_i_abbr, 'i_[M-H]-sn2-H2O', r_lyso_w2_i)
    #                         if rank_mode is True:
    #                             lipid_abbr_df.set_value(_i_abbr, '[M-H]-sn2-H2O',
    #                                                     weight_dct['[M-H]-sn2-H2O'] * (10 - _rank_lw_sn2) / 10)
    #                         else:
    #                             lipid_abbr_df.set_value(_i_abbr, '[M-H]-sn2-H2O',
    #                                                     weight_dct['[M-H]-sn2-H2O'] * r_lyso_w2_i * 0.01)
    #
    #             lipid_abbr_df['Score'] = lipid_abbr_df[weight_type_lst].sum(axis=1, numeric_only=True)
    #             match_reporter = 1
    #
    #
    #         else:
    #             print('!!!!!! NO FA identified =====>--> Skip >>> >>>')
    #         match_info_dct = {'MATCH_INFO': match_reporter, 'SCORE_INFO': lipid_abbr_df, 'FA_INFO': fa_ident_df,
    #                           'LYSO_INFO': lyso_ident_df, 'LYSO_W_INFO': lyso_w_ident_df, 'MG_W_INFO' : mg_w_ident_df}
    #         return match_info_dct

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

    # def get_fa_signals(self, lipid_info_dct, charge_type, mz_lib, ms2_df):
    #
    #     ident_df_dct = {}
    #     ident_checker = 0
    #
    #     sn1_fa = lipid_info_dct['SN1']
    #     sn2_fa = lipid_info_dct['SN2']
    #
    #     sn1_mz = self.fa_def_df.loc[self.fa_def_df['FA'] == sn1_fa, '[M-H]-'].values[0]
    #
    #     sn2_mz = self.fa_def_df.loc[self.fa_def_df['FA'] == sn2_fa, '[M-H]-'].values[0]
    #
    #     print('sn1_mz', sn1_mz)
    #     print('sn2_mz', sn2_mz)
    #
    #     # if PC OCP with COOH on sn2, the charge will be [M-H]-
    #     pc_chg_lst = ['[M+HCOO]-', '[M+CH3COO]-', '[M+FA]-', '[M+OAc]-']
    #     got_lpp_pr = False
    #     if (mz_lib, charge_type) in self.pr_info_lst:
    #         got_lpp_pr = True
    #     else:
    #         if charge_type in pc_chg_lst:
    #             if (mz_lib, '[M-H]-') in self.pr_info_lst:
    #                 got_lpp_pr = True
    #                 charge_type = '[M-H]-'
    #             else:
    #                 pass
    #         else:
    #             pass
    #     # End PC charge check
    #
    #     if got_lpp_pr:
    #         print('got PR in list', mz_lib, charge_type)
    #         pr_query_dct = self.pr_query_dct[mz_lib]
    #         if sn1_mz in pr_query_dct.keys():
    #             pass
    #         else:
    #             print('sn1_mz not in dict, try to reduce decimals -->', sn1_mz, round(sn1_mz, 6))
    #             if round(sn1_mz, 6) in pr_query_dct.keys():
    #                 sn1_mz = round(sn1_mz, 6)
    #                 print('found with new sn1_mz')
    #             else:
    #                 print('Not found with new sn1_mz')
    #         if sn2_mz in pr_query_dct.keys():
    #             pass
    #         else:
    #             print('sn2_mz not in dict, try to reduce decimals -->', sn2_mz, round(sn2_mz, 6))
    #             if round(sn2_mz, 6) in pr_query_dct.keys():
    #                 sn2_mz = round(sn2_mz, 6)
    #                 print('found with new sn2_mz')
    #             else:
    #                 print('Not found with new sn2_mz')
    #
    #         if sn1_mz in pr_query_dct.keys():
    #             _fa_dct = pr_query_dct[sn1_mz]
    #             print(_fa_dct)
    #             for _frag_type in ['sn1', '[M-H]-sn1', '[M-H]-sn1-H2O']:
    #
    #                 if _frag_type == 'sn1':
    #                     _frag_mz = _fa_dct['sn']
    #                     _frag_df = ms2_df.query(_fa_dct['[M-H]-_query'])
    #                 elif _frag_type == '[M-H]-sn1':
    #                     _frag_mz = _fa_dct['[M-H]-sn']
    #                     _frag_df = ms2_df.query(_fa_dct['[M-H]-sn_query'])
    #                 elif _frag_type == '[M-H]-sn1-H2O':
    #                     _frag_mz = _fa_dct['[M-H]-sn-H2O']
    #                     _frag_df = ms2_df.query(_fa_dct['[M-H]-sn-H2O_query'])
    #
    #                 if _frag_df.shape[0] > 0:
    #                     _frag_df.loc[:, 'ppm'] = 1e6 * (_frag_df['mz'] - _frag_mz) / _frag_mz
    #                     _frag_df.loc[:, 'ppm_abs'] = _frag_df['ppm'].abs()
    #
    #                     if _frag_df.shape[0] > 1:
    #                         _frag_i_df = _frag_df.sort_values(by='i', ascending=False).head(1)
    #                         _frag_ppm_df = _frag_df.sort_values(by='ppm_abs').head(1)
    #                         _frag_df = _frag_i_df.copy()
    #                         if _frag_ppm_df['i'].tolist() == _frag_i_df['i'].tolist():
    #                             pass
    #                         else:
    #                             _frag_df = _frag_i_df.append(_frag_ppm_df)
    #                         # convert df to dict
    #                         ident_df_dct[_frag_type] = _frag_df.iloc[0, :].to_dict()
    #
    #                         if _frag_type == 'sn1':
    #                             ident_df_dct[_frag_type]['Proposed_structures'] = sn1_fa
    #                             ident_checker += 1
    #                         else:
    #                             ident_df_dct[_frag_type]['Proposed_structures'] = _frag_type
    #                     elif _frag_df.shape[0] == 1:
    #                         # convert df to dict
    #                         ident_df_dct[_frag_type] = _frag_df.iloc[0, :].to_dict()
    #
    #                         if _frag_type == 'sn1':
    #                             ident_df_dct[_frag_type]['Proposed_structures'] = sn1_fa
    #                             ident_checker += 1
    #                         else:
    #                             ident_df_dct[_frag_type]['Proposed_structures'] = _frag_type
    #
    #         if sn2_mz in pr_query_dct.keys():
    #             _fa_dct = pr_query_dct[sn2_mz]
    #             print(_fa_dct)
    #             for _frag_type in ['sn2', '[M-H]-sn2', '[M-H]-sn2-H2O']:
    #
    #                 if _frag_type == 'sn2':
    #                     _frag_mz = _fa_dct['sn']
    #                     _frag_df = ms2_df.query(_fa_dct['[M-H]-_query'])
    #                 elif _frag_type == '[M-H]-sn2':
    #                     _frag_mz = _fa_dct['[M-H]-sn']
    #                     _frag_df = ms2_df.query(_fa_dct['[M-H]-sn_query'])
    #                 elif _frag_type == '[M-H]-sn2-H2O':
    #                     _frag_mz = _fa_dct['[M-H]-sn-H2O']
    #                     _frag_df = ms2_df.query(_fa_dct['[M-H]-sn-H2O_query'])
    #
    #                 if _frag_df.shape[0] > 0:
    #                     _frag_df.loc[:, 'ppm'] = 1e6 * (_frag_df['mz'] - _frag_mz) / _frag_mz
    #                     _frag_df.loc[:, 'ppm_abs'] = _frag_df['ppm'].abs()
    #
    #                     if _frag_df.shape[0] > 1:
    #                         _frag_i_df = _frag_df.sort_values(by='i', ascending=False).head(1)
    #                         _frag_ppm_df = _frag_df.sort_values(by='ppm_abs').head(1)
    #                         _frag_df = _frag_i_df.copy()
    #                         if _frag_ppm_df['i'].tolist() == _frag_i_df['i'].tolist():
    #                             pass
    #                         else:
    #                             _frag_df = _frag_i_df.append(_frag_ppm_df)
    #                         # convert df to dict
    #                         ident_df_dct[_frag_type] = _frag_df.iloc[0, :].to_dict()
    #
    #                         if _frag_type == 'sn2':
    #                             ident_df_dct[_frag_type]['Proposed_structures'] = sn2_fa
    #                             ident_checker += 1
    #                         else:
    #                             ident_df_dct[_frag_type]['Proposed_structures'] = _frag_type
    #
    #                     elif _frag_df.shape[0] == 1:
    #                         # convert df to dict
    #                         ident_df_dct[_frag_type] = _frag_df.iloc[0, :].to_dict()
    #                         if _frag_type == 'sn2':
    #                             ident_df_dct[_frag_type]['Proposed_structures'] = sn2_fa
    #                             ident_checker += 1
    #                         else:
    #                             ident_df_dct[_frag_type]['Proposed_structures'] = _frag_type
    #
    #     return ident_df_dct, ident_checker

    def get_rankscore(self, master_info_df, abbr_bulk, charge, mz_lib, ms2_df, lipid_type):

        lite_info_df = master_info_df.query('Bulk_ABBR == "%s"' % abbr_bulk)

        query_type_lst = []

        if charge == '[M-H]-':
            if lipid_type in ['PA', 'PC', 'PE', 'PG', 'PI', 'PS']:
                query_type_lst = ['SN1_[FA-H]-', 'SN2_[FA-H]-', '[LPL(SN1)-H2O-H]-', '[LPL(SN1)-H]-',
                                  '[LPL(SN2)-H2O-H]-', '[LPL(SN2)-H]-']

        q_results_df = pd.DataFrame()

        for _idx, _info_se in lite_info_df.iterrows():
            for _q in query_type_lst:
                if _q in self.weight_type_lst:
                    _q_str = _info_se[_q + '_Q']
                    print(_info_se[_q + '_ABBR'], _q, _q_str)
                    _q_tmp_df = ms2_df.query(_q_str)
                    if _q_tmp_df.shape[1] > 0:
                        _q_tmp_df['lib_mz'] = _info_se[_q + '_MZ']
                        _q_tmp_df['obs_ppm'] = 1e6 * (_q_tmp_df['mz'] - _q_tmp_df['lib_mz']) / _q_tmp_df['lib_mz']
                        _q_tmp_df['obs_ppm'] = _q_tmp_df['obs_ppm'].astype(int)
                        _q_tmp_df['obs_ppm_abs'] = _q_tmp_df['obs_ppm'].abs()
                        _q_tmp_df['obs_abbr'] = _info_se[_q + '_ABBR']
                        _q_tmp_df['obs_type'] = _q
                        _q_tmp_df['score_group'] = self.weight_dct[_q]['Group']
                        q_results_df = q_results_df.append(_q_tmp_df)

        # print(query_results_dct)
        q_results_df = q_results_df.sort_values(by=['obs_type', 'obs_ppm_abs', 'i'],
                                                ascending=[True, True, False])
        print(q_results_df)

        if len(self.weight_group_lst) > 1:
            for _g in self.weight_group_lst:
                _sub_query_results_df = q_results_df[q_results_df['score_group'] == _g]
                _sub_query_results_df = _sub_query_results_df.sort_values(by=['obs_ppm_abs', 'i'],
                                                                          ascending=[True, False])
                _sub_query_results_df = _sub_query_results_df.drop_duplicates(subset=['obs_abbr', 'obs_type'],
                                                                              keep='first')
                print(_sub_query_results_df)

        #         else:
        #             print('!! No structure related signals found !!')
        #             match_info_dct = {'MATCH_INFO': matched_checker, 'Rank_score': 0.0}
        #     else:
        #         print('!! No structure related signals found !!')
        #         match_info_dct = {'MATCH_INFO': matched_checker, 'Rank_score': 0.0}
        # else:
        #     print('!! No structure related signals found !!')
        #
        #     match_info_dct = {'MATCH_INFO': matched_checker, 'Rank_score': 0.0}
        matched_checker = 0
        match_info_dct = {'MATCH_INFO': matched_checker, 'Rank_score': 0.0}

        return match_info_dct, matched_checker


if __name__ == '__main__':
    pass

# pre_Proposed.Propose_pre_str()
