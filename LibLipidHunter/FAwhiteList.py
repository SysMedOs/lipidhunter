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

from IsotopeHunter import  IsotopeHunter
import pandas as pd


class FA_list(object):
    def __init__(self):
        self.periodic_table_dct = {'H': [1.0078250321, 0.999885],
                                   'D': [2.0141017780, 0.0001157],
                                   'C': [12.0, 0.9893],
                                   'N': [14.0030740052, 0.99632],
                                   'O': [15.9949146221, 0.99757],
                                   'Na': [22.98976967, 1.0],
                                   'P': [30.97376151, 1.0],
                                   'S': [31.97207069, 0.9493],
                                   'K': [38.9637069, 0.932581]
                                   }

    def elemental_composition (self, link_t, c_t, db_t):
        if link_t == 'A':
            element_dict = { 'C' : c_t, 'H' : 2*c_t - 2*db_t, 'O' : 2}
            elemental_comp = 'C%dH%dO2' % (c_t, 2*c_t - 2*db_t)
        elif link_t == 'O':
            element_dict = {'C' : c_t, 'H' : 2*c_t - 2*db_t + 2, 'O': 1}
            elemental_comp = 'C%dH%dO1' % (c_t, 2*c_t - 2*db_t + 2)
        elif link_t == 'P':
            element_dict = {'C' : c_t, 'H' : 2*c_t - 2*db_t, 'O': 1}
            elemental_comp = 'C%dH%dO1' % (c_t, 2*c_t - 2*db_t)

        return element_dict, elemental_comp

    def M2_isotope (self):
        pass

    def pos_geo (self, elem_dct, mode):
        print ('Come on know')
        info_dct={}
        elem_dct['H'] = elem_dct['H'] - 1
        elem_dct['O'] = elem_dct['O'] - 1
        pos_value = IsotopeHunter().get_mono_mz(elem_dct)
        info_dct['pos_value'] = pos_value
        elem_dct['H'] = elem_dct['H'] - 1
        m_water = IsotopeHunter().get_mono_mz(elem_dct)
        info_dct['m_water'] = m_water
        elem_dct['C'] = elem_dct['C'] + 3
        elem_dct['H'] = elem_dct['H'] + 7
        elem_dct['O'] = elem_dct['O'] + 2
        mg_value = IsotopeHunter().get_mono_mz(elem_dct)
        info_dct['mg_value'] = mg_value
        if mode in ['[M+Na]+']:
            elem_dct['C'] = elem_dct['C'] - 3
            elem_dct['H'] = elem_dct['H'] - 6
            elem_dct['Na'] = 1
        m_Na = IsotopeHunter().get_mono_mz(elem_dct)
        info_dct['m_Na'] = m_Na
        return info_dct




    def negative_column (self, elem_dct):
        info_dct={}
        print elem_dct
        print ('hellooo')
        elem_dct['H'] = elem_dct['H'] - 1
        neg_value = IsotopeHunter().get_mono_mz(elem_dct)
        print neg_value
        elem_dct['O'] = elem_dct['O'] -1
        elem_dct['H'] = elem_dct['H'] - 2
        print elem_dct
        neg_value_w = IsotopeHunter().get_mono_mz(elem_dct)
        print neg_value_w
        elem_dct['H'] = elem_dct['H'] + 1
        m_water = IsotopeHunter().get_mono_mz(elem_dct)
        print elem_dct
        print m_water
        info_dct['neg_value'] = neg_value
        info_dct['neg_value_w'] = neg_value_w
        info_dct['m_water'] = m_water
        print info_dct
        print ('well')
        return info_dct