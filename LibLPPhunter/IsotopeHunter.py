# -*- coding: utf-8 -*-
# Copyright 2015-2016 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

from __future__ import division
from __future__ import print_function


import re

import pandas as pd
from scipy import stats


class IsotopeHunter(object):
    def __init__(self):
        # iupac '97
        self.periodic_table_dct = {'H': [(1.0078250321, 0.999885), (2.0141017780, 0.0001157)],
                                   'D': [(2.0141017780, 0.0001157)],
                                   'C': [(12.0, 0.9893), (13.0033548378, 0.0107)],
                                   'N': [(14.0030740052, 0.99632), (15.0001088984, 0.00368)],
                                   'O': [(15.9949146221, 0.99757), (16.99913150, 0.00038), (17.9991604, 0.00205)],
                                   'Na': [(22.98976967, 1.0)],
                                   'P': [(30.97376151, 1.0)],
                                   'S': [(31.97207069, 0.9493), (32.97145850, 0.0076),
                                         (33.96786683, 0.0429), (35.96708088, 0.0002)],
                                   'K': [(38.9637069, 0.932581), (39.96399867, 0.000117), (40.96182597, 0.067302)],
                                   }


    def get_elements(self, formula):
        elem_dct = {}
        elem_key_lst = self.periodic_table_dct.keys()
        tmp_formula = formula

        elem_lst = re.findall('[A-Z][a-z]*[0-9]*', formula)
        for _e in range(0, len(elem_lst)):

            _elem = re.findall('[A-Z][a-z]*', elem_lst[_e])
            _elem_count = re.findall('[0-9]+', elem_lst[_e])
            if len(_elem_count) == 0:
                _elem_count = 1
            else:
                _elem_count = sum([int(x) for x in _elem_count])
            if _elem[0] in elem_dct.keys():
                elem_dct[_elem[0]] += _elem_count
            else:
                elem_dct[_elem[0]] = _elem_count

        return elem_dct

    def get_mono_mz(self, elem_dct):

        mono_mz = 0.0
        for _elem in elem_dct.keys():
            mono_mz += elem_dct[_elem] * self.periodic_table_dct[_elem][0][0]

        return mono_mz

    def get_isotope_mz(self, elem_dct, isotope_number=2):

        if isotope_number <= 5:
            isotope_count_lst = range(1, isotope_number + 1)
        else:
            isotope_count_lst = [1, 2]

        mono_mz = self.get_mono_mz(elem_dct)

        # consider C only
        # delta_13c = self.periodic_table_dct['C'][1][0] - self.periodic_table_dct['C'][0][0]
        # ration_13c = self.periodic_table_dct['C'][1][0] - self.periodic_table_dct['C'][0][0]

        r_12c = self.periodic_table_dct['C'][0][1]
        r_13c = self.periodic_table_dct['C'][1][1]
        r_1h = self.periodic_table_dct['H'][0][1]
        r_2h = self.periodic_table_dct['H'][1][1]
        # r_14n = self.periodic_table_dct['N'][0][1]
        # r_15n = self.periodic_table_dct['N'][1][1]

        c_count = elem_dct['C']
        h_count = elem_dct['C']
        delta_13c = 1.0033548378
        ration_13c12c = 0.011
        isotope_mz_lst = [mono_mz]
        for _i_count in isotope_count_lst:
            _isotope_mz = mono_mz + delta_13c * _i_count
            isotope_mz_lst.append(_isotope_mz)

        # consider C only
        isotope_pattern = stats.binom.pmf(range(0, isotope_number + 1), c_count, ration_13c12c)
        m0_i = isotope_pattern[0]
        isotope_pattern = [100 * x / m0_i for x in isotope_pattern]

        isotope_distribution_df = pd.DataFrame(data={'mz': isotope_mz_lst, 'ratio': isotope_pattern})
        isotope_distribution_df['ratio'] = isotope_distribution_df['ratio'].astype(int)
        print(isotope_distribution_df)
        return isotope_distribution_df

    def check_isotope(self, formula, spec_df, isotope_number=2):

        elem_dct = self.get_elements(formula)
        # mono_mz = self.get_mono_mz(elem_dct)
        isotope_pattern_dct = self.get_isotope_mz(elem_dct, isotope_number=2)

        return isotope_pattern_dct


if __name__ == '__main__':
    f = 'C39H67NO8P'  # PE(34:5)
    f_lst = [f, f + 'K', f + 'Na', f + 'NH4', f + 'S', f + 'D']
    usr_spec_df = pd.DataFrame()
    isohunter = IsotopeHunter()
    for _f in f_lst:
        print(_f)
        isotope_pattern_dct = isohunter.check_isotope(_f, usr_spec_df, isotope_number=2)

        print(isotope_pattern_dct)
