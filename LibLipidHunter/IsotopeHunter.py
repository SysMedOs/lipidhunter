# -*- coding: utf-8 -*-
# Copyright 2015-2017 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

from __future__ import division
from __future__ import print_function


import re

import pandas as pd
from scipy import stats
from numpy.polynomial.polynomial import Polynomial


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

    def get_isotope_mz(self, elem_dct, isotope_number=2, only_c=False):

        if isotope_number <= 5:
            isotope_count_lst = range(1, isotope_number + 1)
        else:
            isotope_count_lst = [1, 2]

        mono_mz = self.get_mono_mz(elem_dct)

        # consider C only
        # delta_13c = self.periodic_table_dct['C'][1][0] - self.periodic_table_dct['C'][0][0]
        # ration_13c = self.periodic_table_dct['C'][1][0] - self.periodic_table_dct['C'][0][0]

        # r_12c = self.periodic_table_dct['C'][0][1]
        # r_13c = self.periodic_table_dct['C'][1][1]
        # r_1h = self.periodic_table_dct['H'][0][1]
        # r_2h = self.periodic_table_dct['H'][1][1]
        # r_14n = self.periodic_table_dct['N'][0][1]
        # r_15n = self.periodic_table_dct['N'][1][1]

        c_count = elem_dct['C']

        delta_13c = 1.0033548378
        ration_13c12c = 0.011
        isotope_mz_lst = [mono_mz]
        for _i_count in isotope_count_lst:
            _isotope_mz = mono_mz + delta_13c * _i_count
            isotope_mz_lst.append(_isotope_mz)

        # calc distribution by selected algorithms
        if only_c is True:
            # consider C only --> binomial expansion 3x faster
            isotope_pattern = stats.binom.pmf(range(0, isotope_number + 1), c_count, ration_13c12c)

        else:
            try:
                # consider more elements --> binomial/polynomial and McLaurin expansion [doi:10.1038/nmeth.3393]
                h_count = elem_dct['H']
                o_count = elem_dct['O']

                c_ploy = Polynomial((0.9893, 0.0107))
                h_ploy = Polynomial((0.999885, 0.0001157))
                o_ploy = Polynomial((0.99757, 0.00038, 0.00205))

                isotope_pattern_calc = c_ploy ** c_count * h_ploy ** h_count * o_ploy ** o_count

                if 'N' in elem_dct.keys():
                    n_count = elem_dct['N']
                    if n_count > 0:
                        n_ploy = Polynomial((0.99632, 0.00368))
                        isotope_pattern_calc *= n_ploy

                if 'S' in elem_dct.keys():
                    s_count = elem_dct['S']
                    if s_count > 0:
                        s_ploy = Polynomial((0.9493, 0.0076, 0.0429, 0.0002))
                        isotope_pattern_calc *= s_ploy

                if 'K' in elem_dct.keys():
                    k_count = elem_dct['K']
                    if k_count > 0:
                        k_ploy = Polynomial((0.932581, 0.000117, 0.067302))
                        isotope_pattern_calc *= k_ploy
                isotope_pattern = list(isotope_pattern_calc.coef)[:3]
            except ValueError:
                # consider C only --> binomial expansion 3x faster
                isotope_pattern = stats.binom.pmf(range(0, isotope_number + 1), c_count, ration_13c12c)

        m0_i = isotope_pattern[0]
        isotope_pattern = [x / m0_i for x in isotope_pattern]

        isotope_distribution_df = pd.DataFrame(data={'mz': isotope_mz_lst, 'ratio': isotope_pattern})
        isotope_distribution_df = isotope_distribution_df.round({'ratio': 2})
        # print(isotope_distribution_df)
        return isotope_distribution_df

    def get_isotope_score(self, ms1_pr_mz, ms1_pr_i, formula, spec_df, isotope_number=2,
                          ms1_precision=50e-6, pattern_tolerance=5, only_c=False):

        mz_delta = ms1_pr_mz * ms1_precision
        delta_13c = 1.0033548378
        i_df = spec_df.query('%f <= mz <= %f' % (ms1_pr_mz - delta_13c - mz_delta, ms1_pr_mz - delta_13c + mz_delta))

        if i_df.shape[0] > 0:
            max_pre_m_i = i_df['i'].max()
            if ms1_pr_i > max_pre_m_i:
                elem_dct = self.get_elements(formula)
                # mono_mz = self.get_mono_mz(elem_dct)
                isotope_pattern_df = self.get_isotope_mz(elem_dct, isotope_number=2)

                isotope_checker_dct = {}
                isotope_score_delta = 0
                for _i, _se in isotope_pattern_df.iterrows():

                    if _i > 0:
                        _mz = _se['mz']
                        _ratio = _se['ratio']
                        _mz_delta = _mz * ms1_precision
                        _i_df = spec_df.query('%f <= mz <= %f' % (_mz - _mz_delta, _mz + _mz_delta))
                        _i_info_dct = {'theo_mz': _mz, 'theo_i': ms1_pr_i * _ratio, 'theo_ratio': _ratio}
                        if _i_df.shape[0] > 0:
                            _i_df = _i_df.sort_values(by='i', ascending=False).head(1)
                            _i_max = _i_df['i'].tolist()[0]
                            _mz_max = _i_df['mz'].tolist()[0]
                            _i_info_dct['obs_i'] = _i_max
                            _i_info_dct['obs_mz'] = _mz_max
                        else:
                            _i_max = 0.0
                            _i_info_dct['obs_i'] = 0
                            _i_info_dct['obs_mz'] = 0
                        _i_r = _i_max / ms1_pr_i
                        _i_info_dct['obs_ratio'] = _i_r
                        isotope_score_delta += abs(_i_r - _ratio)
                        isotope_checker_dct[_i] = _i_info_dct

                isotope_score = 100 * (1 - isotope_score_delta)
                print('isotope_pattern_df', isotope_pattern_df)

            else:
                print('MS1 PR is an isotope !!!!!!')
                isotope_score = 0
                isotope_checker_dct = {}
        else:
            elem_dct = self.get_elements(formula)
            # mono_mz = self.get_mono_mz(elem_dct)
            isotope_pattern_df = self.get_isotope_mz(elem_dct, isotope_number=2, only_c=only_c)

            isotope_checker_dct = {}
            isotope_score_delta = 0
            for _i, _se in isotope_pattern_df.iterrows():

                if _i > 0:
                    _mz = _se['mz']
                    _ratio = _se['ratio']
                    _mz_delta = _mz * ms1_precision
                    _i_df = spec_df.query('%f <= mz <= %f' % (_mz - _mz_delta, _mz + _mz_delta))
                    _i_info_dct = {'theo_mz': _mz, 'theo_i': ms1_pr_i * _ratio, 'theo_ratio': _ratio}
                    if _i_df.shape[0] > 0:
                        _i_df = _i_df.sort_values(by='i', ascending=False).head(1)
                        _i_max = _i_df['i'].tolist()[0]
                        _mz_max = _i_df['mz'].tolist()[0]
                        _i_info_dct['obs_i'] = _i_max
                        _i_info_dct['obs_mz'] = _mz_max
                    else:
                        _i_max = 0.0
                        _i_info_dct['obs_i'] = 0
                        _i_info_dct['obs_mz'] = 0
                    _i_r = _i_max / ms1_pr_i
                    _i_info_dct['obs_ratio'] = _i_r
                    isotope_score_delta += abs(_i_r - _ratio)
                    isotope_checker_dct[_i] = _i_info_dct

            isotope_score = 100 * (1 - isotope_score_delta)
            print('isotope_pattern_df', isotope_pattern_df)

        return isotope_score, isotope_checker_dct


if __name__ == '__main__':
    # f = 'C39H67NO8P'  # PE(34:5)
    # f_lst = [f, f + 'K', f + 'Na', f + 'NH4', f + 'S', f + 'D']
    # usr_spec_df = pd.DataFrame()
    # iso_hunter = IsotopeHunter()
    # for _f in f_lst:
    #     print(_f)
    #     isotope_pattern_dct = iso_hunter.get_elements(_f)
    #
    #     print(isotope_pattern_dct)

    f = 'C41H71NO7P'  # PE
    # f_lst = [f, f + 'H+']
    f_lst = [f]
    usr_spec_df = pd.DataFrame()
    iso_hunter = IsotopeHunter()
    for _f in f_lst:
        print(_f)
        isotope_pattern_dct = iso_hunter.get_elements(_f)
        isotope_distribute = iso_hunter.get_isotope_mz(isotope_pattern_dct, only_c=False)
        print(isotope_distribute)
        isotope_distribute = iso_hunter.get_isotope_mz(isotope_pattern_dct, only_c=True)
        print(isotope_distribute)
