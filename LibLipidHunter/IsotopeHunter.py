# -*- coding: utf-8 -*-
# Copyright 2016-2017 LPP team, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of LipidHunter.
# For more info please contact:
#     LPP team oxlpp@bbz.uni-leipzig.de
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de
#     Developer Georgia Angelidou georgia.angelidou@uni-leipzig.de

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

    def get_isotope_mz(self, elem_dct, only_c=False, isotope_number=2):

        # calc M+0 --> M+2
        if isotope_number == 2:
            isotope_count_lst = [1, 2]
        else:
            isotope_count_lst = range(1, isotope_number + 1)

        mono_mz = self.get_mono_mz(elem_dct)

        # consider C only
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
                isotope_pattern = list(isotope_pattern_calc.coef)[:isotope_number + 1]
            except ValueError:
                print('==>Elements error --> change to 13C mode for this compound -->')
                # consider C only --> binomial expansion 3x faster
                isotope_pattern = stats.binom.pmf(range(0, isotope_number + 1), c_count, ration_13c12c)

        m0_i = isotope_pattern[0]
        isotope_pattern = [x / m0_i for x in isotope_pattern]

        isotope_distribution_df = pd.DataFrame(data={'mz': isotope_mz_lst, 'ratio': isotope_pattern})
        isotope_distribution_df = isotope_distribution_df.round({'ratio': 2})

        return isotope_distribution_df

    @staticmethod
    def calc_isotope_score(isotope_pattern_df, spec_df, ms1_precision, ms1_pr_i, deconv=[]):

        isotope_checker_dct = {}
        isotope_score_delta = 0
        isotope_m1_score_delta = 0
        m2_i = 0
        theo_i_lst = []

        isotope_score = 0
        isotope_m1_score = 0

        if ms1_pr_i > 0:
            for _i, _se in isotope_pattern_df.iterrows():

                if len(deconv) == 3:
                    _base_i = deconv[_i]
                else:
                    _base_i = 0

                # [M+0] has _i == 0
                _mz = _se['mz']
                _ratio = _se['ratio']
                _mz_delta = _mz * ms1_precision
                _i_df = spec_df.query('%f <= mz <= %f' % (_mz - _mz_delta, _mz + _mz_delta))
                if _i < 2:
                    theo_i = ms1_pr_i * _ratio + _base_i
                else:
                    theo_i = ms1_pr_i * _ratio
                _i_info_dct = {'theo_mz': _mz, 'theo_i': theo_i, 'theo_ratio': _ratio}
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

                theo_i_lst.append(theo_i)

                _i_r = _i_max / ms1_pr_i
                _i_info_dct['obs_ratio'] = _i_r
                isotope_checker_dct[_i] = _i_info_dct

                if _i > 0:
                    isotope_score_delta += abs(_i_r - _ratio)

                    if _i == 1:
                        isotope_m1_score_delta += abs(_i_r - _ratio)
                    else:
                        pass
                    if _i == 2:
                        m2_i = _i_max

            isotope_score = 100 * (1 - isotope_score_delta)
            isotope_m1_score = 100 * (1 - isotope_m1_score_delta)

        isotope_calc_dct = {'isotope_checker_dct': isotope_checker_dct, 'isotope_score': isotope_score,
                            'isotope_m1_score': isotope_m1_score, 'm2_i': m2_i, 'theo_i_lst': theo_i_lst}

        return isotope_calc_dct

    def get_deconvolution(self, elem_dct, spec_df, mz_delta, base_i, only_c=False):

        base_m1_i = 0
        base_m2_i = 0
        base_m3_i = 0

        m_pre1_isotope_pattern_df = self.get_isotope_mz(elem_dct, only_c=only_c, isotope_number=3)
        m_pre1_mz = m_pre1_isotope_pattern_df.get_value(0, 'mz')
        pre_i_df = spec_df.query('%f <= mz <= %f' %
                                 (m_pre1_mz - mz_delta, m_pre1_mz + mz_delta))

        if pre_i_df.shape[0] > 0:
            max_m_pre1_i = pre_i_df['i'].max()
            max_m_pre1_i -= base_i
            if max_m_pre1_i > 0:
                base_m1_i += max_m_pre1_i * m_pre1_isotope_pattern_df.get_value(1, 'ratio')
                base_m2_i += max_m_pre1_i * m_pre1_isotope_pattern_df.get_value(2, 'ratio')
                base_m3_i += max_m_pre1_i * m_pre1_isotope_pattern_df.get_value(3, 'ratio')

        return base_m1_i, base_m2_i, base_m3_i

    def get_isotope_score(self, ms1_pr_mz, ms1_pr_i, formula, spec_df, isotope_number=2,
                          ms1_precision=50e-6, pattern_tolerance=5, only_c=False, score_filter=75, decon=True):

        mz_delta = ms1_pr_mz * ms1_precision
        delta_13c = 1.0033548378

        if decon is True:
            deconv_elem_dct = self.get_elements(formula)
            # M-2
            deconv_elem_dct['H'] += -2
            base_i = 0
            pre2_base_m1_i, pre2_base_m2_i, pre2_base_m3_i = self.get_deconvolution(deconv_elem_dct, spec_df,
                                                                                    mz_delta, base_i, only_c=only_c)

            # M+0
            # m0_base_abs = pre2_base_m2_i + pre1_base_m1_i
            m0_base_abs = pre2_base_m2_i
            base_m1_i, base_m2_i, base_m3_i = self.get_deconvolution(self.get_elements(formula), spec_df, mz_delta,
                                                                     m0_base_abs, only_c=only_c)
            # M+1
            m1_base_abs = pre2_base_m3_i

            # M+2
            m2_base_abs = base_m2_i

            # M+3
            m3_base_abs = base_m3_i

            print('Deconvolution_i_abs_corrections: [M+0] %.1f, [M+1] %.1f, [M+2] %.1f, [M+3] %.1f'
                  % (m0_base_abs, m1_base_abs, m2_base_abs, m3_base_abs))
            deconv_lst = [m0_base_abs, m1_base_abs, m2_base_abs, m3_base_abs]
        else:
            m0_base_abs = 0
            m1_base_abs = 0
            m2_base_abs = 0
            m3_base_abs = 0

            deconv_lst = [m0_base_abs, m1_base_abs, m2_base_abs, m3_base_abs]

        i_df = spec_df.query('%f <= mz <= %f' % (ms1_pr_mz - delta_13c - mz_delta, ms1_pr_mz - delta_13c + mz_delta))

        if i_df.shape[0] > 0:
            max_pre_m_i = i_df['i'].max()
            if ms1_pr_i > max_pre_m_i:
                elem_dct = self.get_elements(formula)
                mono_mz = self.get_mono_mz(elem_dct)
                if abs((ms1_pr_mz - mono_mz)) <= ms1_precision * ms1_pr_mz:
                    isotope_pattern_df = self.get_isotope_mz(elem_dct, only_c=only_c)

                    ms1_pr_i -= m0_base_abs
                    m0_deconv_lst = [m0_base_abs, m1_base_abs, m2_base_abs]
                    isotope_calc_dct = self.calc_isotope_score(isotope_pattern_df, spec_df,
                                                               ms1_precision, ms1_pr_i, deconv=m0_deconv_lst)

                    isotope_checker_dct = isotope_calc_dct['isotope_checker_dct']
                    isotope_score = isotope_calc_dct['isotope_score']
                    isotope_m1_score = isotope_calc_dct['isotope_m1_score']
                    m2_i = isotope_calc_dct['m2_i']

                    m2_checker_dct = {}
                    m2_score = 0

                    if isotope_score < score_filter:
                        # check if M+2 is potential M+0 of M+H2
                        # M+H2 elements

                        m2_elem_dct = self.get_elements(formula + 'H2')
                        m2_isotope_pattern_df = self.get_isotope_mz(m2_elem_dct, only_c=only_c)
                        m2_i -= m2_base_abs
                        m2_deconv_lst = [m2_base_abs, m3_base_abs, 0]
                        m2_calc_dct = self.calc_isotope_score(m2_isotope_pattern_df, spec_df,
                                                              ms1_precision, m2_i, deconv=m2_deconv_lst)

                        m2_checker_dct = m2_calc_dct['isotope_checker_dct']
                        # use M+1 only
                        m2_score = m2_calc_dct['isotope_m1_score']
                        m2_m2i = m2_calc_dct['m2_i']

                        if m2_score > 0:
                            pass
                        else:
                            m2_score = 0

                        if 2 in m2_checker_dct.keys():
                            del m2_checker_dct[2]

                        print('M+2-> M+4 has isotope score for [M+H2]: %.1f' % m2_score)
                        if m2_score >= 60 and isotope_m1_score >= score_filter:
                            isotope_score = isotope_m1_score
                            del isotope_checker_dct[2]
                        else:
                            pass
                else:
                    print('!! MS1 PR m/z not fit to Formula check bulk identification !!!!!!')
                    isotope_score = 0
                    isotope_checker_dct = {}
                    m2_checker_dct = {}
                    m2_score = 0

            else:
                print('MS1 PR is an isotope !!!!!!')
                isotope_score = 0
                isotope_checker_dct = {}
                m2_checker_dct = {}
                m2_score = 0
        else:
            isotope_score = 0
            isotope_checker_dct = {}
            m2_checker_dct = {}
            m2_score = 0

        isotope_score_info_dct = {'isotope_score': isotope_score, 'isotope_checker_dct': isotope_checker_dct,
                                  'm2_score': m2_score, 'm2_checker_dct': m2_checker_dct,
                                  'deconv_lst': deconv_lst}

        return isotope_score_info_dct


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
        isotope_distribute = isotope_distribute.head(3)
        print(isotope_distribute)
        print(isotope_distribute.get_value(0, 'mz'))
