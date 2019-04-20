# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
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
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de
#     Developer Georgia Angelidou georgia.angelidou@uni-leipzig.de

from __future__ import division
from __future__ import print_function

import re

from numpy.polynomial.polynomial import Polynomial
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
        elem_key_lst = list(self.periodic_table_dct.keys())
        tmp_formula = formula.strip('+')
        tmp_formula = tmp_formula.strip('-')

        elem_lst = re.findall('[A-Z][a-z]*[0-9]*', formula)
        for _e in range(0, len(elem_lst)):

            _elem = re.findall('[A-Z][a-z]*', elem_lst[_e])
            _elem_count = re.findall('[0-9]+', elem_lst[_e])
            if len(_elem_count) == 0:
                _elem_count = 1
            else:
                _elem_count = sum([int(x) for x in _elem_count])
            if _elem[0] in list(elem_dct.keys()):
                elem_dct[_elem[0]] += _elem_count
            else:
                elem_dct[_elem[0]] = _elem_count

        return elem_dct

    def get_mono_mz(self, elem_dct={}):
        mono_mz = 0.0
        for _elem in list(elem_dct.keys()):
            mono_mz += elem_dct[_elem] * self.periodic_table_dct[_elem][0][0]

        return mono_mz

    def get_isotope_mz(self, elem_dct, only_c=False, isotope_number=2):

        # calc M+0 --> M+2
        if isotope_number == 2:
            isotope_count_lst = [1, 2]
        else:
            isotope_count_lst = list(range(1, isotope_number + 1))

        # Calclulates the elemenatl mass from the elemental composition
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
            isotope_pattern = stats.binom.pmf(list(range(0, isotope_number + 1)), c_count, ration_13c12c)

        else:
            try:
                # consider more elements --> binomial/polynomial and McLaurin expansion [doi:10.1038/nmeth.3393]
                h_count = elem_dct['H']
                o_count = elem_dct['O']

                c_ploy = Polynomial((0.9893, 0.0107))
                h_ploy = Polynomial((0.999885, 0.0001157))
                o_ploy = Polynomial((0.99757, 0.00038, 0.00205))

                isotope_pattern_calc = c_ploy ** c_count * h_ploy ** h_count * o_ploy ** o_count

                if 'N' in list(elem_dct.keys()):
                    n_count = elem_dct['N']
                    if n_count > 0:
                        n_ploy = Polynomial((0.99632, 0.00368))
                        isotope_pattern_calc *= n_ploy

                if 'S' in list(elem_dct.keys()):
                    s_count = elem_dct['S']
                    if s_count > 0:
                        s_ploy = Polynomial((0.9493, 0.0076, 0.0429, 0.0002))
                        isotope_pattern_calc *= s_ploy

                if 'K' in list(elem_dct.keys()):
                    k_count = elem_dct['K']
                    if k_count > 0:
                        k_ploy = Polynomial((0.932581, 0.000117, 0.067302))
                        isotope_pattern_calc *= k_ploy
                isotope_pattern = list(isotope_pattern_calc.coef)[:isotope_number + 1]
            except ValueError:
                # print(_e)
                print('[INFO] --> Too large to use full elements for isotope pattern ... '
                      'use 13C only mode for this compound...')
                # consider C only --> binomial expansion 3x faster
                isotope_pattern = stats.binom.pmf(list(range(0, isotope_number + 1)), c_count, ration_13c12c)

        m0_i = isotope_pattern[0]
        isotope_pattern = [x / m0_i for x in isotope_pattern]

        isotope_distribution_df = pd.DataFrame(data={'mz': isotope_mz_lst, 'ratio': isotope_pattern})
        isotope_distribution_df = isotope_distribution_df.round({'ratio': 2})

        return isotope_distribution_df

    @staticmethod
    def peak_top_checker(ms1_pr_mz, spec_df, core_count=1, ms1_precision=50e-6):

        peak_top = False
        # top_obs_i = 0

        if ms1_precision <= 100e-6:
            top_precision = min(ms1_precision * 5, 100e-6)
        elif 100e-6 < ms1_precision <= 200e-6:
            top_precision = min(ms1_precision * 5, 200e-6)
        else:
            top_precision = min(ms1_precision * 2, 500e-6)
        pr_delta = ms1_pr_mz * ms1_precision

        ms_pr_df = spec_df.query('%f <= mz <= %f' % (ms1_pr_mz - pr_delta, ms1_pr_mz + pr_delta))
        if not ms_pr_df.empty:
            _i_df = ms_pr_df.sort_values(by='i', ascending=False)
            pr_obs_i = _i_df['i'].values.tolist()[0]
            pr_obs_mz = _i_df['mz'].values.tolist()[0]
            top_delta = pr_obs_mz * top_precision
            peak_top_df = spec_df.query('%f <= mz <= %f' % (ms1_pr_mz - top_delta, ms1_pr_mz + top_delta))
            top_obs_i = max(peak_top_df['i'].values.tolist())
        else:
            pr_obs_i = 0
            pr_obs_mz = 0
            top_obs_i = 0
            top_delta = 0

        if pr_obs_i >= 0.75 * top_obs_i > 0:
            peak_top = True
            # print(core_count,
            #       '[PASSED]MS1 PR m/z {pr} is peak top in +/- {delta} m/z, PR i {obs_i}, max i {max_i}'
            #       .format(pr=pr_obs_mz, delta=top_delta, obs_i=pr_obs_i, max_i=top_obs_i))
        else:
            print(core_count,
                  '[Warning] MS1 PR m/z {pr} is not peak top in +/- {delta} m/z, PR i {obs_i}, max i {max_i}'
                  .format(pr=pr_obs_mz, delta=top_delta, obs_i=pr_obs_i, max_i=top_obs_i))
        return peak_top, top_obs_i

    def calc_isotope_score(self, isotope_pattern_df, spec_df, ms1_precision, ms1_pr_i,
                           core_count=1, deconv=[], mode='m', score_filter=75):

        isotope_checker_dct = {}
        isotope_score_delta = 0
        isotope_m1_score_delta = 0
        m2_i = 0
        theo_i_lst = []

        isotope_score = 0
        isotope_m1_score = 0
        # obs_pr_mz = 0
        # ms1_theo_mz = isotope_pattern_df.at[0, 'mz']

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

                if not _i_df.empty:
                    _i_df = _i_df.sort_values(by='i', ascending=False).head(1)
                    _i_max = _i_df['i'].max()
                    _mz_max = _i_df['mz'].max()
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

        try:
            theo_pr_mz = isotope_checker_dct[0]['theo_mz']
            obs_pr_mz = isotope_checker_dct[0]['obs_mz']
            obs_pr_i = isotope_checker_dct[0]['obs_i']
        except Exception as _e:
            if mode == 'm':
                print(core_count, '[Exception] !!! Cannot get Theoretical or observed PR m/z ...', _e)
            else:
                print(core_count, '[Exception] !!! M+2 is NOT potential M+0 of M+2H ...', _e)
            theo_pr_mz = 0
            obs_pr_mz = 0
            obs_pr_i = 0

        if mode == 'm':
            if score_filter > 0:
                peak_top, top_obs_i = self.peak_top_checker(theo_pr_mz, spec_df,
                                                            core_count=core_count, ms1_precision=ms1_precision)
                if peak_top is True:
                    isotope_calc_dct = {'isotope_checker_dct': isotope_checker_dct, 'isotope_score': isotope_score,
                                        'isotope_m1_score': isotope_m1_score, 'theo_i_lst': theo_i_lst,
                                        'obs_pr_mz': obs_pr_mz, 'obs_pr_i': obs_pr_i}
                else:
                    isotope_calc_dct = {'isotope_checker_dct': {}, 'isotope_score': 0,
                                        'isotope_m1_score': 0, 'theo_i_lst': [],
                                        'obs_pr_mz': 0, 'obs_pr_i': 0}
            else:
                isotope_calc_dct = {'isotope_checker_dct': isotope_checker_dct, 'isotope_score': isotope_score,
                                    'isotope_m1_score': isotope_m1_score,  'theo_i_lst': theo_i_lst,
                                    'obs_pr_mz': obs_pr_mz, 'obs_pr_i': obs_pr_i}
        else:
            isotope_calc_dct = {'isotope_checker_dct': isotope_checker_dct, 'isotope_score': isotope_score,
                                'isotope_m1_score': isotope_m1_score, 'theo_i_lst': theo_i_lst,
                                'obs_pr_mz': obs_pr_mz, 'obs_pr_i': obs_pr_i}

        # print('mode', mode)
        # print(isotope_calc_dct)

        return isotope_calc_dct

    def get_deconvolution(self, spec_df, mz_delta, base_i, elem_dct={}, only_c=False):

        """
        :param elem_dct: (dict)
        :param spec_df:
        :param mz_delta:
        :param base_i:
        :param only_c:
        :return:
        """

        base_m1_i = 0
        base_m2_i = 0
        base_m3_i = 0

        m_pre1_isotope_pattern_df = self.get_isotope_mz(elem_dct, only_c=only_c, isotope_number=3)
        m_pre1_mz = m_pre1_isotope_pattern_df.at[0, 'mz']
        pre_i_df = spec_df.query('%f <= mz <= %f' %
                                 (m_pre1_mz - mz_delta, m_pre1_mz + mz_delta))

        if not pre_i_df.empty:
            max_m_pre1_i = pre_i_df['i'].max()
            max_m_pre1_i -= base_i
            if max_m_pre1_i > 0:
                base_m1_i += max_m_pre1_i * m_pre1_isotope_pattern_df.at[1, 'ratio']
                base_m2_i += max_m_pre1_i * m_pre1_isotope_pattern_df.at[2, 'ratio']
                base_m3_i += max_m_pre1_i * m_pre1_isotope_pattern_df.at[3, 'ratio']

        return base_m1_i, base_m2_i, base_m3_i

    def get_isotope_score(self, ms1_pr_mz, ms1_pr_i, formula, spec_df, core_count,
                          ms1_precision=50e-6, only_c=False, score_filter=75, decon=True, exp_mode='LC-MS',
                          isotope_number=2, pattern_tolerance=5):

        mz_delta = ms1_pr_mz * ms1_precision
        delta_13c = 1.0033548378
        obs_pr_mz = 0
        obs_pr_i = 0

        if exp_mode == 'Shotgun':
            pseudo_pr_check = 0
        else:
            pseudo_pr_check = 1

        if decon is True:
            deconv_elem_dct = self.get_elements(formula)
            # M-2
            deconv_elem_dct['H'] += -2
            base_i = 0
            pre2_base_m1_i, pre2_base_m2_i, pre2_base_m3_i = self.get_deconvolution(spec_df, mz_delta, base_i,
                                                                                    deconv_elem_dct, only_c=only_c)

            # M+0
            # m0_base_abs = pre2_base_m2_i + pre1_base_m1_i
            m0_base_abs = pre2_base_m2_i
            base_m1_i, base_m2_i, base_m3_i = self.get_deconvolution(spec_df, mz_delta, m0_base_abs,
                                                                     self.get_elements(formula), only_c=only_c)
            # M+1
            m1_base_abs = pre2_base_m3_i

            # M+2
            m2_base_abs = base_m2_i

            # M+3
            m3_base_abs = base_m3_i

            # print(core_count, 'Deconvolution_i_abs_corrections: [M+0] %.1f, [M+1] %.1f, [M+2] %.1f, [M+3] %.1f'
            #       % (m0_base_abs, m1_base_abs, m2_base_abs, m3_base_abs))
            deconv_lst = [m0_base_abs, m1_base_abs, m2_base_abs, m3_base_abs]
        else:
            m0_base_abs = 0
            m1_base_abs = 0
            m2_base_abs = 0
            m3_base_abs = 0

            deconv_lst = [m0_base_abs, m1_base_abs, m2_base_abs, m3_base_abs]

        i_df = spec_df.query('%f <= mz <= %f' % (ms1_pr_mz - delta_13c - mz_delta, ms1_pr_mz - delta_13c + mz_delta))

        if not i_df.empty:
            max_pre_m_i = i_df['i'].max()
            peak_top, top_obs_i = self.peak_top_checker(ms1_pr_mz - delta_13c, spec_df,
                                                        core_count=core_count, ms1_precision=ms1_precision)
            if peak_top is True:
                pass
            else:
                max_pre_m_i = top_obs_i
        else:
            # print(core_count, '... No pseudo-precursor peaks ... ')
            max_pre_m_i = 0

        if ms1_pr_i > max_pre_m_i or pseudo_pr_check == 0:
            elem_dct = self.get_elements(formula)
            mono_mz = self.get_mono_mz(elem_dct)
            if abs((ms1_pr_mz - mono_mz)) <= ms1_precision * ms1_pr_mz:
                isotope_pattern_df = self.get_isotope_mz(elem_dct, only_c=only_c)
                # print(isotope_pattern_df)
                ms1_pr_i -= m0_base_abs
                m0_deconv_lst = [m0_base_abs, m1_base_abs, m2_base_abs]
                isotope_calc_dct = self.calc_isotope_score(isotope_pattern_df, spec_df, ms1_precision, ms1_pr_i,
                                                           core_count=core_count, deconv=m0_deconv_lst, mode='m',
                                                           score_filter=score_filter)

                isotope_checker_dct = isotope_calc_dct['isotope_checker_dct']
                isotope_score = isotope_calc_dct['isotope_score']
                isotope_m1_score = isotope_calc_dct['isotope_m1_score']
                # m2_i = isotope_calc_dct['m2_i']
                m2_checker_dct = {}
                m2_score = 0

                obs_pr_mz = isotope_calc_dct['obs_pr_mz']
                obs_pr_i = isotope_calc_dct['obs_pr_i']

                if isotope_score <= score_filter <= isotope_m1_score:
                    print(core_count, '[INFO] --> check if M+2 is potential M+0 of M+2H ...')

                    # check if M+2 is potential M+0 of M+2H
                    # M+H2 elements

                    m2_elem_dct = self.get_elements(formula + 'H2')
                    m2_isotope_pattern_df = self.get_isotope_mz(m2_elem_dct, only_c=only_c)
                    # get exact M+2 pr i, especially for ppm < 10
                    m2_mz_df = m2_isotope_pattern_df.head(1)
                    m2_mz = m2_mz_df['mz'].tolist()[0]
                    m2_i_df = spec_df.query('%f <= mz <= %f' % (m2_mz * (1 - ms1_precision),
                                                                m2_mz * (1 + ms1_precision)))
                    m2_i = m2_i_df['i'].max()
                    m2_i -= m2_base_abs
                    m2_deconv_lst = [m2_base_abs, m3_base_abs, 0]
                    m2_calc_dct = self.calc_isotope_score(m2_isotope_pattern_df, spec_df, ms1_precision, m2_i,
                                                          core_count=core_count, deconv=m2_deconv_lst, mode='m+2',
                                                          score_filter=score_filter)

                    m2_checker_dct = m2_calc_dct['isotope_checker_dct']
                    # use M+1 only
                    m2_score = m2_calc_dct['isotope_m1_score']
                    # m2_m2i = m2_calc_dct['m2_i']

                    if m2_score > 0:
                        pass
                    else:
                        m2_score = 0

                    if 2 in list(m2_checker_dct.keys()):
                        del m2_checker_dct[2]

                    if m2_score >= max(score_filter * 0.75, 60):  # set lower filter for M+2H
                        print(core_count, '[INFO] --> M+2 ~ M+4 has isotope score for [M+H2]: %.1f' % m2_score)
                        isotope_score = isotope_m1_score
                        del isotope_checker_dct[2]['obs_i']
                        del isotope_checker_dct[2]['obs_mz']
                    else:
                        print(core_count, '[WARNING] !!! M+2 ~ M+4 has isotope score for [M+H2]: %.1f' % m2_score)
                        isotope_score = isotope_m1_score
                        m2_checker_dct = {}
                        m2_score = 0
                        del isotope_checker_dct[2]['obs_i']
                        del isotope_checker_dct[2]['obs_mz']

            else:
                print(core_count, '[WARNING] !!! MS1 PR m/z not fit to Formula check bulk identification !!!')
                isotope_score = 0
                isotope_checker_dct = {}
                m2_checker_dct = {}
                m2_score = 0

        else:
            print(core_count, '[WARNING] !!! MS1 PR is an isotope !!!')
            isotope_score = 0
            isotope_checker_dct = {}
            m2_checker_dct = {}
            m2_score = 0

        isotope_score_info_dct = {'isotope_score': isotope_score, 'isotope_checker_dct': isotope_checker_dct,
                                  'm2_score': m2_score, 'm2_checker_dct': m2_checker_dct,
                                  'deconv_lst': deconv_lst, 'obs_pr_mz': obs_pr_mz, 'obs_pr_i': obs_pr_i}

        return isotope_score_info_dct

    def get_isotope_fragments(self, ms1_pr_mz, ms1_pr_i, formula, spec_df, core_count,
                              ms1_precision=50e-6, only_c=False,
                              decon=True, exp_mode='LC-MS'):
        # TODO (georgia.angelidou@uni-leipzig.de): Need to check the reason why we do not get any output
        mz_delta = ms1_pr_mz * ms1_precision
        delta_13c = 1.0033548378

        if exp_mode == 'Shotgun':
            pseudo_pr_check = 0
        else:
            pseudo_pr_check = 1
        if decon is True:
            deconv_elem_dct = self.get_elements(formula)

            deconv_elem_dct['H'] += -2
            base_i = 0
            pre2_base_m1_i, pre2_base_m2_i, pre2_base_m3_i = self.get_deconvolution(spec_df, mz_delta, base_i,
                                                                                    deconv_elem_dct, only_c=only_c)
            m0_base_abs = pre2_base_m2_i
            base_m1_i, base_m2_i, base_m3_i = self.get_deconvolution(spec_df, mz_delta, m0_base_abs,
                                                                     self.get_elements(formula), only_c=only_c)

            m1_base_abs = pre2_base_m3_i

            # M+2
            m2_base_abs = base_m2_i

            # M+3
            m3_base_abs = base_m3_i

            print(core_count, '[INFO] --> Deconvolution_i_abs_corrections: '
                              '[M+0] %.1f, [M+1] %.1f, [M+2] %.1f, [M+3] %.1f'
                  % (m0_base_abs, m1_base_abs, m2_base_abs, m3_base_abs))

            deconv_lst = [m0_base_abs, m1_base_abs, m2_base_abs, m3_base_abs]
        else:
            m0_base_abs = 0
            m1_base_abs = 0
            m2_base_abs = 0
            m3_base_abs = 0

            deconv_lst = [m0_base_abs, m1_base_abs, m2_base_abs, m3_base_abs]

        i_df = spec_df.query('%f <= mz <= %f' % (ms1_pr_mz - delta_13c - mz_delta, ms1_pr_mz - delta_13c + mz_delta))
        isotope_flag = 0
        if not i_df.empty:
            max_pre_m_i = i_df['i'].max()

            if ms1_pr_i > max_pre_m_i or pseudo_pr_check == 0:
                elem_dct = self.get_elements(formula)
                mono_mz = self.get_mono_mz(elem_dct)
                if abs((ms1_pr_mz - mono_mz)) <= ms1_precision * ms1_pr_mz:
                    isotope_flag = 0
                else:
                    print('[WARNING] !!! MS1 PR m/z not fit to Formula check bulk identification !!!')
            else:
                print('[WARNING] !!! MS1 PR is an isotope !!!')
                isotope_flag = 1
        else:
            pass
        return isotope_flag

    def get_isotope_fragments_sec(self, ms1_pr_mz, ms1_pr_i, formula, spec_df, core_count, ms1_precision=50e-6,
                                  only_c=False, decon=True, exp_mode='LC-MS'):

        mz_delta = ms1_pr_mz * ms1_precision
        delta_13c = 1.0033548378

        if exp_mode == 'Shotgun':
            pseudo_pr_check = 0
        else:
            pseudo_pr_check = 1
        if decon is True:
            deconv_elem_dct = self.get_elements(formula)

            deconv_elem_dct['H'] += -2
            base_i = 0
            pre2_base_m1_i, pre2_base_m2_i, pre2_base_m3_i = self.get_deconvolution(spec_df, mz_delta, base_i,
                                                                                    deconv_elem_dct, only_c=only_c)
            m0_base_abs = pre2_base_m2_i
            base_m1_i, base_m2_i, base_m3_i = self.get_deconvolution(spec_df, mz_delta, m0_base_abs,
                                                                     self.get_elements(formula), only_c=only_c)

            m1_base_abs = pre2_base_m3_i

            # M+2
            m2_base_abs = base_m2_i

            # M+3
            m3_base_abs = base_m3_i

            print(core_count, '[INFO] --> Deconvolution_i_abs_corrections: '
                              '[M+0] %.1f, [M+1] %.1f, [M+2] %.1f, [M+3] %.1f'
                  % (m0_base_abs, m1_base_abs, m2_base_abs, m3_base_abs))

            deconv_lst = [m0_base_abs, m1_base_abs, m2_base_abs, m3_base_abs]
        else:
            m0_base_abs = 0
            m1_base_abs = 0
            m2_base_abs = 0
            m3_base_abs = 0

            deconv_lst = [m0_base_abs, m1_base_abs, m2_base_abs, m3_base_abs]

        i_df = spec_df.query('%f <= mz <= %f' % (ms1_pr_mz - delta_13c - mz_delta, ms1_pr_mz - delta_13c + mz_delta))
        isotope_flag = 0
        if not i_df.empty:
            max_pre_m_i = i_df['i'].max()

            if ms1_pr_i > max_pre_m_i or pseudo_pr_check == 0:
                elem_dct = self.get_elements(formula)
                mono_mz = self.get_mono_mz(elem_dct)
                if abs((ms1_pr_mz - mono_mz)) <= ms1_precision * ms1_pr_mz:
                    isotope_flag = 0
                else:
                    print(core_count, '[WARNING] !!! MS1 PR m/z not fit to Formula check bulk identification !!!')
            else:
                print(core_count, '[WARNING] !!! MS1 PR is an isotope !!!')
                isotope_flag = 1
        else:
            pass
        return isotope_flag


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

    # f = 'C41H71NO7P'  # PE
    f = 'C51H92O6Na'
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
        print(isotope_distribute.at[0, 'mz'])
