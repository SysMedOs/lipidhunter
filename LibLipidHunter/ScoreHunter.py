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

from __future__ import division
from __future__ import print_function

import re

import pandas as pd

try:
    from LibLipidHunter.IsotopeHunter import IsotopeHunter
    from LibLipidHunter.PanelPlotter import plot_spectra
except ImportError:  # for python 2.7.14
    from IsotopeHunter import IsotopeHunter
    from PanelPlotter import plot_spectra


def get_specific_peaks(key_frag_dct, mz_lib, ms2_df, hg_ms2_ppm=100, vendor='waters', exp_mode='LC-MS'):
    ms2_max_i = ms2_df['i'].max()
    ms2_precision = hg_ms2_ppm * 0.000001
    target_frag_df = key_frag_dct['target_frag_df']
    target_nl_df = key_frag_dct['target_nl_df']
    other_frag_df = key_frag_dct['other_frag_df']
    other_nl_df = key_frag_dct['other_nl_df']

    _target_frag_df = pd.DataFrame()
    _target_nl_df = pd.DataFrame()
    _other_frag_df = pd.DataFrame()
    _other_nl_df = pd.DataFrame()

    for _i, _frag_se in target_frag_df.iterrows():

        _frag_mz = _frag_se['EXACTMASS']
        _frag_class = _frag_se['CLASS']
        _frag_label = _frag_se['LABEL']
        _frag_mz_query_code = '%f <= mz <= %f' % (_frag_mz * (1 - ms2_precision), _frag_mz * (1 + ms2_precision))

        _frag_df = ms2_df.query(_frag_mz_query_code)

        if _frag_df.shape[0] > 0:
            _frag_df = _frag_df.sort_values(by='i', ascending=False)
            _frag_df.loc[:, 'CLASS'] = _frag_class
            _frag_df.loc[:, 'LABEL'] = _frag_label
            _frag_df.loc[:, _frag_label] = 100 * _frag_df['i'] / ms2_max_i
            _target_frag_df = _target_frag_df.append(_frag_df.head(1))

    if other_frag_df.shape[0] > 0:
        for _i, _frag_se in other_frag_df.iterrows():

            _frag_mz = _frag_se['EXACTMASS']
            _frag_class = _frag_se['CLASS']
            _frag_label = _frag_se['LABEL']
            _frag_mz_query_code = '%f <= mz <= %f' % (_frag_mz * (1 - ms2_precision), _frag_mz * (1 + ms2_precision))

            _frag_df = ms2_df.query(_frag_mz_query_code)

            if _frag_df.shape[0] > 0:
                _frag_df = _frag_df.sort_values(by='i', ascending=False)
                _frag_df.loc[:, 'CLASS'] = _frag_class
                _frag_df.loc[:, 'LABEL'] = _frag_label
                _frag_df.loc[:, _frag_label] = 100 * _frag_df['i'] / ms2_max_i
                _other_frag_df = _other_frag_df.append(_frag_df.head(1))

    for _i, _nl_se in target_nl_df.iterrows():

        _nl_mz = _nl_se['EXACTMASS']
        _nl_class = _nl_se['CLASS']
        _nl_label = _nl_se['LABEL']
        _nl_mz_query_code = '%f <= mz <= %f' % ((mz_lib - _nl_mz) * (1 - ms2_precision),
                                                (mz_lib - _nl_mz) * (1 + ms2_precision))

        _nl_df = ms2_df.query(_nl_mz_query_code)

        if _nl_df.shape[0] > 0:
            _nl_df = _nl_df.sort_values(by='i', ascending=False)
            _nl_df.loc[:, 'CLASS'] = _nl_class
            _nl_df.loc[:, 'LABEL'] = _nl_label
            _nl_df.loc[:, _nl_label] = 100 * _nl_df['i'] / ms2_max_i
            _target_nl_df = _target_nl_df.append(_nl_df.head(1))

    for _i, _nl_se in other_nl_df.iterrows():

        _nl_mz = _nl_se['EXACTMASS']
        _nl_class = _nl_se['CLASS']
        _nl_label = _nl_se['LABEL']
        _nl_mz_query_code = '%f <= mz <= %f' % ((mz_lib - _nl_mz) * (1 - ms2_precision),
                                                (mz_lib - _nl_mz) * (1 + ms2_precision))
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


def get_all_fa_frag(fa_df, ms2_df, lipid_type):
    obs_peaks_df = pd.DataFrame()
    bp_i = ms2_df['i'].max()
    if lipid_type in ['TG']:
        fa_mode = '[FA-H2O+H]+'
    else:
        fa_mode = '[FA-H]-'
    # find all possible FA
    for _idx, _fa_se in fa_df.iterrows():
        # print(_fa_se)
        _q_str = _fa_se['[FA-H]-_Q']
        _q_tmp_df = ms2_df.query(_q_str).copy()
        _q_tmp_df.is_copy = False
        if _q_tmp_df.shape[0] > 0:
            _q_tmp_df.loc[:, 'lib_mz'] = _fa_se['%s_MZ' % fa_mode]
            _q_tmp_df.loc[:, 'obs_mz'] = _q_tmp_df['mz']
            _q_tmp_df.loc[:, 'obs_i_r'] = 100 * _q_tmp_df['i'] / bp_i
            _q_tmp_df.loc[:, 'obs_ppm'] = 1e6 * (_q_tmp_df['mz'] - _q_tmp_df['lib_mz']) / _q_tmp_df['lib_mz']
            _q_tmp_df.loc[:, 'obs_ppm'] = _q_tmp_df['obs_ppm'].astype(int)
            _q_tmp_df.loc[:, 'obs_ppm_abs'] = _q_tmp_df['obs_ppm'].abs()
            _q_tmp_df.loc[:, 'obs_abbr'] = _fa_se['%s_ABBR' % fa_mode]
            _q_tmp_df.loc[:, 'obs_label'] = _q_tmp_df['lib_mz']
            _q_tmp_df = _q_tmp_df.round({'obs_mz': 4, 'obs_i_r': 1, 'obs_label': 2})
            _q_tmp_df.loc[:, 'obs_label'] = _q_tmp_df['obs_label'].astype(str)
            obs_peaks_df = obs_peaks_df.append(_q_tmp_df)

    if obs_peaks_df.shape[0] > 0:
        obs_peaks_df.sort_values(by=['obs_abbr', 'i', 'obs_ppm_abs'], ascending=[False, False, True], inplace=True)
        obs_peaks_df.drop_duplicates(subset=['obs_abbr'], keep='first', inplace=True)
        obs_peaks_df.sort_values(by=['i', 'obs_ppm_abs'], ascending=[False, True], inplace=True)
        obs_peaks_df.reset_index(inplace=True, drop=True)
        obs_peaks_df['obs_rank'] = obs_peaks_df.index + 1
    else:
        print('Warning: get_all_fa_frag report no FRAG peak found !!!')
    return obs_peaks_df.head(10)


def get_all_fa_nl(fa_df, ms2_df, lyso_type_lst, lipid_type='LPL'):
    dg_fa_rgx = re.compile(r'\[M-(.*)\+H\]+')
    obs_peaks_df = pd.DataFrame()
    bp_i = ms2_df['i'].max()

    for _idx, _fa_se in fa_df.iterrows():
        for lyso_typ in lyso_type_lst:
            _q_str = _fa_se['%s_Q' % lyso_typ]
            _q_tmp_df = ms2_df.query(_q_str).copy()
            _q_tmp_df.is_copy = False
            if _q_tmp_df.shape[0] > 0:
                _q_tmp_df.loc[:, 'lib_mz'] = _fa_se['%s_MZ' % lyso_typ]
                _q_tmp_df.loc[:, 'obs_mz'] = _q_tmp_df['mz']
                _q_tmp_df.loc[:, 'obs_i_r'] = 100 * _q_tmp_df['i'] / bp_i
                _q_tmp_df.loc[:, 'obs_ppm'] = 1e6 * (_q_tmp_df['mz'] - _q_tmp_df['lib_mz']) / _q_tmp_df['lib_mz']
                _q_tmp_df.loc[:, 'obs_ppm'] = _q_tmp_df['obs_ppm'].astype(int)
                _q_tmp_df.loc[:, 'obs_ppm_abs'] = _q_tmp_df['obs_ppm'].abs()
                _q_tmp_df.loc[:, 'obs_abbr'] = _fa_se['%s_ABBR' % lyso_typ]
                _q_tmp_df.loc[:, 'obs_label'] = _q_tmp_df['lib_mz']
                _q_tmp_df = _q_tmp_df.round({'obs_mz': 4, 'obs_i_r': 1, 'obs_label': 2})
                _q_tmp_df.loc[:, 'obs_label'] = _q_tmp_df['obs_label'].astype(str)
                if lipid_type in ['TG'] and re.match(dg_fa_rgx, _fa_se['%s_ABBR' % lyso_typ]):
                    _fa_abbr_match = re.match(dg_fa_rgx, _fa_se['%s_ABBR' % lyso_typ])
                    _fa_abbr_lst = _fa_abbr_match.groups()
                    _q_tmp_df.loc[:, 'fa_abbr'] = _fa_abbr_lst[0]
                else:
                    _q_tmp_df.loc[:, 'fa_abbr'] = _idx

                obs_peaks_df = obs_peaks_df.append(_q_tmp_df)

    if obs_peaks_df.shape[0] > 0:
        obs_peaks_df.sort_values(by=['obs_abbr', 'i', 'obs_ppm_abs'], ascending=[False, False, True], inplace=True)
        obs_peaks_df.drop_duplicates(subset=['obs_abbr'], keep='first', inplace=True)
        obs_peaks_df.sort_values(by=['i', 'obs_ppm_abs'], ascending=[False, True], inplace=True)
        obs_peaks_df.reset_index(inplace=True, drop=True)
        obs_peaks_df['obs_rank'] = obs_peaks_df.index + 1
    else:
        print('Warning: get_all_fa_frag report no NL peak found !!!')

    return obs_peaks_df.head(10)


def get_all_frag(frag_calc_df, ms2_df, lipid_class, frag_type_lst):
    obs_peaks_df = pd.DataFrame()
    bp_i = ms2_df['i'].max()

    for _idx, _fa_se in frag_calc_df.iterrows():
        for frag_type in frag_type_lst:
            _q_str = _fa_se['%s_Q' % frag_type]
            _q_tmp_df = ms2_df.query(_q_str).copy()
            _q_tmp_df.is_copy = False
            if _q_tmp_df.share[0] > 0:
                _q_tmp_df.loc[:, 'lib_mz'] = _fa_se['%s_MZ' % frag_type]
                _q_tmp_df.loc[:, 'obs_mz'] = _q_tmp_df['mz']
                _q_tmp_df.loc[:, 'obs_i_r'] = 100 * _q_tmp_df['i'] / bp_i
                _q_tmp_df.loc[:, 'obs_ppm'] = 1e6 * (_q_tmp_df['mz'] - _q_tmp_df['lib_mz']) / _q_tmp_df['lib_mz']
                _q_tmp_df.loc[:, 'obs_ppm'] = _q_tmp_df['obs_ppm'].astype(int)
                _q_tmp_df.loc[:, 'obs_ppm_abs'] = _q_tmp_df['obs_ppm'].abs()
                _q_tmp_df.loc[:, 'obs_abbr'] = _fa_se['%s_ABBR' % frag_type]
                _q_tmp_df.loc[:, 'obs_label'] = _q_tmp_df['lib_mz']
                _q_tmp_df = _q_tmp_df.round({'obs_mz': 4, 'obs_i_r': 1, 'obs_label': 2})
                _q_tmp_df.loc[:, 'obs_label'] = _q_tmp_df['obs_label'].astype(str)
                obs_peaks_df = obs_peaks_df.append(_q_tmp_df)
    if obs_peaks_df.shape[0] > 0:
        obs_peaks_df.sort_values(by=['obs_abbr', 'i', 'obs_ppm_abs'], ascending=[False, False, True], inplace=True)
        obs_peaks_df.drop_duplicates(subset=['obs_abbr'], keep='first', inplace=True)
        obs_peaks_df.sort_values(by=['i', 'obs_ppm_abs'], ascending=[False, True], inplace=True)
        obs_peaks_df.reset_index(inplace=True, drop=True)
        obs_peaks_df['obs_rank'] = obs_peaks_df.index + 1
    else:
        print('Warning')

    return obs_peaks_df.head(10)


def get_rankscore(fa_df, master_info_df, abbr_bulk, charge, ms2_df, _ms2_idx, lipid_type, weight_dct, core_count,
                  rankscore_filter=27.5, all_sn=True):
    lite_info_df = master_info_df.query('BULK_ABBR == "%s" and spec_index == %f' % (abbr_bulk, _ms2_idx))

    lite_info_df.is_copy = False
    lite_info_df['RANK_SCORE'] = 0
    lite_info_df['sn1_found'] = 0
    lite_info_df['sn2_found'] = 0

    # TODO (georgia.angelidou@uni-leipzig.de): Check and rename unused dicts
    ident_peak_dct = {}
    ident_peak_dct2 = {}
    ident_peak_dct3 = {'discrete_abbr': [], 'obs_label': [], 'i': [], 'mz': [], 'obs_abbr': [], 'obs_rank_type': [],
                       'obs_rank': []}
    ident_peak_multIndex = [[], []]
    obs_dct = {}
    frag_lst_dg = []
    if lipid_type in ['PA', 'PE', 'PG', 'PI', 'PS'] and charge == '[M-H]-':
        frag_lst = ['[L%s-H]-' % lipid_type, '[L%s-H2O-H]-' % lipid_type]
        frag_lst_fa = ['[FA-H]-']
    elif lipid_type in ['PC'] and charge in ['[M+HCOO]-', '[M+CH3COO]-']:
        frag_lst = ['[L%s-H]-' % lipid_type, '[L%s-H2O-H]-' % lipid_type]
        frag_lst_fa = ['[FA-H]-']
    elif lipid_type in ['TG']:
        frag_lst_fa = ['[FA-H2O+H]+']
        frag_lst = ['[MG-H2O+H]+']
        frag_lst_dg = ['[M-(SN1-H2O)+H]+', '[M-(SN2-H2O)+H]+', '[M-(SN3-H2O)+H]+', '[M-(SN1)+H]+', '[M-(SN2)+H]+',
                       '[M-(SN3)+H]+']
    else:
        # TODO (zhixu.ni@uni-leipzig.de): consider more situations here
        frag_lst = ['[L%s-H]-' % lipid_type, '[L%s-H2O-H]-' % lipid_type]
        frag_lst_fa = ['[FA-H]-']

    # TODO (georgia.angelidou@uni-leipzig.de): put the MG and the FA in the same dataframe and the DG in another one
    # this way will avoid to many unecessary dataframes
    obs_fa_frag_df = get_all_fa_nl(fa_df, ms2_df, frag_lst_fa, lipid_type)
    obs_fa_nl_df = get_all_fa_nl(fa_df, ms2_df, frag_lst, lipid_type)
    obs_dg_frag_df = pd.DataFrame()
    # obs_fa_frag_df = get_all_fa_frag(fa_df, ms2_df, lipid_type)
    # obs_fa_frag_df = get_all_frag(fa_df, ms2_df, lipid_type, ['[FA-H]-'])
    # obs_fa_frag_df = get_all_frag(fa_df, ms2_df, lipid_type, ['[L%s-H]-' % lipid_type, '[L%s-H2O-H]-' % lipid_type])
    if obs_fa_frag_df.shape[0] + obs_fa_nl_df.shape[0] > 0 and lipid_type in ['PA', 'PE', 'PG', 'PI', 'PS', 'PC']:
        if charge == '[M-H]-':
            obs_dct = {'[FA-H]-': [obs_fa_frag_df, ['SN1_[FA-H]-', 'SN2_[FA-H]-']],
                       '[L%s-H]-' % lipid_type: [obs_fa_nl_df, ['[LPL(SN1)-H]-', '[LPL(SN2)-H]-']],
                       '[L%s-H2O-H]-' % lipid_type: [obs_fa_nl_df, ['[LPL(SN1)-H2O-H]-', '[LPL(SN2)-H2O-H]-']]}
        elif lipid_type in ['PC'] and charge in ['[M+HCOO]-', '[M+CH3COO]-']:
            frag_lst = ['[L%s-H]-' % lipid_type, '[L%s-H2O-H]-' % lipid_type]
            frag_lst_fa = ['[FA-H]-']
            obs_fa_frag_df = get_all_fa_nl(fa_df, ms2_df, frag_lst_fa, lipid_type)
            obs_fa_nl_df = get_all_fa_nl(fa_df, ms2_df, frag_lst, lipid_type)
            obs_dct = {'[FA-H]-': [obs_fa_frag_df, ['SN1_[FA-H]-', 'SN2_[FA-H]-']],
                       '[L%s-H]-' % lipid_type: [obs_fa_nl_df, ['[LPL(SN1)-H]-', '[LPL(SN2)-H]-']],
                       '[L%s-H2O-H]-' % lipid_type: [obs_fa_nl_df, ['[LPL(SN1)-H2O-H]-', '[LPL(SN2)-H2O-H]-']]}

        elif charge == '[M+H]+':
            pass
            # TODO(zhixu.ni@uni-leipzig.de): add support to positive mode
        else:
            pass
    elif lipid_type in ['TG']:
        lite_info_df['sn3_found'] = 0
        obs_dg_frag_df = get_all_fa_nl(lite_info_df, ms2_df, frag_lst_dg, lipid_type)
        # obs_dct = {'[FA-H2O+H]+': [obs_fa_frag_df, ['SN1_[FA-H2O+H]+', 'SN2_[FA-H2O+H]+', 'SN3_[FA-H2O+H]+']],
        #            '[MG-H2O+H]+': [obs_fa_nl_df, ['[MG(SN1)-H2O+H]+', '[MG(SN2)-H2O+H]+', '[MG(SN3)-H2O+H]+']],
        #            '[M-FA+H]+': [obs_dg_frag_df, ['[M-(SN1)+H]+', '[M-(SN2)+H]+', '[M-(SN3)+H]+']],
        #            '[M-(FA-H2O)+H]+': [obs_dg_frag_df,
        #                                ['[M-(SN1-H2O)+H]+', '[M-(SN2-H2O)+H]+', '[M-(SN3-H2O)+H]+']]}
        #####
        # there is an error with the [M-(SN1-H2O)+H]+ temporary deactivated
        obs_dct = {'[FA-H2O+H]+': [obs_fa_frag_df, ['SN1_[FA-H2O+H]+', 'SN2_[FA-H2O+H]+', 'SN3_[FA-H2O+H]+']],
                   '[MG-H2O+H]+': [obs_fa_nl_df, ['[MG(SN1)-H2O+H]+', '[MG(SN2)-H2O+H]+', '[MG(SN3)-H2O+H]+']],
                   '[M-FA+H]+': [obs_dg_frag_df, ['[M-(SN1)+H]+', '[M-(SN2)+H]+', '[M-(SN3)+H]+']]}
    else:
        print(core_count, 'Warning: No informative peak found !!!')

    if len(list(obs_dct.keys())) > 0:

        for obs_type in list(obs_dct.keys()):

            _obs_df = pd.DataFrame(obs_dct[obs_type][0])
            _obs_lst = obs_dct[obs_type][1]
            _obs_drop_idx = []

            for _obs in _obs_lst:

                # default will change in a letter step
                lite_info_df.loc[:, '%s_RANK' % _obs] = 10  # set to Rank 10 +1 , so the score will be 0
                lite_info_df.loc[:, '%s_WEIGHT' % _obs] = weight_dct['%s' % _obs]['Weight']

                _obs_df2 = pd.DataFrame(_obs_df)
                if _obs_df2.shape[0] == 0:
                    _obs_df2 = pd.DataFrame({'i': [], 'mz': [], 'lib_mz': [], 'obs_mz': [], 'obs_i_r': [],
                                             'obs_ppm': [], 'obs_ppm_abs': [], 'obs_abbr': [], 'obs_label': [],
                                             'fa_abbr': [], 'obs_rank': []})

                # Take by one each observation and then check for the the different position of the observed fragment
                for _idx, _lite_se in lite_info_df.iterrows():

                    _lipid_abbr_comp_dct = {}
                    _abbr = _lite_se['%s_ABBR' % _obs]
                    _lipid_abbr = _lite_se['DISCRETE_ABBR']
                    print(core_count, _lipid_abbr)
                    # Part that check the intensities and change the score
                    # TODO (georgia.angelidou@uni-leipzig.de): need to be supported also for phospholipids
                    if lipid_type in ['TG']:
                        _sn1_abbr = _lite_se['SN1_ABBR']
                        _sn2_abbr = _lite_se['SN2_ABBR']
                        _sn3_abbr = _lite_se['SN3_ABBR']

                        for _sn_abbr in [_sn1_abbr, _sn2_abbr, _sn3_abbr]:
                            if _sn_abbr in list(_lipid_abbr_comp_dct.keys()):
                                _lipid_abbr_comp_dct[_sn_abbr] = _lipid_abbr_comp_dct[_sn_abbr] + 1
                            else:
                                _lipid_abbr_comp_dct[_sn_abbr] = 1

                        if len(_lipid_abbr_comp_dct) == 2:

                            _fa_abbr = list(_lipid_abbr_comp_dct.keys())[list(_lipid_abbr_comp_dct.values()).index(2)]
                            if _fa_abbr in list(_obs_df2['fa_abbr']):
                                _position_df = _obs_df2.loc[_obs_df2['fa_abbr'] == _fa_abbr]
                                # _obs_df2.set_value(_position_df.index[0], 'i',
                                #                    _obs_df2.loc[_position_df.index[0], 'i'] / 2)
                                _obs_df2.at[_position_df.index[0], 'i'] = _obs_df2.loc[_position_df.index[0], 'i'] / 2
                                _obs_df2.sort_values(by=['i', 'obs_ppm_abs'], ascending=[False, True], inplace=True)
                                _obs_df2.reset_index(inplace=True, drop=True)
                                _obs_df2['obs_rank'] = _obs_df2.index + 1

                        elif len(_lipid_abbr_comp_dct) == 1:

                            _fa_abbr = list(_lipid_abbr_comp_dct.keys())[list(_lipid_abbr_comp_dct.values()).index(3)]
                            if _fa_abbr in list(_obs_df2['fa_abbr']):
                                _position_df = _obs_df2.loc[_obs_df2['fa_abbr'] == _fa_abbr]
                                # _obs_df2.set_value(_position_df.index[0], 'i',
                                #                    _obs_df2.loc[_position_df.index[0], 'i'] / 3)
                                _obs_df2.at[_position_df.index[0], 'i'] = _obs_df2.loc[_position_df.index[0], 'i'] / 2
                                _obs_df2.sort_values(by=['i', 'obs_ppm_abs'], ascending=[False, True], inplace=True)
                                _obs_df2.reset_index(inplace=True, drop=True)
                                _obs_df2['obs_rank'] = _obs_df2.index + 1
                        else:
                            pass

                        _sn_total_count = 0
                        _sn1_count = lite_info_df.loc[_idx, 'sn1_found']
                        if _sn1_abbr in list(_obs_df2['fa_abbr']) and _sn1_count == 0:
                            # lite_info_df.set_value(_idx, 'sn1_found', 1)
                            lite_info_df.at[_idx, 'sn1_found'] = 1
                            _sn_total_count = _sn_total_count + 1
                        _sn2_count = lite_info_df.loc[_idx, 'sn2_found']
                        if _sn2_abbr in list(_obs_df2['fa_abbr']) and _sn2_count == 0:
                            # lite_info_df.set_value(_idx, 'sn2_found', 1)
                            lite_info_df.at[_idx, 'sn2_found'] = 1
                            _sn_total_count = _sn_total_count + 1
                        _sn3_count = lite_info_df.loc[_idx, 'sn3_found']
                        if _sn3_abbr in list(_obs_df2['fa_abbr']) and _sn3_count == 0:
                            # lite_info_df.set_value(_idx, 'sn3_found', 1)
                            lite_info_df.at[_idx, 'sn3_found'] = 1
                            _sn_total_count = _sn_total_count + 1
                        if all_sn is True and _sn_total_count < 3:
                            # lite_info_df.set_value(_idx, 'sn3_found', 0)
                            lite_info_df.at[_idx, 'sn3_found'] = 0
                            # lite_info_df.set_value(_idx, 'sn2_found', 0)
                            lite_info_df.at[_idx, 'sn2_found'] = 0
                            # lite_info_df.set_value(_idx, 'sn1_found', 0)
                            lite_info_df.at[_idx, 'sn1_found'] = 0
                        else:
                            pass

                    else:
                        # TODO (georgia.angelidou@uni-leipzig.ge): to the above also for the phospholipids
                        pass

                    try:
                        if _abbr in _obs_df['obs_abbr'].values:
                            _rank_idx = _obs_df.loc[_obs_df['obs_abbr'] == _abbr].index[0]
                            _i = _obs_df.loc[_rank_idx, 'i']
                            _i_r = _obs_df.loc[_rank_idx, 'obs_i_r']
                            _mz = _obs_df.loc[_rank_idx, 'mz']
                            _label = _obs_df.loc[_rank_idx, 'obs_label']
                        else:
                            _rank_idx = 10
                            _i = 0
                            _i_r = 0
                            _mz = 0
                            _label = ''
                    except (IndexError, KeyError):
                        _rank_idx = 10
                        _i = 0
                        _i_r = 0
                        _mz = 0
                        _label = ''

                    if _rank_idx <= 10 and _i > 0:
                        lite_info_df.at[_idx, '%s_RANK' % _obs] = _rank_idx
                        lite_info_df.at[_idx, '%s_i' % _obs] = _i
                        lite_info_df.at[_idx, '{obs}_i_per'.format(obs=_obs)] = _i_r
                        # print(_abbr, _rank_idx, _i)
                        # ident_peak_dct[_abbr] = {'discrete_abbr': _lipid_abbr, 'obs_label': _label, 'i': _i,
                        #                          'mz': _mz, 'obs_abbr': _abbr, 'obs_rank_type': '%s_RANK' % _obs,
                        #                          'obs_rank': _rank_idx}
                        ident_peak_multIndex[0].append(_abbr)
                        ident_peak_multIndex[1].append(_lipid_abbr)

                        ident_peak_dct3['discrete_abbr'].append(_lipid_abbr)
                        ident_peak_dct3['obs_label'].append(_label)
                        ident_peak_dct3['i'].append(_i)
                        ident_peak_dct3['mz'].append(_mz)
                        ident_peak_dct3['obs_abbr'].append(_abbr)
                        ident_peak_dct3['obs_rank_type'].append('%s_RANK' % _obs)
                        ident_peak_dct3['obs_rank'].append(_rank_idx)

                        print(core_count, _abbr)
                        _obs_drop_idx.append(_rank_idx)

                    else:
                        pass
                        # print(_obs, _abbr, 'Not Found!')

                lite_info_df.loc[:, '%s_SCORE' % _obs] = ((10 - lite_info_df['%s_RANK' % _obs]) * 0.1
                                                          * lite_info_df['%s_WEIGHT' % _obs])
                lite_info_df.loc[:, 'RANK_SCORE'] += lite_info_df['%s_SCORE' % _obs]

            _obs_drop_idx = list(set(_obs_drop_idx))
            obs_dct[obs_type].append(_obs_drop_idx)

    if lipid_type in ['TG']:
        # TODO (georgia.angelidou@uni-leipzig.de): keep only the 2 decimenals
        lite_info_df['RANK_SCORE'] = round(lite_info_df['RANK_SCORE'] * (
                (lite_info_df['sn1_found'] + lite_info_df['sn2_found'] + lite_info_df['sn3_found']) / 3), 2)

    print(core_count, 'lite_info')

    if lite_info_df.shape[0] > 0:

        lite_info_df = lite_info_df[lite_info_df['RANK_SCORE'] >= rankscore_filter]
        if lipid_type in ['PA', 'PE', 'PG', 'PI', 'PS'] and charge == '[M-H]-':
            lite_info_df.loc[:, 'ident_rank'] = lite_info_df['SN1_[FA-H]-_RANK'] + lite_info_df['SN2_[FA-H]-_RANK']
            lite_info_df.sort_values(by=['RANK_SCORE', 'ident_rank'], ascending=[False, True], inplace=True)
        elif lipid_type in ['TG', 'DG', 'MG'] and charge == '[M+NH4]+':
            # lite_info_df.loc[:, 'ident_rank'] = lite_info_df['SN1_[FA-H2O+H]+_RANK'] + lite_info_df[
            #     'SN2_[FA-H2O+H]+_RANK'] + lite_info_df['SN3_[FA-H2O+H]+_RANK']
            # lite_info_df.sort_values(by=['RANK_SCORE', 'ident_rank'], ascending=[False, True], inplace=True)
            # print ('Need the TG')
            # works also without TG
            pass
            # TODO(zhixu.ni@uni-leipzig.de): @Georgia add TG here please :)
        lite_info_df.reset_index(drop=True, inplace=True)
        # ident_peak_df = pd.DataFrame(ident_peak_dct).T
        tuples = list(zip(*ident_peak_multIndex))
        ident_peak_index = pd.MultiIndex.from_tuples(tuples, names=['fragment_abbr', 'discrete_abbr'])
        ident_peak_df2 = pd.DataFrame(ident_peak_dct3, index=ident_peak_index)
    else:
        ident_peak_df = pd.DataFrame()
        ident_peak_df2 = pd.DataFrame()
    if lite_info_df.shape[0] > 0 and ident_peak_df2.shape[0] > 0:
        matched_checker = 1
        checked_abbr_lst = lite_info_df['DISCRETE_ABBR'].values.tolist()
        # ident_peak_df = ident_peak_df[ident_peak_df['discrete_abbr'].isin(checked_abbr_lst)]
        ident_peak_df2 = ident_peak_df2[ident_peak_df2['discrete_abbr'].isin(checked_abbr_lst)]
        # ident_peak_df.sort_values(by='mz', inplace=True)
        # ident_peak_df.reset_index(drop=True, inplace=True)
        ident_peak_df2.sort_values(by='mz', inplace=True)
        ident_peak_df2.reset_index(drop=True, inplace=True)

    else:
        matched_checker = 0

    # obs_info_dct = {'INFO': lite_info_df, 'OBS_FA': obs_fa_frag_df, 'OBS_LYSO': obs_dg_frag_df,
    #                 'IDENT': ident_peak_df, 'IDENT2': ident_peak_df2}
    if lipid_type in ['PA', 'PC', 'PE', 'PG', 'PI', 'PS']:
        obs_info_dct = {'INFO': lite_info_df, 'OBS_FA': obs_fa_frag_df, 'OBS_LYSO': obs_fa_nl_df,
                        'IDENT': ident_peak_df2}
    else:
        # TODO(georgia.angelidou@uni-leipzig.de): Reason why in the current output doesnt give DG
        # probably put in the OBS_FA the MG an free FA and OBS_LYSO put the FA loss
        # obs_info_dct = {'INFO': lite_info_df, 'OBS_FA': obs_fa_frag_df, 'OBS_LYSO': obs_fa_nl_df,
        #                 'IDENT': ident_peak_df, 'IDENT2' : ident_peak_df2}
        obs_fa_frag_df = obs_fa_frag_df.append(obs_fa_nl_df)
        if obs_fa_frag_df.shape[0] > 0:
            obs_fa_frag_df.sort_values(by='i', inplace=True)
            obs_fa_frag_df.reset_index(drop=True, inplace=True)

        obs_info_dct = {'INFO': lite_info_df, 'OBS_FA': obs_fa_frag_df, 'OBS_LYSO': obs_dg_frag_df,
                        'IDENT': ident_peak_df2}
    return matched_checker, obs_info_dct


def get_lipid_info(param_dct, fa_df, checked_info_df, checked_info_groups, core_list, usr_weight_df,
                   key_frag_dct, core_spec_dct, xic_dct, core_count, save_fig=True, os_type='windows', queue=None):
    core_count = 'Core_#%i' % core_count

    usr_lipid_type = param_dct['lipid_type']
    charge_mode = param_dct['charge_mode']
    output_folder = param_dct['img_output_folder_str']
    usr_ms2_th = param_dct['ms2_th']
    # usr_hg_th = param_dct['hg_th']
    usr_ms1_precision = param_dct['ms_ppm'] * 1e-6
    # usr_ms2_ppm = param_dct['ms2_ppm']
    usr_hg_ppm = param_dct['hg_ppm']
    usr_ms2_info_th_p = param_dct['ms2_infopeak_threshold']
    usr_hg_info_th_p = param_dct['ms2_hginfopeak_threshold']

    usr_isotope_score_filter = param_dct['isotope_score_filter']
    usr_rankscore_filter = param_dct['rank_score_filter']

    img_typ = param_dct['img_type']
    img_dpi = param_dct['img_dpi']
    usr_vendor = param_dct['vendor']
    exp_mode = param_dct['experiment_mode']
    usr_fast_isotope = param_dct['fast_isotope']
    usr_tag_all_sn = param_dct['tag_all_sn']

    hunter_start_time_str = param_dct['hunter_start_time']
    isotope_hunter = IsotopeHunter()

    # score_calc = ScoreGenerator(param_dct, usr_weight_df, usr_key_frag_df, usr_lipid_type,
    #                             checked_info_df, ion_charge=charge_mode, ms2_ppm=usr_ms2_ppm)

    usr_weight_dct = usr_weight_df.to_dict(orient='index')

    tmp_df = pd.DataFrame()

    img_plt_lst = []

    for group_key in core_list:
        _subgroup_df = checked_info_groups.get_group(group_key)
        _usr_abbr_bulk_lst = list(set(_subgroup_df['BULK_ABBR'].values.tolist()))
        usr_spec_info_dct = core_spec_dct[group_key]
        _samemz_se = _subgroup_df.iloc[0, :].squeeze()  # compress df to se for lipids with same bulk structures
        _usr_ms2_rt = _samemz_se['scan_time']
        _usr_ms2_dda_rank = _samemz_se['DDA_rank']
        _usr_ms2_scan_id = _samemz_se['scan_number']
        _usr_mz_lib = _samemz_se['Lib_mz']
        _usr_formula = _samemz_se['FORMULA']
        _usr_charge = _samemz_se['Ion']
        _usr_formula_charged = _samemz_se['%s_FORMULA' % _usr_charge]
        _usr_ms2_pr_mz = _samemz_se['Lib_mz']
        _ms1_pr_i = usr_spec_info_dct['ms1_i']
        _ms1_pr_mz = usr_spec_info_dct['ms1_mz']
        _ms1_df = usr_spec_info_dct['ms1_df']
        _ms2_df = usr_spec_info_dct['ms2_df']
        _ms2_idx = usr_spec_info_dct['_ms2_spec_idx']

        print(core_count, _usr_ms2_rt, _ms1_pr_mz, _usr_formula_charged)

        # use the max threshold from abs & relative intensity settings
        if 'i' in _ms2_df.columns.values.tolist():
            _ms2_max_i = _ms2_df['i'].max()
            ms2_threshold = max(usr_ms2_th, _ms2_max_i * usr_ms2_info_th_p)
            ms2_hg_threshold = max(usr_ms2_th, _ms2_max_i * usr_hg_info_th_p)
            _score_ms2_df = _ms2_df.query('i > %f' % ms2_threshold)
            if ms2_hg_threshold != ms2_threshold:
                _score_ms2_hg_df = _ms2_df.query('i > %f' % ms2_hg_threshold)
            else:
                _score_ms2_hg_df = _score_ms2_df
        else:
            _ms2_df = pd.DataFrame()
            _score_ms2_df = pd.DataFrame()
            _score_ms2_hg_df = pd.DataFrame()
        if _ms1_pr_mz > 0.0 and _ms1_df.shape[0] > 0 and _ms2_df.shape[0] > 0 and _ms1_pr_i > 0.0:

            print(core_count, '>>> >>> >>> >>> Best PR on MS1: %f' % _ms1_pr_mz)

            isotope_score_info_dct = isotope_hunter.get_isotope_score(_ms1_pr_mz, _ms1_pr_i,
                                                                      _usr_formula_charged, _ms1_df, core_count,
                                                                      ms1_precision=usr_ms1_precision,
                                                                      isotope_number=2,
                                                                      only_c=usr_fast_isotope,
                                                                      score_filter=usr_isotope_score_filter)

            isotope_score = isotope_score_info_dct['isotope_score']

            print(core_count, 'isotope_score: %f' % isotope_score)
            if isotope_score >= usr_isotope_score_filter:
                print(core_count, '>>> isotope_check PASSED! >>> >>> >>>')
                print(core_count, '>>> >>> >>> >>> Entry Info >>> >>> >>> >>> ')
                _samemz_se.at['MS1_obs_mz'] = _ms1_pr_mz
                _exact_ppm = 1e6 * (_ms1_pr_mz - _usr_mz_lib) / _usr_mz_lib
                _samemz_se.at['ppm'] = _exact_ppm
                _samemz_se.at['abs_ppm'] = abs(_exact_ppm)
                print(core_count, 'Proposed_bulk_structure can be:', _usr_abbr_bulk_lst)
                for _usr_abbr_bulk in _usr_abbr_bulk_lst:
                    print(core_count, 'Now check_proposed_structure:', _usr_abbr_bulk)

                    matched_checker, obs_info_dct = get_rankscore(fa_df, checked_info_df, _usr_abbr_bulk, charge_mode,
                                                                  _score_ms2_df, _ms2_idx, usr_lipid_type,
                                                                  usr_weight_dct, core_count,
                                                                  rankscore_filter=usr_rankscore_filter,
                                                                  all_sn=usr_tag_all_sn)

                    obs_info_df = obs_info_dct['INFO']
                    rank_score = obs_info_df['RANK_SCORE'].values.tolist()

                    if matched_checker > 0:
                        if len(key_frag_dct) > 0:
                            specific_dct = get_specific_peaks(key_frag_dct, _usr_mz_lib, _score_ms2_hg_df,
                                                              hg_ms2_ppm=usr_hg_ppm, vendor=usr_vendor,
                                                              exp_mode=exp_mode)
                        else:
                            specific_dct = {}

                        print(core_count, 'Rank_score: --> passed', rank_score)
                        print(core_count, obs_info_df[['BULK_ABBR', 'DISCRETE_ABBR', 'RANK_SCORE', 'scan_time']])

                        # format abbr. for file names
                        _save_abbr_bulk = _usr_abbr_bulk
                        _save_abbr_bulk = _save_abbr_bulk.replace(r'(', r'[')
                        _save_abbr_bulk = _save_abbr_bulk.replace(r')', r']')
                        _save_abbr_bulk = _save_abbr_bulk.replace(r'<', r'[')
                        _save_abbr_bulk = _save_abbr_bulk.replace(r'>', r']')
                        _save_abbr_bulk = _save_abbr_bulk.replace(r':', r'-')
                        _save_abbr_bulk = _save_abbr_bulk.replace(r'@', r'-')
                        _save_abbr_bulk = _save_abbr_bulk.replace('\\', r'_')
                        _save_abbr_bulk = _save_abbr_bulk.replace(r'/', r'_')

                        img_name_core = ('/%.4f_rt%.3f_DDAtop%.0f_scan%.0f_%s.%s'
                                         % (_usr_ms2_pr_mz, _usr_ms2_rt, _usr_ms2_dda_rank,
                                            _usr_ms2_scan_id, _save_abbr_bulk, img_typ)
                                         )
                        img_name = (output_folder +
                                    r'/LipidHunter_Results_Figures_%s'
                                    % hunter_start_time_str + img_name_core)

                        print(core_count, '==> check for output -->')

                        obs_info_df['Proposed_structures'] = _usr_abbr_bulk
                        obs_info_df['Bulk_identification'] = _usr_abbr_bulk
                        # obs_info_df['Discrete_identification'] = _usr_abbr_bulk
                        obs_info_df['Formula_neutral'] = _usr_formula
                        obs_info_df['Formula_ion'] = _usr_formula_charged
                        obs_info_df['Charge'] = _usr_charge
                        obs_info_df['MS1_obs_mz'] = _ms1_pr_mz
                        obs_info_df['MS1_obs_i'] = '%.2e' % float(_ms1_pr_i)
                        obs_info_df['Lib_mz'] = _usr_mz_lib
                        obs_info_df['MS2_scan_time'] = _usr_ms2_rt
                        obs_info_df['DDA#'] = _usr_ms2_dda_rank
                        obs_info_df['MS2_PR_mz'] = _usr_ms2_pr_mz
                        obs_info_df['Scan#'] = _usr_ms2_scan_id
                        # obs_info_df['#Specific_peaks'] = (target_frag_count + target_nl_count)
                        # obs_info_df['#Contaminated_peaks'] = (other_frag_count + other_nl_count)
                        obs_info_df['ppm'] = _exact_ppm

                        obs_info_df['img_name'] = img_name_core[1:]

                        tmp_df = tmp_df.append(obs_info_df)
                        if save_fig is True:
                            img_param_dct = {'abbr': _usr_abbr_bulk, 'mz_se': _samemz_se, 'xic_dct': xic_dct,
                                             'ident_info_dct': obs_info_dct, 'spec_info_dct': usr_spec_info_dct,
                                             'isotope_score_info_dct': isotope_score_info_dct,
                                             'specific_dct': specific_dct, 'formula_charged': _usr_formula_charged,
                                             'charge': _usr_charge, 'save_img_as': img_name}

                            img_plt_lst.append(img_param_dct.copy())
                        else:
                            pass

    if tmp_df.shape[0] > 0:
        print(core_count, 'Size of the identified LPP_df %i, %i' % (tmp_df.shape[0], tmp_df.shape[1]))
        tmp_df.reset_index(drop=True, inplace=True)
        tmp_df.index += 1

    else:
        print(core_count, '!! Size of the identified LPP_df == 0')
        tmp_df = 'error'
    if save_fig is True:
        print('img_plt_lst', len(img_plt_lst))
    else:
        img_plt_lst = []

    r_lst = (tmp_df, img_plt_lst)

    if os_type == 'linux_multi':
        queue.put(r_lst)
    else:
        return r_lst
