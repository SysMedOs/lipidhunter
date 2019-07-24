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

import logging
import os
import sys

import pandas as pd
import configparser


class HunterConfig:

    def __init__(self, file_groups, log_level=logging.DEBUG):

        logging.basicConfig(format='%(asctime)s-%(levelname)s - %(message)s',
                            datefmt='%b-%d@%H:%M:%S',
                            level=log_level)
        self.logger = logging.getLogger('log')

        self.opt_dct = {
            'lipid_class': 'lipid_class',
            'charge': 'charge_mode',
            'mzml': 'mzml_path_str',
            'xlsx': 'xlsx_output_path_str',
        }

        file_df = pd.read_excel(file_groups, index_col=0)
        self.file_dct = file_df.to_dict(orient='index')
        self.logger.debug(self.file_dct)

        self.file_abbr_lst = []
        self.abbr_group_dct = {}
        self.group_abbr_dct = {}
        for f in self.file_dct:
            self.file_abbr_lst.append(self.file_dct[f]['ABBR'])
            self.abbr_group_dct[self.file_dct[f]['ABBR']] = self.file_dct[f]['GROUP']
            if self.file_dct[f]['GROUP'] in self.group_abbr_dct:
                self.group_abbr_dct[self.file_dct[f]['GROUP']].append(self.file_dct[f]['ABBR'])
            else:
                self.group_abbr_dct[self.file_dct[f]['GROUP']] = [self.file_dct[f]['ABBR']]
        self.headers = ['lipid_class', 'Proposed_structures', 'DISCRETE_ABBR', 'Formula_neutral', 'Formula_ion',
                        'Charge', 'Lib_mz', 'AVG_RANK_SCORE', 'AVG_ISOTOPE_SCORE', 'AVG_ppm', 'AVG_RT', 'MIN_RT',
                        'MAX_RT', 'IDENT_FA', 'IDENT_HG', 'IDENT_COUNT']
        self.headers.extend(sorted(list(self.group_abbr_dct.keys())))
        self.headers.extend(sorted(list(self.abbr_group_dct.keys())))

    def check_file(self, file_path):
        is_file = False
        try:
            if os.path.isfile(file_path):
                is_file = True
                self.logger.debug(f'Load file: {file_path}')
            else:
                self.logger.error(f'Failed to load file: {file_path}')
        except IOError:
            self.logger.error(f'Failed to load file: {file_path}')

        return is_file

    def read_cfg(self, cfg):

        config = configparser.ConfigParser()
        config.read(cfg)

        if config.has_section('settings'):
            usr_cfg = 'settings'
        elif config.has_section('parameters'):
            usr_cfg = 'parameters'
        else:
            if config.has_section('default'):
                usr_cfg = 'default'
            else:
                usr_cfg = None

        cfg_dct = {}

        if usr_cfg is not None:
            options = config.options(usr_cfg)
            if self.opt_dct['xlsx'] in options:
                xlsx_path = config.get(usr_cfg, self.opt_dct['xlsx'])

                if self.check_file(xlsx_path):
                    self.logger.debug(xlsx_path)
                    cfg_dct['xlsx_path'] = xlsx_path
                    cfg_dct['xlsx'] = os.path.splitext(os.path.basename(xlsx_path))[0]

                    if self.opt_dct['lipid_class'] in options:
                        cfg_dct['lipid_class'] = config.get(usr_cfg, self.opt_dct['lipid_class'])

                    if self.opt_dct['charge'] in options:
                        cfg_dct['charge'] = config.get(usr_cfg, self.opt_dct['charge'])

                    if self.opt_dct['mzml'] in options:
                        mzml_path = config.get(usr_cfg, self.opt_dct['mzml'])
                        cfg_dct['mzml'] = os.path.splitext(os.path.basename(mzml_path))[0]

                else:
                    self.logger.error(f'Failed to load file: {xlsx_path}')

        if 'charge' in cfg_dct and 'xlsx' in cfg_dct:
            self.logger.info(cfg_dct)
            return cfg_dct

        else:
            return None

    def load_batch_cfg(self, cfg_lst, merge_table, rank_score=40, isotope_score=80):

        sum_cfg_dct = {}

        for cfg in cfg_lst:
            cfg_dct = self.read_cfg(cfg)
            try:
                sum_cfg_dct[cfg_dct['xlsx']] = cfg_dct
            except Exception as err:
                self.logger.error(err)

        cfg_df = pd.DataFrame(sum_cfg_dct).T

        sum_df = self.merge_xlsx(cfg_df)
        self.merge_features(sum_df, merge_table, rank_score=rank_score, isotope_score=isotope_score)

    def merge_xlsx(self, cfg_df):

        sum_df = pd.DataFrame()

        for i, r in cfg_df.iterrows():
            mzml = r['mzml']
            tmp_df = pd.read_excel(r['xlsx_path'])
            tmp_df['xlsx'] = r['xlsx']
            tmp_df['mzml'] = mzml
            if mzml in self.file_dct:
                mzml_abbr = self.file_dct[mzml]['ABBR']
                tmp_df[mzml_abbr] = 1
                tmp_df[self.file_dct[mzml]['GROUP']] = 1
                tmp_df['abs_ppm'] = tmp_df['ppm'].abs()

            tmp_df['lipid_class'] = r['lipid_class']

            sum_df = sum_df.append(tmp_df, sort=False)

            del tmp_df

        sum_df.sort_values(by=['Lib_mz', 'MS1_obs_mz', 'RANK_SCORE', 'ISOTOPE_SCORE', 'lipid_class'], inplace=True)
        sum_df.reset_index(drop=True, inplace=True)
        print(sum_df.head())

        return sum_df

    def merge_features(self, sum_df, merge_table, rank_score=40, isotope_score=80):

        unique_df = self.unique_features(sum_df, rank_score=rank_score, isotope_score=isotope_score)
        unique_df.to_excel(merge_table)

    def unique_features(self, sum_df, rank_score=40, isotope_score=80):

        unique_df = pd.DataFrame()

        unique_discrete_lst = sum_df['DISCRETE_ABBR'].unique().tolist()
        for discrete in unique_discrete_lst:
            tmp_df = sum_df[sum_df['DISCRETE_ABBR'] == discrete]
            tmp1_df = tmp_df.query(f'RANK_SCORE >= {rank_score} and ISOTOPE_SCORE >= {isotope_score}')
            tmp2_df = tmp1_df.sort_values(by=['RANK_SCORE', 'ISOTOPE_SCORE', 'abs_ppm'],
                                          ascending=[False, False, True])
            r_df = tmp2_df.head(1)  # type: pd.DataFrame()
            if not r_df.empty:
                self.logger.debug(tmp2_df[['RANK_SCORE', 'ISOTOPE_SCORE', 'abs_ppm']])
                r_df.at[:, 'AVG_RANK_SCORE'] = tmp2_df['RANK_SCORE'].mean()
                r_df.at[:, 'AVG_ISOTOPE_SCORE'] = tmp2_df['ISOTOPE_SCORE'].mean()
                r_df.at[:, 'AVG_ppm'] = tmp2_df['abs_ppm'].mean()
                r_df.at[:, 'AVG_RT'] = tmp2_df['MS2_scan_time'].mean()
                r_df.at[:, 'MIN_RT'] = tmp2_df['MS2_scan_time'].min()
                r_df.at[:, 'MAX_RT'] = tmp2_df['MS2_scan_time'].max()
                r_df.at[:, 'IDENT_FA'] = tmp2_df['#Observed_FA'].max()
                try:
                    r_df.at[:, 'IDENT_HG'] = tmp2_df['#Specific_peaks'].max()
                except KeyError:
                    r_df.at[:, 'IDENT_HG'] = 0
                for file in self.file_abbr_lst:
                    if tmp2_df[file].max() > 0:
                        r_df.at[:, file] = 1
                #         if self.abbr_group_dct[file] in r_df.columns.tolist():
                #             r_df[self.abbr_group_dct[file]] += 1
                #         else:
                #             r_df.at[:, self.abbr_group_dct[file]] = 1
                for gp in self.group_abbr_dct:
                    r_df.at[:, gp] = r_df[self.group_abbr_dct[gp]].sum(axis=1)
                unique_df = unique_df.append(r_df)

        unique_df.at[:, 'IDENT_COUNT'] = unique_df[list(self.group_abbr_dct.keys())].sum(axis=1)
        pre1_output_unique_df = unique_df[self.headers]
        decimal_lst = {'AVG_RANK_SCORE': 2, 'AVG_ISOTOPE_SCORE': 2, 'AVG_ppm': 2, 'AVG_RT': 2, 'MIN_RT': 2, 'MAX_RT': 2}
        pre2_output_unique_df = pre1_output_unique_df.round(decimal_lst)
        output_unique_df = pre2_output_unique_df.sort_values(by=['lipid_class', 'Lib_mz', 'Proposed_structures'])
        output_unique_df.reset_index(drop=True, inplace=True)

        return output_unique_df


if __name__ == '__main__':
    cfg_folder = r'D:\AdipoAtlas_Project\HunterCfg'
    usr_file_groups = r'../temp/file_groups.xlsx'

    # cfg_folder = r'D:\AdipoAtlas_Project\HunterCfg\IT_PL'
    # usr_file_groups = r'../temp/file_groups_PL_IT.xlsx'

    # cfg_folder = r'D:\AdipoAtlas_Project\HunterCfg\AquireX_TGDG'
    # usr_file_groups = r'../temp/file_groups_AquireX.xlsx'

    # cfg_folder = r'D:\AdipoAtlas_Project\HunterCfg\TGDG'
    # usr_file_groups = r'../temp/file_groups_TGDG.xlsx'

    usr_cfg_lst = [os.path.join(cfg_folder, f) for f in os.listdir(cfg_folder)
                   if os.path.isfile(os.path.join(cfg_folder, f))]

    hunter_cfg = HunterConfig(usr_file_groups)
    hunter_cfg.load_batch_cfg(usr_cfg_lst, merge_table=r'../temp/unique_PL_QEx.xlsx',
                              rank_score=30, isotope_score=80)

    print('FIN')
