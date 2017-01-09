# -*- coding: utf-8 -*-
# Copyright 2015-2016 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

import re

import pandas as pd


class ScoreGenerator:
    def __init__(self, fa_def_df, weight_df):
        self.fa_def_df = fa_def_df
        self.weight_df = weight_df

    def get_fa_search(self, abbr, charge_type, mz_lib, ms2_df, ms2_precision=500e-6, ms2_threshold=100):

        fa_ident_df = pd.DataFrame()
        lyso_ident_df = pd.DataFrame()
        lyso_w_ident_df = pd.DataFrame()

        if charge_type in ['[M-H]-', '[M+HCOO]-']:
            charge_mode = 'NEG'
        elif charge_type in ['[M+H]+', '[M+Na]+', '[M+NH4]+', '[M+K]+']:
            charge_mode = 'POS'
        else:
            charge_mode = 'NEG'

        if abbr[:2] in ['PA', 'PC', 'PE', 'PG', 'PI', 'PS', 'SM']:
            lipid_type = 'PL'
        elif abbr[:2] in ['TA', 'TG', 'DA', 'DG', 'MA', 'MG']:
            lipid_type = 'GL'
        else:
            lipid_type = 'PL'

        if lipid_type == 'PL' and charge_mode == 'NEG':
            fa_chk_df = self.fa_def_df[['FA', 'Link', 'C', 'DB', 'mass', '[M-H]-', 'NL-H2O']]
            fa_chk_df = fa_chk_df.rename(columns={'[M-H]-': 'sn', 'mass': 'NL'})
            fa_chk_df['M-sn'] = mz_lib - fa_chk_df['NL-H2O']
            fa_chk_df['M-(sn-H2O)'] = mz_lib - fa_chk_df['NL']

            for _i, _fa_se in fa_chk_df.iterrows():

                _fa_abbr = _fa_se['FA']

                for _frag_type in ['sn', 'M-sn', 'M-(sn-H2O)']:
                    _frag_mz = _fa_se[_frag_type]
                    _frag_mz_low = _frag_mz - _frag_mz * ms2_precision
                    _frag_mz_high = _frag_mz + _frag_mz * ms2_precision
                    _frag_mz_query_code = '%f <= mz <= %f' % (_frag_mz_low, _frag_mz_high)

                    _frag_df = ms2_df.query(_frag_mz_query_code)

                    if _frag_df.shape[0] == 1:
                        _frag_df.loc[:, 'ppm'] = 1e6 * (_frag_df['mz'] - _frag_mz) / _frag_mz
                        _frag_df['ppm_abs'] = _frag_df['ppm'].abs()
                        _frag_df['FA'] = _fa_abbr

                    if _frag_df.shape[0] > 1:
                        _frag_df.loc[:, 'ppm'] = 1e6 * (_frag_df['mz'] - _frag_mz) / _frag_mz
                        _frag_df['ppm_abs'] = _frag_df['ppm'].abs()
                        _frag_df['FA'] = _fa_abbr
                        _frag_i_df = _frag_df.sort_values(by='i', ascending=False).head(1)
                        _frag_ppm_df = _frag_df.sort_values(by='ppm_abs').head(1)
                        _frag_df = _frag_i_df
                        if _frag_ppm_df['i'].tolist() == _frag_i_df['i'].tolist():
                            pass
                        else:
                            _frag_df = _frag_df.append(_frag_ppm_df)

                    if _frag_type == 'sn':
                        fa_ident_df = fa_ident_df.append(_frag_df)
                    elif _frag_type == 'M-sn':
                        lyso_ident_df = lyso_ident_df.append(_frag_df)
                    elif _frag_type == 'M-(sn-H2O)':
                        lyso_w_ident_df = lyso_ident_df.append(_frag_df)

        fa_ident_df = fa_ident_df.query('i > %f' % ms2_threshold)
        lyso_ident_df = lyso_ident_df.query('i > %f' % ms2_threshold)
        lyso_w_ident_df = lyso_w_ident_df.query('i > %f' % ms2_threshold)
        fa_ident_df = fa_ident_df[['FA', 'mz', 'i', 'ppm']].reset_index(drop=True)
        lyso_ident_df = fa_ident_df[['FA', 'mz', 'i', 'ppm']].reset_index(drop=True)
        lyso_w_ident_df = fa_ident_df[['FA', 'mz', 'i', 'ppm']].reset_index(drop=True)
        fa_ident_df = fa_ident_df.sort_values(by='i', ascending=False).head(10)
        lyso_ident_df = lyso_ident_df.sort_values(by='i', ascending=False).head(5)
        lyso_w_ident_df = lyso_w_ident_df.sort_values(by='i', ascending=False).head(5)

        return fa_ident_df, lyso_ident_df, lyso_w_ident_df

    def get_structure(self, abbr):

        lipid_abbr_lst = []
        lipid_sn1_lst = []
        lipid_sn2_lst = []
        db_sn1_lst = []
        db_sn2_lst = []

        pl_checker = re.compile(r'(P[ACEGSI])([(])(.*)([)])')
        pip_checker = re.compile(r'(PIP)([(])(.*)([)])')
        tg_checker = re.compile(r'(TG)([(])(.*)([)])')
        fa_checker = re.compile(r'(\d{1,2})([:])(\d)')
        fa_o_checker = re.compile(r'(O-)(\d{1,2})([:])(\d)')
        fa_p_checker = re.compile(r'(P-)(\d{1,2})([:])(\d)')

        # Check PL Type
        _pl_typ = ''
        bulk_fa_typ = ''
        bulk_fa_linker = ''
        bulk_fa_c = 0
        bulk_fa_db = 0

        if pl_checker.match(abbr):
            print ('PL')
            pl_re_chk = pl_checker.match(abbr)
            pl_typ_lst = pl_re_chk.groups()
            _pl_typ = pl_typ_lst[0]
            bulk_fa_typ = pl_typ_lst[2]
        if pip_checker.match(abbr):
            print ('PIP')
            pip_re_chk = pip_checker.match(abbr)
            pip_typ_lst = pip_re_chk.groups()
            _pl_typ = pip_typ_lst[0]
            bulk_fa_typ = pip_typ_lst[2]
        if tg_checker.match(abbr):
            print ('TG')
            pip_re_chk = pip_checker.match(abbr)
            pip_typ_lst = pip_re_chk.groups()
            _pl_typ = pip_typ_lst[0]
            bulk_fa_typ = pip_typ_lst[2]

        print(bulk_fa_typ)

        if fa_checker.match(bulk_fa_typ):
            bulk_fa_linker = 'A-A-'
            fa_chk = fa_checker.match(bulk_fa_typ)
            bulk_fa_lst = fa_chk.groups()
            bulk_fa_c = bulk_fa_lst[0]
            bulk_fa_db = bulk_fa_lst[2]
        elif fa_o_checker.match(bulk_fa_typ):
            bulk_fa_linker = 'O-A-'
            fa_chk = fa_o_checker.match(bulk_fa_typ)
            bulk_fa_lst = fa_chk.groups()
            bulk_fa_c = bulk_fa_lst[1]
            bulk_fa_db = bulk_fa_lst[3]
        elif fa_p_checker.match(bulk_fa_typ):
            bulk_fa_linker = 'P-A-'
            fa_chk = fa_p_checker.match(bulk_fa_typ)
            bulk_fa_lst = fa_chk.groups()
            bulk_fa_c = bulk_fa_lst[1]
            bulk_fa_db = bulk_fa_lst[3]

        bulk_fa_c = int(bulk_fa_c)
        bulk_fa_db = int(bulk_fa_db)

        for _i, _fa_se in self.fa_def_df.iterrows():
            # FA, Link, C, DB
            _fa_abbr = _fa_se['FA']
            _fa_link = _fa_se['Link']
            _fa_c = _fa_se['C']
            _fa_db = _fa_se['DB']

            if _fa_db <= bulk_fa_db and _fa_c <= bulk_fa_c:
                if _fa_link == bulk_fa_linker[0:1]:
                    _rest_fa_link = bulk_fa_linker[2]
                    _rest_fa_c = bulk_fa_c - _fa_c
                    _rest_fa_db = bulk_fa_db - _fa_db

                    _rest_fa_df = self.fa_def_df.query('Link == "%s" and C == %i and DB == %i'
                                                       % (_rest_fa_link, _rest_fa_c, _rest_fa_db)
                                                       )

                    if _rest_fa_df.shape[0] == 1:
                        # if _fa_link == 'A-':
                        #     _fa_link = ''
                        # else:
                        #     pass
                        _rest_fa_abbr = _rest_fa_df['FA'].tolist()[0]
                        lipid_abbr = '%s(%s_%s)' % (_pl_typ, _fa_abbr, _rest_fa_abbr)

                        lipid_abbr_lst.append(lipid_abbr)
                        lipid_sn1_lst.append(_fa_abbr)
                        lipid_sn2_lst.append(_rest_fa_abbr)
                        db_sn1_lst.append(_fa_db)
                        db_sn2_lst.append(_rest_fa_db)

        lipid_abbr_df = pd.DataFrame(data={'Lipid_abbr': lipid_abbr_lst, 'sn1_abbr': lipid_sn1_lst,
                                           'sn2_abbr': lipid_sn2_lst, 'sn1_DB': db_sn1_lst, 'sn2_DB': db_sn2_lst})

        lipid_abbr_df = lipid_abbr_df.query('sn1_DB <=sn2_DB')
        lipid_abbr_df = lipid_abbr_df[['Lipid_abbr', 'sn1_abbr', 'sn2_abbr']]

        return lipid_abbr_df

    def get_match(self, abbr, charge_type, mz_lib, ms2_df, ms2_precision=500e-6, ms2_threshold=100):

        fa_ident_df, lyso_ident_df, lyso_w_ident_df = self.get_fa_search(abbr, charge_type, mz_lib, ms2_df,
                                                                         ms2_precision=ms2_precision,
                                                                         ms2_threshold=ms2_threshold
                                                                         )

        print(fa_ident_df)
        print(lyso_ident_df)
        print(lyso_w_ident_df)

        lipid_abbr_df = self.get_structure(abbr)

        weight_type_lst = ['sn1', 'sn2', 'M-sn1', 'M-sn2', 'M-(sn1-H2O)', 'M-(sn2-H2O)']
        weight_dct = {}
        for _type in weight_type_lst:
            lipid_abbr_df[_type] = 0
        # lipid_abbr_df['Score'] = 0
        if fa_ident_df.shape[0] > 0:
            fa_ident_lst = fa_ident_df['FA'].tolist()

            try:
                lyso_ident_lst = lyso_ident_df['FA'].tolist()
            except KeyError:
                lyso_ident_lst = []
            try:
                lyso_w_ident_lst = lyso_w_ident_df['FA'].tolist()
            except KeyError:
                lyso_w_ident_lst = []

            self.weight_df['mz'] = 0.0
            for _i, _weight_se in self.weight_df.iterrows():
                _type = _weight_se['Type']
                _weight = _weight_se['Weight']
                # _pr_h = _weight_se['[M-H]-']
                # _pr_formate = _weight_se['[M+HCOO]-']
                # _fa = _weight_se['FA']
                # _water = _weight_se['H2O']
                weight_dct[_type] = _weight

            for _i_abbr, _abbr_se in lipid_abbr_df.iterrows():
                # _pl_abbr = _abbr_se['Lipid_abbr']
                _sn1_abbr = _abbr_se['sn1_abbr']
                _sn2_abbr = _abbr_se['sn2_abbr']

                if _sn1_abbr in fa_ident_lst:
                    lipid_abbr_df.set_value(_i_abbr, 'sn1',
                                            weight_dct['sn1'] * (10 - fa_ident_lst.index(_sn1_abbr)) / 10)
                if _sn2_abbr in fa_ident_lst:
                    lipid_abbr_df.set_value(_i_abbr, 'sn2',
                                            weight_dct['sn2'] * (10 - fa_ident_lst.index(_sn2_abbr)) / 10)
                if _sn1_abbr in lyso_ident_lst:
                    lipid_abbr_df.set_value(_i_abbr, 'M-sn1',
                                            weight_dct['M-sn1'] * (10 - lyso_ident_lst.index(_sn1_abbr)) / 10)
                if _sn2_abbr in lyso_ident_lst:
                    lipid_abbr_df.set_value(_i_abbr, 'M-sn2',
                                            weight_dct['M-sn2'] * (10 - lyso_ident_lst.index(_sn2_abbr)) / 10)
                if _sn1_abbr in lyso_w_ident_lst:
                    lipid_abbr_df.set_value(_i_abbr, 'M-(sn1-H2O)',
                                            weight_dct['M-(sn1-H2O)'] * (10 - lyso_w_ident_lst.index(_sn1_abbr)) / 10)
                if _sn2_abbr in lyso_w_ident_lst:
                    lipid_abbr_df.set_value(_i_abbr, 'M-(sn2-H2O)',
                                            weight_dct['M-(sn2-H2O)'] * (10 - lyso_w_ident_lst.index(_sn2_abbr)) / 10)

                lipid_abbr_df['Score'] = lipid_abbr_df[weight_type_lst].sum(axis=1, numeric_only=True)
        else:
            print('!!!!!! NO FA identified =====>--> Skip >>> >>>')
        return lipid_abbr_df, fa_ident_df, lyso_ident_df, lyso_w_ident_df


if __name__ == '__main__':
    fa_list_csv = r'D:\LPPhunter\FA_list.csv'
    score_cfg = r'D:\LPPhunter\Score_cfg.xlsx'

    usr_fa_def_df = pd.read_csv(fa_list_csv)
    usr_weight_df = pd.read_excel(score_cfg)
    usr_fa_def_df['C'] = usr_fa_def_df['C'].astype(int)
    usr_fa_def_df['DB'] = usr_fa_def_df['DB'].astype(int)

    # for usr_abbr in ['PC(O-36:3)', 'PC(P-38:4)', 'PC(40:5)']:
    #     struc_df = get_structure(usr_fa_def_df, usr_abbr)
    #     print(struc_df)

    usr_lp = 'PC(40:5)'
    score_calc = ScoreGenerator(usr_fa_def_df, usr_weight_df)
    struc_df = score_calc.get_structure(usr_lp)
    print(struc_df)

    usr_fa_ident_df = usr_fa_def_df.head(8)
    usr_lyso_ident_df = usr_fa_def_df.head(4)
    usr_lyso_w_ident_df = usr_fa_def_df.head(2)

    # usr_weight_df = score_calc.get_match(usr_lp, '[M+HCOO]-', mz_lib, ms2_df, ms2_precision=500e-6)
    # print(usr_weight_df)

    # fa_indicator = FAindicator(fa_list_csv)
    # lyso_indicator = Lyso_indicator(fa_list_csv)
    #
    # (fa_df, fa_info_dct) = fa_indicator.indicate(ms2_df)
    #
    # if fa_info_dct is not None:
    #     # print fa_info
    #     fa_count = len(fa_info_dct['mz'])
    #     fa_count_lst = range(fa_count)
    #     _fa_name_lst = fa_info_dct['name']
    #     _fa_mz_lst = fa_info_dct['mz']
    #     _fa_i_lst = fa_info_dct['i']
    #     _fa_ppm_lst = fa_info_dct['ppm']
    #     _fa_delta_lst = fa_info_dct['D']
    #
    #     fa_info_df = pd.DataFrame(fa_info_dct, columns=['name', 'mz', 'i', 'ppm'])
    #     fa_info_df = fa_info_df.sort_values(by='i', ascending=False)
    #     fa_info_df = fa_info_df.round({'mz': 4})
    #     # fa_info_df.index = range(1, len(fa_info['mz']) + 1)
    #     table_info_df = fa_info_df
    #
    #     # print table_info_df.head(5)
    #     old_idx_lst = table_info_df.index.tolist()
    #     table_info_df.index = range(1, len(old_idx_lst) + 1)
    #
    #     col_labels = table_info_df.columns.tolist()
    #     # col_labels = mz_info_df.head().tolist()
    #     row_labels = table_info_df.index.tolist()
    #     table_vals = map(list, table_info_df.values)
    #     try:
    #         # the rectangle is where I want to place the table
    #         fa_table = msms_pic.table(cellText=table_vals, rowLabels=row_labels,
    #                                   colWidths=[.10] * len(col_labels),
    #                                   colLabels=col_labels, loc='upper center')
    #         fa_table.set_fontsize(5)
    #         _auto_ident_chker += 1
    #     except IndexError:
    #         pass
    #         # table_props = the_table.properties()
    #         # table_cells = table_props['child_artists']
    #         # for cell in table_cells:
    #         #     cell.set_height(0.12)
    #
    # # get Lyso identification
    # (lyso_df, lyso_info_dct) = lyso_indicator.indicate(ms2_df, pr_mz, PLtype=_usr_pl_class)
    #
    # if lyso_info_dct is not None:
    #     # print lyso_info
    #     lyso_count = len(lyso_info_dct['mz'])
    #     lyso_count_lst = range(lyso_count)
    #     _lyso_name_lst = lyso_info_dct['name']
    #     _lyso_mz_lst = lyso_info_dct['mz']
    #     _lyso_i_lst = lyso_info_dct['i']
    #     _lyso_ppm_lst = lyso_info_dct['ppm']
    #     _lyso_delta_lst = lyso_info_dct['D']
    #
    #     lyso_info_df = pd.DataFrame(lyso_info_dct, columns=['name', 'mz', 'i', 'ppm'])
    #     lyso_info_df = lyso_info_df.sort_values(by='i', ascending=False)
    #     lyso_info_df = lyso_info_df.round({'mz': 4})
    #
    #     # print table_info_df.head(5)
    #     old_idx_lst = lyso_info_df.index.tolist()
    #     lyso_info_df.index = range(1, len(old_idx_lst) + 1)
    #
    #     col_labels = lyso_info_df.columns.tolist()
    #     # col_labels = mz_info_df.head().tolist()
    #     row_labels = lyso_info_df.index.tolist()
    #     table_vals = map(list, lyso_info_df.values)
    #     # plot lyso table
    #     try:
    #         # the rectangle is where I want to place the table
    #         lyso_table = msms_high_pic.table(cellText=table_vals, rowLabels=row_labels,
    #                                          colWidths=[.22, .1, .1, .1],
    #                                          colLabels=col_labels, loc='upper center')
    #         lyso_table.set_fontsize(7)
    #         _auto_ident_chker += 1
    #     except IndexError:
    #         pass
    #
    #     # get assignment
    #     ident_struct = AssignStructure()
    #     _match_fa_df = pd.DataFrame()
    #     usr_std_fa_df = pd.read_csv(fa_list_csv)
    #     usr_std_fa_df['C'].astype(int)
    #     usr_std_fa_df['DB'].astype(int)
    #     for _ident_idx, _ident_row in usr_std_fa_df.iterrows():
    #         pre_ident_fa_df = usr_std_fa_df
    #         pre_ident_fa_df['abs'] = 0
    #         tmp_fa = pre_ident_fa_df.ix[_ident_idx]['FA']
    #         if tmp_fa in fa_info_dct['fa']:
    #             _tmp_fa_idx = fa_info_dct['fa'].index(tmp_fa)
    #             _tmp_i = fa_info_dct['abs'][_tmp_fa_idx]
    #             pre_ident_fa_df.set_value(_ident_idx, 'abs', _tmp_i)
    #             _match_fa_df = _match_fa_df.append(pre_ident_fa_df.ix[_ident_idx])
    #     _match_lyso_df = pd.DataFrame()
    #     for _ident_idx, _ident_row in usr_std_fa_df.iterrows():
    #         pre_lyso_fa_df = usr_std_fa_df
    #         pre_lyso_fa_df['type'] = ''
    #         pre_lyso_fa_df['abs'] = 0
    #         tmp_lyso = pre_lyso_fa_df.ix[_ident_idx]['FA']
    #         lyso_info_zip_lst = zip(lyso_info_dct['type'], lyso_info_dct['fa'])
    #         if ('Lyso-H2O', tmp_lyso) in lyso_info_zip_lst:
    #             # _tmp_lyso_idx = lyso_info_dct['fa'].index(tmp_lyso)
    #             # _tmp_type = lyso_info_dct['type'][_tmp_lyso_idx]
    #             _tmp_lyso_idx = lyso_info_zip_lst.index(('Lyso-H2O', tmp_lyso))
    #             _tmp_i = lyso_info_dct['abs'][_tmp_lyso_idx]
    #             pre_lyso_fa_df.set_value(_ident_idx, 'type', 'Lyso-H2O')
    #             pre_lyso_fa_df.set_value(_ident_idx, 'abs', _tmp_i)
    #             _match_lyso_df = _match_lyso_df.append(pre_lyso_fa_df.ix[_ident_idx])
    #         if ('Lyso', tmp_lyso) in lyso_info_zip_lst:
    #             _tmp_lyso_idx = lyso_info_zip_lst.index(('Lyso', tmp_lyso))
    #             _tmp_i = lyso_info_dct['abs'][_tmp_lyso_idx]
    #             pre_lyso_fa_df.set_value(_ident_idx, 'type', 'Lyso')
    #             pre_lyso_fa_df.set_value(_ident_idx, 'abs', _tmp_i)
    #             _match_lyso_df = _match_lyso_df.append(pre_lyso_fa_df.ix[_ident_idx])
    #
    #     print ('_match_fa_df', _match_fa_df.shape, '_match_lyso_df', _match_lyso_df.shape)
    #
    #     _ident_df = ident_struct.check(abbr_id, _match_fa_df,
    #                                    _match_lyso_df, usr_std_fa_df)
    #     print ('_ident_df')
    #     print (_ident_df)
    #
    #     if _ident_df.shape[0] > 0:
    #         # print fa_info
    #
    #         _ident_table_df = _ident_df.loc[:, ['Abbr', 'Score']]
    #         _ident_table_df = _ident_table_df.sort_values(by=['Score', 'Abbr'],
    #                                                       ascending=[False, True])
    #         # print table_info_df.head(5)
    #         old_idx_lst = _ident_table_df.index.tolist()
    #         _ident_table_df.index = range(1, len(old_idx_lst) + 1)
    #
    #         col_labels = _ident_table_df.columns.tolist()
    #         # col_labels = mz_info_df.head().tolist()
    #         row_labels = _ident_table_df.index.tolist()
    #         table_vals = map(list, _ident_table_df.values)
    #
    #         try:
    #             # the rectangle is where I want to place the table
    #             fa_table = ms_pic.table(cellText=table_vals, rowLabels=row_labels,
    #                                     colWidths=[.2] * len(col_labels),
    #                                     colLabels=col_labels, loc='upper left')
    #             fa_table.set_fontsize(6)
    #             _auto_ident_chker += 1
    #         except:
    #             pass
