# -*- coding: utf-8 -*-
# Copyright 2014-16 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release.
# For more info please contact zhixu.ni@uni-leipzig.de

import re
import pandas as pd


class AssignStructure:

    def fa_c_sum(self, fa_df):
        return fa_df['C'] + fa_df['_add_c']

    def fa_db_sum(self, fa_df):
        return fa_df['DB'] + fa_df['_add_db']

    def fa_link_sum(self, fa_df):
        return fa_df['Link'] + fa_df['_add_link']

    def fa_code_sum(self, fa_df):
        return fa_df['Sum_Link'] + fa_df['_add_sep1'] + fa_df['Sum_C'] + fa_df['_add_sep1'] + fa_df['Sum_DB']

    def fa_abbr_sum(self, fa_df):
        return (fa_df['Sum_Link'] + fa_df['_add_sep1'] +
                fa_df['_base_c_str'] + fa_df['_add_sep2'] + fa_df['_base_db_str'] + fa_df['_add_sep3'] +
                fa_df['_add_c_str'] + fa_df['_add_sep2'] + fa_df['_add_db_str'])

    def fa_assign_sum(self, fa_df):
        return (fa_df['_base_c_str'] + fa_df['_add_sep2'] + fa_df['_base_db_str'] + fa_df['_add_sep3'] +
                fa_df['_add_c_str'] + fa_df['_add_sep2'] + fa_df['_add_db_str'])
    def fa_assign_op_sum(self, fa_df):
        return (fa_df['Link'] + fa_df['_add_sep1'] +
                fa_df['_base_c_str'] + fa_df['_add_sep2'] + fa_df['_base_db_str'] + fa_df['_add_sep3'] +
                fa_df['_add_c_str'] + fa_df['_add_sep2'] + fa_df['_add_db_str'])

    def check(self, usr_bulk, fa_df, lyso_df):

        fa_lst = fa_df.index.tolist()
        lyso_lst = lyso_df.index.tolist()

        pl_checker = re.compile(r'(P[ACEGS])([(])(.*)([)])')
        pip_checker = re.compile(r'(PIP)([(])(.*)([)])')
        fa_checker = re.compile(r'(\d{1,2})([:])(\d)')
        fa_o_checker = re.compile(r'(O-)(\d{1,2})([:])(\d)')
        fa_p_checker = re.compile(r'(P-)(\d{1,2})([:])(\d)')

        # Check PL Type
        _pl_typ = ''
        _fa_typ = ''
        if pl_checker.match(usr_bulk):
            print ('PL')
            pl_re_chk = pl_checker.match(usr_bulk)
            pl_typ_lst = pl_re_chk.groups()
            _pl_typ = pl_typ_lst[0]
            _fa_typ = pl_typ_lst[2]
        if pip_checker.match(usr_bulk):
            print ('PIP')
            pip_re_chk = pip_checker.match(usr_bulk)
            pip_typ_lst = pip_re_chk.groups()
            _pl_typ = pip_typ_lst[0]
            _fa_typ = pip_typ_lst[2]

        # check FA composition
        _fa_link = 'PA'
        _fa_sum_num = 0
        _db_sum_num = 0
        if fa_checker.match(_fa_typ):
            fa_re_chk = fa_checker.match(_fa_typ)
            fa_typ_lst = fa_re_chk.groups()
            _fa_link = 'AA'
            _fa_sum_num = int(fa_typ_lst[0])
            _db_sum_num = int(fa_typ_lst[2])
        elif fa_o_checker.match(_fa_typ):
            fa_o_re_chk = fa_o_checker.match(_fa_typ)
            fa_typ_lst = fa_o_re_chk.groups()
            _fa_link = 'OA'
            _fa_sum_num = int(fa_typ_lst[1])
            _db_sum_num = int(fa_typ_lst[3])
        elif fa_p_checker.match(_fa_typ):
            fa_p_re_chk = fa_p_checker.match(_fa_typ)
            fa_typ_lst = fa_p_re_chk.groups()
            _fa_link = 'PA'
            _fa_sum_num = int(fa_typ_lst[1])
            _db_sum_num = int(fa_typ_lst[3])

        usr_pl_cmp_lst = [_pl_typ, _fa_link, _fa_sum_num, _db_sum_num]
        print (usr_pl_cmp_lst)
        usr_pl_fa_lst = [_fa_link, str(_fa_sum_num), str(_db_sum_num)]
        usr_pl_fa_code = '-'.join(usr_pl_fa_lst)
        print (usr_pl_fa_code)

        _match_fa_df = pd.DataFrame()

        for idx, row in fa_df.iterrows():
            _add_c = fa_df.ix[idx]['C']
            _add_db = fa_df.ix[idx]['DB']
            _add_link = fa_df.ix[idx]['Link']
            _chk_fa_df = fa_df
            _chk_fa_df['_add_c'] = _add_c
            _chk_fa_df['_add_db'] = _add_db
            _chk_fa_df['_add_link'] = _add_link
            _chk_fa_df['_add_sep1'] = '-'
            _chk_fa_df['_add_sep2'] = ':'
            _chk_fa_df['_add_sep3'] = '/'
            _chk_fa_df['_base_c_str'] = _chk_fa_df['C'].astype(str)
            _chk_fa_df['_base_db_str'] = _chk_fa_df['DB'].astype(str)
            _chk_fa_df['_add_c_str'] = str(_add_c)
            _chk_fa_df['_add_db_str'] = str(_add_db)
            _chk_fa_df['Sum_C'] = _chk_fa_df.apply(self.fa_c_sum, axis=1).astype(str)
            _chk_fa_df['Sum_DB'] = _chk_fa_df.apply(self.fa_db_sum, axis=1).astype(str)
            _chk_fa_df['Sum_Link'] = _chk_fa_df.apply(self.fa_link_sum, axis=1)
            _chk_fa_df['Sum_Code'] = _chk_fa_df.apply(self.fa_code_sum, axis=1)
            _chk_fa_df['Sum_abbr'] = _chk_fa_df.apply(self.fa_abbr_sum, axis=1)

            # print (_chk_fa_df.head())
            _tmp_fa_df = _chk_fa_df[_chk_fa_df['Sum_Code'] == usr_pl_fa_code]
            # print(_match_fa_df.head())
            if _tmp_fa_df.shape[0] > 0:
                for _tmp_idx, _tmp_row in _tmp_fa_df.iterrows():
                    if _tmp_fa_df.ix[_tmp_idx]['Link'] in ['O', 'P']:
                        _tmp_fa_df['Abbr'] = _chk_fa_df.apply(self.fa_assign_op_sum, axis=1)
                    else:
                        _tmp_fa_df['Abbr'] = _chk_fa_df.apply(self.fa_assign_sum, axis=1)
                print (_tmp_fa_df['Sum_abbr'].tolist())
                _match_fa_df = _match_fa_df.append(_tmp_fa_df)
            else:
                print ('Not Match!')
        return _match_fa_df

# usr_csv = r'D:\LPPhunter\FA_list.csv'
# chk = AssignStructure()
# fadf = pd.read_csv(usr_csv)
#
# chk.check('PG(O-36:5)', fadf, fadf)
