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
        return (fa_df['_add_pl'] + fa_df['_add_sep0'] +
                fa_df['_base_c_str'] + fa_df['_add_sep2'] + fa_df['_base_db_str'] + fa_df['_add_sep3'] +
                fa_df['_add_c_str'] + fa_df['_add_sep2'] + fa_df['_add_db_str'] + fa_df['_add_sep00'])
    def fa_assign_op_sum(self, fa_df):
        return (fa_df['_add_pl'] + fa_df['_add_sep0'] + fa_df['Link'] + fa_df['_add_sep1'] +
                fa_df['_base_c_str'] + fa_df['_add_sep2'] + fa_df['_base_db_str'] + fa_df['_add_sep3'] +
                fa_df['_add_c_str'] + fa_df['_add_sep2'] + fa_df['_add_db_str'] + fa_df['_add_sep00'])

    def check(self, usr_bulk, fa_df, lyso_df, usr_df):

        print ('Calc Scores!!!')

        try:
            fa_df = fa_df.sort_values(by='abs', ascending=False)
            fa_lst = fa_df['FA'].tolist()
        except KeyError:
            fa_lst = []
        try:
            lyso_df = lyso_df.sort_values(by='abs', ascending=False)
            lyso_lst = lyso_df['FA'].tolist()
            lyso_typ_lst = lyso_df['type'].tolist()
            lyso_zip_lst = zip(lyso_lst, lyso_typ_lst)

        except KeyError:
            lyso_zip_lst = [('no_lyso', 'no_type')]

        print('fa_lst', fa_lst)
        print('lyso_zip_lst', lyso_zip_lst)

        pl_checker = re.compile(r'(P[ACEGSI])([(])(.*)([)])')
        pip_checker = re.compile(r'(PIP)([(])(.*)([)])')
        tg_checker = re.compile(r'(TG)([(])(.*)([)])')
        fa_checker = re.compile(r'(\d{1,2})([:])(\d)')
        fa_o_checker = re.compile(r'(O-)(\d{1,2})([:])(\d)')
        fa_p_checker = re.compile(r'(P-)(\d{1,2})([:])(\d)')

        # Check PL Type
        _pl_typ = ''
        _pre_fa_typ_lst = ['']
        if pl_checker.match(usr_bulk):
            print ('PL')
            pl_re_chk = pl_checker.match(usr_bulk)
            pl_typ_lst = pl_re_chk.groups()
            _pl_typ = pl_typ_lst[0]
            _pre_fa_typ_lst = pl_typ_lst[2].split('/')
        if pip_checker.match(usr_bulk):
            print ('PIP')
            pip_re_chk = pip_checker.match(usr_bulk)
            pip_typ_lst = pip_re_chk.groups()
            _pl_typ = pip_typ_lst[0]
            _pre_fa_typ_lst = pip_typ_lst[2].split('/')
        if tg_checker.match(usr_bulk):
            print ('TG')
            pip_re_chk = pip_checker.match(usr_bulk)
            pip_typ_lst = pip_re_chk.groups()
            _pl_typ = pip_typ_lst[0]
            _pre_fa_typ_lst = pip_typ_lst[2].split('/')

        print ('_pl_typ', _pl_typ)
        # check FA composition
        _fa_link = ''
        _fa_sum_num = 0
        _db_sum_num = 0
        for _fa_typ in _pre_fa_typ_lst:
            if fa_checker.match(_fa_typ):
                fa_re_chk = fa_checker.match(_fa_typ)
                fa_typ_lst = fa_re_chk.groups()
                _fa_link += 'A'
                _fa_sum_num += int(fa_typ_lst[0])
                _db_sum_num += int(fa_typ_lst[2])
            elif fa_o_checker.match(_fa_typ):
                fa_o_re_chk = fa_o_checker.match(_fa_typ)
                fa_typ_lst = fa_o_re_chk.groups()
                _fa_link += 'O'
                _fa_sum_num += int(fa_typ_lst[1])
                _db_sum_num += int(fa_typ_lst[3])
            elif fa_p_checker.match(_fa_typ):
                fa_p_re_chk = fa_p_checker.match(_fa_typ)
                fa_typ_lst = fa_p_re_chk.groups()
                _fa_link += 'P'
                _fa_sum_num += int(fa_typ_lst[1])
                _db_sum_num += int(fa_typ_lst[3])

        usr_pl_cmp_lst = [_pl_typ, _fa_link, _fa_sum_num, _db_sum_num]
        # print (usr_pl_cmp_lst)
        print('_fa_link', _fa_link)

        if _fa_link in ['OA', 'PA', 'AA']:
            usr_pl_fa_lst = [_fa_link, str(_fa_sum_num), str(_db_sum_num)]
            usr_pl_fa_code = '-'.join(usr_pl_fa_lst)
        elif _fa_link in ['O', 'P']:
            _fa_link += 'A'
            usr_pl_fa_lst = [_fa_link, str(_fa_sum_num), str(_db_sum_num)]
            usr_pl_fa_code = '-'.join(usr_pl_fa_lst)
        elif _fa_link in ['A', '']:
            _fa_link = 'AA'
            usr_pl_fa_lst = [_fa_link, str(_fa_sum_num), str(_db_sum_num)]
            usr_pl_fa_code = '-'.join(usr_pl_fa_lst)
        else:
            _fa_link = 'AA'
            usr_pl_fa_lst = [_fa_link, str(_fa_sum_num), str(_db_sum_num)]
            usr_pl_fa_code = '-'.join(usr_pl_fa_lst)
        print ('usr_pl_fa_code', usr_pl_fa_code)

        _match_fa_df = pd.DataFrame()

        for idx, row in fa_df.iterrows():
            _add_fa = fa_df.ix[idx]['FA']
            _add_c = int(fa_df.ix[idx]['C'])
            _add_db = int(fa_df.ix[idx]['DB'])
            _add_link = fa_df.ix[idx]['Link']
            _chk_fa_df = usr_df
            _chk_fa_df['_add_fa'] = _add_fa
            _chk_fa_df['_add_pl'] = _pl_typ
            _chk_fa_df['_add_c'] = _add_c
            _chk_fa_df['_add_db'] = _add_db
            _chk_fa_df['_add_link'] = _add_link
            _chk_fa_df['_add_sep1'] = '-'
            _chk_fa_df['_add_sep2'] = ':'
            _chk_fa_df['_add_sep3'] = '/'
            _chk_fa_df['_add_sep0'] = '('
            _chk_fa_df['_add_sep00'] = ')'
            _chk_fa_df['_base_c_str'] = _chk_fa_df['C'].astype(str)
            _chk_fa_df['_base_db_str'] = _chk_fa_df['DB'].astype(str)
            _chk_fa_df['_add_c_str'] = str(_add_c)
            _chk_fa_df['_add_db_str'] = str(_add_db)
            _chk_fa_df['Sum_C'] = _chk_fa_df.apply(self.fa_c_sum, axis=1).astype(str)
            _chk_fa_df['Sum_DB'] = _chk_fa_df.apply(self.fa_db_sum, axis=1).astype(str)
            _chk_fa_df['Sum_Link'] = _chk_fa_df.apply(self.fa_link_sum, axis=1)
            _chk_fa_df['Sum_Code'] = _chk_fa_df.apply(self.fa_code_sum, axis=1)
            _chk_fa_df['Sum_abbr'] = _chk_fa_df.apply(self.fa_abbr_sum, axis=1)

            # print(_chk_fa_df)
            # print (_chk_fa_df.head())
            _tmp_fa_df = _chk_fa_df[_chk_fa_df['Sum_Code'] == usr_pl_fa_code]
            print('_tmp_fa_df', _tmp_fa_df)
            # print(_match_fa_df.head())
            if _tmp_fa_df.shape[0] > 0:
                for _tmp_idx, _tmp_row in _tmp_fa_df.iterrows():
                    if _tmp_fa_df.ix[_tmp_idx]['Link'] in ['O', 'P']:
                        _tmp_fa_df['Abbr'] = _chk_fa_df.apply(self.fa_assign_op_sum, axis=1)
                    else:
                        _tmp_fa_df['Abbr'] = _chk_fa_df.apply(self.fa_assign_sum, axis=1)

                    _tmp_fa_df['Score'] = 0.0
                    _tmp_score_lst = []
                    _tmp_add_fa = _tmp_fa_df.ix[_tmp_idx]['_add_fa']
                    _tmp_base_fa = _tmp_fa_df.ix[_tmp_idx]['FA']
                    if _tmp_add_fa in fa_lst:
                        _tmp_add_fa_idx = fa_lst.index(_tmp_add_fa)
                        if _tmp_add_fa_idx < 10:
                            _tmp_add_fa_factor = (10 - _tmp_add_fa_idx) * 0.1
                            _tmp_score_lst.append(25 * _tmp_add_fa_factor)
                    else:
                        pass
                    if (_tmp_add_fa, 'Lyso-H2O') in lyso_zip_lst:
                        _tmp_add_lysow_idx = lyso_zip_lst.index((_tmp_add_fa, 'Lyso-H2O'))
                        if _tmp_add_lysow_idx < 10:
                            _tmp_add_lysow_factor = (10 - _tmp_add_lysow_idx) * 0.1
                            _tmp_score_lst.append(15 * _tmp_add_lysow_factor)
                    else:
                        pass
                    if (_tmp_add_fa, 'Lyso') in lyso_zip_lst:
                        _tmp_add_lyso_idx = lyso_zip_lst.index((_tmp_add_fa, 'Lyso'))
                        if _tmp_add_lyso_idx < 10:
                            _tmp_add_lyso_factor = (10 - _tmp_add_lyso_idx) * 0.1
                            _tmp_score_lst.append(10 * _tmp_add_lyso_factor)
                    else:
                        pass
                    if _tmp_base_fa in fa_lst:
                        _tmp_base_fa_idx = fa_lst.index(_tmp_base_fa)
                        if _tmp_base_fa_idx < 10:
                            _tmp_base_fa_factor = (10 - _tmp_base_fa_idx) * 0.1
                            _tmp_score_lst.append(25 * _tmp_base_fa_factor)
                    else:
                        pass
                    # if (_tmp_base_fa, 'Lyso-H2O') in lyso_zip_lst:
                    #     _tmp_base_lysow_idx = lyso_zip_lst.index((_tmp_base_fa, 'Lyso-H2O'))
                    #     if _tmp_base_lysow_idx < 10:
                    #         _tmp_base_lysow_factor = (10 - _tmp_base_lysow_idx) * 0.1
                    #         _tmp_score_lst.append(15 * _tmp_base_lysow_factor)
                    # else:
                    #     pass
                    # if (_tmp_base_fa, 'Lyso') in lyso_zip_lst:
                    #     _tmp_base_lyso_idx = lyso_zip_lst.index((_tmp_base_fa, 'Lyso'))
                    #     if _tmp_base_lyso_idx < 10:
                    #         _tmp_base_lyso_factor = (10 - _tmp_base_lyso_idx) * 0.1
                    #         _tmp_score_lst.append(10 * _tmp_base_lyso_factor)
                    # else:
                    #     pass
                    _tmp_score = sum(_tmp_score_lst)
                    print ('_tmp_score', _tmp_score)
                    _tmp_fa_df.set_value(_tmp_idx, 'Score', _tmp_score)
                # print (_tmp_fa_df['Sum_abbr'].tolist())

                _match_fa_df = _match_fa_df.append(_tmp_fa_df)
            else:
                # print ('Not Match!')
                pass
        return _match_fa_df

# usr_csv = r'D:\LPPhunter\FA_list.csv'
# chk = AssignStructure()
# fadf = pd.read_csv(usr_csv)
#
# chk.check('PG(O-36:5)', fadf, fadf)
