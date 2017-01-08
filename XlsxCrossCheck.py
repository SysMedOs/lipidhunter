# -*- coding: utf-8 -*-
# Copyright 2014-16 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release.
# For more info please contact zhixu.ni@uni-leipzig.de
from __future__ import print_function
import pandas as pd

in_a_f = r'D:\LPPhunter\images\70min_SIN_II\PI\Marked_PI_70min_mgfHunter_3.xlsx'
in_b_f = r'D:\LPPhunter\images\C16h\Manual_Checked_PI_C16h_mgfHunter_3.xlsx'
out_f = r'D:\LPPhunter\images\PI_CrossComp_C16h_70minSIN_Same.xlsx'
out_c_f = r'D:\LPPhunter\images\PI_CrossComp_C16h_70minSIN_only_c.xlsx'
out_70_f = r'D:\LPPhunter\images\PI_CrossComp_C16h_70minSIN_only_70.xlsx'

c16h_df = pd.read_excel(in_b_f, sheetname='Marked')
sin70min_df = pd.read_excel(in_a_f, sheetname='Marked')

print(c16h_df.head())
print(sin70min_df.head())

same_df = pd.DataFrame()
only_c_df = pd.DataFrame()
only_70_df = pd.DataFrame()
for c_idx, c_row in c16h_df.iterrows():

    _pl_class = c16h_df.get_value(c_idx, 'Class')
    _pl_abbr = c16h_df.get_value(c_idx, 'Abbreviation')
    _pl_assi = c16h_df.get_value(c_idx, 'Assignment')
    _pl_rt = c16h_df.get_value(c_idx, 'rt')

    _query_str = 'Class == "%s" and Abbreviation == "%s" and Assignment == "%s"' % (_pl_class, _pl_abbr, _pl_assi)

    _tmp_df = sin70min_df.query(_query_str)
    if _tmp_df.shape[0] > 0:
        print('Found', _query_str)
        _tmp_df['rt_c'] = _pl_rt
        _tmp_df['Sample2'] = 'C16h_II'
        same_df = same_df.append(_tmp_df)
    else:
        print('Unique', _query_str)
        only_c_df = only_c_df.append(c_row)

print('check for 70 min ------------------------->')

for s_idx, s_row in sin70min_df.iterrows():
    _pl_class = sin70min_df.get_value(s_idx, 'Class')
    _pl_abbr = sin70min_df.get_value(s_idx, 'Abbreviation')
    _pl_assi = sin70min_df.get_value(s_idx, 'Assignment')

    _query_str = 'Class == "%s" and Abbreviation == "%s" and Assignment == "%s"' % (_pl_class, _pl_abbr, _pl_assi)

    _tmp_df = same_df.query(_query_str)
    if _tmp_df.shape[0] > 0:
        print('Found', _query_str)
        pass
    else:
        print('Unique', _query_str)
        only_70_df = only_70_df.append(s_row)

same_df.to_excel(out_f)
only_c_df.to_excel(out_c_f)
only_70_df.to_excel(out_70_f)
print('fin!')
