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

from operator import itemgetter
import os
import shutil

import pandas as pd


class LogPageCreator(object):
    def __init__(self, output_folder, start_time, params):
        print(os.getcwd())
        self.output_folder = output_folder
        self.output_img_folder = output_folder + r'/LipidHunter_Results_Figures_%s' % start_time
        self.main_page = output_folder + r'/LipidHunter_Results_%s.html' % start_time
        self.logo = r'LipidHunter_Results_Figures_%s/LipidHunter.ico' % start_time
        _image_lst_page = r'LipidHunter_Results_Figures_%s/LipidHunter_Results_Figures_list.html' % start_time
        _params_lst_page = r'LipidHunter_Results_Figures_%s/LipidHunter_Params_list.html' % start_time
        _idx_lst_page = r'LipidHunter_Results_Figures_%s/LipidHunter_Identification_list.html' % start_time
        self.image_lst_page = self.output_img_folder + r'/LipidHunter_Results_Figures_list.html'
        self.params_lst_page = self.output_img_folder + r'/LipidHunter_Params_list.html'
        self.idx_lst_page = self.output_img_folder + r'/LipidHunter_Identification_list.html'
        self.cfg_sum_page = self.output_img_folder + r'/LipidHunter_Configuration_Summary.html'

        self.lipid_class = params['lipid_class']
        hunter_folder = params['hunter_folder']

        if params['rank_score'] is True:
            score_mode = ''
        else:
            score_mode = '(Intensity mode)'

        if params['fast_isotope'] is True:
            isotope_score_mode = '(Fast mode)'
        else:
            isotope_score_mode = ''

        with open(self.main_page, 'w') as _m_page:
            # Merge sub pages into one report page, set _params_lst_page width 345, height 345
            m_info_lst = ['<html>\n', '<link rel="icon" href="', self.logo,
                          '" type="image/x-icon"/>\n<title>LipidHunter_Results ', start_time,
                          '</title>\n<frameset cols="345,*">\n<frameset rows="345,*">\n',
                          '<frame src="', _params_lst_page, '" frameborder="0" >\n',
                          '<frame src="', _idx_lst_page, '" frameborder="0" >\n</frameset>\n',
                          '<frame src="', _image_lst_page, '"name ="results_frame">\n</frameset>\n</html>\n']
            _m_page.write(''.join(m_info_lst))

        with open(self.image_lst_page, 'w') as _img_page:
            _img_page.write('''
                            <html>\n<body>\n<style type="text/css">\n
                            p {margin-left: 16px; text-decoration: none; font-family: sans-serif;}\n
                            h3 {font-size:20px; margin-left: 16px; text-decoration: none; font-family: sans-serif;}\n
                            body {font-family: sans-serif;}\n
                            table, th, td {font-size:14px;text-align: center; font-family: sans-serif;}\n
                            th{background-color:#0066B2;color:white; margin:center;}\n
                            tr:nth-child(odd){background-color: #B1D3EC;}\n
                            tr:nth-child(even){background-color: #7C94A5;}\n
                            a:link {text-decoration:none; color:black} a:hover{text-decoration:underline; color:black;}
                            a:visited {text-decoration:none; color:black;}\n</style>\n''')

        with open(self.cfg_sum_page, 'w') as _cfg_page:

            cfg_sum_dct = {'fawhitelist_path_str': 'FA white list:',
                           'mzml_path_str': 'mzML file used:',
                           'img_output_folder_str': 'Output images saved in folder:',
                           'xlsx_output_path_str': 'Summary Output table saved as:',
                           'lipid_specific_cfg': 'PL specific ions defined in:',
                           'score_cfg': 'Weight factor defined in:',
                           'hunter_start_time': 'Run start from:',
                           'vendor': 'Instrument vendor:', 'experiment_mode': 'Experiment type:',
                           'lipid_class': 'Lipid class:', 'charge_mode': 'Precursor charge mode:',
                           'rt_start': 'RT start:', 'rt_end': 'RT end:',
                           'mz_start': r'm/z start:', 'mz_end': r'm/z end:', 'score_filter': 'Score >=',
                           'rank_score': 'Rank score mode:', 'rank_score_filter': 'Rank score >=',
                           'isotope_score_filter': 'Isotope score >=', 'fast_isotope': 'Fast Isotope score mode:',
                           'ms_th': 'Threshold @ MS level =', 'ms_ppm': 'ppm @ MS level <=',
                           'ms_max': 'For precursor intensity <=', 'pr_window': 'Precursor selection window +/-',
                           'dda_top': 'DDA Top =', 'ms2_th': 'Threshold @ MS/MS level =',
                           'ms2_ppm': 'ppm @ MS/MS level <=', 'ms2_infopeak_threshold': 'Consider MS/MS with i (%) >=',
                           'hg_th': 'For PL specific peaks Consider MS/MS with i (abs) >=',
                           'hg_ppm': 'For PL specific peaks Consider MS/MS ppm <=',
                           'ms2_hginfopeak_threshold': 'For PL specific peaks Consider MS/MS with i (%) >=',
                           'hunter_folder': 'LipidHunter program folder',
                           'core_number': 'Run with max CPU core number =',
                           'max_ram': 'Run with max RAM (GB) =',
                           'img_type': 'Save image format:', 'img_dpi': 'Save image with dpi =',
                           'tag_all_sn': 'Prefer all sn identified for TAGs:'}

            cfg_sum_lst = ['lipid_class', 'charge_mode', 'vendor', 'experiment_mode', 'hunter_start_time',
                           'img_output_folder_str', 'xlsx_output_path_str',
                           'mzml_path_str', 'fawhitelist_path_str', 'lipid_specific_cfg', 'score_cfg',
                           'rt_start', 'rt_end', 'mz_start', 'mz_end', 'dda_top',
                           'hunter_folder', 'isotope_score_filter', 'fast_isotope',
                           'score_filter', 'rank_score', 'rank_score_filter',
                           'ms_ppm', 'ms_th', 'pr_window',
                           'ms2_ppm', 'ms2_th', 'ms2_infopeak_threshold',
                           'hg_ppm', 'hg_th', 'ms2_hginfopeak_threshold',
                           'core_number', 'max_ram', 'img_type', 'img_dpi', 'ms_max', 'tag_all_sn']

            param_key_lst = list(params.keys())
            disp_cfg_lst = []

            for _key in cfg_sum_lst:
                if _key in param_key_lst:
                    disp_cfg_lst.append('<li><strong>{k}  </strong>{v}</li>\n'.
                                        format(k=cfg_sum_dct[_key], v=params[_key]))

            params_li_str = ''.join(disp_cfg_lst)

            cfg_template = '''
                           <html>
                           <style type="text/css">
                           p {margin-left: 20px; text-decoration: none; font-family: sans-serif;}
                           body {font-family: sans-serif;}\
                           table{width:100%;}
                           table, th, td {font-size:14px;text-align: center; font-family: sans-serif;}
                           th{background-color:#0066B2;color:white; margin:center;}
                           a:link {text-decoration:none} a:hover{text-decoration:underline }
                           ul
                           </style>
                           <body>
                           <h3><img src="LipidHunter.ico" height=30/>  LipidHunter</h3>
                           <hr> <h3>Parameters:</h3>\n<ul>\n
                           '''
            cfg_template += params_li_str
            cfg_template += '</ul>\n</body>\n</html>\n'
            _cfg_page.write(cfg_template)

        with open(self.params_lst_page, 'w') as _params_page:
            _params_page.write('''
                                <html>
                                <style type="text/css">
                                p {margin-left: 20px; text-decoration: none; font-family: sans-serif;}
                                body {font-family: sans-serif;}
                                table{width:100%s;}
                                table, th, td {font-size:14px;text-align: center; font-family: sans-serif;}
                                th{background-color:#0066B2;color:white; margin:center;}
                                a:link {text-decoration:none} a:hover{text-decoration:underline }
                                ul {font-size:14px; width: 260px;}
                                </style>
                                <body>
                                <h3><img src="LipidHunter.ico" height=30/>LipidHunter</h3><font size="1">
                                <hr> 
                                <h3>Parameters:</h3>
                                <ul>
                                <li>File Name: %s</li>
                                <li>Start time: %s</li>
                                <li>Mode: %s %s</li>
                                <li><i>m/z</i> range: %.1f - %.1f <i>m/z</i></li>
                                <li>RT range: %.1f - %.1f min</li>
                                <li>MS1 Threshold: %i</li>
                                <li>MS2 Threshold: %i</li>
                                <li>MS1 ppm: %i</li>
                                <li>MS2 ppm: %i</li>
                                <li>LipidHunter score > %.1f %s</li>
                                <li>Isotope score > %.1f %s</li>
                                </ul>
                                <h4><a href ="LipidHunter_Configuration_Summary.html" target ="_blank">
                                View all parameters...</a><h4><hr>
                                </body>
                                </html>
                                ''' % ('%', os.path.basename(params['mzml_path_str']), params['hunter_start_time'],
                                       self.lipid_class, params['charge_mode'],
                                       params['mz_start'], params['mz_end'], params['rt_start'], params['rt_end'],
                                       params['ms_th'], params['ms2_th'], params['ms_ppm'], params['ms2_ppm'],
                                       params['rank_score_filter'], score_mode,
                                       params['isotope_score_filter'], isotope_score_mode))
        with open(self.idx_lst_page, 'w') as _idx_page:
            _idx_page.write('''
                            <html>
                            <style type="text/css">
                            p {margin-left: 20px; text-decoration: none; font-family: sans-serif;}
                            body {background-color: #B1D3EC;font-family: sans-serif;}
                            table{width:100%s;}
                            table, th, td {font-size:14px;text-align: center; font-family: sans-serif;}
                            th{background-color:#0066B2;color:white;margin:center;}
                            tr:nth-child(even){background-color: #7C94A5;}
                            a:link {text-decoration:none; color:black}a:hover{text-decoration:underline; color:black;}
                            a:visited {text-decoration:none; color:black;}
                            </style>
                            <body>
                            <h3>Lipid identification list:</h3><font size="1">
                            <table>
                                <thead>
                                <tr style="text-align: center;">
                                <th>ID#</th>
                                <th> MS1_obs_mz </th>
                                <th>RT(min)</th>
                                <th>Discrete</th>
                                <th>Score</th>
                                </tr>
                                </thead>
                            <tbody>
                            ''' % '%')
        try:
            shutil.copy('%s\LipidHunter.ico' % hunter_folder, self.output_img_folder)
        except IOError:
            pass

    def add_all_info(self, ident_info_df):

        with open(self.image_lst_page, 'a') as img_page:
            with open(self.idx_lst_page, 'a') as idx_page:

                _log_info_df = ident_info_df
                _log_info_df.is_copy = False
                _log_info_df['MS1_log_mz'] = _log_info_df['MS1_obs_mz'].round(1)
                _log_info_df = _log_info_df.sort_values(by=['MS1_log_mz', 'Proposed_structures', 'MS2_scan_time',
                                                            'RANK_SCORE'], ascending=[True, True, True, False])
                _log_info_df.reset_index(drop=True, inplace=True)
                _log_info_df.index += 1
                _log_info_groups = _log_info_df.groupby(['MS1_log_mz', 'Proposed_structures', 'Charge',
                                                         'MS2_scan_time'])
                _log_info_groups_key_lst = list(_log_info_groups.groups.keys())
                _log_info_groups_key_lst = sorted(_log_info_groups_key_lst, key=itemgetter(0, 1, 3))
                # _log_info_groups_key_lst = sorted(_log_info_groups_key_lst, key=lambda x: x[0])

                for _idx in range(len(_log_info_groups_key_lst)):
                    _subgroup_df = _log_info_groups.get_group(_log_info_groups_key_lst[_idx])

                    img_path = str(_subgroup_df['img_name'].values.tolist()[0])
                    ms1_pr_mz = _subgroup_df['MS1_obs_mz'].values.tolist()[0]
                    ms2_rt = _subgroup_df['MS2_scan_time'].values.tolist()[0]
                    dda = _subgroup_df['DDA#'].values.tolist()[0]
                    ms2_scan_id = _subgroup_df['Scan#'].values.tolist()[0]
                    ident_abbr = str(_subgroup_df['Proposed_structures'].values.tolist()[0])
                    try:
                        ident_abbr = ident_abbr.replace('<', '&lt;')
                        ident_abbr = ident_abbr.replace('>', '&gt;')
                    except AttributeError:
                        pass
                    score = _subgroup_df['RANK_SCORE'].values.tolist()[0]
                    formula_ion = _subgroup_df['Formula_ion'].values.tolist()[0]
                    charge = _subgroup_df['Charge'].values.tolist()[0]

                    # convert info df to html table code
                    plot_df_cols = []
                    if self.lipid_class in ['PA', 'PC', 'PE', 'PG', 'PI', 'PIP', 'PS']:
                        plot_df_cols = ['Proposed_structures', 'DISCRETE_ABBR', 'RANK_SCORE',
                                        'FA1_[FA-H]-_i_per', 'FA2_[FA-H]-_i_per',
                                        '[LPL(FA1)-H]-_i_per', '[LPL(FA2)-H]-_i_per',
                                        '[LPL(FA1)-H2O-H]-_i_per', '[LPL(FA2)-H2O-H]-_i_per']
                        peak_info_df = pd.DataFrame(_subgroup_df, columns=plot_df_cols)
                        peak_info_df.rename(
                            columns={'FA1_[FA-H]-_i_per': 'FA1_[FA-H]-_i (%)', 'FA2_[FA-H]-_i_per': 'FA2_[FA-H]-_i (%)',
                                     '[LPL(FA1)-H]-_i_per': '[LPL(FA1)-H]-_i (%)',
                                     '[LPL(FA2)-H]-_i_per': '[LPL(FA2)-H]-_i (%)',
                                     '[LPL(FA1)-H2O-H]-_i_per': '[LPL(FA1)-H2O-H]-_i (%)',
                                     '[LPL(FA2)-H2O-H]-_i_per': '[LPL(FA2)-H2O-H]-_i (%)'}, inplace=True)
                    elif self.lipid_class in ['TG', 'TAG', 'MG', 'MAG'] and charge in ['[M+H]+', '[M+NH4]+']:
                        plot_df_cols = ['Proposed_structures', 'DISCRETE_ABBR', 'RANK_SCORE', 'FA1_[FA-H2O+H]+_i_per',
                                        'FA2_[FA-H2O+H]+_i_per', 'FA3_[FA-H2O+H]+_i_per', '[M-(FA1)+H]+_i_per',
                                        '[M-(FA2)+H]+_i_per', '[M-(FA3)+H]+_i_per', '[MG(FA1)-H2O+H]+_i_per',
                                        '[MG(FA2)-H2O+H]+_i_per', '[MG(FA3)-H2O+H]+_i_per']
                        peak_info_df = pd.DataFrame(_subgroup_df, columns=plot_df_cols)
                        peak_info_df.rename(
                            columns={'FA1_[FA-H2O+H]+_i_per': 'FA1_[FA-H2O+H]+_i (%)',
                                     'FA2_[FA-H2O+H]+_i_per': 'FA2_[FA-H2O+H]+_i (%)',
                                     'FA3_[FA-H2O+H]+_i_per': 'FA3_[FA-H2O+H]+_i (%)',
                                     '[M-(FA1)+H]+_i_per': '[M-(FA1)+H]+_i (%)',
                                     '[M-(FA2)+H]+_i_per': '[M-(FA2)+H]+_i (%)',
                                     '[M-(FA3)+H]+_i_per': '[M-(FA3)+H]+_i (%)',
                                     '[MG(FA1)-H2O+H]+_i_per': '[MG(FA1)-H2O+H]+_i (%)',
                                     '[MG(FA2)-H2O+H]+_i_per': '[MG(FA2)-H2O+H]+_i (%)',
                                     '[MG(FA3)-H2O+H]+_i_per': '[MG(FA3)-H2O+H]+_i (%)'}, inplace=True)
                    elif self.lipid_class in ['TG'] and charge in ['[M+Na]+']:
                        plot_df_cols = ['Proposed_structures', 'DISCRETE_ABBR', 'RANK_SCORE', 'FA1_[FA-H2O+H]+_i_per',
                                        'FA2_[FA-H2O+H]+_i_per', 'FA3_[FA-H2O+H]+_i_per', '[M-(FA1)+Na]+_i_per',
                                        '[M-(FA2)+Na]+_i_per', '[M-(FA3)+Na]+_i_per', '[MG(FA1)-H2O+H]+_i_per',
                                        '[MG(FA2)-H2O+H]+_i_per', '[MG(FA3)-H2O+H]+_i_per']
                        peak_info_df = pd.DataFrame(_subgroup_df, columns=plot_df_cols)
                        peak_info_df.rename(
                            columns={'FA1_[FA-H2O+H]+_i_per': 'FA1_[FA-H2O+H]+_i (%)',
                                     'FA2_[FA-H2O+H]+_i_per': 'FA2_[FA-H2O+H]+_i (%)',
                                     'FA3_[FA-H2O+H]+_i_per': 'FA3_[FA-H2O+H]+_i (%)',
                                     '[M-(FA1)+Na]+_i_per': '[M-(FA1)+Na]+_i (%)',
                                     '[M-(FA2)+Na]+_i_per': '[M-(FA2)+Na]+_i (%)',
                                     '[M-(FA3)+Na]+_i_per': '[M-(FA3)+Na]+_i (%)',
                                     '[MG(FA1)-H2O+H]+_i_per': '[MG(FA1)-H2O+H]+_i (%)',
                                     '[MG(FA2)-H2O+H]+_i_per': '[MG(FA2)-H2O+H]+_i (%)',
                                     '[MG(FA3)-H2O+H]+_i_per': '[MG(FA3)-H2O+H]+_i (%)'}, inplace=True)
                    elif self.lipid_class in ['DG'] and charge in ['[M+H]+', '[M+NH4]+']:
                        plot_df_cols = ['Proposed_structures', 'DISCRETE_ABBR', 'RANK_SCORE', 'FA1_[FA-H2O+H]+_i_per',
                                        'FA2_[FA-H2O+H]+_i_per', '[MG(FA1)-H2O+H]+_i_per', '[MG(FA2)-H2O+H]+_i_per']
                        peak_info_df = pd.DataFrame(_subgroup_df, columns=plot_df_cols)
                        peak_info_df.rename(
                            columns={'FA1_[FA-H2O+H]+_i_per': 'FA1_[FA-H2O+H]+_i (%)',
                                     'FA2_[FA-H2O+H]+_i_per': 'FA2_[FA-H2O+H]+_i (%)',
                                     '[M-(FA1)+H]+_i_per': '[M-(FA1)+H]+_i (%)',
                                     '[M-(FA2)+H]+_i_per': '[M-(FA2)+H]+_i (%)',
                                     '[MG(FA1)-H2O+H]+_i_per': '[MG(FA1)-H2O+H]+_i (%)',
                                     '[MG(FA2)-H2O+H]+_i_per': '[MG(FA2)-H2O+H]+_i (%)'}, inplace=True)
                    else:
                        plot_df_cols = ['Proposed_structures', 'DISCRETE_ABBR', 'RANK_SCORE',
                                        'FA1_[FA-H]-i_per', 'FA2_[FA-H]-i_per',
                                        '[LPL(FA1)-H]-i_per', '[LPL(FA2)-H]-i_per',
                                        '[LPL(FA1)-H2O-H]-i_per', '[LPL(FA2)-H2O-H]-i_per']
                        peak_info_df = pd.DataFrame(_subgroup_df, columns=plot_df_cols)

                    try:
                        table_buf_code = peak_info_df.to_html(float_format='%.1f', border=0, index=False)
                    except TypeError:
                        table_buf_code = peak_info_df.to_html(index=False)
                    table_buf_code = table_buf_code.replace('NaN', '')

                    _idx += 1  # set start from 0 to 1

                    img_title_str = ('{mz}_RT{rt:.3}_DDArank{dda}_Scan{scan}_{ident}_{f}_{chg}_score{score}'
                                     .format(mz='%.4f' % ms1_pr_mz, rt=ms2_rt, dda=dda, scan=ms2_scan_id,
                                             ident=ident_abbr, score=score, f=formula_ion, chg=charge))
                    img_info_lst = ['<a name="', '%i' % _idx, '"><h3>', '<a href="', img_path, '" target="blank">',
                                    img_title_str, '</a></h3></a>', '<a href="', img_path, '" target="blank">',
                                    '<img src="', img_path, '" height="800" /></a>', table_buf_code, '\n<hr>\n']
                    img_page.write(''.join(img_info_lst))

                    idx_str = ('''
                            <tr>\n<td>
                            <a href ="LipidHunter_Results_Figures_list.html#{id}" target ="results_frame">{id}
                            </td>\n<td>
                            <a href ="LipidHunter_Results_Figures_list.html#{id}" target ="results_frame">{mz}
                            </td>\n<td>
                            <a href ="LipidHunter_Results_Figures_list.html#{id}" target ="results_frame">{rt}
                            </td>\n<td>
                            <a href ="LipidHunter_Results_Figures_list.html#{id}" target ="results_frame">{ident}
                            </td>\n<td>
                            <a href ="LipidHunter_Results_Figures_list.html#{id}" target ="results_frame">{score}
                            </td>\n</tr>\n
                            '''.format(id='%i' % _idx, mz='%.4f' % ms1_pr_mz, rt='%.1f' % ms2_rt,
                                       ident=ident_abbr, score=score))
                    idx_page.write(idx_str)

            print('==> info added to report html -->')

    def close_page(self):
        with open(self.main_page, 'a') as _m_page:
            _m_page.write('\n</body></html>\n')

        with open(self.image_lst_page, 'a') as _img_page:
            _img_page.write('\n</body></html>\n')

        with open(self.idx_lst_page, 'a') as _idx_page:
            _idx_page.write('\n</tbody>\n</table>\n</body></html>\n')
