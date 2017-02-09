# -*- coding: utf-8 -*-
# Copyright 2016-2017 LPP team, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of LipidHunter.
# For more info please contact:
#     LPP team oxlpp@bbz.uni-leipzig.de
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de
#     Developer Georgia Angelidou georgia.angelidou@uni-leipzig.de

import os
import shutil
import pandas as pd


class LogPageCreator(object):
    def __init__(self, output_folder, start_time, params):
        print(os.getcwd())
        self.output_folder = output_folder
        self.output_img_folder = output_folder + r'\LipidHunter_Results_Figures_%s' % start_time
        self.main_page = output_folder + r'\LipidHunter_Results_%s.html' % start_time
        self.logo = r'LipidHunter_Results_Figures_%s\LipidHunter.ico' % start_time
        _image_lst_page = r'LipidHunter_Results_Figures_%s\LipidHunter_Results_Figures_list.html' % start_time
        _params_lst_page = r'LipidHunter_Results_Figures_%s\LipidHunter_Params_list.html' % start_time
        _idx_lst_page = r'LipidHunter_Results_Figures_%s\LipidHunter_Identification_list.html' % start_time
        self.image_lst_page = self.output_img_folder + r'\LipidHunter_Results_Figures_list.html'
        self.params_lst_page = self.output_img_folder + r'\LipidHunter_Params_list.html'
        self.idx_lst_page = self.output_img_folder + r'\LipidHunter_Identification_list.html'

        self.lipid_type = params['lipid_type']
        hunter_folder = params['hunter_folder']

        if params['rank_score'] is True:
            score_mode = '(Rank mode)'
        else:
            score_mode = '(Intensity mode)'

        if params['fast_isotope'] is True:
            isotope_score_mode = '(Fast mode)'
        else:
            isotope_score_mode = ''

        with open(self.main_page, 'w') as _m_page:
            m_info_lst = ['<html>\n', '<link rel="icon" href="', self.logo, '" type="image/x-icon"/>\n'
                          '<title>LipidHunter_Results ', start_time,
                          '</title>\n<frameset cols="350,*">\n<frameset rows="390,*">\n',
                          '<frame src="', _params_lst_page, '" frameborder="0" >\n',
                          '<frame src="', _idx_lst_page, '" frameborder="0" >\n</frameset>\n',
                          '<frame src="', _image_lst_page, '"name ="results_frame">\n</frameset>\n</html>\n'
                          ]
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

        with open(self.params_lst_page, 'w') as _params_page:
            _params_page.write('''
                                <html>\n<body>\n<style type="text/css">\n
                                p {margin-left: 20px; text-decoration: none; font-family: sans-serif;}\n
                                body {font-family: sans-serif;}\n table{width:100%s;}\n
                                table, th, td {font-size:14px;text-align: center; font-family: sans-serif;}\n
                                th{background-color:#0066B2;color:white; margin:center;}\n
                                a:link {text-decoration:none} a:hover{text-decoration:underline }\n
                                ul {font-size:14px; width: 260px;}\n </style>\n
                                <h3><img src="LipidHunter.ico" height=30/>  LipidHunter</h3><font size="1">\n
                                <hr> <h4>Parameters:</h4>\n<ul>\n
                                <li>Start time: %s</li>\n<li>Mode: %s %s</li>\n
                                <li><i>m/z</i> range: %.1f - %.1f <i>m/z</i></li>\n<li>RT range: %.1f - %.1f min</li>\n
                                <li>MS1 Threshold: %i</li>\n<li>MS2 Threshold: %i</li>\n
                                <li>MS1 ppm: %i</li>\n<li>MS2 ppm: %i</li>\n
                                <li>LipidHunter score > %.1f %s</li>\n<li>Isotope score > %.1f %s</li>\n
                                </ul>\n<hr>\n<h4>Lipid identification list:</h4><font size="1">\n<table>\n<thead>\n
                                <tr style="text-align: center;">\n
                                <th>id#</th>\n<th> Exact Mass </th>\n<th>RT(min)</th>\n<th>Bulk</th>\n<th>Score</th>\n
                                </tr>\n</thead>\n</table>\n</body>\n</html>\n
                                ''' % ('%', params['hunter_start_time'], self.lipid_type, params['charge_mode'],
                                       params['mz_start'], params['mz_end'], params['rt_start'], params['rt_end'],
                                       params['ms_th'], params['ms2_th'], params['ms_ppm'], params['ms2_ppm'],
                                       params['score_filter'], score_mode,
                                       params['isotope_score_filter'], isotope_score_mode))

        with open(self.idx_lst_page, 'w') as _idx_page:
            _idx_page.write('''
                            <html>\n<body>\n<style type="text/css">\n
                            p {margin-left: 20px; text-decoration: none; font-family: sans-serif;}\n
                            body {background-color: #B1D3EC;font-family: sans-serif;}\n table{width:100%s;}\n
                            table, th, td {font-size:14px;text-align: center; font-family: sans-serif;}\n
                            th{background-color:#0066B2;color:white;margin:center;}\n
                            tr:nth-child(even){background-color: #7C94A5;}\n
                            a:link {text-decoration:none; color:black}a:hover{text-decoration:underline; color:black;}
                            a:visited {text-decoration:none; color:black;}\n
                            </style>\n<table>\n<tbody>
                            ''' % '%')
        try:
            shutil.copy('%s\LipidHunter.ico' % hunter_folder, self.output_img_folder)
        except IOError:
            pass

    def add_info(self, img_name, ident_idx, ident_info_df):

        print('try to add identification to report html')

        img_path = img_name[1:]
        ident_idx = str(ident_idx)

        ms2_pr_mz = ident_info_df.get_value(1, r'MS2_PR_mz')
        ms2_rt = ident_info_df.get_value(1, 'MS2_scan_time')
        dda = ident_info_df.get_value(1, 'DDA#')
        ms2_scan_id = ident_info_df.get_value(1, 'Scan#')
        abbr_bulk = ident_info_df.get_value(1, 'Bulk_identification')
        ident_abbr = ident_info_df.get_value(1, 'Proposed_structures')
        score = ident_info_df.get_value(1, 'Score')
        formula_ion = ident_info_df.get_value(1, 'Formula_ion')
        charge = ident_info_df.get_value(1, 'Charge')

        with open(self.image_lst_page, 'a') as img_page:
            # convert info df to html table code
            if self.lipid_type in ['PA', 'PC', 'PE', 'PG', 'PI', 'PIP', 'PS']:
                plot_df_cols = ['Proposed_structures', 'Score', 'i_sn1', 'i_sn2', 'i_M-sn1', 'i_M-sn2',
                                'i_M-(sn1-H2O)', 'i_M-(sn2-H2O)']
            elif self.lipid_type in ['TG', 'TAG', 'DG', 'DAG', 'MG', 'MAG']:
                plot_df_cols = ['Proposed_structures', 'Score', 'i_sn1', 'i_sn2', 'i_sn3',
                                'i_M-sn1', 'i_M-sn2', 'i_M-sn3',
                                'i_M-(sn1+sn2)', 'i_M-(sn1+sn3)', 'i_M-(sn2+sn3)']
            else:
                plot_df_cols = ['Proposed_structures', 'Score', 'i_sn1', 'i_sn2', 'i_M-sn1', 'i_M-sn2',
                                'i_M-(sn1-H2O)', 'i_M-(sn2-H2O)']
            ident_col = ident_info_df.columns.tolist()

            for _col in plot_df_cols:
                if _col not in ident_col:
                    ident_info_df.loc[:, _col] = ''
            table_buf_code = ident_info_df.to_html(columns=plot_df_cols, float_format='%.1f', border=0)
            table_buf_code = table_buf_code.replace('NaN', '')
            img_title_str = ('{mz}_RT{rt:.3}_DDAtop{dda}_Scan{scan}_{bulk}_{ident}_{f}_{chg}_score{score}'
                             .format(mz=ms2_pr_mz, rt=ms2_rt, dda=dda, scan=ms2_scan_id, bulk=abbr_bulk,
                                     ident=ident_abbr, score=score, f=formula_ion, chg=charge))
            img_info_lst = ['<a name="', ident_idx, '"><h3>', '<a href="', img_path, '" target="blank">', img_title_str,
                            '</a></h3></a>', '<a href="', img_path, '" target="blank">',
                            '<img src="', img_path, '" height="720" /></a>', table_buf_code, '\n<hr>\n']
            img_page.write(''.join(img_info_lst))

        with open(self.idx_lst_page, 'a') as idx_page:

            idx_str = ('''
                        <tr>\n
                        <td><a href ="LipidHunter_Results_Figures_list.html#{id}" target ="results_frame">{id}</td>\n
                        <td><a href ="LipidHunter_Results_Figures_list.html#{id}" target ="results_frame">{mz}</td>\n
                        <td><a href ="LipidHunter_Results_Figures_list.html#{id}" target ="results_frame">{rt}</td>\n
                        <td><a href ="LipidHunter_Results_Figures_list.html#{id}" target ="results_frame">{bulk}</td>\n
                        <td><a href ="LipidHunter_Results_Figures_list.html#{id}" target ="results_frame">{score}</td>\n
                        </tr>\n
                        '''.format(id=ident_idx, mz='%.4f' % ms2_pr_mz, rt='%.1f' % ms2_rt,
                                   bulk=abbr_bulk, score=score))
            idx_page.write(idx_str)

        print('==> info added to report html -->')

    def close_page(self):
        with open(self.main_page, 'a') as _m_page:
            _m_page.write('\n</body></html>\n')

        with open(self.image_lst_page, 'a') as _img_page:
            _img_page.write('\n</body></html>\n')

        with open(self.idx_lst_page, 'a') as _idx_page:
            _idx_page.write('\n</tbody>\n</table>\n</body></html>\n')
