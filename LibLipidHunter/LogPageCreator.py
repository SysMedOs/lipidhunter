# -*- coding: utf-8 -*-
# Copyright 2015-2017 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de


class LogPageCreator(object):

    def __init__(self, output_folder):
        self.output_folder = output_folder
        self.main_page = output_folder + r'\LipidHunter_Results.html'
        _image_lst_page = r'lipid_images\lipid_images_list.html'
        _idx_lst_page = r'lipid_images\lipid_index_list.html'
        self.image_lst_page = self.output_folder + r'\lipid_images\lipid_images_list.html'
        self.idx_lst_page = self.output_folder + r'\lipid_images\lipid_index_list.html'

        with open(self.main_page, 'w') as _m_page:
            m_info_lst = ['<html><frameset cols="350,*"><frame src="', _idx_lst_page, '" ><frame src="',
                          _image_lst_page, '"name ="showframe"></frameset></html>\n'
                          ]
            _m_page.write(''.join(m_info_lst))

        with open(self.image_lst_page, 'w') as _img_page:
            _img_page.write('<html><body>')

        with open(self.idx_lst_page, 'w') as _idx_page:
            _idx_page.write('<html><body><h4>List of output images</h4><font size="1">')

    def add_info(self, img_name, ident_idx):
        img_info = img_name[1:-4]
        img_path = img_name[1:]
        ident_idx = str(ident_idx)
        with open(self.image_lst_page, 'a') as img_page:
            img_info_lst = ['<a name="', ident_idx, '"><h2>', '<a href="', img_path, '" target="blank">', img_info,
                            '</a></h2></a>', '<a href="', img_path, '" target="blank">',
                            '<img src="', img_path, '" height="720" /></a>\n']
            img_page.write(''.join(img_info_lst))

        with open(self.idx_lst_page, 'a') as idx_page:
            idx_info_lst = ['<a href ="lipid_images_list.html#', ident_idx, '" target ="showframe">', img_info, '</a>\n']
            idx_page.write(''.join(idx_info_lst))

    def close_page(self):
        with open(self.main_page, 'a') as _m_page:
            _m_page.write('\n</body></html>\n')

        with open(self.image_lst_page, 'a') as _img_page:
            _img_page.write('\n</body></html>\n')

        with open(self.idx_lst_page, 'a') as _idx_page:
            _idx_page.write('\n</body></html>\n')
