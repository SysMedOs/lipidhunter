# -*- coding: utf-8 -*-
# Copyright 2014-16 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig          
# The software is currently  under development and is not ready to be released. 
# A suitable license will be chosen before the official release.               
# For more info please contact zhixu.ni@uni-leipzig.de

import xml.etree.cElementTree as ET


spectra_obj = ET.ElementTree(file='orbi.mzML')

spectra_obj.getroot()

