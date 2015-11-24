# -*- coding: utf-8 -*-
# Copyright 2014-16 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig          
# The software is currently  under development and is not ready to be released. 
# A suitable license will be chosen before the official release.
# For more info please contact zhixu.ni@uni-leipzig.de

from base64 import b64decode as b64dec
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

import zlib
import pymzml
import pymzml.spec
import pymzml.run


class XIC(object):

    def __init__(self, mzml):

        self.mzml_obj = pymzml.run.Reader(mzml, MS1_Precision=20e-6, MSn_Precision=200e-6)


    def find_mz(self, mz):

        # timeDependentIntensities = []
        '''
        ['total ion current', 'filter string', 'PY:0000000', 'id', 'MS:1000285', 'MS:1000512', 'MS:1000511',
        'MS:1000515', 'MS:1000514', 'scan start time', 'intensity array', 'defaultArrayLength', 'm/z array',
        'MS:1000128', 'ms level', 'MS:1000574', 'BinaryArrayOrder', 'profile spectrum', '32-bit float',
        'zlib compression', 'MS:1000016', 'encodedData', 'MS:1000521']

        'zlib compression': ['zlib compression', 'zlib compression']
        'MS:1000574': ['zlib compression', 'zlib compression']


        ['MS:1000511', 'MS:1000576', 'MS:1000016', 'encodedData', 'MS:1000515', None, 'MS:1000523',
        'BinaryArrayOrder', 'MS:1000514', 'defaultArrayLength', 'MS:1000127', 'PY:0000000', 'id', 'MS:1000285']

        :param mz:
        :return:
        '''
        '''
        :param mz:
        :return:
        '''
        for spectrum in self.mzml_obj:
            # print spectrum.keys()
            if spectrum['MS:1000511'] == 1:  # 'MS:1000511' == 'ms level'

                # print spectrum['encodedData']
                # Add 32-bit/64-bit encoding to mzML
                try:
                    print spectrum['BinaryArrayOrder']
                # [('encoding', '32-bit float'), ('compression', 'zlib'), ('arrayType', 'mz'),
                # ('encoding', '32-bit float'), ('compression', 'zlib'), ('arrayType', 'i')]
                    if spectrum['BinaryArrayOrder'] == []:
                        spectrum['BinaryArrayOrder'] = [('encoding', '32-bit float'),
                                                        ('compression', 'no'), ('arrayType', 'mz'),
                                                        ('encoding', '32-bit float'), ('compression', 'no'),
                                                        ('arrayType', 'i')]
                        print r"spectrum['BinaryArrayOrder'] existed"
                        print spectrum['BinaryArrayOrder']
                    else:
                        print r"spectrum['BinaryArrayOrder'] existed"
                except KeyError:
                    spectrum['BinaryArrayOrder'] = [('encoding', '32-bit float'),
                                                    ('compression', 'no'), ('arrayType', 'mz'),
                                                    ('encoding', '32-bit float'), ('compression', 'no'),
                                                    ('arrayType', 'i')]
                    print 'BinaryArrayOrder Added!'

                print spectrum.keys()
                for _m, _i in spectrum.highestPeaks(300):
                    print _m, _i


                break

