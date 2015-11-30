# -*- coding: utf-8 -*-
# Copyright 2014-16 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig          
# The software is currently  under development and is not ready to be released. 
# A suitable license will be chosen before the official release.               
# For more info please contact zhixu.ni@uni-leipzig.de

import pymzml


def check_encode(mzml):

    # The waters.raw converted mzML by proteowizard may not contain encode info.
    # The '32-bit float' or '64-bit float' will be add to allow following process
    # :param spectrum: the spectrum form pymzML
    # :return: the spectrum with encode info inside

    mzml_obj = pymzml.run.Reader(mzml, MS1_Precision=20e-6, MSn_Precision=200e-6)
    print 'try 32-bit/64-bit float'

    # take 1 scan from mzML
    spectrum = mzml_obj[1]

    # check if 32-bit/64-bit encoding to mzML
    # [('encoding', '32-bit float'), ('compression', 'zlib'), ('arrayType', 'mz'),
    # ('encoding', '32-bit float'), ('compression', 'zlib'), ('arrayType', 'i')]

    if 'BinaryArrayOrder' in spectrum.keys():
        if len(spectrum['BinaryArrayOrder']) == 6:
            print spectrum['BinaryArrayOrder'][0][1]
            return spectrum['BinaryArrayOrder'][0][1]
        else:
            try:
                # print 'try 32-bit float'
                spectrum['BinaryArrayOrder'] = [('encoding', '32-bit float'),
                                                ('compression', 'no'), ('arrayType', 'mz'),
                                                ('encoding', '32-bit float'), ('compression', 'no'),
                                                ('arrayType', 'i')]
                if spectrum.highestPeaks(1)[0][1] > 1:
                    # print 'assume:', spectrum['BinaryArrayOrder'][0][1]
                    return spectrum['BinaryArrayOrder'][0][1]

            except:
                # print 'try 64-bit float'
                spectrum['BinaryArrayOrder'] = [('encoding', '64-bit float'),
                                            ('compression', 'no'), ('arrayType', 'mz'),
                                            ('encoding', '64-bit float'), ('compression', 'no'),
                                            ('arrayType', 'i')]
                if spectrum.highestPeaks(1)[0][1] > 1:
                    # print 'assume:', spectrum['BinaryArrayOrder'][0][1]
                    return spectrum['BinaryArrayOrder'][0][1]

    else:
        try:
            # print 'try 64-bit float'
            spectrum['BinaryArrayOrder'] = [('encoding', '64-bit float'),
                                            ('compression', 'no'), ('arrayType', 'mz'),
                                            ('encoding', '64-bit float'), ('compression', 'no'),
                                            ('arrayType', 'i')]
            if spectrum.highestPeaks(1)[0][1] > 1:
                # print 'assume:', spectrum['BinaryArrayOrder'][0][1]
                return spectrum['BinaryArrayOrder'][0][1]

        except:
            # print 'try 32-bit float'
            spectrum['BinaryArrayOrder'] = [('encoding', '32-bit float'),
                                        ('compression', 'no'), ('arrayType', 'mz'),
                                        ('encoding', '32-bit float'), ('compression', 'no'),
                                        ('arrayType', 'i')]
            if spectrum.highestPeaks(1)[0][1] > 1:
                # print 'assume:', spectrum['BinaryArrayOrder'][0][1]
                return spectrum['BinaryArrayOrder'][0][1]


def add_encode_info(spectrum, encode_type):

    if encode_type in ['32-bit float', '64-bit float']:
        pass
    elif encode_type in ['32-bit', '64-bit']:
        if encode_type == '32-bit':
            encode_type = '32-bit float'
        elif encode_type == '64-bit':
            encode_type = '64-bit float'
    elif encode_type in [32, 64]:
        if encode_type == 32:
            encode_type = '32-bit float'
        elif encode_type == 64:
            encode_type = '64-bit float'

    # Add 32-bit/64-bit encoding to mzML
    try:
        # print spectrum['BinaryArrayOrder']
        # [('encoding', '32-bit float'), ('compression', 'zlib'), ('arrayType', 'mz'),
        # ('encoding', '32-bit float'), ('compression', 'zlib'), ('arrayType', 'i')]

        if len(spectrum['BinaryArrayOrder']) == 6:
            pass
        else:
            spectrum['BinaryArrayOrder'] = [('encoding', encode_type),
                                            ('compression', 'no'), ('arrayType', 'mz'),
                                            ('encoding', encode_type), ('compression', 'no'),
                                            ('arrayType', 'i')]

    except KeyError:
        spectrum['BinaryArrayOrder'] = [('encoding', encode_type),
                                        ('compression', 'no'), ('arrayType', 'mz'),
                                        ('encoding', encode_type), ('compression', 'no'),
                                        ('arrayType', 'i')]
        print 'BinaryArrayOrder Added!'

    return spectrum
