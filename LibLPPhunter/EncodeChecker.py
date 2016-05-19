# -*- coding: utf-8 -*-
# Copyright 2014-16 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig          
# The software is currently  under development and is not ready to be released. 
# A suitable license will be chosen before the official release.               
# For more info please contact zhixu.ni@uni-leipzig.de

import pymzml


def check_encode(mzml):
    """
    This is to solve the bug? from Proteowizard 3.0.9xxx
    The mzML converted from Waters .raw files are without encoder type info.
    This is to test if the file is 32bit or 64bit and rerun the encoder.

    :param mzml: file name of the .mzML file
    :type mzml: str
    :return _encoder: the type of encoder. Can be  ``32-bit float`` or ``64-bit float``
    :type _encoder: str
    """

    # The waters.raw converted mzML by proteowizard may not contain encode info.
    # The '32-bit float' or '64-bit float' will be add to allow following process
    # :param spectrum: the spectrum form pymzML
    # :return: the spectrum with encode info inside

    mzml_obj = pymzml.run.Reader(mzml, MS1_Precision=20e-6, MSn_Precision=200e-6)
    print 'try 32-bit/64-bit float'

    # take 1 scan from mzML
    spectrum = mzml_obj[100]

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
                    _encoder = spectrum['BinaryArrayOrder'][0][1]
                    return _encoder

            except:
                # print 'try 64-bit float'
                spectrum['BinaryArrayOrder'] = [('encoding', '64-bit float'),
                                            ('compression', 'no'), ('arrayType', 'mz'),
                                            ('encoding', '64-bit float'), ('compression', 'no'),
                                            ('arrayType', 'i')]
                if spectrum.highestPeaks(1)[0][1] > 1:
                    # print 'assume:', spectrum['BinaryArrayOrder'][0][1]
                    _encoder = spectrum['BinaryArrayOrder'][0][1]
                    return _encoder

    else:
        try:
            # print 'try 64-bit float'
            spectrum['BinaryArrayOrder'] = [('encoding', '64-bit float'),
                                            ('compression', 'no'), ('arrayType', 'mz'),
                                            ('encoding', '64-bit float'), ('compression', 'no'),
                                            ('arrayType', 'i')]
            if spectrum.highestPeaks(1)[0][1] > 1:
                # print 'assume:', spectrum['BinaryArrayOrder'][0][1]
                _encoder = spectrum['BinaryArrayOrder'][0][1]
                return _encoder

        except:
            # print 'try 32-bit float'
            spectrum['BinaryArrayOrder'] = [('encoding', '32-bit float'),
                                        ('compression', 'no'), ('arrayType', 'mz'),
                                        ('encoding', '32-bit float'), ('compression', 'no'),
                                        ('arrayType', 'i')]
            if spectrum.highestPeaks(1)[0][1] > 1:
                # print 'assume:', spectrum['BinaryArrayOrder'][0][1]
                _encoder = spectrum['BinaryArrayOrder'][0][1]
                return _encoder


def add_encode_info(spectrum, encode_type):
    """
    This is to decode a spectrum from a given encoder type.

    :param spectrum: ``mzML spectrum`` a spectrum from .mzML
    :param encode_type: ``string`` the type of encoder. Can be  ``32-bit float`` or ``64-bit float``
    :type encode_type: str
    :return spectrum: ``mzML spectrum`` the spectrum with correct encoder info inside
    """

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
