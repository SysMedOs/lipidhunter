# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
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
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de
#     Developer Georgia Angelidou georgia.angelidou@uni-leipzig.de

import re


class ElemCalc:
    def __init__(self):

        pa_hg_elem = {"C": 0, "H": 3, "O": 4, "P": 1, "N": 0}
        pc_hg_elem = {"C": 5, "H": 14, "O": 4, "P": 1, "N": 1}
        pe_hg_elem = {"C": 2, "H": 8, "O": 4, "P": 1, "N": 1}
        pg_hg_elem = {"C": 3, "H": 9, "O": 6, "P": 1, "N": 0}
        pi_hg_elem = {"C": 6, "H": 13, "O": 9, "P": 1, "N": 0}
        pip_hg_elem = {"C": 6, "H": 14, "O": 12, "P": 2, "N": 0}
        ps_hg_elem = {"C": 3, "H": 8, "O": 6, "P": 1, "N": 1}
        tg_hg_elem = {"C": 0, "H": 0, "O": 0, "P": 0, "N": 0}
        fa_hg_elem = {"C": 0, "H": 0, "O": 0, "P": 0, "N": 0}

        self.lipid_hg_elem_dct = {
            "PA": pa_hg_elem,
            "PC": pc_hg_elem,
            "PE": pe_hg_elem,
            "PG": pg_hg_elem,
            "PI": pi_hg_elem,
            "PS": ps_hg_elem,
            "PIP": pip_hg_elem,
            "LPA": pa_hg_elem,
            "LPC": pc_hg_elem,
            "LPE": pe_hg_elem,
            "LPG": pg_hg_elem,
            "LPI": pi_hg_elem,
            "LPS": ps_hg_elem,
            "LPIP": pip_hg_elem,
            "TG": tg_hg_elem,
            "FA": fa_hg_elem,
            "DG": tg_hg_elem,
        }

        self.glycerol_bone_elem_dct = {"C": 3, "H": 2}
        self.link_o_elem_dct = {"O": -1, "H": 2}
        self.link_p_elem_dct = {"O": -1}

        self.periodic_table_dct = {
            "H": [(1.0078250321, 0.999885), (2.0141017780, 0.0001157)],
            "D": [(2.0141017780, 0.0001157)],
            "C": [(12.0, 0.9893), (13.0033548378, 0.0107)],
            "N": [(14.0030740052, 0.99632), (15.0001088984, 0.00368)],
            "O": [
                (15.9949146221, 0.99757),
                (16.99913150, 0.00038),
                (17.9991604, 0.00205),
            ],
            "Na": [(22.98976967, 1.0)],
            "P": [(30.97376151, 1.0)],
            "S": [
                (31.97207069, 0.9493),
                (32.97145850, 0.0076),
                (33.96786683, 0.0429),
                (35.96708088, 0.0002),
            ],
            "K": [
                (38.9637069, 0.932581),
                (39.96399867, 0.000117),
                (40.96182597, 0.067302),
            ],
        }

    @staticmethod
    def decode_abbr(abbr):

        pl_checker = re.compile(r"(P[ACEGSI])([(])(.*)([)])")
        lpl_checker = re.compile(r"(LP[ACEGSI])([(])(.*)([)])")
        pip_checker = re.compile(r"(PIP)([(])(.*)([)])")
        tg_checker = re.compile(r"(TG)([(])(.*)([)])")
        dg_checker = re.compile(r"(DG)([(])(.*)([)])")
        fa_checker = re.compile(r"(FA)(\d{1,2})([:])(\d{1,2})")
        fa_short_checker = re.compile(r"(\d{1,2})([:])(\d{1,2})")
        fa_o_checker = re.compile(r"(O-)(\d{1,2})([:])(\d)")
        fa_p_checker = re.compile(r"(P-)(\d{1,2})([:])(\d)")

        # Check PL Type
        _pl_typ = ""
        bulk_fa_typ = ""
        bulk_fa_linker = ""
        bulk_fa_c = 0
        bulk_fa_db = 0
        lyso_fa_linker_dct = {"fa1": "", "fa2": ""}

        # TODO(georgia.angelidou@uni-leipzig.de): consideration to put se elif statement
        if pl_checker.match(abbr):
            # print('PL')
            pl_re_chk = pl_checker.match(abbr)
            pl_typ_lst = pl_re_chk.groups()
            _pl_typ = pl_typ_lst[0]
            bulk_fa_typ = pl_typ_lst[2]
        if lpl_checker.match(abbr):
            # print('PL')
            lpl_re_chk = lpl_checker.match(abbr)
            lpl_typ_lst = lpl_re_chk.groups()
            _pl_typ = lpl_typ_lst[0]
            bulk_fa_typ = lpl_typ_lst[2]
        if pip_checker.match(abbr):
            # print('PIP')
            pip_re_chk = pip_checker.match(abbr)
            pip_typ_lst = pip_re_chk.groups()
            _pl_typ = pip_typ_lst[0]
            bulk_fa_typ = pip_typ_lst[2]
        if tg_checker.match(abbr):
            # print('TG')
            tg_re_chk = tg_checker.match(abbr)
            tg_typ_lst = tg_re_chk.groups()
            _pl_typ = tg_typ_lst[0]
            bulk_fa_typ = tg_typ_lst[2]
        if dg_checker.match(abbr):
            dg_re_chk = dg_checker.match(abbr)
            dg_typ_lst = dg_re_chk.groups()
            _pl_typ = dg_typ_lst[0]
            bulk_fa_typ = dg_typ_lst[2]
        if fa_checker.match(abbr):
            # print('FA')
            _pl_typ = "FA"
            bulk_fa_typ = abbr
            fa_chk = fa_checker.match(abbr)
            bulk_fa_lst = fa_chk.groups()
            bulk_fa_c = bulk_fa_lst[1]
            bulk_fa_db = bulk_fa_lst[3]
            bulk_fa_linker = "A-"
            lyso_fa_linker_dct = {"A": ""}
        if fa_short_checker.match(abbr):
            # print('FA')
            _pl_typ = "FA"
            bulk_fa_typ = abbr
            fa_chk = fa_checker.match(abbr)
            bulk_fa_lst = fa_chk.groups()
            bulk_fa_c = bulk_fa_lst[0]
            bulk_fa_db = bulk_fa_lst[2]
            bulk_fa_linker = "A-"
            lyso_fa_linker_dct = {"A": ""}
        if fa_o_checker.match(abbr):
            # print('FA')
            _pl_typ = "FA"
            bulk_fa_typ = abbr
            fa_chk = fa_o_checker.match(abbr)
            bulk_fa_lst = fa_chk.groups()
            bulk_fa_c = bulk_fa_lst[1]
            bulk_fa_db = bulk_fa_lst[3]
            bulk_fa_linker = "O-"
            lyso_fa_linker_dct = {"O": ""}
        if fa_p_checker.match(abbr):
            # print('FA')
            _pl_typ = "FA"
            bulk_fa_typ = abbr
            fa_chk = fa_p_checker.match(abbr)
            bulk_fa_lst = fa_chk.groups()
            bulk_fa_c = bulk_fa_lst[1]
            bulk_fa_db = bulk_fa_lst[3]
            bulk_fa_linker = "P-"
            lyso_fa_linker_dct = {"P": ""}

        if _pl_typ in ["PL", "PA", "PC", "PE", "PG", "PI", "PIP", "PS"]:
            if fa_short_checker.match(bulk_fa_typ):
                bulk_fa_linker = "A-A-"
                lyso_fa_linker_dct = {"A": ""}
                fa_chk = fa_short_checker.match(bulk_fa_typ)
                bulk_fa_lst = fa_chk.groups()
                bulk_fa_c = bulk_fa_lst[0]
                bulk_fa_db = bulk_fa_lst[2]
            # elif fa_short_checker.match(bulk_fa_typ):
            #     bulk_fa_linker = ''
            #     lyso_fa_linker_dct = {'A': ''}
            #     fa_chk = fa_short_checker.match(bulk_fa_typ)
            #     bulk_fa_lst = fa_chk.groups()
            #     bulk_fa_c = bulk_fa_lst[0]
            #     bulk_fa_db = bulk_fa_lst[2]
            elif fa_o_checker.match(bulk_fa_typ):
                bulk_fa_linker = "O-A-"
                lyso_fa_linker_dct = {
                    "O": "",
                    "A": "O-",
                }  # link of the other sn after NL of this sn
                fa_chk = fa_o_checker.match(bulk_fa_typ)
                bulk_fa_lst = fa_chk.groups()
                bulk_fa_c = bulk_fa_lst[1]
                bulk_fa_db = bulk_fa_lst[3]
            elif fa_p_checker.match(bulk_fa_typ):
                bulk_fa_linker = "P-A-"
                lyso_fa_linker_dct = {
                    "P": "",
                    "A": "P-",
                }  # link of the other sn after NL of this sn
                fa_chk = fa_p_checker.match(bulk_fa_typ)
                bulk_fa_lst = fa_chk.groups()
                bulk_fa_c = bulk_fa_lst[1]
                bulk_fa_db = bulk_fa_lst[3]

        elif _pl_typ in ["LPL", "LPA", "LPC", "LPE", "LPG", "LPI", "LPIP", "LPS"]:
            if fa_short_checker.match(bulk_fa_typ):
                bulk_fa_linker = "A-"
                lyso_fa_linker_dct = {"A": ""}
                fa_chk = fa_short_checker.match(bulk_fa_typ)
                bulk_fa_lst = fa_chk.groups()
                bulk_fa_c = bulk_fa_lst[0]
                bulk_fa_db = bulk_fa_lst[2]
            elif fa_o_checker.match(bulk_fa_typ):
                bulk_fa_linker = "O-"
                lyso_fa_linker_dct = {
                    "O": ""
                }  # link of the other sn after NL of this sn
                fa_chk = fa_o_checker.match(bulk_fa_typ)
                bulk_fa_lst = fa_chk.groups()
                bulk_fa_c = bulk_fa_lst[1]
                bulk_fa_db = bulk_fa_lst[3]
            elif fa_p_checker.match(bulk_fa_typ):
                bulk_fa_linker = "P-"
                lyso_fa_linker_dct = {
                    "P": ""
                }  # link of the other sn after NL of this sn
                fa_chk = fa_p_checker.match(bulk_fa_typ)
                bulk_fa_lst = fa_chk.groups()
                bulk_fa_c = bulk_fa_lst[1]
                bulk_fa_db = bulk_fa_lst[3]

        elif _pl_typ in ["TG"]:
            if fa_short_checker.match(bulk_fa_typ):
                bulk_fa_linker = "A-A-A-"
                lyso_fa_linker_dct = {"A": ""}
                fa_chk = fa_short_checker.match(bulk_fa_typ)
                bulk_fa_lst = fa_chk.groups()
                bulk_fa_c = bulk_fa_lst[0]
                bulk_fa_db = bulk_fa_lst[2]
            elif fa_o_checker.match(bulk_fa_typ):
                bulk_fa_linker = "O-A-A-"
                lyso_fa_linker_dct = {
                    "O": "",
                    "A": "O-",
                }  # link of the other sn after NL of this sn
                fa_chk = fa_o_checker.match(bulk_fa_typ)
                bulk_fa_lst = fa_chk.groups()
                bulk_fa_c = bulk_fa_lst[1]
                bulk_fa_db = bulk_fa_lst[3]
            elif fa_p_checker.match(bulk_fa_typ):
                bulk_fa_linker = "P-A-A-"
                lyso_fa_linker_dct = {
                    "P": "",
                    "A": "P-",
                }  # link of the other sn after NL of this sn
                fa_chk = fa_p_checker.match(bulk_fa_typ)
                bulk_fa_lst = fa_chk.groups()
                bulk_fa_c = bulk_fa_lst[1]
                bulk_fa_db = bulk_fa_lst[3]
        elif _pl_typ in ["DG"]:
            if fa_short_checker.match(bulk_fa_typ):
                bulk_fa_linker = "A-A-"
                lyso_fa_linker_dct = {"A": ""}
                fa_chk = fa_short_checker.match(bulk_fa_typ)
                bulk_fa_lst = fa_chk.groups()
                bulk_fa_c = bulk_fa_lst[0]
                bulk_fa_db = bulk_fa_lst[2]
            elif fa_o_checker.match(bulk_fa_typ):
                bulk_fa_linker = "O-A-"
                lyso_fa_linker_dct = {
                    "O": "",
                    "A": "O-",
                }  # link of the other sn after NL of this sn
                fa_chk = fa_o_checker.match(bulk_fa_typ)
                bulk_fa_lst = fa_chk.groups()
                bulk_fa_c = bulk_fa_lst[1]
                bulk_fa_db = bulk_fa_lst[3]
            elif fa_p_checker.match(bulk_fa_typ):
                bulk_fa_linker = "P-A-"
                lyso_fa_linker_dct = {
                    "P": "",
                    "A": "P-",
                }  # link of the other sn after NL of this sn
                fa_chk = fa_p_checker.match(bulk_fa_typ)
                bulk_fa_lst = fa_chk.groups()
                bulk_fa_c = bulk_fa_lst[1]
                bulk_fa_db = bulk_fa_lst[3]

        bulk_fa_c = int(bulk_fa_c)
        bulk_fa_db = int(bulk_fa_db)

        lipid_info_dct = {
            "TYPE": _pl_typ,
            "LINK": bulk_fa_linker,
            "C": bulk_fa_c,
            "DB": bulk_fa_db,
            "LYSO_LINK": lyso_fa_linker_dct,
        }

        return lipid_info_dct

    def get_neutral_elem(self, abbr):

        usr_lipid_info_dct = self.decode_abbr(abbr)

        lipid_type = usr_lipid_info_dct["TYPE"]

        if lipid_type in list(self.lipid_hg_elem_dct.keys()):
            if lipid_type in ["FA"]:
                # print(abbr)
                tmp_lipid_elem_dct = {
                    "C": usr_lipid_info_dct["C"],
                    "O": 2,
                    "H": (usr_lipid_info_dct["C"] * 2 - usr_lipid_info_dct["DB"] * 2),
                }
                if usr_lipid_info_dct["LINK"] == "":
                    pass
                elif usr_lipid_info_dct["LINK"] == "O-":
                    tmp_lipid_elem_dct["O"] += -1
                    tmp_lipid_elem_dct["H"] += 2
                elif usr_lipid_info_dct["LINK"] == "P-":
                    tmp_lipid_elem_dct["O"] += -1
                else:
                    pass

                return tmp_lipid_elem_dct

            if lipid_type in ["PA", "PC", "PE", "PG", "PI", "PIP", "PS"]:
                tmp_lipid_elem_dct = self.lipid_hg_elem_dct[
                    usr_lipid_info_dct["TYPE"]
                ].copy()
                tmp_lipid_elem_dct["O"] += 4
                tmp_lipid_elem_dct["C"] += (
                    self.glycerol_bone_elem_dct["C"] + usr_lipid_info_dct["C"]
                )
                tmp_lipid_elem_dct["H"] += (
                    self.glycerol_bone_elem_dct["H"]
                    + usr_lipid_info_dct["C"] * 2
                    - usr_lipid_info_dct["DB"] * 2
                )  # DBE = DB + 2xC=O from FA
                if usr_lipid_info_dct["LINK"] == "O-A-":
                    tmp_lipid_elem_dct["O"] += -1
                    tmp_lipid_elem_dct["H"] += 2
                elif usr_lipid_info_dct["LINK"] == "P-A-":
                    tmp_lipid_elem_dct["O"] += -1
                else:
                    pass

                return tmp_lipid_elem_dct

            elif lipid_type in ["LPA", "LPC", "LPE", "LPG", "LPI", "LPIP", "LPS"]:
                tmp_lipid_elem_dct = self.lipid_hg_elem_dct[
                    usr_lipid_info_dct["TYPE"]
                ].copy()
                tmp_lipid_elem_dct["O"] += 3
                tmp_lipid_elem_dct["C"] += (
                    self.glycerol_bone_elem_dct["C"] + usr_lipid_info_dct["C"]
                )
                tmp_lipid_elem_dct["H"] += (
                    self.glycerol_bone_elem_dct["H"]
                    + 2
                    + usr_lipid_info_dct["C"] * 2
                    - usr_lipid_info_dct["DB"] * 2
                )  # DBE = DB + 2xC=O from FA

                if usr_lipid_info_dct["LINK"] == "O-":
                    tmp_lipid_elem_dct["O"] += -1
                    tmp_lipid_elem_dct["H"] += 2
                elif usr_lipid_info_dct["LINK"] == "P-":
                    tmp_lipid_elem_dct["O"] += -1
                else:
                    pass

                return tmp_lipid_elem_dct

            elif lipid_type in ["TG"]:
                tmp_lipid_elem_dct = self.lipid_hg_elem_dct[
                    usr_lipid_info_dct["TYPE"]
                ].copy()
                tmp_lipid_elem_dct["O"] += 6
                tmp_lipid_elem_dct["C"] += (
                    self.glycerol_bone_elem_dct["C"] + usr_lipid_info_dct["C"]
                )
                tmp_lipid_elem_dct["H"] += (
                    self.glycerol_bone_elem_dct["H"]
                    + usr_lipid_info_dct["C"] * 2
                    - usr_lipid_info_dct["DB"] * 2
                )  # DBE = DB + 2xC=O from FA
                if usr_lipid_info_dct["LINK"] == "O-A-A-":
                    tmp_lipid_elem_dct["O"] += -1
                    tmp_lipid_elem_dct["H"] += 2
                elif usr_lipid_info_dct["LINK"] == "P-A-A-":
                    tmp_lipid_elem_dct["O"] += -1
                else:
                    pass

                return tmp_lipid_elem_dct

            elif lipid_type in ["DG"]:
                tmp_lipid_elem_dct = self.lipid_hg_elem_dct[
                    usr_lipid_info_dct["TYPE"]
                ].copy()
                tmp_lipid_elem_dct["O"] += 5
                tmp_lipid_elem_dct["C"] += (
                    self.glycerol_bone_elem_dct["C"] + usr_lipid_info_dct["C"]
                )
                tmp_lipid_elem_dct["H"] += (
                    self.glycerol_bone_elem_dct["H"]
                    + usr_lipid_info_dct["C"] * 2
                    - usr_lipid_info_dct["DB"] * 2
                    + 2
                )

                return tmp_lipid_elem_dct
            else:
                return {"C": 0, "H": 0, "O": 0, "P": 0}
        else:
            return {"C": 0, "H": 0, "O": 0, "P": 0}

    def get_charged_elem(self, abbr, charge="[M-H]-"):

        lipid_elem_dct = self.get_neutral_elem(abbr)
        if charge == "[M-H]-":
            lipid_elem_dct["H"] += -1
        elif charge == "[M+HCOO]-" or charge == "[M+FA-H]-":
            lipid_elem_dct["H"] += 1
            lipid_elem_dct["C"] += 1
            lipid_elem_dct["O"] += 2
        elif charge == "[M+CH3COO]-":
            lipid_elem_dct["H"] += 3
            lipid_elem_dct["C"] += 2
            lipid_elem_dct["O"] += 2
        elif charge == "[M+OAc]-":
            lipid_elem_dct["H"] += 3
            lipid_elem_dct["C"] += 2
            lipid_elem_dct["O"] += 2
        elif charge == "[M+H]+":
            lipid_elem_dct["H"] += 1
        elif charge == "[M+NH4]+":
            if "N" in list(lipid_elem_dct.keys()):
                lipid_elem_dct["N"] += 1
            else:
                lipid_elem_dct["N"] = 1
            lipid_elem_dct["H"] += 4
        elif charge == "[M+Na]+":
            lipid_elem_dct["Na"] = 1

        return lipid_elem_dct

    def get_formula(self, abbr, charge=""):

        if charge in ["neutral", "Neutral", "", None]:

            elem_dct = self.get_neutral_elem(abbr)
        else:
            elem_dct = self.get_charged_elem(abbr, charge=charge)

        formula_str = "C{c}H{h}".format(c=elem_dct["C"], h=elem_dct["H"])

        if "N" in list(elem_dct.keys()):
            if elem_dct["N"] == 1:
                formula_str += "N"
            elif elem_dct["N"] > 1:
                formula_str += "N%i" % elem_dct["N"]

        if "O" in list(elem_dct.keys()):
            if elem_dct["O"] == 1:
                formula_str += "O"
            elif elem_dct["O"] > 1:
                formula_str += "O%i" % elem_dct["O"]

        if "P" in list(elem_dct.keys()):
            if elem_dct["P"] == 1:
                formula_str += "P"
            elif elem_dct["P"] > 1:
                formula_str += "P%i" % elem_dct["P"]

        if "Na" in list(elem_dct.keys()):
            if elem_dct["Na"] == 1:
                formula_str += "Na"
            elif elem_dct["Na"] > 1:
                formula_str += "Na%i" % elem_dct["Na"]

        if "K" in list(elem_dct.keys()):
            if elem_dct["K"] == 1:
                formula_str += "K"
            elif elem_dct["K"] > 1:
                formula_str += "K%i" % elem_dct["K"]

        if charge in ["neutral", "Neutral", "", None]:
            pass
        elif charge in ["[M-H]-", "[M+HCOO]-"]:
            formula_str += "-"
        elif charge in ["[M+H]+", "[M+NH4]+", "[M+Na]+"]:
            formula_str += "+"
        # print ('lets see if you manage to get out from this one')
        return formula_str, elem_dct

    def get_exactmass(self, elem_dct):

        mono_mz = 0.0
        for _elem in list(elem_dct.keys()):
            mono_mz += elem_dct[_elem] * self.periodic_table_dct[_elem][0][0]

        return round(mono_mz, 6)

    # Function to calculate the possible precursor of [M+NH4]+ for TG and DG
    # Current step is working for TG
    # This is needed to understand if the precursors is true [M+H]+ or [M+Na]+ (given the users option)
    def get_NH3_pos_mode(self, charge, mz_pr, amm_elem_dct):
        mz_NH3_pr_H = (
            ""
        )  # The corespond mz of [M+NH4]+ if the given precursor corresponds to the [M+H]+
        mz_NH3_pr_Na = (
            ""
        )  # The corespond mz of [M+NH4]+ if the given precursor corresponds to the [M+Na]+
        if charge in ["[M+H]+"]:
            # Problem need to calculate the theoritical one and not according to the experimental identification
            # mz_NH3_pr_H = mz_pr - self.periodic_table_dct['H'][0][0] + (4 * self.periodic_table_dct['H'][0][0]) + \
            #               self.periodic_table_dct['N'][0][0]
            mz_NH3_pr_H = (
                amm_elem_dct["C"] * self.periodic_table_dct["C"][0][0]
                + amm_elem_dct["H"] * self.periodic_table_dct["H"][0][0]
                + amm_elem_dct["O"] * self.periodic_table_dct["O"][0][0]
                + self.periodic_table_dct["N"][0][0]
            )
            # If this precursor corresponds to the [M+Na]+ then to calculate the bulk identification which it will be different from the above
            # We do the following calculations
            # First remove the charge and the atoms which form the glycerol bond (Na, C3, H5)
            # C3H5Na = 64.028895
            mz_pr_Na = (
                mz_pr
                - self.periodic_table_dct["Na"][0][0]
                - (self.periodic_table_dct["C"][0][0] * 3)
                - (self.periodic_table_dct["H"][0][0] * 5)
            )
            # Second Step: For TG Remove the 6O and the first C from the 3 FA and the last C with the 3 H from each
            # O = 15.994915, C = 12, CH3 = 15.023475
            # 6O, 3xC, 3xCH3 => C6H9O6 = 177.039915
            mz_pr_Na = (
                mz_pr_Na
                - (self.periodic_table_dct["O"][0][0] * 6)
                - (self.periodic_table_dct["C"][0][0] * 6)
                - (self.periodic_table_dct["H"][0][0] * 9)
            )
            mz_pr_Na = int(mz_pr_Na)
            db_count_Na = 0
            while mz_pr_Na % 14 > 1:
                db_count_Na = db_count_Na + 1
                mz_pr_Na = mz_pr_Na - 26
            c_count_Na = int(mz_pr_Na / 14) + 6 + db_count_Na * 2
            tg_abbt_bulk_Na = "TG(" + str(c_count_Na) + ":" + str(db_count_Na) + ")"
            _mz_Na_formula, _mz_Na_elemdct = self.get_formula(
                tg_abbt_bulk_Na, "[M+Na]+"
            )
            mz_NH3_pr_Na = (
                _mz_Na_elemdct["C"] * self.periodic_table_dct["C"][0][0]
                + ((_mz_Na_elemdct["H"] + 4) * self.periodic_table_dct["H"][0][0])
                + _mz_Na_elemdct["O"] * self.periodic_table_dct["O"][0][0]
                + self.periodic_table_dct["N"][0][0]
            )
            # Third step: the rest of the elements correspond to the CH2 * x and if DB to the (CHx2) * y (For 1DB == 2xCH)
            # CH2 = 14.015650, CH = 13.007825
            # Need to create a loop and make all the numbers integers
            # elif charge in ['[M+Na]+']:
            #     mz_NH3_pr_Na = mz_pr - self.periodic_table_dct['Na'][0][0] + (4 * self.periodic_table_dct['H'][0][0]) + \
            #                    self.periodic_table_dct['N'][0][0]

            # TODO(georgia.angelidou@uni-leipzig.de): _mz_H_formula is used before ref, fix here!
            # Temp add

            mz_NH4_H_form = (
                "C" + str(amm_elem_dct["C"]) + "H" + str(amm_elem_dct["H"]) + "O6N"
            )

            # mz_NH4_Na_form = _mz_H_formula
            mz_NH4_Na_form = (
                "C"
                + str(_mz_Na_elemdct["C"])
                + "H"
                + str((_mz_Na_elemdct["H"] + 3))
                + "O6N"
            )

        elif charge in ["[M+Na]+"]:
            mz_NH3_pr_Na = (
                amm_elem_dct["C"] * self.periodic_table_dct["C"][0][0]
                + amm_elem_dct["H"] * self.periodic_table_dct["H"][0][0]
                + amm_elem_dct["O"] * self.periodic_table_dct["O"][0][0]
                + self.periodic_table_dct["N"][0][0]
            )
            print(mz_pr)
            C5H3 = (
                self.periodic_table_dct["H"][0][0] * 6
                + self.periodic_table_dct["C"][0][0] * 3
            )
            rest_sampl = (
                self.periodic_table_dct["O"][0][0] * 6
                + self.periodic_table_dct["H"][0][0] * 9
                + self.periodic_table_dct["C"][0][0] * 6
            )
            mz_pr_H = mz_pr - C5H3 - rest_sampl
            mz_pr_H = int(mz_pr_H)
            db_count_H = 0
            while mz_pr_H % 14 > 1:
                db_count_H = db_count_H + 1
                mz_pr_H = mz_pr_H - 26
            c_count_H = int(mz_pr_H / 14) + 6 + db_count_H * 2
            tg_abbt_bulk_H = "TG(" + str(c_count_H) + ":" + str(db_count_H) + ")"
            _mz_H_formula, _mz_H_elemdct = self.get_formula(tg_abbt_bulk_H, "[M+H]+")
            mz_NH3_pr_H = (
                _mz_H_elemdct["C"] * self.periodic_table_dct["C"][0][0]
                + (_mz_H_elemdct["H"] + 3) * self.periodic_table_dct["H"][0][0]
                + _mz_H_elemdct["O"] * self.periodic_table_dct["O"][0][0]
                + self.periodic_table_dct["N"][0][0]
            )
            mz_NH4_Na_form = (
                "C" + str(amm_elem_dct["C"]) + "H" + str(amm_elem_dct["H"]) + "O6N"
            )
            mz_NH4_H_form = (
                "C"
                + str(_mz_H_elemdct["C"])
                + "H"
                + str((_mz_H_elemdct["H"] + 3))
                + "O6N"
            )

        return (mz_NH3_pr_H, mz_NH4_H_form, mz_NH3_pr_Na, mz_NH4_Na_form)


if __name__ == "__main__":

    usr_bulk_abbr_lst = [
        # 'TG(P-48:2)',
        # 'PC(O-36:3)',
        # 'PC(P-36:3)',
        # 'PC(36:3)',
        "LPC(20:3)"
    ]
    charge_lst = ["[M+NH4]+", "[M-H]-", "[M+HCOO]-", "[M+OAc]-"]
    # usr_bulk_abbr_lst = ['PC(36:3)', 'PC(O-36:3)', 'PC(P-36:3)']
    # charge_lst = ['', '[M-H]-', '[M+HCOO]-', '[M+OAc]-']

    abbr2formula = ElemCalc()

    for usr_abbr in usr_bulk_abbr_lst:
        for _charge in charge_lst:
            usr_formula, usr_elem_dct = abbr2formula.get_formula(
                usr_abbr, charge=_charge
            )
            print((usr_abbr, _charge))
            print(usr_elem_dct)
            print(usr_formula)
