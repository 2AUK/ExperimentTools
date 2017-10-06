#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 12:24:18 2017

@author: AbdullahAhmad
"""

import numpy as np

class PDBManipulator(object):
    """
    General class for manipulating PDB files, because having to do it everytime
    is really annoying. Mainly just extending with functionality as needed.
    """
    def __init__(self, path):
        self.path = path
        self.name = ""
        self.sturcture = []

    def readprotein(self):
        """
        PDB file reader. Need a way to read in any PDB file, no matter how crap
        the formatting is. Going to base this on the PDB format standard given
        here:
        ftp://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_A4.pdf
        """
        with open(self.path, 'r') as input_file:
            self.name = input_file.name
            for line in input_file:
                elements = line.split()
                if len(line.strip()) != 0:
                    ident = elements[0]
                    print ident
        