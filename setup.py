# -*- coding: utf-8 -*-
"""
Created on Sun Jan 01 01:50:20 2012

@author: Pavlo Shchelokovskyy
"""

from cx_Freeze import setup, Executable

setup(
        name = "pores",
        version = "0.1",
        description = "find pores",
        executables = [Executable("pores.py")])
