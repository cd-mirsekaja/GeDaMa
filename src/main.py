#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Ronja Rösner

This module brings the other modules together.

"""
import os, sys
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_dir)

from src.mainInterface import MainInterface
from src.setup import NAME, VERSION

# turn the imported program info into strings
program_name=str(NAME[0])
program_version=str(VERSION[0])

# function for constructing the application
def main():
	main_window=MainInterface(program_name, program_version)
	main_window.focus_set()
	main_window.mainloop()


if __name__ == "__main__":
	main()