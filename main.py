#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Ronja Rösner

This module brings the other modules together.

"""

# import custom functions for constructing interface
from .mainInterface import MainInterface

# import the program name and version from the setup file
from .setup import NAME, VERSION

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