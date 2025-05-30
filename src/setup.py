"""
This is a setup.py script generated by py2applet

Usage:
	python setup.py py2app
"""

from setuptools import setup
import os

NAME = ['GeDaMa']
VERSION = ['0.2.1']
APP = ['main.py']
DATA_FILES = [
	'species_database.db'
]
OPTIONS = {
	"iconfile": "images/icon.icns",
	'argv_emulation': True,
	"packages": ['pandas','pygbif']
	}

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR = os.path.dirname(SCRIPT_DIR)
CONFIG_FILE = f"{PARENT_DIR}/data/config.ini"
DB_FILE = f"{PARENT_DIR}/data/species_database.db"


if __name__=='__main__':
	setup(
	    app=APP,
	    data_files=DATA_FILES,
	    options={'py2app': OPTIONS},
	    setup_requires=['py2app'],
	)