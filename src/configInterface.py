#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Ronja Rösner


"""

#import tkinter for managing GUI
import tkinter as tk
from tkinter import ttk
from configparser import ConfigParser
import os

from .setup import CONFIG_FILE

CONFIG_ADDON_VERSION = "0.0.0"
CONFIG_PARSE = ConfigParser()
CONFIG_KEYS = {
	"API_OPTIONS": {
		"NCBI_API_KEY": "",
		"NCBI_API_MAIL": ""
	}
}

def makeConfigFile(replace: bool=False):
	"""
	Creates a config file for the program.
	"""
	if not os.path.exists(CONFIG_FILE) or replace:
		with open(CONFIG_FILE, 'w') as configfile:
			CONFIG_PARSE.read_dict(CONFIG_KEYS)
			CONFIG_PARSE.write(configfile)
			CONFIG_PARSE.clear()

def getFullConfig():
	with open(CONFIG_FILE, 'r') as file:
		file_content = file.read()
		return file_content

def getSingleConfig(key: str):
	CONFIG_PARSE.read(CONFIG_FILE)
	parsed_data = CONFIG_PARSE.items(key)
	parsed_data = dict(parsed_data)
	CONFIG_PARSE.clear()
	return parsed_data

class OptionsInterface(tk.Toplevel):
	
	def resizeWindow(self, x: int, y: int, min: bool=True, max: bool=True):
		self.geometry(f'{x}x{y}')
		self.minsize(x,y) if min else None
		self.maxsize(x,y) if max else None
	
	def __init__(self, comment: str=""):
		super().__init__()
		self.title(f"Edit Config File {comment}")
		self.resizeWindow(380, 410)
		
		main_frame=tk.Frame(self)
		main_frame.pack(fill='both',expand=True)
		
		makeConfigFile()
		WindowContent(main_frame)


class WindowContent(tk.Frame):
	
	def __init__(self,main_frame):
		super().__init__(main_frame)

		self.option_frame=ttk.Frame(main_frame,height=200,width=900, relief='flat')
		self.option_frame.columnconfigure(0, weight=0)
		self.option_frame.columnconfigure(1, weight=1)
		self.option_frame.columnconfigure(2, weight=1)

		self.option_frame.pack(side='top',fill='x',padx=5,pady=10,ipadx=10, ipady=10)

		self.optionsArea()

	def optionsArea(self):
		self.text_options = ttk.Frame(self.option_frame)
		self.button_row = ttk.Frame(self.option_frame)
		self.text_options.grid(column=0, row=0, padx=10, pady=5, sticky='nsew')
		self.button_row.grid(column=0, row=1, padx=10, pady=5, sticky='ew')
		
		def _place_text():
			self.config_text = tk.Text(self.text_options, wrap='word', height=10, width=10)
			self.config_text.insert('1.0', getFullConfig())
			self.config_text.config(state='disabled')
			self.config_text.grid(column=0, row=0, columnspan=2, padx=10, pady=10, sticky='ew')

		def _update_text():
			self.config_text.config(state='normal')
			self.config_text.delete('1.0', tk.END)
			with open(CONFIG_FILE, 'r') as file:
				file_content = file.read()
				self.config_text.insert(1.0, file_content)
			self.config_text.config(state='disabled')
		
		def _update_entries():
			for main_key in self.entries:
				config_data = getSingleConfig(main_key)
				for sub_key in self.entries[main_key]:
					entry = self.entries[main_key][sub_key]
					config_item = config_data[sub_key.lower()]
					entry.delete(0, tk.END)
					entry.insert(0, config_item)

		def _place_entries():
			entry_dict = {}
			main_index = 1
			for main_key in CONFIG_KEYS:
				main_label = ttk.Label(self.text_options, text=main_key, font=("Arial", 12, "bold"))
				main_label.grid(column=0, row=main_index, columnspan=2, padx=5, pady=5, sticky='ew')
				main_content = getSingleConfig(main_key)
				entry_dict[main_key] = {}
				
				for sub_index, key in enumerate(CONFIG_KEYS[main_key], start=main_index+1):
					ttk.Label(self.text_options, text=key).grid(column=0, row=sub_index, padx=10, pady=5, sticky='w')
					entry = ttk.Entry(self.text_options, textvariable=CONFIG_KEYS[main_key][key])
					entry.insert(0, main_content[key.lower()])
					entry.grid(column=1, row=sub_index, padx=10, pady=5, sticky='ew')
					entry_dict[main_key][key] = entry
				
				main_index = main_index + len(CONFIG_KEYS[main_key].items())+1

			return entry_dict
		
		def _save_entries():
			CONFIG_PARSE.read(CONFIG_FILE)
			for main_key in self.entries:
				config_key = CONFIG_PARSE[main_key]
				for sub_key in self.entries[main_key]:
					config_key[sub_key] = self.entries[main_key][sub_key].get()
			
			with open(CONFIG_FILE, 'w') as file:
				CONFIG_PARSE.write(file)
			
			CONFIG_PARSE.clear()
			_update_text()
		
		def _reset_config():
			makeConfigFile(replace=True)
			_update_text()
			_update_entries()
		
		_place_text()
		self.entries = _place_entries()
		
		self.save_button = ttk.Button(self.button_row, text="Save & Update", command=lambda: _save_entries())
		self.reset_button = ttk.Button(self.button_row, text="Reset Config", command=lambda: _reset_config())

		for widget in self.button_row.winfo_children():
			if 'button' in str(widget):
				widget.pack(side = 'left', fill='x', padx=10, pady=2, expand=True)

if __name__ == "__main__":
	makeConfigFile(replace=False)

	options_window=OptionsInterface("[only for testing]")
	options_window.focus_set()
	options_window.mainloop()