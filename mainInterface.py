#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Ronja Rösner

This module constructs the main interface for the program.

"""

# import various packages
import os, shutil, threading
#import tkinter for managing GUI
import tkinter as tk
from tkinter import ttk
from tkinter.filedialog import askopenfilename, asksaveasfilename
from tkinter import messagebox
# import custom methods
from .createDatabase import CreateDatabase, DB_FILE
from .downloadSpeciesData import getScientificNames
from .databaseInterface import DatabaseInterface
from .optionsInterface import OptionsInterface, CONFIG_FILE, makeConfigFile

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

class DatabaseMakerInterface(tk.Toplevel):
	def resizeWindow(self, x: int, y: int, min: bool=True, max: bool=True):
		self.geometry(f'{x}x{y}')
		self.minsize(x,y) if min else None
		self.maxsize(x,y) if max else None
	
	def __init__(self, app_title: str, app_version: str=""):
		super().__init__()
		self.title(f'{app_title} {app_version}')
		self.resizeWindow(1000, 550)

		self.main_frame = ttk.Frame(self)
		self.main_frame.pack(fill='both', expand=True)
		
		WindowContent(self.main_frame)

class MainInterface(tk.Tk):
	def resizeWindow(self, x: int, y: int, min: bool=True, max: bool=True):
		self.geometry(f'{x}x{y}')
		self.minsize(x,y) if min else None
		self.maxsize(x,y) if max else None
	
	def __init__(self, app_title: str, app_version: str=""):
		super().__init__()
		self.title(f'{app_title} {app_version}')
		self.resizeWindow(1000, 550)

		self.main_frame = ttk.Frame(self)
		self.main_frame.pack(fill='both', expand=True)
		
		WindowContent(self.main_frame)


class WindowContent(tk.Frame):
	
	def __init__(self, main_frame):
		super().__init__(main_frame)
		self.log_frame = ttk.LabelFrame(main_frame, text = "Log", border = 2, relief = "flat")
		self.log_frame.pack(side = "right", padx = 20, pady = 20, fill = "both")

		self.text_frame = ttk.LabelFrame(main_frame, text = f"Enter Data:", border = 0, relief = "flat")
		self.text_frame.pack(side = "right", padx = 20, pady = 20, fill = "y")

		self.options_frame = ttk.LabelFrame(main_frame, text = "Options", border = 2, relief = "solid")
		self.options_frame.pack(side = "left", padx = 20, pady = 20, fill = "both")

		self.current_thread = None
		self.stop_event = threading.Event()

		self.optionsArea()
		self.speciesNameInput()
		self.logOutput()
	
	def optionsArea(self):
		self.general_options_frame = ttk.Frame(self.options_frame)
		self.general_options_frame.pack(side = "top", fill = "x", pady = 10)

		def _switch_state(widget):
			if self.save_db_onoff.get() == 1:
				widget.configure(state="normal")
			else:
				widget.configure(state="disabled")
		
		self.button_frame_right = ttk.Frame(self.options_frame)
		self.button_frame_right.pack(side = "bottom", fill = "x", pady = 10)

		self.save_db_onoff = tk.IntVar()
		ttk.Checkbutton(self.general_options_frame, text = "Save output to database file", variable = self.save_db_onoff, onvalue = 1, offvalue = 0, command=lambda: _switch_state(append_button))
		self.save_db_onoff.set(1)

		# when appending and not overwriting, the IDX gets reset. Needs to be fixed.
		self.db_append_onoff = tk.IntVar()
		append_button=ttk.Checkbutton(self.general_options_frame, text = "Append output to existing database", variable = self.db_append_onoff, onvalue = 1, offvalue = 0)
		self.db_append_onoff.set(0)

		self.itemsort_onoff = tk.IntVar()
		ttk.Checkbutton(self.general_options_frame, text = "Sort output alphabetically", variable = self.itemsort_onoff, onvalue = 1, offvalue = 0)
		self.itemsort_onoff.set(1)

		ttk.Separator(self.general_options_frame)

		self.ncbi_search_onoff = tk.IntVar()
		ttk.Checkbutton(self.general_options_frame, text = "Search for NCBI Accession Numbers", variable = self.ncbi_search_onoff, onvalue = 1, offvalue = 0)
		self.ncbi_search_onoff.set(1)

		ttk.Separator(self.general_options_frame)

		# widgets for choosing the input method
		def _select_input():
			self.text_frame.configure(text = f"Enter {self.inputmethod_selection.get()}:")
			self.ncbi_search_onoff.set(0) if self.inputmethod_selection.get() == "Accession Numbers" else self.ncbi_search_onoff.set(1)
		
		ttk.Label(self.general_options_frame, text = "Select input method:")
		input_list = [
			"Accession Numbers",
			"Scientific Names"
		]
		self.inputmethod_selection = tk.StringVar()
		self.inputmethod_selection.set(input_list[1])
		_select_input()
		for input in input_list:
			ttk.Radiobutton(self.general_options_frame, text = input, variable = self.inputmethod_selection, value = input, command = lambda: _select_input())

		ttk.Separator(self.general_options_frame)

		# widgets for choosing the item delimiter
		ttk.Label(self.general_options_frame, text = "Set item delimiter:")
		delim_tuple = (
			"Linebreak",
			",",
			";"
		)
		self.delim_selector = tk.StringVar()
		choose_delimiter = ttk.Combobox(self.general_options_frame, textvariable = self.delim_selector)
		choose_delimiter["values"] = delim_tuple
		choose_delimiter.current(0)

		# place all general options widgets on the interface
		for widget in self.general_options_frame.winfo_children():
			if ".!frame.!labelframe3.!frame.!label" in str(widget):
				widget.pack(fill = "x", padx = 20, pady = 0)
			elif ".!frame.!labelframe3.!frame.!checkbutton" in str(widget):
				widget.pack(fill = "x", padx = 20, pady = 2.5)
			elif ".!frame.!labelframe3.!frame.!radiobutton" in str(widget):
				widget.pack(fill = "x", padx = 20, pady = 2.5)
			else:
				widget.pack(fill = "x", padx = 20, pady = 5)

		# set list with possible text formats
		textfile_types=[
			("Simple Text Files", "*.txt"),
			("Complex Text Files", "*.rtf"),
			("Markdown Files", ".md")
			]

		# function for opening the options window
		self.options_window_open=False
		def _view_config():
				if not os.path.isfile(CONFIG_FILE):
					messagebox.showwarning("File missing", "No config file found, empty file will be initialized.")
					makeConfigFile()

				# function for destroying the window after it has been closed
				def onTableClose():
					self.options_window.destroy()
					self.options_window_open=False
				
				if not self.options_window_open:
					self.options_window=OptionsInterface()
					self.options_window.protocol('WM_DELETE_WINDOW',lambda: onTableClose())
					self.options_window_open=True
				else:
					self.options_window.focus_set()
		
		# function for uploading a file to the data folder
		def _upload_file():
			filepath=askopenfilename(filetypes=textfile_types)
			if not filepath:
				return
			
			shutil.copy(filepath, os.path.join(SCRIPT_DIR, "data/"))

		# function for uploading a text file into the text field
		def _open_file():
			filepath=askopenfilename(filetypes=textfile_types)
			if not filepath:
				return
			
			self.text_field.delete(1.0, tk.END)
			with open(filepath,"r") as file:
				content = file.read()
				self.text_field.insert(tk.END,content)
		
		# function for saving the current text field to a file
		def _save_file():
			filepath=asksaveasfilename(filetypes=textfile_types)
			if not filepath:
				return
			
			with open(filepath,"w") as file:
				content=self.text_field.get(1.0,tk.END)
				file.write(content)

		# function for resetting the text field
		def _reset():
			self.text_field.delete(1.0, tk.END)


		ttk.Separator(self.button_frame_right)
		ttk.Button(self.button_frame_right, text = "Advanced Options", command = lambda: _view_config())
		ttk.Separator(self.button_frame_right)
		#ttk.Button(self.button_frame_right, text = "Upload Accession Numbers", command = lambda: _upload_file())
		ttk.Button(self.button_frame_right, text = "Upload Text File", command = lambda: _open_file())
		ttk.Button(self.button_frame_right, text = "Save Text File", command = lambda: _save_file())
		ttk.Button(self.button_frame_right, text = "Remove Text", command = lambda: _reset())

		for item in self.button_frame_right.winfo_children():
			item.pack(fill = "x", padx = 20, pady = 2.5)

	def speciesNameInput(self):

		self.button_frame_middle = ttk.Frame(self.text_frame)
		self.button_frame_middle.pack(side = "bottom", fill = "x", pady = 10)

		self.text_field = tk.Text(self.text_frame, width = 30, height = 27)

		def _confirm():
			_cancel()
			self.stop_event.clear()

			# display confirmation overlay
			conti = messagebox.askyesno("Continue?", "Do you want to start a search with your input?")
			if not conti:
				return

			delim = self.delim_selector.get() if self.delim_selector.get() != "Linebreak" else "\n"
			input_list = self.text_field.get(1.0, tk.END)
			input_list = [s for s in input_list.split(delim) if s]
			input_list = [item.strip() for item in input_list]
			input_list = sorted(input_list) if self.itemsort_onoff.get() == 1 else input_list

			def _background_task():
				if self.inputmethod_selection.get() == "Scientific Names":
					sciNameList = input_list
					accessionNumberList = []
				elif self.inputmethod_selection.get() == "Accession Numbers":
					self.updateLog("\n=== Now starting search for scientific names ===")
					accessionNumberList = input_list
					sciNameList = getScientificNames(input_list, log_function=self.updateLog, stop_event=self.stop_event)["all"]
					if self.stop_event.is_set():
						self.updateLog("\nSearch cancelled during name retrieval.\n")
						return
	
				out_dict = {
					"sciNameList": sciNameList,
					"accessionNumberList": accessionNumberList,
					"getMissing": self.ncbi_search_onoff.get(),
					"append_data": self.db_append_onoff.get()
				}
				sql = CreateDatabase(out_dict, self.updateLog, stop_event=self.stop_event)

				if self.save_db_onoff.get() == 1:
					self.updateLog("\n=== Now starting general data download ===")
					sql.makeSQL()
				else:
					self.updateLog("\n=== Now starting general data download ===")
					data_dict = sql.makeDataframes()
					for key in data_dict:
						self.updateLog(f"{key}:\n{data_dict[key]}")

				self.updateLog("\nSearch completed!\n")
				self.updateLog("---------------------------\n")
			
			# start the search in the background so that the interface does not freeze
			self.current_thread = threading.Thread(target = _background_task)
			self.current_thread.daemon = True
			self.current_thread.start()

		def _cancel():
			if self.current_thread and self.current_thread.is_alive():
				self.stop_event.set()
				self.current_thread.join(timeout=0.5)
				self.updateLog("\n--- Search cancelled. ---\n")

		def _export_database():
			if not os.path.isfile(DB_FILE):
				messagebox.showwarning("File missing", "No database file found!")
				return

			filepath=asksaveasfilename(filetypes=[("Database File","*.db")])
			if not filepath:
				return
			
			shutil.copyfile(DB_FILE, filepath)
		
		# function for opening the database viewer
		self.table_window_open=False
		def _view_database():
				if not os.path.isfile(DB_FILE):
					messagebox.showwarning("File missing", "No database file found!")
					return

				# function for destroying the window after it has been closed
				def onTableClose():
					self.table_window.destroy()
					self.table_window_open=False
				
				if not self.table_window_open:
					self.table_window=DatabaseInterface()
					self.table_window.protocol('WM_DELETE_WINDOW',lambda: onTableClose())
					self.table_window_open=True
				else:
					self.table_window.focus_set()

		# create buttons, listed in reverse order of appearance
		ttk.Button(self.button_frame_middle, text = "Export Current Database", command = lambda: _export_database())
		ttk.Button(self.button_frame_middle, text = "View Current Database", command=lambda: _view_database())
		ttk.Button(self.button_frame_middle, text = "Cancel Search", command=lambda: _cancel())
		ttk.Button(self.button_frame_middle, text = "Confirm Search", command=lambda: _confirm())

		self.text_field.pack(fill = "y", pady = 5, padx = 5)
		# place widgets in the text frame
		for widget in self.button_frame_middle.winfo_children():
			if "!button" in str(widget):
				widget.pack(fill = "x", padx = 10, pady = 0.1, side = "bottom")
			else:
				widget.pack()
	
	def logOutput(self):
		self.log_field = tk.Text(self.log_frame, width = 50, height = 37, state = "disabled", cursor = "cross")
		self.log_field.pack(fill = "y")
	
	def updateLog(self, message):
		self.log_field.configure(state = "normal")
		self.log_field.insert(tk.END, message + "\n")
		self.log_field.see(tk.END)  # Scroll to the end
		self.log_field.configure(state = "disabled")





if __name__ == "__main__":
	main_window = MainInterface("Genome Database Maker", "[only for testing]")
	main_window.focus_set()
	main_window.mainloop()
