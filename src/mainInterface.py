#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Ronja Rösner

This module constructs the main interface for the program.

"""

# import various packages
import os, shutil, threading, numpy
#import tkinter for managing GUI
import tkinter as tk
from tkinter import ttk
from tkinter.filedialog import askopenfilename, asksaveasfilename
from tkinter import messagebox

from .createDatabase import CreateDatabase
from .databaseInterface import DatabaseInterface
from .configInterface import OptionsInterface, makeConfigFile
from .setup import CONFIG_FILE, DB_FILE

class DatabaseMakerInterface(tk.Toplevel):
	def resizeWindow(self, x: int, y: int, min: bool=True, max: bool=True):
		self.geometry(f'{x}x{y}')
		self.minsize(x,y) if min else None
		self.maxsize(x,y) if max else None
	
	def __init__(self, db_file: str, app_title: str, app_version: str=""):
		super().__init__()
		self.title(f'{app_title} {app_version}')
		self.resizeWindow(1220, 580)

		self.main_frame = ttk.Frame(self)
		self.main_frame.pack(fill='both', expand=True)
		
		WindowContent(self.main_frame, module=True, db_file=db_file)

class MainInterface(tk.Tk):
	def resizeWindow(self, x: int, y: int, min: bool=True, max: bool=True):
		self.geometry(f'{x}x{y}')
		self.minsize(x,y) if min else None
		self.maxsize(x,y) if max else None
	
	def __init__(self, app_title: str, app_version: str=""):
		super().__init__()
		self.title(f'{app_title} {app_version}')
		self.resizeWindow(1220, 580, max=False)

		self.main_frame = ttk.Frame(self)
		self.main_frame.pack(fill='both', expand=True)
		
		WindowContent(self.main_frame)


class WindowContent(tk.Frame):
	
	def __init__(self, main_frame, module: bool=False, db_file: str=DB_FILE):
		super().__init__(main_frame)
		self.options_frame = ttk.LabelFrame(main_frame, text = "Options", border = 2, relief = "solid")
		self.options_frame.pack(side = "left", padx = 20, pady = 20, fill = "both")

		self.text_frame = ttk.LabelFrame(main_frame, text = f"Enter Data:", border = 0, relief = "flat")
		self.text_frame.pack(side = "left", padx = 20, pady = 20, fill = "y")

		self.log_frame = ttk.LabelFrame(main_frame, text = "Log", border = 2, relief = "flat")
		self.log_frame.pack(side = "left", padx = 20, pady = 20, fill = "both")

		self.current_thread = None
		self.stop_event = threading.Event()

		self.db_file = db_file

		self.optionsArea()
		self.speciesNameInput()
		self.logOutput()
	
	def optionsArea(self):
		self.general_options_frame = ttk.Frame(self.options_frame)
		self.general_options_frame.pack(side = "top", fill = "x", pady = 10)

		def _switch_state(widget, variable, trigger=0):
			if variable.get() == trigger:
				widget.configure(state="disabled")
			else:
				widget.configure(state="normal")
		
		self.delim_selector = tk.StringVar()
		self.delim_selector.set(";")

		self.button_frame_left = ttk.Frame(self.options_frame)
		self.button_frame_left.pack(side = "bottom", fill = "x", pady = 10)

		# widgets for choosing output saving parameters
		def _output_options():
			self.output_opt_frame = ttk.Frame(self.general_options_frame)

			self.save_db_onoff = tk.IntVar()
			ttk.Checkbutton(self.output_opt_frame, text = "Save output to database file", variable = self.save_db_onoff, onvalue = 1, offvalue = 0,
					command=lambda: _switch_state(append_check, self.save_db_onoff))
			self.save_db_onoff.set(1)

			self.db_append_onoff = tk.IntVar()
			append_check=ttk.Checkbutton(self.output_opt_frame, text = "Append output to existing database", variable = self.db_append_onoff, onvalue = 1, offvalue = 0)
			self.db_append_onoff.set(1)

			"""
			self.itemsort_onoff = tk.IntVar()
			ttk.Checkbutton(self.general_options_frame, text = "Sort output alphabetically", variable = self.itemsort_onoff, onvalue = 1, offvalue = 0)
			self.itemsort_onoff.set(1)
			"""

			# widgets for sorting the output
			ttk.Label(self.output_opt_frame, text = "Sort output by:")
			sorting_options = [
				None,
				'Accession Numbers',
				'Scientific Names',
				]
			self.sorting_selection = tk.StringVar()
			self.sorting_selection.set(sorting_options[0])
			tk.OptionMenu(self.output_opt_frame, self.sorting_selection, *sorting_options)

		# widgets for choosing additional search parameters
		def _search_options():
			self.search_opt_frame = ttk.Frame(self.general_options_frame)
			self.ncbi_search_onoff = tk.IntVar()
			self.ncbi_search_button = ttk.Checkbutton(self.search_opt_frame, text = "Search for NCBI Accession Numbers",
					variable = self.ncbi_search_onoff, onvalue = 1, offvalue = 0)
			self.ncbi_search_onoff.set(1)
		
		# widgets for choosing the input method
		def _input_options():
			def _select_input():
				[_switch_state(radio, self.inputmethod_selection, trigger="Both_On") for radio in self.radiobuttons.values()]
				sel_input = self.inputmethod_selection.get()
				if sel_input == "Both_On":
					self.sorting_selection.set(None)
					self.text_frame.configure(text = f"Enter values as Number{self.delim_selector.get()}Name")
				elif sel_input == "Both_Off":
					self.inputmethod_selection.set(self.input_list[0])
					self.sorting_selection.set(self.input_list[0])
					self.text_frame.configure(text = f"Enter list of {self.inputmethod_selection.get()}")
				else:
					self.text_frame.configure(text = f"Enter list of {self.inputmethod_selection.get()}")
					self.sorting_selection.set(sel_input)
				_switch_state(self.ncbi_search_button, self.inputmethod_selection, trigger="Accession Numbers")
				[_switch_state(combo, self.inputmethod_selection, trigger="Both_Off") for combo in self.comboboxes.values()]
				
			self.input_opt_frame = ttk.Labelframe(self.general_options_frame, text = "Select input method:")
			self.input_list = [
				"Accession Numbers",
				"Scientific Names",
				"Both_On",
				"Both_Off"
			]
			
			self.inputmethod_selection = tk.StringVar()
			self.inputmethod_selection.set(self.input_list[2])
			self.radiobuttons, self.comboboxes = {}, {}
			ttk.Checkbutton(self.input_opt_frame, text = "Both values", variable=self.inputmethod_selection, onvalue="Both_On", offvalue="Both_Off",
					command=lambda: _select_input())
			
			self.radiobuttons[self.input_list[0]] = ttk.Radiobutton(self.input_opt_frame, text = self.input_list[0],
					variable = self.inputmethod_selection, value = self.input_list[0], command = lambda: _select_input())
			self.radiobuttons[self.input_list[1]] = ttk.Radiobutton(self.input_opt_frame, text = self.input_list[1],
					variable = self.inputmethod_selection, value = self.input_list[1], command = lambda: _select_input())
			_select_input()

			ttk.Separator(self.input_opt_frame)

			# widgets for choosing the item delimiter
			ttk.Label(self.input_opt_frame, text = "Set item delimiter:")
			delim_tuple = (
				";",
				",",
				"|"
			)
			self.comboboxes["delimiter"] = ttk.Combobox(self.input_opt_frame, textvariable = self.delim_selector)
			self.comboboxes["delimiter"]["values"] = delim_tuple
			self.comboboxes["delimiter"].current(0)


		# call the functions for widget creation
		_output_options()
		ttk.Separator(self.general_options_frame)
		_search_options()
		ttk.Separator(self.general_options_frame)
		_input_options()


		# place all general options widgets on the interface
		for widget in self.general_options_frame.winfo_children():
			if widget.winfo_class() == "TFrame" or widget.winfo_class() == "TLabelframe":
				widget.pack(fill="both", padx=10, pady=5, expand=True)
				for sub_widget in widget.winfo_children():
					if sub_widget.winfo_class() == "TLabel":
						sub_widget.pack(fill = "x", padx = 10, pady = 2.5, side="left")
					else:
						sub_widget.pack(fill = "x", padx = 10, pady = 2.5)
			else:
				widget.pack(fill = "x", padx = 10, pady = 5)

		# set list with possible text formats
		textfile_types=[
			("Simple Text Files", "*.txt"),
			("Complex Text Files", "*.rtf"),
			("Comma Separated Values", "*.csv"),
			("Markdown Files", ".md")
			]

		# function for opening the options window
		self.config_window_open=False
		def _view_config():
				if not os.path.isfile(CONFIG_FILE):
					messagebox.showwarning("File missing", "No config found, empty file will be initialized.")
					makeConfigFile()

				# function for destroying the window after it has been closed
				def _onConfigClose():
					self.config_window.destroy()
					self.config_window_open=False
				
				if not self.config_window_open:
					self.config_window=OptionsInterface()
					self.config_window.protocol('WM_DELETE_WINDOW',lambda: _onConfigClose())
					self.config_window_open=True
				else:
					self.config_window.focus_set()
		
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

		ttk.Separator(self.button_frame_left)
		ttk.Button(self.button_frame_left, text = "Change Config", command = lambda: _view_config())
		ttk.Separator(self.button_frame_left)
		ttk.Button(self.button_frame_left, text = "Upload Text", command = lambda: _open_file())
		ttk.Button(self.button_frame_left, text = "Save Text", command = lambda: _save_file())
		ttk.Button(self.button_frame_left, text = "Remove Text", command = lambda: _reset())

		for item in self.button_frame_left.winfo_children():
			item.pack(fill = "x", padx = 10, pady = 1)

	def speciesNameInput(self):

		self.button_frame_middle = ttk.Frame(self.text_frame)
		self.button_frame_middle.pack(side = "bottom", fill = "x", pady = 10)

		self.text_field = tk.Text(self.text_frame, width = 40, height = 27)

		def _confirm():
			_cancel()
			self.stop_event.clear()

			# display confirmation overlay
			conti = messagebox.askyesno("Continue?", "Do you want to start a search with your input?")
			if not conti:
				return

			delim = self.delim_selector.get()
			input_list = self.text_field.get(1.0, tk.END)
			input_list = [s for s in input_list.split("\n") if s]
			input_list = [item.strip() for item in input_list]
			
			def _background_task():
				input_sorting = self.sorting_selection.get()
				selected_input = self.inputmethod_selection.get()
				if selected_input == "Scientific Names":
					accessionNumberList = []
					sciNameList = sorted(input_list) if input_sorting == "Scientific Names" else input_list
				elif selected_input == "Accession Numbers":
					accessionNumberList = sorted(input_list) if input_sorting == "Accession Numbers" else input_list
					sciNameList = []
				elif selected_input == "Both_On":
					
					accessionNumberList = [item.split(delim)[0] for item in input_list]
					sciNameList = [item.split(delim)[1] for item in input_list]

					if input_sorting == "Accession Numbers":
						indices = numpy.argsort(accessionNumberList)
						accessionNumberList = [accessionNumberList[i] for i in indices]
						sciNameList = [sciNameList[i] for i in indices]
					elif input_sorting == "Scientific Names":
						indices = numpy.argsort(sciNameList)
						accessionNumberList = [accessionNumberList[i] for i in indices]
						sciNameList = [sciNameList[i] for i in indices]

				out_dict = {
					"sciNameList": sciNameList,
					"accessionNumberList": accessionNumberList,
					"getMissing": self.ncbi_search_onoff.get(),
					"append_data": self.db_append_onoff.get()
				}
				sql = CreateDatabase(out_dict, db_file=self.db_file, log_function=self.updateLog, stop_event=self.stop_event)

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
			if not os.path.isfile(self.db_file):
				messagebox.showwarning("File missing", "No database file found!")
				return

			filepath=asksaveasfilename(filetypes=[("Database File","*.db")])
			if not filepath:
				return
			
			shutil.copyfile(self.db_file, filepath)
		
		# function for opening the database viewer
		self.database_window_open=False
		def _view_database():
				if not os.path.isfile(self.db_file):
					messagebox.showwarning("File missing", "No database file found!")
					return

				# function for destroying the window after it has been closed
				def _onDatabaseClose():
					self.database_window.destroy()
					self.database_window_open=False
				
				if not self.database_window_open:
					self.database_window=DatabaseInterface(db_file=self.db_file)
					self.database_window.protocol('WM_DELETE_WINDOW',lambda: _onDatabaseClose())
					self.database_window_open=True
				else:
					self.database_window.focus_set()

		# create buttons, listed in reverse order of appearance
		ttk.Button(self.button_frame_middle, text = "Confirm Search", command=lambda: _confirm())
		ttk.Button(self.button_frame_middle, text = "Cancel Search", command=lambda: _cancel())
		
		ttk.Button(self.button_frame_middle, text = "View Current Database", command=lambda: _view_database())
		ttk.Button(self.button_frame_middle, text = "Export Current Database", command = lambda: _export_database())

		self.text_field.pack(fill = "y", pady = 5, padx = 10)
		ttk.Separator(self.button_frame_middle).grid(column=1, row=0, columnspan=2, sticky="ew", padx=5, pady=10)
		# place middle buttons
		button_idx = 1
		for button in self.button_frame_middle.winfo_children():
			if "!button" in str(button):
				if button_idx<=2:
					button.grid(column=button_idx, row=1, padx=10, pady=0.1, sticky="ew")
					button_idx = button_idx+1
				elif button_idx>=3:
					button.grid(column=1, row=button_idx-1, padx=10, pady=0.1, columnspan=2, sticky="ew")
					button_idx = button_idx+1
	
	def logOutput(self):
		self.log_field = tk.Text(self.log_frame, width = 75, height = 40, state = "disabled", cursor = "cross")
		self.log_field.pack(fill = "y")
	
	def updateLog(self, message, end: str="\n"):
		self.log_field.configure(state = "normal")
		self.log_field.insert(tk.END, f"{message}{end}")
		self.log_field.see(tk.END)  # Scroll to the end
		self.log_field.configure(state = "disabled")


if __name__ == "__main__":
	main_window = MainInterface("Taxonda", "[only for testing]")
	main_window.focus_set()
	main_window.mainloop()
