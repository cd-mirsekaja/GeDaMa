#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Ronja Rösner


"""

 
#import tkinter for managing GUI
import tkinter as tk
from tkinter import ttk
from tkinter.filedialog import asksaveasfilename
import os, sqlite3

from .createDatabase import column_names, DB_FILE

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TABLE_ADDON_VERSION = "0.1.1"

class DatabaseInterface(tk.Toplevel):
	
	def resizeWindow(self, x: int, y: int, min: bool=True, max: bool=True):
		self.geometry(f'{x}x{y}')
		self.minsize(x,y) if min else None
		self.maxsize(x,y) if max else None
	
	def __init__(self, comment: str=""):
		super().__init__()
		self.title(f"Database Viewer {comment}")
		self.resizeWindow(830, 900)
		
		main_frame=tk.Frame(self)
		main_frame.pack(fill='both',expand=True)
		
		WindowContent(main_frame)


class WindowContent(tk.Frame):
	
	def __init__(self,main_frame):
		super().__init__(main_frame)
		
		self.table_frame=tk.LabelFrame(main_frame,text="Table",height=540,width=780)
		self.table_frame.pack_propagate(False)
		self.option_frame=tk.LabelFrame(main_frame,text="Options",height=200,width=780)
		self.option_frame.columnconfigure(0, weight=0)
		self.option_frame.columnconfigure(1, weight=1)
		self.option_frame.columnconfigure(2, weight=1)

		self.option_frame.grid(column=1,row=1,sticky='wne',padx=25,pady=10,ipadx=10, ipady=10)
		self.table_frame.grid(column=1,row=2,sticky='enw',padx=25,pady=10)

		db_conn = sqlite3.connect(DB_FILE)
		self.cursor = db_conn.cursor()

		self.optionsArea()
		self.tableView()

	def updateTableView(self):
		selection=self.table_selector.get()
		query=f"SELECT * FROM {selection}"
		self.cursor.execute(query)
		rows=self.cursor.fetchall()
		self.tax_tree.config(columns=(column_names[selection]))
		for i, col in enumerate(column_names[selection]):
			self.tax_tree.heading(i, text=col)

			# set the column width to the longest item in the column
			max_width = len(col) * 12  # start with the header width
			for row in rows:
				cell_value = str(row[i])
				cell_width = len(cell_value) * 8
				if cell_width > max_width:
					max_width = cell_width

			self.tax_tree.column(i, minwidth=0, width=max_width, stretch=False)


		self.tax_tree.delete(*self.tax_tree.get_children())
		for row in rows:
			self.tax_tree.insert("", tk.END, values=row)
		
		self.query_field.delete(0, tk.END)
		self.query_field.insert(0, f"SELECT * FROM {selection}")
	
	def tableView(self):

		tree_scroll_x = ttk.Scrollbar(self.table_frame, orient='horizontal')
		tree_scroll_y = ttk.Scrollbar(self.table_frame,orient='vertical')

		tree_scroll_x.pack(side='bottom',fill='x')
		tree_scroll_y.pack(side='left', fill='y')
		
		self.tax_tree = ttk.Treeview(self.table_frame, xscrollcommand=tree_scroll_x.set, yscrollcommand=tree_scroll_y.set,  show="headings",height=20)
		
		self.tax_tree.pack(side='top',fill='both',expand=True)


		tree_scroll_x.config(command=self.tax_tree.xview)
		tree_scroll_y.config(command=self.tax_tree.yview)

		self.updateTableView()	
		
	def optionsArea(self):
		def executeQuery(query):
			self.cursor.execute(query)
			rows=self.cursor.fetchall()
			
			self.query_output_field.delete(1.0, tk.END)
			self.query_output_field.insert(tk.END, rows)


		tk.Label(self.option_frame,text="Choose Table").grid(column=1,row=0,sticky='nw',padx=10,pady=5)

		self.table_selector=tk.StringVar()
		self.table_selector.set("taxonomy")
		table_menu=tk.OptionMenu(self.option_frame, self.table_selector, *column_names.keys(), command=lambda event: self.updateTableView())
		table_menu.grid(column=1,row=1,sticky='nwe',padx=10,columnspan=2)

		ttk.Separator(self.option_frame,orient='horizontal').grid(column=1,row=2,sticky='we',padx=10,pady=5,columnspan=3)

		tk.Label(self.option_frame,text="Enter SQL Query").grid(column=1,row=3,sticky='nw',padx=10)
		self.query_field=tk.Entry(self.option_frame,width=50)
		self.query_field.grid(column=1,row=4,sticky='nwe',padx=10)

		ttk.Button(self.option_frame,text="Execute Query",command=lambda: executeQuery(self.query_field.get())).grid(column=2,row=4,sticky='nwe',padx=10)

		
		y_scroll_text=ttk.Scrollbar(self.option_frame,orient='vertical')
		y_scroll_text.grid(column=0,row=5,sticky='nes')
		
		self.query_output_field=tk.Text(self.option_frame,height=10)
		self.query_output_field.config(yscrollcommand=y_scroll_text.set)
		self.query_output_field.grid(column=1,row=5,sticky='nwes',padx=5,columnspan=2)
		
		y_scroll_text.config(command=self.query_output_field.yview)

		def switchSelection(selector, selection):
			selector.set(selection)
			self.updateTableView()

		# bind keyboard shortcuts for switching between selections
		self.option_frame.bind_all("<Command-Key-1>", lambda event: switchSelection(self.table_selector, "taxonomy"))
		self.option_frame.bind_all("<Command-Key-2>", lambda event: switchSelection(self.table_selector, "habitats"))
		self.option_frame.bind_all("<Command-Key-3>", lambda event: switchSelection(self.table_selector, "ids"))
		self.option_frame.bind_all("<Return>", lambda event: executeQuery(self.query_field.get()))



if __name__ == "__main__":
	table_window=DatabaseInterface("[only for testing]")
	table_window.focus_set()
	table_window.mainloop()

