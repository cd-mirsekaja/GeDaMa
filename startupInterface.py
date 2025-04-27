#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Ronja Rösner


"""

 
#import tkinter for managing GUI
import tkinter as tk
from tkinter import ttk
from tkinter.filedialog import asksaveasfilename
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

class StartupInterface(tk.Toplevel):
	
	def resizeWindow(self, x: int, y: int):
		self.geometry(f'{x}x{y}')
		self.minsize(x,y)
		self.maxsize(x,y)
	
	def __init__(self, comment: str=""):
		super().__init__()
		self.title(f"Welcome!")
		self.resizeWindow(830, 900)
		
		main_frame=tk.Frame(self)
		main_frame.pack(fill='both',expand=True)
		
		WindowContent(main_frame)

class WindowContent(tk.Frame):
	def __init__(self,main_frame):
		super().__init__(main_frame)
		
