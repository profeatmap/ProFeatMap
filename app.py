# -*- coding: utf-8 -*-
# python 3.9

# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

# Modules to install

# Reading xlsx and odf files
# openpyxl 3.0.7 for reading excel files as DataFrames (should not be necessary anymore)
# odfpy 1.4.1 for opening ods files as DataFrames

# Dash
import dash # 2.0.0
from dash.dash import no_update
from dash import dcc # dash 2.0.0
from dash import html # dash 2.0.0
from dash import dash_table # 5.0.0 included in dash 2.0.0
from dash.dependencies import Input, Output, State, MATCH, ALL
import numpy as np # 1.19.1 included in dash 2.0.0
import pandas as pd # 1.4.1

# Dash Bootstrap Components (GUI)
import dash_bootstrap_components as dbc # dash-bootstrap-components 1.0.3

# Downloading files
# dash-extensions : https://stackoverflow.com/questions/61784556/download-csv-file-in-dash/61814287
from dash_extensions import Download # 0.0.60
from dash_extensions.snippets import send_data_frame
from dash_extensions.snippets import send_bytes

# Base Python
import random as rd
import math
import base64
import re
import io
import subprocess
import os
import copy
import time
import itertools
from datetime import datetime

# Downloading data from Uniprot
import urllib.request # urllib3 1.25.9

# Map drawing
from PIL import Image, ImageDraw, ImageFont # pillow 7.2.0

# Colormaps handling
import matplotlib as plt # 3.5.1
import pylab as pl

# Creating xlsx files
import xlsxwriter # 1.2.9

# Generating unique ID
import uuid # 1.30



def b64_image(image_filename):
	"""
		DESCRIPTION: used to display images in dash application
	"""
	with open(image_filename, 'rb') as f:
		image = f.read()

	return 'data:image/png;base64,' + base64.b64encode(image).decode('utf-8')

def getUniprotDBData(protein_list) :
	"""
		DESCRIPTION: Reads actual state of the uniprot "database" file as dict
	"""
	file = None
	uniprot_db = {}
	# Try to open file
	try :
		for line in open('uniprot_db.csv', 'r') :
			sline = line.strip('\n').split('\t')
			if sline[0] != 'Uniprot_code' :
				uniprot_code = sline[0]
				latest_download = sline[1]
				download_nb = int(sline[2])
				uniprot_db[uniprot_code] = {'latest' : latest_download, 'nb_dl' : download_nb}
	# Creates an empty database files if not file was found
	except :
		file = open('uniprot_db.csv', 'w')
		file.write('Uniprot_code\tlatest_download\tdownload_nb\n')
		file.close()

	return uniprot_db

def updateUniprotDBData(uniprot_db) :
	"""
		DESCRIPTION: Updates uniprot "database" file with latests downloads and download number
	"""
	output = open('uniprot_db.csv', 'w')
	output.write('Uniprot_code\tlatest_download\tdownload_nb\n')
	for key, values in uniprot_db.items() :
		output.write(key+'\t'+values['latest']+'\t'+str(values['nb_dl'])+'\n')

	output.close()

def updateTargetsUniprotDBData(uniprot_db, protein_list, outdated_thres=3) :
	"""
		DESCRIPTION: Return list of all uniprot files that needs to be replaced.
			Need to be replaced because to old (>3 months by default) or forced update by the user
	"""

	# Search for all proteins not already present in the 'database'
	new_protein_list = list(set(protein_list) - set(list(uniprot_db.keys())))
	# Get the list of all protein where information is already in the 'database'
	already_downloaded_data = list(set(protein_list).intersection(list(uniprot_db.keys())))

	# Initializing the list containing all proteins that should be newly downloaded or updated.
	to_update_protein_list = new_protein_list

	# Get todays date (year, month and day)
	today_date = datetime.today().strftime('%Y-%m-%d')
	year, month, day = [int(elem) for elem in today_date.split('-')]

	# For all proteins already in database
	for protein in already_downloaded_data :
		# Get latest date
		db_year, db_month, db_day = [int(elem) for elem in uniprot_db[protein]['latest'].split('-')]
		# Check if date in database is older then defined threshold
		if year * month - db_year * db_month >= outdated_thres :
			# If existing file is older, add to the download list
			to_update_protein_list.append(protein)

	return to_update_protein_list

def parse_content(content, filename):
	"""
		DESCRIPTION: Function that opens files delivered by drag and drops as DataFrames.
	"""

	# Get the content of the file as string
	content_type, content_string = content.split(',')

	# Decode files for coded formats (.xls, .xlsx, .ods)
	decoded = base64.b64decode(content_string)
	#Get the file extension
	file_extension = filename.split('.')[-1]

	try:
		# Non coded file formats
		if file_extension in ['csv', 'tsv', 'txt', 'tab'] :
			# Assume that the user uploaded a CSV file
			file_str = decoded.decode('utf-8')

			if '\t' in file_str :
				print('Tab separator detected')
				df = pd.read_csv(io.StringIO(file_str), sep='\t')
			elif ';' in file_str :
				print('; detected')
				file_str = file_str.replace(',', '.') # here we assume that ; csv files have decimal separators as commas
				df = pd.read_csv(io.StringIO(file_str), sep=';')
			else :
				print('Comma separator detected')
				df = pd.read_csv(io.StringIO(file_str), sep=',')

		# Encoded file formats
		elif file_extension == 'xlsx' :
			# Assume that the user uploaded an excel file
			# Opening of xlsx files isn't working properly of different OS. Two different approaches are tryed here.
			try :
				print('Trying to open without openpyxl')
				df = pd.read_excel(io.BytesIO(decoded), engine=None) # Uses the default engine
			except :
				print('Trying to open with openpyxl')
				df = pd.read_excel(io.BytesIO(decoded), engine='openpyxl') # Uses the outdated openpyxl engine

		elif file_extension == 'xls' :
			# Assume that the user uploaded an excel file
			df = pd.read_excel(io.BytesIO(decoded))

		elif file_extension == 'ods' :
			# Assume that the user uploaded an ods file
			df = pd.read_excel(io.BytesIO(decoded), engine='odf')

	# If loading of file fails
	except Exception as e:
		print(e)
		return html.Div([
			'There was an error processing this file.'
		])

	return df

def extractData(base_protein_list, df, FT_order_list, case=None):
	"""
		DESCRIPTION: Extracts information for all proteins in a list (base_protein_list) from a pandas FataFrame (df) and orders the features according to the order list (FT_order_lsit)
			FT_order_list can be changed by the user to avoid overlaps of features as the the order in the dict is used in the map drawing process.
	"""

	# Initialize 
	data_dict = {}
	feature_index = 0

	# Initialize dictionaries for each protein
	for protein_name in base_protein_list :
		data_dict[protein_name] = {}

	# Fill dictionaries
	for index, row in df.iterrows():
		protein_name = row['protein']

		# Extract data when no values are given
		if case == 'noCase' :
			data_dict[protein_name][feature_index] = [str(row[elem]) if str(row[elem]) != 'nan' else '' for elem in ('feature_type', 'feature', 'start', 'length', 'intensity')]
		# Extract data for a specific case (for values)
		else :
			data_dict[protein_name][feature_index] = [str(row[elem]) if str(row[elem]) != 'nan' else '' for elem in ('feature_type', 'feature', 'start', 'length', case)]

		feature_index += 1

	# Sort data given the FT_order_list
	new_dict = sortDataDict(data_dict, FT_order_list)

	return new_dict

def sortDataDict(data_dict, FT_order_list) :
	"""
		DESCRIPTION: From an exsisting data_dict, creates a new dict sorted in the order defined in the FT_order_list.
			The order is used in the map drawing process.
	"""

	# Initialize new dict
	new_dict = {}

	# If the FT_order_list if not given, do not change the order.
	if FT_order_list == [''] :
		new_dict = data_dict
	else :
		# Initialize dict for each protein
		for protein in data_dict :
			new_dict[protein] = {}

		# For each feature type in FT_order_list, add data first
		for FT in FT_order_list :
			for protein in data_dict :
				for FT_ID, data in data_dict[protein].items() :
					if data_dict[protein][FT_ID][0] == FT :
						new_dict[protein][FT_ID] = data

		# Complete data for all feature type not in FT_order_list
		for protein in data_dict :
			for FT_ID, data in data_dict[protein].items() :
				if data[0] not in FT_order_list :
					new_dict[protein][FT_ID] = data

	return  new_dict

def extractLength(protein_length_df) :
	"""
		DESCRIPTION: Extracts the length for each protein
	"""

	# Initialize dict containing protein lengths
	prot_length_dict = {}

	# Add information of the DataFrame to a dict
	for index, row in protein_length_df.iterrows():
		prot_length_dict[row['protein']] = int(row['total_length'])

	return prot_length_dict

def extractProtCut(protein_cut_df) :
	"""
		DESCRIPTION: Transforms the protein cut DataFrame to a simple dict used for the map drawing process
	"""

	# Initialize dict containing protein cut regions
	prot_cut_dict = {}

	# Transform the DataFrame to a dict used in the map creation process
	for index, row in protein_cut_df.iterrows():
		if row['protein'] != '' and row['length'] != '' and row['start'] != ''	:
			prot_cut_dict[row['protein']] = (int(row['length']), int(row['start']))

	return prot_cut_dict

def getDefault(shape_df, height, draw_coverage, draw_2nd_struct, draw_disorder, draw_mod_res, draw_region) :
	"""
		DESCRIPTION: Function rassembling all default reprensation (shapes and colors) for all Feature parameters. These are used only if nothing is specifically defined by the user.
	"""

	# Initialize shape and colors DataFrame index from the actual size of the DataFrame
	shape_df_i = len(shape_df.index)

	# If 3D structure coverage selected
	if draw_coverage == 'True' :
		if 'PDB' not in list(shape_df['feature']) :
			if height != None :
				shape_df.loc[shape_df_i] = ['PDB', 'line', 'down', str(int(height * 1.5)), None, 'autumn_r', None, None, None, None, None]
			else :
				# Use default height if None is given
				shape_df.loc[shape_df_i] = ['PDB', 'line', 'down', str(int(20 * 1.5)), None, 'autumn_r', None, None, None, None, None]
			shape_df_i += 1

	# If secondary structure selected
	if draw_2nd_struct == 'True' :
		# Add Helix default if None is defined
		if 'Helix' not in list(shape_df['feature']) :
			shape_df.loc[shape_df_i] = ['Helix', 'rectangle', None, None, None, None, None, '#FF4040', None, None, None]
			shape_df_i += 1
		# Add Strand default if None is defined
		if 'Strand' not in list(shape_df['feature']) :
			shape_df.loc[shape_df_i] = ['Strand', 'arrow', 'right', None, None, None, None, 'skyblue', None, None, None]
			shape_df_i += 1
		# Add Turn default if None is defined
		if 'Turn' not in list(shape_df['feature']) :
			shape_df.loc[shape_df_i] = ['Turn', 'rhombus', None, None, None, None, None, 'gold', None, None, None]
			shape_df_i += 1

	# If disorder is selected :
	if draw_disorder == 'True':
		# Add default representation of disorder if None is defined
		if 'Disordered' not in list(shape_df['feature']) :
			shape_df.loc[shape_df_i] = ['Disordered', 'rectangle', None, '0', 'red', None, None, None, None, None, None]
			shape_df_i += 1

	# If modified residues is selected
	if draw_mod_res == 'True' :
		# Add default representation of modified residues if None is defined
		if 'MOD_RES' not in list(shape_df['feature']) :
			shape_df.loc[shape_df_i] = ['MOD_RES', 'mark', 'up', None, 'red', None, None, None, None, None, None]
			shape_df_i += 1
		# Add default representation of carbohydrates if None is defined
		if 'CARBOHYD' not in list(shape_df['feature']) :
			shape_df.loc[shape_df_i] = ['CARBOHYD', 'mark', 'up', None, 'blue', None, None, None, None, None, None]
			shape_df_i += 1
		# Add default representation of lipids if None is defined
		if 'LIPID' not in list(shape_df['feature']) :
			shape_df.loc[shape_df_i] = ['LIPID', 'mark', 'up', None, 'green', None, None, None, None, None, None]
			shape_df_i += 1

	# If compotion biased region is selected
	if draw_region != [] and draw_region != None and draw_region != 'None' :
		# Add default representation of composition biased region if None is defined
		if 'COMPBIAS' not in list(shape_df['feature']) :
			if height != None :
				shape_df.loc[shape_df_i] = ['COMPBIAS', 'l', 'down-left', str(int(height/4)), None, None, None, None, None, None, None]
			else :
				# Use default height if None is given
				shape_df.loc[shape_df_i] = ['COMPBIAS', 'l', 'down-left', str(int(20/4)), None, None, None, None, None, None, None]

			shape_df_i += 1

	return shape_df

def getShapeDict(shape_df, height, draw_coverage, draw_2nd_struct, draw_disorder, draw_mod_res, draw_region, uniform_contour=False):
	"""
		DESCRIPTION: Transforms the DataFrame containing all shapes and colors information into a dict used in the map drawing process.
	"""

	# Initialize shape and color dict
	shape_dict = {}

	# For each row in the shape and color DataFrame
	for index, row in shape_df.iterrows():
		# Initialize dict for the feature
		shape_dict[row['feature']] = {}
		# If a shape is defined add it to the dict in lowercase
		if row['shape'] != np.nan and row['shape'] != '' and row['shape'] != None :
			shape_dict[row['feature']]['shape'] = str(row['shape']).lower()
		# If an orientation is defined, add it to the dict in lowercase
		if row['orientation'] != np.nan and row['orientation'] != '' and row['orientation'] != None :
			shape_dict[row['feature']]['orientation'] = str(row['orientation']).lower()
		# If a height is defined, add the value to the dict
		if row['height'] != np.nan and row['height'] != '' and row['height'] != None :
			try :
				shape_dict[row['feature']]['height'] = int(float(row['height']))
			except : pass
		# If a contour color is defined, add it to the dict in lowercase
		if isinstance(row['contour_color'], str) and row['contour_color'] != '' : 
			shape_dict[row['feature']]['contour_color'] = str(row['contour_color']).lower()
		# If a contour colormap is defined, add it to the dict. No lowercase here because maptplotlib colormaps are case sensitive.
		if isinstance(row['contour_colormap'], str) and row['contour_colormap'] != '' :
			shape_dict[row['feature']]['contour_colormap'] = str(row['contour_colormap'])
		# If a threshold for the contour colormap is defined, cast it to float and add it to the dict
		if row['contour_threshold'] != np.nan and row['contour_threshold'] != '' and row['contour_threshold'] != None :
			try :
				shape_dict[row['feature']]['contour_threshold'] = float(row['contour_threshold'])
			except : pass
		# If a color is defined, add it to the dict in lowercase
		if row['color'] != np.nan and row['color'] != '' and row['color'] != None :
			shape_dict[row['feature']]['color'] = str(row['color']).lower()
		# If a colormap is defined, add it to the dict. No lowercase here because maptplotlib colormaps are case sensitive.
		if row['colormap'] != np.nan and row['colormap'] != '' and row['colormap'] != None :
			shape_dict[row['feature']]['colormap'] = str(row['colormap'])
		# If a threshold for the colormap is defined, cast it to float and add it to the dict
		if row['threshold'] != np.nan and row['threshold'] != '' and row['threshold'] != None :
			try :
				shape_dict[row['feature']]['threshold'] = float(row['threshold'])
			except : pass
		# If a pensize is defined, cast the value to an int and add it to the dict
		if row['pensize'] != np.nan and row['pensize'] != '' and row['pensize'] != None :
			try :
				shape_dict[row['feature']]['pensize'] = int(row['pensize'])
			except : pass

		# If the consistent color option is selected, change the contour color according to the color inside the shape.
		if uniform_contour and 'contour_color' not in shape_dict[row['feature']] and 'color' in shape_dict[row['feature']] :
			shape_dict[row['feature']]['contour_color'] = shape_dict[row['feature']]['color']

		# If the consistent color option is selected, change the contour colormap according to the colormap inside the shape.
		if uniform_contour and 'contour_colormap' not in shape_dict[row['feature']] and 'colormap' in shape_dict[row['feature']] :
			shape_dict[row['feature']]['contour_colormap'] = shape_dict[row['feature']]['colormap']

	return shape_dict

def drawContour(draw, vertex_list, contour_color, pensize, complete=True) :
	"""
		DESCRIPTION: Function that draw the contour of the shape. Borders are added to the input vertex_list to apply the joint='curve' option on the whole contour.
	"""

	# Make a copy of the vertex list
	outline_vertex_list = vertex_list.copy()
	if complete :
		# Add the first two elements of the vertex list to the end off the new vertex list if "complete" option is selected.
		outline_vertex_list.append(vertex_list[0])
		outline_vertex_list.append(vertex_list[1])
	# Draw a line according to the newly built vertex list
	draw.line(outline_vertex_list, fill=contour_color, width=pensize, joint='curve')

def drawForm(draw, x_pos, y_pos, size, shape='None', color='white', height=80, orientation='', contour_color='black', pensize=3, region=False) :
	"""
		DESCRIPTION: Draws the shape according to the shape ID and the additionnal orientation if needed.
			The function contains all shapes available to the user.
	"""

	# If a shape is defined
	if shape != 'None' :

		# Put the shape name in lowercase
		shape = shape.lower()

		# Horizontal gap
		if shape == 'gap' and orientation == 'horizontal' :
			# Vertex list containing all points to define the contour of the shape
			vertex_list = [
				(x_pos, y_pos - (height/2)),
				(x_pos + size, y_pos - (height/2)),
				(x_pos + size, y_pos + (height/2)),
				(x_pos, y_pos + (height/2))
			]
			# Fill the defined polygone with the fill color
			draw.polygon(vertex_list, fill=color)
			# Draw contour with the contour color
			drawContour(draw, [vertex_list[0], vertex_list[1]], contour_color, pensize, complete=False)
			drawContour(draw, [vertex_list[2], vertex_list[3]], contour_color, pensize, complete=False)

		# Vertical gap
		elif shape == 'gap' and orientation == 'vertical' :
			vertex_list = [
				(x_pos, y_pos - (height/2)),
				(x_pos + size, y_pos - (height/2)),
				(x_pos + size, y_pos + (height/2)),
				(x_pos, y_pos + (height/2))
			]
			draw.polygon(vertex_list, fill=color)
			drawContour(draw, [vertex_list[1], vertex_list[2]], contour_color, pensize, complete=False)
			drawContour(draw, [vertex_list[3], vertex_list[0]], contour_color, pensize, complete=False)

		# Lower side of rectangle linked on the left side. Creates a gap in the protein line.
		elif shape == 'l' and (orientation == 'down-left' or orientation == 'left-down') :
			vertex_list = [
				(x_pos, y_pos - (height/2)),
				(x_pos + size, y_pos - (height/2)),
				(x_pos + size, y_pos + (height/2)),
				(x_pos, y_pos + (height/2))
			]
			draw.polygon(vertex_list, fill=color)
			drawContour(draw, [(x_pos + size, y_pos + (height/2)), (x_pos, y_pos + (height/2)), (x_pos, y_pos)], contour_color, pensize, complete=False)

		# Upper side of rectangle linked on the left side. Creates a gap in the protein line.
		elif shape == 'l' and (orientation == 'up-left' or orientation == 'left-up') :
			vertex_list = [
				(x_pos, y_pos - (height/2)),
				(x_pos + size, y_pos - (height/2)),
				(x_pos + size, y_pos + (height/2)),
				(x_pos, y_pos + (height/2))
			]
			draw.polygon(vertex_list, fill=color)
			drawContour(draw, [(x_pos, y_pos), (x_pos, y_pos - (height/2)), (x_pos + size, y_pos - (height/2))], contour_color, pensize, complete=False)

		# Lower side of rectangle linked on the right side. Creates a gap in the protein line.
		elif shape == 'l' and (orientation == 'down-right' or orientation == 'right-down') :
			vertex_list = [
				(x_pos, y_pos - (height/2)),
				(x_pos + size, y_pos - (height/2)),
				(x_pos + size, y_pos + (height/2)),
				(x_pos, y_pos + (height/2))
			]
			draw.polygon(vertex_list, fill=color)
			drawContour(draw, [(x_pos, y_pos + (height/2)), (x_pos + size, y_pos + (height/2)), (x_pos + size, y_pos)], contour_color, pensize, complete=False)

		# Upper side of rectangle linked on the right side. Creates a gap in the protein line.
		elif shape == 'l' and (orientation == 'up-right' or orientation == 'right-up') :
			vertex_list = [
				(x_pos, y_pos - (height/2)),
				(x_pos + size, y_pos - (height/2)),
				(x_pos + size, y_pos + (height/2)),
				(x_pos, y_pos + (height/2))
			]
			draw.polygon(vertex_list, fill=color)
			drawContour(draw, [(x_pos, y_pos - (height/2)), (x_pos + size, y_pos - (height/2)), (x_pos + size, y_pos)], contour_color, pensize, complete=False)

		# Upper half of an isocele triangle pointing to the left
		elif shape == 'semitriangle' and (orientation == 'up-left' or orientation == 'left-up') :
			vertex_list = [
				(x_pos + size,  y_pos - (height/2)),
				(x_pos + size, y_pos),
				(x_pos, y_pos)
			]
			draw.polygon(vertex_list, fill=color, outline=contour_color)
			drawContour(draw, vertex_list, contour_color, pensize)

		# Lower half of an isocele triangle pointing to the left
		elif shape == 'semitriangle' and (orientation == 'down-left' or orientation == 'left-down') :
			vertex_list = [
				(x_pos + size,  y_pos + (height/2)),
				(x_pos + size, y_pos),
				(x_pos, y_pos)
			]
			draw.polygon(vertex_list, fill=color, outline=contour_color)
			drawContour(draw, vertex_list, contour_color, pensize)

		# Upper half of an isocele triangle pointing to the right
		elif shape == 'semitriangle' and (orientation == 'up-right' or orientation == 'right-up') :
			vertex_list = [
				(x_pos,  y_pos - (height/2)),
				(x_pos + size, y_pos),
				(x_pos, y_pos)
			]
			draw.polygon(vertex_list, fill=color, outline=contour_color)
			drawContour(draw, vertex_list, contour_color, pensize)

		# Lower half of an isocele triangle pointing to the right
		elif shape == 'semitriangle' and (orientation == 'down-right' or orientation == 'right-down') :
			vertex_list = [
				(x_pos,  y_pos + (height/2)),
				(x_pos + size, y_pos),
				(x_pos, y_pos)
			]
			draw.polygon(vertex_list, fill=color, outline=contour_color)
			drawContour(draw, vertex_list, contour_color, pensize)

		# Triangle pointing upwards which base is on the protein
		elif shape == 'semitriangle' and orientation == 'up' :
			vertex_list = [
				(x_pos + (size/2),  y_pos - (height/2)),
				(x_pos + size, y_pos),
				(x_pos, y_pos)
			]
			draw.polygon(vertex_list, fill=color, outline=contour_color)
			drawContour(draw, vertex_list, contour_color, pensize)

		# Triangle pointing downwards which base is on the protein
		elif shape == 'semitriangle' and orientation == 'down' :
			vertex_list = [
				(x_pos + (size/2),  y_pos + (height/2)),
				(x_pos + size, y_pos),
				(x_pos, y_pos)
			]
			draw.polygon(vertex_list, outline=contour_color)
			drawContour(draw, vertex_list, contour_color, pensize)

		# Rectangle on the upper side of the protein, with height halved compared to standard rectangle
		elif shape == 'semirectangle' and orientation == 'up' :
			if region == True :
				color = 'white'
			vertex_list = [
				(x_pos + size, y_pos),
				(x_pos + size, y_pos - (height/2)),
				(x_pos, y_pos - (height/2)),
				(x_pos, y_pos)
			]
			if region == True :
				draw.polygon(vertex_list, fill=color)
				drawContour(draw, vertex_list, contour_color, pensize, complete=False)
			else :
				draw.polygon(vertex_list, fill=color, outline=contour_color)
				drawContour(draw, vertex_list, contour_color, pensize)

		# Rectangle on the lower side of the protein, with height halved compared to standard rectangle
		elif shape == 'semirectangle' and orientation == 'down' :
			if region == True :
				color = 'white'
			vertex_list = [
				(x_pos + size, y_pos),
				(x_pos + size, y_pos + (height/2)),
				(x_pos, y_pos + (height/2)),
				(x_pos, y_pos)
			]
			if region == True :
				draw.polygon(vertex_list, fill=color)
				drawContour(draw, vertex_list, contour_color, pensize, complete=False)
			else :
				draw.polygon(vertex_list, fill=color, outline=contour_color)
				drawContour(draw, vertex_list, contour_color, pensize)

		# Parallel line above the protein line
		elif shape == 'line' and orientation == 'up' :
			vertex_list = [
				(x_pos, y_pos - (height/2)),
				(x_pos + size, y_pos - (height/2))
			]
			draw.polygon(vertex_list, outline=contour_color)
			drawContour(draw, vertex_list, contour_color, pensize)

		# Parallel line below the protein line
		elif shape == 'line' and orientation == 'down' :
			vertex_list = [
				(x_pos, y_pos + (height/2)),
				(x_pos + size, y_pos + (height/2))
			]
			draw.polygon(vertex_list, outline=contour_color)
			drawContour(draw, vertex_list, contour_color, pensize)

		# Hexagon where 2 of the sides are perpendicular to the protein line
		elif shape == 'hexagon' and orientation == 'vertical' :
			vertex_list = [
				(x_pos, y_pos - (height/4)),
				(x_pos + (size/2), y_pos - (height/2)),
				(x_pos + size, y_pos - (height/4)),
				(x_pos + size, y_pos + (height/4)),
				(x_pos + (size/2), y_pos + (height/2)),
				(x_pos, y_pos + (height/4)),
				(x_pos, y_pos)
			]
			draw.polygon(vertex_list, fill=color, outline=contour_color)
			drawContour(draw, vertex_list, contour_color, pensize)

		# Upper half of hexagon where 2 of the sides are perpendicular to the protein line.
		elif shape == 'semihexagon' and (orientation == 'vertical-up' or orientation == 'up-vertical') :
			vertex_list = [
				(x_pos, y_pos - (height/4)),
				(x_pos + (size/2), y_pos - (height/2)),
				(x_pos + size, y_pos - (height/4)),
				(x_pos + size, y_pos),
				(x_pos, y_pos)
			]
			draw.polygon(vertex_list, fill=color, outline=contour_color)
			drawContour(draw, vertex_list, contour_color, pensize)

		# Lower half of hexagon where 2 of the sides are perpendicular to the protein line.
		elif shape == 'semihexagon' and (orientation == 'vertical-down' or orientation == 'down-vertical') :
			vertex_list = [
				(x_pos, y_pos + (height/4)),
				(x_pos + (size/2), y_pos + (height/2)),
				(x_pos + size, y_pos + (height/4)),
				(x_pos + size, y_pos),
				(x_pos, y_pos)
			]
			draw.polygon(vertex_list, fill=color, outline=contour_color)
			drawContour(draw, vertex_list, contour_color, pensize)

		# Hexagon where 2 of the sides are parallel to the protein line
		elif shape == 'hexagon' and orientation == 'horizontal' :
			vertex_list = [
				(x_pos + (size/4), y_pos - (height/2)),
				(x_pos + (3*size/4), y_pos - (height/2)),
				(x_pos + size, y_pos),
				(x_pos + (3*size/4), y_pos + (height/2)),
				(x_pos + (size/4), y_pos + (height/2)),
				(x_pos, y_pos)
			]
			draw.polygon(vertex_list, fill=color, outline=contour_color)
			drawContour(draw, vertex_list, contour_color, pensize)

		# Upper half of hexagon where 2 of the sides are parallel to the protein line.
		elif shape == 'semihexagon' and (orientation == 'horizontal-up' or orientation == 'up-horizontal') :
			vertex_list = [
				(x_pos + (size/4), y_pos - (height/2)),
				(x_pos + (3*size/4), y_pos - (height/2)),
				(x_pos + size, y_pos),
				(x_pos, y_pos)
			]
			draw.polygon(vertex_list, fill=color, outline=contour_color)
			drawContour(draw, vertex_list, contour_color, pensize)

		# Lower half of hexagon where 2 of the sides are parallel to the protein line.
		elif shape == 'semihexagon' and (orientation == 'horizontal-down' or orientation == 'down-horizontal') :
			vertex_list = [
				(x_pos + (size/4), y_pos + (height/2)),
				(x_pos + (3*size/4), y_pos + (height/2)),
				(x_pos + size, y_pos),
				(x_pos, y_pos)
			]
			draw.polygon(vertex_list, fill=color, outline=contour_color)
			drawContour(draw, vertex_list, contour_color, pensize)

		# Rhombus
		elif shape == 'rhombus' :
			vertex_list = [
				(x_pos + (size/2), y_pos - (height/2)),
				(x_pos + size, y_pos),
				(x_pos + (size/2), y_pos + (height/2)),
				(x_pos, y_pos)
			]
			draw.polygon(vertex_list, fill=color, outline=contour_color)
			drawContour(draw, vertex_list, contour_color, pensize)

		# Upper half of rhombus. (same as triangle up)
		elif shape == 'semirhombus' and orientation == 'up' :
			vertex_list = [
				(x_pos + (size/2), y_pos - (height/2)),
				(x_pos + size, y_pos),
				(x_pos, y_pos)
			]
			draw.polygon(vertex_list, fill=color, outline=contour_color)
			drawContour(draw, vertex_list, contour_color, pensize)

		# Lower half of rhombus. (same as triangle down)
		elif shape == 'semirhombus' and orientation == 'down' :
			vertex_list = [
				(x_pos + (size/2), y_pos + (height/2)),
				(x_pos + size, y_pos),
				(x_pos, y_pos)
			]
			draw.polygon(vertex_list, fill=color, outline=contour_color)
			drawContour(draw, vertex_list, contour_color, pensize)

		# Rectangle
		elif shape == 'rectangle' :
			vertex_list = [
				(x_pos, y_pos - (height/2)),
				(x_pos + size, y_pos - (height/2)),
				(x_pos + size, y_pos + (height/2)),
				(x_pos, y_pos + (height/2))
			]
			draw.polygon(vertex_list, fill=color, outline=contour_color)
			drawContour(draw, vertex_list, contour_color, pensize)

		# Arrow pointing to the right
		elif shape == 'arrow' and orientation == 'right' :
			vertex_list = [
				(x_pos, y_pos - (height/2)),
				(x_pos + (size/2), y_pos - (height/2)),
				(x_pos + size, y_pos),
				(x_pos + (size/2), y_pos + (height/2)),
				(x_pos, y_pos + (height/2)),
				(x_pos, y_pos)
			]
			draw.polygon(vertex_list, fill=color, outline=contour_color)
			drawContour(draw, vertex_list, contour_color, pensize)

		# Upper half of arrow pointing to the right
		elif shape == 'semiarrow' and (orientation == 'right-up' or orientation == 'up-right') :
			vertex_list = [
				(x_pos, y_pos - (height/2)),
				(x_pos + (size/2), y_pos - (height/2)),
				(x_pos + size, y_pos),
				(x_pos, y_pos)
			]
			draw.polygon(vertex_list, fill=color, outline=contour_color)
			drawContour(draw, vertex_list, contour_color, pensize)

		# Lower half of arrow pointing to the right
		elif shape == 'semiarrow' and (orientation == 'right-down' or orientation == 'down-right') :
			vertex_list = [
				(x_pos, y_pos + (height/2)),
				(x_pos + (size/2), y_pos + (height/2)),
				(x_pos + size, y_pos),
				(x_pos, y_pos)
			]
			draw.polygon(vertex_list, fill=color, outline=contour_color)
			drawContour(draw, vertex_list, contour_color, pensize)

		# Arrow pointing to the left
		elif shape == 'arrow' and orientation == 'left' :
			vertex_list = [
				(x_pos + (size/2), y_pos - (height/2)),
				(x_pos + size, y_pos - (height/2)),
				(x_pos + size, y_pos + (height/2)),
				(x_pos + (size/2), y_pos + (height/2)),
				(x_pos, y_pos)
			]
			draw.polygon(vertex_list, fill=color, outline=contour_color)
			drawContour(draw, vertex_list, contour_color, pensize)

		# Upper half of arrow pointing to the left
		elif shape == 'semiarrow' and (orientation == 'left-up' or orientation == 'up-left') :
			vertex_list = [
				(x_pos + (size/2), y_pos - (height/2)),
				(x_pos + size, y_pos - (height/2)),
				(x_pos + size, y_pos),
				(x_pos, y_pos)
			]
			draw.polygon(vertex_list, fill=color, outline=contour_color)
			drawContour(draw, vertex_list, contour_color, pensize)

		# Lower half of arrow pointing to the left
		elif shape == 'semiarrow' and (orientation == 'left-down' or orientation == 'down-left') :
			vertex_list = [
				(x_pos + (size/2), y_pos + (height/2)),
				(x_pos + size, y_pos + (height/2)),
				(x_pos + size, y_pos),
				(x_pos, y_pos)
			]
			draw.polygon(vertex_list, fill=color, outline=contour_color)
			drawContour(draw, vertex_list, contour_color, pensize)

		# Triangle pointing to the left
		elif shape == 'triangle' and orientation == 'left' :
			vertex_list = [
				(x_pos + size, y_pos - (height/2)),
				(x_pos + size, y_pos + (height/2)),
				(x_pos, y_pos)
			]
			draw.polygon(vertex_list, fill=color, outline=contour_color)
			drawContour(draw, vertex_list, contour_color, pensize)

		# Triangle pointing to the right
		elif shape == 'triangle' and orientation == 'right' :
			vertex_list = [
				(x_pos, y_pos - (height/2)),
				(x_pos + size, y_pos),
				(x_pos, y_pos + (height/2))
			]
			draw.polygon(vertex_list, fill=color, outline=contour_color)
			drawContour(draw, vertex_list, contour_color, pensize)

		# Small line detached and perpendicular to the protein at the beginning of the feature limit. Pointing upwards
		elif shape == 'mark' and orientation == 'up' :
			vertex_list = [
				(x_pos, y_pos - (height/2)),
				(x_pos, y_pos - (3*height/4)),
			]
			draw.polygon(vertex_list, outline=contour_color)
			drawContour(draw, vertex_list, contour_color, pensize)

		# Small line detached and perpendicular to the protein at the beginning of the feature limit. Pointing downwards
		elif shape == 'mark' and orientation == 'down' :
			vertex_list = [
				(x_pos, y_pos + (height/2)),
				(x_pos, y_pos + (3*height/4)),
			]
			draw.polygon(vertex_list, outline=contour_color)
			drawContour(draw, vertex_list, contour_color, pensize)

		# Bar perpendicular to the protein at the beginning of the feature limit.
		elif shape == 'bar' :
			vertex_list = [
				(x_pos, y_pos - (height/2)),
				(x_pos, y_pos + (height/2)),
			]
			draw.polygon(vertex_list, outline=contour_color)
			drawContour(draw, vertex_list, contour_color, pensize)

		else :
			print('Invalid shape and/or orientation!')
			print(shape, orientation)

def drawDiscontinuity(draw, x_pos, y_pos, size=20, height=50, pensize=0) :
	"""
		DESCRIPTION: Draws a discontinuity (-//-) of the line representing the protein.
	"""
	# Recolor the fraction of the base line in white
	draw.line([(x_pos+size*0.25, y_pos), (x_pos+size*0.75, y_pos)], fill='white', width=pensize)

	# Draw the discontinuity
	draw.line([(x_pos, y_pos + (height/2)), (x_pos + (size/2), y_pos - (height/2))], fill='black', width=pensize)
	draw.line([(x_pos + (size/2), y_pos + (height/2)), (x_pos + size, y_pos - (height/2))], fill='black', width=pensize)

def takeSecond(elem):
	"""
		DESCRIPTION: Return the 2nd element of the object
	"""
	return elem[1]

def proteinListSorting(protein_list, data_dict, protein_length_dict=None, sorting=None, focus=None, threshold=0.20, selection=False) :
	"""
		DESCRIPTION: Sorts the protein list
			'abc': Alphabetical sorting 
			'value': Sorts proteins in descending order of values on a given feature (focus). If multiple values are in a protein, the highest will be used as protein representative.
			'feature_number_distance': Sorts the proteins using the number of common features, the number of features not in common and protein length. It will aggregate protein will similar protein contents.
		WARNING: 'feature_number_distance' sorting option is of n² complexity. Protein list of more than 1k protein will make the sorting lasting more than 5-10mins.
	"""

	# Not sorted protein list
	if sorting == None :
		pass
	# Alphabetically sorted protein list
	elif sorting == 'abc' :
		protein_list = sorted(protein_list)
   	# Sorting by numerical value
	elif sorting == 'value' and focus != None :
		new_prot_list = []
		for protein in protein_list :
			# Search for the highest value of all the focused features
			BI_max = -1
			try :
				BI_max = max([float(data_dict[protein][j][4]) for j in data_dict[protein] if data_dict[protein][j][1] == focus])
			except :
				BI_max = -1
			# Select only the protein with a max value above the given threshold
			if BI_max >= threshold :
				new_prot_list.append((protein, BI_max))
		# Sort proteins by descending max values
		new_prot_list.sort(key = takeSecond, reverse = True)
		protein_list = [elem[0] for elem in new_prot_list]

	# Sorting by feature content
	elif sorting == 'feature_number_distance' :
		# Get all features present in protein list
		feature_list = []
		for protein, protein_data in data_dict.items():
			for feature_id, feature_data in protein_data.items():
				FT_type = feature_data[0]
				feature = feature_data[1]
				if FT_type in ['DOMAIN', 'REPEAT', 'MOTIF', 'REGION'] :
					feature_list.append(feature)

		# Get unique feature list and sort
		feature_list = list(set(feature_list))
		feature_list.sort()

		# Initialise dictionnary for each protein
		feature_number_dict = {}
		for protein, protein_data in data_dict.items():
			feature_number_dict[protein] = {}
			for feature in feature_list :
				feature_number_dict[protein][feature] = 0

		# Parse data_dict and count each number of each feature
		for protein, protein_data in data_dict.items():
			for feature_id, feature_data in protein_data.items():
				FT_type = feature_data[0]
				feature = feature_data[1]
				if FT_type in ['DOMAIN', 'REPEAT', 'MOTIF', 'REGION'] :
					feature_number_dict[protein][feature] += 1

		# Create DataFrame from dict containing number of each feature
		df = pd.DataFrame.from_dict(feature_number_dict, orient='index')


		# Create similarity and dissimilarity dataframe
		disim_dist_dict = {}
		sim_dist_dict = {}
		# For each protein in the list
		for protein in protein_list :
			# Initialize similarity and dissimilary dict for each protein
			disim_dist_dict[protein] = {}
			sim_dist_dict[protein] = {}
			# Get data for reference protein
			protein_data = list(df.loc[protein])
			# For each protein to compare with
			for vs_protein in protein_list :
				# Get data for protein to compare with
				vs_protein_data = list(df.loc[vs_protein])
				# Initialize similarity and dissimilarity dict for compared protein to 0
				disim_dist_dict[protein][vs_protein] = 0
				sim_dist_dict[protein][vs_protein] = 0
				# For each feature found in the list of protein
				for item in zip(protein_data, vs_protein_data) :
					# Get the number of the given feature
					val_protein, val_vs_protein = item
					# Add to the dissimilarity score the absolute value of the difference in number of the given feature
					disim_dist_dict[protein][vs_protein] += abs(val_protein - val_vs_protein)
					# If there is more than one feature for the two proteins
					if val_protein > 0 and val_vs_protein > 0 :
						# Add the number of common features to the similiarity score
						sim_dist_dict[protein][vs_protein] += min(val_protein, val_vs_protein)

		# Create DataFrames from the similarity and dissimilarity dict
		disim_dist_df = pd.DataFrame.from_dict(disim_dist_dict)
		sim_dist_df = pd.DataFrame.from_dict(sim_dist_dict)

		# Setting diagonal to zeros
		disim_dist_df.values[[np.arange(disim_dist_df.shape[0])]*2] = -1
		sim_dist_df.values[[np.arange(sim_dist_df.shape[0])]*2] = -1

		# Get max similarity and dissimilarity scores
		max_sim = max(sim_dist_df.max())
		max_disim = max(disim_dist_df.max())

		# Create all combination of pairs of proteins
		protein_comb = list(set(list(itertools.combinations(protein_list, 2))))
		protein_comb.sort()

		## Creating a list of priorities to link the proteins based on the content
		# In order of priority : Similarity (descending) > Dissimilarity (ascending) > Sequence size difference (ascending)

		# Initializing step number and list of priorities
		step = 0
		floor_list = []

		# From highest to lowest similarity score
		for sim in [i for i in range(0,max_sim+1)][::-1] :
			# Initialize list of links with same similarity score
			disim_list = []
			# For each combination of pairs of proteins
			for comb in protein_comb :
				# If similarity score corresponds to the searched similarity score
				if sim_dist_df[comb[0]][comb[1]] == sim :
					# Add combination to the list (protein1, protein2, dissimilarity score, absolute difference of sequence length)
					disim_list.append((comb[0], comb[1], disim_dist_df[comb[0]][comb[1]], abs(protein_length_dict[comb[0]] - protein_length_dict[comb[1]])))

			# Sort the list of links based on the difference of sequence length (ascending)
			disim_list = sorted(disim_list, key=lambda tup: tup[3])
			# Sort the list of links based on the dissimilarity score (ascending)
			disim_list = sorted(disim_list, key=lambda tup: tup[2])

			# For each link in the link list
			for comb in disim_list :
				# Append it to the global priority list
				floor_list.append((step, (comb[0], comb[1])))
				step += 1

		# Initialize dict containing the number of links existing for each protein
		link_dict = {}
		for protein in protein_list :
			link_dict[protein] = 0

		## Linking algorithm
		# The objective of this part is to create a list of protein, where proteins are aggregated based on their feature content. 

		serie_dict = {}
		# For each possible link from highest to lowest priority
		for floor in floor_list :
			# Get the 2 proteins inplicated in the potential link
			protein_1 = floor[1][0]
			protein_2 = floor[1][1]
			# Check if proteins are already linked or not
			out_1 = isUsed(serie_dict, protein_1)
			out_2 = isUsed(serie_dict, protein_2)
			# In case the two protein have existing links :
			if out_1 != False and out_2 != False :
				# If the two proteins have only one link
				if link_dict[protein_1] == 1 and link_dict[protein_2] == 1 :
					# Get to position in priority list of the 2 proteins
					floor_id_1, elem_nb_1 = out_1
					floor_id_2, elem_nb_2 = out_2
					if floor_id_1 != floor_id_2 :
						# Get minimum and maximum postion of the 2 proteins in the priority list
						min_floor = min(floor_id_1, floor_id_2)
						max_floor = max(floor_id_1, floor_id_2)
						# Add a link to each protein
						link_dict[protein_1] += 1
						link_dict[protein_2] += 1
						# Create link between the two existing lists depending on the orientation of the lists
						if elem_nb_1 == 0 and elem_nb_2 == 0 :
							# Flipping list containing protein 1 and append list containing protein 2
							serie_dict[min_floor] = serie_dict[floor_id_1][::-1] + serie_dict[floor_id_2]
						elif elem_nb_1 == 0 and elem_nb_2 != 0 :
							# Append to list containing protein 2 list containing protein 1
							serie_dict[min_floor] = serie_dict[floor_id_2] +  serie_dict[floor_id_1]
						elif elem_nb_1 != 0 and elem_nb_2 == 0 :
							# Append to list containing protein 1 list containing protein 2
							serie_dict[min_floor] = serie_dict[floor_id_1] + serie_dict[floor_id_2]
						elif elem_nb_1 != 0 and elem_nb_2 != 0 :
							# Add to list containging protein 1 flipped list containing protein 2
							serie_dict[min_floor] = serie_dict[floor_id_1] + serie_dict[floor_id_2][::-1]
						del serie_dict[max_floor]

			# In case only protein 2 has already a link 
			elif out_1 != False :
				# Get to position in priority list of the not linked protein
				floor_id, elem_nb = out_1
				floor_protein_list = serie_dict[floor_id]
				# Protein to link at the start of the existing linked list
				if elem_nb == 0 :
					# Create link between isolated protein and linked list of proteins by the protein 1 end
					serie_dict[floor_id] = [protein_2] + floor_protein_list
					# Add a link to each protein
					link_dict[protein_1] += 1
					link_dict[protein_2] += 1
				# Protein to link at the end of the existing linked list
				elif elem_nb == len(floor_protein_list)-1 :
					# Create link between isolated protein and linked list of proteins by the protein 1 end
					serie_dict[floor_id] = floor_protein_list + [protein_2]
					# Add a link to each protein
					link_dict[protein_1] += 1
					link_dict[protein_2] += 1

			# In case only protein 1 has already one link
			elif out_2 != False :
				# Get to position in priority list of the not linked protein
				floor_id, elem_nb = out_2
				floor_protein_list = serie_dict[floor_id]
				# Protein to link at the start of the existing linked list
				if elem_nb == 0 :
					# Create link between isolated protein and linked list of proteins by the protein 2 end
					serie_dict[floor_id] = [protein_1] + floor_protein_list
					# Add a link to each protein
					link_dict[protein_1] += 1
					link_dict[protein_2] += 1
				# Protein to link at the end of the existing linked list
				elif elem_nb == len(floor_protein_list)-1 :
					# Create link between isolated protein and linked list of proteins by the protein 2 end
					serie_dict[floor_id] = floor_protein_list + [protein_1]
					# Add a link to each protein
					link_dict[protein_1] += 1
					link_dict[protein_2] += 1

			# In case none of the two protein have exisiting links
			else :
				# Link the two proteins
				serie_dict[floor[0]] = [protein_1, protein_2]
				# Add a link to each protein
				link_dict[protein_1] += 1
				link_dict[protein_2] += 1

		protein_list = serie_dict[0]

		## Reorder the feature count DataFrame in case you want to output it
		#df = df.reindex(protein_list)

	return protein_list

def isUsed(serie_dict, protein) :
	"""
		DESCRIPTION: Get the position of the priority list if protein is used and the position in the list
	"""
	for floor_id, floor in serie_dict.items() :
		for i in range(len(floor)) :
			if floor[i] == protein :
				return (floor_id, i)
	return False

def drawFigure(protein_list_df, data_df, shape_df, protein_cut_df, protein_length_df, pdb_coverage_list, case, length_factor=0.2, height=800, text_size=5, sorting=None, focus=None, threshold=0.2, biased_region_text_size=4, draw_region=[], draw_coverage='False', draw_2nd_struct='False', draw_disorder='False', draw_mod_res='False', draw_protein_length='False', uniform=False, pensize=0, FT_order_list=['']) :
	"""
		DESCRIPTION: Function handling all the map drawing process
	"""

	# Dict where keys are category of composition biased regions and values the associated letter
	comp_bias_dict = {'Acidic residues' : 'A', 'Basic residues' : 'B', 'Basic and acidic residues' : 'C', 'Polar residues' : 'H', 'Pro residues' : 'P'}

	# Get the protein list
	base_protein_list = list(protein_list_df['protein'])

	print('drawFigure')
	print('case', case)

	# Get data as dict
	data_dict = extractData(base_protein_list, data_df, FT_order_list, case=case)

	# Get dict with all protein lengths
	length_dict = extractLength(protein_length_df)

	# Sort protein list
	protein_list = proteinListSorting(base_protein_list, data_dict, protein_length_dict=length_dict, sorting=sorting, focus=focus, threshold = threshold)
	sorted_protein_list = protein_list.copy()

	# Get number of proteins
	max_length = len(protein_list)
	print('max length', max_length)

	#Chooses a text_size depending on the number of displayed proteins, if None is given.
	if text_size == None :
		text_size = int(-0.1 * max_length) + 18

	# Get cut regions if any as dict
	prot_cut_dict = extractProtCut(protein_cut_df)

	# If no pensize is defined, put the default value
	if pensize == 0 :
		pensize = 3

	# Reorder columns of saved data
	shape_df = shape_df[['feature', 'shape', 'orientation', 'height', 'contour_color', 'contour_colormap', 'contour_threshold', 'color', 'colormap', 'threshold', 'pensize']]
	# Replace all empty cells by None
	shape_df = shape_df.replace({np.nan: None})

	# Get all default shapes and color parameters for selected parameters
	shape_df = getDefault(shape_df, height, draw_coverage, draw_2nd_struct, draw_disorder, draw_mod_res, draw_region)

	# Transform shape DataFrame to dict
	shape_dict = getShapeDict(shape_df, height, draw_coverage, draw_2nd_struct, draw_disorder, draw_mod_res, draw_region, uniform_contour=uniform)
	feature_shape_dict = copy.deepcopy(shape_dict)

	# Define space between 2 proteins based on the height value
	space_between_protein = height*2

	# Create empty image with the adequate height based on the number of proteins
	im = Image.new('RGB', (4000, space_between_protein * max_length), 'white')
	# Initialize the draw Object used to draw on the image
	draw = ImageDraw.Draw(im)
	# Load the text font
	text_font = ImageFont.truetype('arial-bold.ttf', text_size)

	print('start of createFigure')

	# Check if the protein list is not empty
	if protein_list != [] :

		# Get the maximum text width of the protein names. Used to align names to the right side.
		max_text_w = max([draw.textsize(protein, font=text_font)[0] for protein in protein_list if protein != None])

		# For each protein in the list.
		for i in range(max_length) :
			# Get name of the  protein if there is one
			protein = None
			try :
				protein  = protein_list[i]
			except :
				protein = None

			# If there is a protein
			if protein != None :

				# Get position of the protein name zone
				plot_x_pos = 20
				plot_y_pos = height + space_between_protein * i

				# Get height and width of protein name
				w, h = draw.textsize(protein, font=text_font)
				# Draw text
				draw.text((10+max_text_w-w+plot_x_pos, plot_y_pos-(text_size/3)), protein, (0,0,0), font=text_font)

				# Get position of the protein
				plot_x_pos = 50+max_text_w+plot_x_pos
				plot_y_pos = height + space_between_protein * i

				# Get cut region postion and size if defined
				remove = 0
				if protein in prot_cut_dict :
					remove = prot_cut_dict[protein][0]

				# Draw a line for the protein corresponding to its length
				draw.line([(plot_x_pos, plot_y_pos), (plot_x_pos + (length_dict[protein] - remove) * length_factor, plot_y_pos)], fill='black', width=pensize, joint='curve')
				
				# If draw_protein_length option is activated, draw protein size at the end of the protein
				if draw_protein_length == 'True' :
					draw.text(((50 + plot_x_pos + (length_dict[protein] - remove) * length_factor, plot_y_pos-(text_size/3))), str(length_dict[protein]), (0,0,0), font=text_font)

				# If the protein cut regions are defined, draw discontinuity at the corresponding place 
				if protein in prot_cut_dict :
					drawDiscontinuity(draw, plot_x_pos+prot_cut_dict[protein][1] * length_factor, plot_y_pos, height * length_factor, height = height, pensize = pensize)

				# For each feature in the protein
				for j in data_dict[protein] :

					# Check if position of feature is a valid number (float or int)
					if data_dict[protein][j][2].replace('.0', '').isdigit():

						data_dict[protein][j][2] = int(float(data_dict[protein][j][2]))

						# Switch indicating if the feature should be drawn or not
						do_draw = True
						# Variable containing the size of the cut region
						cut = 0
						# If the protein contains a cut region
						if protein in prot_cut_dict :
							# If feature position is after cut region, get the value to remove to shift position of feature
							if int(data_dict[protein][j][2]) > (prot_cut_dict[protein][1] + prot_cut_dict[protein][0]) :
								cut =  prot_cut_dict[protein][0]
							# If the feature position is in a cut region, do not draw it
							if int(data_dict[protein][j][2]) > prot_cut_dict[protein][1] and int(data_dict[protein][j][2]) < (prot_cut_dict[protein][1] + prot_cut_dict[protein][0]) :
								do_draw = False

						# If the feature should be drawn
						if do_draw == True :
							# print(do_draw)

							# Get the type of feature
							FT_type = data_dict[protein][j][0]
							# Get the feature name
							feature = data_dict[protein][j][1]

							# If the feature is a DOMAIN or REPEAT
							if FT_type in ['DOMAIN', 'REPEAT'] :

								# Get feature information
								feature_data = data_dict[protein][j]


								# For each feature, draw the corresponding shape
								if feature in shape_dict :

									# Get the defined shape information
									dom_dict = feature_shape_dict[feature_data[1]]

									# Add the rescaled size for the feature
									dom_dict['size'] = int(float(feature_data[3])) * length_factor

									# If no height is defined, use the general height value
									if 'height' not in dom_dict :
										dom_dict['height'] = height

									# If no pensize is defined, use the protein thickness 
									if 'pensize' not in dom_dict :
										dom_dict['pensize'] = pensize

									# Get numerical value if any
									BI = None
									if feature_data[4] != 'None' and feature_data[4] != '' :
										BI = float(feature_data[4])

									# Draw the shape corresponding to the feature
									drawShape(draw, plot_x_pos+(int(feature_data[2]) - cut) * length_factor, plot_y_pos, dom_dict, BI)

								# If Others category is defined in the shape and color table
								elif 'Others' in shape_dict :
									# And if feature not in the rest of the shape and color table
									if feature not in shape_dict :
										feature_data = data_dict[protein][j]

										dom_dict = feature_shape_dict['Others']

										# Add the rescaled size for the feature
										dom_dict['size'] = int(float(feature_data[3])) * length_factor

										if 'height' not in dom_dict :
											dom_dict['height'] = height

										# If no pensize is defined, use the protein thickness 
										if 'pensize' not in dom_dict :
											dom_dict['pensize'] = pensize

										# If no color is defined, take white.
										if 'color' not in dom_dict or dom_dict['color'] == 'nan' :
											dom_dict['color'] = 'white'

										BI = None

										drawShape(draw, plot_x_pos+(int(feature_data[2]) - cut) * length_factor, plot_y_pos, dom_dict, BI)


							#Filter out the compositon biased regions (Rich or Poly aa)
							elif FT_type == 'COMPBIAS' :

								if FT_type in shape_dict :
									dom_dict = feature_shape_dict['COMPBIAS']
									if 'height' not in dom_dict :
										dom_dict['height'] = height
									if 'pensize' not in dom_dict :
										dom_dict['pensize'] = pensize

									feature_data = data_dict[protein][j]
									dom_dict['size'] = int(float(feature_data[3])) * length_factor

									if 'A' in draw_region and feature == 'Acidic residues' :
										drawShape(draw, plot_x_pos+(int(feature_data[2]) - cut) * length_factor, plot_y_pos, dom_dict, biased_region_text='A', biased_region_text_size=biased_region_text_size)
									elif 'B' in draw_region and feature == 'Basic residues' :
										drawShape(draw, plot_x_pos+(int(feature_data[2]) - cut) * length_factor, plot_y_pos, dom_dict, biased_region_text='B', biased_region_text_size=biased_region_text_size)
									elif 'C' in draw_region and feature == 'Basic and acidic residues' :
										drawShape(draw, plot_x_pos+(int(feature_data[2]) - cut) * length_factor, plot_y_pos, dom_dict, biased_region_text='C', biased_region_text_size=biased_region_text_size)
									elif 'H' in draw_region and feature == 'Polar residues' :
										drawShape(draw, plot_x_pos+(int(feature_data[2]) - cut) * length_factor, plot_y_pos, dom_dict, biased_region_text='H', biased_region_text_size=biased_region_text_size)
									elif 'P' in draw_region and feature == 'Pro residues' :
										drawShape(draw, plot_x_pos+(int(feature_data[2]) - cut) * length_factor, plot_y_pos, dom_dict, biased_region_text='P', biased_region_text_size=biased_region_text_size)

							elif FT_type == 'VARIANT' :
								if feature in shape_dict :

									feature_data = data_dict[protein][j]

									dom_dict = feature_shape_dict[feature_data[1]]

									dom_dict['size'] = 1 * length_factor

									if 'height' not in dom_dict :
										dom_dict['height'] = height

									if 'pensize' not in dom_dict :
										dom_dict['pensize'] = pensize

									BI = None

									if feature_data[4] != 'None' and feature_data[4] != '' :
										BI = float(feature_data[4])

									drawShape(draw, plot_x_pos+(int(feature_data[2]) - cut) * length_factor, plot_y_pos, dom_dict, BI)

							elif FT_type in ['MOD_RES', 'LIPID', 'CARBOHYD'] :

								feature_data = data_dict[protein][j]
								must_draw = False

								# Check if feature is indicated in feature shape
								if feature in shape_dict :
									dom_dict = feature_shape_dict[feature_data[1]]

									if 'height' not in dom_dict :
										dom_dict['height'] = height

									if 'pensize' not in dom_dict :
										dom_dict['pensize'] = pensize

									must_draw = True

								# Check if general feature type is indicated instead
								elif FT_type in shape_dict :
									dom_dict = feature_shape_dict[FT_type]

									if 'height' not in dom_dict :
										dom_dict['height'] = height

									if 'pensize' not in dom_dict :
										dom_dict['pensize'] = pensize

									must_draw = True


								# Use default, only if modified residues toggle is activated
								elif draw_mod_res == 'True' :
									#Default
									dom_dict = {}
									dom_dict['pensize'] = pensize

									# if FT_type == 'MOD_RES' :
									# 	dom_dict['shape'] = 'mark'
									# 	dom_dict['orientation'] = 'up'
									# 	dom_dict['contour_color'] = 'red'
									# 	dom_dict['color'] = 'red'
									# 	dom_dict['height'] = height
									# elif FT_type == 'LIPID' :
									# 	dom_dict['shape'] = 'mark'
									# 	dom_dict['orientation'] = 'up'
									# 	dom_dict['contour_color'] = 'green'
									# 	dom_dict['color'] = 'green'
									# 	dom_dict['height'] = height
									# elif FT_type == 'CARBOHYD' :
									# 	dom_dict['shape'] = 'mark'
									# 	dom_dict['orientation'] = 'up'
									# 	dom_dict['contour_color'] = 'blue'
									# 	dom_dict['color'] = 'blue'
									# 	dom_dict['height'] = height

									must_draw = True

								if must_draw :

									dom_dict['size'] = 1 * length_factor

									BI = None

									if feature_data[4] != 'None' and feature_data[4] != '' :
										BI = float(feature_data[4])

									drawShape(draw, plot_x_pos+(int(feature_data[2]) - cut) * length_factor, plot_y_pos, dom_dict, BI)

							elif FT_type == 'MOTIF' :
								feature_data = data_dict[protein][j]
								if feature in shape_dict :

									dom_dict = feature_shape_dict[feature_data[1]]

									dom_dict['size'] = 1 * length_factor

									if 'height' not in dom_dict :
										dom_dict['height'] = height

									if 'pensize' not in dom_dict :
										dom_dict['pensize'] = pensize

									BI = None

									if feature_data[4] != 'None' and feature_data[4] != '' :
										BI = float(feature_data[4])

									drawShape(draw, plot_x_pos+(int(feature_data[2]) - cut) * length_factor, plot_y_pos, dom_dict, BI)

							elif FT_type in ['HELIX', 'STRAND', 'TURN'] :
								feature_data = data_dict[protein][j]
								if feature in shape_dict :
									#If shapes defined in feature_shape file
									dom_dict = feature_shape_dict[feature_data[1]]
									if 'height' not in dom_dict :
										dom_dict['height'] = height

									if 'pensize' not in dom_dict :
										dom_dict['pensize'] = pensize

								# else :
								# 	#Use default if not
								# 	dom_dict = {}
								# 	dom_dict['pensize'] = pensize

								# 	if FT_type == 'HELIX' :
								# 		dom_dict['shape'] = 'rectangle'
								# 		dom_dict['contour_color'] = 'black'
								# 		dom_dict['color'] = '#FF4040'
								# 		dom_dict['height'] = height
								# 	elif FT_type == 'STRAND' :
								# 		dom_dict['shape'] = 'arrow'
								# 		dom_dict['orientation'] = 'right'
								# 		dom_dict['contour_color'] = 'black'
								# 		dom_dict['color'] = 'skyblue'
								# 		dom_dict['height'] = height
								# 	elif FT_type == 'TURN' :
								# 		dom_dict['shape'] = 'rhombus'
								# 		dom_dict['contour_color'] = 'black'
								# 		dom_dict['color'] = 'gold'
								# 		dom_dict['height'] = height

									dom_dict['size'] = int(float(feature_data[3])) * length_factor

									feature_data = data_dict[protein][j]

									BI = None

									if feature_data[4] != 'None' and feature_data[4] != '' :
										BI = float(feature_data[4])

									drawShape(draw, plot_x_pos+(int(feature_data[2]) - cut) * length_factor, plot_y_pos, dom_dict, BI)

							elif FT_type == 'ZN_FING' :
								feature_data = data_dict[protein][j]
								if feature in shape_dict :
									dom_dict = feature_shape_dict[feature_data[1]]
									if 'height' not in dom_dict :
										dom_dict['height'] = height

									if 'pensize' not in dom_dict :
										dom_dict['pensize'] = pensize

									dom_dict['size'] = int(float(feature_data[3])) * length_factor

									feature_data = data_dict[protein][j]

									BI = None

									if feature_data[4] != 'None' and feature_data[4] != '' :
										BI = float(feature_data[4])

									drawShape(draw, plot_x_pos+(int(feature_data[2]) - cut) * length_factor, plot_y_pos, dom_dict, BI)

							elif FT_type == 'PDB' :
								feature_data = data_dict[protein][j]
								if feature in shape_dict :
									dom_dict = feature_shape_dict[feature_data[1]]
									if 'height' not in dom_dict :
										dom_dict['height'] = height

									if 'pensize' not in dom_dict :
										dom_dict['pensize'] = pensize

								# else :
								# 	#Use default
								# 	dom_dict = {}
								# 	dom_dict['pensize'] = pensize
								# 	dom_dict['shape'] = 'line'
								# 	dom_dict['orientation'] = 'down'
								# 	dom_dict['contour_colormap'] = 'autumn_r'
								# 	dom_dict['height'] = height * 1.5


									dom_dict['size'] = int(float(feature_data[3])) * length_factor

									feature_data = data_dict[protein][j]

									BI = None

									if feature_data[4] != 'None' and feature_data[4] != '' :
										BI = float(feature_data[4])

									drawShape(draw, plot_x_pos+(int(feature_data[2]) - cut) * length_factor, plot_y_pos, dom_dict, BI)

							elif FT_type == 'DISORDER' :
								feature_data = data_dict[protein][j]

								if feature in shape_dict :
									dom_dict = feature_shape_dict[feature_data[1]]
									if 'height' not in dom_dict :
										dom_dict['height'] = height

									if 'pensize' not in dom_dict :
										dom_dict['pensize'] = pensize
									# else :
									# 	#Use default
									# 	dom_dict = {}
									# 	dom_dict['pensize'] = pensize
									# 	dom_dict['shape'] = 'rectangle'
									# 	dom_dict['contour_color'] = 'red'
									# 	dom_dict['height'] = 0

									dom_dict['size'] = int(float(feature_data[3])) * length_factor

									feature_data = data_dict[protein][j]

									BI = None

									if feature_data[4] != 'None' and feature_data[4] != '' :
										BI = float(feature_data[4])

									drawShape(draw, plot_x_pos+(int(feature_data[2]) - cut) * length_factor, plot_y_pos, dom_dict, BI)

							elif FT_type in ['REGION', 'DNA_BIND'] :

								#For each feature, draw the corresponding shape
								if feature in shape_dict :

									feature_data = data_dict[protein][j]

									dom_dict = feature_shape_dict[feature_data[1]]

									dom_dict['size'] = int(float(feature_data[3])) * length_factor

									if 'height' not in dom_dict :
										dom_dict['height'] = height

									if 'pensize' not in dom_dict :
										dom_dict['pensize'] = pensize

									BI = None

									if feature_data[4] != 'None' and feature_data[4] != '' :
										BI = float(feature_data[4])

									drawShape(draw, plot_x_pos+(int(feature_data[2]) - cut) * length_factor, plot_y_pos, dom_dict, BI)

								elif 'Others_R' in shape_dict :
									if feature not in shape_dict :
										feature_data = data_dict[protein][j]

										dom_dict = feature_shape_dict['Others_R']

										dom_dict['size'] = int(float(feature_data[3])) * length_factor

										if 'height' not in dom_dict :
											dom_dict['height'] = height

										if 'pensize' not in dom_dict :
											dom_dict['pensize'] = pensize

										BI = None

										if feature_data[4] != 'None' and feature_data[4] != '' :
											BI = float(feature_data[4])

										drawShape(draw, plot_x_pos+(int(feature_data[2]) - cut) * length_factor, plot_y_pos, dom_dict, BI)


		print('end of createFigure')
		# print(im.size)

		uuid_im = str(uuid.uuid4())

		im.save(os.path.join('Figures', uuid_im+'.png'), dpi=(600,600))

		children = html.Img(src=b64_image(os.path.join('Figures', uuid_im+'.png')), style={'width' : '100%', 'display' : 'flex', 'align-items' : 'center', 'justify-content' : 'center'})

		uuid_legend = drawLegend(pensize, height, text_size, shape_dict, length_factor, draw_region, biased_region_text_size, pdb_coverage_list)

		legend_children = html.Img(src=b64_image(os.path.join('Figures', uuid_legend+'.png')), style={'width' : '100%', 'display' : 'flex', 'align-items' : 'center', 'justify-content' : 'center'})

		return children, legend_children, uuid_im, uuid_legend, sorted_protein_list

	return no_update, no_update, no_update, no_update, no_update

def drawLegend(pensize, height, text_size, shape_dict, length_factor, draw_region, biased_region_text_size, pdb_coverage_list) :
	"""
		DESCRIPTION: Function handling all the legend drawing process
	"""

	# Calculate the height in pixels between two protein (defaut, 2 times the protein feature height)
	space_between_protein = height * 2

	shape_list = list(shape_dict.keys())

	im = Image.new('RGB', (4000, space_between_protein * len(shape_list)), 'white')
	draw = ImageDraw.Draw(im)
	text_font = ImageFont.truetype('arial-bold.ttf', text_size)

	len_treshold_list = len(pdb_coverage_list)
	value_list = [99/100 if i == 100 else i/100 for i in range(int(1/len_treshold_list*100), 101, int(1/len_treshold_list*100))]
	# print(pdb_coverage_list)
	# print(value_list)

	max_text_w = max([draw.textsize(shape, font=text_font)[0] for shape in shape_list if shape != None])

	for i in range(len(shape_list)) :
		#Get name of the protein if there is one
		shape = None
		try :
			shape = shape_list[i]
		except :
			shape = None

		#If there is a protein
		if shape != None :
			dom_dict = shape_dict[shape]
			if shape == 'Others' :
				shape = 'Others*'

			plot_x_pos = 20
			plot_y_pos = height + space_between_protein * i

			w, h = draw.textsize(shape, font=text_font)
			draw.text((10+max_text_w-w+plot_x_pos, plot_y_pos-(text_size/3)), shape, (0,0,0), font=text_font)

			draw.line([(50+max_text_w+plot_x_pos, plot_y_pos), (50+max_text_w+plot_x_pos + 200* length_factor, plot_y_pos)], fill='black', width=pensize, joint='curve')

			dom_dict['size'] = 100 * length_factor

			if 'height' not in dom_dict :
				dom_dict['height'] = height

			if 'pensize' not in dom_dict :
				dom_dict['pensize'] = pensize

			if shape == 'COMPBIAS' :
				drawShape(draw, 50+max_text_w+plot_x_pos + 50*length_factor, plot_y_pos, dom_dict, 0.5, biased_region_text='H', biased_region_text_size=biased_region_text_size)
			else :
				drawShape(draw, 50+max_text_w+plot_x_pos + 50*length_factor, plot_y_pos, dom_dict, 0.5)

			# if draw_region != [] :
			# 	drawRegion(draw, 50+max_text_w+plot_x_pos + 50*length_factor, plot_y_pos, 'H', biased_region_text_size, data_dict, protein, j, height, cut, length_factor, pensize)


	# Add colormaps on figure
	for i in range(len(shape_list)) :
		#Get name of the protein if there is one
		shape = None
		try :
			shape = shape_list[i]
		except :
			shape = None

		#If there is a protein
		if shape != None :
			dom_dict = shape_dict[shape]

			plot_x_pos = 20
			plot_y_pos = int(height/2) + space_between_protein * i

			if ('colormap' in dom_dict and dom_dict['colormap']!= 'None' and dom_dict['colormap']!= np.nan) or ('contour_colormap' in dom_dict and dom_dict['contour_colormap'] != 'None' and dom_dict['contour_colormap']!= np.nan) :

				colormap = None
				contour_colormap = None

				if ('colormap' in dom_dict and dom_dict['colormap']!= 'None' and dom_dict['colormap']!= np.nan) :
					colormap = dom_dict['colormap']
				if ('contour_colormap' in dom_dict and dom_dict['contour_colormap']!= 'None' and dom_dict['contour_colormap']!= np.nan) :
					contour_colormap = dom_dict['contour_colormap']

				# print(colormap, contour_colormap)
				nb_colormap = len([cmap for cmap in [colormap, contour_colormap] if cmap != None])
				# print('nb_colormap', nb_colormap)
				if colormap == contour_colormap :
					nb_colormap = 1

				for j in range(nb_colormap) :
					cmap_list = [cmap for cmap in [colormap, contour_colormap] if cmap != None]
					try :
						cmap_image = Image.open(os.path.join('Colormaps', cmap_list[j]+'.png'))

						cmap_image = cmap_image.resize((height*8, int(height/2)))

						im.paste(cmap_image, (50+max_text_w+plot_x_pos + 200* length_factor + (height*8+50)*j + 50, plot_y_pos))
					except : pass

	draw = ImageDraw.Draw(im)

	# Add sclaes on colormaps
	for i in range(len(shape_list)) :
		#Get name of the protein if there is one
		shape = None
		try :
			shape = shape_list[i]
		except :
			shape = None

		#If there is a protein
		if shape != None :
			dom_dict = shape_dict[shape]

			plot_x_pos = 20
			plot_y_pos = height + space_between_protein * i

			if ('colormap' in dom_dict and dom_dict['colormap']!= 'None' and dom_dict['colormap']!= np.nan) or ('contour_colormap' in dom_dict and dom_dict['contour_colormap'] != 'None' and dom_dict['contour_colormap']!= np.nan) :

				colormap = None
				contour_colormap = None

				if ('colormap' in dom_dict and dom_dict['colormap']!= 'None' and dom_dict['colormap']!= np.nan) :
					colormap = dom_dict['colormap']
				if ('contour_colormap' in dom_dict and dom_dict['contour_colormap']!= 'None' and dom_dict['contour_colormap']!= np.nan) :
					contour_colormap = dom_dict['contour_colormap']

				# print(colormap, contour_colormap)
				nb_colormap = len([cmap for cmap in [colormap, contour_colormap] if cmap != None])
				# print('nb_colormap', nb_colormap)
				if colormap == contour_colormap :
					nb_colormap = 1

				for j in range(nb_colormap) :

					draw.line([(50+max_text_w+plot_x_pos + int(200* length_factor) + (int(height*8)+50)*j + 50, plot_y_pos), (50+max_text_w+plot_x_pos + int(200* length_factor) + (height*8+50)*j + 50 + height*8, plot_y_pos)], fill='black', width=pensize, joint='curve')

					for tick in range(50+max_text_w+plot_x_pos + int(200* length_factor) + (int(height*8)+50)*j + 50, 50+max_text_w+plot_x_pos + int(200* length_factor) + (height*8+50)*j + 50 + height*8 + 1, int(height*8/10)) :
						draw.line([(tick, plot_y_pos), (tick, plot_y_pos+int(height/8))], fill='black', width=pensize, joint='curve')

	uuid_legend = str(uuid.uuid4())
	im.save(os.path.join('Figures', uuid_legend+'.png'), 'PNG')

	return uuid_legend

def drawShape(draw, x_pos, y_pos, dom_dict, BI=None, biased_region_text=None, biased_region_text_size=None) :
	"""
		DESCRIPTION: Function handling the drawing of a single feature
	"""

	dom_dict_copy = copy.deepcopy(dom_dict)

	# print(dom_dict_copy)

	colormap = ''
	threshold = 0.0

	contour_colormap = ''
	contour_threshold = 0.0

	if 'color' in dom_dict_copy :
		color = dom_dict_copy['color']
		if color == '' or color == 'nan' or color == 'none':
			dom_dict_copy['color'] = 'white'
	if 'contour_color' in dom_dict_copy :
		contour_color = dom_dict_copy['contour_color']
		if contour_color == '' or contour_color == 'nan' or contour_color == 'none':
			dom_dict_copy['contour_color'] = 'black'

	if 'colormap' in dom_dict :
		colormap = dom_dict_copy.pop('colormap')
	if 'threshold' in dom_dict :
		threshold = dom_dict_copy.pop('threshold')

	if 'contour_colormap' in dom_dict :
		contour_colormap = dom_dict_copy.pop('contour_colormap')
	if 'contour_threshold' in dom_dict :
		contour_threshold = dom_dict_copy.pop('contour_threshold')

	#Set the color from the colormap with the BI values. If no value is given (BI = -1), the shape will be colored grey
	if BI != None and colormap != '' and colormap != 'nan' and colormap != 'None':
		# print(BI, type(BI))
		if BI >= threshold :
			cmap = pl.cm.get_cmap(colormap)
			# print('colormap', colormap)
			color = plt.colors.rgb2hex(cmap(BI)[:3])
			# print('color', contour_colormap)
			dom_dict_copy['color'] = color
		if BI == -1 :
			dom_dict_copy['color'] = 'grey' #Color for missing values

	if BI != None and contour_colormap != '' and contour_colormap != 'nan' and contour_colormap != 'None':
		if BI >= contour_threshold :
			cmap = pl.cm.get_cmap(contour_colormap)
			contour_color = plt.colors.rgb2hex(cmap(BI)[:3])
			dom_dict_copy['contour_color'] = contour_color
		elif BI == -1 :
			dom_dict_copy['contour_color'] = 'black' #Color for missing value


	drawForm(draw, x_pos, y_pos, **dom_dict_copy)

	if biased_region_text != None and biased_region_text_size != None :

		text_font = ImageFont.truetype('arial-bold.ttf', biased_region_text_size)
		draw.text((x_pos, y_pos-biased_region_text_size), biased_region_text, (0,0,0), font=text_font)

def findData(case_list, feature_list, intensity_data, db_PDB, features_df):
	"""
		DESCRIPTION: Merging extracted data and numerical values to make the final DataFrame for the map drawing process.
	"""

	output_df = None

	if case_list != ['noCase'] :
		if case_list == ['PDB'] :
			data_dict = loadIntensities(case_list, feature_list, intensity_data, db_PDB)
		else :
			data_dict = loadIntensities(case_list, feature_list, intensity_data, db_PDB, with_feature_name=True)

		output_dict_list = []
		header_list = ['protein', 'feature_type', 'feature', 'start', 'length']+case_list

		for index, row in features_df.iterrows():
			feature = row['feature']
			if feature in feature_list :
				protein = row['protein']
				start = row['start']
				value = '-1'

				value_list = []
				for case in case_list :
					if protein in data_dict[feature][case] :
						data_list = data_dict[feature][case][protein]

						if case_list == ['PDB'] :
							res_list = [elem[1] for elem in data_list if start == elem[0]]
						else :
							res_list = [elem[2] for elem in data_list if start == elem[1] and feature == elem[0]]

						if len(res_list) > 0 :
							value = str(res_list[0])
							if res_list[0] == np.nan : value == str(-1)

							value_list.append(value)

				output_dict_list.append({x[0] : x[1] for x in zip(header_list, [protein, row['feature_type'], feature, start, row['length']] + value_list)})

			else :
				output_dict_list.append({x[0] : x[1] for x in zip(header_list, [row['protein'], row['feature_type'], row['feature'], row['start'], row['length']] + [np.nan for i in range(len(case_list))])})

		print('-----------------------------')
		output_df = pd.DataFrame.from_dict(output_dict_list)
		# print(output_df)

	else :

		features_df['noCase'] = features_df['value']

	return output_df

def loadIntensities(case_list, feature_list, intensity_data, db_PDB, with_feature_name=False) :
	"""
		DESCRIPTION: Reading numerical values of the 3D structure coverage and user (optionnal) and tranforming it into a dict.
	"""

	data_dict = {}
	# print(case_list, feature_list)
	for feature in feature_list :
		data_dict[feature] = {}
		for case in case_list :
			data_dict[feature][case] = {}

		df = None

		if feature == 'PDB' :
			df = pd.read_json(db_PDB, orient='split')
		else :
			df = intensity_data

		for index, row in df.iterrows():
			protein = row['protein']
			start = row['start']
			for case in case_list :

				value = row[case]

				if with_feature_name :
					try :
						data_dict[feature][case][protein].append((row['feature'], int(start), float(value)))
					except :
						data_dict[feature][case][protein] = [(row['feature'], int(start), float(value))]
				else :
					try :
						data_dict[feature][case][protein].append((int(start), float(value)))
					except :
						data_dict[feature][case][protein] = [(int(start), float(value))]

	# print(data_dict)

	return data_dict

def autoFeatChoice(occurrence_df, protein_nb=None, auto_toggle='more than x', auto_select_threshold=None) :
	"""
		DESCRIPTION: Choosing shapes and colors for features based on the feature occurrence (only DOMAIN, REPEAT, REGION and MOTIF). Most common features (domains and repeats, >=10 and >= 50 respectively in the human proteome) use a default shape and color. Other feature ahve a random shape and color selected.
			Two options available:
				x most represented: Selecting the x features that have the most occurrences in the extraced data. Default x: math.sqrt(number of proteins in the list)
				more than x occurrence : Selecting all features that appear at least x number of time in the protein list. Default x: 2
			x: Can also be defined by the user ('auto_select_threshold')
	"""

	if protein_nb != None :
		protein_nb = int(protein_nb.split(': ')[1])
	# print('Protein number : ', protein_nb)

	contour_color_list = ['orange', 'red', 'violet', 'purple', 'blue', 'skyblue', 'green', 'gray']
	fill_color_list = ['gold', 'violet', 'purple', 'skyblue', 'lightgreen', 'black', 'gray']
	shape_list = [('arrow', 'left'), ('arrow', 'right'), ('triangle', 'left'), ('triangle', 'right'), ('rhombus', ''), ('hexagon', 'horizontal'), ('hexagon', 'vertical'), ('rectangle', ''), ('semirectangle', 'up'), ('semirectangle', 'down')]

	data_dict = {}
	feature_type = ['DOMAIN', 'CHAIN', 'INIT_MET', 'PEPTIDE', 'ZN_FING', 'DNA_BIND', 'REGION', 'ACT_SITE', 'METAL', 'SITE', 'LIPID', 'HELIX', 'STRAND', 'TURN', 'CONFLICT', 'CARBOHYD', 'BINDING', 'MOTIF', 'MOD_RES', 'COMPBIAS', 'REPEAT', 'VARIANT']
	for feature in feature_type :
		data_dict[feature] = {}

	column_names = ['feature', 'shape', 'orientation', 'height', 'contour_color', 'contour_colormap', 'contour_threshold', 'color', 'colormap', 'threshold', 'pensize']

	shape_df = pd.DataFrame(columns=column_names)

	occurrence_df = occurrence_df.sort_values(by='occurrence', ascending=False)
	# print(occurrence_df)

	# print(auto_toggle, auto_select_threshold)
	if auto_toggle == [] or auto_toggle == None :
		# x Most represented
		if auto_select_threshold == None or auto_select_threshold.strip() == '':
			# No nb given, select default value (sqrt(protein_nb))
			auto_select_threshold = int(protein_nb**(1/2))
		else :
			# Take given value
			auto_select_threshold = int(auto_select_threshold)

		nb_feature = 0
		shape_id = 0
		for index, row in occurrence_df.iterrows():
			if row['feature_type'] in ['DOMAIN', 'REPEAT', 'MOTIF', 'REGION'] :
				if nb_feature < auto_select_threshold :
					feature_name = row['feature']
					if feature_name != 'Disordered' :
						shape, orientation = rd.choice(shape_list)
						contour_color = rd.choice(contour_color_list)
						color = rd.choice(fill_color_list)
						if row['feature_type'] == 'MOTIF' :
							shape_df.loc[shape_id] = [feature_name, 'bar', None, None, contour_color, None, None, contour_color, None, None, 5]
						else :
							shape_df.loc[shape_id] = [feature_name, shape, orientation, None, None, None, None, color, None, None, None]
						shape_id += 1

						nb_feature += 1

		shape_df.loc[shape_id] = ['Others', 'rectangle', None, None, 'black', None, None, 'white', None, None, None]

	else :
		# More than x occurrence
		if auto_select_threshold == None or auto_select_threshold.strip() == '':
			# No nb given, select default value (sqrt(protein_nb))
			auto_select_threshold = 2
		else :
			# Take given value
			auto_select_threshold = int(auto_select_threshold)

		shape_id = 0
		for index, row in occurrence_df.iterrows():
			if row['feature_type'] in ['DOMAIN', 'REPEAT', 'MOTIF', 'REGION'] and row['occurrence'] >= auto_select_threshold :
				feature_name = row['feature']
				if feature_name != 'Disordered' :
					shape, orientation = rd.choice(shape_list)
					contour_color = rd.choice(contour_color_list)
					color = rd.choice(fill_color_list)
					if row['feature_type'] == 'MOTIF' :
						shape_df.loc[shape_id] = [feature_name, 'bar', None, None, contour_color, None, None, contour_color, None, None, 5]
					else :
						shape_df.loc[shape_id] = [feature_name, shape, orientation, None, None, None, None, color, None, None, None]
					shape_id += 1

		shape_df.loc[shape_id] = ['Others', 'rectangle', None, None, 'black', None, None, 'white', None, None, None]

		print('auto feat choice')
		# print(shape_df)


	# Replaces most represented features (DOMAINS and REPEATS) by their predefined colors
	most_common_feature_list = list(most_repr_feat_shape_df.index)
	column_list = list(shape_df)

	final_feature_dict_list = []

	for index, row in shape_df.iterrows():
		feature = row['feature']
		if feature in most_common_feature_list :
			single_shape_dict = {}
			single_shape_dict['feature'] = feature

			single_shape_dict = {**single_shape_dict, **most_repr_feat_shape_df.loc[feature].to_dict()}

			final_feature_dict_list.append(single_shape_dict)
		else :
			final_feature_dict_list.append({column : row[column] for column in column_list})

	shape_df = pd.DataFrame.from_dict(final_feature_dict_list)


	return shape_df

def prepareMapParameters(length_factor, height, text_size, sorting_dropdown, focus_input, threshold, biased_region_text_size, uniform_toggles, draw_region_dropdown, coverage_toggles, second_struct_toggles, disorder_toggles, mod_res_toggles, protein_length_toggles, pensize, FT_order_list) :

	if length_factor == None :
		length_factor = 1

	if height == None :
		height = 20

	if text_size == None :
		text_size = 30

	if biased_region_text_size == None :
		biased_region_text_size = 20

	draw_region = draw_region_dropdown

	switch_toggles = coverage_toggles + second_struct_toggles + disorder_toggles + mod_res_toggles + protein_length_toggles

	if 'coverage' in switch_toggles :
		draw_coverage = 'True'
	else :
		draw_coverage = 'False'

	if '2nd_struct' in switch_toggles :
		draw_2nd_struct = 'True'
	else :
		draw_2nd_struct = 'False'

	if 'disorder' in switch_toggles :
		draw_disorder = 'True'
	else :
		draw_disorder = 'False'

	if 'mod_res' in switch_toggles :
		draw_mod_res = 'True'
	else :
		draw_mod_res = 'False'

	if 'show_length' in switch_toggles :
		draw_protein_length = 'True'
	else :
		draw_protein_length = 'False'

	if pensize == None :
		pensize = 0

	if threshold == None :
		threshold = 0

	if sorting_dropdown == 'None' :
		sorting_dropdown = None

	if focus_input == 'None' or focus_input == None or focus_input == '' :
		focus_input = None

	parameter_dict = {}

	parameter_dict['length_factor'] = length_factor
	parameter_dict['height'] = height
	parameter_dict['text_size'] = text_size
	parameter_dict['sorting'] = sorting_dropdown
	parameter_dict['focus'] = focus_input
	parameter_dict['threshold'] = threshold
	parameter_dict['biased_region_text_size'] = biased_region_text_size
	parameter_dict['uniform'] = uniform_toggles
	parameter_dict['draw_region'] = draw_region
	parameter_dict['draw_coverage'] = draw_coverage
	parameter_dict['draw_2nd_struct'] = draw_2nd_struct
	parameter_dict['draw_disorder'] = draw_disorder
	parameter_dict['draw_mod_res'] = draw_mod_res
	parameter_dict['draw_protein_length'] = draw_protein_length
	parameter_dict['pensize'] = pensize
	parameter_dict['FT_order_list'] = FT_order_list

	return parameter_dict

def getMissingProteinNames(protein_list_df) :
	full_entry_list = list(protein_list_df['code'])

	protein_number_dict = {}
	protein_name_list = []
	for entry in full_entry_list :
		f = open(os.path.join('Uniprot_files', entry + '.txt'))
		all_lines = f.readlines()
		f.close()

		for line in all_lines :
			sline = line.strip().split('{')[0].split()

			if sline[0] == 'ID' :
				protein_name = sline[1]

				if protein_name not in protein_number_dict :
					protein_number_dict[protein_name] = 1
				else :
					protein_number_dict[protein_name] += 1
					protein_name = protein_name+'-dup_'+str(protein_number_dict[protein_name]-1)

				protein_name_list.append(protein_name)
				break
	return protein_name_list

def proteinDataGathering(complete_protein_list_df, dl_latest=None, uniprot_availability_test=False) :

	if not os.path.exists('Uniprot_files') :
		print('Creating Uniprot_files folder')
		os.mkdir('Uniprot_files')

	protein_list_df = pd.DataFrame()
	if uniprot_availability_test == False :
		if 'Entry' in list(complete_protein_list_df) :
			protein_list_df['code'] = complete_protein_list_df['Entry']
		elif 'code' in list(complete_protein_list_df) :
			protein_list_df['code'] = complete_protein_list_df['code']
		if 'Entry name' in list(complete_protein_list_df)  :
			protein_list_df['protein'] = complete_protein_list_df['Entry name']
		elif 'Entry Name' in list(complete_protein_list_df)  :
			protein_list_df['protein'] = complete_protein_list_df['Entry Name']
		elif 'protein' in list(complete_protein_list_df) :
			protein_list_df['protein'] = complete_protein_list_df['protein']

		print(protein_list_df)

		protein_list = list(protein_list_df['code'])
		protein_list = [protein_name.strip() for protein_name in protein_list if protein_name.strip() != '']

	else :
		protein_list = complete_protein_list_df

	if protein_list != []:
		unique_protein_list = list(set(protein_list))

		uniprot_db = getUniprotDBData(protein_list)

		if dl_latest == ['latest'] :
			to_add_protein_list = updateTargetsUniprotDBData(uniprot_db, unique_protein_list, outdated_thres=0)
		else :
			to_add_protein_list = updateTargetsUniprotDBData(uniprot_db, unique_protein_list, outdated_thres=3)

		# Clear all files and db file
		# file = open('uniprot_db.csv', 'w')
		# file.write('Uniprot_code\tlatest_download\tdownload_nb\n')
		# file.close()
		# out_dir_list = os.listdir('Uniprot_files')
		# for item in out_dir_list:
		# 	if item.endswith(".txt"):
		# 		os.remove(os.path.join('Uniprot_files', item))

		nb_Uniprot = len(set(to_add_protein_list))
		i_Uniprot = 1

		not_dled_df = pd.DataFrame(columns=['Unrecognized code'])
		not_dled_df_i = 0

		#For each protein
		for protein in to_add_protein_list :
			perc_finished = round(i_Uniprot/nb_Uniprot*100, 1)
			print(i_Uniprot, '/', nb_Uniprot, '->', perc_finished, '%')
			i_Uniprot += 1
			file = protein.strip() + ".txt"
			print("Downloading " + file)

			if protein not in uniprot_db :
				uniprot_db[protein] = {'latest' : datetime.today().strftime('%Y-%m-%d'), 'nb_dl' : 1}
			else :
				uniprot_db[protein] = {'latest' : datetime.today().strftime('%Y-%m-%d'), 'nb_dl' : uniprot_db[protein]['nb_dl'] + 1}

			url = 'https://www.uniprot.org/uniprot/' + file
			try :
				urllib.request.urlretrieve(url, os.path.join('Uniprot_files', file)) #Make the request for Uniprot file
			except :
				not_dled_df.loc[not_dled_df_i] = [protein]
				not_dled_df_i += 1

		updateUniprotDBData(uniprot_db)

		obsolete_df = pd.DataFrame(columns=['Obsolete code'])
		obsolete_df_i = 0
		for protein in protein_list :
			# print(protein)
			# print(list(not_dled_df['Unrecognized code']))
			if protein not in list(not_dled_df['Unrecognized code']) :
				file = protein.strip() + ".txt"
				if os.path.getsize(os.path.join('Uniprot_files', file)) == 0 :
					obsolete_df.loc[obsolete_df_i] = [protein]
					obsolete_df_i += 1

		# HAS to check of obsolete codes even if file already exists

		# print(obsolete_df)
		# print(not_dled_df)

		obsolete_data = obsolete_df.to_json(date_format='iso', orient='split')
		unrecognized_data = not_dled_df.to_json(date_format='iso', orient='split')
		if obsolete_df_i != 0 or not_dled_df_i != 0 :
			return obsolete_data, unrecognized_data, {'display': 'block'}, b64_image(os.path.join('GUI', 'fail'+'.png')), None, None
		else :
			if uniprot_availability_test == False :
				missing_protein_names = False
				if 'protein' not in list(protein_list_df) :
					print('No protein names detected')
					missing_protein_names = True

					protein_name_list = getMissingProteinNames(protein_list_df)

					protein_list_df['protein'] = protein_name_list
					print(protein_list_df)

				return obsolete_data, unrecognized_data, {'display': 'none'}, b64_image(os.path.join('GUI', 'success'+'.png')), missing_protein_names, protein_list_df
			else :
				return None, None, {'display': 'none'}, b64_image(os.path.join('GUI', 'success'+'.png')), None, None

	else :
		print('No proteins in protein list')
		return no_update, no_update, no_update, b64_image(os.path.join('GUI', 'fail'+'.png'))

def featureExtraction(protein_list_df, modification_df=None, pdb_threshold_list=[1,2,3,5,10]) :
	if not isinstance(modification_df, pd.DataFrame) :
		modification_df = pd.DataFrame()

	exception_dict = {}
	exception_dict['remove'] = {}
	exception_dict['add'] = {}

	for index, row in modification_df.iterrows():
		if row['ex_type'] == '-' or row['ex_type'] == 'remove' :
			exception_dict['remove'][index] = {}
			exception_dict['remove'][index]['protein'] = row['protein']
			exception_dict['remove'][index]['feature_type'] = row['feature_type']
			exception_dict['remove'][index]['feature'] = row['feature']
		elif row['ex_type'] == '+' or row['ex_type'] == 'add' :
			exception_dict['add'][index] = {}
			exception_dict['add'][index]['protein'] = row['protein']
			exception_dict['add'][index]['feature_type'] = row['feature_type']
			exception_dict['add'][index]['feature'] = row['feature']
			exception_dict['add'][index]['start'] = int(row['start'])
			exception_dict['add'][index]['length'] = int(row['length'])

	# print(exception_dict)

	uni_protein_dict = {}

	# print(df)

	full_entry_list = list(protein_list_df['code'])
	full_base_protein_list = list(protein_list_df['protein'])

	code_dict = {}
	for i in range(len(full_entry_list)) :
		code_dict[full_base_protein_list[i]] = full_entry_list[i]
		code_dict[full_entry_list[i]] = full_base_protein_list[i]

	entry_list = []
	base_protein_list = []

	for i in range(max(len(full_entry_list), len(full_base_protein_list))) :
		if full_base_protein_list[i].strip() != '' and full_entry_list[i].strip() != '' :
			uni_protein_dict[full_entry_list[i]] = full_base_protein_list[i]
			uni_protein_dict[full_base_protein_list[i]] = full_entry_list[i]
			entry_list.append(full_entry_list[i].strip())
			base_protein_list.append(full_base_protein_list[i].strip())

	# print(uni_protein_dict)

	file_list = [os.path.join('Uniprot_files', uniprot_code+'.txt') for uniprot_code in list(set(list(uni_protein_dict.keys())))]

	protein_length_dict_list = []
	protein_length_dict = {}

	pdb_list_dict_list = []

	feature_df_dict_list = []

	base_protein_list_len = len(base_protein_list)
	current_protein_id = 1

	for protein in base_protein_list :

		print('Extracting:', str(current_protein_id)+'/'+str(base_protein_list_len))
		#Read file
		f = open(os.path.join('Uniprot_files', uni_protein_dict[protein] + '.txt'))
		all_lines = f.readlines()
		f.close()

		#Get name of corresponding protein
		#protein = uni_protein_dict[os.path.split(file)[1][0:-4]]

		name = 'ERROR'
		feature_length = None
		FT_type = None
		wline = ""
		is_exception = False
		get_name = False
		feature_length = None
		feature_start = None
		feature_full_name = None
		variant = ''

		for line in all_lines :
			sline = line.strip().split('{')[0].split()

			if len(sline) >= 5 and sline[0] == 'DR' and sline[1].replace(';', '') == 'PDB' :
				# print(sline)
				PDB = sline[2].replace(';', '')
				Method = sline[3].replace(';', '')
				Resolution = sline[4].replace(';', '')
				Chains = 'Unknown'
				StartEnd = 'Unknown'
				Start = 'Unknown'
				End = 'Unknown'
				try :
					Chains, StartEnd = sline[-1].replace('.', '').split('=')
					Start, End = StartEnd.split('-')
					Start = int(Start)
					End = int(End)
					Resolution = float(Resolution)
				except : pass

				pdb_list_dict_list.append({'code' : code_dict[protein], 'protein' : protein, 'PDB' : PDB, 'Method' : Method, 'Resolution' : Resolution, 'Chains' : Chains, 'Start' : Start, 'End' : End})

			if get_name == False :
				name = 'ERROR'
				feature_length = None
				FT_type = None
				wline = ""
				is_exception = False
				feature_length = None
				feature_start = None
				feature_full_name = None
				variant = ''

			if get_name == True :
				feature_full_name = 'None'
				try :
					feature_full_name = line.strip('\n').strip(" ").split('/note=')[1].strip('"')
				except : pass

				if FT_type in ['DOMAIN', 'CHAIN', 'INIT_MET', 'PEPTIDE', 'ZN_FING', 'DNA_BIND', 'REGION', 'ACT_SITE', 'METAL', 'SITE', 'LIPID', 'HELIX', 'STRAND', 'TURN', 'CARBOHYD', 'BINDING'] : # == 'DOMAIN' or FT_type == 'REGION' or FT_type == 'DNA_BIND' or FT_type == 'ZN_FING' :
					full_feature_name = ' '.join(feature_full_name.split(' ')).replace('.', '').replace(';', '')
					if len(full_feature_name.split(' ')) > 1 :
						try :
							int(full_feature_name.split(' ')[-1])
							name = ' '.join(full_feature_name.split(' ')[:-1])
						except :
							name = full_feature_name
						if 'EGF-like' in name and 'calcium-binding' in name :
							name_elem_list = full_feature_name.split(' ')
							name = name_elem_list[0] + ' ' + name_elem_list[-1]
					else :
						name = full_feature_name

					if FT_type == 'REGION' and full_feature_name == 'Disordered' :
						FT_type = 'DISORDER'

					if FT_type == 'ZN_FING' : name = 'ZN_FING'

					#For secondary structure, put according name
					if FT_type == 'HELIX' : name = 'Helix'
					elif FT_type == 'STRAND' : name = 'Strand'
					elif FT_type == 'TURN' : name = 'Turn'


				elif FT_type == 'MOTIF' :
					full_feature_name = ' '.join(feature_full_name.split(' ')).split('{')[0].split('.')[0]
					try :
						int(full_feature_name.split('_')[-1])
						name = ' '.join(full_feature_name.split('_')[:-1])
					except :
						name = full_feature_name
				elif FT_type == 'MOD_RES' :
					name = ' '.join(feature_full_name.split(' ')).split(';')[0].split('.')[0]
				elif FT_type == 'COMPBIAS' :
					name = ' '.join(feature_full_name.split(' ')).split('.')[0]
				elif FT_type == 'REPEAT' :
					full_feature_name = feature_full_name.replace('.', '').replace(' ', '_')
					if len(full_feature_name.split('_')) > 1 :
						try :
							int(full_feature_name.split('_')[-1])
							name = ' '.join(full_feature_name.split('_')[:-1])
						except :
							name = full_feature_name
					else :
						try :
							int(full_feature_name)
							name = 'Repeated element of protein'
						except :
							name = full_feature_name

				elif FT_type == 'VARIANT' or FT_type == 'CONFLICT':
					name = 'None'
					if FT_type == 'VARIANT' : name = 'Variant'
					elif FT_type == 'CONFLICT' : name = 'Conflict'
					variant = ' '.join(feature_full_name.split(' ')).replace('.', '').split('(')[0].replace('_', '')
					feature_length = 1


				get_name = False

				for elem in exception_dict['remove'] :
					ex_protein = exception_dict['remove'][elem]['protein']
					ex_FT_type = exception_dict['remove'][elem]['feature_type']
					ex_name = exception_dict['remove'][elem]['feature']
					if protein == ex_protein and FT_type == ex_FT_type and name == ex_name :
						# print("#################### REMOVED #####################")
						# print(protein + "\t" + FT_type + "\t" + name + "\t" + str(feature_start) + "\t" + str(feature_length) + "\t" + "None" + "\n")
						is_exception = True

				if is_exception == False :
					if variant == '' :
						feature_df_dict_list.append({'code' : code_dict[protein], 'protein' : protein, 'feature_type' : FT_type, 'feature' : name, 'start' : feature_start, 'length' : feature_length, 'intensity' : np.nan, 'variant' : np.nan})
					else :
						feature_df_dict_list.append({'code' : code_dict[protein], 'protein' : protein, 'feature_type' : FT_type, 'feature' : name, 'start' : feature_start, 'length' : feature_length, 'intensity' : np.nan, 'variant' : variant})

			elif len(sline) >= 2 :
				FT_type = sline[1]
				#If the line contains information about a domain/composition biased region/phosphosite add it to file (with start position and length)
				if sline[0] == "FT" and FT_type in ['DOMAIN', 'CHAIN', 'INIT_MET', 'PEPTIDE', 'ZN_FING', 'DNA_BIND', 'REGION', 'ACT_SITE', 'METAL', 'SITE', 'LIPID', 'HELIX', 'STRAND', 'TURN', 'CONFLICT', 'CARBOHYD', 'BINDING', 'MOTIF', 'MOD_RES', 'COMPBIAS', 'REPEAT', 'VARIANT'] :
					#print(sline)
					ft_elem = sline[2].replace('<','').replace('>','').split('..')
					nb_elem = len(ft_elem)
					if nb_elem == 1 :
						feature_start = sline[2]
						feature_length = 1
					else :
						try :
							feature_start, feature_stop = [int(elem) for elem in ft_elem]
							feature_length = feature_stop - feature_start + 1 #int(sline[3]) - int(sline[2]) + 1
						except :
							pass

					get_name = True

			#Get length of the full protein
			if len(sline) != 0 and sline[0] == "SQ" and len(sline) > 1 :
				protein_length_dict_list.append({'code' : code_dict[protein], 'protein' : protein, 'total_length' : int(sline[2])})

				protein_length_dict[protein] = int(sline[2])

		current_protein_id += 1

	for protein in base_protein_list :
		for exception in exception_dict['add'] :
			if exception_dict['add'][exception]['protein'] == protein :
				# print("##################### ADDED ######################")
				# print(protein + "\t" + exception_dict['add'][exception]['feature_type'] + "\t" + exception_dict['add'][exception]['feature'] + "\t" + str(exception_dict['add'][exception]['start']) + "\t" + str(exception_dict['add'][exception]['length']) + "\t" + "None" + "\n")
				feature_df_dict_list.append({'code' : code_dict[protein], 'protein' : protein, 'feature_type' : exception_dict['add'][exception]['feature_type'], 'feature' : exception_dict['add'][exception]['feature'], 'start' : exception_dict['add'][exception]['start'], 'length' : exception_dict['add'][exception]['length'], 'intensity' : np.nan, 'variant' : np.nan})

	features_df = pd.DataFrame.from_dict(feature_df_dict_list)

	protein_length_df = pd.DataFrame.from_dict(protein_length_dict_list)

	protein_occurrence_dict_list = []

	feature_list = list(set(list(features_df['feature'])))

	name_dict = {}
	feature_dict = {}
	for feature in feature_list :
		name_dict[feature] = 0

	for index, row in features_df.iterrows():
		name_dict[row['feature']] += 1
		feature_dict[row['feature']] = row['feature_type']

	for feature in name_dict :
		protein_occurrence_dict_list.append({'feature_type' : feature_dict[feature], 'feature' : feature, 'occurrence' : name_dict[feature]})

	protein_occurrence_df = pd.DataFrame.from_dict(protein_occurrence_dict_list)

	pdb_list_df = pd.DataFrame.from_dict(pdb_list_dict_list)

	seq_dict = {}

	for protein in base_protein_list :
		seq = ''
		for line in open(os.path.join('Uniprot_files', uni_protein_dict[protein] + '.txt')) :
			sline = line.strip('\n').split()
			if (len(sline) > 1 and len(sline[0]) > 2) or (len(sline) == 1 and sline[0] != '//') :
				seq += ''.join(sline)
		seq_dict[protein] = seq

	protein_seq_str = ''

	for protein in seq_dict :
		protein_seq_str += '>' + protein + '\n'
		for start in range(0, len(seq_dict[protein]), 60) :
			try :
				protein_seq_str += seq_dict[protein][start:start+60] + '\n'
			except :
				protein_seq_str += seq_dict[protein][start:] + '\n'

	coverage_dict_list = []

	knowledge_dict = {}
	for protein in base_protein_list :
		knowledge_dict[protein] = [0 for i in range(protein_length_dict[protein])]

	protein_i = 0
	for index, row in pdb_list_df.iterrows():
		protein = row['protein']
		protein_i += 1
		# print(protein_i)
		if protein in base_protein_list :
			try :
				start = int(row['Start'])
				end = int(row['End'])
			except :
				start = 1
				for elem in protein_length_dict_list :
					if elem['protein'] == protein :
						end = int(elem['total_length'])

			for i in range(start, end) :
				knowledge_dict[protein][i] += 1

	for protein, coverage in knowledge_dict.items() :
		# print(protein, coverage)
		last_value = 0
		start = 1
		end = 1
		for i in range(1, len(coverage)) :
			last_value = coverage[i-1]
			value = coverage[i]

			if value != last_value :
				end = i
				coverage_dict_list.append({'protein' : protein, 'start' : start, 'end' : end, 'length' : end-start+1, 'coverage' : last_value})
				start = i+1

			elif i == len(coverage)-1 :
				end = i+1
				coverage_dict_list.append({'protein' : protein, 'start' : start, 'end' : end, 'length' : end-start+1, 'coverage' : last_value})

			if i == len(coverage)-1 and coverage[-1] != coverage[-2] and coverage[-1] != 0 :
				start = len(coverage)
				end = len(coverage)
				last_value = coverage[-1]
				coverage_dict_list.append({'protein' : protein, 'start' : start, 'end' : end, 'length' : end-start+1, 'coverage' : last_value})

	# print(coverage_df)
	coverage_df = pd.DataFrame.from_dict(coverage_dict_list)

	db_coverage_dict_list = []

	pdb_threshold_list = [int(elem) for elem in pdb_threshold_list]
	len_treshold_list = len(pdb_threshold_list)

	for index, row in coverage_df.iterrows():
		if row['protein'] != 'protein' :
			coverage = int(row['coverage'])
			val_list = [99/100 if i == 100 else i/100 for i in range(int(1/len_treshold_list*100), 101, int(1/len_treshold_list*100))]
			value = 0
			for i in range(len_treshold_list) :
				if coverage >= pdb_threshold_list[i] :
					value = val_list[i]

			db_coverage_dict_list.append({'protein' : row['protein'], 'feature' : np.nan, 'start' : row['start'], 'PDB' : value})

			if value != 0 :
				feature_df_dict_list.append({'code' : code_dict[row['protein']], 'protein' : row['protein'], 'feature_type' : 'PDB', 'feature' : 'PDB', 'start' : row['start'], 'length' : row['length'], 'intensity' : str(value), 'variant' : np.nan})

	db_coverage_df = pd.DataFrame.from_dict(db_coverage_dict_list)

	features_df = pd.DataFrame.from_dict(feature_df_dict_list)

	return protein_length_df, protein_occurrence_df, pdb_list_df, coverage_df, db_coverage_df, features_df, protein_seq_str


external_stylesheets = ['https://cdn.jsdelivr.net/npm/bootswatch@4.5.2/dist/darkly/bootstrap.min.css']

storage_type='local'

# Create Figure folder if it doesn't exist
if not os.path.exists('Figures') :
	print('Creating Figures folder')
	os.mkdir('Figures')

# Read most represented features xlsx file
try :
	print('Trying to open without openpyxl')
	most_repr_feat_shape_df = pd.read_excel('most_repr_feature_shape.xlsx', index_col=0, engine=None) # Uses the default engine
except :
	print('Trying to open with openpyxl')
	most_repr_feat_shape_df = pd.read_excel('most_repr_feature_shape.xlsx', index_col=0, engine='openpyxl') # Uses the outdated openpyxl engine

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)


logo_ratio = (1709/390)
tooltip_font_size = '16px'

app.layout = html.Div([
	# Data storage
	dcc.Store(id='Step_1', storage_type=storage_type),
	dcc.Store(id='file_download_prog_perc', storage_type=storage_type),
	# Step 2 Input Store
	dcc.Store(id='Step_2_exceptions', storage_type=storage_type),
	# Step 2_to_3  Store
	dcc.Store(id='2_to_3_features', storage_type=storage_type),
	dcc.Store(id='2_to_3_protein_seq', storage_type=storage_type),
	dcc.Store(id='2_to_3_disorder', storage_type=storage_type),
	dcc.Store(id='2_to_3_pdb_list', storage_type=storage_type),
	dcc.Store(id='2_to_3_pdb_coverage', storage_type=storage_type),
	dcc.Store(id='2_to_3_pdb_db_coverage', storage_type=storage_type),
	dcc.Store(id='2_to_3_pdb_threshold_list', storage_type=storage_type),
	dcc.Store(id='obsolete_list', storage_type=storage_type),
	dcc.Store(id='unrecognized_list', storage_type=storage_type),

	dcc.Store(id='feature_extract', storage_type=storage_type),
	dcc.Store(id='feature_regex_data', storage_type=storage_type),

	# Step 2_to_4 Store
	dcc.Store(id='2_to_4_protein_length', storage_type=storage_type),
	# Step 3 Store
	dcc.Store(id='Step_3_db_PDB', storage_type=storage_type),
	dcc.Store(id='Step_3_db_user', storage_type=storage_type),

	dcc.Store(id='Step_3_feature_shape', storage_type=storage_type),
	dcc.Store(id='Step_3_protein_cut', storage_type=storage_type),

	dcc.Store(id='Step_4_pdb_coverage_list', storage_type=storage_type),
	dcc.Store(id='Step_4_pdb_coverage_steps', storage_type=storage_type),
	dcc.Store(id='Step_4_parameters', storage_type=storage_type),

	dcc.Store(id='Sorted_protein_list', storage_type=storage_type),
	dcc.Store(id='Figure_code', storage_type=storage_type),
	dcc.Store(id='Legend_code', storage_type=storage_type),

	dcc.Store(id='QR_Figure_code', storage_type=storage_type),
	dcc.Store(id='QR_Legend_code', storage_type=storage_type),

	# Step 4 Store
	dcc.Store(id='final_store', storage_type=storage_type),
	dcc.Store(id='Step_4_feature_occurrence', storage_type=storage_type),

	dcc.Store(id='last_seed', storage_type=storage_type),
	dcc.Store(id='last_help', storage_type='memory'),

	html.Div(id='intensity_store_div'),
	html.Div(id='placeholder_div'),


	
	dbc.Row(
		[
			html.Div(children=[
				html.Div(html.Img(src=b64_image(os.path.join('GUI', 'ProFeatMap_invert'+'.png')), style={'width':str(int(logo_ratio*10))+'vh', 'height':'10vh'}), style={'position':'absolute', 'margin-top':'2vh', 'margin-left':'4vh'}),
				html.Div([
						html.Img(src=b64_image(os.path.join('GUI', 'question'+'.png')), style={'width':'5vh', 'height':'5vh', 'margin-top':'4vh'})
					],
					style={'position':'absolute', 'left':'61vw'}
				),
				html.Div([
						html.H3('User guide', style={'position':'absolute', 'width':'12vw', 'background-color':'rgb(223, 170, 24)', 'color':'black', 'font-family':'courier', 'font-weight':'bold', 'text-decoration':'underline', 'margin-top':'5vh', 'text-align':'center'}),
						dbc.Button(id='user_guide_btn', style={'position':'absolute', 'background-color':'transparent', 'width':'12vw', 'height':'4vh', 'margin-top':'5vh'}),
						Download(id='user_guide_download'),
					],
					style={'position':'absolute', 'left':'65vw'}
				),
				html.Div([
						html.H4('ProFeatMap v1.0.0', style={'width':'20vw', 'text-align':'left', 'margin-top':'2vh'}),
						html.H4('For bug report and feedback:', style={'width':'30vw', 'text-align':'left', 'margin-top':'1vh'}),
						html.H4('profeatmap@gmail.com', style={'width':'30vw', 'text-align':'left', 'margin-top':'0vh'}),
					],
					style={'position':'absolute', 'left':'80vw'}
				),

				],
				style={'position':'relative', 'height':'10vh'},
			),
		]
	),
	html.Hr(),
	# Quick run
	dbc.Row(
		[
			html.Div(children=[
				html.Div(html.Img(src=b64_image(os.path.join('GUI', 'red'+'.png')), style={'width':'100vw', 'height':'8vh'}), style={'position':'absolute', 'left':0}),
				html.Div(
					id='quick_run_btn',
					children=[
						html.H2('Quick run guide: Click on me!', style={'text-align':'center', 'width':'100vw', 'margin-top':'1vh'}),
					],
					style={'position':'absolute', 'left':0}
				),
				],
				style={'position':'relative', 'height':'8vh'}
			),
		]
	),
	html.Hr(),
	dbc.Collapse(
		id='quick_guide_collapse',
		children=[
			dbc.Row(
				[
					dbc.Col(
						[
							dbc.Row([
								dbc.Col(
									html.Div(html.Img(id='quick_run_btn_2', src=b64_image(os.path.join('GUI', 'red'+'.png')), style={'width':'100%', 'height':'35vh'})),
									width={'size': 6},
								),
								]
							)
						],
						width={'size': 1},
					),
					dbc.Col(
						[
							dbc.Row([
								dbc.Col([
									html.H4('1', style={'text-align':'right'}),
								],
								width={'size': 1},
								),
								dbc.Col([
									html.H4('Upload a protein list in Step 1 section in the Drag and Drop zone. (A protein “List example” file can be downloaded if needed)', style={'text-align':'left'}),
								],
								width={'size': 9},
								),
							]),
							dbc.Row([
								dbc.Col([
									html.H5('', style={'text-align':'left', 'margin-left':'5vh'}),
								],
								width={'size': 10, 'offset': 1},
								),
							]),
							dbc.Row([
								dbc.Col([
									html.H4('2', style={'text-align':'right'}),
								],
								width={'size': 1},
								),
								dbc.Col([
									html.H4('Click on the Normal run button and wait until finished.', style={'text-align':'left'}),
								],
								width={'size': 9},
								),
							]),
							dbc.Row([
								dbc.Col([
									html.H5('', style={'text-align':'left', 'margin-left':'5vh'}),
								],
								width={'size':10, 'offset': 1},
								),
							]),
							dbc.Row([
								dbc.Col([
									html.H4('3', style={'text-align':'right'}),
								],
								width={'size': 1},
								),
								dbc.Col([
									html.H4('Go to Step 2 section and click the Run button.', style={'text-align':'left'}),
								],
								width={'size': 9},
								),
							]),
							dbc.Row([
								dbc.Col([
									html.H5('', style={'text-align':'left', 'margin-left':'5vh'}),
								],
								width={'size': 10, 'offset': 1},
								),
							]),
							dbc.Row([
								dbc.Col([
									html.H4('4', style={'text-align':'right'}),
								],
								width={'size': 1},
								),
								dbc.Col([
									html.H4('Skip Step 3 and go directly to Step 4. Under the Shapes and colors table, click on the AUTOMATIC feature selection and save the changes.', style={'text-align':'left'}),
								],
								width={'size': 9},
								),
							]),
							dbc.Row([
								dbc.Col([
									html.H5('', style={'text-align':'left', 'margin-left':'5vh'}),
								],
								width={'size': 10, 'offset': 1},
								),
							]),
							dbc.Row([
								dbc.Col([
									html.H4('5', style={'text-align':'right'}),
								],
								width={'size': 1},
								),
								dbc.Col([
									html.H4('Click the Run button in Step 4. The map will be shown once generated.', style={'text-align':'left'}),
								],
								width={'size': 9},
								),
							]),
							dbc.Row([
								dbc.Col([
									html.H5('', style={'text-align':'left', 'margin-left':'5vh'}),
								],
								width={'size': 10, 'offset': 1},
								),
							]),
							dbc.Row([
								dbc.Col([
									html.H4('Map Only', style={'text-align':'right', 'color':'rgb(223, 170, 24)'}),
								],
								width={'size': 1},
								),
								dbc.Col([
									html.H4('After uploading a list of protein, click the Map only button.', style={'text-align':'left'}),
								],
								width={'size': 9},
								),
							]),
							dbc.Row(
								style={'height':'4vh'}
							),
							dbc.Row([
								dbc.Col([
									html.H5('Tip: You can close the Quick run guide by clicking again on the red banner.', style={'text-align':'left'}),
								],
								width={'size': 10, 'offset': 1},
								),
							]),
						],
						width={'size': 10},
					),
				]
			),
			html.Hr(),
		],
	),
	# Step 1
	dbc.Row(
		[
			html.Div(
				id='step_1_collapse_btn',
				children=[
				html.Div(html.Img(src=b64_image(os.path.join('GUI', 'lavender'+'.png')), style={'width':'100vw', 'height':'8vh'}), style={'position':'absolute', 'left':0}),
				html.Div([
					html.H2('Step 1: Protein data gathering', style={'text-align':'center', 'width':'100vw', 'margin-top':'1vh'}),
					],
					style={'position':'absolute', 'left':0}
				),
				],
				style={'position':'relative', 'height':'8vh'}
			),
		]
	),
	dbc.Collapse(
		id='step_1_collapse',
		children=[
			html.Hr(),
			dbc.Row(
				[
					dbc.Col(
						[
							dbc.Row(
								dbc.Col(
									html.Div(html.Img(src=b64_image(os.path.join('GUI', 'lavender'+'.png')), style={'width':'100%', 'height':'8vh'})),
									width={'size': 6},
								),
							)
						],
						width={'size': 1},
					),
					dbc.Col([
							html.H4('This first step needs a list of proteins, with their Uniprot accession code and a unique protein name.'),
							html.H4('ProFeatMap will retrieve data from Uniprot based on these codes.'),
							html.H4('The name given for the proteins will be used as identifier and will appear on the final map.'),
						],
						width={'size': 10},
					),
				]
			),
			html.Hr(),
			dbc.Row(
				[
					dbc.Col(
						[
							dbc.Row(
								dbc.Col(
									html.Div(html.Img(src=b64_image(os.path.join('GUI', 'lavender'+'.png')), style={'width':'100%', 'height':'8vh'})),
									width={'size': 6},
								),
							)
						],
						width={'size': 1},
					),
					dbc.Col(
						[
							html.H3('Protein list', style={'text-align':'center'}),
						],
						width={'size': 10},
					),
				]
			),
			dbc.Row(
				[
					dbc.Col(
						[
							dbc.Row(
								dbc.Col(
									html.Div(html.Img(id='step_1_height', src=b64_image(os.path.join('GUI', 'lavender'+'.png')), style={'width':'100%', 'height':'28vh'})),
									width={'size': 6},
								),
							)
						],
						width={'size': 1},
					),
					dbc.Col([
						dbc.Row(
							[
								html.H4('Uniprot website', style={'margin-top':'1vh'}),
								dbc.Col(
									children=[
										html.Div(
											children=[
												html.Div(html.Img(src=b64_image(os.path.join('GUI', 'link'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
												html.Div([
													html.A(
														dbc.Button(style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}),
														href='https://www.uniprot.org/',
														target='_blank'
													)
													],
													style={'position':'absolute'}
												),
											],
											style={'position':'relative'})
										],
									width={'size': 2},
								),
								dbc.Col(
									children=[
										html.Div(
											id='uniprot_start_state',
											children=[
												html.Img(id='uniprot_end_state', src=b64_image(os.path.join('GUI', 'fail'+'.png')), style={'width':'4vh', 'height':'4vh'}),
												dbc.Tooltip(
													'If cross on red background, Uniprot can temporarly not be reached by ProFeatMap. Please retry later. Reload page to update status.',
													target='uniprot_end_state',
													placement='bottom',
													style={'font-size':tooltip_font_size}
												)
											],
										),
									],
									width={'size': 1},
								),
							],
						),
						dbc.Row(
							style={'height':'2vh'}
						),
						dbc.Row(
							[
								html.H4(id='entry_nb', children='Proteins:'),
							],
						),
						dbc.Row(
							[
								html.H4(id='unique_entry_nb', children='Unique proteins:'),
							],
						),
						dbc.Row(
							[
								html.H4(id='unique_organism'),
							],
						),
						dbc.Row(
							[
								html.Div(
									id='duplicate_entry_names_div',
								),
							],
						),
						],
						width={'size': 2}
					),
					dbc.Col([
						dbc.Row(
							style={'height':'2vh'}
						),
						dbc.Row(
							[
								dbc.Col([
									html.Div(
										id='protein_list_dnd_div',
										children=dcc.Upload(
											id='protein_list_dnd',
											children=html.Div([
												'Drag and Drop or ',
												html.A('Select Files'),
												dbc.Tooltip(
													'Click to select file or drag and drop file directly in box (.xlsx, .ods or .csv)',
													target='protein_list_dnd_div',
													placement='bottom',
													style={'font-size':tooltip_font_size}
												)
											]),
											style={'height': '60px','lineHeight': '60px','borderWidth': '1px','borderStyle': 'dashed','borderRadius': '5px','textAlign': 'center','margin': '10px'},
										),
									),
								],
								),
							]
						),
						dbc.Row(
							[
								dbc.Col([
									html.Div(id='protein_list_datatable_div', children=dash_table.DataTable(id='protein_list_datatable')),
								],
								),
							]
						),
						dbc.Row(
							style={'height':'1vh'}
						),
						dbc.Row(
							[
								dbc.Col(
									id='step_1_add_row_div',
									children=[
										dbc.Row(
											[
												dbc.Col([
													html.Div(
														children=[
															html.Div(html.Img(src=b64_image(os.path.join('GUI', 'add'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															html.Div(dbc.Button(id='protein_list_editing_rows_btn', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															dbc.Tooltip(
																	'Add a row to the datatable',
																	target='step_1_add_row_div',
																	placement='bottom',
																	style={'font-size':tooltip_font_size}
															)
														],
														style={'position':'relative'})
													],
												),
											]
										)
									],
									width={'size': 1, 'offset':1}
								),
								dbc.Col(
									id='step_1_clear_div',
									children=[
										dbc.Row(
											[
												dbc.Col([
													html.Div(
														children=[
															html.Div(html.Img(src=b64_image(os.path.join('GUI', 'clear'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															html.Div(dbc.Button(id='protein_list_clear', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															dbc.Tooltip(
																	'Clear datatable content',
																	target='step_1_clear_div',
																	placement='bottom',
																	style={'font-size':tooltip_font_size}
															)
														],
														style={'position':'relative'})
													],
												),
											]
										)
									],
									width={'size': 1, 'offset':2}
								),
								dbc.Col(
									id='step_1_save_div',
									children=[
										dbc.Row(
											[
												dbc.Col([
													html.Div(
														children=[
															html.Div(id='step_1_save_img', children=[html.Img(src=b64_image(os.path.join('GUI', 'save'+'.png')), style={'width':'4vh', 'height':'4vh'})], style={'position':'absolute'}),
															html.Div(dbc.Button(id='protein_list_changes_save', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															dbc.Tooltip(
																	'Save datatable changes',
																	target='step_1_save_div',
																	placement='bottom',
																	style={'font-size':tooltip_font_size}
															)
														],
														style={'position':'relative'})
													],
												),
											]
										)
									],
									width={'size': 1, 'offset':2}
								),
								dbc.Col(
									id='step_1_download_div',
									children=[
										dbc.Row(
											[
												dbc.Col([
													html.Div(
														children=[
															html.Div(html.Img(src=b64_image(os.path.join('GUI', 'load_white'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															html.Div(dbc.Button(id='protein_list_download_button', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															dbc.Tooltip(
																	'Download datatable content as excel file',
																	target='step_1_download_div',
																	placement='bottom',
																	style={'font-size':tooltip_font_size}
															)
														],
														style={'position':'relative', 'height':'4vh'})
													],
												),
											]
										)
									],
									width={'size': 1, 'offset':2}
								),
								Download(id='protein_list_download'),
							]
						),
						],
						width={'size': 6}
					),

					dbc.Col(
						children=[
							dbc.Row(
								id='step_1_options1_div',
								children=[
									html.H4('List example', style={'margin-top':'1vh'}),
									dbc.Col(
										children=[
											html.Div(
												children=[
													html.Div(html.Img(src=b64_image(os.path.join('GUI', 'load_white'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
													html.Div(dbc.Button(id='protein_list_example_download_button', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
													dbc.Tooltip(
															'Download an example of protein list with TRIO_HUMAN',
															target='protein_list_example_download_button',
															placement='bottom',
															style={'font-size':tooltip_font_size}
													),
													Download(id='protein_list_example_download'),
												],
												# style={'position':'relative'}
												)
											],
											width={'size': 2},
									),
									dbc.Col(
										children=[
											html.Div(
												children=[
													html.Div(html.Img(src=b64_image(os.path.join('GUI', 'load_white'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
													html.Div(dbc.Button(id='protein_list_example_download_button_2', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
													dbc.Tooltip(
															'Download an example of protein list of PDZ domain containing proteins',
															target='protein_list_example_download_button_2',
															placement='bottom',
															style={'font-size':tooltip_font_size}
													),
													Download(id='protein_list_example_download_2'),
												],
												# style={'position':'relative'}
												)
											],
											width={'size': 2},
									),
								]
							),
							dbc.Row(
								style={'height':'2vh'}
							),							
							dbc.Row(
								id='step_1_options2_div',
								children=[
									html.H4('Remove organism', style={'margin-top':'1vh'}),
									dbc.Col(
										children=[
											html.Div(
												children=[
													html.Div(html.Img(src=b64_image(os.path.join('GUI', 'option_run'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
													html.Div(dbc.Button(id='remove_organism_button', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
													dbc.Tooltip(
															'Remove organism tag of current protein list',
															target='step_1_options2_div',
															placement='bottom',
															style={'font-size':tooltip_font_size}
														)
												],
												style={'position':'relative'})
											],
									),
								]
							),							
						],
						width={'size': 3}
					),
				]
			),
			dbc.Row(
				style={'height':'2vh'}
			),	
			dbc.Row(
				[
					dbc.Col(
						[
							dbc.Row(
								dbc.Col(
									html.Div(html.Img(src=b64_image(os.path.join('GUI', 'lavender'+'.png')), style={'width':'100%', 'height':'6vh'})),
									width={'size': 6},
								),
							)
						],
						width={'size': 1},
					),
					dbc.Col([
						html.H4(id='step_1_dowload_latest', children='Download latest data', style={'text-align':'right'}),
						dbc.Tooltip(
							'Select to force download of latest available data from the Uniprot database',
							target='step_1_dowload_latest',
							placement='bottom',
							style={'font-size':tooltip_font_size}
						),
						],
						width={'size': 2},
					),
					dbc.Col([
						html.Div(
							children=[
								dbc.Checklist(
									id='dl_latest_toggle',
									options=[
										{'label': '', 'value': 'latest'},
									],
									value=[],
								),
							],
						),
						],
						width={'size': 1},
					),
					dbc.Col(
						id='UQR_button_div',
						children=[
							dbc.Row(
								[
									html.H3('Map only', style={'text-align':'center', 'margin-top':'0vh'}),
									dbc.Col([
										html.Div(
											children=[
												html.Div(html.Img(src=b64_image(os.path.join('GUI', 'quick_run'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
												html.Div(dbc.Button(id='UQR_button', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
												dbc.Tooltip(
														'Use default parameters except for map parameters to generate map directly',
														target='UQR_button_div',
														placement='bottom',
														style={'font-size':tooltip_font_size}
												)
											],
											style={'position':'relative'})
										],
									),
								]
							),
							dcc.Loading(
								id='step_1_UQR_loading',
								children=[html.Div(id='step_UQR_placeholder')],
								type='default',
								color='#DFAA18',
							),
						],
						width={'size': 2}
					),
					dbc.Col(
						id='step_1_run_div',
						children=[
						dbc.Row(
							[
								html.H3('Normal run', style={'text-align':'center', 'margin-top':'0vh'}),
								dbc.Col([
									html.Div(
										children=[
											html.Div(html.Img(src=b64_image(os.path.join('GUI', 'run'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
											html.Div(dbc.Button(id='step_1_launch_button', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
											dbc.Tooltip(
													'Start protein data gathering',
													target='step_1_run_div',
													placement='bottom',
													style={'font-size':tooltip_font_size}
											)
										],
										style={'position':'relative'})
									],
								),
							]
						),
						dcc.Loading(
							id='step_1_loading',
							children=[html.Div(id='step_1_placeholder')],
							type='default',
							color='#94A3BC',
						),
						],
						width={'size': 2}
					),
					dbc.Col(
						children=[
							html.Div(id='step_1_start_state', children=html.Img(id='step_1_end_state', src=b64_image(os.path.join('GUI', 'neutral'+'.png')), style={'width':'4vh', 'height':'4vh'})),
						],
						width={'size': 1},
					),
				]
			),
			dbc.Row(
				id='obsolete_code_div',
				children=
				[
					dbc.Col([
							html.Img(src=b64_image(os.path.join('GUI', 'warning'+'.png')), style={'width':'6vh', 'height':'6vh'})
						],
						width={'size': 1, 'offset':2}
					),
					dbc.Col([
							html.H4(children='Some Uniprot codes are found to be obsolete.', style={'color': '#a50104', 'margin-top':'1vh'}),
							html.H4(children='Please remove or replace them, then rerun Step 1.', style={'color': '#a50104'}),
							html.H4(children='Download the concerned list with the button below.')
						],
						width={'size': 6, 'offset':3}
					),
					dbc.Col([
						html.Div(
							children=[
								html.Div(children=[html.Img(src=b64_image(os.path.join('GUI', 'load_white'+'.png')), style={'width':'4vh', 'height':'4vh'})], style={'position':'absolute'}),
								html.Div(dbc.Button(id='obsolete_code_download_btn', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
								dbc.Tooltip(
									'Download list of proteins that ProFeatMap could not download',
									target='obsolete_code_download_btn',
									placement='bottom',
									style={'font-size':tooltip_font_size}
								),
								Download(id='obsolete_code_download'),
							],
							style={'position':'relative', 'height':'4vh'})
						],
						width={'size': 2, 'offset':5}
					),
				],
				style={'display': 'none'}
			),
		]
	),
	dbc.Collapse(
		id='quick_run_collapse',
		children=[
			dbc.Col(
				children=[
					dbc.Row(
						[
							dbc.Col([
								html.Div(
									children=[
										html.Div(children=[html.Img(src=b64_image(os.path.join('GUI', 'hide'+'.png')), style={'width':'4vh', 'height':'4vh'})], style={'position':'absolute'}),
										html.Div(dbc.Button(id='hide_quick_run_btn', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
									],
									style={'position':'relative', 'height':'4vh'})
								],
							),
							dbc.Col([
								html.Div(
									children=[
										html.Div(children=[html.Img(src=b64_image(os.path.join('GUI', 'load_white'+'.png')), style={'width':'4vh', 'height':'4vh'})], style={'position':'absolute'}),
										html.Div(dbc.Button(id='quick_run_result_download', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
										dbc.Tooltip(
											'Download Map and Legend',
											target='quick_run_result_download',
											placement='bottom',
											style={'font-size':tooltip_font_size}
										),
										Download(id='QR_figure_download'),
										Download(id='QR_legend_download'),
									],
									style={'position':'relative', 'height':'4vh'})
								],
							),
						]
					)
				],
				width={'size': 1, 'offset': 1}
			),
			html.Hr(),
			dbc.Row(
				[
					dbc.Col([
							html.Div(
								id='UQR_result_img',
							)
						],
						width={'size': 10, 'offset':1}
					),
				]
			),
			html.Hr(),
			dbc.Row(
				[
					dbc.Col([
							html.Div(
								id='UQR_legend_img',
							)
						],
						width={'size': 10, 'offset':1}
					),
				]
			),
		],
	),
	html.Hr(),
	# Step 2
	dbc.Row(
		[
			html.Div(
				id='step_2_collapse_btn',
				children=[
				html.Div(html.Img(src=b64_image(os.path.join('GUI', 'almond'+'.png')), style={'width':'100vw', 'height':'8vh'}), style={'position':'absolute', 'left':0}),
				html.Div([
					html.H2('Step 2: Feature extraction', style={'text-align':'center', 'width':'100vw', 'margin-top':'1vh', 'color':'black'}),
					],
					style={'position':'absolute', 'left':0}
				),
				],
				style={'position':'relative', 'height':'8vh'}
			),
		]
	),
	dbc.Collapse(
		id='step_2_collapse',
		children=[
			html.Hr(),
			dbc.Row(
				[
					dbc.Col(
						[
							dbc.Row(
								dbc.Col(
									html.Div(html.Img(src=b64_image(os.path.join('GUI', 'almond'+'.png')), style={'width':'100%', 'height':'8vh'})),
									width={'size': 6},
								),
							)
						],
						width={'size': 1},
					),
					dbc.Col([
							html.H4('Once data is properly gathered from the Uniprot database, ProFeatMap will proceed to the feature extraction.'),
							html.H4('Results of the extraction can be downloaded as an Excel file, with the Download button.'),
							html.H4('Modification of the extraction can be made by providing an exception file.'),
						],
						width={'size': 10},
					),
				]
			),
			html.Hr(),
			dbc.Row(
				[
					dbc.Col(
						[
							dbc.Row(
								dbc.Col(
									html.Div(html.Img(src=b64_image(os.path.join('GUI', 'almond'+'.png')), style={'width':'100%', 'height':'8vh'})),
									width={'size': 6},
								),
							)
						],
						width={'size': 1},
					),
					dbc.Col(
						[
							html.H3('Modifications (optional)', style={'text-align':'center'}),
						],
						width={'size': 10},
					),
				]
			),
			dbc.Row(
				[
					dbc.Col(
						[
							dbc.Row(
								dbc.Col(
									html.Div(html.Img(id='step_2_height', src=b64_image(os.path.join('GUI', 'almond'+'.png')), style={'width':'100%', 'height':'87vh'})),
									width={'size': 6},
								),
							)
						],
						width={'size': 1},
					),
					dbc.Col([
						],
						width={'size': 2}
					),
					dbc.Col([
						dbc.Row(
							[
								dbc.Col([
									html.Div(
										id='exception_dnd_div',
										children=dcc.Upload(
											id='exception_dnd',
											children=html.Div([
												'Drag and Drop or ',
												html.A('Select Files')
											]),
											style={'height': '60px','lineHeight': '60px','borderWidth': '1px','borderStyle': 'dashed','borderRadius': '5px','textAlign': 'center','margin': '10px'},
										),
									),
								],
								),
							]
						),
						dbc.Row(
							[
								dbc.Col([
									html.Div(id='exception_datatable_div', children=dash_table.DataTable(id='exception_datatable')),
								],
								),
							]
						),
						dbc.Row(
							style={'height':'1vh'}
						),
						dbc.Row(
							[
								dbc.Col(
									id='step_2_add_row_div',
									children=[
										dbc.Row(
											[
												dbc.Col([
													html.Div(
														children=[
															html.Div(html.Img(src=b64_image(os.path.join('GUI', 'add'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															html.Div(dbc.Button(id='exception_editing_rows_btn', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															dbc.Tooltip(
																	'Add a row to the datatable',
																	target='step_2_add_row_div',
																	placement='bottom',
																	style={'font-size':tooltip_font_size}
															)
														],
														style={'position':'relative'})
													],
												),
											]
										)
									],
									width={'size': 1, 'offset':1}
								),
								dbc.Col(
									id='step_2_clear_div',
									children=[
										dbc.Row(
											[
												dbc.Col([
													html.Div(
														children=[
															html.Div(html.Img(src=b64_image(os.path.join('GUI', 'clear'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															html.Div(dbc.Button(id='exception_clear', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															dbc.Tooltip(
																	'Clear datatable content',
																	target='step_2_clear_div',
																	placement='bottom',
																	style={'font-size':tooltip_font_size}
															)
														],
														style={'position':'relative'})
													],
												),
											]
										)
									],
									width={'size': 1, 'offset':2}
								),
								dbc.Col(
									id='step_2_save_div',
									children=[
										dbc.Row(
											[
												dbc.Col([
													html.Div(
														children=[
															html.Div(id='step_2_save_img', children=[html.Img(src=b64_image(os.path.join('GUI', 'save'+'.png')), style={'width':'4vh', 'height':'4vh'})], style={'position':'absolute'}),
															html.Div(dbc.Button(id='exception_changes_save', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															dbc.Tooltip(
																	'Save datatable changes',
																	target='step_2_save_div',
																	placement='bottom',
																	style={'font-size':tooltip_font_size}
															)
														],
														style={'position':'relative'})
													],
												),
											]
										)
									],
									width={'size': 1, 'offset':2}
								),
								dbc.Col(
									id='step_2_download_div',
									children=[
										dbc.Row(
											[
												dbc.Col([
													html.Div(
														children=[
															html.Div(html.Img(src=b64_image(os.path.join('GUI', 'load_white'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															html.Div(dbc.Button(id='exception_download_button', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															dbc.Tooltip(
																	'Download datatable content as excel file',
																	target='step_2_download_div',
																	placement='bottom',
																	style={'font-size':tooltip_font_size}
															)
														],
														style={'position':'relative', 'height':'4vh'})
													],
												),
											]
										)
									],
									width={'size': 1, 'offset':2}
								),
								Download(id='exception_download'),
							]
						),
						],
						width={'size': 6}
					),
					dbc.Col(
						children=[
							dbc.Row(
								id='step_2_options1_div',
								children=[
									html.H4('Modifications example', style={'margin-top':'1vh'}),
									dbc.Col(
										children=[
											html.Div(
												children=[
													html.Div(html.Img(src=b64_image(os.path.join('GUI', 'load_white'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
													html.Div(dbc.Button(id='modifications_example_download_button', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
													dbc.Tooltip(
															'Download an example of modifications for TRIO_HUMAN',
															target='modifications_example_download_button',
															placement='bottom',
															style={'font-size':tooltip_font_size}
													),
													Download(id='modifications_example_download'),
												],
												# style={'position':'relative'}
												)
											],
											width={'size': 2},
									),
									dbc.Col(
										children=[
											html.Div(
												children=[
													html.Div(html.Img(src=b64_image(os.path.join('GUI', 'load_white'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
													html.Div(dbc.Button(id='modifications_example_download_button_2', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
													dbc.Tooltip(
															'Download an example of modifications for PDZ domain containing proteins',
															target='modifications_example_download_button_2',
															placement='bottom',
															style={'font-size':tooltip_font_size}
													),
													Download(id='modifications_example_download_2'),
												],
												# style={'position':'relative'}
												)
											],
											width={'size': 2},
									),
								]
							),						
						],
						width={'size': 3}
					),
				]
			),
			dbc.Row(
				style={'height':'2vh'}
			),	
			dbc.Row(
				[
					dbc.Col(
						[
							dbc.Row(
								dbc.Col(
									html.Div(html.Img(src=b64_image(os.path.join('GUI', 'almond'+'.png')), style={'width':'100%', 'height':'8vh'})),
									width={'size': 6},
								),
							)
						],
						width={'size': 1},
					),
					dbc.Col(
						[
							dbc.Row(
								dbc.Col(
									html.H3('Structural coverage extraction', style={'text-align':'center'}),
								),
							)
						],
						width={'size': 10},
					),
				]
			),
			dbc.Row(
				[
					dbc.Col(
						[
							dbc.Row(
								dbc.Col(
									html.Div(html.Img(src=b64_image(os.path.join('GUI', 'almond'+'.png')), style={'width':'100%', 'height':'10vh'})),
									width={'size': 6},
								),
							)
						],
						width={'size': 1},
					),
					dbc.Col([
						html.H4('Colormap thresholds', style={'text-align':'right'}),
					],
					width={'size': 2}
					),
					dbc.Col([
						dcc.RangeSlider(
							id='pdb_threshold_list',
							min=0,
							max=50,
							value=[1,2,3,5,10],
							marks={1: '1', 2: '2', 3: '3', 5: '5', 10: '10', 15: '15', 20: '20', 25: '25', 30: '30', 50: '50'},
							pushable=1,
						)
					],
					width={'size': 5}
					),
					dbc.Col([
						html.H4('Step number', style={'text-align':'right'}),
					],
					width={'size': 1}
					),
					dbc.Col([
						dcc.Input(
							id='nb_steps_pdb_threshold',
							type='number',
							min=1,
							max=10,
							value=5,
							style={'width':'100%'}
						)
					],
					width={'size': 1}
					),
					dbc.Col(
						children=[
							dbc.Row(
								[
									dbc.Col([
										html.Div(
											children=[
												html.Div(children=[html.Img(src=b64_image(os.path.join('GUI', 'default'+'.png')), style={'width':'4vh', 'height':'4vh'})], style={'position':'absolute'}),
												html.Div(dbc.Button(id='pdb_threshold_default_btn', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
												dbc.Tooltip(
														'Reset to default parameters',
														target='pdb_threshold_default_btn',
														placement='bottom',
														style={'font-size':tooltip_font_size}
												)
											],
											style={'position':'relative'})
										],
									),
								]
							)
						],
						width={'size': 1}
					),
				]
			),
			dbc.Row(
				[
					dbc.Col(
						[
							dbc.Row(
								dbc.Col(
									html.Div(html.Img(src=b64_image(os.path.join('GUI', 'almond'+'.png')), style={'width':'100%', 'height':'6vh'})),
									width={'size': 6},
								),
							)
						],
						width={'size': 1},
					),
					dbc.Col([
						html.H4(id='extract_only_txt', children='Extract only', style={'text-align':'right'}),
						Download(id='direct_download'),
						dbc.Tooltip(
								'Download extracted data at the end of the step directly. Activating this option does not allow to generate maps or use feature/sequence extraction tools',
								target='extract_only_txt',
								placement='bottom',
								style={'font-size':tooltip_font_size}
						),
					],
					width={'size': 2},
					),
					dbc.Col([
						dbc.Checklist(
							id='big_data_toggle',
							options=[
								{'label': '', 'value': 'big_data'},
							],
							value=[],
						),
					],
					width={'size': 1},
					),
					dbc.Col(
						id='step_2_run_div',
						children=[
						dbc.Row(
							[
								html.H3('Run', style={'text-align':'center', 'margin-top':'0vh'}),
								dbc.Col([
									html.Div(
										children=[
											html.Div(html.Img(src=b64_image(os.path.join('GUI', 'run'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
											html.Div(dbc.Button(id='step_2_launch_button', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
											dbc.Tooltip(
													'Start feature extraction',
													target='step_2_run_div',
													placement='bottom',
													style={'font-size':tooltip_font_size}
												)
										],
										style={'position':'relative'})
									],
								),
							]
						),
						dcc.Loading(
							id='step_2_loading',
							children=[html.Div(id='step_2_placeholder')],
							type='default',
							color='#94A3BC',
						),
						],
						width={'size': 1, 'offset':2}
					),
					dbc.Col(
						children=[
							html.Div(id='step_2_start_state', children=html.Img(id='step_2_end_state', src=b64_image(os.path.join('GUI', 'neutral'+'.png')), style={'width':'4vh', 'height':'4vh'})),
						],
						width={'size': 1},
					),
				]
			),
			html.Hr(),
			dbc.Collapse(
				id='step_2_extraction_collapse',
				children=[
					dbc.Row(
						[
							dbc.Col(
								[
									dbc.Row(
										dbc.Col(
											html.Div(html.Img(src=b64_image(os.path.join('GUI', 'almond'+'.png')), style={'width':'100%', 'height':'8vh'})),
											width={'size': 6},
										),
									)
								],
								width={'size': 1},
							),
							dbc.Col(
								[
									dbc.Row(
										id='step_2_extraction_download_div',
										children=[
											html.H3('Extracted data'),
											dbc.Col(
												children=[
													dbc.Row(
														[
															dbc.Col([
																html.Div(
																	children=[
																		html.Div(html.Img(src=b64_image(os.path.join('GUI', 'load_white'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
																		html.Div(dbc.Button(id='protein_length_download_btn', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
																		dbc.Tooltip(
																				'Download extracted data',
																				target='step_2_extraction_download_div',
																				placement='bottom',
																				style={'font-size':tooltip_font_size}
																		),
																		Download(id='protein_length_download'),
																	],
																	style={'position':'relative', 'height':'4vh'})
																],
															),
														]
													),
												],
												width={'size': 1}
											),
										],
									)
								],
								width={'size': 2, 'offset':4},
							),
						]
					),
					dbc.Row(
						[
							dbc.Col(
								[
									dbc.Row(
										dbc.Col(
											html.Div(html.Img(src=b64_image(os.path.join('GUI', 'almond'+'.png')), style={'width':'100%', 'height':'8vh'})),
											width={'size': 6},
										),
									)
								],
								width={'size': 1},
							),
							dbc.Col(
								[
									dbc.Row(
										dbc.Col(
											html.H3('Feature sequence extraction parameters', style={'text-align':'center'}),
										),
									)
								],
								width={'size': 10},
							),
						]
					),
					dbc.Row(
						[
							dbc.Col([
								html.H5('Feature name', style={'text-align':'right'}),
							],
							width={'size': 2, 'offset': 1}
							),
							dbc.Col([
								dcc.Input(
									id='feature_extract_name',
									style={'width':'100%'}
								)
							],
							width={'size': 2}
							),
							dbc.Col([
								html.H5('N-ter ext', style={'text-align':'right'}),
							],
							width={'size': 1}
							),
							dbc.Col([
								dcc.Input(
									id='feature_extract_nter_ext',
									value='0',
									style={'width':'100%'}
								)
							],
							width={'size': 1}
							),
							dbc.Col([
								html.H5('C-ter ext', style={'text-align':'right'}),
							],
							width={'size': 1}
							),
							dbc.Col([
								dcc.Input(
									id='feature_extract_cter_ext',
									value='0',
									style={'width':'100%'}
								)
							],
							width={'size': 1}
							),
							dbc.Col(
								children=[
									html.Div(
										children=[
											html.Div(html.Img(src=b64_image(os.path.join('GUI', 'option_run'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
											html.Div(dbc.Button(id='feature_extract_launch_button', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
											dbc.Tooltip(
													'Extract sequences of all features found in your protein list',
													target='feature_extract_launch_button',
													placement='bottom',
													style={'font-size':tooltip_font_size}
											),
											dcc.Loading(
												id='feature_extract_loading',
												children=[html.Div(id='feature_extract_placeholder')],
												type='default',
												color='#DFAA18',
											),
										],
										style={'position':'relative'})
									],
								width={'size': 1}
							),
							dbc.Col(
								children=[
									html.Div(
										children=[			
											html.Div(html.Img(src=b64_image(os.path.join('GUI', 'load_white'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
											html.Div(dbc.Button(id='feature_extract_download_button', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
											dbc.Tooltip(
													'Download extracted sequences',
													target='feature_extract_download_button',
													placement='bottom',
													style={'font-size':tooltip_font_size}
											),
											Download(id='feature_extract_download'),
										],
										style={'position':'relative'})
									],
								width={'size': 1}
							),
						]
					),



					dbc.Row(
						[
							dbc.Col(
								[
									dbc.Row(
										dbc.Col(
											html.Div(html.Img(src=b64_image(os.path.join('GUI', 'almond'+'.png')), style={'width':'100%', 'height':'8vh'})),
											width={'size': 6},
										),
									)
								],
								width={'size': 1},
							),
							dbc.Col(
								[
									dbc.Row(
										dbc.Col(
											html.H3('Feature/Motif search by regular expression', style={'text-align':'center'}),
										),
									)
								],
								width={'size': 10},
							),
						]
					),
					dbc.Row(
						[
							dbc.Col([
								html.H5('Feature regular expression', style={'text-align':'right'}),
							],
							width={'size': 2, 'offset': 1}
							),
							dbc.Col([
								dcc.Input(
									id='feature_regex',
									style={'width':'100%'}
								)
							],
							width={'size': 2}
							),
							dbc.Col(
								children=[
									html.Div(
										children=[
											html.Div(html.Img(src=b64_image(os.path.join('GUI', 'link'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
											html.Div([
												html.A(
													dbc.Button(id='regex_link', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}),
													href='https://pythex.org/',
													target='_blank'
												),
												dbc.Tooltip(
														'Create and test your regular expression here',
														target='regex_link',
														placement='bottom',
														style={'font-size':tooltip_font_size}
												),
												],
												style={'position':'absolute'}
											),
										],
										style={'position':'relative'})
									],
								width={'size': 1}
							),
							dbc.Col([
								html.H5('Feature name', style={'text-align':'right'}),
							],
							width={'size': 1}
							),
							dbc.Col([
								dcc.Input(
									id='feature_regex_name',
									style={'width':'100%'}
								)
							],
							width={'size': 2}
							),

							dbc.Col(
								children=[
									html.Div(
										children=[
											html.Div(html.Img(src=b64_image(os.path.join('GUI', 'option_run'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
											html.Div(dbc.Button(id='feature_regex_launch_button', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
											dbc.Tooltip(
													'Extract domain position matching the regular expression',
													target='feature_regex_launch_button',
													placement='bottom',
													style={'font-size':tooltip_font_size}
											),
											dcc.Loading(
												id='feature_regex_loading',
												children=[html.Div(id='feature_regex_placeholder')],
												type='default',
												color='#DFAA18',
											),
										],
										style={'position':'relative'})
									],
								width={'size': 1}
							),
							dbc.Col(
								children=[
									html.Div(
										children=[			
											html.Div(html.Img(src=b64_image(os.path.join('GUI', 'load_white'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
											html.Div(dbc.Button(id='feature_regex_download_button', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
											dbc.Tooltip(
													'Download position of matching regular expressions as modification file format',
													target='feature_regex_download_button',
													placement='bottom',
													style={'font-size':tooltip_font_size}
											),
											Download(id='feature_regex_download'),
										],
										style={'position':'relative'})
									],
								width={'size': 1}
							),
						]
					),
				]
			),
		],
	),
	html.Hr(),
	# Step 3
	dbc.Row(
		[
			html.Div(
				id='step_3_collapse_btn',
				children=[
				html.Div(html.Img(src=b64_image(os.path.join('GUI', 'gold'+'.png')), style={'width':'100vw', 'height':'8vh'}), style={'position':'absolute', 'left':0}),
				html.Div([
					html.H2('Step 3: Numerical values addition', style={'text-align':'center', 'width':'100vw', 'margin-top':'1vh', 'color':'black'}),
					],
					style={'position':'absolute', 'left':0}
				),
				],
				style={'position':'relative', 'height':'8vh'}
			),
		]
	),
	dbc.Collapse(
		id='step_3_collapse',
		children=[
			html.Hr(),
			dbc.Row(
				[
					dbc.Col(
						[
							dbc.Row(
								dbc.Col(
									html.Div(html.Img(src=b64_image(os.path.join('GUI', 'gold'+'.png')), style={'width':'100%', 'height':'8vh'})),
									width={'size': 6},
								),
							)
						],
						width={'size': 1},
					),
					dbc.Col([
							html.H3('This step is optional.'),
							html.H4('Unless you have numerical values you want to display on one or more feature.'),
							html.H4('If so, you will have to provide a file containing these values.'),
						],
						width={'size': 10},
					),
				]
			),
			html.Hr(),
			dbc.Row(
				[
					dbc.Col(
						[
							dbc.Row(
								dbc.Col(
									html.Div(html.Img(src=b64_image(os.path.join('GUI', 'gold'+'.png')), style={'width':'100%', 'height':'8vh'})),
									width={'size': 6},
								),
							)
						],
						width={'size': 1},
					),
					dbc.Col(
						[
							html.H3('Numerical values (optional)', style={'text-align':'center'}),
						],
						width={'size': 10},
					),
				]
			),

			dbc.Row(
				[
					dbc.Col(
						[
							dbc.Row(
								dbc.Col(
									html.Div(html.Img(id='step_3_height', src=b64_image(os.path.join('GUI', 'gold'+'.png')), style={'width':'100%', 'height':'87vh'})),
									width={'size': 6},
								),
							)
						],
						width={'size': 1},
					),
					dbc.Col([
						],
						width={'size': 2}
					),
					dbc.Col([
						dbc.Row(
							[
								dbc.Col([
									html.Div(
										id='values_dnd_div',
										children=dcc.Upload(
											id='values_dnd',
											children=html.Div([
												'Drag and Drop or ',
												html.A('Select Files')
											]),
											style={'height': '60px','lineHeight': '60px','borderWidth': '1px','borderStyle': 'dashed','borderRadius': '5px','textAlign': 'center','margin': '10px'},
										),
									),
								],
								),
							]
						),
						dbc.Row(
							[
								dbc.Col([
									html.Div(id='values_datatable_div', children=dash_table.DataTable(id='values_datatable')),
								],
								),
							]
						),
						dbc.Row(
							style={'height':'1vh'}
						),
						#############
						dbc.Row(
							[
								dbc.Col(
									id='step_3_add_row_div',
									children=[
										dbc.Row(
											[
												dbc.Col([
													html.Div(
														children=[
															html.Div(html.Img(src=b64_image(os.path.join('GUI', 'add'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															html.Div(dbc.Button(id='values_editing_rows_btn', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															dbc.Tooltip(
																	'Add a row to the datatable',
																	target='step_3_add_row_div',
																	placement='bottom',
																	style={'font-size':tooltip_font_size}
															)
														],
														style={'position':'relative'})
													],
												),
											]
										)
									],
									width={'size': 1, 'offset':1}
								),
								dbc.Col(
									id='step_3_clear_div',
									children=[
										dbc.Row(
											[
												dbc.Col([
													html.Div(
														children=[
															html.Div(html.Img(src=b64_image(os.path.join('GUI', 'clear'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															html.Div(dbc.Button(id='values_clear', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															dbc.Tooltip(
																	'Clear datatable content',
																	target='step_3_clear_div',
																	placement='bottom',
																	style={'font-size':tooltip_font_size}
															)
														],
														style={'position':'relative'})
													],
												),
											]
										)
									],
									width={'size': 1, 'offset':2}
								),
								dbc.Col(
									id='step_3_save_div',
									children=[
										dbc.Row(
											[
												dbc.Col([
													html.Div(
														children=[
															html.Div(id='step_3_save_img', children=[html.Img(src=b64_image(os.path.join('GUI', 'save'+'.png')), style={'width':'4vh', 'height':'4vh'})], style={'position':'absolute'}),
															html.Div(dbc.Button(id='values_changes_save', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															dbc.Tooltip(
																	'Save datatable changes',
																	target='step_3_save_div',
																	placement='bottom',
																	style={'font-size':tooltip_font_size}
															)
														],
														style={'position':'relative'})
													],
												),
											]
										)
									],
									width={'size': 1, 'offset':2}
								),
								dbc.Col(
									id='step_3_download_div',
									children=[
										dbc.Row(
											[
												dbc.Col([
													html.Div(
														children=[
															html.Div(html.Img(src=b64_image(os.path.join('GUI', 'load_white'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															html.Div(dbc.Button(id='values_download_button', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															dbc.Tooltip(
																	'Download datatable content as excel file',
																	target='step_3_download_div',
																	placement='bottom',
																	style={'font-size':tooltip_font_size}
															)
														],
														style={'position':'relative', 'height':'4vh'})
													],
												),
											]
										)
									],
									width={'size': 1, 'offset':2}
								),
								Download(id='values_download'),
							]
						),
						],
						width={'size': 6}
					),
					dbc.Col(
						children=[
							dbc.Row(
								id='step_3_options1_div',
								children=[
									html.H4('Numerical values example', style={'margin-top':'1vh'}),
									dbc.Col(
										children=[
											html.Div(
												children=[
													html.Div(html.Img(src=b64_image(os.path.join('GUI', 'load_white'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
													html.Div(dbc.Button(id='values_example_download_button', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
													dbc.Tooltip(
															'Download an example of numerical values file for TRIO_HUMAN',
															target='values_example_download_button',
															placement='bottom',
															style={'font-size':tooltip_font_size}
													),
													Download(id='values_example_download'),
												],
												# style={'position':'relative'}
												)
											],
											width={'size': 2},
									),
									dbc.Col(
										children=[
											html.Div(
												children=[
													html.Div(html.Img(src=b64_image(os.path.join('GUI', 'load_white'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
													html.Div(dbc.Button(id='values_example_download_button_2', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
													dbc.Tooltip(
															'Download an example of numerical values file for PDZ domain containing proteins',
															target='values_example_download_button_2',
															placement='bottom',
															style={'font-size':tooltip_font_size}
													),
													Download(id='values_example_download_2'),
												],
												# style={'position':'relative'}
												)
											],
											width={'size': 2},
									),
								]
							),						
						],
						width={'size': 3}
					),
				]
			),
			dbc.Row(
				style={'height':'2vh'}
			),			
			dbc.Row(
				[
					dbc.Col(
						[
							dbc.Row(
								dbc.Col(
									html.Div(html.Img(src=b64_image(os.path.join('GUI', 'gold'+'.png')), style={'width':'100%', 'height':'6vh'})),
									width={'size': 6},
								),
							)
						],
						width={'size': 1},
					),
					dbc.Col(
						id='step_3_run_div',
						children=[
						dbc.Row(
							[
								html.H3('Run', style={'text-align':'center', 'margin-top':'0vh'}),
								dbc.Col([
									html.Div(
										children=[
											html.Div(html.Img(src=b64_image(os.path.join('GUI', 'run'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
											html.Div(dbc.Button(id='step_3_launch_button', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
											dbc.Tooltip(
													'Start addition of values (optional)',
													target='step_3_run_div',
													placement='bottom',
													style={'font-size':tooltip_font_size}
												)
										],
										style={'position':'relative'})
									],
								),
							]
						),
						dcc.Loading(
							id='step_3_loading',
							children=[html.Div(id='step_3_placeholder')],
							type='default',
							color='#94A3BC',
						),
						],
						width={'size': 1, 'offset':5}
					),
					dbc.Col(
						children=[
							html.Div(id='step_3_start_state', children=html.Img(id='step_3_end_state', src=b64_image(os.path.join('GUI', 'neutral'+'.png')), style={'width':'4vh', 'height':'4vh'})),
						],
						width={'size': 1},
					),
				]
			),
		]
	),

	html.Hr(),
	# Step 4
	dbc.Row(
		[
			html.Div(
				id='step_4_collapse_btn',
				children=[
				html.Div(html.Img(src=b64_image(os.path.join('GUI', 'steel_blue'+'.png')), style={'width':'100vw', 'height':'8vh'}), style={'position':'absolute', 'left':0}),
				html.Div([
					html.H2('Step 4: Map creation', style={'text-align':'center', 'width':'100vw', 'margin-top':'1vh', 'color':'black'}),
					],
					style={'position':'absolute', 'left':0}
				),
				],
				style={'position':'relative', 'height':'8vh'}
			),
		]
	),
	dbc.Collapse(
		id='step_4_collapse',
		children=[
			html.Hr(),
			dbc.Row(
				[
					dbc.Col(
						[
							dbc.Row(
								dbc.Col(
									html.Div(html.Img(src=b64_image(os.path.join('GUI', 'steel_blue'+'.png')), style={'width':'100%', 'height':'8vh'})),
									width={'size': 6},
								),
							)
						],
						width={'size': 1},
					),
					dbc.Col([
							html.H4('This last step will create a map and the according legend.'),
							html.H4('You have to specify a shape for each feature you want to appear on the map in the Feature shapes section.'),
							html.H4('You can let ProFeatMap choose for you by using the Automatic feature select button.'),
							html.H4('You can then change the general figure and feature parameters to explore your protein list.'),
						],
						width={'size': 10},
					),
				]
			),
			html.Hr(),
			dbc.Row(
				[
					dbc.Col(
						[
							dbc.Row(
								dbc.Col(
									html.Div(html.Img(src=b64_image(os.path.join('GUI', 'steel_blue'+'.png')), style={'width':'100%', 'height':'8vh'})),
									width={'size': 6},
								),
							)
						],
						width={'size': 1},
					),
					dbc.Col(
						[
							html.H3('Shapes and colors', style={'text-align':'center'}),
						],
						width={'size': 6},
					),
					dbc.Col(
						[
							html.H3('Cut regions (optional)', style={'text-align':'center'}),
						],
						width={'size': 3, 'offset':1},
					),
				]
			),

			dbc.Row(
				[
					dbc.Col([
						dbc.Row(
							[
								dbc.Col(
									[
										dbc.Row(
											dbc.Col(
												html.Div(html.Img(src=b64_image(os.path.join('GUI', 'steel_blue'+'.png')), style={'width':'100%', 'height':'10vh'})),
												width={'size': 6},
											),
										)
									],
									width={'size': 1},
								),
								dbc.Col([
									html.Div(
										id='feature_shape_dnd_div',
										children=dcc.Upload(
											id='feature_shape_dnd',
											children=html.Div([
												'Drag and Drop or ',
												html.A('Select Files')
											]),
											style={'height': '60px','lineHeight': '60px','borderWidth': '1px','borderStyle': 'dashed','borderRadius': '5px','textAlign': 'center','margin': '10px'},
										),
									),
								],
								width={'size': 6},
								),
								dbc.Col([
									html.Div(
										id='protein_cut_dnd_div',
										children=dcc.Upload(
											id='protein_cut_dnd',
											children=html.Div([
												'Drag and Drop or ',
												html.A('Select Files')
											]),
											style={'height': '60px','lineHeight': '60px','borderWidth': '1px','borderStyle': 'dashed','borderRadius': '5px','textAlign': 'center','margin': '10px'},
										),
									),
								],
								width={'size': 3, 'offset' : 1}
								),
							]
						),
					],

					),
				]
			),
			dbc.Row(
				[
					dbc.Col([
						dbc.Row(
							[
								dbc.Col(
									[
										dbc.Row(
											dbc.Col(
												html.Div(html.Img(src=b64_image(os.path.join('GUI', 'steel_blue'+'.png')), style={'width':'100%', 'height':'16vh'})),
												width={'size': 6},
											),
										)
									],
									width={'size': 1},
								),
								dbc.Col([
									dbc.Row(
										[
											dbc.Col(
												id='step_4_auto_feat_select_div',
												children=[
												dbc.Row(
													[
														html.H4(id='auto_feature_select_text', children='AUTOMATIC feature selection', style={'text-align':'center', 'margin-top':'0vh', 'background-color':'rgb(223, 170, 24)', 'color':'black', 'font-family':'courier', 'font-weight':'bold'}),
														dbc.Col([
															html.Div(
																children=[
																	html.Div(html.Img(id='auto_feature_select_img', src=b64_image(os.path.join('GUI', 'auto'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
																	html.Div(dbc.Button(id='auto_feature_select', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
																	dbc.Tooltip(
																			'Selects for you either the x most represented features or features with more than x occurrences in your proteins and assign a shape and color to each.',
																			target='step_4_auto_feat_select_div',
																			placement='bottom',
																			style={'font-size':tooltip_font_size}
																		)
																],
																style={'position':'relative'})
															],
														),
													]
												)
												],
												width={'size': 7, 'offset':1}
											),
											dbc.Col([
												html.H5('Lock seed', style={'text-align':'right'}),
											],
											width={'size': 3},
											),
											dbc.Col([
												dbc.Checklist(
													id='lock_seed',
													options=[
														{'label': '', 'value': 'lock_seed'},
													],
													value=[],
												),
											],
											width={'size': 1},
											),
										]
									),
									dbc.Row(
										style={'height':'1vh'}
									),
									dbc.Row(
										[	
											dbc.Col([
													dbc.Input(
														id='auto_select_threshold_1',
														placeholder='x',
														style={'width':'100%'},
													),
												],
												width={'size': 1, 'offset':2}
											),								
											dbc.Col([
													html.H5('most represented features', style={'text-align':'center', 'margin-top':'1vh'}),
												],
												width={'size': 6}
											),
											dbc.Col([
													dbc.Checklist(
														id='auto_toggle_1',
														options=[
															{'label': '', 'value': 'most_represented'},
														],
														value=[],
														style={'text-align':'center', 'margin-top':'1vh'}
													),
												],
												width={'size': 1}
											),
										]
									),
									dbc.Row(
										[	
											dbc.Col([
													html.H4('Or', style={'text-align':'right', 'margin-top':'1vh'}),
												],
												width={'size': 1, 'offset':1}
											),
											dbc.Col([
													html.H5('More than', style={'text-align':'center', 'margin-top':'1vh'}),
												],
												width={'size': 2}
											),
											dbc.Col([
													dbc.Input(
														id='auto_select_threshold_2',
														placeholder='x',
														style={'width':'100%'}
													),
												],
												width={'size': 1}
											),
											dbc.Col([
													html.H5('feature occurences', style={'text-align':'center', 'margin-top':'1vh'}),
												],
												width={'size': 4}
											),
											dbc.Col([
													dbc.Checklist(
														id='auto_toggle_2',
														options=[
															{'label': '', 'value': 'most_represented'},
														],
														value=['most_represented'],
														style={'text-align':'center', 'margin-top':'1vh'}
													),
												],
												width={'size': 1}
											),
										],
									),
									dbc.Row(
										style={'height':'1vh'}
									),
									],
									width={'size': 6},
								),
							],
						),
					]),
				]
			),

			# html.Hr(),

			dbc.Row(
				[
					dbc.Col(
						[
							dbc.Row(
								dbc.Col(
									html.Div(html.Img(id='step_4_height', src=b64_image(os.path.join('GUI', 'steel_blue'+'.png')), style={'width':'100%', 'height':'80vh'})),
									width={'size': 6},
								),
							)
						],
						width={'size': 1},
					),
					dbc.Col([
						dbc.Row(
							[
								dbc.Col([
									html.Div(id='feature_shape_datatable_div', children=dash_table.DataTable(id='feature_shape_datatable')),
								],
								),
							]
						),
						dbc.Row(
							style={'height':'1vh'}
						),
						dbc.Row(
							[
								dbc.Col(
									id='step_4_protein_shape_add_row_div',
									children=[
										dbc.Row(
											[
												dbc.Col([
													html.Div(
														children=[
															html.Div(html.Img(src=b64_image(os.path.join('GUI', 'add'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															html.Div(dbc.Button(id='feature_shape_editing_rows_btn', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															dbc.Tooltip(
																	'Add a row to the datatable',
																	target='step_4_protein_shape_add_row_div',
																	placement='bottom',
																	style={'font-size':tooltip_font_size}
															)
														],
														style={'position':'relative'})
													],
												),
											]
										)
									],
									width={'size': 1, 'offset':1}
								),
								dbc.Col(
									id='step_4_protein_shape_clear_div',
									children=[
										dbc.Row(
											[
												dbc.Col([
													html.Div(
														children=[
															html.Div(html.Img(src=b64_image(os.path.join('GUI', 'clear'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															html.Div(dbc.Button(id='feature_shape_clear', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															dbc.Tooltip(
																	'Clear datatable content',
																	target='step_4_protein_shape_clear_div',
																	placement='bottom',
																	style={'font-size':tooltip_font_size}
															)
														],
														style={'position':'relative'})
													],
												),
											]
										)
									],
									width={'size': 1, 'offset':2}
								),
								dbc.Col(
									id='step_4_protein_shape_save_div',
									children=[
										dbc.Row(
											[
												dbc.Col([
													html.Div(
														children=[
															html.Div(id='step_4_protein_shape_save_img', children=[html.Img(src=b64_image(os.path.join('GUI', 'save'+'.png')), style={'width':'4vh', 'height':'4vh'})], style={'position':'absolute'}),
															html.Div(dbc.Button(id='feature_shape_changes_save', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															dbc.Tooltip(
																	'Save datatable changes',
																	target='step_4_protein_shape_save_div',
																	placement='bottom',
																	style={'font-size':tooltip_font_size}
															)
														],
														style={'position':'relative'})
													],
												),
											]
										)
									],
									width={'size': 1, 'offset':2}
								),
								dbc.Col(
									id='step_4_protein_shape_download_div',
									children=[
										dbc.Row(
											[
												dbc.Col([
													html.Div(
														children=[
															html.Div(html.Img(src=b64_image(os.path.join('GUI', 'load_white'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															html.Div(dbc.Button(id='feature_shape_download_button', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															dbc.Tooltip(
																	'Download datatable content as excel file',
																	target='step_4_protein_shape_download_div',
																	placement='bottom',
																	style={'font-size':tooltip_font_size}
															)
														],
														style={'position':'relative', 'height':'4vh'})
													],
												),
											]
										)
									],
									width={'size': 1, 'offset':2}
								),
								Download(id='feature_shape_download'),
							]
						),
						dbc.Row(
							style={'height':'1vh'}
						),
						dbc.Row(
							[
							dbc.Col(
								children=[
								dbc.Row(
									[
										html.H4('Shapes', style={'text-align':'center'}),
										dbc.Col([
											html.Div(
												children=[
													html.Div(html.Img(src=b64_image(os.path.join('GUI', 'show'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
													html.Div(dbc.Button(id='step_4_shapes_show_button', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
												],
												style={'position':'relative'})
											],
										),
									]
								)
								],
								width={'size': 3, 'offset':1}
							),
							dbc.Col(
								children=[
								dbc.Row(
									[
										html.H4('Colors', style={'text-align':'center'}),
										dbc.Col([
											html.Div(
												children=[
													html.Div(html.Img(src=b64_image(os.path.join('GUI', 'show'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
													html.Div(dbc.Button(id='step_4_colors_show_button', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
												],
												style={'position':'relative'})
											],
										),
									]
								)
								],
								width={'size': 3}
							),
							dbc.Col(
								children=[
								dbc.Row(
									[
										html.H4('Colormaps', style={'text-align':'center'}),
										dbc.Col([
											html.Div(
												children=[
													html.Div(html.Img(src=b64_image(os.path.join('GUI', 'show'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
													html.Div(dbc.Button(id='step_4_colormaps_show_button', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
												],
												style={'position':'relative'})
											],
										),
									]
								)
								],
								width={'size': 3}
							),
							]
						),
						dbc.Row(
							style={'height':'1vh'}
						),
						dbc.Collapse(
							id='shape_and_color_collapse',
							children=[
								dbc.Row(
									[
										dbc.Col([
											html.Div(
												children=[
													html.Div(html.Img(src=b64_image(os.path.join('GUI', 'hide'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
													html.Div(dbc.Button(id='step_4_shapes_and_colors_hide_button', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
												],
												style={'position':'relative'})
											],
										),
										dbc.Col(
											[
												html.Div(html.Img(id='shapes_colors_colormaps_img', style={'width':'100%', 'height':'100%'})),
											],
											width={'size': 12},
										),
									]
								),
							],
						),
						],
						width={'size': 6}

					),


					dbc.Col([
						dbc.Row(
							[
								dbc.Col([
									html.Div(id='protein_cut_datatable_div', children=dash_table.DataTable(id='protein_cut_datatable')),
								],
								),
							]
						),
						dbc.Row(
							style={'height':'1vh'}
						),
						dbc.Row(
							[
								dbc.Col(
									id='step_4_protein_cut_add_row_div',
									children=[
										dbc.Row(
											[
												dbc.Col([
													html.Div(
														children=[
															html.Div(html.Img(src=b64_image(os.path.join('GUI', 'add'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															html.Div(dbc.Button(id='protein_cut_editing_rows_btn', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															dbc.Tooltip(
																	'Add a row to the datatable',
																	target='step_4_protein_cut_add_row_div',
																	placement='bottom',
																	style={'font-size':tooltip_font_size}
															)
														],
														style={'position':'relative'})
													],
												),
											]
										)
									],
									width={'size': 1, 'offset':1}
								),
								dbc.Col(
									id='step_4_protein_cut_clear_div',
									children=[
										dbc.Row(
											[
												dbc.Col([
													html.Div(
														children=[
															html.Div(html.Img(src=b64_image(os.path.join('GUI', 'clear'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															html.Div(dbc.Button(id='protein_cut_clear', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															dbc.Tooltip(
																	'Clear datatable content',
																	target='step_4_protein_cut_clear_div',
																	placement='bottom',
																	style={'font-size':tooltip_font_size}
															)
														],
														style={'position':'relative'})
													],
												),
											]
										)
									],
									width={'size': 1, 'offset':2}
								),
								dbc.Col(
									id='step_4_protein_cut_save_div',
									children=[
										dbc.Row(
											[
												dbc.Col([
													html.Div(
														children=[
															html.Div(id='step_4_save_img', children=[html.Img(src=b64_image(os.path.join('GUI', 'save'+'.png')), style={'width':'4vh', 'height':'4vh'})], style={'position':'absolute'}),
															html.Div(dbc.Button(id='protein_cut_changes_save', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															dbc.Tooltip(
																	'Save datatable changes',
																	target='step_4_protein_cut_save_div',
																	placement='bottom',
																	style={'font-size':tooltip_font_size}
															)
														],
														style={'position':'relative'})
													],
												),
											]
										)
									],
									width={'size': 1, 'offset':2}
								),
								dbc.Col(
									id='step_4_protein_cut_download_div',
									children=[
										dbc.Row(
											[
												dbc.Col([
													html.Div(
														children=[
															html.Div(html.Img(src=b64_image(os.path.join('GUI', 'load_white'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															html.Div(dbc.Button(id='protein_cut_download_button', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															dbc.Tooltip(
																	'Download datatable content as excel file',
																	target='step_4_protein_cut_download_div',
																	placement='bottom',
																	style={'font-size':tooltip_font_size}
															)
														],
														style={'position':'relative', 'height':'4vh'})
													],
												),
											]
										)
									],
									width={'size': 1, 'offset':2}
								),
								Download(id='protein_cut_download'),
							],
						),
						dbc.Row(
							style={'height':'10vh'}
						),

						dbc.Row(
							id='step_4_options1_div',
							children=[
								html.H4('Shapes and colors example', style={'margin-top':'1vh'}),
								dbc.Col(
									children=[
										html.Div(
											children=[
												html.Div(html.Img(src=b64_image(os.path.join('GUI', 'load_white'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
												html.Div(dbc.Button(id='shapes_and_colors_example_download_button', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
												dbc.Tooltip(
														'Download an example of shapes and colors file for TRIO_HUMAN',
														target='shapes_and_colors_example_download_button',
														placement='bottom',
														style={'font-size':tooltip_font_size}
												),
												Download(id='shapes_and_colors_example_download'),
											],
											# style={'position':'relative'}
											)
										],
										width={'size': 2},
								),
								dbc.Col(
									children=[
										html.Div(
											children=[
												html.Div(html.Img(src=b64_image(os.path.join('GUI', 'load_white'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
												html.Div(dbc.Button(id='shapes_and_colors_example_download_button_2', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
												dbc.Tooltip(
														'Download an example of shapes and colors file for PDZ domain containing proteins',
														target='shapes_and_colors_example_download_button_2',
														placement='bottom',
														style={'font-size':tooltip_font_size}
												),
												Download(id='shapes_and_colors_example_download_2'),
											],
											# style={'position':'relative'}
											)
										],
										width={'size': 2},
								),
							],
						),
						dbc.Row(
							id='step_4_options2_div',
							children=[
								html.H4('Cut regions example', style={'margin-top':'1vh'}),
								dbc.Col(
									children=[
										html.Div(
											children=[
												html.Div(html.Img(src=b64_image(os.path.join('GUI', 'load_white'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
												html.Div(dbc.Button(id='cut_regions_example_download_button', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
												dbc.Tooltip(
														'Download an example of cut regions file for TRIO_HUMAN',
														target='cut_regions_example_download_button',
														placement='bottom',
														style={'font-size':tooltip_font_size}
												),
												Download(id='cut_regions_example_download'),
											],
											# style={'position':'relative'}
											)
										],
										width={'size': 2},
								),
								dbc.Col(
									children=[
										html.Div(
											children=[
												html.Div(html.Img(src=b64_image(os.path.join('GUI', 'load_white'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
												html.Div(dbc.Button(id='cut_regions_example_download_button_2', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
												dbc.Tooltip(
														'Download an example of cut regions file for PDZ domain containing proteins',
														target='cut_regions_example_download_button_2',
														placement='bottom',
														style={'font-size':tooltip_font_size}
												),
												Download(id='cut_regions_example_download_2'),
											],
											# style={'position':'relative'}
											)
										],
										width={'size': 2},
								),
							],
						),

						],
						width={'size': 3, 'offset' : 1}
					)
				]
			),
			dbc.Row(
				[
					dbc.Col(
						[
							dbc.Row([
								dbc.Col(
									html.Div(html.Img(src=b64_image(os.path.join('GUI', 'steel_blue'+'.png')), style={'width':'100%', 'height':'6vh'})),
									width={'size': 6},
								),
								]
							)
						],
						width={'size': 1},
					),
				]
			),
			# html.Hr(),
			dbc.Row(
				[
					dbc.Col(
						[
							dbc.Row(
								dbc.Col(
									html.Div(html.Img(src=b64_image(os.path.join('GUI', 'steel_blue'+'.png')), style={'width':'100%', 'height':'8vh'})),
									width={'size': 6},
								),
							)
						],
						width={'size': 1},
					),
					dbc.Col(
						children=[html.H3('Map parameters', style={'text-align':'center',})],
						width={'size': 10},
					),
				]
			),
			dbc.Row(
				[
					dbc.Col(
						[
							dbc.Row(
								dbc.Col(
									html.Div(html.Img(src=b64_image(os.path.join('GUI', 'steel_blue'+'.png')), style={'width':'100%', 'height':'50vh'})),
									width={'size': 6},
								),
							)
						],
						width={'size': 1},
					),
					dbc.Col([
						dbc.Row(
							[
								dbc.Col([
									html.H4('Case to draw', style={'text-align':'right',}),
								],
								width={'size': 3},
								),

								dbc.Col([
									dcc.Dropdown(
										id='case_dropdown',
										style={'color':'black', 'width':'100%'},
										placeholder='Case',
										clearable=True
									),
								],
								width={'size': 3},
								),
								dbc.Col([
									html.H4('Focus on', style={'text-align':'right',}),
								],
								width={'size': 3},
								),

								dbc.Col([
									dcc.Input(
										#Have to be flaot (will need casting of value given in input as float). Revove arrows ?
										id='focus_input',
										type='text',
										placeholder='enter feature name here',
										style={'width':'100%'}
									),
								],
								width={'size': 3},
								),
							]
						),
						dbc.Row(
							style={'height':'1vh'}
						),
						dbc.Row(
							[
								dbc.Col([
									html.H4('Sorting', style={'text-align':'right'}),
								],
								width={'size': 3},
								),
								dbc.Col([
									dcc.Dropdown(
										id='sorting_dropdown',
										options=[{"label": x, "value": x} for x in ['None', 'abc', 'value', 'feature_number_distance']],
										value='feature_number_distance',
										style={'color':'black', 'width':'100%'},
										clearable=False
									),
								],
								width={'size': 3},
								),
								dbc.Col([
									html.H4('Threshold', style={'text-align':'right'}),
								],
								width={'size': 3},
								),
								dbc.Col([
									dcc.Input(
										#Have to be flaot (will need casting of value given in input as float). Revove arrows ?
										id='threshold',
										type='number',
										placeholder='default: None',
										style={'width':'100%'}
									),
								],
								width={'size': 3},
								),
							],
						),
						dbc.Row(
							style={'height':'2vh'}
						),
						dbc.Row(
							[
								dbc.Col([
									html.H4('General figure parameters', style={'text-align':'center'}),
									dbc.Row(
										[
											dbc.Col([
												html.H5('Horizontal stretch factor', style={'text-align':'right'}),
											],
											width={'size': 6},
											),
											dbc.Col([
												dcc.Input(
													#Have to be float (will need casting of value given in input as float). Revove arrows ?
													id='length_factor',
													type='number',
													placeholder='default: 1',
													style={'width':'100%'}
												),
											],
											width={'size': 6},
											),
										]
									),
									dbc.Row(
										[
											dbc.Col([
												html.H5('Vertical height', style={'text-align':'right'}),
											],
											width={'size': 6},
											),
											dbc.Col([
												dcc.Input(
													id='height',
													type='number',
													placeholder='default: 20',
													style={'width':'100%'}
												),
											],
											width={'size': 6},
											),
										]
									),
									dbc.Row(
										[
											dbc.Col([
												html.H5('Protein thickness', style={'text-align':'right'}),
											],
											width={'size': 6},
											),
											dbc.Col([
												dcc.Input(
													id='pensize',
													type='number',
													placeholder='default: 3',
													style={'width':'100%'}
												),
											],
											width={'size': 6},
											),
										]
									),
									dbc.Row(
										[
											dbc.Col([
												html.H5('Protein name size', style={'text-align':'right'}),
											],
											width={'size': 6},
											),
											dbc.Col([
												dcc.Input(
													id='text_size',
													type='number',
													placeholder='default: 30',
													style={'width':'100%'}
												),
											],
											width={'size': 6},
											),
										]
									),
									dbc.Row(
										[
											dbc.Col([
												html.H5('Biased regions text size', style={'text-align':'right'}),
											],
											width={'size': 6},
											),
											dbc.Col([
												dcc.Input(
													id='biased_region_text_size',
													type='number',
													placeholder='default: 20',
													style={'width':'100%'}
												),
											],
											width={'size': 6},
											),
										]
									),
									dbc.Row(
										[
											dbc.Col([
												html.H5('Show protein length', style={'text-align':'right'}),
											],
											width={'size': 6},
											),
											dbc.Col([
												dbc.Checklist(
													id='protein_length_toggles',
													options=[
														{'label': '', 'value': 'show_length'},
													],
													value=['show_length'],
												),
											],
											width={'size': 6},
											),
										]
									),
									dbc.Row(
										[
											dbc.Col([
												html.H5('Consistent shape fill/contour', style={'text-align':'right'}),
											],
											width={'size': 6},
											),
											dbc.Col([
												dbc.Checklist(
													id='uniform_toggles',
													options=[
														{'label': '', 'value': 'uniform'},
													],
												),
											],
											width={'size': 6},
											),
										]
									),
								],
								width={'size': 6},
								),
								dbc.Col([
									html.H4('Feature parameters', style={'text-align':'center'}),
									dbc.Row(
										[
											dbc.Col([
												html.H5('3D structure coverage', style={'text-align':'right'}),
											],
											width={'size': 6},
											),
											dbc.Col([
												dbc.Checklist(
													id='coverage_toggles',
													options=[
														{'label': '', 'value': 'coverage'},
													],
													value=[],
												),
											],
											width={'size': 6},
											),
										]
									),
									dbc.Row(
										[
											dbc.Col([
												html.H5('Secondary structure', style={'text-align':'right'}),
											],
											width={'size': 6},
											),
											dbc.Col([
												dbc.Checklist(
													id='2nd_struct_toggles',
													options=[
														{'label': '', 'value': '2nd_struct'},
													],
													value=[],
												),
											],
											width={'size': 6},
											),
										]
									),
									dbc.Row(
										[
											dbc.Col([
												html.H5('Disorder', style={'text-align':'right'}),
											],
											width={'size': 6},
											),
											dbc.Col([
												dbc.Checklist(
													id='disorder_toggles',
													options=[
														{'label': '', 'value': 'disorder'},
													],
													value=['disorder'],
												),
											],
											width={'size': 6},
											),
										]
									),
									dbc.Row(
										[
											dbc.Col([
												html.H5('Modified residues', style={'text-align':'right'}),
											],
											width={'size': 6},
											),
											dbc.Col([
												dbc.Checklist(
													id='mod_res_toggles',
													options=[
														{'label': '', 'value': 'mod_res'},
													],
													value=[],
												),
											],
											width={'size': 6},
											),
										]
									),
									dbc.Row(
										[
											dbc.Col([
												html.H5('Composition biased regions', style={'text-align':'right'}),
											],
											width={'size': 6},
											),
											dbc.Col([
												dcc.Dropdown(
													id='draw_region_dropdown',
													options=[
													{"label": '(A) Acidic residues', "value": 'A'},
													{"label": '(B) Basic residues', "value": 'B'},
													{"label": '(C) Basic and acidic residues', "value": 'C'},
													{"label": '(H) Polar residues', "value": 'H'},
													{"label": '(P) Pro residues', "value": 'P'}
													],
													value='None',
													style={'color':'black', 'width':'100%'},
													clearable=True,
													multi=True
												),
											],
											width={'size': 6},
											),
										],
									),
									dbc.Row(
										style={'height':'1vh'}
									),									
									dbc.Row(
										[									
											dbc.Col(
												children=[
												dbc.Row(
													[
														dbc.Col([
															html.H5('Feature parameters defaults', style={'text-align':'right'}),
															],
															width={'size': 6},
														),
														dbc.Col([
															html.Div(
																children=[
																	html.Div(html.Img(src=b64_image(os.path.join('GUI', 'option_run'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
																	html.Div(dbc.Button(id='default_feat_param', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
																],
																style={'position':'relative'})
															],
															width={'size': 6},
														),
													]
												)
												],
												width={'size': 12}
											),
										]
									),
								],
								width={'size': 6},
								),
							],
						),
						dbc.Row(
							style={'height':'2vh'}
						),
						dbc.Row(
							dbc.Col([
								html.H4('Order of feature drawing', style={'text-align':'center'}),
							],
							width={'size': 12},
							),
						),
						dcc.Dropdown(
							id='FT_order_list',
							options=[{"label": x, "value": x} for x in 'DISORDER,DOMAIN,CHAIN,INIT_MET,PEPTIDE,ZN_FING,DNA_BIND,REGION,ACT_SITE,METAL,SITE,LIPID,HELIX,STRAND,TURN,CONFLICT,CARBOHYD,BINDING, MOTIF,MOD_RES,COMPBIAS,REPEAT,VARIANT,PDB'.split(',')],
							value=[x for x in 'DISORDER,DOMAIN,CHAIN,INIT_MET,PEPTIDE,ZN_FING,DNA_BIND,REGION,ACT_SITE,METAL,SITE,LIPID,HELIX,STRAND,TURN,CONFLICT,CARBOHYD,BINDING, MOTIF,MOD_RES,COMPBIAS,REPEAT,VARIANT,PDB'.split(',')],
							style={'color':'black', 'width':'100%'},
							multi=True,
						),
						dbc.Row(
							style={'height':'2vh'}
						),
						dbc.Row(
							[
								html.H4('Map parameters example', style={'margin-top':'1vh'}),
								dbc.Col(
									children=[
										html.Div(
											children=[
												html.Div(html.Img(src=b64_image(os.path.join('GUI', 'load_white'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
												html.Div(dbc.Button(id='map_parameters_example_download_button', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
												dbc.Tooltip(
														'Download an example of map parameters file for TRIO_HUMAN',
														target='map_parameters_example_download_button',
														placement='bottom',
														style={'font-size':tooltip_font_size}
												),
												Download(id='map_parameters_example_download'),
											],
											# style={'position':'relative'}
											)
										],
										width={'size': 1},
								),
								dbc.Col(
									children=[
										html.Div(
											children=[
												html.Div(html.Img(src=b64_image(os.path.join('GUI', 'load_white'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
												html.Div(dbc.Button(id='map_parameters_example_download_button_2', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
												dbc.Tooltip(
														'Download an example of map parameters file for PDZ domain containing proteins',
														target='map_parameters_example_download_button_2',
														placement='bottom',
														style={'font-size':tooltip_font_size}
												),
												Download(id='map_parameters_example_download_2'),
											],
											# style={'position':'relative'}
											)
										],
										width={'size': 1},
								),
								dbc.Col(
									id='step_4_map_parameters_clear_div',
									children=[
										dbc.Row(
											[
												dbc.Col([
													html.Div(
														children=[
															html.Div(html.Img(src=b64_image(os.path.join('GUI', 'default'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															html.Div(dbc.Button(id='parameters_clear', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															dbc.Tooltip(
																	'Reset to default parameters',
																	target='step_4_protein_cut_clear_div',
																	placement='bottom',
																	style={'font-size':tooltip_font_size}
															)
														],
														style={'position':'relative'})
													],
												),
											]
										)
									],
									width={'size': 1, 'offset':1}
								),
								dbc.Col(
									id='step_4_map_parameters_save_div',
									children=[
										dbc.Row(
											[
												dbc.Col([
													html.Div(
														children=[
															html.Div(children=[html.Img(src=b64_image(os.path.join('GUI', 'save'+'.png')), style={'width':'4vh', 'height':'4vh'})], style={'position':'absolute'}),
															html.Div(dbc.Button(id='parameters_changes_save', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															dbc.Tooltip(
																	'Save parameters changes. Saving is not needed to impact the map creation.',
																	target='step_4_map_parameters_save_div',
																	placement='bottom',
																	style={'font-size':tooltip_font_size}
															)
														],
														style={'position':'relative'})
													],
												),
											]
										)
									],
									width={'size': 1, 'offset':1}
								),
								dbc.Col(
									id='step_4_map_parameters_download_div',
									children=[
										dbc.Row(
											[
												dbc.Col([
													html.Div(
														children=[
															html.Div(html.Img(src=b64_image(os.path.join('GUI', 'load_white'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															html.Div(dbc.Button(id='parameters_download_button', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
															dbc.Tooltip(
																	'Download datatable content as excel file',
																	target='step_4_map_parameters_download_div',
																	placement='bottom',
																	style={'font-size':tooltip_font_size}
															)
														],
														style={'position':'relative', 'height':'4vh'})
													],
												),
											]
										)
									],
									width={'size': 1, 'offset':1}
								),
								Download(id='parameters_download'),
							]
						),
					],
					width={'size': 8, 'offset':1},
					),
				]
			),
			# html.Hr(),
			dbc.Row(
				style={'height':'2vh'}
			),	
			dbc.Row(
				[
					dbc.Col(
						[
							dbc.Row(
								dbc.Col(
									html.Div(html.Img(src=b64_image(os.path.join('GUI', 'steel_blue'+'.png')), style={'width':'100%', 'height':'6vh'})),
									width={'size': 6},
								),
							)
						],
						width={'size': 1},
					),
					dbc.Col(
						id='step_4_run_div',
						children=[
						dbc.Row(
							[
								html.H3('Run', style={'text-align':'center', 'margin-top':'0vh'}),
								dbc.Col([
									html.Div(
										children=[
											html.Div(html.Img(src=b64_image(os.path.join('GUI', 'run'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
											html.Div(dbc.Button(id='step_4_launch_button', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
											dbc.Tooltip(
													'Start map creation',
													target='step_4_run_div',
													placement='bottom',
													style={'font-size':tooltip_font_size}
												)
										],
										style={'position':'relative'})
									],
								),
							]
						),
						dcc.Loading(
							id='step_4_loading',
							children=[html.Div(id='step_4_placeholder')],
							type='default',
							color='#94A3BC',
						),
						],
						width={'size': 1, 'offset':5}
					),


					dbc.Col(
						children=[
							html.Div(id='step_4_start_state', children=html.Img(id='step_4_end_state', src=b64_image(os.path.join('GUI', 'neutral'+'.png')), style={'width':'4vh', 'height':'4vh'})),
						],
						width={'size': 1},
					),
				]
			),
			html.Hr(),
		]
	),
	dbc.Row(
		style={'height':'2vh'}
	),
	dbc.Row(
		[
			dbc.Col(
				html.Div(html.H4('Latest sorted protein list', style={'text-align':'right', 'width':'100%'}), style={'margin-top':'1vh'}),
				width={'size': 2, 'offset': 1},
			),
			dbc.Col(
				id='step_4_sorting_download_div',
				children=[
					dbc.Row(
						[
							dbc.Col([
								html.Div(
									children=[
										html.Div(html.Img(src=b64_image(os.path.join('GUI', 'load_white'+'.png')), style={'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
										html.Div(dbc.Button(id='protein_sorting_download_button', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
										dbc.Tooltip(
												'Download latest run sorted protein list',
												target='step_4_sorting_download_div',
												placement='bottom',
												style={'font-size':tooltip_font_size}
										),
										Download(id='protein_sorting_download'),
									],
									style={'position':'relative', 'height':'4vh'})
								],
							),
						]
					)
				],
				width={'size': 1}
			),
			dbc.Col(
				html.Div(html.H4('Map and Legend', style={'text-align':'right', 'width':'100%'}), style={'margin-top':'1vh'}),
				width={'size': 2, 'offset': 1},
			),
			dbc.Col(
				id='step_4_figure_download_div',
				children=[
					dbc.Row(
						[
							dbc.Col([
								html.Div(
									children=[
										html.Div(children=[html.Img(src=b64_image(os.path.join('GUI', 'load_white'+'.png')), style={'width':'4vh', 'height':'4vh'})], style={'position':'absolute'}),
										html.Div(dbc.Button(id='normal_run_result_download', style={'background-color':'transparent', 'width':'4vh', 'height':'4vh'}), style={'position':'absolute'}),
										dbc.Tooltip(
											'Download Map and Legend',
											target='step_4_figure_download_div',
											placement='bottom',
											style={'font-size':tooltip_font_size}
										),
										Download(id='NR_figure_download'),
										Download(id='NR_legend_download'),
									],
									style={'position':'relative', 'height':'4vh'})
								],
							),
						]
					)
				],
				width={'size': 1}
			),
		],
	),
	dbc.Row(
		[
			dbc.Col([
					html.Div(
						id='result_img',
					)
				],
				width={'size': 10, 'offset':1}
			),
			Download(id='map_download'),
		]
	),
	html.Hr(),
	dbc.Row(
		[
			dbc.Col([
					html.Div(
						id='legend_img',
					)
				],
				width={'size': 10, 'offset':1}
			),
			Download(id='legend_download'),
		]
	),


	html.Div(id='placeholder'),
	html.Div(id='placeholder_parameters'),

])

@app.callback(
	Output('quick_guide_collapse', 'is_open'),
	Input('quick_run_btn', 'n_clicks'),
	Input('quick_run_btn_2', 'n_clicks'),
	State('quick_guide_collapse', 'is_open'),
	)
def toggle_collapse(n, n2, is_open):
	if n:
		return not is_open
	return is_open

@app.callback(
	Output('step_1_collapse', 'is_open'),
	Input('step_1_collapse_btn', 'n_clicks'),
	State('step_1_collapse', 'is_open'),
	)
def toggle_collapse(n, is_open):
	if n:
		return not is_open
	return is_open

@app.callback(
	Output('step_2_collapse', 'is_open'),
	Input('step_2_collapse_btn', 'n_clicks'),
	State('step_2_collapse', 'is_open'),
	)
def toggle_collapse(n, is_open):
	if n:
		return not is_open
	return is_open

@app.callback(
	Output('step_3_collapse', 'is_open'),
	Input('step_3_collapse_btn', 'n_clicks'),
	State('step_3_collapse', 'is_open'),
	)
def toggle_collapse(n, is_open):
	if n:
		return not is_open
	return is_open

@app.callback(
	Output('step_4_collapse', 'is_open'),
	Input('step_4_collapse_btn', 'n_clicks'),
	State('step_4_collapse', 'is_open'),
	)
def toggle_collapse(n, is_open):
	if n:
		return not is_open
	return is_open

@app.callback(
	Output('shape_and_color_collapse', 'is_open'),
	Output('shapes_colors_colormaps_img', 'src'),
	Output('last_help', 'data'),

	Input('step_4_shapes_show_button', 'n_clicks'),
	Input('step_4_colors_show_button', 'n_clicks'),
	Input('step_4_colormaps_show_button', 'n_clicks'),
	Input('step_4_shapes_and_colors_hide_button', 'n_clicks'),

	State('shape_and_color_collapse', 'is_open'),
	State('last_help', 'data'),

	prevent_initial_call=True
	)
def toggle_collapse(shapes_show, colors_show, colormaps_show, shapes_and_colors_hide, is_open, last_help):
	ctx = dash.callback_context
	triggered = ctx.triggered[0]['prop_id'].split('.')[0]
	print("TRIGERED_load_param: "+triggered)

	if triggered == 'step_4_shapes_show_button' :
		if is_open and last_help == 'step_4_shapes_show_button':
			return False, no_update, no_update
		else :
			src = b64_image(os.path.join('GUI', 'shapes_examples'+'.png'))
			return True, src, 'step_4_shapes_show_button'

	elif triggered == 'step_4_colors_show_button' :
		if is_open and last_help == 'step_4_colors_show_button':
			return False, no_update, no_update
		else :
			src = b64_image(os.path.join('GUI', 'colors_examples'+'.png'))
			return True, src, 'step_4_colors_show_button'

	elif triggered == 'step_4_colormaps_show_button' :
		if is_open and last_help == 'step_4_colormaps_show_button':
			return False, no_update, no_update
		else :
			src = b64_image(os.path.join('GUI', 'colormaps_examples'+'.png'))
			return True, src, 'step_4_colormaps_show_button'

@app.callback(

	Output('sorting_dropdown', 'value'), # Sorting
	Output('focus_input', 'value'), # Focus on
	Output('threshold', 'value'), # Threshold

	#General figure parameters
	Output('length_factor', 'value'), # Horizontal stretch factor
	Output('height', 'value'), # Vertical stretch factor
	Output('pensize', 'value'), # Protein thickness
	Output('text_size', 'value'), # Protein name size
	Output('biased_region_text_size', 'value'), # Biased regions text size
	Output('protein_length_toggles', 'value'), # Show protein length
	Output('uniform_toggles', 'value'), # Uniform shape fill/contour

	# Feature parameters
	Output('coverage_toggles', 'value'), # 3D structure coverage
	Output('2nd_struct_toggles', 'value'), # Secondary structure
	Output('disorder_toggles', 'value'), # Disorder
	Output('mod_res_toggles', 'value'), # Modified residues
	Output('draw_region_dropdown', 'value'), # Composition biased regions

	Output('FT_order_list', 'value'), # Order of feature drawing

	Input('placeholder_parameters', 'n_clicks'),
	Input('parameters_clear', 'n_clicks'),

	State('Step_4_parameters', 'data'),
	)
def download_uniprot_files(n, clear, parameters):
	ctx = dash.callback_context
	triggered = ctx.triggered[0]['prop_id'].split('.')[0]
	print("TRIGERED_load_param: "+triggered)

	if triggered == '' and parameters != None :
		sorting = parameters['Sorting']
		focus = parameters['Focus on']
		threshold = parameters['Threshold']
		if threshold != None :
			threshold = float(str(threshold).replace(',', '.'))

		length_factor = parameters['Horizontal stretch factor']
		height_factor = parameters['Vertical stretch factor']
		pensize = parameters['Protein thickness']
		text_size = parameters['Protein name size']
		biased_region_text_size = parameters['Biased regions text size']
		protein_length = parameters['Show protein length']
		if protein_length == None :
			protein_length = []
		uniform = parameters['Uniform shape fill/contour']
		if uniform == None :
			uniform = []

		coverage = parameters['3D structure coverage']
		if coverage == None :
			coverage = []
		secondary_struct = parameters['Secondary structure']
		if secondary_struct == None :
			secondary_struct = []
		disorder = parameters['Disorder']
		if disorder == None :
			disorder = []
		mod_res = parameters['Modified residues']
		if mod_res == None :
			mod_res = []

		draw_region = parameters['Composition biased regions']

		FT_order_list = parameters['Order of feature drawing']
		if FT_order_list != None :
			FT_order_list.split(', ')

		# print(parameters)

		return sorting, focus, threshold, length_factor, height_factor, pensize, text_size, biased_region_text_size, protein_length, uniform, coverage, secondary_struct, disorder, mod_res, draw_region, FT_order_list.split(', ')


	elif triggered == 'parameters_clear' :

		sorting = 'feature_number_distance'
		focus = None
		threshold = None

		length_factor = None
		height_factor = None
		pensize = None
		text_size = None
		biased_region_text_size = None
		protein_length = ['show_length']
		uniform = []

		coverage = []
		secondary_struct = []
		disorder = ['disorder']
		mod_res = []

		draw_region = []

		FT_order_list = [x for x in 'DISORDER,DOMAIN,CHAIN,INIT_MET,PEPTIDE,ZN_FING,DNA_BIND,REGION,ACT_SITE,METAL,SITE,LIPID,HELIX,STRAND,TURN,CONFLICT,CARBOHYD,BINDING, MOTIF,MOD_RES,COMPBIAS,REPEAT,VARIANT,PDB'.split(',')]

		return sorting, focus, threshold, length_factor, height_factor, pensize, text_size, biased_region_text_size, protein_length, uniform, coverage, secondary_struct, disorder, mod_res, draw_region, FT_order_list

	else :
		return [no_update for i in range(16)]

@app.callback(
	Output('protein_sorting_download', 'data'),

	Input('protein_sorting_download_button', 'n_clicks'),

	State('Sorted_protein_list', 'data'),
	State('sorting_dropdown', 'value'),
	State('Step_1', 'data'),

	prevent_initial_call=True
	)
def download_uniprot_files(button_click, sorted_list, sorting, data):
	if sorted_list == None :
		sorted_list = []

	protein_list_df = pd.read_json(data, orient='split')

	df = protein_list_df.set_index('protein', drop=False)
	df = df.loc[sorted_list]

	return send_data_frame(df.to_excel, filename=sorting+'.xlsx', index=False)

@app.callback(
	Output('parameters_download', 'data'),

	Input('parameters_download_button', 'n_clicks'),

	State('Step_4_parameters', 'data'),

	prevent_initial_call=True
	)
def download_uniprot_files(button_click, parameters):

	# print(parameters)

	df = pd.DataFrame(columns=['Parameter', 'Value'])
	df_i = 0

	if parameters['Case to draw'] == None :
		df.loc[df_i] = ['Case to draw', None]
		df_i += 1
	elif parameters['Case to draw'] != None :
		df.loc[df_i] = ['Case to draw', parameters['Case to draw']]
		df_i += 1
	if parameters['Sorting'] == None :
		df.loc[df_i] = ['Sorting', None]
		df_i += 1
	elif parameters['Sorting'] != None :
		df.loc[df_i] = ['Sorting', parameters['Sorting']]
		df_i += 1
	if parameters['Focus on'] == None :
		df.loc[df_i] = ['Focus on', None]
		df_i += 1
	elif parameters['Focus on'] != None :
		df.loc[df_i] = ['Focus on', str(parameters['Focus on'])]
		df_i += 1
	if parameters['Threshold'] == None :
		df.loc[df_i] = ['Threshold', None]
		df_i += 1
	elif parameters['Threshold'] != None :
		df.loc[df_i] = ['Threshold', str(parameters['Threshold'])]
		df_i += 1

	df.loc[df_i] = [None, None]
	df_i += 1

	df.loc[df_i] = ['General figure parameters', None]
	df_i += 1

	if parameters['Horizontal stretch factor'] == None :
		df.loc[df_i] = ['Horizontal stretch factor', None]
		df_i += 1
	elif parameters['Horizontal stretch factor'] != None :
		df.loc[df_i] = ['Horizontal stretch factor', str(parameters['Horizontal stretch factor'])]
		df_i += 1
	if parameters['Vertical stretch factor'] == None :
		df.loc[df_i] = ['Vertical stretch factor', None]
		df_i += 1
	elif parameters['Vertical stretch factor'] != None :
		df.loc[df_i] = ['Vertical stretch factor', str(parameters['Vertical stretch factor'])]
		df_i += 1
	if parameters['Protein thickness'] == None :
		df.loc[df_i] = ['Protein thickness', None]
		df_i += 1
	elif parameters['Protein thickness'] != None :
		df.loc[df_i] = ['Protein thickness', str(parameters['Protein thickness'])]
		df_i += 1
	if parameters['Protein name size'] == None :
		df.loc[df_i] = ['Protein name size', None]
		df_i += 1
	elif parameters['Protein name size'] != None :
		df.loc[df_i] = ['Protein name size', str(parameters['Protein name size'])]
		df_i += 1
	if parameters['Biased regions text size'] == None :
		df.loc[df_i] = ['Biased regions text size', None]
		df_i += 1
	elif parameters['Biased regions text size'] != None :
		df.loc[df_i] = ['Biased regions text size', str(parameters['Biased regions text size'])]
		df_i += 1
	if parameters['Show protein length'] == None or parameters['Show protein length'] == [] :
		df.loc[df_i] = ['Show protein length', 'No']
		df_i += 1
	else :
		df.loc[df_i] = ['Show protein length', 'Yes']
		df_i += 1
	if parameters['Uniform shape fill/contour'] == None or parameters['Uniform shape fill/contour'] == [] :
		df.loc[df_i] = ['Uniform shape fill/contour', 'No']
		df_i += 1
	else :
		df.loc[df_i] = ['Uniform shape fill/contour', 'Yes']
		df_i += 1

	df.loc[df_i] = [None, None]
	df_i += 1

	df.loc[df_i] = ['Feature parameters', None]
	df_i += 1

	if parameters['3D structure coverage'] == None or parameters['3D structure coverage'] == [] :
		df.loc[df_i] = ['3D structure coverage', 'No']
		df_i += 1
	else :
		df.loc[df_i] = ['3D structure coverage', 'Yes']
		df_i += 1
	if parameters['Secondary structure'] == None or parameters['Secondary structure'] == [] :
		df.loc[df_i] = ['Secondary structure', 'No']
		df_i += 1
	else :
		df.loc[df_i] = ['Secondary structure', 'Yes']
		df_i += 1
	if parameters['Disorder'] == None or parameters['Disorder'] == [] :
		df.loc[df_i] = ['Disorder', 'No']
		df_i += 1
	else :
		df.loc[df_i] = ['Disorder', 'Yes']
		df_i += 1
	if parameters['Modified residues'] == None or parameters['Modified residues'] == [] :
		df.loc[df_i] = ['Modified residues', 'No']
		df_i += 1
	else :
		df.loc[df_i] = ['Modified residues', 'Yes']
		df_i += 1
	if parameters['Composition biased regions'] == None or parameters['Composition biased regions'] == [] :
		df.loc[df_i] = ['Composition biased regions', None]
		df_i += 1
	else :
		df.loc[df_i] = ['Composition biased regions', ', '.join(parameters['Composition biased regions'])]
		df_i += 1

	df.loc[df_i] = [None, None]
	df_i += 1

	if parameters['Order of feature drawing'] == None or parameters['Order of feature drawing'] == [] :
		df.loc[df_i] = ['Order of feature drawing', None]
		df_i += 1
	else :
		df.loc[df_i] = ['Order of feature drawing', parameters['Order of feature drawing']]
		df_i += 1

	# print(df)

	return send_data_frame(df.to_excel, filename='parameters.xlsx', index=False, header=False)

# Save parameters
@app.callback(
	Output('Step_4_parameters', 'data'),

	Input('parameters_changes_save', 'n_clicks'),

	State('case_dropdown', 'value'), # Case to draw
	State('sorting_dropdown', 'value'), # Sorting
	State('focus_input', 'value'), # Focus on
	State('threshold', 'value'), # Threshold

	#General figure parameters
	State('length_factor', 'value'), # Horizontal stretch factor
	State('height', 'value'), # Vertical stretch factor
	State('pensize', 'value'), # Protein thickness
	State('text_size', 'value'), # Protein name size
	State('biased_region_text_size', 'value'), # Biased regions text size
	State('protein_length_toggles', 'value'), # Show protein length
	State('uniform_toggles', 'value'), # Uniform shape fill/contour

	# Feature parameters
	State('coverage_toggles', 'value'), # 3D structure coverage
	State('2nd_struct_toggles', 'value'), # Secondary structure
	State('disorder_toggles', 'value'), # Disorder
	State('mod_res_toggles', 'value'), # Modified residues
	State('draw_region_dropdown', 'value'), # Composition biased regions

	State('FT_order_list', 'value'), # Order of feature drawing

	prevent_initial_call=True
	)
def download_uniprot_files(
	button_press,

	case,
	sorting,
	focus,
	threshold,

	length_factor,
	height_factor,
	pensize,
	text_size,
	biased_region_text_size,
	protein_length,
	uniform,

	coverage,
	secondary_struct,
	disorder,
	mod_res,
	draw_region,

	FT_order_list
	):

	ctx = dash.callback_context
	triggered = ctx.triggered[0]['prop_id'].split('.')[0]
	print("TRIGERED: "+triggered)

	if case == [] :
		case = None
	if case != None :
		''.join(case)
	if sorting == [] :
		sorting = None
	if sorting != None :
		''.join(sorting)
	if focus == [] :
		focus = None
	if focus != None :
		''.join(focus)
	if threshold == [] :
		threshold = None
	if threshold != None :
		str(threshold)

	if length_factor == [] :
		length_factor = None
	if length_factor != None :
		str(length_factor)
	if height_factor == [] :
		height_factor = None
	if height_factor != None :
		str(height_factor)
	if pensize == [] :
		pensize = None
	if pensize != None :
		str(pensize)
	if text_size == [] :
		text_size = None
	if text_size != None :
		str(text_size)
	if biased_region_text_size == [] :
		biased_region_text_size = None
	if biased_region_text_size != None :
		str(biased_region_text_size)
	if protein_length == [] :
		protein_length = None
	if protein_length != None :
		''.join(protein_length)
	if uniform == [] :
		uniform = None
	if uniform != None :
		''.join(uniform)

	if coverage == [] :
		coverage = None
	if coverage != None :
		''.join(coverage)
	if secondary_struct == [] :
		secondary_struct = None
	if secondary_struct != None :
		''.join(secondary_struct)
	if disorder == [] :
		disorder = None
	if disorder != None :
		''.join(disorder)
	if mod_res == [] :
		mod_res = None
	if mod_res != None :
		''.join(mod_res)
	if draw_region == [] :
		draw_region = None
	if draw_region != None :
		''.join(draw_region)

	param_dict = {
	'Case to draw' : case,
	'Sorting' : sorting,
	'Focus on' : focus,
	'Threshold' : threshold,

	'Horizontal stretch factor' : length_factor,
	'Vertical stretch factor' : height_factor,
	'Protein thickness' : pensize,
	'Protein name size' : text_size,
	'Biased regions text size' : biased_region_text_size,
	'Show protein length' : protein_length,
	'Uniform shape fill/contour' : uniform,

	'3D structure coverage' : coverage,
	'Secondary structure' : secondary_struct,
	'Disorder' : disorder,
	'Modified residues' : mod_res,
	'Composition biased regions' : draw_region,

	'Order of feature drawing' : ', '.join(FT_order_list)
	}

	# print(param_dict)

	return param_dict

@app.callback(
	Output('user_guide_download', 'data'),
	Input('user_guide_btn', 'n_clicks'),
	prevent_initial_call=True,
	)
def update_output(placeholder):
	return dcc.send_file(os.path.join('User_guide_v1.0.0'+'.pdf'))


# Protein list example download
@app.callback(
	Output('protein_list_example_download', 'data'),
	Input('protein_list_example_download_button', 'n_clicks'),
	prevent_initial_call=True,
	)
def update_output(placeholder):
	return dcc.send_file(os.path.join('Examples', 'TRIO_protein_list'+'.xlsx'))

@app.callback(
	Output('protein_list_example_download_2', 'data'),
	Input('protein_list_example_download_button_2', 'n_clicks'),
	prevent_initial_call=True,
	)
def update_output(placeholder):
	return dcc.send_file(os.path.join('Examples', 'TRIO_protein_list'+'.xlsx'))

# Modifications example download
@app.callback(
	Output('modifications_example_download', 'data'),
	Input('modifications_example_download_button', 'n_clicks'),
	prevent_initial_call=True,
	)
def update_output(placeholder):
	return dcc.send_file(os.path.join('Examples', 'TRIO_modifications'+'.xlsx'))

@app.callback(
	Output('modifications_example_download_2', 'data'),
	Input('modifications_example_download_button_2', 'n_clicks'),
	prevent_initial_call=True,
	)
def update_output(placeholder):
	return dcc.send_file(os.path.join('Examples', 'TRIO_modifications'+'.xlsx'))

# Numerical values example downloads
@app.callback(
	Output('values_example_download', 'data'),
	Input('values_example_download_button', 'n_clicks'),
	prevent_initial_call=True,
	)
def update_output(placeholder):
	return dcc.send_file(os.path.join('Examples', 'TRIO_numerical_values'+'.xlsx'))

@app.callback(
	Output('values_example_download_2', 'data'),
	Input('values_example_download_button_2', 'n_clicks'),
	prevent_initial_call=True,
	)
def update_output(placeholder):
	return dcc.send_file(os.path.join('Examples', 'TRIO_numerical_values'+'.xlsx'))

# Shapes and colors example downloads
@app.callback(
	Output('shapes_and_colors_example_download', 'data'),
	Input('shapes_and_colors_example_download_button', 'n_clicks'),
	prevent_initial_call=True,
	)
def update_output(placeholder):
	return dcc.send_file(os.path.join('Examples', 'TRIO_shapes_and_colors'+'.xlsx'))

@app.callback(
	Output('shapes_and_colors_example_download_2', 'data'),
	Input('shapes_and_colors_example_download_button_2', 'n_clicks'),
	prevent_initial_call=True,
	)
def update_output(placeholder):
	return dcc.send_file(os.path.join('Examples', 'TRIO_shapes_and_colors'+'.xlsx'))

# Cur regions example downloads
@app.callback(
	Output('cut_regions_example_download', 'data'),
	Input('cut_regions_example_download_button', 'n_clicks'),
	prevent_initial_call=True,
	)
def update_output(placeholder):
	return dcc.send_file(os.path.join('Examples', 'TRIO_cut_regions'+'.xlsx'))

@app.callback(
	Output('cut_regions_example_download_2', 'data'),
	Input('cut_regions_example_download_button_2', 'n_clicks'),
	prevent_initial_call=True,
	)
def update_output(placeholder):
	return dcc.send_file(os.path.join('Examples', 'TRIO_cut_regions'+'.xlsx'))


@app.callback(
	Output('uniprot_end_state', 'src'),
	Input('uniprot_start_state', 'children'),
	)
def update_output(placeholder):
	print('Checking Uniprot availability')
	obsolete_data, unrecognized_data, show_obsolete, state_src, missing_protein_names, protein_list_df = proteinDataGathering(['P04637'], dl_latest=True, uniprot_availability_test=True)
	return state_src

# Protein list datatable management
@app.callback(


	Output('protein_list_datatable_div', 'children'),
	Output('protein_list_datatable', 'data'),
	Output('step_1_height', 'style'),
	Output('Step_1', 'data'),
	Output('protein_list_dnd_div', 'children'),
	Output('protein_list_download', 'data'),
	Output('duplicate_entry_names_div', 'children'),
	Output('step_1_save_img', 'children'),

	Output('step_1_placeholder', 'children'),
	Output('obsolete_list', 'data'),
	Output('unrecognized_list', 'data'),
	Output('obsolete_code_div', 'style'),
	Output('step_1_end_state', 'src'),

	Input('Step_1', 'data'),
	Input('protein_list_editing_rows_btn', 'n_clicks'),
	Input('protein_list_clear', 'n_clicks'),
	Input('protein_list_dnd', 'contents'),
	Input('protein_list_changes_save', 'n_clicks'),
	Input('Step_1', 'clear_data'),
	Input('protein_list_download_button', 'n_clicks'),
	Input('remove_organism_button', 'n_clicks'),
	Input('protein_list_datatable_div', 'n_clicks'),

	Input('step_1_start_state', 'children'),

	State('protein_list_datatable', 'data'),
	State('protein_list_datatable', 'columns'),
	State('protein_list_datatable', 'data'),
	State('protein_list_dnd', 'filename'),
	State('Step_1', 'data'),

	State('dl_latest_toggle', 'value'),
)
def update_output(
	data,
	n_clicks,
	clear_btn,
	content,
	save_changes,
	clear,
	download_click,
	remove_organism,
	datatable_click,

	button_press,

	rows,
	columns,
	protein_list_data,
	filename,
	protein_list_stored,

	dl_latest
	):
	ctx = dash.callback_context
	triggered = ctx.triggered[0]['prop_id'].split('.')[0]
	print("TRIGERED: "+triggered)

	data = protein_list_stored

	# print('-----------------------')
	# print('data', data)
	# print('rows', rows)

	children = None
	df = pd.DataFrame()

	saved = [html.Img(src=b64_image(os.path.join('GUI', 'save'+'.png')), style={'width':'4vh', 'height':'4vh'})]
	need_save = [html.Img(src=b64_image(os.path.join('GUI', 'need_save'+'.png')), style={'width':'4vh', 'height':'4vh'})]

	if triggered == 'protein_list_editing_rows_btn' :
		if n_clicks > 0:
			rows.append({c['id']: '' for c in columns})

			return no_update, rows, {'width':'100%', 'height':str(int(min(len(rows),25)*1.7)+30)+'vh'}, no_update, no_update, no_update, no_update, need_save, no_update, no_update, no_update, no_update, no_update

	elif triggered == 'Step_1' :

		if data == None :

			children = html.Div(dash_table.DataTable(id='protein_list_datatable'))

		else :
			df = pd.read_json(data, orient='split')

			children = html.Div(
				dash_table.DataTable(
					id='protein_list_datatable',
					data=df.to_dict('records'),
					columns=[{'name': i, 'id': i} for i in df.columns],
					style_header={
						'backgroundColor': 'rgb(125, 111, 134)',
						'fontWeight': 'bold',
						'color': 'white'
					},
					style_cell={
						'backgroundColor': 'rgb(34, 34, 34)',
						'color': 'white',
						'minWidth': '150px', 'width': '150px', 'maxWidth': '150px',
						'overflow': 'hidden',
						'textOverflow': 'ellipsis',
					},
					page_size=20,
					editable=True,
					row_deletable=True,
					fixed_columns={'headers': True, 'data': 2},
					style_table={'minWidth': '100%', 'maxWidth': '100%'},
				),
			)

		return children, no_update, {'width':'100%', 'height':str(int(min(len(df.index),25)*1.7)+30)+'vh'}, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update

	elif triggered == 'protein_list_clear' :

		column_names = ['code', 'protein']

		data_table_data = pd.DataFrame(columns=column_names)

		header_color = 'rgb(125, 111, 134)'
		header_font_color = 'white'

		children = html.Div(
			dash_table.DataTable(
				id='protein_list_datatable',
				data=[{elem:'' for elem in column_names}],
				columns=[{'name': i, 'id': i} for i in data_table_data.columns],
				style_header={
					'backgroundColor': header_color,
					'fontWeight': 'bold',
					'color': header_font_color
				},
				style_cell={
					'backgroundColor': 'rgb(34,34,34)',
					'color': 'white',
					'minWidth': '150px', 'width': '150px', 'maxWidth': '150px',
					'overflow': 'hidden',
					'textOverflow': 'ellipsis',
				},
				page_size=20,
				editable=True,
				row_deletable=True,
				fixed_columns={'headers': True, 'data': 2},
				style_table={'minWidth': '100%', 'maxWidth': '100%'},
			),
		)

		return children, no_update, {'width':'100%', 'height':str(int(2*1.7)+30)+'vh'}, no_update, no_update, no_update, no_update, need_save, no_update, no_update, no_update, no_update, no_update

	elif triggered == 'protein_list_changes_save' :

		modified_df = pd.DataFrame.from_dict(protein_list_data)

		if len(list(modified_df['protein'])) != len(list(set(list(modified_df['protein'])))) :
			# if duplicate names
			children = [
				html.Img(src=b64_image(os.path.join('GUI', 'warning'+'.png')), style={'width':'6vh', 'height':'6vh'}),
				html.H4('Duplicate Entry names detected ! Please remove duplicates.', style={'color': '#a50104'})
			]

			return no_update, no_update, no_update, modified_df.to_json(date_format='iso', orient='split'), no_update, no_update, children, saved, no_update, no_update, no_update, no_update, no_update

		else :
			return no_update, no_update, no_update, modified_df.to_json(date_format='iso', orient='split'), no_update, no_update, [], saved, no_update, no_update, no_update, no_update, no_update

	elif triggered == 'protein_list_dnd' :
		# Upload of same file twice does not work : https://github.com/plotly/dash-core-components/issues/816
		# Problem solved by regenerating after each upload a new drag and drop Div (de La Rache method)

		children = dcc.Upload(id='protein_list_dnd',children=html.Div(['Drag and Drop or ', html.A('Select Files')]),style={'height': '60px','lineHeight': '60px','borderWidth': '1px','borderStyle': 'dashed','borderRadius': '5px','textAlign': 'center','margin': '10px'},)

		if content is not None:
			complete_protein_list_df = parse_content(content, filename)

			complete_protein_list_df = complete_protein_list_df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

			protein_list_df = pd.DataFrame()
			if 'Entry' in list(complete_protein_list_df) :
				protein_list_df['code'] = complete_protein_list_df['Entry']
			elif 'code' in list(complete_protein_list_df) :
				protein_list_df['code'] = complete_protein_list_df['code']
			if 'Entry name' in list(complete_protein_list_df) :
				protein_list_df['protein'] = complete_protein_list_df['Entry name']
			if 'Entry Name' in list(complete_protein_list_df) :
				protein_list_df['protein'] = complete_protein_list_df['Entry Name']
			elif 'protein' in list(complete_protein_list_df) :
				protein_list_df['protein'] = complete_protein_list_df['protein']

			print(protein_list_df)
			if 'code' not in list(protein_list_df) :
				# Assumes that they are only Uniprot codes in the list
				print('no code column name found')
				first_entry = list(complete_protein_list_df)[0]
				other_entries = list(complete_protein_list_df[first_entry])

				entry_list = [first_entry] + other_entries
				entry_list = [elem.strip() for elem in entry_list if elem.strip() != '']

				protein_list_df = pd.DataFrame()
				protein_list_df['code'] = entry_list

			datatable_children = html.Div(
				dash_table.DataTable(
					id='protein_list_datatable',
					data=protein_list_df.to_dict('records'),
					columns=[{'name': i, 'id': i} for i in protein_list_df.columns],
					style_header={
						'backgroundColor': 'rgb(125, 111, 134)',
						'fontWeight': 'bold',
						'color': 'white'
					},
					style_cell={
						'backgroundColor': 'rgb(34, 34, 34)',
						'color': 'white',
						'minWidth': '150px', 'width': '150px', 'maxWidth': '150px',
						'overflow': 'hidden',
						'textOverflow': 'ellipsis',
					},
					page_size=20,
					editable=True,
					row_deletable=True,
					fixed_columns={'headers': True, 'data': 2},
					style_table={'minWidth': '100%', 'maxWidth': '100%'},
				),
			)

			if 'protein' in list(protein_list_df) and  len(list(protein_list_df['protein'])) != len(list(set(list(protein_list_df['protein'])))) :
				duplicate_children = [
					html.Img(src=b64_image(os.path.join('GUI', 'warning'+'.png')), style={'width':'6vh', 'height':'6vh'}),
					html.H4('Duplicate Entry names detected ! Please remove duplicates.', style={'color': '#a50104'})
				]
				return datatable_children, no_update, no_update, protein_list_df.to_json(date_format='iso', orient='split'), children, no_update, duplicate_children, no_update, no_update, no_update, no_update, no_update, no_update
			else :
				return datatable_children, no_update, no_update, protein_list_df.to_json(date_format='iso', orient='split'), children, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update
		else :
			return no_update, no_update, no_update, no_update, children, no_update, [], no_update, no_update, no_update, no_update, no_update, no_update

	elif triggered == 'protein_list_download_button' :

		if protein_list_stored == None or protein_list_stored == '{"columns":[],"index":[],"data":[]}' :
			df = pd.DataFrame(columns=['code', 'protein'])
		else :
			df = pd.read_json(protein_list_stored, orient='split')
			df = df[['code', 'protein']]

		return no_update, no_update, no_update, no_update, no_update, send_data_frame(df.to_excel, filename='protein_list.xlsx', index=False), no_update, no_update, no_update, no_update, no_update, no_update, no_update

	elif triggered == 'remove_organism_button' :

		df = pd.read_json(data, orient='split')
		df = df[['code', 'protein']]

		entry_name_list = list(df['protein'])
		organism_list = [entry.split('_')[-1] for entry in entry_name_list]
		if len(list(set(organism_list))) == 1 :
			name_list = ['_'.join(entry.split('_')[0:-1]) for entry in entry_name_list]
			df['protein'] = name_list

			datatable_children = html.Div(
				dash_table.DataTable(
					id='protein_list_datatable',
					data=df.to_dict('records'),
					columns=[{'name': i, 'id': i} for i in df.columns],
					style_header={
						'backgroundColor': 'rgb(125, 111, 134)',
						'fontWeight': 'bold',
						'color': 'white'
					},
					style_cell={
						'backgroundColor': 'rgb(34, 34, 34)',
						'color': 'white',
						'minWidth': '150px', 'width': '150px', 'maxWidth': '150px',
						'overflow': 'hidden',
						'textOverflow': 'ellipsis',
					},
					page_size=20,
					editable=True,
					row_deletable=True,
					fixed_columns={'headers': True, 'data': 2},
					style_table={'minWidth': '100%', 'maxWidth': '100%'},
				),
			)

			return datatable_children, no_update, no_update, no_update, no_update, no_update, no_update, need_save, no_update, no_update, no_update, no_update, no_update

		else :
			print('Not same organism for each entry')

			return no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update

	elif triggered == 'step_1_start_state' :
		print('step_1 start')
		start_time = time.time()

		complete_protein_list_df = pd.read_json(data, orient='split')

		obsolete_data, unrecognized_data, show_obsolete, state_src, missing_protein_names, protein_list_df = proteinDataGathering(complete_protein_list_df, dl_latest=dl_latest)

		end_time = time.time()
		print('Step 1 execution time:', str(end_time - start_time))

		if missing_protein_names :
			datatable_children = html.Div(
				dash_table.DataTable(
					id='protein_list_datatable',
					data=protein_list_df.to_dict('records'),
					columns=[{'name': i, 'id': i} for i in protein_list_df.columns],
					style_header={
						'backgroundColor': 'rgb(125, 111, 134)',
						'fontWeight': 'bold',
						'color': 'white'
					},
					style_cell={
						'backgroundColor': 'rgb(34, 34, 34)',
						'color': 'white',
						'minWidth': '150px', 'width': '150px', 'maxWidth': '150px',
						'overflow': 'hidden',
						'textOverflow': 'ellipsis',
					},
					page_size=20,
					editable=True,
					row_deletable=True,
					fixed_columns={'headers': True, 'data': 2},
					style_table={'minWidth': '100%', 'maxWidth': '100%'},
				),
			)

			return datatable_children, no_update, no_update, protein_list_df.to_json(date_format='iso', orient='split'), no_update, no_update, no_update, no_update, no_update, obsolete_data, unrecognized_data, show_obsolete, state_src
		else :
			return no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, obsolete_data, unrecognized_data, show_obsolete, state_src

	elif datatable_click != None and datatable_click > 0 :
		return no_update, no_update, no_update, no_update, no_update, no_update, no_update, need_save, no_update, no_update, no_update, no_update, no_update

	else :
		return no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update

# Step 1 enter
@app.callback(
	Output('step_1_start_state', 'children'),
	Input('step_1_launch_button', 'n_clicks'),
	prevent_initial_call=True
	)
def download_uniprot_files(button_press):
	return html.Img(id='step_1_end_state', src=b64_image(os.path.join('GUI', 'fail'+'.png')), style={'width':'4vh', 'height':'4vh'})

# Exception datatable management
@app.callback(
	Output('exception_datatable_div', 'children'),
	Output('exception_datatable', 'data'),
	Output('step_2_height', 'style'),
	Output('Step_2_exceptions', 'data'),
	Output('exception_dnd_div', 'children'),
	Output('exception_download', 'data'),
	Output('step_2_save_img', 'children'),

	Input('Step_2_exceptions', 'data'),
	Input('exception_editing_rows_btn', 'n_clicks'),
	Input('exception_clear', 'n_clicks'),
	Input('exception_dnd', 'contents'),
	Input('exception_changes_save', 'n_clicks'),
	Input('Step_2_exceptions', 'clear_data'),
	Input('exception_download_button', 'n_clicks'),
	Input('exception_datatable_div', 'n_clicks'),

	State('exception_datatable', 'data'),
	State('exception_datatable', 'columns'),
	State('exception_datatable', 'data'),
	State('exception_dnd', 'filename'),
	State('Step_2_exceptions', 'data'),
)
def update_output(
	data,
	n_clicks,
	clear_btn,
	exception_content,
	exception_save_changes,
	exception_clear,
	exception_download_click,
	datatable_click,
	rows,
	columns,
	exception_data,
	exception_filename,
	exception_stored
	):
	ctx = dash.callback_context
	triggered = ctx.triggered[0]['prop_id'].split('.')[0]
	print("TRIGERED: "+triggered)

	# print('-----------------------')
	# print('data', data)
	# print('rows', rows)

	children = None
	df = pd.DataFrame()

	saved = [html.Img(src=b64_image(os.path.join('GUI', 'save'+'.png')), style={'width':'4vh', 'height':'4vh'})]
	need_save = [html.Img(src=b64_image(os.path.join('GUI', 'need_save'+'.png')), style={'width':'4vh', 'height':'4vh'})]

	if triggered == 'exception_editing_rows_btn' :
		if n_clicks > 0:
			rows.append({c['id']: '' for c in columns})

			return no_update, rows, {'width':'100%', 'height':str(int(min(len(rows),25)*1.7)+30)+'vh'}, no_update, no_update, no_update, need_save

	elif triggered == 'Step_2_exceptions' :

		if data == None :

			children = html.Div(dash_table.DataTable(id='exception_datatable'))

		else :
			df = pd.read_json(data, orient='split')
			df = df[['ex_type', 'protein', 'feature_type', 'feature', 'start', 'length']]

			children = html.Div(
				dash_table.DataTable(
					id='exception_datatable',
					data=df.to_dict('records'),
					columns=[{'name': i, 'id': i} for i in df.columns],
					style_header={
						'backgroundColor': 'rgb(231, 215, 193)',
						'fontWeight': 'bold',
						'color': 'black'
					},
					style_cell={
						'backgroundColor': 'rgb(34, 34, 34)',
						'color': 'white',
						'minWidth': '150px', 'width': '150px', 'maxWidth': '150px',
						'overflow': 'hidden',
						'textOverflow': 'ellipsis',
					},
					page_size=20,
					editable=True,
					row_deletable=True,
					fixed_columns={'headers': True, 'data': 0},
					style_table={'minWidth': '100%', 'maxWidth': '100%'},
				),
			)

		return children, no_update, {'width':'100%', 'height':str(int(min(len(df.index),25)*1.7)+30)+'vh'}, no_update, no_update, no_update, no_update

	elif triggered == 'exception_clear' :

		column_names = ['ex_type', 'protein', 'feature_type', 'feature', 'start', 'length']

		data_table_data = pd.DataFrame(columns=column_names)

		header_color = 'rgb(231, 215, 193)'
		header_font_color = 'black'

		children = html.Div(
			dash_table.DataTable(
				id='exception_datatable',
				data=[{elem:'' for elem in column_names}],
				columns=[{'name': i, 'id': i} for i in data_table_data.columns],
				style_header={
					'backgroundColor': header_color,
					'fontWeight': 'bold',
					'color': header_font_color
				},
				style_cell={
					'backgroundColor': 'rgb(34,34,34)',
					'color': 'white',
					'minWidth': '150px', 'width': '150px', 'maxWidth': '150px',
					'overflow': 'hidden',
					'textOverflow': 'ellipsis',
				},
				page_size=20,
				editable=True,
				row_deletable=True,
				fixed_columns={'headers': True, 'data': 0},
				style_table={'minWidth': '100%', 'maxWidth': '100%'},
			),
		)

		return children, no_update, {'width':'100%', 'height':str(int(2*1.7)+30)+'vh'}, no_update, no_update, no_update, need_save

	elif triggered == 'exception_changes_save' :

		df = pd.DataFrame.from_dict(exception_data)

		return no_update, no_update, no_update, df.to_json(date_format='iso', orient='split'), no_update, no_update, saved

	elif triggered == 'exception_dnd' :
		# Upload of same file twice does not work : https://github.com/plotly/dash-core-components/issues/816
		# Problem solved by regenerating after each upload a new drag and drop Div (de La Rache method)

		children = dcc.Upload(id='exception_dnd',children=html.Div(['Drag and Drop or ', html.A('Select Files')]),style={'height': '60px','lineHeight': '60px','borderWidth': '1px','borderStyle': 'dashed','borderRadius': '5px','textAlign': 'center','margin': '10px'},)

		if exception_content is not None:
			# print(exception_filename)
			df = parse_content(exception_content, exception_filename)
			df = df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

			datatable_children = html.Div(
				dash_table.DataTable(
					id='exception_datatable',
					data=df.to_dict('records'),
					columns=[{'name': i, 'id': i} for i in df.columns],
					style_header={
						'backgroundColor': 'rgb(231, 215, 193)',
						'fontWeight': 'bold',
						'color': 'black'
					},
					style_cell={
						'backgroundColor': 'rgb(34, 34, 34)',
						'color': 'white',
						'minWidth': '150px', 'width': '150px', 'maxWidth': '150px',
						'overflow': 'hidden',
						'textOverflow': 'ellipsis',
					},
					page_size=20,
					editable=True,
					row_deletable=True,
					fixed_columns={'headers': True, 'data': 0},
					style_table={'minWidth': '100%', 'maxWidth': '100%'},
				),
			)

			return datatable_children, no_update, no_update,df.to_json(date_format='iso', orient='split'), children, no_update, no_update
		else :
			return datatable_children, no_update, no_update,no_update, children, no_update, no_update

	elif triggered == 'exception_download_button' :

		if exception_stored == None or exception_stored == '{"columns":[],"index":[],"data":[]}' :
			df = pd.DataFrame(columns=['ex_type', 'protein', 'feature_type', 'feature', 'start', 'length'])
		else :
			df = pd.read_json(exception_stored, orient='split')
			df = df[['ex_type', 'protein', 'feature_type', 'feature', 'start', 'length']]

		return no_update, no_update, no_update,no_update, no_update, send_data_frame(df.to_excel, filename='modifications.xlsx', index=False), no_update

	elif datatable_click != None and datatable_click > 0 :
		return no_update, no_update, no_update, no_update, no_update, no_update, need_save

	else :
		return no_update, no_update, no_update, no_update, no_update, no_update, no_update

# Values datatable management
@app.callback(
	Output('values_datatable_div', 'children'),
	Output('values_datatable', 'data'),
	Output('step_3_height', 'style'),
	Output('Step_3_db_user', 'data'),
	Output('values_dnd_div', 'children'),
	Output('values_download', 'data'),
	Output('step_3_save_img', 'children'),

	Input('Step_3_db_user', 'data'),
	Input('values_editing_rows_btn', 'n_clicks'),
	Input('values_clear', 'n_clicks'),
	Input('values_dnd', 'contents'),
	Input('values_changes_save', 'n_clicks'),
	Input('Step_3_db_user', 'clear_data'),
	Input('values_download_button', 'n_clicks'),
	Input('values_datatable_div', 'n_clicks'),

	State('values_datatable', 'data'),
	State('values_datatable', 'columns'),
	State('values_datatable', 'data'),
	State('values_dnd', 'filename'),
	State('Step_3_db_user', 'data'),
)
def update_output(
	data,
	n_clicks,
	clear_btn,
	values_content,
	values_save_changes,
	values_clear,
	values_download_click,
	datatable_click,
	rows,
	columns,
	values_data,
	values_filename,
	values_stored
	):
	ctx = dash.callback_context
	triggered = ctx.triggered[0]['prop_id'].split('.')[0]
	print("TRIGERED: "+triggered)

	# print('-----------------------')
	# print('data', data)
	# print('rows', rows)

	children = None
	df = pd.DataFrame()

	saved = [html.Img(src=b64_image(os.path.join('GUI', 'save'+'.png')), style={'width':'4vh', 'height':'4vh'})]
	need_save = [html.Img(src=b64_image(os.path.join('GUI', 'need_save'+'.png')), style={'width':'4vh', 'height':'4vh'})]

	if triggered == 'values_editing_rows_btn' :
		if n_clicks > 0:
			rows.append({c['id']: '' for c in columns})

			return no_update, rows, {'width':'100%', 'height':str(int(min(len(rows),25)*1.7)+30)+'vh'}, no_update, no_update, no_update, need_save

	elif triggered == 'Step_3_db_user' :

		if data == None :

			children = html.Div(dash_table.DataTable(id='values_datatable'))

		else :
			df = pd.read_json(data, orient='split')
			# df = df[['ex_type', 'protein', 'feature_type', 'feature', 'start', 'length']]

			children = html.Div(
				dash_table.DataTable(
					id='values_datatable',
					data=df.to_dict('records'),
					columns=[{'name': i, 'id': i} for i in df.columns],
					style_header={
						'backgroundColor': 'rgb(223, 170, 24)',
						'fontWeight': 'bold',
						'color': 'black'
					},
					style_cell={
						'backgroundColor': 'rgb(34, 34, 34)',
						'color': 'white',
						'minWidth': '150px', 'width': '150px', 'maxWidth': '150px',
						'overflow': 'hidden',
						'textOverflow': 'ellipsis',
					},
					page_size=20,
					editable=True,
					row_deletable=True,
					fixed_columns={'headers': True, 'data': 0},
					style_table={'minWidth': '100%', 'maxWidth': '100%'},
				),
			)

		return children, no_update, {'width':'100%', 'height':str(int(min(len(df.index),25)*1.7)+30)+'vh'}, no_update, no_update, no_update, no_update

	elif triggered == 'values_clear' :

		column_names = ['protein', 'feature', 'start']

		data_table_data = pd.DataFrame(columns=column_names)

		header_color = 'rgb(223, 170, 24)'
		header_font_color = 'black'

		children = html.Div(
			dash_table.DataTable(
				id='values_datatable',
				data=[{elem:'' for elem in column_names}],
				columns=[{'name': i, 'id': i} for i in data_table_data.columns],
				style_header={
					'backgroundColor': header_color,
					'fontWeight': 'bold',
					'color': header_font_color
				},
				style_cell={
					'backgroundColor': 'rgb(34,34,34)',
					'color': 'white',
					'minWidth': '150px', 'width': '150px', 'maxWidth': '150px',
					'overflow': 'hidden',
					'textOverflow': 'ellipsis',
				},
				page_size=20,
				editable=True,
				row_deletable=True,
				fixed_columns={'headers': True, 'data': 0},
				style_table={'minWidth': '100%', 'maxWidth': '100%'},
			),
		)

		return children, no_update, {'width':'100%', 'height':str(int(2*1.7)+30)+'vh'}, no_update, no_update, no_update, need_save

	elif triggered == 'values_changes_save' :

		df = pd.DataFrame.from_dict(values_data)

		return no_update, no_update, no_update, df.to_json(date_format='iso', orient='split'), no_update, no_update, saved

	elif triggered == 'values_dnd' :
		# Upload of same file twice does not work : https://github.com/plotly/dash-core-components/issues/816
		# Problem solved by regenerating after each upload a new drag and drop Div (de La Rache method)

		children = dcc.Upload(id='values_dnd',children=html.Div(['Drag and Drop or ', html.A('Select Files')]),style={'height': '60px','lineHeight': '60px','borderWidth': '1px','borderStyle': 'dashed','borderRadius': '5px','textAlign': 'center','margin': '10px'},)

		if values_content is not None:
			# print(values_filename)
			df = parse_content(values_content, values_filename)
			df = df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

			datatable_children = html.Div(
				dash_table.DataTable(
					id='values_datatable',
					data=df.to_dict('records'),
					columns=[{'name': i, 'id': i} for i in df.columns],
					style_header={
						'backgroundColor': 'rgb(223, 170, 24)',
						'fontWeight': 'bold',
						'color': 'black'
					},
					style_cell={
						'backgroundColor': 'rgb(34, 34, 34)',
						'color': 'white',
						'minWidth': '150px', 'width': '150px', 'maxWidth': '150px',
						'overflow': 'hidden',
						'textOverflow': 'ellipsis',
					},
					page_size=20,
					editable=True,
					row_deletable=True,
					fixed_columns={'headers': True, 'data': 0},
					style_table={'minWidth': '100%', 'maxWidth': '100%'},
				),
			)

			return datatable_children, no_update, no_update,df.to_json(date_format='iso', orient='split'), children, no_update, no_update
		else :
			return datatable_children, no_update, no_update,no_update, children, no_update, no_update

	elif triggered == 'values_download_button' :

		if values_stored == None or values_stored == '{"columns":[],"index":[],"data":[]}' :
			df = pd.DataFrame(columns=['protein', 'feature', 'start'])
		else :
			df = pd.read_json(values_stored, orient='split')
			# df = df[['ex_type', 'protein', 'feature_type', 'feature', 'start', 'length']]

		return no_update, no_update, no_update,no_update, no_update, send_data_frame(df.to_excel, filename='values.xlsx', index=False), no_update

	elif datatable_click != None and datatable_click > 0 :
		return no_update, no_update, no_update, no_update, no_update, no_update, need_save

	else :
		return no_update, no_update, no_update, no_update, no_update, no_update, no_update


# Feature shape datatable management
@app.callback(
	Output('feature_shape_datatable_div', 'children'),
	Output('feature_shape_datatable', 'data'),
	Output('step_4_height', 'style'),
	Output('Step_3_feature_shape', 'data'),
	Output('feature_shape_dnd_div', 'children'),
	Output('feature_shape_download', 'data'),
	Output('step_4_protein_shape_save_img', 'children'),
	Output('last_seed', 'data'),
	Output('auto_feature_select_img', 'src'),
	Output('auto_feature_select_text', 'style'),

	Input('Step_3_feature_shape', 'data'),
	Input('feature_shape_editing_rows_btn', 'n_clicks'),
	Input('feature_shape_clear', 'n_clicks'),
	Input('auto_feature_select', 'n_clicks'),
	Input('feature_shape_dnd', 'contents'),
	Input('feature_shape_changes_save', 'n_clicks'),
	Input('Step_3_feature_shape', 'clear_data'),
	Input('feature_shape_download_button', 'n_clicks'),
	Input('feature_shape_datatable_div', 'n_clicks'),
	Input('default_feat_param', 'n_clicks'),

	State('feature_shape_datatable', 'data'),
	State('feature_shape_datatable', 'columns'),
	State('Step_4_feature_occurrence', 'data'),
	State('auto_toggle_2', 'value'),
	State('auto_select_threshold_2', 'value'),
	State('unique_entry_nb', 'children'),
	State('feature_shape_datatable', 'data'),
	State('feature_shape_dnd', 'filename'),
	State('Step_3_feature_shape', 'data'),

	State('height', 'value'),

	State('coverage_toggles', 'value'),
	State('2nd_struct_toggles', 'value'),
	State('disorder_toggles', 'value'),
	State('mod_res_toggles', 'value'),
	State('draw_region_dropdown', 'value'),

	State('last_seed', 'data'),
	State('lock_seed', 'value')
	)
def update_output(
	data,
	n_clicks,
	clear_btn,
	autofeatchoice,
	feature_shape_content,
	feature_shape_save_changes,
	feature_shape_clear,
	feature_shape_download_click,
	datatable_click,
	default_feat_param_click,
	rows,
	columns,
	feature_occurrence_data,
	auto_toggle,
	auto_select_threshold,
	protein_nb,
	feature_shape_data,
	feature_shape_filename,
	feature_shape_stored,

	height,
	draw_coverage,
	draw_2nd_struct,
	draw_disorder,
	draw_mod_res,
	draw_region,

	last_seed,
	lock_seed
	):
	ctx = dash.callback_context
	triggered = ctx.triggered[0]['prop_id'].split('.')[0]
	print("TRIGERED: "+triggered)

	children = None
	df = pd.DataFrame()

	saved = [html.Img(src=b64_image(os.path.join('GUI', 'save'+'.png')), style={'width':'4vh', 'height':'4vh'})]
	need_save = [html.Img(src=b64_image(os.path.join('GUI', 'need_save'+'.png')), style={'width':'4vh', 'height':'4vh'})]

	autoed = b64_image(os.path.join('GUI', 'auto'+'.png'))
	need_auto = b64_image(os.path.join('GUI', 'need_auto'+'.png'))

	autoed_text = {'text-align':'center', 'margin-top':'1vh', 'background-color':'rgb(34, 34, 34)', 'color':'white', 'font-family':'courier', 'font-weight':'bold'}
	need_auto_text = {'text-align':'center', 'margin-top':'1vh', 'background-color':'rgb(223, 170, 24)', 'color':'black', 'font-family':'courier', 'font-weight':'bold'}

	if triggered == 'feature_shape_editing_rows_btn' :
		if n_clicks > 0:
			rows.append({c['id']: '' for c in columns})

			return no_update, rows, {'width':'100%', 'height':str(int(min(len(rows),25)*1.7)+30)+'vh'}, no_update, no_update, no_update, need_save, no_update, no_update, no_update

	elif triggered == 'Step_3_feature_shape' :
		if data == None :

			children = html.Div(dash_table.DataTable(id='feature_shape_datatable'))
			autoed = need_auto
			autoed_text = need_auto_text

		else :
			df = pd.read_json(data, orient='split')
			df = df[['feature', 'shape', 'orientation', 'height', 'contour_color', 'contour_colormap', 'contour_threshold', 'color', 'colormap', 'threshold', 'pensize']]

			if sum([1 for elem in list(df['feature']) if elem != None and elem != '']) == 0 or sum([1 for elem in list(df['shape']) if elem != None and elem != '']) == 0: 
				autoed = need_auto
				autoed_text = need_auto_text

			children = html.Div(
				dash_table.DataTable(
					id='feature_shape_datatable',
					data=df.to_dict('records'),
					columns=[{'name': i, 'id': i} for i in df.columns],
					style_header={
						'backgroundColor': 'rgba(148,163,188)',
						'fontWeight': 'bold',
						'color': 'black'
					},
					style_cell={
						'backgroundColor': 'rgb(34, 34, 34)',
						'color': 'white',
						'minWidth': '150px', 'width': '150px', 'maxWidth': '150px',
						'overflow': 'hidden',
						'textOverflow': 'ellipsis',
					},
					page_size=20,
					editable=True,
					row_deletable=True,
					fixed_columns={'headers': True, 'data': 1},
					style_table={'minWidth': '100%', 'maxWidth': '100%'},
				),
			)

		return children, no_update, {'width':'100%', 'height':str(int(min(len(df.index),25)*1.7)+30)+'vh'}, no_update, no_update, no_update, no_update, no_update, autoed, autoed_text

	elif triggered == 'feature_shape_clear' :

		column_names = ['feature', 'shape', 'orientation', 'height', 'contour_color', 'contour_colormap', 'contour_threshold', 'color', 'colormap', 'threshold', 'pensize']

		data_table_data = pd.DataFrame(columns=column_names)

		children = html.Div(
			dash_table.DataTable(
				id='feature_shape_datatable',
				data=[{elem:'' for elem in column_names}],
				columns=[{'name': i, 'id': i} for i in data_table_data.columns],
				style_header={
					'backgroundColor': 'rgba(148,163,188)',
					'fontWeight': 'bold',
					'color': 'black'
				},
				style_cell={
					'backgroundColor': 'rgb(34,34,34)',
					'color': 'white',
					'minWidth': '150px', 'width': '150px', 'maxWidth': '150px',
					'overflow': 'hidden',
					'textOverflow': 'ellipsis',
				},
				page_size=20,
				editable=True,
				row_deletable=True,
				fixed_columns={'headers': True, 'data': 1},
				style_table={'minWidth': '100%', 'maxWidth': '100%'},
			),
		)

		return children, no_update, {'width':'100%', 'height':str(int(2*1.7)+30)+'vh'}, no_update, no_update, no_update, need_save, no_update, need_auto, need_auto_text

	elif triggered == 'auto_feature_select' :
		print('AUTO FEAT CHOICE')
		####feature_occurrence_data

		occurrence_df = pd.read_json(feature_occurrence_data, orient='split')
		occurrence_df = occurrence_df[['feature_type', 'feature', 'occurrence']]

		if lock_seed == ['lock_seed'] :
			new_seed = last_seed
		else :
			new_seed = rd.randint(1,999999)
		rd.seed(new_seed)

		print(new_seed)

		shape_df = autoFeatChoice(occurrence_df, protein_nb, auto_toggle=auto_toggle, auto_select_threshold=auto_select_threshold)
		print('--------------------------')
		# print(shape_df)

		children = html.Div(
			dash_table.DataTable(
				id='feature_shape_datatable',
				data=shape_df.to_dict('records'),
				columns=[{'name': i, 'id': i} for i in shape_df.columns],
				style_header={
					'backgroundColor': 'rgba(148,163,188)',
					'fontWeight': 'bold',
					'color': 'black'
				},
				style_cell={
					'backgroundColor': 'rgb(34, 34, 34)',
					'color': 'white',
					'minWidth': '150px', 'width': '150px', 'maxWidth': '150px',
					'overflow': 'hidden',
					'textOverflow': 'ellipsis',
				},
				page_size=20,
				editable=True,
				row_deletable=True,
				fixed_columns={'headers': True, 'data': 1},
				style_table={'minWidth': '100%', 'maxWidth': '100%'},
			),
		)

		return children, no_update, {'width':'100%', 'height':str(int(min(len(df.index),25)*1.7)+30)+'vh'}, no_update, no_update, no_update, need_save, new_seed, autoed, autoed_text

	elif triggered == 'feature_shape_changes_save' :

		# print(feature_shape_data)
		df = pd.DataFrame.from_dict(feature_shape_data)
		df = df[['feature', 'shape', 'orientation', 'height', 'contour_color', 'contour_colormap', 'contour_threshold', 'color', 'colormap', 'threshold', 'pensize']]
		# print(df)
		if sum([1 for elem in list(df['feature']) if elem != None and elem != '']) == 0 or sum([1 for elem in list(df['shape']) if elem != None and elem != '']) == 0 : 
			autoed = need_auto
			autoed_text = need_auto_text

		return no_update, no_update, no_update, df.to_json(date_format='iso', orient='split'), no_update, no_update, saved, no_update, autoed, autoed_text

	elif triggered == 'feature_shape_dnd' :
		# Upload of same file twice does not work : https://github.com/plotly/dash-core-components/issues/816
		# Problem solved by regenerating after each upload a new drag and drop Div (de La Rache method)

		children = dcc.Upload(id='feature_shape_dnd',children=html.Div(['Drag and Drop or ', html.A('Select Files')]),style={'height': '60px','lineHeight': '60px','borderWidth': '1px','borderStyle': 'dashed','borderRadius': '5px','textAlign': 'center','margin': '10px'},)

		if feature_shape_content is not None:
			# print(feature_shape_filename)
			df = parse_content(feature_shape_content, feature_shape_filename)
			df = df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

			datatable_children = html.Div(
				dash_table.DataTable(
					id='feature_shape_datatable',
					data=df.to_dict('records'),
					columns=[{'name': i, 'id': i} for i in df.columns],
					style_header={
						'backgroundColor': 'rgba(148,163,188)',
						'fontWeight': 'bold',
						'color': 'black'
					},
					style_cell={
						'backgroundColor': 'rgb(34, 34, 34)',
						'color': 'white',
						'minWidth': '150px', 'width': '150px', 'maxWidth': '150px',
						'overflow': 'hidden',
						'textOverflow': 'ellipsis',
					},
					page_size=20,
					editable=True,
					row_deletable=True,
					fixed_columns={'headers': True, 'data': 1},
					style_table={'minWidth': '100%', 'maxWidth': '100%'},
				),
			)

			return datatable_children, no_update, no_update, df.to_json(date_format='iso', orient='split'), children, no_update, no_update, no_update, no_update, no_update
		else :
			return datatable_children, no_update, no_update, no_update, children, no_update, no_update, no_update, no_update, no_update

	elif triggered == 'feature_shape_download_button' :

		if feature_shape_stored == None or feature_shape_stored == '{"columns":[],"index":[],"data":[]}' :
			df = pd.DataFrame(columns=['feature', 'shape', 'orientation', 'height', 'contour_color', 'contour_colormap', 'contour_threshold', 'color', 'colormap', 'threshold', 'pensize'])
		else :
			df = pd.read_json(feature_shape_stored, orient='split')
			df = df[['feature', 'shape', 'orientation', 'height', 'contour_color', 'contour_colormap', 'contour_threshold', 'color', 'colormap', 'threshold', 'pensize']]

		return no_update, no_update, no_update, no_update, no_update, send_data_frame(df.to_excel, filename='shapes_and_colors.xlsx', index=False), no_update, no_update, no_update, no_update

	elif triggered == 'default_feat_param' :

		if feature_shape_stored == None or feature_shape_stored == '{"columns":[],"index":[],"data":[]}' :
			df = pd.DataFrame(columns=['feature', 'shape', 'orientation', 'height', 'contour_color', 'contour_colormap', 'contour_threshold', 'color', 'colormap', 'threshold', 'pensize'])
		else :
			df = pd.read_json(feature_shape_stored, orient='split')
			df = df[['feature', 'shape', 'orientation', 'height', 'contour_color', 'contour_colormap', 'contour_threshold', 'color', 'colormap', 'threshold', 'pensize']]
			df = df.replace({np.nan: None})

		# print(df)

		switch_toggles = draw_coverage + draw_2nd_struct + draw_disorder + draw_mod_res

		draw_coverage = 'False'
		draw_2nd_struct = 'False'
		draw_disorder = 'False'
		draw_mod_res = 'False'

		if 'coverage' in switch_toggles :
			draw_coverage = 'True'
		else :
			draw_coverage = 'False'

		if '2nd_struct' in switch_toggles :
			draw_2nd_struct = 'True'
		else :
			draw_2nd_struct = 'False'

		if 'disorder' in switch_toggles :
			draw_disorder = 'True'
		else :
			draw_disorder = 'False'

		if 'mod_res' in switch_toggles :
			draw_mod_res = 'True'
		else :
			draw_mod_res = 'False'

		df = getDefault(df, height, draw_coverage, draw_2nd_struct, draw_disorder, draw_mod_res, draw_region)

		children = html.Div(
			dash_table.DataTable(
				id='feature_shape_datatable',
				data=df.to_dict('records'),
				columns=[{'name': i, 'id': i} for i in df.columns],
				style_header={
					'backgroundColor': 'rgba(148,163,188)',
					'fontWeight': 'bold',
					'color': 'black'
				},
				style_cell={
					'backgroundColor': 'rgb(34, 34, 34)',
					'color': 'white',
					'minWidth': '150px', 'width': '150px', 'maxWidth': '150px',
					'overflow': 'hidden',
					'textOverflow': 'ellipsis',
				},
				page_size=20,
				editable=True,
				row_deletable=True,
				fixed_columns={'headers': True, 'data': 1},
				style_table={'minWidth': '100%', 'maxWidth': '100%'},
			),
		)

		return children, no_update, {'width':'100%', 'height':str(int(min(len(df.index),25)*1.7)+30)+'vh'}, no_update, no_update, no_update, need_save, no_update, no_update, no_update

	elif datatable_click != None and datatable_click > 0 :
		return no_update, no_update, no_update, no_update, no_update, no_update, need_save, no_update, no_update, no_update

	else :
		return no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update


# Protein cut datatable management
@app.callback(
	Output('protein_cut_datatable_div', 'children'),
	Output('protein_cut_datatable', 'data'),
	Output('Step_3_protein_cut', 'data'),
	Output('protein_cut_dnd_div', 'children'),
	Output('protein_cut_download', 'data'),
	Output('step_4_save_img', 'children'),

	Input('Step_3_protein_cut', 'data'),
	Input('protein_cut_editing_rows_btn', 'n_clicks'),
	Input('protein_cut_clear', 'n_clicks'),
	Input('protein_cut_download_button', 'n_clicks'),
	Input('protein_cut_dnd', 'contents'),
	Input('protein_cut_changes_save', 'n_clicks'),
	Input('Step_3_protein_cut', 'clear_data'),
	Input('protein_cut_datatable_div', 'n_clicks'),

	State('protein_cut_datatable', 'data'),
	State('protein_cut_datatable', 'columns'),
	State('protein_cut_dnd', 'filename'),
	State('Step_3_protein_cut', 'data'),
)
def update_output(
	data,
	n_clicks,
	clear_btn,
	download_btn,
	protein_cut_content,
	protein_cut_save_changes,
	protein_cut_clear,
	datatable_click,
	rows,
	columns,
	protein_cut_filename,
	protein_cut_stored
	):

	ctx = dash.callback_context
	triggered = ctx.triggered[0]['prop_id'].split('.')[0]
	print("TRIGERED_protein_cut: "+triggered)

	# print('-----------------------')
	# print('data', data)
	# print('rows', rows)

	# print('---', datatable_click)

	children = None
	df = pd.DataFrame()

	saved = [html.Img(src=b64_image(os.path.join('GUI', 'save'+'.png')), style={'width':'4vh', 'height':'4vh'})]
	need_save = [html.Img(src=b64_image(os.path.join('GUI', 'need_save'+'.png')), style={'width':'4vh', 'height':'4vh'})]

	if triggered == 'protein_cut_editing_rows_btn' :
		if n_clicks > 0:
			rows.append({c['id']: '' for c in columns})

			return no_update, rows, no_update, no_update, no_update, need_save

	elif triggered == 'Step_3_protein_cut' :

		if data == None :

			children = html.Div(dash_table.DataTable(id='protein_cut_datatable'))

		else :
			df = pd.read_json(data, orient='split')
			df = df[['protein', 'start', 'length']]

			children = html.Div(
				dash_table.DataTable(
					id='protein_cut_datatable',
					data=df.to_dict('records'),
					columns=[{'name': i, 'id': i} for i in df.columns],
					style_header={
						'backgroundColor': 'rgba(148,163,188)',
						'fontWeight': 'bold',
						'color': 'black'
					},
					style_cell={
						'backgroundColor': 'rgba(34,34,34)',
						'color': 'white',
						'minWidth': '150px', 'width': '150px', 'maxWidth': '150px',
						'overflow': 'hidden',
						'textOverflow': 'ellipsis',
					},
					page_size=20,
					editable=True,
					row_deletable=True,
					fixed_columns={'headers': True, 'data': 1},
					style_table={'minWidth': '100%', 'maxWidth': '100%'},
				),
			)

		return children, no_update, no_update, no_update, no_update, no_update

	elif triggered == 'protein_cut_clear' :

		column_names = ['protein', 'start', 'length']

		data_table_data = pd.DataFrame(columns=column_names)

		header_color = 'rgba(148,163,188)'
		header_font_color = 'black'

		children = html.Div(
			dash_table.DataTable(
				id='protein_cut_datatable',
				data=[{elem:'' for elem in column_names}],
				columns=[{'name': i, 'id': i} for i in data_table_data.columns],
				style_header={
					'backgroundColor': header_color,
					'fontWeight': 'bold',
					'color': header_font_color
				},
				style_cell={
					'backgroundColor': 'rgb(34,34,34)',
					'color': 'white',
					'minWidth': '150px', 'width': '150px', 'maxWidth': '150px',
					'overflow': 'hidden',
					'textOverflow': 'ellipsis',
				},
				page_size=20,
				editable=True,
				row_deletable=True,
				fixed_columns={'headers': True, 'data': 1},
				style_table={'minWidth': '100%', 'maxWidth': '100%'},
			),
		)

		return children, no_update, no_update, no_update, no_update, need_save

	elif triggered == 'protein_cut_changes_save' :

		df = pd.DataFrame.from_dict(rows)

		return no_update, no_update, df.to_json(date_format='iso', orient='split'), no_update, no_update, saved

	elif triggered == 'protein_cut_dnd' :
		# Upload of same file twice does not work : https://github.com/plotly/dash-core-components/issues/816
		# Problem solved by regenerating after each upload a new drag and drop Div (de La Rache method)

		children = dcc.Upload(id='protein_cut_dnd',children=html.Div(['Drag and Drop or ', html.A('Select Files')]),style={'height': '60px','lineHeight': '60px','borderWidth': '1px','borderStyle': 'dashed','borderRadius': '5px','textAlign': 'center','margin': '10px'},)

		if protein_cut_content is not None:
			# print(protein_cut_filename)
			df = parse_content(protein_cut_content, protein_cut_filename)
			df = df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

			datatable_children = html.Div(
				dash_table.DataTable(
					id='protein_cut_datatable',
					data=df.to_dict('records'),
					columns=[{'name': i, 'id': i} for i in df.columns],
					style_header={
						'backgroundColor': 'rgba(148,163,188)',
						'fontWeight': 'bold',
						'color': 'black'
					},
					style_cell={
						'backgroundColor': 'rgba(34,34,34)',
						'color': 'white',
						'minWidth': '150px', 'width': '150px', 'maxWidth': '150px',
						'overflow': 'hidden',
						'textOverflow': 'ellipsis',
					},
					page_size=20,
					editable=True,
					row_deletable=True,
					fixed_columns={'headers': True, 'data': 1},
					style_table={'minWidth': '100%', 'maxWidth': '100%'},
				),
			)

			return datatable_children, no_update, df.to_json(date_format='iso', orient='split'), children, no_update, no_update
		else :
			return datatable_children, no_update, no_update, children, no_update, no_update

	elif triggered == 'protein_cut_download_button' :

		if protein_cut_stored == None or protein_cut_stored == '{"columns":[],"index":[],"data":[]}' :
			df = pd.DataFrame(columns=['protein', 'start', 'length'])
		else :
			df = pd.read_json(protein_cut_stored, orient='split')
			df = df[['protein', 'start', 'length']]

		return no_update, no_update, no_update, no_update, send_data_frame(df.to_excel, filename='cut_regions.xlsx', index=False), no_update

	elif datatable_click != None and datatable_click > 0 :
		return no_update, no_update, no_update, no_update, no_update, need_save

	else :
		return no_update, no_update, no_update, no_update, no_update, no_update

@app.callback(
	Output('NR_figure_download', 'data'),
	Output('NR_legend_download', 'data'),
	Input('normal_run_result_download', 'n_clicks'),
	State('Figure_code', 'data'),
	State('Legend_code', 'data'),
	prevent_initial_call=True
	)
def update_output(n_clicks, Figure_code, Legend_code) :

	return dcc.send_file(os.path.join('Figures', Figure_code+'.png')), dcc.send_file(os.path.join('Figures', Legend_code+'.png'))

@app.callback(
	Output('QR_figure_download', 'data'),
	Output('QR_legend_download', 'data'),
	Input('quick_run_result_download', 'n_clicks'),
	State('QR_Figure_code', 'data'),
	State('QR_Legend_code', 'data'),
	prevent_initial_call=True
	)
def update_output(n_clicks, QR_Figure_code, QR_Legend_code) :

	return dcc.send_file(os.path.join('Figures', QR_Figure_code+'.png')), dcc.send_file(os.path.join('Figures', QR_Legend_code+'.png'))

@app.callback(
	Output('obsolete_code_download', 'data'),
	Input('obsolete_code_download_btn', 'n_clicks'),
	State('obsolete_list', 'data'),
	State('unrecognized_list', 'data'),
	prevent_initial_call=True
	)
def update_output(n_clicks, obsolete_list, unrecognized_list) :


	writer = pd.ExcelWriter('Uniprot_code_errors.xlsx', engine='xlsxwriter')

	df = pd.read_json(obsolete_list, orient='split')
	# print(df)
	df.to_excel(writer, sheet_name='obsolete', index=False)

	df = pd.read_json(unrecognized_list, orient='split')
	# print(df)
	df.to_excel(writer, sheet_name='unrecognized', index=False)

	writer.save()

	# return send_data_frame(to_xlsx, filename='pandas_multiple.xlsx')
	return dcc.send_file('Uniprot_code_errors.xlsx')

@app.callback(
	Output('feature_extract_placeholder', 'children'),
	Output('feature_extract', 'data'),
	Input('feature_extract_launch_button', 'n_clicks'),
	State('feature_extract_name', 'value'),
	State('feature_extract_nter_ext', 'value'),
	State('feature_extract_cter_ext', 'value'),
	State('2_to_3_features', 'data'),
	State('2_to_3_protein_seq', 'data'),
	prevent_initial_call=True
	)
def update_output(n_clicks, feature_name, nter_ext, cter_ext, feature_data, protein_seq_data):
	ctx = dash.callback_context
	triggered = ctx.triggered[0]['prop_id'].split('.')[0]
	print("TRIGERED: "+triggered)

	feature_name = feature_name.lower()
	# feature_name = feature_name.replace(' ', '_')

	if nter_ext.strip() != '' and nter_ext != None :
		nter_ext = int(nter_ext)

	if cter_ext.strip() != '' and cter_ext != None :
		cter_ext = int(cter_ext)

	if feature_name.strip() != '' and feature_name != None :
		feature_name = feature_name.strip()

		feature_df = pd.read_json(feature_data, orient='split')
		feature_df = feature_df[['protein', 'feature_type', 'feature', 'start', 'length', 'intensity', 'variant']]

		protein_seq = protein_seq_data

		# print(protein_seq)

		seq_dict = {}

		protein_name = ''
		for line in protein_seq.split('\n') :
			# print(line)
			if line != '' :
				if line[0] == '>' :
					protein_name = line[1:]
					seq_dict[protein_name] = ''
				else :
					seq_dict[protein_name] += line

		output_str = ''
		last_prot = ''
		i = 1

		# print(list(feature_df))

		for index, row in feature_df.iterrows():
			if row['feature'].lower() == feature_name :
				protein = row['protein']
				start = int(row['start'])
				length = int(row['length'])
				end = start + length - 1
				start_plus_ext = start-nter_ext-1
				end_plus_ext = end+cter_ext

				if start_plus_ext < 0 :
					start_plus_ext = 0
				if end_plus_ext < 0 :
					start_plus_ext = 0
					end_plus_ext = 0
				if start_plus_ext > len(seq_dict[protein]) :
					start_plus_ext = len(seq_dict[protein])
					end_plus_ext = len(seq_dict[protein])
				if end_plus_ext > len(seq_dict[protein]) :
					end_plus_ext = len(seq_dict[protein])

				seq = seq_dict[protein][start_plus_ext:end_plus_ext]
				# print(protein, seq)
				if protein == last_prot :
					i += 1
				else :
					i = 1
					last_prot = protein
				output_str += '>' + protein + '-' + str(i) + '\n'
				output_str += seq + '\n'

		# print(output_str)
		return no_update, output_str

	else :
		return no_update, no_update


@app.callback(
	Output('feature_regex_placeholder', 'children'),
	Output('feature_regex_data', 'data'),
	Input('feature_regex_launch_button', 'n_clicks'),
	State('feature_regex', 'value'),
	State('feature_regex_name', 'value'),
	State('2_to_3_protein_seq', 'data'),
	prevent_initial_call=True
	)
def update_output(n_clicks, feature_regex, feature_regex_name, protein_seq):
	ctx = dash.callback_context
	triggered = ctx.triggered[0]['prop_id'].split('.')[0]
	print("TRIGERED: "+triggered)

	if feature_regex != None and feature_regex.strip() != '' and feature_regex_name != None and feature_regex_name.strip() != '' :
		# print('regex', feature_regex, feature_regex_name)

		# print(protein_seq)

		seq_dict = {}
		protein_name = ''
		for line in protein_seq.split('\n') :
			if line != '' :
				if line[0] == '>' :
					protein_name = line[1:]
					seq_dict[protein_name] = ''
				else :
					seq_dict[protein_name] += line

		# print(seq_dict)

		df = pd.DataFrame(columns=['ex_type', 'protein', 'feature_type', 'feature', 'start', 'length'])
		df_id = 0

		for protein, sequence in seq_dict.items() :
			for m in re.finditer(feature_regex, sequence):
				# print(protein, m.start(), m.group())
				df.loc[df_id] = ['+', protein, 'DOMAIN', feature_regex_name, m.start()+1, m.span()[1]-m.start()] # len = end-start here
				df_id += 1

		# print(df)
		regex_data = df.to_json(date_format='iso', orient='split')

		return no_update, regex_data

	return no_update, no_update

@app.callback(
	Output('feature_regex_download', 'data'),
	Input('feature_regex_download_button', 'n_clicks'),
	State('feature_regex_name', 'value'),
	State('feature_regex_data', 'data'),
	prevent_initial_call=True
	)
def update_output(n_clicks, feature_regex_name, feature_regex_data) :

	if feature_regex_name != None and feature_regex_data != None and feature_regex_name.strip() != '' :
		feature_regex_name = feature_regex_name.strip()
		df = pd.read_json(feature_regex_data, orient='split')

		return send_data_frame(df.to_excel, filename=feature_regex_name+'.xlsx', index=False)

	else :
		return no_update

@app.callback(
	Output('feature_extract_download', 'data'),
	Input('feature_extract_download_button', 'n_clicks'),
	State('feature_extract_name', 'value'),
	State('feature_extract', 'data'),
	prevent_initial_call=True
	)
def update_output(n_clicks, feature_name, feature_extract_data) :

	if feature_name != None and feature_extract_data != None and feature_name.strip() != '' :
		feature_name = feature_name.strip()

		return dict(content=feature_extract_data, filename=feature_name+'.fa')

	else :
		return no_update

@app.callback(
	Output('protein_length_download', 'data'),
	Input('protein_length_download_btn', 'n_clicks'),
	State('2_to_4_protein_length', 'data'),
	State('2_to_3_features', 'data'),
	State('2_to_3_pdb_list', 'data'),
	State('Step_4_feature_occurrence', 'data'),
	State('2_to_3_protein_seq', 'data'),
	prevent_initial_call=True
	)
def update_output(n_clicks, protein_length_data, features_data, pdb_list_data, feature_occurrence_data, protein_seq) :

	if protein_length_data == None :
		return no_update
	else :
		writer = pd.ExcelWriter('extracted_data.xlsx', engine='xlsxwriter')

		df = pd.read_json(features_data, orient='split')
		df = df[['code', 'protein', 'feature_type', 'feature', 'start', 'length', 'intensity', 'variant']]
		df.to_excel(writer, sheet_name='feature_list', index=False)

		df = pd.read_json(feature_occurrence_data, orient='split')
		df = df[['feature_type', 'feature', 'occurrence']]
		df = df.sort_values(by='occurrence', ascending=False)
		df.to_excel(writer, sheet_name='feature_occurrence', index=False)

		df = pd.read_json(pdb_list_data, orient='split')
		df = df[['code', 'protein', 'PDB', 'Method', 'Resolution', 'Chains', 'Start', 'End']]
		df.to_excel(writer, sheet_name='pdb_list', index=False)

		df = pd.read_json(protein_length_data, orient='split')
		df = df[['code', 'protein', 'total_length']]
		df.to_excel(writer, sheet_name='protein_length', index=False)

		seq_dict = {}

		protein_name = ''
		for line in protein_seq.split('\n') :
			if line != '' :
				if line[0] == '>' :
					protein_name = line[1:]
					seq_dict[protein_name] = ''
				else :
					seq_dict[protein_name] += line

		df = pd.DataFrame(columns=['sequence'])
		df_i = 0
		for protein, sequence in seq_dict.items() :
			df.loc[df_i] = '>'+protein
			df_i += 1
			df.loc[df_i] = sequence
			df_i += 1

		df.to_excel(writer, sheet_name='protein_sequence', index=False, header=False)


		writer.save()

	# return send_data_frame(to_xlsx, filename='pandas_multiple.xlsx')
	return dcc.send_file('extracted_data.xlsx')

@app.callback(
	Output('pdb_threshold_list', 'value'),
	Output('nb_steps_pdb_threshold', 'value'),
	Output('2_to_3_pdb_threshold_list', 'data'),
	Output('Step_4_pdb_coverage_list', 'data'),
	Output('Step_4_pdb_coverage_steps', 'data'),

	Input('nb_steps_pdb_threshold', 'value'),
	Input('pdb_threshold_default_btn', 'n_clicks'),
	Input('pdb_threshold_list', 'value'),

	State('Step_4_pdb_coverage_list', 'data'),
	State('Step_4_pdb_coverage_steps', 'data'),
	)
def update_output(nb_steps, default, threshold_list, threshold_list_saved, nb_steps_saved):
	ctx = dash.callback_context
	triggered = ctx.triggered[0]['prop_id'].split('.')[0]
	print("TRIGERED_coverage: "+triggered)
	print(threshold_list_saved, nb_steps_saved)

	marks_list = [1,2,3,5,10,15,20,25,30,50]

	if triggered == 'pdb_threshold_default_btn' :
		return [marks_list[i] for i in range(5)], 5, [marks_list[i] for i in range(5)], [marks_list[i] for i in range(5)], 5

	if triggered == 'nb_steps_pdb_threshold' :
		return [marks_list[i] for i in range(nb_steps)], no_update, [marks_list[i] for i in range(nb_steps)], [marks_list[i] for i in range(nb_steps)], nb_steps

	if triggered == 'pdb_threshold_list' :
		return no_update, no_update, threshold_list, threshold_list, no_update

	else :
		if nb_steps_saved == None :
			nb_steps_saved = nb_steps
		if threshold_list_saved == None :
			threshold_list_saved = threshold_list

		return threshold_list_saved, nb_steps_saved, [marks_list[i] for i in range(5)], threshold_list_saved, nb_steps_saved

@app.callback(
	Output('entry_nb', 'children'),
	Output('unique_entry_nb', 'children'),
	Input('Step_1', 'data'),
	)
def update_output(data):
	if data == None :
		return no_update, no_update
	else :
		df = pd.read_json(data, orient='split')
		if not df.empty :
			full_list = list(df['code'])
			full_list = [protein_name.strip() for protein_name in full_list if protein_name.strip() != '']
			full_list_nb = 'Proteins: '+str(len(full_list))
			unique_list_nb = 'Unique proteins: '+str(len(list(set(full_list))))

			return full_list_nb, unique_list_nb

		else :
			full_list_nb = 'Proteins: '+str(0)
			unique_list_nb = 'Unique proteins: '+str(0)

			return full_list_nb, unique_list_nb

# Step 2 and 3 enter
@app.callback(
	Output('step_2_start_state', 'children'),
	Output('step_3_start_state', 'children'),

	Input('step_2_launch_button', 'n_clicks'),
	Input('step_3_launch_button', 'n_clicks'),

	prevent_initial_call=True
	)
def download_uniprot_files(button_press, button_press_2):

	ctx = dash.callback_context
	triggered = ctx.triggered[0]['prop_id'].split('.')[0]
	print("TRIGERED: "+triggered)

	if triggered == 'step_2_launch_button' :

		return html.Img(id='step_2_end_state', src=b64_image(os.path.join('GUI', 'fail'+'.png')), style={'width':'4vh', 'height':'4vh'}), no_update

	if triggered == 'step_3_launch_button' :

		return no_update, html.Img(id='step_3_end_state', src=b64_image(os.path.join('GUI', 'fail'+'.png')), style={'width':'4vh', 'height':'4vh'})

# Step 2 and 3 main
@app.callback(
	Output('2_to_3_features', 'data'),
	Output('2_to_4_protein_length', 'data'),
	Output('Step_4_feature_occurrence', 'data'),
	Output('2_to_3_pdb_list', 'data'),
	Output('2_to_3_pdb_coverage', 'data'),
	Output('Step_3_db_PDB', 'data'),
	Output('2_to_3_protein_seq', 'data'),
	Output('step_2_placeholder', 'children'),
	Output('step_2_end_state', 'src'),

	Output('final_store', 'data'),
	Output('step_3_placeholder', 'children'),
	Output('step_3_end_state', 'src'),

	Output('direct_download', 'data'),
	Output('step_2_extraction_collapse', 'is_open'),

	Input('step_2_start_state', 'children'),

	Input('step_3_start_state', 'children'),

	State('Step_1', 'data'),
	State('Step_2_exceptions', 'data'),
	State('2_to_3_pdb_threshold_list', 'data'),

	State('Step_3_db_user', 'data'),
	State('Step_3_db_PDB', 'data'),
	State('2_to_3_features', 'data'),

	State('big_data_toggle', 'value'),

	State('step_2_extraction_collapse', 'is_open'),

	prevent_initial_call=True
	)
def download_uniprot_files(
	button_press,
	start,

	protein_list_data,
	exception_data,
	pdb_threshold_list,

	intensity_data,
	db_PDB,
	features,

	big_data,
	collapse_state
	):
	ctx = dash.callback_context
	triggered = ctx.triggered[0]['prop_id'].split('.')[0]
	print("TRIGERED: "+triggered)

	if triggered == 'step_2_start_state' :

		start_time = time.time()

		if exception_data != None :
			modification_df = pd.read_json(exception_data, orient='split')
			modification_df = modification_df[['ex_type', 'protein', 'feature_type', 'feature', 'start', 'length']]
		else :
			modification_df = pd.DataFrame(columns=['ex_type', 'protein', 'feature_type', 'feature', 'start', 'length'])

		protein_list_df = pd.read_json(protein_list_data, orient='split')


		protein_length_df, protein_occurrence_df, pdb_list_df, coverage_df, db_coverage_df, features_df, protein_seq_str = featureExtraction(protein_list_df, modification_df=modification_df, pdb_threshold_list=pdb_threshold_list)

		protein_length_json = protein_length_df.to_json(date_format='iso', orient='split')
		protein_occurrence_json = protein_occurrence_df.to_json(date_format='iso', orient='split')
		pdb_list_json = pdb_list_df.to_json(date_format='iso', orient='split')
		coverage_json = coverage_df.to_json(date_format='iso', orient='split')
		db_coverage_json = db_coverage_df.to_json(date_format='iso', orient='split')
		features_json = features_df.to_json(date_format='iso', orient='split')

		send_big_data = no_update
		if big_data == ['big_data'] :

			# If user has a big list, which extracted data can not be saved localy on this browser, he can select this option to download directly extracted data.
			# When this option is chosen, it will not be possilbe to create a map.
			# A theoretical limit of the browser must be around 5k proteins. 

			writer = pd.ExcelWriter('extracted_data.xlsx', engine='xlsxwriter')

			features_df = features_df[['code', 'protein', 'feature_type', 'feature', 'start', 'length', 'intensity', 'variant']]
			features_df.to_excel(writer, sheet_name='feature_list', index=False)

			protein_occurrence_df = protein_occurrence_df[['feature_type', 'feature', 'occurrence']]
			protein_occurrence_df = protein_occurrence_df.sort_values(by='occurrence', ascending=False)
			protein_occurrence_df.to_excel(writer, sheet_name='feature_occurrence', index=False)

			pdb_list_df = pdb_list_df[['code', 'protein', 'PDB', 'Method', 'Resolution', 'Chains', 'Start', 'End']]
			pdb_list_df.to_excel(writer, sheet_name='pdb_list', index=False)

			protein_length_df = protein_length_df[['code', 'protein', 'total_length']]
			protein_length_df.to_excel(writer, sheet_name='protein_length', index=False)

			seq_dict = {}

			protein_name = ''
			for line in protein_seq_str.split('\n') :
				if line != '' :
					if line[0] == '>' :
						protein_name = line[1:]
						seq_dict[protein_name] = ''
					else :
						seq_dict[protein_name] += line

			sequence_dict_list = []

			for protein, sequence in seq_dict.items() :
				sequence_dict_list.append({'sequence':'>'+protein})
				sequence_dict_list.append({'sequence':sequence})

			sequence_df = pd.DataFrame.from_dict(sequence_dict_list)

			sequence_df.to_excel(writer, sheet_name='protein_sequence', index=False, header=False)

			writer.save()

			send_big_data = dcc.send_file('extracted_data.xlsx')

		end_time = time.time()
		print('Step 2 execution time:', str(end_time - start_time))

		return features_json, protein_length_json, protein_occurrence_json, pdb_list_json, coverage_json, db_coverage_json, protein_seq_str, no_update, b64_image(os.path.join('GUI', 'success'+'.png')), no_update, no_update, no_update, send_big_data, True

	if triggered == 'step_3_start_state' :

		start_time = time.time()

		db_user_df = pd.read_json(intensity_data, orient='split')
		# print(db_user_df)
		feature_list = list(set(list(db_user_df['feature'])))
		feature_list.sort()
		clean_feature_list = []
		for feature in feature_list :
			if feature.strip() != '' :
				clean_feature_list.append(feature.strip())
		# print(clean_feature_list)

		case_list = []
		for column in db_user_df.columns :
			if column not in ['protein', 'feature', 'start'] :
				case_list.append(column)
		# print(case_list)
		if case_list == [] :
			case_list = ['noCase']

		if clean_feature_list != [] :

			features_df = pd.read_json(features, orient='split')
			df = findData(case_list, clean_feature_list, db_user_df, db_PDB, features_df)

			output_PDB = findData(['PDB'], ['PDB'], db_user_df, db_PDB, features_df)

			children = []

			pdb_df = output_PDB

			without_PDB_dict_list = []
			header_list = ['protein', 'feature_type', 'feature', 'start', 'length']+case_list

			for index, row in df.iterrows():
				if row['feature'] != 'PDB' :
					without_PDB_dict_list.append({x[0] : x[1] for x in zip(header_list, list(row))})
					# print({x[0] : x[1] for x in zip(header_list, list(row))})
			for index, row in pdb_df.iterrows():
				if row['feature'] == 'PDB' :
					new_row = [row[elem] for elem in ['protein', 'feature_type', 'feature', 'start', 'length']+['PDB' for i in range(len(case_list))]]
					without_PDB_dict_list.append({x[0] : x[1] for x in zip(header_list, new_row)})

			without_PDB_df = pd.DataFrame.from_dict(without_PDB_dict_list)


			children = without_PDB_df.to_json(date_format='iso', orient='split')

			# without_PDB_df.to_csv(elem[0]+'.csv', index=False)

			end_time = time.time()
			print('Step 3 execution time:', str(end_time - start_time))

			return no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, children, no_update, b64_image(os.path.join('GUI', 'success'+'.png')), no_update, False

		else :

			return no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, None, no_update, b64_image(os.path.join('GUI', 'success'+'.png')), no_update, False

# Step 4 case dropdown update
@app.callback(
	Output('case_dropdown', 'options'),
	Output('case_dropdown', 'value'),
	Input('final_store', 'data'),
	)
def download_uniprot_files(final_store):
	if final_store != None :
		df = pd.read_json(final_store, orient='split')
		columns_list = list(df)
		# print(columns_list)

		for elem in ['protein', 'feature_type', 'feature', 'start', 'length'] :
			columns_list.remove(elem)
		columns_list.sort()
		# print(columns_list)

		if columns_list == ['noCase'] :
			options = None
			value = None
		else :
			options = [{"label": x, "value": x} for x in columns_list]
			value = columns_list[0]

		return options, value

	else :
		return no_update, no_update

# Step 4 enter
@app.callback(
	Output('step_4_start_state', 'children'),
	Input('step_4_launch_button', 'n_clicks'),
	prevent_initial_call=True
	)
def download_uniprot_files(button_press):

	return html.Img(id='step_4_end_state', src=b64_image(os.path.join('GUI', 'fail'+'.png')), style={'width':'4vh', 'height':'4vh'})

# Step 4 main
@app.callback(
	Output('result_img', 'children'),
	Output('legend_img', 'children'),
	Output('step_4_placeholder', 'children'),
	Output('step_4_end_state', 'src'),

	Output('Sorted_protein_list', 'data'),

	Output('Figure_code', 'data'),
	Output('Legend_code', 'data'),

	Output('map_download', 'data'),
	Output('legend_download', 'data'),

	Input('step_4_start_state', 'children'),

	State('final_store', 'data'),

	State('case_dropdown', 'value'),
	State('length_factor', 'value'),
	State('height', 'value'),
	State('text_size', 'value'),
	State('biased_region_text_size', 'value'),
	State('draw_region_dropdown', 'value'),

	State('coverage_toggles', 'value'),
	State('2nd_struct_toggles', 'value'),
	State('disorder_toggles', 'value'),
	State('mod_res_toggles', 'value'),

	State('protein_length_toggles', 'value'),
	State('uniform_toggles', 'value'),

	State('pensize', 'value'),
	State('sorting_dropdown', 'value'),
	State('focus_input', 'value'),
	State('threshold', 'value'),
	State('FT_order_list', 'value'),

	State('2_to_3_features', 'data'),
	State('Step_3_feature_shape', 'data'),
	State('Step_3_protein_cut', 'data'),
	State('2_to_4_protein_length', 'data'),
	State('Step_1', 'data'),
	State('Step_4_pdb_coverage_list', 'data'),

	State('Figure_code', 'data'),
	State('Legend_code', 'data'),
	)
def blabla(
	button_click,

	final_store,

	case_dropdown,
	length_factor,
	height,
	text_size,
	biased_region_text_size,
	draw_region_dropdown,
	coverage_toggles,
	second_struct_toggles,
	disorder_toggles,
	mod_res_toggles,
	protein_length_toggles,
	uniform_toggles,

	pensize,
	sorting_dropdown,
	focus_input,
	threshold,
	FT_order_list,

	features_list_data,
	feature_shape,
	protein_cut_data,
	protein_length_data,
	protein_list_data,
	pdb_coverage_list,

	old_figure_code,
	old_legend_code
	) :

	ctx = dash.callback_context
	triggered = ctx.triggered[0]['prop_id'].split('.')[0]
	print("TRIGERED: "+triggered)

	if triggered == 'step_4_start_state' :

		start_time = time.time()

		parameter_dict = prepareMapParameters(
			length_factor, 
			height, 
			text_size,
			sorting_dropdown, 
			focus_input, 
			threshold, 
			biased_region_text_size, 
			uniform_toggles, 
			draw_region_dropdown, 
			coverage_toggles, 
			second_struct_toggles, 
			disorder_toggles, 
			mod_res_toggles, 
			protein_length_toggles, 
			pensize, 
			FT_order_list)

		if case_dropdown == None or case_dropdown == 'None' :
			case_dropdown = 'noCase'

		print('----------------Step_4----------------')
		protein_list_df = pd.read_json(protein_list_data, orient='split')

		if final_store != None and case_dropdown != 'noCase' :
			data_df = pd.read_json(final_store, orient='split')
		else :
			data_df = pd.read_json(features_list_data, orient='split')

		# Make a DataFrame out of the stored protein length data (json)
		protein_length_df = pd.read_json(protein_length_data, orient='split')

		# If data is store, make a DataFrame out of it
		if protein_cut_data != None :
			protein_cut_df = pd.read_json(protein_cut_data, orient='split')
		# If no stored data is found, create empty template DataFrame
		else :
			protein_cut_df = pd.DataFrame(columns=['protein', 'start', 'length'])

		# Make DataFrame from stored feature shape an color data
		shape_df = pd.read_json(feature_shape, orient='split')

		im_children, legend_children, uuid_im, uuid_legend, sorted_protein_list = drawFigure(
			protein_list_df,
			data_df,
			shape_df,
			protein_cut_df,
			protein_length_df,
			pdb_coverage_list,
			case_dropdown,
			**parameter_dict
		)

		end_time = time.time()
		print('Step 4 execution time:', str(end_time - start_time))

		# Removing old figures if they exist
		if old_figure_code != None and os.path.exists(os.path.join('Figures', old_figure_code+'.png')) :
			os.remove(os.path.join('Figures', old_figure_code+'.png'))

		if old_legend_code != None and os.path.exists(os.path.join('Figures', old_legend_code+'.png')) :
			os.remove(os.path.join('Figures', old_legend_code+'.png'))

		return im_children, legend_children, no_update, b64_image(os.path.join('GUI', 'success'+'.png')), sorted_protein_list, uuid_im, uuid_legend, no_update, no_update #dcc.send_file(os.path.join('Figures', uuid_im+'.png')), dcc.send_file(os.path.join('Figures', uuid_legend+'.png'))

	else :

		im_children = html.Div()
		legend_children = html.Div()

		if old_figure_code != None and os.path.exists(os.path.join('Figures', old_figure_code+'.png')) :
			im_children = html.Img(src=b64_image(os.path.join('Figures', old_figure_code+'.png')), style={'width' : '100%', 'display' : 'flex', 'align-items' : 'center', 'justify-content' : 'center'})

		if old_legend_code != None and os.path.exists(os.path.join('Figures', old_legend_code+'.png')) :
			legend_children = html.Img(src=b64_image(os.path.join('Figures', old_legend_code+'.png')), style={'width' : '100%', 'display' : 'flex', 'align-items' : 'center', 'justify-content' : 'center'})

		return im_children, legend_children, no_update, no_update, no_update, no_update, no_update, no_update, no_update

@app.callback(
	Output('auto_toggle_1', 'value'),
	Output('auto_toggle_2', 'value'),
	Input('auto_toggle_1', 'value'),
	Input('auto_toggle_2', 'value'),
	prevent_initial_call=True
	)
def main(x1, x2) :
	ctx = dash.callback_context
	triggered = ctx.triggered[0]['prop_id'].split('.')[0]
	print("TRIGERED: "+triggered)
	choice_list = [[],['most_represented']]

	if triggered == 'auto_toggle_1' :
		return no_update, choice_list[choice_list.index(x1)-1]
	elif triggered == 'auto_toggle_2' :
		return choice_list[choice_list.index(x2)-1], no_update

@app.callback(
	Output('auto_select_threshold_1', 'value'),
	Output('auto_select_threshold_2', 'value'),
	Input('auto_select_threshold_1', 'value'),
	Input('auto_select_threshold_2', 'value'),
	prevent_initial_call=True
	)
def main(x1, x2) :
	ctx = dash.callback_context
	triggered = ctx.triggered[0]['prop_id'].split('.')[0]
	print("TRIGERED: "+triggered)

	if triggered == 'auto_select_threshold_1' :
		return no_update, x1
	elif triggered == 'auto_select_threshold_2' :
		return x2, no_update


# Ultra Quick Run
@app.callback(
	Output('QR_Figure_code', 'data'),
	Output('QR_Legend_code', 'data'),

	Output('UQR_result_img', 'children'),
	Output('UQR_legend_img', 'children'),
	Output('quick_run_collapse', 'is_open'),
	Output('step_UQR_placeholder', 'children'),

	Input('UQR_button', 'n_clicks'),
	Input('hide_quick_run_btn', 'n_clicks'),

	State('dl_latest_toggle', 'value'),

	State('Step_1', 'data'),
	State('unique_entry_nb', 'children'),

	State('length_factor', 'value'),
	State('height', 'value'),
	State('text_size', 'value'),
	State('biased_region_text_size', 'value'),
	State('draw_region_dropdown', 'value'),

	State('coverage_toggles', 'value'),
	State('2nd_struct_toggles', 'value'),
	State('disorder_toggles', 'value'),
	State('mod_res_toggles', 'value'),

	State('protein_length_toggles', 'value'),
	State('uniform_toggles', 'value'),

	State('pensize', 'value'),
	State('sorting_dropdown', 'value'),
	State('FT_order_list', 'value'),

	State('quick_run_collapse', 'is_open'),

	prevent_initial_call=True
	)
def download_uniprot_files(
	button_press,
	hide_quick_run_btn,

	latest,

	protein_list_data,
	protein_nb,

	length_factor,
	height,
	text_size,
	biased_region_text_size,
	draw_region_dropdown,
	coverage_toggles,
	second_struct_toggles,
	disorder_toggles,
	mod_res_toggles,
	protein_length_toggles,
	uniform_toggles,

	pensize,
	sorting_dropdown,
	FT_order_list,

	quick_run_collapse
	):

	ctx = dash.callback_context
	triggered = ctx.triggered[0]['prop_id'].split('.')[0]
	print("TRIGERED_load_param: "+triggered)

	if triggered == 'UQR_button' :	

		print('UQR initiated...')
		complete_protein_list_df = pd.read_json(protein_list_data, orient='split')
		missing_protein_names = False
		if 'protein' not in list(complete_protein_list_df) :
			print('No protein names detected')
			missing_protein_names = True
			protein_name_list = getMissingProteinNames(complete_protein_list_df)
			complete_protein_list_df['protein'] = protein_name_list
		print('Entering: Protein data gathering...')
		proteinDataGathering(complete_protein_list_df, dl_latest=latest)
		print('Quiting: Protein data gathering')
		print('Entering: Data extraction...')
		protein_length_df, protein_occurrence_df, pdb_list_df, coverage_df, db_coverage_df, features_df, protein_seq_str = featureExtraction(complete_protein_list_df)
		print('Quiting: Data extraction...')
		print('Entering: Auto feature choice...')
		shape_df = autoFeatChoice(protein_occurrence_df)
		print(shape_df)
		print('Quiting: Auto feature choice...')
		print('Entering: Map creation...')

		parameter_dict = prepareMapParameters(
			length_factor, 
			height, 
			text_size,
			sorting_dropdown, 
			None, 
			None, 
			biased_region_text_size, 
			uniform_toggles, 
			draw_region_dropdown, 
			coverage_toggles, 
			second_struct_toggles, 
			disorder_toggles, 
			mod_res_toggles, 
			protein_length_toggles, 
			pensize, 
			FT_order_list)

		im_children, legend_children, uuid_im, uuid_legend, sorted_protein_list = drawFigure(
			complete_protein_list_df,
			features_df,
			shape_df,
			pd.DataFrame(columns=['protein', 'start', 'length']),
			protein_length_df,
			[1,2,3,5,10],
			'noCase', #case_dropdown
			**parameter_dict
		)
		print('Quiting: Map creation...')

		return uuid_im, uuid_legend, im_children, legend_children, True, no_update

	else :
		return no_update, no_update, no_update, no_update, False, no_update



if __name__ == '__main__':
	app.run_server(debug=True)



