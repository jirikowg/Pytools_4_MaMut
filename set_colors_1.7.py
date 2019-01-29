'''
THis program allows you to assign a customized color code to cell lineage data in MAMuT-XML format.
It requires that information on cell cenerations was previously aquired and written to the MaMuT-File
using "lineage_from_mamut"
'''


from lxml import etree as ET
from io import StringIO, BytesIO
import os, sys
import glob
import math
from colour import Color


# Create dictionary with basic colors mapped to initials
col_dict = {}
col_dict['r'] = -65536
col_dict['g'] = -16711936
col_dict['b'] = -16776961
col_dict['y'] = -256
col_dict['o'] = -23296
col_dict['c'] = -16711681
col_dict['m'] = -65281
col_dict['bl'] = -16777216
col_dict['w'] = -1
col_dict['gr'] = -8421505

rgb_dict = {}
rgb_dict['r'] = [1,0,0]
rgb_dict['g'] = [0,1,0]
rgb_dict['b'] = [0,0,1]
rgb_dict['y'] = [1,1,0]
rgb_dict['o'] = [1,0.6,0]
rgb_dict['c'] = [0,1,1]
rgb_dict['m'] = [1,0,1]
rgb_dict['bl'] = [0,0,0]
rgb_dict['w'] = [1,1,1]
rgb_dict['gr'] = [0.5,0.5,0.5]


# Create a set of base colors:
bas_cols = ['r','g','b','y','o','c','m','bl','w','gr']

print(col_dict)
	
# calculates color value in WIN-format from either base colors or rgb-values
def col_trans (cols_slin, col_dict):
	for i in range(len(cols_slin)):
		if cols_slin[i][1] in bas_cols:
			cd = cols_slin[i][1]
			cols_slin[i][1] = col_dict[cd]
		else:
			s_col = cols_slin[i][1].split(',')
			B = int(s_col[2])
			G = int(s_col[1])
			R = int(s_col[0])
			# RGB to Windows(decimal)
			WIN = (-16777216 + (R*65536) + (G*256) + B)
			print('WIN:',WIN)

			cols_slin[i][1] = WIN
	return cols_slin


# Function that calculates base colors or rgb entries (0-255) to rgb-format used by the "Colour"-module (r/255,g/255,b/255)

def col_2_rgb(cols,cols_n): 						# The argument "cols" is a list of colors in either base color format or (r,g,b) format
	for i in range (len(cols)):
		if cols[i] in bas_cols:
			cd = cols[i]
			cols_n.append(rgb_dict[cd])
		else:
			cd = cols[i].split(',')
			B = int(cd[2])/255
			G = int(cd[1])/255
			R = int(cd[0])/255
			cd = [R,G,B]
			cols_n.append(cd)
	return cols,cols_n



# Calculate Code color by generation

# This function takes user-specified root and leaf color, calculates a range of colors from them with the length of the total generatinons in the track.
# Converts the colors to WIN-format writes them to a list called l_win[]. 
col_range_win = []
def create_color_range (cols,gens,col_range_win):			
	cols_n = []
	col_2_rgb(cols,cols_n)				# funcion
	rot = Color(rgb=(cols_n[0][0],cols_n[0][1],cols_n[0][2]))
	lef = Color(rgb=(cols_n[1][0],cols_n[1][1],cols_n[1][2]))
	ge = int(gens)
	col_range = list(rot.range_to(lef,ge))		
	print('col_range:',col_range)

	# extract r, g and b values from PI-color-format
	for i in range(len(col_range)):
		print('(col_range[i]',col_range[i])
		l_red = int(float(col_range[i].red)*255)
		print ('l_red:',l_red)
		l_green = int(float(col_range[i].green)*255)
		print ('l_green:',l_green)
		l_blue = int(float(col_range[i].blue)*255)
		print ('l_blue:',l_blue)
		col_win_i = (-16777216 + (l_red*65536) + (l_green*256) + l_blue)
		col_range_win.append(col_win_i)
		print('col_range_win:',col_range_win)
	return col_range_win
#!!! There can be eerors when calculating a range including cyan AND blue (both calculated as blue)



#Version of "create_color_range" modified to alculate color range within sublineages. 
def create_color_range_sublin (gens,col_range_win,cols_slin_gen,cols_slin_gen_ran):			
	for i in range (len(cols_slin_gen)):
		cols = []
		cols_n = []
		slin_gen_ran = []
		col_range_win = []
		slin_gen_ran.append(cols_slin_gen[i][0])
		cols.append(cols_slin_gen[i][1])
		cols.append(cols_slin_gen[i][2])
		col_2_rgb(cols,cols_n)	
		rot = Color(rgb=(cols_n[0][0],cols_n[0][1],cols_n[0][2]))
		lef = Color(rgb=(cols_n[1][0],cols_n[1][1],cols_n[1][2]))
		ge = int(gens)+2										#!!! Added 2 generations
		col_range = list(rot.range_to(lef,ge))	
	
		# extract r, g and b values from PI-color-format and calculate WIN-format
		for i in range(len(col_range)):
			l_red = int(float(col_range[i].red)*255)
			l_green = int(float(col_range[i].green)*255)
			l_blue = int(float(col_range[i].blue)*255)
			col_win_i = (-16777216 + (l_red*65536) + (l_green*256) + l_blue)
			col_range_win.append(col_win_i)
		#print('col_range_win:',col_range_win)		
		slin_gen_ran.append(col_range_win)
		cols_slin_gen_ran.append(slin_gen_ran)
		#print('cols_slin_gen_ran',cols_slin_gen_ran)
	return cols_slin_gen_ran
	#print('cols_slin_gen_ran:',cols_slin_gen_ran)



# Import xml file:

get_file = input('enter name of mamut-xml (file without suffix) for which you whant to set the colors. >')
get_file_suff = get_file + '.xml'
export_mamut = ET.parse(get_file_suff)

root = export_mamut.getroot()



#Funktion that enables user to apply color of spots and edges in tree to direct targes spots e.g. after interpolation with MaMut. 
def adapt_colors ():
	print('adapt colors ...')
	# get number of spots 
	n_spots = 0
	for Spot in root.iter('Spot'):
		n_spots =+ 1
	print ('n_spots:',n_spots)

	# get Edge List:
	edge_list = []
	for Edge in root.iter('Edge'):
		edge = []
		edge.append(int(Edge.attrib['SPOT_SOURCE_ID']))
		edge.append(int(Edge.attrib['SPOT_TARGET_ID']))
		edge_list.append(edge)
	# get mit_spots
	source_list=[]
	for Edge in root.iter('Edge'):
		SPOT_SOURCE_ID = Edge.attrib['SPOT_SOURCE_ID']
		source_list.append(int(SPOT_SOURCE_ID))
	mit_spots = [] 		
	sorted_source_list = sorted(source_list)
	#print('sorted_source_list', sorted_source_list)
	for i in range (len(sorted_source_list)):
		if sorted_source_list[i]==sorted_source_list[i-1]:
			mit_spots.append(sorted_source_list[i])
	print('coloring spots downward ...')
	spot_color_downward (mit_spots,edge_list,n_spots)
	
def spot_color_downward (mit_spots,edge_list,n_spots): 	
	print('function stadted: spot_color_downward')
	counter = 0
	for Spot in root.iter('Spot'):
		Spot_dict = Spot.attrib
		if 'MANUAL_COLOR' in Spot_dict and int(Spot_dict['ID']) not in mit_spots:
			col = Spot_dict['MANUAL_COLOR']
			print('loop 1')
			for i in range (len(edge_list)):
				if edge_list[i][0] == int(Spot_dict ['ID']):
					target = edge_list[i][1]
					print('loop 2')
					for j in range (len(edge_list)):
						if edge_list[j][0] == target:
							print('loop 3')
							for Spot in root.iter('Spot'):
								Spot_dict = Spot.attrib
								if Spot_dict ['ID'] == str(target):
									Spot_dict['MANUAL_COLOR'] = col
									counter =+1
									print('spot_color_downward ',((counter/n_spots)*100),' done')
									spot_color_downward (mit_spots,edge_list,n_spots)


#col_adapt = input('Write colors from colored spots to target uncolored spots (e.g. after interpolation with MaMut) (y/n) ? >')
#if col_adapt == 'y':
#	adapt_colors ()




print ('Analyzing lineage ...')
# Get maximum number of cell generations from xml file
gens_tot = 0
for Spot in root.iter('Spot'):
	Spot_dict = Spot.attrib
	#print ('Spot_dict:',Spot_dict)
	if 'cell_generation' in Spot_dict:
		#print('cell_generation:',Spot_dict['cell_generation'],'spot_ID:',Spot_dict['ID'])
		ge = (Spot_dict['cell_generation'])
		int(ge)
		if 'cell_generation' in Spot_dict and int(ge) > gens_tot:
			gens_tot = int(ge)
		if 'cell_generation' not in Spot_dict:
			print('no generation info for spot',Spot_dict['ID'])

# Get minimum number of cell generations from xml file
gens_min = gens_tot
for Spot in root.iter('Spot'):
	Spot_dict = Spot.attrib
	if 'cell_generation' in Spot_dict:
		ge = int(Spot_dict['cell_generation'])
		if ge < gens_min: 
			gens_min = ge

# Number of generations passed within the recorded lineage:
gens = gens_tot - gens_min + 1

print ('gens_tot:',gens_tot)
print ('gens_min:',gens_min)
print ('generations passed in lineage: ', gens)

# Get sublineages from xml file
Sublins = []						# List of sublineage names as strings.
for Spot in root.iter('Spot'):
	Spot_dict = Spot.attrib
	if 'sublineage_name' in Spot_dict and Spot_dict['sublineage_name'] not in Sublins:
		Sublins.append(Spot_dict['sublineage_name'])



print('Longest track has ',gens_tot,' generations')
print('Found sublineages:',Sublins)


# Set background color

b = input ('Enter background color (color given to lineages that are not sublineages). To leave colors unchanged type "n" > ')
if b != 'n':
	if b in bas_cols:
		bc = col_dict[b]
	else:
		s_col = b.split(',')
		B = int(s_col[2])
		G = int(s_col[1])
		R = int(s_col[0])
		# RGB to Windows(decimal)
		WIN = (-16777216 + (R*65536) + (G*256) + B)
		bc = WIN
	print('Background color:',bc)
	# Apply background color to file:
	for Spot in root.iter('Spot'):
		Spot_dict = Spot.attrib
		Spot_dict['MANUAL_COLOR'] = str(bc)
else:
	bc = 'n'


col_gen = input('Color lineage by generations? (y/n) >')
col_slin = input('Color lineage by sublineages? (y/n) >')		

cols_slin = []						# Nested list mapping sublineage names to assigned colors

if col_slin == 'y' and col_gen != 'y':
	print ('To choose a color from a small set of basic colors type "r" (red), "g" (green), "b" (blue)')
	print ('"y" (yellow), "o" (orange), "c" (cyan), "m" (magenta), bl" (black), "w" (white), "gr" (grey)')
	print ('To specify a clolor manually enter the values (0-255) for each channel, separated by commas')
	print('E.g. for "purple" enter "128,0,128"')

	for i in range (len(Sublins)):
		col = []
		print('name of sublineage:',Sublins[i])
		c = input('Choose color for this sublineage >')
		col.append(Sublins[i])
		col.append(c)
		cols_slin.append(col)
	print('cols_slin:',cols_slin)

	col_trans (cols_slin, col_dict)

	print('cols_slin_new:',cols_slin)



def update_xml_colors_sublin (cols_slin,b,bc):
	for Spot in root.iter('Spot'):
		Spot_dict = Spot.attrib
		if 'sublineage_name' in Spot_dict:
			for i in range(len(cols_slin)):
				if ['sublineage_name'] == cols_slin [i][0]:
					del Spot_dict ['MANUAL_COLOR']
					Spot_dict ['MANUAL_COLOR'] = cols_slin [i][1]				
		else:
			if b != 'n' and bc !='n':
				Spot_dict ['MANUAL_COLOR'] = str(bc)
update_xml_colors_sublin(cols_slin,b,bc)	


### !!! new and not yet working ... 
# calgulate color ranges within sublineages:
cols_slin_gen = []
if col_slin == 'y' and col_gen == 'y':
	print ('To choose a color from a small set of basic colors type "r" (red), "g" (green), "b" (blue)')
	print ('"y" (yellow), "o" (orange), "c" (cyan), "m" (magenta), bl" (black), "w" (white), "gr" (grey)')
	print ('To specify a clolor manually enter the values (0-255) for each channel, separated by commas')
	print('E.g. for "purple" enter "128,0,128"')

	for i in range (len(Sublins)):
		col = []
		print('name of sublineage:',Sublins[i])
		cr = input('Choose root color for this sublineage >')
		cl = input('Choose leef color for this sublineage >')
		col.append(Sublins[i])
		col.append(cr)
		col.append(cl)
		cols_slin_gen.append(col)
	print('cols_slin_gen:',cols_slin_gen)

# create color range for each sublineage:

cols_slin_gen_ran = []
create_color_range_sublin (gens,col_range_win,cols_slin_gen,cols_slin_gen_ran)
print('cols_slin_gen_ran',cols_slin_gen_ran)

def update_xml_colors_sublin_gen (cols_slin_gen_ran,b,bc):
	for Spot in root.iter('Spot'):
		Spot_dict = Spot.attrib
		if 'sublineage_name' in Spot_dict and 'cell_generation' in Spot_dict:
			#print ('found Spot')
			for i in range(len(cols_slin_gen_ran)):
				if Spot_dict['sublineage_name'] == str(cols_slin_gen_ran [i][0]):
					#print ('sublineage_name:',Spot_dict['sublineage_name'])
					c = int (Spot_dict['cell_generation'])
					col_list = cols_slin_gen_ran [i][1]
					if c == 9:
						print('c:',c)
						print('col_list',col_list)
						#print('col_list[9]',col_list[9])
					#print('cols_slin_gen_ran',cols_slin_gen_ran)		
					for j in range (len(col_list)):
						if j == c:
							Spot_dict ['MANUAL_COLOR'] = str(col_list[j])
		else:
			print ('else')
			if b != 'n' and bc !='n':
				Spot_dict ['MANUAL_COLOR'] = str(bc)
update_xml_colors_sublin_gen (cols_slin_gen_ran,b,bc)	




# Inquire parameters for color coding lineage by generation:
if col_gen == 'y' and col_slin != 'y':
	csl = 0
	for Feature in root.iter('Feature'):
		Feat_dict = Feature.attrib
		if Feat_dict['feature'] == 'sublineage_name':
			csl = 1
			break
	if csl == 1:
		col_gen_sl = input('Color-code generations in selected sublineages (not yet implemented) (y/n)? >') #not yet implemented
	col_gen_tot = input('Color-code generations in entire Lineage (y/n)?')

	print ('To choose a color from a small set of basic colors type "r" (red), "g" (green), "b" (blue)')
	print ('"y" (yellow), "o" (orange), "c" (cyan), "m" (magenta), bl" (black), "w" (white), "gr" (grey)')
	print ('To specify a clolor manually enter the values (0-255) for each channel, separated by commas')
	print('E.g. for "purple" enter "128,0,128"')

	if col_gen_tot == 'y':
		ro = input('enter root color >')
		le = input('enter leaf color >')
		cols = [ro,le]
		create_color_range (cols,gens,col_range_win)	#!!! There can be eerors when calculating a range including cyan AND blue (both calculated as blue)
		print('col_range_win:',col_range_win)

		for Spot in root.iter('Spot'):
			Spot_dict = Spot.attrib
			#if 'MANUAL_COLOR' in Spot_dict:				# This means that all existing colors are deleted, including background color.
			#	del Spot_dict ['MANUAL_COLOR']
			if 'cell_generation' in Spot_dict:
				g = int(Spot_dict['cell_generation'])
				#print('generation:',g)		# To DO later: Above add code that checks if cell generations are specified and tells the user to run "lineage from mamut" if not.  
				for i in range (len(col_range_win)):
					if i == g - gens_min:
						c_i = col_range_win[i]
						Spot_dict['MANUAL_COLOR'] = str(col_range_win[i])
						# print('col_range_win[i]:',col_range_win[i])	
			 	
			# print('Spot_dict[MANUAL_COLOR]',Spot_dict['MANUAL_COLOR'])


# For entire xml-object set edge-color to the color of the target spot. 
for Spot in root.iter('Spot'):
	Spot_dict = Spot.attrib
	if 'MANUAL_COLOR' in Spot_dict:
		for Edge in root.iter('Edge'):
			Edge_dict = Edge.attrib
			if Edge_dict['SPOT_TARGET_ID'] == Spot_dict['ID']:
				Edge_dict['MANUAL_COLOR'] = Spot_dict['MANUAL_COLOR']



# Remove spot-and edge- Information from root.object (leaving only changed color). 
# This is necessary if changes ti this file in MaMut are to be saved
def remove_new_data ():
	for Spot in root.iter('Spot'):
		Spot_dict = Spot.attrib
		if 'cell_ID' in Spot_dict:
			del Spot_dict ['cell_ID']
		if 'cell_generation' in Spot_dict:
			del Spot_dict ['cell_generation']
		if 'sublineage_name' in Spot_dict:
			del Spot_dict ['sublineage_name']

#!!! Does not work yet
#	for Feature in root.iter ('Feature'):
#		F_dict = Feature.attrib
#		print (F_dict)
#		if F_dict ['feature'] == 'cell_ID':
#			SpotFeatures.remove (Feature)
#		if F_dict ['feature'] == 'cell_generation':
#			SpotFeatures.remove (Feature)
#		if F_dict ['feature'] == 'sublineage_name':
#			SpotFeatures.remove (Feature)

rem_new_data = input('Remove spot-and edge- Information from root.object (leaving only changed color). This is necessary if changes ti this file in MaMut are to be saved. (y/n) >')
if rem_new_data == 'y':
	remove_new_data()


# Test colorcode:
'''
for Spot in root.iter('Spot'):
	Spot_dict = Spot.attrib
	print (Spot_dict['cell_generation'])
	print (Spot_dict['MANUAL_COLOR'])
	print ('Spot_ID:',Spot_dict['ID'])
	for i in range (len(col_range_win)):
		if i == int(Spot_dict['cell_generation']) - gens_min: 
			c_i = col_range_win[i]
	if str(c_i) == Spot_dict['MANUAL_COLOR']:
		print('correct')
	if str(c_i) != Spot_dict['MANUAL_COLOR']:
		print('WRONG!')
'''

#print('cols_n:',cols_n)

export_name = get_file +'_color.xml'
export_mamut.write(export_name)
print ('saved file ',export_name)

#To Do: Plot Color code for generations, sublineages, etc. 



