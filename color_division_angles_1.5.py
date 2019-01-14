# THIS PROGRAM COLOR-CODES ORTHOGONALITY OF CELL DIVISIONS IN MAMUT CELL LINEAGES.
# THE PROGRAM CALCULATES THE ANGLE BETWEEN TWO SUBSEQUEND MITOSES AND APPLIES A USER DEFINED COLOR CODE TO THE 
# BRANCHES ESXTENDING FROM THE MOTHER SPOT TO THE TWO DAUGHTER SPOTS. 
# THE LINEAGE TREE CAN BE MODIFIED BY REMOVING ALL SPOTS BETWEEN DAUGHTER CELLS AND THE FOLLOWING MITOSES. 
# THIS CAN BE USEFUL BECAUSE EXCESSIVE CELL MOVEMENT BETWEEN MITOSES CAN HAVE AN IMPACT ON VISIBILITY OF THE CELL TRACKS. 

'''
TO DO:
Modify code so that saving file in MaMuT works.

PLOT COLOR RANGE

FIX OPTION TO SIMPLIFY BRANCHES BETWEEN MITOSES

Make program faster!

Make scaling optional

'''

import math
import numpy as np
from lxml import etree as ET
from io import StringIO, BytesIO
from colour import Color

get_file = input('enter name of mamut-xml file without suffix. >')
get_file_suff = get_file + '.xml'
export_mamut = ET.parse(get_file_suff)
root = export_mamut.getroot()

### get time interval
for ImageData in root.iter('ImageData'):
	time_interval = float(ImageData.attrib['timeinterval'])
print('time_interval:',time_interval)


print('')
print('')
print('')
print('Caution! Divisions must have exactly 2 daughter cells!')
print('Caution! The program uses the first offspring of a mitosis to calculate it´s angle!')
print('Caution! In the current version calculated angle information will only be valid if the lineage data are in scale (e.g. to voxel size)!')
print('The h5xml file should be provided in the working directory so scaling information can be read!')
print('')
print('')
print('')

'''
TO DO: ALLOW INCORPORATION OF MULTIPLE FRAMES AFTER MITOSIS TO CALCULATE MEAN DAUGHTER CELLS. 
THIS CAN BE NECESSARY BECAUSE CELLS CAN ROTATE DURING MITOSIS, ALTERING THE ANGLES.
CALCULATING MEAN DAUGHTER SPOTS MAY REQUIRE PREVIOUS INTERPOLATION OF SPOTS
'''

# Get IDs of spots at (a) beginning of cells, (b) mitoses and (c) end of cells. 
# This is done by comparing lists: the list of sources, targets and spots. 
source_list=[]
target_list=[]
for Edge in root.iter('Edge'):
	SPOT_SOURCE_ID = Edge.attrib['SPOT_SOURCE_ID']
	source_list.append(int(SPOT_SOURCE_ID))
	SPOT_TARGET_ID = Edge.attrib['SPOT_TARGET_ID']
	target_list.append(int(SPOT_TARGET_ID))


spot_list =[]
start_spots = [] 
end_spots = []
frame_dict = {}
for Spot in root.iter('Spot'):
	Spot_ID = int(Spot.attrib['ID'])
	spot_frame = int(Spot.attrib['FRAME'])
	spot_list.append(Spot_ID)
	frame_dict[Spot_ID] = spot_frame
	if Spot_ID not in target_list:
		start_spots.append(Spot_ID)
	if Spot_ID not in source_list:
		end_spots.append(Spot_ID)
#print('frame_dict:',frame_dict)

# Get mitosis spots (last spot before  division) by iterating throug the source list. If the spot exists in 
# the source list TWICE it is a mitosis spot.
mit_spots = []					#List of mother spots at mitoses
mit_list = []					#List of mitoses [[mother spot, first daughter spot, second daughter spot], [,,],[,,],.... ]	
edge_list = []					#Nested list mapping source ID and Target ID for each edge. 
rev_edge_dict = {}			# Maps target_ID (Key) to SOurce_ID (value) for each edge.
dau_spots = []					#List of Spots which are direct daughters of a mitosis.

for Edge in root.iter('Edge'):
	edge = []
	edge.append(int(Edge.attrib['SPOT_SOURCE_ID']))
	edge.append(int(Edge.attrib['SPOT_TARGET_ID']))
	edge_list.append(edge)
	rev_edge_dict[int(Edge.attrib['SPOT_TARGET_ID'])] = int(Edge.attrib['SPOT_SOURCE_ID'])


sorted_edge_list = sorted(edge_list)

mit = []
for i in range (len(sorted_edge_list)):										
	if sorted_edge_list[i][0] == sorted_edge_list[i-1][0]:
		mit_spots.append(sorted_edge_list[i][0])
		dau_spots.append(sorted_edge_list[i][1])
		dau_spots.append(sorted_edge_list[i-1][1])
		mit.append(sorted_edge_list[i][0])
		mit.append(sorted_edge_list[i][1])
		mit.append(sorted_edge_list[i-1][1])
		mit_list.append(mit)
		mit = []


# List of all spots that mark beginning or end of a cell.
cell_end_spots = mit_spots + start_spots + end_spots 	
#print ('cell_end_spots:',cell_end_spots)
#print('mit_spots:',mit_spots)

#mit_end_spots = mit_spots + end_spots 			# New in "lineage_from_mamut" 2.4.1, this list makes "conect_edges_up" faster.
#print ('mit_end_spots:',mit_end_spots)
	

n_spots = len(spot_list)
n_edges = len(edge_list)
#print('n_spots:',n_spots)
#print ('n_edges:',n_edges)

# Get the source- and target ID for the end spots of each cell (end_spots, mit_spots) 
# creates a nested list called cell in which source- and target-IDS of the spots on the branch are saved. 
cells = []			 # List of spot IDs on one branch, including start- and end spot. 
print('collecting cells ...')
def connect_edges_up (edge_list,cell_end_spots,mit_spots):
	for i in range(len(edge_list)):
		#print(int((i/n_edges)*100),'%')								# Counter
		if edge_list[i][1] in cell_end_spots:
			start = edge_list [i][0]
			end = edge_list [i][1]
			cell = [end,start]
			#print('cell:',cell)
			#print('grand_mit:',grand_mit)
			extend_edges (edge_list,start,end,cell,cells,mit_spots,start_spots)
	return start, end, cell, cells

def extend_edges (edge_list,start,end,cell,cells,mit_spots,start_spots):	# connetcs the edges by iterating upward over the branch until a amitosis-spot or start-spot is reached. 
	for j in range(len(edge_list)):
		if cell in cells:
			#print('cell in cells')
			break
		if start in cell_end_spots:
			cells.append(cell)
			#print('cell:',cell)
			break
		if edge_list[j][1] == start:
			start = edge_list [j][0]
			cell.append(start)
			extend_edges (edge_list,start,end,cell,cells,mit_spots,start_spots)
	return start, end, cell, cells

connect_edges_up (edge_list,cell_end_spots,mit_spots)
#print('cells:',cells)

print('Making mitosis-dictionaries ...')
rev_mit_dict = {} 		### Dictionary that maps mother spots to daugter spots of mitosis
rev_dau_dict = {}		### Dictionary that maps daughter spots of mitoses to mother spot of the next mitosis
for i in range(len(cells)):
	if cells[i][-1] in mit_spots:
		rev_mit_dict[cells[i][-2]] = cells[i][-1]		#Using daughters as values causes problems: Two identical keys (mother spot ID)
		if cells[i][0] not in end_spots:
			rev_dau_dict[cells[i][0]] = cells[i][-2]
#print('rev_mit_dict:',rev_mit_dict)
#print('rev_dau_dict:',rev_dau_dict)

grand_mit_list = []		# List of mother mitoses and their daughter mitoses [[grandmother, mother1, mother2], [mother, daughter1, daughter2]]
print('collecting mitoses ...')
def make_grand_mit_list(start_spots,cells,mit_spots,mit_list,grand_mit_list):
	for i in range(len (cells)):
		start = cells[i][-1]
		end = cells[i][0]
		grand_mit = []
		if start in mit_spots and end in mit_spots: 
			for j in range(len(mit_list)):
				if start == mit_list[j][0]:
					grand_mit.append(mit_list[j])
			for j in range(len(mit_list)):
				if end == mit_list[j][0]:
					grand_mit.append(mit_list[j])	
			grand_mit_list.append(grand_mit)	
	return grand_mit_list
	#print('grand_mit_list:',grand_mit_list)
make_grand_mit_list(start_spots,cells,mit_spots,mit_list,grand_mit_list)
#print('grand_mit_list:',grand_mit_list)
#print('')


print('creating vectors ...')
grand_vec_list = grand_mit_list 		# 	LIST WHERE DAUGHTER CELL IDS ARE REPLACED BY DAUGHTER CELL XYZ-COORDINATES.
def get_coordinates(grand_mit_list):
	# FOR EACH MITOSIS CALCULATE ANGLE BETWEEN ITSELF AND THE PREVIOUS MITOSIS
	# Coordinates are scaled using voxel size? NO!
	for i in range (len(grand_vec_list)):
		#print('getting angles from grand_mit_list ',int((i/len(grand_mit_list))*100),'%')
		for Spot in root.iter('Spot'):
			if	int(Spot.attrib['ID']) == grand_vec_list[i][0][1]:		
				x = float(Spot.attrib['POSITION_X'])#*scal_x
				y = float(Spot.attrib['POSITION_Y'])#*scal_y
				z = float(Spot.attrib['POSITION_Z'])#*scal_z
				A = [x,y,z]	
				grand_vec_list[i][0][1] = A
			if	int(Spot.attrib['ID']) == grand_vec_list[i][0][2]:		
				x = float(Spot.attrib['POSITION_X'])#*scal_x
				y = float(Spot.attrib['POSITION_Y'])#*scal_y
				z = float(Spot.attrib['POSITION_Z'])#*scal_z
				B = [x,y,z]	
				grand_vec_list[i][0][2] = B
			if	int(Spot.attrib['ID']) == grand_vec_list[i][1][1]:		
				x = float(Spot.attrib['POSITION_X'])#*scal_x
				y = float(Spot.attrib['POSITION_Y'])#*scal_y
				z = float(Spot.attrib['POSITION_Z'])#*scal_z
				Q = [x,y,z]	
				grand_vec_list[i][1][1] = Q
			if	int(Spot.attrib['ID']) == grand_vec_list[i][1][2]:		
				x = float(Spot.attrib['POSITION_X'])
				y = float(Spot.attrib['POSITION_Y'])
				z = float(Spot.attrib['POSITION_Z'])
				R = [x,y,z]	
				grand_vec_list[i][1][2] = R
	#print('grand_vec_list:',grand_vec_list)
	#print('')

	return grand_vec_list
get_coordinates(grand_mit_list)

#Calculate Angles
# To DO: Ommit Mitoses which do not have mother mitoses!
grand_ang_list = []		# Nested list of mitosis mother spot IDs and angles between mitosis and following mitosis. 
mit_ang_dict = {}		# Dictionary that maps mitosis spot ID to angle

print('calculating angles ...')
def get_angle (grand_vec_list,grand_ang_list,mit_ang_dict):
	# Turn coordinate lists into np arrays:
	for i in range(len(grand_vec_list)):
		mit_ang = [grand_vec_list[i][1][0]]		# List [mother cell ID]
		a = np.array(((grand_vec_list[i][0][1][0]),(grand_vec_list[i][0][1][1]),(grand_vec_list[i][0][1][2])))
		b = np.array(((grand_vec_list[i][0][2][0]),(grand_vec_list[i][0][2][1]),(grand_vec_list[i][0][2][2])))
		q = np.array(((grand_vec_list[i][1][1][0]),(grand_vec_list[i][1][1][1]),(grand_vec_list[i][1][1][2])))
		r = np.array(((grand_vec_list[i][1][2][0]),(grand_vec_list[i][1][2][1]),(grand_vec_list[i][1][2][2])))
		#print('a,b,q,r:',a,b,q,r)
		vec_ab = a-b
		vec_qr = q-r
		#print('vec_ab',vec_ab)
		#print('vec_qr',vec_qr)

		#Betrag der Vektoren:
		bet_vec_ab = math.sqrt (vec_ab[0]**2 + vec_ab[1]**2 + vec_ab[2]**2)
		bet_vec_qr = math.sqrt (vec_qr[0]**2 + vec_qr[1]**2 + vec_qr[2]**2)

		dot_alfa = np.dot(vec_ab,vec_qr)
		#print ('dot_alfa:',dot_alfa)
		#print('')

		# angles:
		rad_alfa = math.acos(dot_alfa/(bet_vec_ab*bet_vec_qr))
		alfa = math.degrees (rad_alfa)
		#print('alfa:',alfa,)
		#print('')

		# CALCULATE CORRECTED ANGLE, THAT SPECIFIES DEVIATION FROM 90 DEGREES (EG 95 --> 5)
		if alfa > 90:
			corr_alfa = 180 - alfa
		else: 
			corr_alfa = alfa
		mit_ang.append(int(corr_alfa))
		grand_ang_list.append(mit_ang)

		mit_ang_dict[grand_vec_list[i][1][0]] = int(corr_alfa)		# Maybe this will be useful later
	return grand_ang_list
get_angle (grand_vec_list,grand_ang_list,mit_ang_dict)
#print('grand_ang_list ([mit_spot ID, angle between mitosis andd previous mitosis]')
#print(grand_ang_list)
#print('mit_ang_dict:',mit_ang_dict)

# Add color set:
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


def set_background():
# Set background color for spots and edges
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
			WIN = (-16777216 + (R*65536) + (G*256) + B)			# Ausnahmen hunzufügen die in Biocell2XML probleme gemacht haben. 
			bc = WIN
		print('Background color:',bc)
		# Apply background color to file:
		for Spot in root.iter('Spot'):
			Spot_dict = Spot.attrib
			Spot_dict['MANUAL_COLOR'] = str(bc)
		for Edge in root.iter('Edge'):
			Edge_dict = Edge.attrib
			Edge_dict['MANUAL_COLOR'] = str(bc) 
	else:
		bc = 'n'
set_background()

# Inquire color settings from user
print('')
print('')
print('Cell divisions with at least one preceding division will be color coded in the MaMut_XML-file...')
print('')
print('The color code will be applied to the first edges extending from the mother spot ...')
print('')
print('To define the color range please enter start- and end colors for angles 0° and 90° respectively ...')
print('')
print ('To choose a color from a small set of basic colors type "r" (red), "g" (green), "b" (blue), "y" (yellow), "o" (orange), "c" (cyan), "m" (magenta), "bl" (black), "w" (white), "gr" (grey) ...')
print('')
print ('To specify a color manually enter the R-, G-, B-values (0-255) separated by commas (E.g. for "purple" enter "128,0,128"')
print('')
print('')

a1 = input ('Enter start color > ')
a2 = input ('Enter end color > ')
cols = [a1,a2]

'''col_dict = {}
col_dict['r'] = -65536
col_dict['g'] = -16711936
col_dict['b'] = -16776961
col_dict['y'] = -256
col_dict['o'] = -23296
col_dict['c'] = -16711681
col_dict['m'] = -65281
col_dict['bl' 'w'
'gr'] = -8421505
'''

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
	#print('cols:',cols)
	#print('cols_n:',cols_n)
	return cols,cols_n


#MODIFY:
col_range_win = []
def create_color_range (cols,col_range_win):			
	cols_n = []							# Colors transferred to rgb
	col_2_rgb(cols,cols_n)				# funcion
	c1 = Color(rgb=(cols_n[0][0],cols_n[0][1],cols_n[0][2]))
	c2 = Color(rgb=(cols_n[1][0],cols_n[1][1],cols_n[1][2]))
	global col_range
	col_range = list(c1.range_to(c2,90))		
	#print('col_range:',col_range)

	# extract r, g and b values from PI-color-format
	for i in range(len(col_range)):
		#print('(col_range[i]',col_range[i])
		l_red = int(float(col_range[i].red)*255)
		l_green = int(float(col_range[i].green)*255)
		l_blue = int(float(col_range[i].blue)*255)
		#print ('l_red:',l_red,'l_green:',l_green,'l_blue:',l_blue)
		col_win_i = (-16777216 + (l_red*65536) + (l_green*256) + l_blue)
		col_range_win.append(col_win_i)
	#print('col_range_win:',col_range_win)
	return col_range_win
create_color_range (cols,col_range_win)
#!!! There can be eerors when calculating a range including cyan AND blue (both calculated as blue). Also with other colors (yellow) ...

# color spot by angle
mit_col_dict = {}	# Dictionary of mitosis spot IDs and colors
def make_mit_col_dict(mit_col_dict,grand_ang_list,col_range_win):
	for i in range(len(grand_ang_list)):
		m_ID = grand_ang_list[i][0]
		m_col = col_range_win [int(grand_ang_list[i][1])]
		mit_col_dict[m_ID] = m_col
make_mit_col_dict(mit_col_dict,grand_ang_list,col_range_win)
#print('mit_col_dict:',mit_col_dict)

# Update spot colors in XML:
def update_spot_colors(mit_spots,col_range_win,grand_ang_list,mit_col_dict):
	# TO DO: Set background color here
	for Spot in root.iter('Spot'):
		Spot_dict = Spot.attrib
		if int(Spot_dict['ID']) in mit_col_dict:
			m_ID = int(Spot_dict['ID'])
			m_ind = int(mit_col_dict[m_ID])
			Spot_dict['MANUAL_COLOR'] = str(m_ind)
			#print ('Spot_ID',int(Spot_dict['ID']),'mit_col_dict[m_ID]:',mit_col_dict[m_ID])
update_spot_colors(mit_spots,col_range_win,grand_ang_list,mit_col_dict)

# Update color of edges that originate from color-coeded mitosis spots in root object. 
def update_edge_colors(mit_col_dict):
	# TO DO: Set background color here
	for Edge in root.iter('Edge'):
		Edge_dict = Edge.attrib
		if int(Edge_dict['SPOT_SOURCE_ID']) in mit_col_dict:
			m_ID = int(Edge_dict['SPOT_SOURCE_ID'])
			m_ind = int(mit_col_dict[m_ID])
			Edge_dict['MANUAL_COLOR'] = str(m_ind)
update_edge_colors(mit_col_dict)

print('Simplify lineage tree by removing cell anagenesis between first daughter spots')
simp = input('and next mother spots (y/n)? >')

# !!! DOES NOT WORK YET

def simplify_tree(rev_mit_dict,rev_dau_dict,mit_spots,frame_dict,time_interval):

	#### remove spots that are not involved in mitoses
	for SpotsInFrame in root.iter('SpotsInFrame'):
		for Spot in SpotsInFrame.iter('Spot'):
			if int(Spot.attrib['ID']) not in mit_spots and int(Spot.attrib['ID']) not in rev_mit_dict:
				 SpotsInFrame.remove(Spot)

	### remove empty SpotsInFrame- antries
	for AllSpots in root.iter('AllSpots'):
		for SpotsInFrame in AllSpots.findall('SpotsInFrame'):
			if SpotsInFrame.getchildren() == []:	#### Condition: Element does not have Children
				AllSpots.remove(SpotsInFrame)			
		spot_count = 0
		for Spot in AllSpots.iter('Spot'):
			spot_count += 1

	for AllSpots in root.iter('AllSpots'):		#### Write new number of spots
		#print(AllSpots.attrib)
		AllSpots.attrib['nspots'] = str(spot_count)
		#print('AllSpots.attrib[nspots]:',AllSpots.attrib['nspots'])

	for Track in root.iter('Track'):
		#### Remove edges that are not involved in mitoses. combine anagenesis to one single edge.
		track_spots = []
		edge_count = 0
		for Edge in Track.iter('Edge'):
			source = int(Edge.attrib['SPOT_SOURCE_ID'])		
			target = int(Edge.attrib['SPOT_TARGET_ID'])
			if not target in rev_mit_dict and target not in rev_dau_dict:
				Track.remove(Edge)

		for Edge in Track.iter('Edge'):
			source = int(Edge.attrib['SPOT_SOURCE_ID'])		
			target = int(Edge.attrib['SPOT_TARGET_ID'])
			if target in rev_dau_dict:
				#print('target',target)
				Edge.attrib['SPOT_SOURCE_ID'] = str(rev_dau_dict[target]) ### Edge duration has to be updated here to calculate longest Gap.

			
		#### Analyze modified track #####
		edge_count = 0
		lon_gap = 0
		source_frame = 1000000000000 	### n_frames should not be larger than 1.000 billion
		target_frame = 0
		
		for Edge in Track.iter('Edge'):
			edge_count += 1
			source = int(Edge.attrib['SPOT_SOURCE_ID'])		
			target = int(Edge.attrib['SPOT_TARGET_ID'])
			if source not in track_spots:
				track_spots.append(source) 
			if target not in track_spots:
				track_spots.append(target) 
			if frame_dict[source] <= source_frame:		# Find spot with lowest frame number
				source_frame = frame_dict[source]
			if frame_dict[target] >= target_frame:	# Find spot with highest frame number
				target_frame = frame_dict[target]
			gap = frame_dict[target] - frame_dict[source]
			if gap > lon_gap:
				lon_gap = gap

		Track.attrib['NUMBER_SPOTS'] = str(len(track_spots))
		Track.attrib['NUMBER_GAPS'] = str(edge_count)
		Track.attrib['LONGEST_GAP'] = str(lon_gap)
		Track.attrib['NUMBER_SPLITS'] = str(len([val for val in track_spots if val in mit_spots]))		#Intersection of track spots and mit spots		
		Track.attrib['LONGEST_GAP'] = str(lon_gap)
		Track.attrib['TRACK_START'] = str(source_frame * time_interval)
		Track.attrib['TRACK_STOP'] = str(target_frame * time_interval)
		Track.attrib['TRACK_DURATION'] = str((target_frame * time_interval)-(source_frame * time_interval))
		# !!! Track displacement are not yet updated  !!!!

if simp == 'y':
	simplify_tree(rev_mit_dict,rev_dau_dict,mit_spots,frame_dict,time_interval)

def angles_to_xml(grand_ang_list,mit_ang_dict):
	# Insert angle informatiion as elements in xml tree object:
	for SpotFeatures in root.iter('SpotFeatures'):
		if len(grand_ang_list) > 0:
			angle = ET.Element("Feature", feature="Mitosis_ANGLE", name="Mitosis_ANGLE", shortname="angle", dimension="NONE", isint="true") 
			# ET.tostring(f, pretty_print=True)				# Doesn´t work. xml data is not serialized at inserted position. 
			SpotFeatures.insert(0,angle)

	# Add angle to Spot Attributes of Mitosis Spots in XML
	print('mit_ang_dict:',mit_ang_dict)
	for Spot in root.iter('Spot'):
		ID = int(Spot.attrib['ID'])
		print('ID:',ID)
		if ID in mit_ang_dict:
			print('mit_ID:',ID)
			Spot.attrib ['Mitosis_ANGLE'] = str(mit_ang_dict[ID])
			print('angle:',Spot.attrib ['Mitosis_ANGLE'])
angles_to_xml(grand_ang_list,mit_ang_dict)

export_name = get_file +'_col_ang.xml'
export_mamut.write(export_name)
print ('saved file ',export_name)

# plot data.
plot = input('Plot angle distribution ? (y/n) >')
if plot == 'y':
	import matplotlib.pyplot as plt
	from matplotlib.colors import LinearSegmentedColormap

	ang_list = []
	for i in range(len(grand_ang_list)):
		ang_list.append(grand_ang_list[i][1])
	#print (ang_list)
	
	bins = 90#range[0,90]		#Creates 90 groups for Histogram. Less than 90 is not reasonable, since  90 is also the range applieg for color coding angles above. 
	#make custom colormap: 
	colors = []
	for i in range(len(col_range)):
		colors.append(col_range[i].rgb)
	#colors = col_range.rgb	#Testwise: R->G
	#print('colors:',colors)
	cmap_name = 'cmap_ang'
	#fig, axs = plt.subplots(2, 2, figsize=(6, 9))
	#fig.subplots_adjust(left=0.02, bottom=0.06, right=0.95, top=0.94, wspace=0.05)

	n_bins = []		#List of steps in serialized color map
	for i in range (len(range(bins+1))):
		n_bins.append(i)

	for n_bin in zip(n_bins):
		# Create the colormap
		cm = LinearSegmentedColormap.from_list(cmap_name, colors)

	# Plot histogram.
	n, bins, patches = plt.hist(ang_list, bins, color='green', edgecolor='black',linewidth=0.8)	### Copied from Tutorial, don´t understand Why scaling does not work.   
	bin_centers = 0.5 * (bins[:-1] + bins[1:])

	# scale values to interval [0,1]
	col = bin_centers - min(bin_centers)
	col /= max(col)

	for c, p in zip(col, patches):
		plt.setp(p, 'facecolor',cm(c))
	#print(cm)
	plt.title("Angle frequenzy distribution")
	plt.xlabel("Angle (°)")
	plt.ylabel("Frequency")
	plt.show()


