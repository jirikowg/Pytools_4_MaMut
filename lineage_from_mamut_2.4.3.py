
'''
Tool to generate cell IDs and cell generations from Mamut-xml files. 
a new mamut-file is created where this information is added to the spot attributes. 
'''

from lxml import etree as ET
from io import StringIO, BytesIO

get_file = input('enter name of mamut-xml file without suffix. >')
get_file_suff = get_file + '.xml'
export_mamut = ET.parse(get_file_suff)

root = export_mamut.getroot()

# Get IDs of spots at (a) beginning of cells, (b) mitoses and (c) end of cells. 
# This is done by comparing lists: the list of sourches, targets and spots. 

#def spots_from_xml ():

source_list=[]
target_list=[]
for Edge in root.iter('Edge'):
	SPOT_SOURCE_ID = Edge.attrib['SPOT_SOURCE_ID']
	source_list.append(int(SPOT_SOURCE_ID))
	SPOT_TARGET_ID = Edge.attrib['SPOT_TARGET_ID']
	target_list.append(int(SPOT_TARGET_ID))


spot_list =[]
for Spot in root.iter('Spot'):
	Spot_ID = Spot.attrib['ID']
	spot_list.append(int(Spot_ID))

start_spots = [] 				# Spots which are not listed as targetof an edge and are therefore at the beginning of a track (or single spots). 
for Spot in root.iter('Spot'):
	Spot_ID = int(Spot.attrib['ID'])
	if Spot_ID not in target_list:
		start_spots.append(Spot_ID)
#print ('start_spots:',start_spots)

end_spots = []					# # Spots which are not listed as sources in the edges are at the end of a track (or single spots). 
for Spot in root.iter('Spot'):
	Spot_ID = int(Spot.attrib['ID'])
	if Spot_ID not in source_list:
		end_spots.append(Spot_ID)
#print ('end_spots:',end_spots)

# Get mitosis spots (last spot before  division) by iterating throug the source list. If the spot exists in 
# the source list TWICE it is a mitosis spot.

# !!!!!!! TO Do: allow more than one target cell (allow polytomies). !!!!!!




mit_spots = [] 		

sorted_source_list = sorted(source_list)
#print('sorted_source_list', sorted_source_list)
for i in range (len(sorted_source_list)):
	if sorted_source_list[i]==sorted_source_list[i-1]:
		mit_spots.append(sorted_source_list[i])

# Results:
#print ('start_spots:',start_spots)
#print('mit_spots',mit_spots)
#print ('end_spots:',end_spots)

# List of all spots that mark beginning or end of a cell.
cell_end_spots = mit_spots + start_spots + end_spots 	
#print ('cell_end_spots:',cell_end_spots)

mit_end_spots = mit_spots + end_spots 			# New in 2.4.1, this list makes "conect_edges_up" faster.
print ('mit_end_spots:',mit_end_spots)


# Make nested list that maps source to target for each edge.
edge_list = []
rev_edge_dict = {}			# Maps target_ID (Key) to SOurce_ID (value) for each edge. 
for Edge in root.iter('Edge'):
	edge = []
	edge.append(int(Edge.attrib['SPOT_SOURCE_ID']))
	edge.append(int(Edge.attrib['SPOT_TARGET_ID']))
	edge_list.append(edge)
	rev_edge_dict[int(Edge.attrib['SPOT_TARGET_ID'])] = int(Edge.attrib['SPOT_SOURCE_ID'])
#print('edge_list:',edge_list)
print('rev_edge_dict:',rev_edge_dict)

n_spots = len(spot_list)
n_edges = len(edge_list)
#print('n_spots:',n_spots)
#print ('n_edges:',n_edges)

#spots_from_xml ()

# Function:
# gets the source- and target ID for the end spots of each cell (end_spots, mit_spots) 
# creates a nested list called cell in which source- and target-IDS of the spots on the 
# branch are saved. 



cells = []			 # List of spot IDs on one branch, including start- and end spot. 

def connect_edges_up (edge_list,cell_end_spots):
	for i in range(len(mit_end_spots)):
		#print((i/len(mit_end_spots))*100)								# Counter
		#if edge_list[i][1] in cell_end_spots:
		print('mit_end_spots[i]:',mit_end_spots[i])
		if mit_end_spots [i] not in start_spots:		#### Problematic Case: Fitst spot in track is a mitosis spot!!!
			end = mit_end_spots [i]
			start = rev_edge_dict[end]
			cell = [end,start]
			extend_edges (edge_list,start,end,cell,cells)
	return start, end, cell, cells


# connetcs the edges by iterating upward over the branch until a amitosis-spot or start-spot is reached. 
def extend_edges (edge_list,start,end,cell,cells):
	for j in range(len(edge_list)):
		if cell in cells:
			break
		if start in cell_end_spots:
			cells.append(cell)
			break
		if edge_list[j][1] == start:
			start = edge_list [j][0]
			cell.append(start)

			extend_edges (edge_list,start,end,cell, cells)
	return start, end, cell, cells

connect_edges_up (edge_list,cell_end_spots)
print('cells:',cells)


# Make nested list of cells, assign cell_IDs and calculate cell generation by iterating from 
# each cell to the root of the tree. Also add ID of track-start spot. 

# Format: [[[Spots of cell],generation(gen),cell_ID, start-spot_ID],...,...]


cells_ID_gen = [] 	# nested list of cells, cell generation, cell_ID and start spot of the respective track.
cell_ID = 0

def get_cells (cells,cells_ID_gen,cell_ID,start_spots):
	for i in range (len(cells)):
		print('get_cells:', (i/(len(cells))*100))		# Counter
		gen = 0
		cell_ID = cell_ID +1
		track_start_spot = 0
		cell = [cells[i],gen,cell_ID,track_start_spot]
		start = cells [i][-1]
		if cells[i][-1] in start_spots:
			cell[3] = cells [i][-1]			# Track start spot
			cells_ID_gen.append(cell)
		else:
			extend_cells (cells, cell, gen, start, start_spots)
			cells_ID_gen.append(cell)
	return cells_ID_gen

def extend_cells (cells, cell, gen, start, start_spots):
	for j in range (len(cells)):
		if cells [j][0] == start and cells [j][-1] in start_spots:		# This will cause trouble when there are more than one track. Change !!!!!
			gen = gen +1
			cell[1] = gen
			cell[3] = cells [j][-1]
			break
		if cells [j][0] == start:
			gen = gen +1
			start = cells [j][-1]
			extend_cells (cells, cell, gen, start, start_spots)	
	return cells, cell, start, start_spots

get_cells (cells,cells_ID_gen,cell_ID,start_spots)
print('cells_ID_gen:',cells_ID_gen)

		
# Correct generations with start generation for each track. !!! To Do: Make this optional so users can use only the sublineage option if they like. 
start_spot_gens = []
def get_start_gens (start_spots): 		# Makes nested list of start spots of tracks and respective start generations. 
	#Inquire start generation of tracks from the user
	'''TO DO: Make option: Same genearation for every track. (In case lineage has a large number of tracks).'''
	print('... found',len(start_spots),'tracks.')
	print('Enter the start generations for each track of your lineage.')
	print('e.g. if the lineage starts with the cygote enter "0", if track 1 starts after second division (4 cell stage) enter "2", ...)')
	for i in range (len(start_spots)):
		print('Track has start spot ',start_spots[i])
		start_gen = input('Please enter start generation of this track >')
		start_gen = int(start_gen)
		start_spot_gen = [start_spots[i],start_gen]
		start_spot_gens.append(start_spot_gen)
	return start_spot_gen
get_start_gens (start_spots)


# Correct generation values in list "cells_ID_gen" by adding start-generations. 
def correct_gen(cells_ID_gen,start_spot_gens):
	for i in range (len(start_spot_gens)):
		print('correct_gen:',(i/len(start_spot_gens))*100)						# Counter
		for j in range (len(cells_ID_gen)):
			if cells_ID_gen[j][3] == start_spot_gens[i][0]:
				cells_ID_gen[j][1] = cells_ID_gen[j][1] + start_spot_gens[i][1]
	return cells_ID_gen
correct_gen(cells_ID_gen,start_spot_gens)
print('cells_ID_gen, corrected:',cells_ID_gen)

# Remove overlap between start spots and end spots in list "cells_ID_gen":
# If end_spot of a cell has same ID ans Start spot of daughter cell: Delete start Spot of daughter cell. 
# This will assure that every spot has only one generation value. 

cells_ID_gen_nodup = []
def remove_overlaps (cells_ID_gen,mit_spots,cells_ID_gen_nodup):	
	for i in range(len(cells_ID_gen)):
		print('remove_overlaps:',(i/len(cells_ID_gen))*100)
		if cells_ID_gen[i][0][-1] in start_spots and cells_ID_gen[i][0][0] in end_spots:
			cells_ID_gen_nodup.append(cells_ID_gen[i])
			print ('Cell with startspot and endspot',cells_ID_gen[i][0])
		if cells_ID_gen[i][0][-1] in start_spots and cells_ID_gen[i][0][0] in mit_spots:
			cells_ID_gen_nodup.append(cells_ID_gen[i])
		for j in range(len(cells_ID_gen)):
			if cells_ID_gen[i][0][0] in mit_spots and cells_ID_gen[i][0][0] == cells_ID_gen[j][0][-1] and cells_ID_gen[j] not in cells_ID_gen_nodup:
				nu_cell = cells_ID_gen[j][0][0:-1]
				nu_cell_ID_gen = [nu_cell, cells_ID_gen[j][1],cells_ID_gen[j][2],cells_ID_gen[j][3]]
				cells_ID_gen_nodup.append(nu_cell_ID_gen)
remove_overlaps (cells_ID_gen,mit_spots,cells_ID_gen_nodup)
print ('cells_ID_gen_nodup',cells_ID_gen_nodup)





# 		Assign sublineages 
print ('Would you like to assign identities to sublineages ?') 
print ('e.g. You could assign parts of the lineage to Quadrants A, B, C, or D, the 4d lineage in spiral cleavage, or any other category that you wish. ')
l = input ('specify sublineages  (y/n)? >')
print('l:',l)

def get_sublins (sublins):
	ans = input('Enter a(nother) sublineage (y/n)? >')
	if ans == 'y':
		sublin = []
		name = input('enter name >')
		stid = input('Enter start_spot_ID of sublineage >')
		sublin.append(name)
		sublin.append(stid)
		sublins.append(sublin)
		get_sublins (sublins)


if l == 'y':
	n_tracks = 0
	for Track in root.iter('Track'):
		n_tracks += 1
	print ('found ',n_tracks,' tracks')
	print ('start_spot_IDs of lineage:', start_spots)
	sublins = []
	get_sublins (sublins)
	print('sublins:',sublins)



def get_sublin_spots (starts, sublin_spots,edge_list, end_spots):
	if len(starts) != 0:
		for j in range(len(edge_list)):	
			if len(starts) == 0:
				break
			for k in range(len(starts)):
				if starts[k] in end_spots:
					starts.remove(starts[k])
					break
			if len(starts) == 0:
				break
			for k in range(len(starts)):
				if len(starts) != 0:
					if edge_list [j][0] == (starts [k]):
						starts.append (edge_list [j][1])
						sublin_spots.append (edge_list [j][1])
						for l in range (len(edge_list)):
							if edge_list [l][0] == edge_list [j][0] and edge_list [l][1] != edge_list [j][1]:
								starts.append (edge_list [l][1])
								sublin_spots.append (edge_list [l][1])
						starts.remove(starts[k])
						get_sublin_spots (starts, sublin_spots,edge_list, end_spots)
	return sublin_spots, starts
	print ('sublin_spots:',sublin_spots)

sublin_list = []
def assign_sublineages (sublins, edge_list, sublin_list):
	for i in range(len(sublins)):
		print('assign_sublineages:',(i/len(sublins))*100)
		sublin_spots = [int(sublins[i][1])]
		starts = [int(sublins[i][1])]
		get_sublin_spots (starts, sublin_spots, edge_list, end_spots)
		name_sublin = [sublins[i][0],sublin_spots]
		sublin_list.append(name_sublin)
	return sublin_list
if l == 'y':
	assign_sublineages (sublins, edge_list, sublin_list)

	print('sublin_list:',sublin_list)


# Add Cell ID and Cell generation to spot-attributes in mamut-xml-file

def update_sublin_xml (sublin_list):
	for Spot in root.iter('Spot'):
		Spot_ID = int(Spot.attrib['ID'])
		Spot_dict = Spot.attrib
		for i in range (len(sublin_list)):
			if Spot_ID in sublin_list [i][1]:
				Spot_dict ['sublineage_name'] = str(sublin_list [i][0])	# Adds information to xml
if l == 'y':
	update_sublin_xml (sublin_list)



def update_cell_gen_xml (cells_ID_gen_nodup, get_file):
	for Spot in root.iter('Spot'):
		Spot_ID = int(Spot.attrib['ID'])
		Spot_dict = Spot.attrib
		for i in range (len(cells_ID_gen_nodup)):
			if Spot_ID in cells_ID_gen_nodup [i][0]:
				Spot_dict ['cell_generation'] = str(cells_ID_gen_nodup [i][1])
				Spot_dict ['cell_ID'] = str(cells_ID_gen_nodup [i][2])				# "cell generation" = number of divisions from root.
update_cell_gen_xml (cells_ID_gen_nodup,get_file)

# Update_spot_features_in_xml_file :
# The features to be updated depend on the changes made to the spot attributes, e.g. only cell ID and Cell generation added, or also subleanage, ... 

# Check if 
cids = 0
cgens = 0
slins = 0
for Spot in root.iter('Spot'):
	Spot_dict = Spot.attrib

	if 'cell_ID' in Spot_dict:
		cids += 1
	if 'cell_generation' in Spot_dict:
		cgens +=1
	if 'sublineage_name' in Spot_dict:
		slins +=1	

#print('cids:',cids)
#print('cgens:',cgens)
#print('slins:',slins)

# Insert elements in xml tree object:
for SpotFeatures in root.iter('SpotFeatures'):
	if cids > 0:
		cid = ET.Element("Feature", feature="Cell_ID", name="Cell_ID", shortname="Cell_ID", dimension="NONE", isint="true") 
		# ET.tostring(f, pretty_print=True)				# Doesn´t work. xml data is not serialized at inserted position. 
		SpotFeatures.insert(0,cid)
	if cgens > 0:
		cgen = ET.Element("Feature", feature="cell_generation", name="cell_generation", shortname="cell_generation", dimension="NONE", isint="true") 
		# ET.tostring(f, pretty_print=True)				# Doesn´t work. xml data is not serialized at inserted position. 
		SpotFeatures.insert(0,cgen)
	if slins > 0:
		sln = ET.Element("Feature", feature="sublineage_name", name="sublineage_name", shortname="sublineage_name", dimension="NONE", isint="false") 
		# ET.tostring(f, pretty_print=True)				# Doesn´t work. xml data is not serialized at inserted position. 
		SpotFeatures.insert(0,sln)
	#sln = ET.Element("Feature", feature="sublineage_name", name="sublineage_name", shortname="sublineage_name", dimension="NONE", isint="false") 
	# ET.tostring(f, pretty_print=True)				# Doesn´t work. xml data is not serialized at inserted position. 
	#SpotFeatures.insert(0,sln)

#	Write data to xml
export_name = get_file +'_lin.xml'	
export_mamut.write(export_name)
print ('saved file ',export_name)

