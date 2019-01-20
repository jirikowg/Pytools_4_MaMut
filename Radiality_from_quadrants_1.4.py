from lxml import etree as ET
from io import StringIO, BytesIO
import math
import numpy as np
import xlwt
import matplotlib.pyplot as plt

print('Calculates radiality of quadrant arrangement as described in Loose & Scholtz, Development, Genes, Evolution (in prep)')
print('WARNING: Calculating radiality from quadrants requires an interpolated MaMut lineage')
print('The lineage must have at least 4 Branches')
print('Quadrant identities have to be assigned in a previous step using LINEAGE_FROM_MAMUT.PY and use names A, B, C and D for the quadrants')

get_file = input('enter name of mamut-xml file (without suffix) for which radialityshould be calculated >')
get_file_suff = get_file + '.xml'
export_mamut = ET.parse(get_file_suff)
root = export_mamut.getroot()


def get_coordinates(X_A,Y_A,Z_A,X_B,Y_B,Z_B,X_C,Y_C,Z_C,X_D,Y_D,Z_D):

	# Collect coordinates: 
	if Spot.attrib['sublineage_name'] == 'A':
		xa = float(Spot.attrib['POSITION_X'])
		X_A.append(xa)
		ya = float(Spot.attrib['POSITION_Y'])
		Y_A.append(ya)
		za = float(Spot.attrib['POSITION_Z'])
		Z_A.append(za)
	if Spot.attrib['sublineage_name'] == 'B':
		xb = float(Spot.attrib['POSITION_X'])
		X_B.append(xb)
		yb = float(Spot.attrib['POSITION_Y'])
		Y_B.append(yb)
		zb = float(Spot.attrib['POSITION_Z'])
		Z_B.append(zb)
	if Spot.attrib['sublineage_name'] == 'C':
		xc = float(Spot.attrib['POSITION_X'])
		X_C.append(xc)
		yc = float(Spot.attrib['POSITION_Y'])
		Y_C.append(yc)
		zc = float(Spot.attrib['POSITION_Z'])
		Z_C.append(zc)
	if Spot.attrib['sublineage_name'] == 'D':
		xd = float(Spot.attrib['POSITION_X'])
		X_D.append(xd)
		yd = float(Spot.attrib['POSITION_Y'])
		Y_D.append(yd)
		zd = float(Spot.attrib['POSITION_Z'])
		Z_D.append(zd)
	return 	X_A,Y_A,Z_A,X_B,Y_B,Z_B,X_C,Y_C,Z_C,X_D,Y_D,Z_D

def get_angles(X_A,Y_A,Z_A,X_B,Y_B,Z_B,X_C,Y_C,Z_C,X_D,Y_D,Z_D,frame,ang_ABCD,ang_ADCB,ang_mean,quad_pos):
	
	#Calculate mean quadrant coordinates: 
	mean_X_A = sum(X_A)/len(X_A)
	quad_pos.append(mean_X_A)

	mean_Y_A = sum(Y_A)/len(Y_A)
	quad_pos.append(mean_Y_A)

	mean_Z_A = sum(Z_A)/len(Z_A)
	quad_pos.append(mean_Z_A)

	mean_X_B = sum(X_B)/len(X_B)
	quad_pos.append(mean_X_B)

	mean_Y_B = sum(Y_B)/len(Y_B)
	quad_pos.append(mean_Y_B)

	mean_Z_B = sum(Z_B)/len(Z_B)
	quad_pos.append(mean_Z_B)

	mean_X_C = sum(X_C)/len(X_C)
	quad_pos.append(mean_X_C)

	mean_Y_C = sum(Y_C)/len(Y_C)
	quad_pos.append(mean_Y_C)

	mean_Z_C = sum(Z_C)/len(Z_C)
	quad_pos.append(mean_Z_C)

	mean_X_D = sum(X_D)/len(X_D)
	quad_pos.append(mean_X_D )

	mean_Y_D = sum(Y_D)/len(Y_D)
	quad_pos.append(mean_Y_D)

	mean_Z_D = sum(Z_D)/len(Z_D)
	quad_pos.append(mean_Z_D)


	# Numpy-arrays from Coordinates: 
	A = np.array(((mean_X_A),(mean_Y_A),(mean_Z_A)))
	B = np.array(((mean_X_B),(mean_Y_B),(mean_Z_B)))
	C = np.array(((mean_X_C),(mean_Y_C),(mean_Z_C)))
	D = np.array(((mean_X_D),(mean_Y_D),(mean_Z_D)))

	# Calculate vectors between potential sister-quadrants for two alternative scenarios: 
	AB = B-A
	CD = C-D

	AD = D-A
	CB = C-B

	#Lengh of vectors:
	lAB = math.sqrt (AB[0]**2 + AB[1]**2 + AB[2]**2)
	lCD = math.sqrt (CD[0]**2 + CD[1]**2 + CD[2]**2)
	lAD = math.sqrt (AD[0]**2 + AD[1]**2 + AD[2]**2)
	lCB = math.sqrt (CB[0]**2 + CB[1]**2 + CB[2]**2)

	# Dot products: 
	dot_alfa = np.dot(AB,CD)

	dot_beta = np.dot(AD,CB)

	# angles:
	rad_alfa = math.acos(dot_alfa/(lAB*lCD))
	alfa = math.degrees (rad_alfa)

	rad_beta = math.acos(dot_beta/(lAD*lCB))
	beta = math.degrees (rad_beta)

	mean_ang = (alfa + beta) /2
	#print('alfa:',alfa,'beta:',beta)

	ang_ABCD.append(frame)
	ang_ABCD.append(alfa) 
	
	ang_ADCB.append(frame)
	ang_ADCB.append(beta)
	
	ang_mean.append(frame)
	ang_mean.append(mean_ang)

	#print ('quad_pos:',quad_pos)
		
	#print('ang_ABCD:',ang_ABCD,'ang_ABCD:',ang_ABCD,'ang_mean:',ang_mean)
	return ang_ABCD,ang_ADCB,ang_mean,quad_pos



#    MASTER FUNCTION:

# Make excel spreadsheet: 
radiality_table = xlwt.Workbook()
sheet = radiality_table.add_sheet("radiality_sheet")

sheet.write(0, 0, 'AXm') # row, column, value
sheet.write(0, 1, 'AYm')
sheet.write(0, 2, 'AZm')
sheet.write(0, 3, 'BXm')
sheet.write(0, 4, 'BYm')
sheet.write(0, 5, 'BZm')
sheet.write(0, 6, 'CXm')
sheet.write(0, 7, 'CYm')
sheet.write(0, 8, 'CZm')
sheet.write(0, 9, 'DXm')
sheet.write(0, 10, 'DYm')
sheet.write(0, 11, 'DZm')
sheet.write(0, 13, 'angle ABCD')
sheet.write(0, 14, 'angle ADCB')
sheet.write(0, 15, 'angle mean')

# Nested lists with angle for each frame:
angles_ABCD = []
angles_ADCB = []
angles_mean = []

counter = 0

for SpotsInFrame in root.iter('SpotsInFrame'):
	SpotsInFrame_dict = SpotsInFrame.attrib
	frame = int(SpotsInFrame.attrib['frame'])
	#print('frame:',frame)

	# Lists with frame and angle:

	# Count number of quadrants in frame: 
	quad_list = []
	for Spot in SpotsInFrame.iter('Spot'):
		if 'sublineage_name' in Spot.attrib and Spot.attrib['sublineage_name'] not in quad_list:
			quad_list.append(Spot.attrib['sublineage_name']) 
	
	if len(quad_list) == 4:
		counter += 1
		#print('frame:',frame)
		#print('counter:',counter)

		#print('frame:',frame)
		X_A = []
		Y_A = []
		Z_A = []
		X_B = []
		Y_B = []
		Z_B = []
		X_C = []
		Y_C = []
		Z_C = []
		X_D = []
		Y_D = []
		Z_D = []
		quad_pos = []

		for Spot in SpotsInFrame.iter('Spot'):
			# Collect coordinates for spots of each lineage:
			get_coordinates(X_A,Y_A,Z_A,X_B,Y_B,Z_B,X_C,Y_C,Z_C,X_D,Y_D,Z_D)
		#print('X_A:',X_A)
		
		ang_ABCD = []
		ang_ADCB = []
		ang_mean = []
		quad_pos =[]
		get_angles(X_A,Y_A,Z_A,X_B,Y_B,Z_B,X_C,Y_C,Z_C,X_D,Y_D,Z_D,frame,ang_ABCD,ang_ADCB,ang_mean,quad_pos)
		
		angles_ABCD.append(ang_ABCD)
		angles_ADCB.append(ang_ADCB)
		angles_mean.append(ang_mean)

		# Write quadrant coordinates and angles to excel sheet 
		#print ('quad_pos[0]:',quad_pos[0])
		
		for i in range (len(quad_pos)):
			sheet.write(counter, i, quad_pos[i])

		sheet.write(counter, 13, ang_ABCD[1])
		sheet.write(counter, 14, ang_ADCB[1])
		sheet.write(counter, 15, ang_mean[1])

		#print('ang_ABCD:',ang_ABCD,'ang_ABCD:',ang_ABCD,'ang_mean:',ang_mean)

#print('angles_ABCD:',angles_ABCD)
#print('angles_ADCB:',angles_ADCB)
#print('angles_mean:',angles_mean)



#print('sorted_source_list', sorted_source_list)

# Save excel file: 
radiality_table.save("radiality_table.xls") 

# prepare plot:
frames = []
for i in range (len(angles_ABCD)):
	frames.append(angles_ABCD[i][0])
only_angles_ABCD = []
for i in range(len(angles_ABCD)):
	only_angles_ABCD.append(angles_ABCD[i][1])
only_angles_ADCB = []
for i in range(len(angles_ADCB)):
	only_angles_ADCB.append(angles_ADCB[i][1])
only_angles_mean = []
for i in range(len(angles_mean)):
	only_angles_mean.append(angles_mean[i][1])


# Get mitosis spots and their frames for plotting:
source_list=[]
for Edge in root.iter('Edge'):
	SPOT_SOURCE_ID = Edge.attrib['SPOT_SOURCE_ID']
	source_list.append(int(SPOT_SOURCE_ID))

	
mit_spots = [] 		
sorted_source_list = sorted(source_list)
for i in range (len(sorted_source_list)):
	if sorted_source_list[i]==sorted_source_list[i-1]:
		mit_spots.append(sorted_source_list[i])
print('mit_spots:',mit_spots)

mit_frames = []
for SpotsInFrame in root.iter('SpotsInFrame'):
	SpotsInFrame_dict = SpotsInFrame.attrib
	frame = int(SpotsInFrame.attrib['frame'])
	for Spot in SpotsInFrame.iter('Spot'):
		if int(Spot.attrib['ID']) in mit_spots:
			#print('int(Spot.attrib[ID]) in mit_spots = TRUE')
			mit_frames.append(frame)
print ('mit_frames:',mit_frames)

#List of numbers of mitoses for each frame:
mits_freq = []
for i in range(len(frames)):
	n = 0
	if frames[i] in mit_frames:
		for j in range(len(mit_frames)):
			if mit_frames[j] == frames[i]:
				n += 1
	mits_freq.append(n)
print('mits_freq:',mits_freq)



# plot data.
plt.plot(frames, only_angles_ABCD, 'c',frames, only_angles_ADCB, 'm',frames, only_angles_mean, 'k', mits_freq,'go')
plt.axis([30, 300, 0, 90])
plt.show()











#print('edge_list:',edge_list)

