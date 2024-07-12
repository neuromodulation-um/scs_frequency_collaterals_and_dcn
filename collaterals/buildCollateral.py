# from DC_Cell_Adjusted import *
# from ParentCollateralAdjusted import *
from neuron import h, gui
import numpy as np
import random
import math
from glob import glob
from matplotlib import pyplot
from mpl_toolkits.mplot3d import axes3d, Axes3D
# from Stimulators_filter import *
from scipy.interpolate import griddata
#from findNodePositions import *
import sys
import pickle as pickle
import os

h.load_file('stdrun.hoc')

## These parameters define axon selections. Fibnum is a number 1-25, corresponding to the 25 3D axon collateral geometries
## gna_axnode is the sodium conductance (default: gna_axnode = 3.0)
## bouton_L = length of a node/bouton (default: bouton_L = 1.0)

fibnum = sys.argv[1]

gna_axnode = float(sys.argv[2])

bouton_L = float(sys.argv[3])


if len(sys.argv)  > 0:
	fibnum = sys.argv[1]
else:
	fibnum = '1'
os.chdir('./fiber_' + fibnum)

files = glob('*.csv')
axon_section_list = [None] * len(files)

bouton_list = [None] * len(files)

axon_section_list = [None] * len(files)

for file in files:
	with open(file) as filename:
		currentSection = file.strip('.csv').split('_')[0]

		pos_array = np.loadtxt(filename,delimiter = ',')

		axon_index = int(currentSection)

		first_sec_x = np.array(pos_array[0])
		first_sec_y = np.array(pos_array[1])
		first_sec_z = np.array(pos_array[2])


		if len(first_sec_x) > 50000:
			stepsize = 1000
		elif len(first_sec_x) > 20000:
			stepsize = 500
		elif len(first_sec_x) > 10000:
			stepsize = 100
		elif len(first_sec_x) > 500:
			stepsize = 50
		elif len(first_sec_x) > 200:
			stepsize = 5
		else:
			stepsize = 1


		x_list = first_sec_x[::stepsize]
		y_list = first_sec_y[::stepsize]
		z_list = first_sec_z[::stepsize]

		if (x_list[-1] != first_sec_x[-1]) or (y_list[-1] != first_sec_y[-1]):
			np.append(x_list,first_sec_x[-1])
			np.append(y_list,first_sec_y[-1])
			np.append(z_list,first_sec_z[-1])


		direction_vec = np.array([(x_list[-1] - x_list[-2]), (y_list[-1] - y_list[-2]), (z_list[-1] - z_list[-2])]) # get the direction of the last two points
		direction_vec = direction_vec / np.linalg.norm(direction_vec) ## turn it into a unit vector


		axon_section_list[axon_index] = h.Section(name = 'axon_sec_' + currentSection)

		bouton_list[axon_index] = h.Section(name = 'bouton_sec_' + currentSection)

		for i in range(len(x_list)):
			h.pt3dadd(x_list[i], y_list[i],z_list[i], 1.5, sec = axon_section_list[axon_index])


		h.pt3dadd(x_list[-1], y_list[-1], z_list[-1], bouton_diam, sec = bouton_list[axon_index])
		h.pt3dadd(x_list[-1] + bouton_L * direction_vec[0], y_list[-1] + bouton_L * direction_vec[1], z_list[-1] + bouton_L * direction_vec[2], bouton_diam, sec = bouton_list[axon_index])




for i in range(len(files)):
	bouton_list[i].connect(axon_section_list[i](1))


for file in files:

	currentSection = int(file.strip('.csv').split('_')[0])
	parentSection  = int(file.strip('.csv').split('_')[1])

	if parentSection >= 0:
		axon_section_list[currentSection].connect(bouton_list[parentSection](1))



axon_section_list[0].diam = 3.0
bouton_diam = 2.1
parent_list = [axon_section_list[0]]

## Adjust diam factor change how much the axon thins at branch points
diam_factor = 0.65
node_fact = 0.8

while(len(parent_list) > 0):
	parent = parent_list[0]
	parent_diam = parent.diam
	parent_ref = h.SectionRef(sec = parent)
	child_bout = parent_ref.child[0]
	child_bout.diam = node_fact * parent_diam
	child_bout_ref = h.SectionRef(sec = child_bout)
	children_branches = list(child_bout_ref.child)

	for child_branch in children_branches:
		if parent.diam * diam_factor > 0.5:
			child_branch.diam = parent.diam * diam_factor
		else:
			child_branch.diam = 0.5

	parent_list = parent_list + children_branches
	parent_list.pop(0)



new_axon_list = []

for sec in axon_section_list:
	axon_section_ref = h.SectionRef(sec = sec)
	child_bouton = axon_section_ref.child[0]
	axon_sec_diam = sec.diam

	if (sec.L / (sec.diam * 120.0)) > 1.:
		numSec = int(np.ceil(sec.L / (sec.diam * 120.0)))
		numpts3d = int(h.n3d(sec = sec))

		arc_length_cutoffs = [1.0/numSec * q for q in range(1,numSec)]
		cutoff_indices = []

		for cutoff_pct in arc_length_cutoffs:
			for pt_num in range(numpts3d):
				if ((h.arc3d(pt_num, sec = sec) / sec.L) >= cutoff_pct) and ((h.arc3d(pt_num - 1, sec = sec) / sec.L) < cutoff_pct):
					cutoff_indices.append(pt_num)

		cutoff_indices.append(numpts3d)

		pointsSplit3d = [None for secnum in range(numSec)]
		cutoff_start = 0
		for cutoff_pt_index in range(len(cutoff_indices)):
			cutoff_end = cutoff_indices[cutoff_pt_index]
			pointsSplit3d[cutoff_pt_index] = range(cutoff_start, cutoff_end)

			cutoff_start = cutoff_end - 1


		h.disconnect(sec = child_bouton)
		
		sub_axon_list = [sec]
		new_bouton_list = []
		for i in range(1,numSec):
			sub_axon_list.append(h.Section())
			new_bouton_list.append(h.Section())

		xpts    = [None] * numpts3d
		ypts    = [None] * numpts3d
		zpts    = [None] * numpts3d
		diampts = [None] * numpts3d
		for i in range(numpts3d):
			xpts[i]    = h.x3d(i, sec = sec)
			ypts[i]    = h.y3d(i, sec = sec)
			zpts[i]    = h.z3d(i, sec = sec)
			diampts[i] = h.diam3d(i, sec = sec)

		h.pt3dclear(sec = sec)

		for i in range(numSec):
			pts_section = pointsSplit3d[i]
			for j in pts_section:
				h.pt3dadd(xpts[j], ypts[j], zpts[j], diampts[j], sec = sub_axon_list[i])


			lastPtsDirection = [xpts[pts_section[-1]] - xpts[pts_section[-2]], ypts[pts_section[-1]] - ypts[pts_section[-2]], zpts[pts_section[-1]] - zpts[pts_section[-2]]]
			lastPtsDirectionLength = np.linalg.norm(lastPtsDirection)

			lastPtsDirection = [q/lastPtsDirectionLength for q in lastPtsDirection]

			if i < (numSec - 1):
				h.pt3dadd(0,0,0, node_fact * axon_sec_diam, sec= new_bouton_list[i])
				h.pt3dadd(lastPtsDirection[0] * bouton_L, lastPtsDirection[1] * bouton_L, lastPtsDirection[2] * bouton_L, node_fact * axon_sec_diam, sec = new_bouton_list[i])

				new_bouton_list[i].connect(sub_axon_list[i])

				sub_axon_list[i + 1].connect(new_bouton_list[i])

				bouton_list.append(new_bouton_list[i])
			else:
				child_bouton.connect(sub_axon_list[i])


		for j in range(1, numSec):
			new_axon_list.append(sub_axon_list[j])


axon_section_list = axon_section_list + new_axon_list








for i in range(21):
	bouton_list.insert(0,h.Section(name = 'newbout_' + str(i)))
for i in range(21):
	bouton_list[i].diam = 2.4
	bouton_list[i].L = 1.0
	

for j in range(20):
	axon_section_list.insert(0, h.Section(name = 'newaxon_' + str(j)))
for j in range(20):
	axon_section_list[j].L = 300.0
	axon_section_list[j].diam = 3.0





axon_section_list[20].connect(bouton_list[20])
bouton_list[20].connect(axon_section_list[19])
axon_section_list[19].connect(bouton_list[19])
bouton_list[19].connect(axon_section_list[18])
axon_section_list[18].connect(bouton_list[18])
bouton_list[18].connect(axon_section_list[17])
axon_section_list[17].connect(bouton_list[17])
bouton_list[17].connect(axon_section_list[16])
axon_section_list[16].connect(bouton_list[16])
bouton_list[16].connect(axon_section_list[15])
axon_section_list[15].connect(bouton_list[15])
bouton_list[15].connect(axon_section_list[14])
axon_section_list[14].connect(bouton_list[14])
bouton_list[14].connect(axon_section_list[13])
axon_section_list[13].connect(bouton_list[13])
bouton_list[13].connect(axon_section_list[12])
axon_section_list[12].connect(bouton_list[12])
bouton_list[12].connect(axon_section_list[11])
axon_section_list[11].connect(bouton_list[11])
bouton_list[11].connect(axon_section_list[10])
axon_section_list[10].connect(bouton_list[10])
bouton_list[10].connect(axon_section_list[9])
axon_section_list[9].connect(bouton_list[9])
bouton_list[9].connect(axon_section_list[8])
axon_section_list[8].connect(bouton_list[8])
bouton_list[8].connect(axon_section_list[7])
axon_section_list[7].connect(bouton_list[7])
bouton_list[7].connect(axon_section_list[6])
axon_section_list[6].connect(bouton_list[6])
bouton_list[6].connect(axon_section_list[5])
axon_section_list[5].connect(bouton_list[5])
bouton_list[5].connect(axon_section_list[4])
axon_section_list[4].connect(bouton_list[4])
bouton_list[4].connect(axon_section_list[3])

axon_section_list[3].connect(bouton_list[3])
bouton_list[3].connect(axon_section_list[2])

axon_section_list[2].connect(bouton_list[2])
bouton_list[2].connect(axon_section_list[1])

axon_section_list[1].connect(bouton_list[1])
bouton_list[1].connect(axon_section_list[0])

axon_section_list[0].connect(bouton_list[0])


h.define_shape()


os.chdir('..')

for sec in axon_section_list:
	if sec.L < sec.diam * 40.0:
		sec.L = sec.diam * 40.0





unmyelinatedRa = 100
unmyelinatedcm = 0.7
unmyelinated_sodium_density = 1500 # original 2500
unmyelinated_potassium_density = 1000 # original 500
ct = 0

g_ratio = 0.8
rhoa = 0.7e6
mycm = 0.1
mygm = 0.001
ct = 0
for sec in bouton_list:
	sec.insert('newaxnode') ### Node properties
	sec.el_newaxnode = -90.0 # use -85 for 1.5 conductivity, -90 for 3.0
	sec.gnapbar_newaxnode = .005 # 0.01 # original # 0.0 to fix firing
	sec.gnabar_newaxnode = 3.0 # 3.0 
	sec.nseg = 7
	sec.Ra = rhoa/10000.
	sec.cm = 2.0
	sec.insert('extracellular')
	nodeD = sec.diam
	space_p1 = .002
	Rpn0 = (rhoa*0.01)/( math.pi*((((nodeD/2.)+space_p1)**2.)-((nodeD/2.)**2.)))
	sec.xraxial[0] = Rpn0
	sec.xg[0] = 1e10
	sec.xc[0] = 0
	ct = ct+ 1
ct = 0



for sec in axon_section_list:
	sec.nseg = 21

	fiberD = sec.diam ## Diameter including myelin
	paraD2 = g_ratio * sec.diam
	axonD  = g_ratio * sec.diam
	space_p2 = 0.004 ## periaxonal space width, Maybe change to 0.002?
	nl = math.ceil(10.5 * axonD) ## Number of myelin laminae, from Cohen 2020 Fig 2C
	sec.insert('pas')
	sec.g_pas = 0.0001*paraD2/fiberD ## Perhaps .00001, like McIntyre 2004
	sec.e_pas = -80.0
	sec.cm = 2.*paraD2/fiberD
	sec.Ra = rhoa*(1/(paraD2/fiberD)**2)/10000.
	sec.insert('extracellular')
	# sec.insert('xtra') 
	Rpn2=(rhoa*.01)/(math.pi*((((paraD2/2.)+space_p2)**2.)-((paraD2/2.)**2.)))
	sec.xraxial[0] = Rpn2
	sec.xg[0] = mygm/(nl*2.)
	sec.xc[0] = mycm/(nl*2.)


ct = 0
for sec in bouton_list:
	sref = h.SectionRef(sec = sec)

	if sref.nchild() == 0:
		print(ct)
	ct =  ct + 1


ct = 0
terminal_list = []
for sec in bouton_list:
	sref = h.SectionRef(sec = sec)
	if sref.nchild() == 0:
		terminal_list.append(sec)
	ct = ct + 1


