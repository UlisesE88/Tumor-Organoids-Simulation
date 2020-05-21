
from pyMCDS import pyMCDS
import numpy as np
import pandas as pd

def generate_data(xml_file):
	mcds = pyMCDS(xml_file, '/Users/codeplus/Desktop/Tumor-Organoids-Simulation/PhysiCell/output')

	#Print the total volume of cells
	cell_volumes = mcds.data['discrete_cells']['total_volume']
	total_volume = 0
	for volumes in cell_volumes:
		total_volume += volumes
	print('')
	print('Total volume = ' + str(total_volume))
	print('')

	#Grab information just as position, volume, and the phase of each cell
	mcds.get_cell_variables()
	cell_df = mcds.get_cell_df()
	cell_table = cell_df[['position_x', 'position_y', 'total_volume', 'current_phase', 'migration_speed']]
	print(cell_table)

	oxy_dict = mcds.data['continuum_variables']['oxygen']
	oxy_dict['decay_rate']

	oxy_conc = oxy_dict['data']
	print(oxy_conc.shape)
	print(mcds.get_concentrations_at(x=0., y=0., z=0.))



	#Count the number of each cell type
	print('')
	cell_type = cell_df[['cell_type']]
	d = {}
	for index, row in cell_type.iterrows():
		type_cell = row['cell_type']

		if type_cell not in d:
			d[type_cell] = 0
		d[type_cell] += 1

	d_list = []
	for (type_cell, count) in d.items():
		if type_cell == 1.0:
			d_list.append(['stem_like', count])
		elif type_cell == 2.0:
			d_list.append(['epithelial', count])
		else:
			d_list.append(['mesenchymal', count])
  
	df = pd.DataFrame(d_list, columns = ['Cell Type', 'Count'])
	print(df.head())
	print(" ")

	#Count the number of each phase
	cell_phase = cell_df[['current_phase']]
	p = {}
	for index, row in cell_phase.iterrows():
		type_phase = row['current_phase']
		if type_phase not in p:
			p[type_phase] = 0
		p[type_phase] += 1

	p_list = []
	for (type_phase, count) in p.items():
		if type_phase == 1.0:
			p_list.append(['blue', count])
		elif type_phase == 2.0:
			p_list.append(['green', count])
		else:
			p_list.append(['red', count])
  
	pf = pd.DataFrame(p_list, columns = ['Cell Phase', 'Count'])
	print(pf.head())

if __name__ == '__main__':
	lst_of_files = ['output00000000.xml', 'final.xml']
	
	for xml_file in lst_of_files:
		generate_data(xml_file)
		print('')

		#'output00000001.xml','output00000002.xml',
					#'output00000003.xml','output00000004.xml','output00000005.xml',
					#'output00000006.xml',

'''
ID: unique integer that will identify the cell throughout its lifetime in the simulation
position(_x, _y, _z): floating point positions for the cell in 𝑥, 𝑦, and 𝑧 directions
total_volume: total volume of the cell
cell_type: integer label for the cell as used in PhysiCell
cycle_model: integer label for the cell cycle model as used in PhysiCell
current_phase: integer specification for which phase of the cycle model the cell is currently in
elapsed_time_in_phase: time that cell has been in current phase of cell cycle model
nuclear_volume: volume of cell nucleus
cytoplasmic_volume: volume of cell cytoplasm
fluid_fraction: proportion of the volume due to fliud
calcified_fraction: proportion of volume consisting of calcified material
orientation(_x, _y, _z): direction in which cell is pointing
polarity:
migration_speed: current speed of cell
motility_vector(_x, _y, _z): current direction of movement of cell
migration_bias: coefficient for stochastic movement (higher is “more deterministic”)
motility_bias_direction(_x, _y, _z): direction of movement bias
persistence_time: time in-between direction changes for cell
motility_reserved:
'''
