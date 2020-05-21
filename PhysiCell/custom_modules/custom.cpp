/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"

// declare cell definitions here 

//Cell_Definition civilian; 
//Cell_Definition Thanos; 
//Cell_Definition avenger;

Cell_Definition epithelial; 
Cell_Definition mesenchymal;
Cell_Definition stem_like;

Cycle_Model EMT_model;

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	
	SeedRandom(parameters.ints("random_seed") ); // or specify a seed here 
	
	//Set up the EMT Model
	EMT_model.code = PhysiCell_constants::custom_cycle_model; 
	EMT_model.name = "Red Green Blue";
	// add the phases
	int red_index = EMT_model.add_phase( 0 , "Red" );
	int green_index = EMT_model.add_phase( 1 , "Green" ); 
	int blue_index = EMT_model.add_phase( 2 , "Blue" );
	// set the Blue phase to have division at its end 
	EMT_model.phases[blue_index].division_at_phase_exit = true;
	// set up the phase links
	EMT_model.add_phase_link( red_index , green_index , NULL ); // red -> green 
	EMT_model.add_phase_link( green_index , blue_index , NULL ); // green -> blue 
	EMT_model.add_phase_link( blue_index , red_index , NULL ); // blue -> red
	// set the transition rates 
	EMT_model.transition_rate(red_index,green_index) = 1.0/( 60.0 * 5.0 );
	// mean duration: 5 hours 
	EMT_model.transition_rate(green_index,blue_index) = 1.0/( 60.0 * 8.0 );
	// mean duration: 8 hours
	EMT_model.transition_rate(blue_index,red_index) = 1.0/( 60.0 * 2.0 ); // mean duration: 2 hours
	// display the model 
	EMT_model.display( std::cout );




	// housekeeping 
	
	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	// Name the default cell type 
	
	cell_defaults.type = 0; 
	cell_defaults.name = "tumor cell"; 
	
	// set default cell cycle model 

	cell_defaults.functions.cycle_model = live; 
	
	// set default_cell_functions; 
	
	cell_defaults.functions.update_phenotype = update_cell_and_death_parameters_O2_based;
	//NULL; 
	
	// needed for a 2-D simulation: 
	
	/* grab code from heterogeneity */ 
	
	cell_defaults.functions.set_orientation = up_orientation; 
	cell_defaults.phenotype.geometry.polarity = 1.0;
	cell_defaults.phenotype.motility.restrict_to_2D = true; 
	
	// make sure the defaults are self-consistent. 
	
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.sync_to_functions( cell_defaults.functions ); 

	// set the rate terms in the default phenotype 

	// first find index for a few key variables. 
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );
	int oxygen_substrate_index = microenvironment.find_density_index( "oxygen" ); 

	int G0G1_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::G0G1_phase );
	int S_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::S_phase );

	// initially no necrosis 
	cell_defaults.phenotype.death.rates[necrosis_model_index] = 0.0; 

	// add custom data here, if any 
	
	//cell_defaults.custom_data.add_variable( "health" , "dimensionless", 1.0 );
	
	
	
	// Now, let's define another cell type. 
	// It's best to just copy the default and modify it. 
	
	// make this cell type randomly motile, less adhesive, greater survival, 
	// and less proliferative 
	
	stem_like = cell_defaults; 
	stem_like.type = 1; 
	stem_like.name = "stem_like"; 
	
	// make sure the new cell type has its own reference phenotype
	
	stem_like.parameters.pReference_live_phenotype = &( stem_like.phenotype ); 
	
	// disable motility 
	stem_like.phenotype.motility.is_motile = false; 
	
	// set birth rate
	stem_like.phenotype.cycle.data.transition_rate(0,0) = 
		parameters.doubles( "stem_like_birth_rate" ); // 0.0007; 
	
	// set death rate 
	stem_like.phenotype.death.rates[apoptosis_model_index] = 
		parameters.doubles( "stem_like_death_rate" ); // 0.000007; 
	
	// Set cell-cell adhesion to 5% of other cells 
	stem_like.phenotype.mechanics.cell_cell_adhesion_strength = 2.0;
	stem_like.phenotype.mechanics.cell_cell_repulsion_strength = 10.0;
	
	
	//stem_like.custom_data["strength"] = parameters.doubles("stem_like_strength" ); 
	
	// Thanos setup 
	
	epithelial = cell_defaults; 
	epithelial.type = 2; 
	epithelial.name = "epithelial"; 
	
	epithelial.functions.cycle_model = live; 
	epithelial.parameters.pReference_live_phenotype = &( epithelial.phenotype ); 
		

	// Set cell-cell adhesion to 5% of other cells 
	epithelial.phenotype.mechanics.cell_cell_adhesion_strength *= parameters.doubles( "epithelial_relative_adhesion" );
 	epithelial.phenotype.mechanics.cell_cell_repulsion_strength *= parameters.doubles( "epithelial_relative_repulsion" );

 	//epithelial.phenotype.mechanics.cell_cell_adhesion_strength = 1;
 	//epithelial.phenotype.mechanics.cell_cell_repulsion_strength = 90.0;


	epithelial.phenotype.cycle.data.transition_rate(0,0) = parameters.doubles( "epithelial_birth_rate" ); 
	
	epithelial.phenotype.motility.is_motile = true; 
	epithelial.phenotype.motility.migration_speed = parameters.doubles("epithelial_speed"); 
	epithelial.phenotype.motility.migration_bias = parameters.doubles("epithelial_motility_bias");  
	epithelial.phenotype.motility.persistence_time = 1;
	
	epithelial.functions.update_migration_bias = epithelial_taxis; 
	
	//epithelial.custom_data["strength"] = parameters.doubles("epithelial_strength" ); 
	
	int epithelial_sig_i = microenvironment.find_density_index( "epithelial" );
	epithelial.phenotype.secretion.secretion_rates[epithelial_sig_i] = 10; 
	epithelial.phenotype.secretion.saturation_densities[epithelial_sig_i] = 1; 
	
	epithelial.phenotype.geometry.radius = 7.065;

	epithelial.phenotype.death.rates[apoptosis_model_index] = parameters.doubles( "epithelial_apoptosis_rate" );

	epithelial.functions.custom_cell_rule = epithelial_cell_rule; 

	//epithelial.phenotype.cycle.sync_to_cycle_model( EMT_model );


	// Avenger setup 
	


















	mesenchymal = epithelial; 
	mesenchymal.type = 3; 
	mesenchymal.name = "mesenchymal";

	// mesenchymal.phenotype.mechanics.cell_cell_repulsion_strength = 20.0;
	// mesenchymal.phenotype.mechanics.cell_cell_adhesion_strength = 0;


	// mesenchymal.functions.cycle_model = flow_cytometry_separated_cycle_model; 

	// int M_phase_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::M_phase );
	// int G0G1_phase_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::G0G1_phase );


	// mesenchymal.phenotype.cycle.data.transition_rate(G0G1_phase_index,G0G1_phase_index) = 
	// 	parameters.doubles( "mesenchymal_birth_rate" );
	// mesenchymal.phenotype.death.rates[apoptosis_model_index] = 
	// 	parameters.doubles( "mesenchymal_death_rate" );


	// mesenchymal.phenotype.motility.migration_speed = parameters.doubles("mesenchymal_speed"); 
	// mesenchymal.phenotype.motility.migration_bias = parameters.doubles("mesenchymal_motility_bias"); 
	// mesenchymal.phenotype.motility.persistence_time = 1; 
	
	// mesenchymal.functions.update_migration_bias = mesenchymal_taxis; 	

	// //mesenchymal.custom_data["strength"] = parameters.doubles( "mesenchymal_strength" ); 


	// int mesenchymal_sig_i = microenvironment.find_density_index( "mesenchymal" );
	// mesenchymal.phenotype.secretion.secretion_rates[epithelial_sig_i] = 0; 
	// mesenchymal.phenotype.secretion.saturation_densities[epithelial_sig_i] = 1; 

	// mesenchymal.phenotype.secretion.secretion_rates[mesenchymal_sig_i] = 10; 
	// mesenchymal.phenotype.secretion.saturation_densities[mesenchymal_sig_i] = 1; 

	// Thanos is a big guy
	stem_like.phenotype.volume.multiply_by_ratio( 0.1 ); //μm3
	//epithelial.phenotype.volume.multiply_by_ratio( 0.1 ); 
	//mesenchymal.phenotype.volume.multiply_by_ratio( 0.1 ); 
	
	return; 
}


void epithelial_cell_rule( Cell* pCell , Phenotype& phenotype , double dt )
{
	
	std::vector<Cell*> nearby = pCell->cells_in_my_container(); 
	
	// if at least 2 neighbors, check therehold prob to convert into mesenchymal cell	
	update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);

	std::random_device rd; 
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<> dis(0.0, 1.0);
    float prob = dis(gen);

	if( nearby.size() > 5 && prob > 0.5)
	{
		//delete_cell(pCell -> index);

		//pCell->convert_to_cell_definition(mesenchymal);
		pCell->	phenotype.cycle.sync_to_cycle_model( EMT_model );

		//phenotype.motility.is_motile = false;
		//pCell->functions.update_phenotype = NULL; 
	}
	
	return; 
}


void setup_microenvironment( void )
{
	// set domain parameters 
	
/* now this is in XML 
	default_microenvironment_options.X_range = {-1000, 1000}; 
	default_microenvironment_options.Y_range = {-1000, 1000}; 
	default_microenvironment_options.simulate_2D = true; 
*/
	
	// make sure to override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl; 
		default_microenvironment_options.simulate_2D = true; 
	}
	
/* now this is in XML 	
	// no gradients need for this example 

	default_microenvironment_options.calculate_gradients = false; 
	
	// set Dirichlet conditions 

	default_microenvironment_options.outer_Dirichlet_conditions = true;
	
	// if there are more substrates, resize accordingly 
	std::vector<double> bc_vector( 1 , 38.0 ); // 5% o2
	default_microenvironment_options.Dirichlet_condition_vector = bc_vector;
	
	// set initial conditions 
	default_microenvironment_options.initial_condition_vector = { 38.0 }; 
*/
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	// create some cells near the origin
	
	Cell* pC;
	
	// get domain size 
	int X_nodes = microenvironment.mesh.x_coordinates.size(); 
	double X_left  = microenvironment.mesh.x_coordinates[0] - microenvironment.mesh.dx / 2.0 ; 
	double X_right = microenvironment.mesh.x_coordinates[X_nodes-1] + microenvironment.mesh.dx / 2.0 ; 
	double X_length = X_right - X_left; 
	
	int Y_nodes = microenvironment.mesh.y_coordinates.size(); 
	double Y_left  = microenvironment.mesh.y_coordinates[0] - microenvironment.mesh.dy / 2.0 ; 
	double Y_right = microenvironment.mesh.y_coordinates[Y_nodes-1] + microenvironment.mesh.dy / 2.0 ; 
	double Y_length = Y_right - Y_left; 
	
	int number_of_stem_like = parameters.ints( "stem_like_initial_count" ); 
	for( int n = 0; n < number_of_stem_like ; n++ )
	{
		std::vector<double> position = {0,0,0}; 

		double x = X_left + 0.1*X_length + 0.8*X_length*UniformRandom(); 
		double y = Y_left + 0.1*Y_length + 0.8*Y_length*UniformRandom(); 
		position[0] = x;
		position[1] = y; 
		
		pC = create_cell( stem_like ); 
		pC->assign_position( position );
		
	}
	
	// seed avengers 
	
	int number_of_mesenchymal = parameters.ints( "mesenchymal_initial_count" ); 
	for( int n = 0; n < number_of_mesenchymal ; n++ )
	{
		std::vector<double> position = {0,0,0}; 

		double x = X_left + 0.1*X_length + 0.8*X_length*UniformRandom(); 
		double y = Y_left + 0.1*Y_length + 0.8*Y_length*UniformRandom(); 
		position[0] = x;
		position[1] = y; 
		
		pC = create_cell( mesenchymal ); 
		pC->assign_position( position );
	}

	int number_of_epithelial = parameters.ints( "epithelial_initial_count" ); 
	for( int n = 0; n < number_of_epithelial ; n++ )
	{
		std::vector<double> position = {0,0,0}; 

		double x = X_left + 0.1*X_length + 0.8*X_length*UniformRandom(); 
		double y = Y_left + 0.1*Y_length + 0.8*Y_length*UniformRandom(); 
		position[0] = x;
		position[1] = y; 
		
		pC = create_cell( epithelial ); 
		pC->assign_position( position );
	}




	int num_cells = (*all_cells).size();
	double dt = 0.1; // should query global mechanics time step size
	double startup_time = 5; // 5 minutes enough?
	int startup_steps = (int) ceil( startup_time / dt );
 
	for( int n=0 ; n < startup_steps ; n++ )
	{
  		// update velocities   
  		#pragma omp parallel for
  		for( int i = 0; i < num_cells ; i++ )
 		{
   			Cell* pC = (*all_cells)[i];
   			pC->functions.update_velocity( pC, pC->phenotype , dt );
  		}
 
  		// update positions
  		#pragma omp parallel for
  		for( int i = 0; i < num_cells ; i++ )
  		{
   			(*all_cells)[i]->update_position( dt );
  		}
	}
	return;
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	// start off with black 
	std::vector<std::string> output( 4 , "black" ); 
	
	// stay black if dead 
	if( pCell->phenotype.death.dead == true )
	{ return output; }

	//printf("%d\n", pCell->phenotype.cycle.current_phase_index());


	if(pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::live ){
		std::string color = "mediumpurple"; 

		output[0] = "lavender";  
		output[2] = color; 
		output[3] = color; 
		return output; 
	}


	
	if( pCell->phenotype.cycle.current_phase_index() == 0)
	{
		std::string color = "indigo"; 

		output[0] = "lavender";
		output[2] = color; 
		output[3] = color; 
		return output; 
	}

	if( pCell->phenotype.cycle.current_phase_index() == 1)
	{
		std::string color = "green";

		output[0] = color;
		output[2] = color; 
		output[3] = color; 
		return output; 
	}

	if( pCell->phenotype.cycle.current_phase_index() == 2)
	{
		std::string color = "blue";

		output[0] = color;
		output[2] = color; 
		output[3] = color; 
		return output; 
	}



		
	// stem_like cells are green 
	if( pCell->type == stem_like.type )
	{
		std::string color = "limegreen"; 

		output[0] = "honeydew"; 
		output[2] = color; 
		output[3] = color; 
		return output; 
	}
	
	// epithelial is purple 
	// if( pCell->type == epithelial.type )
	// {
	// 	std::string color = "mediumpurple"; 

	// 	output[0] = "lavender";  
	// 	output[2] = color; 
	// 	output[3] = color; 
	// 	return output; 
	// }	
	
	// mesenchymal are red 
	if( pCell->type == mesenchymal.type )
	{
		std::string color = "red"; 

		output[0] = "mistyrose";  
		output[2] = color; 
		output[3] = color; 
		return output; 
	}	
	return output; 
}

void mesenchymal_taxis( Cell* pCell , Phenotype& phenotype , double dt )
{
	static int mesenchymal_i = microenvironment.find_density_index( "mesenchymal" );

	phenotype.motility.migration_bias_direction = pCell->nearest_gradient(mesenchymal_i); 
	normalize( &(phenotype.motility.migration_bias_direction) ); 

	return; 
}

void epithelial_taxis( Cell* pCell , Phenotype& phenotype , double dt )
{
	static int epithelial_i = microenvironment.find_density_index( "epithelial" );

	phenotype.motility.migration_bias_direction = pCell->nearest_gradient(epithelial_i); 
	normalize( &(phenotype.motility.migration_bias_direction) ); 

	return; 
}
















// keep 
void add_elastic_velocity( Cell* pActingOn, Cell* pAttachedTo , double elastic_constant )
{
	std::vector<double> displacement = pAttachedTo->position - pActingOn->position; 
	
	// dettach cells if too far apart 
	static double max_elastic_displacement = parameters.doubles("max_elastic_displacement");
	static double max_displacement_squared = max_elastic_displacement*max_elastic_displacement; 
	
	// if( norm_squared( displacement ) > max_displacement_squared )
	// {
	// 	dettach_cells( pActingOn , pAttachedTo );
	// 	std::cout << "\t\tDETACH!!!!!" << std::endl; 
	// 	return; 
	// }
	
	axpy( &(pActingOn->velocity) , elastic_constant , displacement ); 
	
	return; 
}

void extra_elastic_attachment_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
	// for( int i=0; i < pCell->state.neighbors.size() ; i++ )
	// {
	// 	add_elastic_velocity( pCell, pCell->state.neighbors[i], pCell->custom_data["elastic coefficient"] ); 
	// }

	return; 
}	

// keep 
void attach_cells( Cell* pCell_1, Cell* pCell_2 )
{
	#pragma omp critical
	{
	
	//Iterate through pCell_1 neighbors and check if pCell_2 is already attached
	bool already_attached = false; 
	for( int i=0 ; i < pCell_1->state.neighbors.size() ; i++ )
	{
		if( pCell_1->state.neighbors[i] == pCell_2 )
		{ already_attached = true; }
	}

	//If pCell_2 is not already attached to pCell_1, then add it to its neighbors vector
	if( already_attached == false )
	{ 
		pCell_1->state.neighbors.push_back( pCell_2 ); 
	}
	
	//Iterate through pCell_2 neighbors and check if pCell_1 is already attached
	already_attached = false; 
	for( int i=0 ; i < pCell_2->state.neighbors.size() ; i++ )
	{
		if( pCell_2->state.neighbors[i] == pCell_1 )
		{ already_attached = true; }
	}

	//If pCell_1 is not already attached to pCell_2, then add it to its neighbors vector
	if( already_attached == false )
	{ 
		pCell_2->state.neighbors.push_back( pCell_1 ); 
	}

	}
	return; 
}

// keep 
void dettach_cells( Cell* pCell_1 , Cell* pCell_2 )
{
	#pragma omp critical
	{
		bool found = false; 
		int i = 0; 
		while( !found && i < pCell_1->state.neighbors.size() )
		{
			// if cell 2 is in cell 1's list, remove it
			if( pCell_1->state.neighbors[i] == pCell_2 )
			{
				int n = pCell_1->state.neighbors.size(); 
				// copy last entry to current position 
				pCell_1->state.neighbors[i] = pCell_1->state.neighbors[n-1]; 
				// shrink by one 
				pCell_1->state.neighbors.pop_back(); 
				found = true; 
			}
			i++; 
		}
	
		found = false; 
		i = 0; 
		while( !found && i < pCell_2->state.neighbors.size() )
		{
			// if cell 1 is in cell 2's list, remove it
			if( pCell_2->state.neighbors[i] == pCell_1 )
			{
				int n = pCell_2->state.neighbors.size(); 
				// copy last entry to current position 
				pCell_2->state.neighbors[i] = pCell_2->state.neighbors[n-1]; 
				// shrink by one 
				pCell_2->state.neighbors.pop_back(); 
				found = true; 
			}
			i++; 
		}

	}
	
	return; 
}

// keep! 
Cell* worker_cell_check_neighbors_for_attachment( Cell* pWorker , double dt )
{
	std::vector<Cell*> nearby = pWorker->cells_in_my_container(); 
	int i = 0; 
	while( i < nearby.size() )
	{
		// don't try to attach to yourself 
		if( nearby[i] != pWorker )
		{
			if( worker_cell_attempt_attachment( pWorker, nearby[i] , dt ) )
			{ return nearby[i]; }
		}
		i++; 
	}
	
	return NULL; 
}

// keep! 
bool worker_cell_attempt_attachment( Cell* pWorker, Cell* pCargo , double dt )
{
	static int receptor_i = pCargo->custom_data.find_variable_index( "receptor" ); 

	static double receptor_threshold = 
		parameters.doubles("attachment_receptor_threshold"); // 0.1; 
	
	static double max_attachment_distance = 
		parameters.doubles("max_attachment_distance"); // 18.0; 
	static double min_attachment_distance = 
		parameters.doubles("min_attachment_distance"); // 14.0; 
	static double attachment_difference = max_attachment_distance - min_attachment_distance; 
	
	if( pCargo->custom_data[receptor_i] > receptor_threshold )
	{
		std::vector<double> displacement = pCargo->position - pWorker->position;
		double distance = norm( displacement ); 
		if( distance > max_attachment_distance )
		{ return false; } 
	
		if( distance < min_attachment_distance )
		{ 
			attach_cells( pWorker, pCargo );
			return true; 
		}
		
		return true; 
	}
	
	return false; 
}

void worker_cell_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
	// if I am dead, don't bother

	if( phenotype.death.dead == true )
	{
		// the cell death functions don't automatically turn off custom functions, 
		// since those are part of mechanics. 
		
		// Let's just fully disable now. 
		pCell->functions.custom_cell_rule = NULL; 
		return; 
	}
	
	// am I searching for cargo? if so, see if I've found it
	if( pCell->state.neighbors.size() == 0 )
	{
		std::vector<Cell*> nearby = pCell->cells_in_my_container(); 
		
		//bool attached = false; // want to limit to one attachment 
		for( int i=0; i < nearby.size(); i++ )
		{
			// if it is expressing the receptor, dock with it 
			if( nearby[i]->custom_data["receptor"] > 0.5 )
			{
				attach_cells( pCell, nearby[i] ); 
				// nearby[i]->custom_data["receptor"] = 0.0; // put into cargo cell rule instead? 
				// nearby[i]->phenotype.secretion.set_all_secretion_to_zero(); // put into cargo rule instead? 
				//attached = true; 
			}
		}
		
	}
	
	return; 
}
