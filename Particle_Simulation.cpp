// Particle_Simulation.cpp : Defines the entry point for the application.
//

/*
	NOTES:

	Ensure particle interaction range does not extend further than one neighbor in the cell array.
	Periodic boundary conditions are assumed for now

	Still need to verify correctness of 3D neighbor list
*/
#pragma once

#include <chrono>
#include <random>
#include <cmath>
#include <fstream>
#include <stdexcept>

//#include <emmintrin.h>
#include <immintrin.h>

#include "Particle_Simulation.h"

using namespace std;

Cell_Structure Cell_Grid[2];

Particle_Structure* Particles;

// simulation variables
double Particle_Diameter = 1;
double Colloid_Diameter = 5; // nominal 5

double Reaction_Radius = 7.5; // nominal 7
double Change_Rate = 0.0003;

double End_Time = 1500.0;
double Delta_t = 0.00003;
double sigma = 1;
double Epsilon_WCA = 10;
double Epsilon_LJ = 3;

int Dimension = 2;
int Number_Particles = 750;
int Number_Colloids = 128;

double Region_Size[3];

double Sigma_Square;
double Sigma_Square_Colloid_Colloid;
double Sigma_Square_Colloid_Particle;
double Sigma_Test_WCA;
double Sigma_Test_WCA_Colloid_Particle;
double Sigma_Test_WCA_Colloid_Colloid;
double Sigma_Test_LJ;
double Epsilon_WCA_x_24;
double Epsilon_LJ_x_24;
double Reaction_Radius_Squared;

int Neighbor_Count;

int main(int argc, char* argv[])
{
	Sigma_Square = pow(Particle_Diameter, 2);
	Sigma_Square_Colloid_Colloid = pow(Colloid_Diameter, 2);
	Sigma_Square_Colloid_Particle = pow((Colloid_Diameter + Particle_Diameter) / 2, 2);
	
	Sigma_Test_WCA = pow(2.0, (2 / 6.0)) * Sigma_Square;
	Sigma_Test_LJ = pow(2.5, 2) * Sigma_Square;

	Sigma_Test_WCA_Colloid_Particle = pow(2.0, (2 / 6.0)) * pow((Colloid_Diameter + Particle_Diameter)/2, 2);
	Sigma_Test_WCA_Colloid_Colloid = pow(2.0, (2 / 6.0)) * pow(Colloid_Diameter, 2);;

	Epsilon_WCA_x_24 = Epsilon_WCA * 24.0;
	Epsilon_LJ_x_24 = Epsilon_LJ * 24.0;

	Reaction_Radius_Squared = Reaction_Radius * Reaction_Radius;
	
	Region_Size[0] = 100;
	Region_Size[1] = 100;

	Cell_Grid[0].Divisions[0] = 20;
	Cell_Grid[0].Divisions[1] = 20;

	Cell_Grid[1].Divisions[0] = 10;
	Cell_Grid[1].Divisions[1] = 10;

	if (Dimension == 2)
	{
		Region_Size[2] = 0;

		Cell_Grid[0].Divisions[2] = 0;
		Cell_Grid[0].Total_Divisions = Cell_Grid[0].Divisions[0] * Cell_Grid[0].Divisions[1];

		Cell_Grid[1].Divisions[2] = 0;
		Cell_Grid[1].Total_Divisions = Cell_Grid[1].Divisions[0] * Cell_Grid[1].Divisions[1];

		Neighbor_Count = 4;
	}
	else
	{
		Region_Size[2] = 100;

		Cell_Grid[0].Divisions[2] = 20;
		Cell_Grid[0].Total_Divisions = Cell_Grid[0].Divisions[0] * Cell_Grid[0].Divisions[1] * Cell_Grid[0].Divisions[2];

		Cell_Grid[1].Divisions[2] = 10;
		Cell_Grid[1].Total_Divisions = Cell_Grid[1].Divisions[0] * Cell_Grid[1].Divisions[1] * Cell_Grid[1].Divisions[2];

		Neighbor_Count = 13;
	}

	Build_Memory();

	cout << "memory allocated" << endl;

	Assign_Particle_Data();

	Initialize_Positions();

	cout << "positions assigned" << endl;

	Build_Neighbor_List(0);
	Build_Neighbor_List(1);

	Build_Neighbor_Offset_List(0);
	Build_Neighbor_Offset_List(1);

	double Current_Time = 0;

	auto start = std::chrono::system_clock::now();

	std::random_device Get_Seed;
	std::mt19937 Random_Normal(Get_Seed());
	std::normal_distribution<double> Random_Force(0, 1);

	std::mt19937 Random_Uniform(Get_Seed());
	std::uniform_real_distribution<double> Random_Reaction(0, 1);

	double chance = 0;

	int steps = 0;

	string iter_val = "";

	if (argc == 2)
	{
		iter_val = argv[1];
	}

	string Position_File_Name = "output/Position" + iter_val + ".txt";
	string MSD_File_Name = "output/MSD" + iter_val + ".txt";

	ofstream Output_File;
	Output_File.open(Position_File_Name, ios::out | ios::trunc);

	ofstream MSD_File;
	MSD_File.open(MSD_File_Name, ios::out | ios::trunc);

	double delta_x = 0;
	double delta_y = 0;
	double delta_z = 0;

	cout << "ini complete" << endl;

	if (Dimension == 2)
	{
		while (Current_Time < End_Time)
		{
			Build_Cell_List_2D(0);
			Build_Cell_List_2D(1);

			for (int i = 0; i < Number_Particles; i++)
			{
				Particles[i].Force_X = 0;
				Particles[i].Force_Y = 0;

				// do reaction
				if (Particles[i].Inside_Reaction)
				{
					if (Particles[i].Type == 0)
					{
						if (Random_Reaction(Random_Uniform) < Change_Rate)
						{
							Particles[i].Type = 1;
						}
					}
				}
				else
				{
					if (Particles[i].Type == 1)
					{
						if (Random_Reaction(Random_Uniform) < Change_Rate)
						{
							Particles[i].Type = 0;
						}
					}
				}

				//Particles[i].Type = 0;

				Particles[i].Inside_Reaction = false;
			}

			for (int i = Number_Particles; i < Number_Particles + Number_Colloids; i++)
			{
				Particles[i].Force_X = 0;
				Particles[i].Force_Y = 0;
			}

			Calculate_Forces_2D();

			for (int i = 0; i < Number_Particles + Number_Colloids; i++)
			{				
				double diff_delta_t = Particles[i].Diffusion * Delta_t;

				double rand_num_1 = Random_Force(Random_Normal);
				double rand_num_2 = Random_Force(Random_Normal);

				// normal
				
				double srqt_diff_delta_t = sqrt(2 * diff_delta_t);
				
				delta_x = (Particles[i].Force_X * diff_delta_t + srqt_diff_delta_t * rand_num_1);
				delta_y = (Particles[i].Force_Y * diff_delta_t + srqt_diff_delta_t * rand_num_2);

				Particles[i].Delta_X += delta_x;
				Particles[i].Delta_Y += delta_y;
				
				Particles[i].Position_X += delta_x;
				Particles[i].Position_Y += delta_y;

				while (Particles[i].Position_X >= Region_Size[0])
				{
					Particles[i].Position_X -= Region_Size[0];
					/*
					if (Particles[i].Position_X >= Region_Size[0])
					{
						cout << "bad X + " << i << endl;
					}
					*/
				}
				while (Particles[i].Position_Y >= Region_Size[1])
				{
					Particles[i].Position_Y -= Region_Size[1];
					/*
					if (Particles[i].Position_Y >= Region_Size[1])
					{
						cout << "bad Y + " << i << endl;
					}
					*/
				}

				while (Particles[i].Position_X < 0)
				{
					Particles[i].Position_X += Region_Size[0];
					/*
					if (Particles[i].Position_X < 0)
					{
						cout << "bad X - " << i << endl;
					}
					*/
				}
				while (Particles[i].Position_Y < 0)
				{
					Particles[i].Position_Y += Region_Size[1];
					/*
					if (Particles[i].Position_Y < 0)
					{
						cout << "bad Y - " << i << endl;
					}
					*/
				}
			}

			Current_Time += Delta_t;

			steps += 1;

			if ((steps % 100000) == 0)
			{
				for (int i = 0; i < Number_Particles + Number_Colloids; i++)
				{
					Output_File << Particles[i].Type << "\t" << Particles[i].Position_X << "\t" << Particles[i].Position_Y << "\t" << Particles[i].Position_Z << "\n";
				}

				for (int i = Number_Particles; i < Number_Particles + Number_Colloids; i++)
				{
					MSD_File << Particles[i].Delta_X << "\t" << Particles[i].Delta_Y << "\n";
				}
			}
		}
	}
	else
	{
		while (Current_Time < End_Time)
		{
			Build_Cell_List_3D(0);
			Build_Cell_List_3D(1);

			for (int i = 0; i < Number_Particles; i++)
			{
				Particles[i].Force_X = 0;
				Particles[i].Force_Y = 0;
				Particles[i].Force_Z = 0;

				// do reaction
				if (Particles[i].Inside_Reaction)
				{
					if (Particles[i].Type == 0)
					{
						if (Random_Reaction(Random_Uniform) < Change_Rate)
						{
							Particles[i].Type = 1;
						}
					}
				}
				else
				{
					if (Particles[i].Type == 1)
					{
						if (Random_Reaction(Random_Uniform) < Change_Rate)
						{
							Particles[i].Type = 0;
						}
					}
				}

				Particles[i].Inside_Reaction = false;
			}

			for (int i = Number_Particles; i < Number_Particles + Number_Colloids; i++)
			{
				Particles[i].Force_X = 0;
				Particles[i].Force_Y = 0;
				Particles[i].Force_Z = 0;
			}

			Calculate_Forces_3D();

			for (int i = 0; i < Number_Particles + Number_Colloids; i++)
			{
				delta_x = (Particles[i].Force_X * Particles[i].Diffusion * Delta_t + sqrt(2 * Particles[i].Diffusion * Delta_t) * Random_Force(Random_Normal));
				delta_y = (Particles[i].Force_Y * Particles[i].Diffusion * Delta_t + sqrt(2 * Particles[i].Diffusion * Delta_t) * Random_Force(Random_Normal));
				delta_z = (Particles[i].Force_Z * Particles[i].Diffusion * Delta_t + sqrt(2 * Particles[i].Diffusion * Delta_t) * Random_Force(Random_Normal));

				Particles[i].Delta_X += delta_x;
				Particles[i].Delta_Y += delta_y;
				Particles[i].Delta_Z += delta_z;

				Particles[i].Position_X += delta_x;
				Particles[i].Position_Y += delta_y;
				Particles[i].Position_Z += delta_z;

				while (Particles[i].Position_X >= Region_Size[0])
				{
					Particles[i].Position_X -= Region_Size[0];
					/*
					if (Particles[i].Position_X >= Region_Size[0])
					{
						cout << "bad X + " << i << endl;
					}
					*/
				}
				while (Particles[i].Position_Y >= Region_Size[1])
				{
					Particles[i].Position_Y -= Region_Size[1];
					/*
					if (Particles[i].Position_Y >= Region_Size[1])
					{
						cout << "bad Y + " << i << endl;
					}
					*/
				}
				while (Particles[i].Position_Z >= Region_Size[2])
				{
					Particles[i].Position_Z -= Region_Size[2];
					/*
					if (Particles[i].Position_Z >= Region_Size[2])
					{
						cout << "bad Z + " << i << endl;
					}
					*/
				}

				while (Particles[i].Position_X < 0)
				{
					Particles[i].Position_X += Region_Size[0];
					/*
					if (Particles[i].Position_X < 0)
					{
						cout << "bad X - " << i << endl;
					}
					*/
				}
				while (Particles[i].Position_Y < 0)
				{
					Particles[i].Position_Y += Region_Size[1];
					/*
					if (Particles[i].Position_Y < 0)
					{
						cout << "bad Y - " << i << endl;
					}
					*/
				}
				while (Particles[i].Position_Z < 0)
				{
					Particles[i].Position_Z += Region_Size[2];
					/*
					if (Particles[i].Position_Z < 0)
					{
						cout << "bad Z - " << i << endl;
					}
					*/
				}
			}

			Current_Time += Delta_t;

			steps += 1;

			if ((steps % 333333) == 0)
			{
				for (int i = 0; i < Number_Particles + Number_Colloids; i++)
				{
					Output_File << Particles[i].Type << "\t" << Particles[i].Position_X << "\t" << Particles[i].Position_Y << "\t" << Particles[i].Position_Z << "\n";
				}

				for (int i = Number_Particles; i < Number_Particles + Number_Colloids; i++)
				{
					MSD_File << Particles[i].Delta_X << "\t" << Particles[i].Delta_Y << "\t" << Particles[i].Delta_Z << "\n";
				}				
			}
		}
	}

	Output_File.close();
	MSD_File.close();

	auto end = std::chrono::system_clock::now();
	auto elapsed = end - start;
	std::cout << elapsed.count() << '\n';

	Clear_Memory();

	int a;

	//cin >> a;

	return 0;
}

void Build_Cell_List_2D(int index)
{
	double X_Divisor, Y_Divisor;

	int Local_Cell_X, Local_Cell_Y;

	int Final_Cell;

	int Div_X = Cell_Grid[index].Divisions[0];
	int Div_Y = Cell_Grid[index].Divisions[1];

	for (int i = 0; i < Cell_Grid[index].Total_Divisions; i++)
	{
		Cell_Grid[index].Particles_Per_Cell[i] = 0;
		Cell_Grid[index].Colloids_Per_Cell[i] = 0;
		Cell_Grid[index].Neighbor_Has_Colloid[i] = false;
	}

	X_Divisor = Region_Size[0] / Div_X;
	Y_Divisor = Region_Size[1] / Div_Y;

	for (int i = 0; i < Number_Particles; i++)
	{
		// assign each particle to a cell
		// keeping track of the counts per cell

		Local_Cell_X = (int)floor(Particles[i].Position_X / X_Divisor);
		Local_Cell_Y = (int)floor(Particles[i].Position_Y / Y_Divisor);

		Local_Cell_Y *= Div_X;

		Final_Cell = Local_Cell_Y + Local_Cell_X;

		Cell_Grid[index].Particle_Cell_List[Final_Cell][Cell_Grid[index].Particles_Per_Cell[Final_Cell]] = i;

		Cell_Grid[index].Particles_Per_Cell[Final_Cell] += 1;
	}

	// assign each colloid to a cell if in the colloid grid
	if (index == 1)
	{
		for (int i = Number_Particles; i < Number_Particles + Number_Colloids; i++)
		{
			// assign each colloid to a cell
			// keeping track of the counts per cell

			Local_Cell_X = (int)floor(Particles[i].Position_X / X_Divisor);
			Local_Cell_Y = (int)floor(Particles[i].Position_Y / Y_Divisor);

			Local_Cell_Y *= Div_X;

			Final_Cell = Local_Cell_Y + Local_Cell_X;

			Cell_Grid[index].Colloid_Cell_List[Final_Cell][Cell_Grid[index].Colloids_Per_Cell[Final_Cell]] = i;

			Cell_Grid[index].Colloids_Per_Cell[Final_Cell] += 1;
		}

		for (int i = 0; i < Cell_Grid[index].Total_Divisions; i++)
		{
			for (int j = 0; j < Neighbor_Count; j++)
			{
				if (Cell_Grid[index].Colloids_Per_Cell[Cell_Grid[index].Neighbor_List[i][j]] > 0)
				{
					Cell_Grid[index].Neighbor_Has_Colloid[i] = true;
				}
			}
		}
	}
}

void Build_Cell_List_3D(int index)
{
	double X_Divisor, Y_Divisor, Z_Divisor;

	int Local_Cell_X, Local_Cell_Y, Local_Cell_Z;

	int Final_Cell;

	int Div_X = Cell_Grid[index].Divisions[0];
	int Div_Y = Cell_Grid[index].Divisions[1];
	int Div_Z = Cell_Grid[index].Divisions[2];
	int Div_Tot = Cell_Grid[index].Total_Divisions;

	for (int i = 0; i < Cell_Grid[index].Total_Divisions; i++)
	{
		Cell_Grid[index].Particles_Per_Cell[i] = 0;
		Cell_Grid[index].Colloids_Per_Cell[i] = 0;
		Cell_Grid[index].Neighbor_Has_Colloid[i] = false;
	}

	int Block_Division = Div_X * Div_Y;

	X_Divisor = Region_Size[0] / Div_X;
	Y_Divisor = Region_Size[1] / Div_Y;
	Z_Divisor = Region_Size[2] / Div_Z;

	for (int i = 0; i < Number_Particles; i++)
	{
		// assign each particle to a cell
		// keeping track of the counts per cell

		Local_Cell_X = (int)floor(Particles[i].Position_X / X_Divisor);
		Local_Cell_Y = (int)floor(Particles[i].Position_Y / Y_Divisor);
		Local_Cell_Z = (int)floor(Particles[i].Position_Z / Z_Divisor);

		Local_Cell_Y *= Div_X;

		Local_Cell_Z *= Block_Division;

		Final_Cell = Local_Cell_Z + Local_Cell_Y + Local_Cell_X;

		Cell_Grid[index].Particle_Cell_List[Final_Cell][Cell_Grid[index].Particles_Per_Cell[Final_Cell]] = i;

		Cell_Grid[index].Particles_Per_Cell[Final_Cell] += 1;
	}

	// assign each colloid to a cell if in the colloid grid
	if (index == 1)
	{
		for (int i = Number_Particles; i < Number_Particles + Number_Colloids; i++)
		{
			// assign each colloid to a cell
			// keeping track of the counts per cell

			Local_Cell_X = (int)floor(Particles[i].Position_X / X_Divisor);
			Local_Cell_Y = (int)floor(Particles[i].Position_Y / Y_Divisor);
			Local_Cell_Z = (int)floor(Particles[i].Position_Z / Z_Divisor);

			Local_Cell_Y *= Div_X;

			Local_Cell_Z *= Block_Division;

			Final_Cell = Local_Cell_Z + Local_Cell_Y + Local_Cell_X;

			Cell_Grid[index].Colloid_Cell_List[Final_Cell][Cell_Grid[index].Colloids_Per_Cell[Final_Cell]] = i;

			Cell_Grid[index].Colloids_Per_Cell[Final_Cell] += 1;
		}

		for (int i = 0; i < Cell_Grid[index].Total_Divisions; i++)
		{
			for (int j = 0; j < Neighbor_Count; j++)
			{
				if (Cell_Grid[index].Colloids_Per_Cell[Cell_Grid[index].Neighbor_List[i][j]] > 0)
				{
					Cell_Grid[index].Neighbor_Has_Colloid[i] = true;
				}
			}
		}
	}
}

void Calculate_Forces_2D()
{
	int index = 1;

	for (int i = 0; i < Cell_Grid[index].Total_Divisions; i++)
	{
		if (Cell_Grid[index].Colloids_Per_Cell[i] != 0)
		{
			// note this only runs if more than 1 colloid is in the cell
			for (int j = 0; j < Cell_Grid[index].Colloids_Per_Cell[i] - 1; j++)
			{
				// calculate forces on colloids in the current cell with each neighbor in this cell
				for (int k = j + 1; k < Cell_Grid[index].Colloids_Per_Cell[i]; k++)
				{
					Calculate_Colloid_Forces_2D(Cell_Grid[index].Colloid_Cell_List[i][j], Cell_Grid[index].Colloid_Cell_List[i][k], 0, 0);
				}
			}

			for (int j = 0; j < Cell_Grid[index].Colloids_Per_Cell[i]; j++)
			{
				// calculate forces on particles in the current cell with each neighbor in this cell
				for (int k = 0; k < Cell_Grid[index].Particles_Per_Cell[i]; k++)
				{
					Calculate_Colloid_Forces_2D(Cell_Grid[index].Colloid_Cell_List[i][j], Cell_Grid[index].Particle_Cell_List[i][k], 0, 0);
				}
			}

			for (int n = 0; n < Neighbor_Count; n++)
			{
				if (Cell_Grid[index].Colloids_Per_Cell[Cell_Grid[index].Neighbor_List[i][n]] != 0)
				{
					for (int j = 0; j < Cell_Grid[index].Colloids_Per_Cell[i]; j++)
					{
						// calculate forces on Colloids in the current cell with each neighboring cell
						for (int k = 0; k < Cell_Grid[index].Colloids_Per_Cell[Cell_Grid[index].Neighbor_List[i][n]]; k++)
						{
							Calculate_Colloid_Forces_2D(Cell_Grid[index].Colloid_Cell_List[i][j], Cell_Grid[index].Colloid_Cell_List[Cell_Grid[index].Neighbor_List[i][n]][k],
								Cell_Grid[index].Neighbor_X_Offset[i][n], Cell_Grid[index].Neighbor_Y_Offset[i][n]);
						}
					}
				}

				if (Cell_Grid[index].Particles_Per_Cell[Cell_Grid[index].Neighbor_List[i][n]] != 0)
				{
					for (int j = 0; j < Cell_Grid[index].Colloids_Per_Cell[i]; j++)
					{
						// calculate forces on particles in the current cell with each neighboring cell
						for (int k = 0; k < Cell_Grid[index].Particles_Per_Cell[Cell_Grid[index].Neighbor_List[i][n]]; k++)
						{
							Calculate_Colloid_Forces_2D(Cell_Grid[index].Colloid_Cell_List[i][j], Cell_Grid[index].Particle_Cell_List[Cell_Grid[index].Neighbor_List[i][n]][k],
								Cell_Grid[index].Neighbor_X_Offset[i][n], Cell_Grid[index].Neighbor_Y_Offset[i][n]);
						}
					}
				}
			}
		}

		if (Cell_Grid[index].Neighbor_Has_Colloid[i])
		{
			for (int n = 0; n < Neighbor_Count; n++)
			{
				if (Cell_Grid[index].Colloids_Per_Cell[Cell_Grid[index].Neighbor_List[i][n]] != 0)
				{
					for (int j = 0; j < Cell_Grid[index].Particles_Per_Cell[i]; j++)
					{
						// calculate forces on Colloids in the current cell with each neighboring cell
						for (int k = 0; k < Cell_Grid[index].Colloids_Per_Cell[Cell_Grid[index].Neighbor_List[i][n]]; k++)
						{
							Calculate_Colloid_Forces_2D(Cell_Grid[index].Particle_Cell_List[i][j], Cell_Grid[index].Colloid_Cell_List[Cell_Grid[index].Neighbor_List[i][n]][k],
								Cell_Grid[index].Neighbor_X_Offset[i][n], Cell_Grid[index].Neighbor_Y_Offset[i][n]);
						}
					}
				}
			}
		}
	}
	
	index = 0;

	for (int i = 0; i < Cell_Grid[index].Total_Divisions; i++)
	{
		if (Cell_Grid[index].Particles_Per_Cell[i] != 0)
		{
			// note this only runs if more than 1 particle is in the cell
			for (int j = 0; j < Cell_Grid[index].Particles_Per_Cell[i] - 1; j++)
			{
				// calculate forces on particles in the current cell with each neighbor in this cell
				for (int k = j + 1; k < Cell_Grid[index].Particles_Per_Cell[i]; k++)
				{
					Calculate_Particle_Forces_2D(Cell_Grid[index].Particle_Cell_List[i][j], Cell_Grid[index].Particle_Cell_List[i][k], 0, 0);
				}
			}

			for (int n = 0; n < Neighbor_Count; n++)
			{
				if (Cell_Grid[index].Particles_Per_Cell[Cell_Grid[index].Neighbor_List[i][n]] != 0)
				{
					for (int j = 0; j < Cell_Grid[index].Particles_Per_Cell[i]; j++)
					{
						// calculate forces on particles in the current cell with each neighboring cell
						for (int k = 0; k < Cell_Grid[index].Particles_Per_Cell[Cell_Grid[index].Neighbor_List[i][n]]; k++)
						{
							Calculate_Particle_Forces_2D(Cell_Grid[index].Particle_Cell_List[i][j], Cell_Grid[index].Particle_Cell_List[Cell_Grid[index].Neighbor_List[i][n]][k],
								Cell_Grid[index].Neighbor_X_Offset[i][n], Cell_Grid[index].Neighbor_Y_Offset[i][n]);
						}
					}
				}
			}
		}
	}
}

void Calculate_Forces_3D()
{
	int index = 1;

	for (int i = 0; i < Cell_Grid[index].Total_Divisions; i++)
	{
		if (Cell_Grid[index].Colloids_Per_Cell[i] != 0)
		{
			// note this only runs if more than 1 colloid is in the cell
			for (int j = 0; j < Cell_Grid[index].Colloids_Per_Cell[i] - 1; j++)
			{
				// calculate forces on colloids in the current cell with each neighbor in this cell
				for (int k = j + 1; k < Cell_Grid[index].Colloids_Per_Cell[i]; k++)
				{
					Calculate_Colloid_Forces_3D(Cell_Grid[index].Colloid_Cell_List[i][j], Cell_Grid[index].Colloid_Cell_List[i][k], 0, 0, 0);
				}
			}

			for (int j = 0; j < Cell_Grid[index].Colloids_Per_Cell[i]; j++)
			{
				// calculate forces on particles in the current cell with each neighbor in this cell
				for (int k = 0; k < Cell_Grid[index].Particles_Per_Cell[i]; k++)
				{
					Calculate_Colloid_Forces_3D(Cell_Grid[index].Colloid_Cell_List[i][j], Cell_Grid[index].Particle_Cell_List[i][k], 0, 0, 0);
				}
			}

			for (int n = 0; n < Neighbor_Count; n++)
			{
				if (Cell_Grid[index].Colloids_Per_Cell[Cell_Grid[index].Neighbor_List[i][n]] != 0)
				{
					for (int j = 0; j < Cell_Grid[index].Colloids_Per_Cell[i]; j++)
					{
						// calculate forces on Colloids in the current cell with each neighboring cell
						for (int k = 0; k < Cell_Grid[index].Colloids_Per_Cell[Cell_Grid[index].Neighbor_List[i][n]]; k++)
						{
							Calculate_Colloid_Forces_3D(Cell_Grid[index].Colloid_Cell_List[i][j], Cell_Grid[index].Colloid_Cell_List[Cell_Grid[index].Neighbor_List[i][n]][k],
								Cell_Grid[index].Neighbor_X_Offset[i][n], Cell_Grid[index].Neighbor_Y_Offset[i][n], Cell_Grid[index].Neighbor_Z_Offset[i][n]);
						}
					}
				}

				if (Cell_Grid[index].Particles_Per_Cell[Cell_Grid[index].Neighbor_List[i][n]] != 0)
				{
					for (int j = 0; j < Cell_Grid[index].Colloids_Per_Cell[i]; j++)
					{
						// calculate forces on particles in the current cell with each neighboring cell
						for (int k = 0; k < Cell_Grid[index].Particles_Per_Cell[Cell_Grid[index].Neighbor_List[i][n]]; k++)
						{
							Calculate_Colloid_Forces_3D(Cell_Grid[index].Colloid_Cell_List[i][j], Cell_Grid[index].Particle_Cell_List[Cell_Grid[index].Neighbor_List[i][n]][k],
								Cell_Grid[index].Neighbor_X_Offset[i][n], Cell_Grid[index].Neighbor_Y_Offset[i][n], Cell_Grid[index].Neighbor_Z_Offset[i][n]);
						}
					}
				}
			}
		}

		if (Cell_Grid[index].Neighbor_Has_Colloid[i])
		{
			for (int n = 0; n < Neighbor_Count; n++)
			{
				if (Cell_Grid[index].Colloids_Per_Cell[Cell_Grid[index].Neighbor_List[i][n]] != 0)
				{
					for (int j = 0; j < Cell_Grid[index].Particles_Per_Cell[i]; j++)
					{
						// calculate forces on Colloids in the current cell with each neighboring cell
						for (int k = 0; k < Cell_Grid[index].Colloids_Per_Cell[Cell_Grid[index].Neighbor_List[i][n]]; k++)
						{
							Calculate_Colloid_Forces_3D(Cell_Grid[index].Particle_Cell_List[i][j], Cell_Grid[index].Colloid_Cell_List[Cell_Grid[index].Neighbor_List[i][n]][k],
								Cell_Grid[index].Neighbor_X_Offset[i][n], Cell_Grid[index].Neighbor_Y_Offset[i][n], Cell_Grid[index].Neighbor_Z_Offset[i][n]);
						}
					}
				}
			}
		}
	}

	index = 0;

	for (int i = 0; i < Cell_Grid[index].Total_Divisions; i++)
	{
		if (Cell_Grid[index].Particles_Per_Cell[i] != 0)
		{
			// note this only runs if more than 1 particle is in the cell
			for (int j = 0; j < Cell_Grid[index].Particles_Per_Cell[i] - 1; j++)
			{
				// calculate forces on particles in the current cell with each neighbor in this cell
				for (int k = j + 1; k < Cell_Grid[index].Particles_Per_Cell[i]; k++)
				{
					Calculate_Particle_Forces_3D(Cell_Grid[index].Particle_Cell_List[i][j], Cell_Grid[index].Particle_Cell_List[i][k], 0, 0, 0);
				}
			}

			for (int n = 0; n < Neighbor_Count; n++)
			{
				if (Cell_Grid[index].Particles_Per_Cell[Cell_Grid[index].Neighbor_List[i][n]] != 0)
				{
					for (int j = 0; j < Cell_Grid[index].Particles_Per_Cell[i]; j++)
					{
						// calculate forces on particles in the current cell with each neighboring cell
						for (int k = 0; k < Cell_Grid[index].Particles_Per_Cell[Cell_Grid[index].Neighbor_List[i][n]]; k++)
						{
							Calculate_Particle_Forces_3D(Cell_Grid[index].Particle_Cell_List[i][j], Cell_Grid[index].Particle_Cell_List[Cell_Grid[index].Neighbor_List[i][n]][k],
								Cell_Grid[index].Neighbor_X_Offset[i][n], Cell_Grid[index].Neighbor_Y_Offset[i][n], Cell_Grid[index].Neighbor_Z_Offset[i][n]);
						}
					}
				}
			}
		}
	}
}

void Calculate_Particle_Forces_2D(int Index_1, int Index_2, double Offset_X, double Offset_Y)
{
	double Vector_X = Particles[Index_2].Position_X - (Particles[Index_1].Position_X + Offset_X);
	double Vector_Y = Particles[Index_2].Position_Y - (Particles[Index_1].Position_Y + Offset_Y);
	
	// normal
	double r_square = Vector_X * Vector_X + Vector_Y * Vector_Y;
	double Temp_Force_X;
	double Temp_Force_Y;


	double base_calc;
	double temp_s;

	if ((Particles[Index_1].Type == 1) && (Particles[Index_2].Type == 1))
	{
		// LJ potential
		double f_lj = 0;

		if (r_square < Sigma_Test_LJ)
		{
			temp_s = Sigma_Square / r_square;
			
			base_calc = temp_s * temp_s * temp_s;
			
			f_lj = Epsilon_LJ_x_24 / r_square * (-2.0 * base_calc * base_calc + base_calc);

			// normal
			Temp_Force_X = Vector_X * f_lj;
			Temp_Force_Y = Vector_Y * f_lj;

			Particles[Index_1].Force_X += Temp_Force_X;
			Particles[Index_2].Force_X -= Temp_Force_X;

			Particles[Index_1].Force_Y += Temp_Force_Y;
			Particles[Index_2].Force_Y -= Temp_Force_Y;
		}
	}
	else
	{
		// WCA potential
		double f_wca = 0;

		if (r_square < Sigma_Test_WCA)
		{
			temp_s = Sigma_Square / r_square;
			
			base_calc = temp_s * temp_s * temp_s;
			
			f_wca = Epsilon_WCA_x_24 / r_square * (-2.0 * base_calc * base_calc + base_calc);

			// normal
			Temp_Force_X = Vector_X * f_wca;
			Temp_Force_Y = Vector_Y * f_wca;

			Particles[Index_1].Force_X += Temp_Force_X;
			Particles[Index_2].Force_X -= Temp_Force_X;

			Particles[Index_1].Force_Y += Temp_Force_Y;
			Particles[Index_2].Force_Y -= Temp_Force_Y;
		}
	}
}

void Calculate_Colloid_Forces_2D(int Index_1, int Index_2, double Offset_X, double Offset_Y)
{
	double Vector_X = Particles[Index_2].Position_X - (Particles[Index_1].Position_X + Offset_X);
	double Vector_Y = Particles[Index_2].Position_Y - (Particles[Index_1].Position_Y + Offset_Y);

	// normal
	double r_square = Vector_X * Vector_X + Vector_Y * Vector_Y;

	// WCA potential

	double f_wca = 0;

	double temp_s;

	bool Particle_Present = false;

	double Test_Sigma = Sigma_Test_WCA_Colloid_Particle;
	double Sigma_Square_Use = Sigma_Square_Colloid_Particle;

	if (Index_1 < Number_Particles)
	{
		if (r_square < Reaction_Radius_Squared)
		{
			Particles[Index_1].Inside_Reaction = true;
		}

		Particle_Present = true;
	}

	if (Index_2 < Number_Particles)
	{
		if (r_square < Reaction_Radius_Squared)
		{
			Particles[Index_2].Inside_Reaction = true;
		}

		Particle_Present = true;
	}

	if (!Particle_Present)
	{
		Test_Sigma = Sigma_Test_WCA_Colloid_Colloid;
		Sigma_Square_Use = Sigma_Square_Colloid_Colloid;
	}

	if (r_square < Test_Sigma)
	{
		temp_s = Sigma_Square_Use / r_square;
		
		double base_calc = temp_s * temp_s * temp_s;
		
		f_wca = Epsilon_WCA_x_24 / r_square * (-2.0 * base_calc * base_calc + base_calc);

		// normal
		double Temp_Force_X = Vector_X * f_wca;
		double Temp_Force_Y = Vector_Y * f_wca;

		Particles[Index_1].Force_X += Temp_Force_X;
		Particles[Index_2].Force_X -= Temp_Force_X;

		Particles[Index_1].Force_Y += Temp_Force_Y;
		Particles[Index_2].Force_Y -= Temp_Force_Y;
	}
}

void Calculate_Particle_Forces_3D(int Index_1, int Index_2, double Offset_X, double Offset_Y, double Offset_Z)
{
	double Vector_X = Particles[Index_2].Position_X - (Particles[Index_1].Position_X + Offset_X);
	double Vector_Y = Particles[Index_2].Position_Y - (Particles[Index_1].Position_Y + Offset_Y);
	double Vector_Z = Particles[Index_2].Position_Z - (Particles[Index_1].Position_Z + Offset_Z);

	double r_square = pow((Vector_X), 2) + pow((Vector_Y), 2) + pow((Vector_Z), 2);

	double Temp_Force_X;
	double Temp_Force_Y;
	double Temp_Force_Z;
	double base_calc;

	if ((Particles[Index_1].Type == 1) && (Particles[Index_2].Type == 1))
	{
		// LJ potential
		double f_lj = 0;

		if (r_square < Sigma_Test_LJ)
		{
			base_calc = pow(Sigma_Square / r_square, 3);

			f_lj = Epsilon_LJ_x_24 / r_square * (-2.0 * pow(base_calc, 2) + base_calc);

			Temp_Force_X = Vector_X * f_lj;
			Temp_Force_Y = Vector_Y * f_lj;
			Temp_Force_Z = Vector_Z * f_lj;

			Particles[Index_1].Force_X += Temp_Force_X;
			Particles[Index_2].Force_X -= Temp_Force_X;

			Particles[Index_1].Force_Y += Temp_Force_Y;
			Particles[Index_2].Force_Y -= Temp_Force_Y;

			Particles[Index_1].Force_Z += Temp_Force_Z;
			Particles[Index_2].Force_Z -= Temp_Force_Z;
		}
	}
	else
	{
		// WCA potential
		double f_wca = 0;

		if (r_square < Sigma_Test_WCA)
		{
			base_calc = pow(Sigma_Square / r_square, 3);

			f_wca = Epsilon_WCA_x_24 / r_square * (-2.0 * pow(base_calc, 2) + base_calc);

			Temp_Force_X = Vector_X * f_wca;
			Temp_Force_Y = Vector_Y * f_wca;
			Temp_Force_Z = Vector_Z * f_wca;

			Particles[Index_1].Force_X += Temp_Force_X;
			Particles[Index_2].Force_X -= Temp_Force_X;

			Particles[Index_1].Force_Y += Temp_Force_Y;
			Particles[Index_2].Force_Y -= Temp_Force_Y;

			Particles[Index_1].Force_Z += Temp_Force_Z;
			Particles[Index_2].Force_Z -= Temp_Force_Z;
		}
	}
}

void Calculate_Colloid_Forces_3D(int Index_1, int Index_2, double Offset_X, double Offset_Y, double Offset_Z)
{
	double Vector_X = Particles[Index_2].Position_X - (Particles[Index_1].Position_X + Offset_X);
	double Vector_Y = Particles[Index_2].Position_Y - (Particles[Index_1].Position_Y + Offset_Y);
	double Vector_Z = Particles[Index_2].Position_Z - (Particles[Index_1].Position_Z + Offset_Z);

	double r_square = pow((Vector_X), 2) + pow((Vector_Y), 2) + pow((Vector_Z), 2);

	// WCA potential

	double f_wca = 0;

	bool Particle_Present = false;

	double Test_Sigma = Sigma_Test_WCA_Colloid_Particle;

	if (Index_1 < Number_Particles)
	{
		if (r_square < Reaction_Radius_Squared)
		{
			Particles[Index_1].Inside_Reaction = true;

			//cout << "1: " << Index_1 << endl;
		}

		Particle_Present = true;


	}

	if (Index_2 < Number_Particles)
	{
		if (r_square < Reaction_Radius_Squared)
		{
			Particles[Index_2].Inside_Reaction = true;

			//cout << "2: " << Index_2 << endl;
		}

		Particle_Present = true;
	}

	if (!Particle_Present)
	{
		Test_Sigma = Sigma_Test_WCA_Colloid_Colloid;
	}

	if (r_square < Test_Sigma)
	{
		double base_calc = pow(Test_Sigma / r_square, 3);

		f_wca = Epsilon_WCA_x_24 / r_square * (-2.0 * pow(base_calc, 2) + base_calc);

		double Temp_Force_X = Vector_X * f_wca;
		double Temp_Force_Y = Vector_Y * f_wca;
		double Temp_Force_Z = Vector_Z * f_wca;

		Particles[Index_1].Force_X += Temp_Force_X;
		Particles[Index_2].Force_X -= Temp_Force_X;

		Particles[Index_1].Force_Y += Temp_Force_Y;
		Particles[Index_2].Force_Y -= Temp_Force_Y;

		Particles[Index_1].Force_Z += Temp_Force_Z;
		Particles[Index_2].Force_Z -= Temp_Force_Z;
	}
}