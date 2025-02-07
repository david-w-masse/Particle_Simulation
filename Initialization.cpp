// Initialization code
#pragma once

#include <chrono>
#include <random>
#include <cmath>
#include <fstream>
#include <stdexcept>

#include "Particle_Simulation.h"

using namespace std;

extern Cell_Structure Cell_Grid[2];

extern Particle_Structure* Particles;

// simulation variables
extern double Particle_Diameter;
extern double Colloid_Diameter;

extern double Reaction_Radius;
extern double Change_Rate;

extern double End_Time;
extern double Delta_t;
extern double sigma;
extern double Epsilon_WCA;
extern double Epsilon_LJ;

extern int Dimension;
extern int Number_Particles;
extern int Number_Colloids;

extern double Region_Size[3];

extern double Sigma_Square;
extern double Sigma_Test_WCA;
extern double Sigma_Test_WCA_Colloid_Particle;
extern double Sigma_Test_WCA_Colloid_Colloid;
extern double Sigma_Test_LJ;
extern double Epsilon_WCA_x_24;
extern double Epsilon_LJ_x_24;
extern double Reaction_Radius_Squared;

extern int Neighbor_Count;

double X_Offsets_2D[9];
double Y_Offsets_2D[9];

void Assign_Particle_Data()
{
	for (int i = 0; i < Number_Particles; i++)
	{
		Particles[i].Position_X = 0;
		Particles[i].Position_Y = 0;
		Particles[i].Position_Z = 0;
		Particles[i].Force_X = 0;
		Particles[i].Force_Y = 0;
		Particles[i].Force_Z = 0;

		Particles[i].Size = Particle_Diameter/2.0;
		Particles[i].Type = 0;

		Particles[i].Effective_Range = 0;
		Particles[i].Diffusion = 1;

		Particles[i].Inside_Reaction = false;
	}

	for (int i = Number_Particles; i < Number_Particles + Number_Colloids; i++)
	{
		Particles[i].Position_X = 0;
		Particles[i].Position_Y = 0;
		Particles[i].Position_Z = 0;
		Particles[i].Force_X = 0;
		Particles[i].Force_Y = 0;
		Particles[i].Force_Z = 0;

		Particles[i].Size = Colloid_Diameter/2.0;
		Particles[i].Type = 2;

		Particles[i].Effective_Range = 5;
		Particles[i].Diffusion = 0.2;

		Particles[i].Inside_Reaction = false;
	}

	X_Offsets_2D[0] = 0;						Y_Offsets_2D[0] = 0;
	X_Offsets_2D[1] = Region_Size[0];			Y_Offsets_2D[1] = 0;
	X_Offsets_2D[2] = -Region_Size[0];			Y_Offsets_2D[2] = 0;
	X_Offsets_2D[3] = 0;						Y_Offsets_2D[3] = Region_Size[1];
	X_Offsets_2D[4] = 0;						Y_Offsets_2D[4] = -Region_Size[1];
	X_Offsets_2D[5] = Region_Size[0];			Y_Offsets_2D[5] = Region_Size[1];
	X_Offsets_2D[6] = Region_Size[0];			Y_Offsets_2D[6] = -Region_Size[1];
	X_Offsets_2D[7] = -Region_Size[0];			Y_Offsets_2D[7] = Region_Size[1];
	X_Offsets_2D[8] = -Region_Size[0];			Y_Offsets_2D[8] = -Region_Size[1];
}

void Initialize_Positions()
{
	std::random_device Get_Seed;
	std::mt19937 Random_Generator(Get_Seed());
	std::uniform_real_distribution<double> Random_X(0, Region_Size[0]);
	std::uniform_real_distribution<double> Random_Y(0, Region_Size[1]);
	std::uniform_real_distribution<double> Random_Z(0, Region_Size[2]);

	bool Good_Position;

	double Temp_X, Temp_Y, Temp_Z;

	double Current_X, Current_Y, Current_Z;

	double Temp_Dist_Min;

	int Tries;

	// Assign positions to larger partricles first
	// may be harder to fit them in later after small particles are placed
	if (Dimension == 2)
	{
		for (int i = Number_Particles; i < Number_Particles + Number_Colloids; i++)
		{
			Good_Position = false;

			Tries = 0;

			while (!Good_Position)
			{
				Temp_X = Random_X(Random_Generator);
				Temp_Y = Random_Y(Random_Generator);

				Good_Position = true;

				for (int j = Number_Particles; j < i; j++)
				{
					Temp_Dist_Min = (Particles[i].Size + Particles[j].Size) * 1.175;

					Current_X = Particles[j].Position_X;
					Current_Y = Particles[j].Position_Y;

					// check center and all boundaries for every particle
					for (int k = 0; k < 9; k++)
					{
						Good_Position &= Determine_Proper_Seperation_2D(Current_X, Current_Y, (Temp_X + X_Offsets_2D[k]), (Temp_Y + Y_Offsets_2D[k]), Temp_Dist_Min);
					}
				}

				if (!Good_Position)
				{
					Tries += 1;

					if ((Tries == 1000000) && !Good_Position)
					{
						throw std::invalid_argument("Error: Cannot place particle.");
					}
				}
			}

			Particles[i].Position_X = Temp_X;
			Particles[i].Position_Y = Temp_Y;
		}

		for (int i = 0; i < Number_Particles; i++)
		{
			Good_Position = false;

			Tries = 0;

			while (!Good_Position)
			{
				Temp_X = Random_X(Random_Generator);
				Temp_Y = Random_Y(Random_Generator);

				Good_Position = true;

				for (int j = 0; j < i; j++)
				{
					Temp_Dist_Min = (Particles[i].Size + Particles[j].Size) * 1.175;

					Current_X = Particles[j].Position_X;
					Current_Y = Particles[j].Position_Y;

					// check center and all boundaries for every particle
					for (int k = 0; k < 9; k++)
					{
						Good_Position &= Determine_Proper_Seperation_2D(Current_X, Current_Y, (Temp_X + X_Offsets_2D[k]), (Temp_Y + Y_Offsets_2D[k]), Temp_Dist_Min);
					}
				}

				for (int j = Number_Particles; j < Number_Particles + Number_Colloids; j++)
				{
					Temp_Dist_Min = (Particles[i].Size + Particles[j].Size) * 1.175;

					Current_X = Particles[j].Position_X;
					Current_Y = Particles[j].Position_Y;

					// check center and all boundaries for every particle
					for (int k = 0; k < 9; k++)
					{
						Good_Position &= Determine_Proper_Seperation_2D(Current_X, Current_Y, (Temp_X + X_Offsets_2D[k]), (Temp_Y + Y_Offsets_2D[k]), Temp_Dist_Min);
					}
				}

				if (!Good_Position)
				{
					Tries += 1;

					if ((Tries == 1000000) && !Good_Position)
					{
						throw std::invalid_argument("Error: Cannot place particle.");
					}
				}
			}

			cout << "particle placed " << i << endl; 

			Particles[i].Position_X = Temp_X;
			Particles[i].Position_Y = Temp_Y;
		}
	}
	else
	{
		for (int i = Number_Particles; i < Number_Particles + Number_Colloids; i++)
		{
			Good_Position = false;

			Tries = 0;

			while (!Good_Position)
			{
				Temp_X = Random_X(Random_Generator);
				Temp_Y = Random_Y(Random_Generator);
				Temp_Z = Random_Y(Random_Generator);

				Good_Position = true;

				for (int j = Number_Particles; j < i; j++)
				{
					Temp_Dist_Min = (Particles[i].Size + Particles[j].Size) * 1.25;

					Current_X = Particles[j].Position_X;
					Current_Y = Particles[j].Position_Y;
					Current_Z = Particles[j].Position_Z;

					// check center and all boundaries for every particle
					for (int k = 0; k < 9; k++)
					{
						Good_Position &= Determine_Proper_Seperation_3D(Current_X, Current_Y, Current_Z,
										(Temp_X + X_Offsets_2D[k]), (Temp_Y + Y_Offsets_2D[k]), Temp_Z, Temp_Dist_Min);
					}

					for (int k = 0; k < 9; k++)
					{
						Good_Position &= Determine_Proper_Seperation_3D(Current_X, Current_Y, Current_Z,
										(Temp_X + X_Offsets_2D[k]), (Temp_Y + Y_Offsets_2D[k]), (Temp_Z + Region_Size[2]), Temp_Dist_Min);
					}

					for (int k = 0; k < 9; k++)
					{
						Good_Position &= Determine_Proper_Seperation_3D(Current_X, Current_Y, Current_Z,
										(Temp_X + X_Offsets_2D[k]), (Temp_Y + Y_Offsets_2D[k]), (Temp_Z - Region_Size[2]), Temp_Dist_Min);
					}
				}

				if (!Good_Position)
				{
					Tries += 1;

					if ((Tries == 1000000) && !Good_Position)
					{
						throw std::invalid_argument("Error: Cannot place particle.");
					}
				}
			}

			Particles[i].Position_X = Temp_X;
			Particles[i].Position_Y = Temp_Y;
			Particles[i].Position_Z = Temp_Z;
		}
		
		for (int i = 0; i < Number_Particles; i++)
		{
			Good_Position = false;

			Tries = 0;

			while (!Good_Position)
			{
				Temp_X = Random_X(Random_Generator);
				Temp_Y = Random_Y(Random_Generator);
				Temp_Z = Random_Y(Random_Generator);

				Good_Position = true;

				for (int j = 0; j < i; j++)
				{
					Temp_Dist_Min = (Particles[i].Size + Particles[j].Size) * 1.25;

					Current_X = Particles[j].Position_X;
					Current_Y = Particles[j].Position_Y;
					Current_Z = Particles[j].Position_Z;

					// check center and all boundaries for every particle
					for (int k = 0; k < 9; k++)
					{
						Good_Position &= Determine_Proper_Seperation_3D(Current_X, Current_Y, Current_Z,
										(Temp_X + X_Offsets_2D[k]), (Temp_Y + Y_Offsets_2D[k]), Temp_Z, Temp_Dist_Min);
					}

					for (int k = 0; k < 9; k++)
					{
						Good_Position &= Determine_Proper_Seperation_3D(Current_X, Current_Y, Current_Z,
										(Temp_X + X_Offsets_2D[k]), (Temp_Y + Y_Offsets_2D[k]), (Temp_Z + Region_Size[2]), Temp_Dist_Min);
					}

					for (int k = 0; k < 9; k++)
					{
						Good_Position &= Determine_Proper_Seperation_3D(Current_X, Current_Y, Current_Z,
										(Temp_X + X_Offsets_2D[k]), (Temp_Y + Y_Offsets_2D[k]), (Temp_Z - Region_Size[2]), Temp_Dist_Min);
					}
				}

				for (int j = Number_Particles; j < Number_Particles + Number_Colloids; j++)
				{
					Temp_Dist_Min = (Particles[i].Size + Particles[j].Size) * 1.25;

					Current_X = Particles[j].Position_X;
					Current_Y = Particles[j].Position_Y;
					Current_Z = Particles[j].Position_Z;

					// check center and all boundaries for every particle
					for (int k = 0; k < 9; k++)
					{
						Good_Position &= Determine_Proper_Seperation_3D(Current_X, Current_Y, Current_Z,
										(Temp_X + X_Offsets_2D[k]), (Temp_Y + Y_Offsets_2D[k]), Temp_Z, Temp_Dist_Min);
					}

					for (int k = 0; k < 9; k++)
					{
						Good_Position &= Determine_Proper_Seperation_3D(Current_X, Current_Y, Current_Z,
										(Temp_X + X_Offsets_2D[k]), (Temp_Y + Y_Offsets_2D[k]), (Temp_Z + Region_Size[2]), Temp_Dist_Min);
					}

					for (int k = 0; k < 9; k++)
					{
						Good_Position &= Determine_Proper_Seperation_3D(Current_X, Current_Y, Current_Z,
										(Temp_X + X_Offsets_2D[k]), (Temp_Y + Y_Offsets_2D[k]), (Temp_Z - Region_Size[2]), Temp_Dist_Min);
					}
				}

				if (!Good_Position)
				{
					Tries += 1;

					if ((Tries == 1000000) && !Good_Position)
					{
						throw std::invalid_argument("Error: Cannot place particle.");
					}
				}
			}

			Particles[i].Position_X = Temp_X;
			Particles[i].Position_Y = Temp_Y;
			Particles[i].Position_Z = Temp_Z;
		}
	}
}

bool Determine_Proper_Seperation_2D(double Particle_X, double Particle_Y, double Test_X, double Test_Y, double Min_Distance)
{
	double Diff_X, Diff_Y, Temp_Dist;

	Diff_X = Particle_X - Test_X;
	Diff_Y = Particle_Y - Test_Y;

	Temp_Dist = sqrt(pow((Diff_X), 2) + pow((Diff_Y), 2));

	if (Temp_Dist < Min_Distance)
	{
		return false;
	}

	return true;
}

bool Determine_Proper_Seperation_3D(double Particle_X, double Particle_Y, double Particle_Z, double Test_X, double Test_Y, double Test_Z, double Min_Distance)
{
	double Diff_X, Diff_Y, Diff_Z, Temp_Dist;

	Diff_X = Particle_X - Test_X;
	Diff_Y = Particle_Y - Test_Y;
	Diff_Z = Particle_Z - Test_Z;

	Temp_Dist = sqrt(pow(Diff_X, 2) + pow(Diff_Y, 2) + pow(Diff_Z, 2));

	if (Temp_Dist < Min_Distance)
	{
		return false;
	}

	return true;
}

void Build_Neighbor_List(int index)
{
	int Base_Block, Neighbor_Index;

	int Div_X = Cell_Grid[index].Divisions[0];
	int Div_Y = Cell_Grid[index].Divisions[1];
	int Div_Z = Cell_Grid[index].Divisions[2];
	int Div_Tot = Cell_Grid[index].Total_Divisions;

	int Block_Size = Div_X * Div_Y;

	if (Dimension == 2)
	{
		for (int i = 0; i < Div_X; i++)
		{
			for (int j = 0; j < Div_Y; j++)
			{
				Neighbor_Index = j * Div_X + i;

				// Right Neighbor
				Cell_Grid[index].Neighbor_List[Neighbor_Index][0] = (j * Div_X + ((i + 1) % Div_X)) % Div_Tot;

				// Upper Right Neighbor
				Cell_Grid[index].Neighbor_List[Neighbor_Index][1] = (((j + 1) % Div_Y) * Div_X + ((i + 1) % Div_X)) % Div_Tot;

				// Upper Neighbor
				Cell_Grid[index].Neighbor_List[Neighbor_Index][2] = (((j + 1) % Div_Y) * Div_X + i) % Div_Tot;

				// Upper Left Neighbor
				if (i == 0)
				{
					Cell_Grid[index].Neighbor_List[Neighbor_Index][3] = (((j + 1) % Div_Y) * Div_X + Div_X - 1) % Div_Tot;
				}
				else
				{
					Cell_Grid[index].Neighbor_List[Neighbor_Index][3] = (((j + 1) % Div_Y) * Div_X + (i - 1)) % Div_Tot;
				}
			}
		}
	}
	else
	{
		for (int i = 0; i < Div_X; i++)
		{
			for (int j = 0; j < Div_Y; j++)
			{
				for (int k = 0; k < Div_Z; k++)
				{
					Base_Block = k * Block_Size;

					Neighbor_Index = Base_Block + j * Div_X + i;

					// Right Neighbor
					Cell_Grid[index].Neighbor_List[Neighbor_Index][0] = Base_Block + ((j * Div_X + ((i + 1) % Div_X)) % Block_Size);

					// Upper Right Neighbor
					Cell_Grid[index].Neighbor_List[Neighbor_Index][1] = Base_Block + ((((j + 1) % Div_Y) * Div_X + ((i + 1) % Div_X)) % Block_Size);

					// Upper Neighbor
					Cell_Grid[index].Neighbor_List[Neighbor_Index][2] = Base_Block + ((((j + 1) % Div_Y) * Div_X + i) % Block_Size);

					// Upper Left Neighbor
					if (i == 0)
					{
						Cell_Grid[index].Neighbor_List[Neighbor_Index][3] = Base_Block + ((((j + 1) % Div_Y) * Div_X + Div_X - 1) % Block_Size);
					}
					else
					{
						Cell_Grid[index].Neighbor_List[Neighbor_Index][3] = Base_Block + ((((j + 1) % Div_Y) * Div_X + (i - 1)) % Block_Size);
					}

					// Forward Right Neighbor
					Cell_Grid[index].Neighbor_List[Neighbor_Index][4] = ((Base_Block + Block_Size) % Div_Tot) + ((j * Div_X + ((i + 1) % Div_X)) % Block_Size);

					// Forward Upper Right Neighbor
					Cell_Grid[index].Neighbor_List[Neighbor_Index][5] = ((Base_Block + Block_Size) % Div_Tot) + ((((j + 1) % Div_Y) * Div_X + ((i + 1) % Div_X)) % Block_Size);

					// Forward Upper Neighbor
					Cell_Grid[index].Neighbor_List[Neighbor_Index][6] = ((Base_Block + Block_Size) % Div_Tot) + ((((j + 1) % Div_Y) * Div_X + i) % Block_Size);

					// Forward Upper Left Neighbor
					if (i == 0)
					{
						Cell_Grid[index].Neighbor_List[Neighbor_Index][7] = ((Base_Block + Block_Size) % Div_Tot) + ((((j + 1) % Div_Y) * Div_X + Div_X - 1) % Block_Size);
					}
					else
					{
						Cell_Grid[index].Neighbor_List[Neighbor_Index][7] = ((Base_Block + Block_Size) % Div_Tot) + ((((j + 1) % Div_Y) * Div_X + (i - 1)) % Block_Size);
					}

					// Forward Left Neighbor
					if (i == 0)
					{
						Cell_Grid[index].Neighbor_List[Neighbor_Index][8] = ((Base_Block + Block_Size) % Div_Tot) + ((j * Div_X + Div_X - 1) % Block_Size);
					}
					else
					{
						Cell_Grid[index].Neighbor_List[Neighbor_Index][8] = ((Base_Block + Block_Size) % Div_Tot) + ((j * Div_X + (i - 1)) % Block_Size);
					}

					if (j == 0)
					{
						// Forward Lower Left Neighbor
						if (i == 0)
						{
							Cell_Grid[index].Neighbor_List[Neighbor_Index][9] = ((Base_Block + Block_Size) % Div_Tot) + ((Block_Size - 1) % Block_Size);
						}
						else
						{
							Cell_Grid[index].Neighbor_List[Neighbor_Index][9] = ((Base_Block + Block_Size) % Div_Tot) + ((Block_Size - Div_X + (i - 1)) % Block_Size);
						}

						// Forward Lower Neighbor
						Cell_Grid[index].Neighbor_List[Neighbor_Index][10] = ((Base_Block + Block_Size) % Div_Tot) + ((Block_Size - Div_X + i) % Block_Size);

						// Forward Lower Right Neighbor
						Cell_Grid[index].Neighbor_List[Neighbor_Index][11] = ((Base_Block + Block_Size) % Div_Tot) + ((Block_Size - Div_X + ((i + 1) % Div_X)) % Block_Size);
					}
					else
					{
						// Forward Lower Left Neighbor
						if (i == 0)
						{
							Cell_Grid[index].Neighbor_List[Neighbor_Index][9] = ((Base_Block + Block_Size) % Div_Tot) + (((j - 1) * Div_X + Div_X - 1) % Block_Size);
						}
						else
						{
							Cell_Grid[index].Neighbor_List[Neighbor_Index][9] = ((Base_Block + Block_Size) % Div_Tot) + (((j - 1) * Div_X + (i - 1)) % Block_Size);
						}

						// Forward Lower Neighbor
						Cell_Grid[index].Neighbor_List[Neighbor_Index][10] = ((Base_Block + Block_Size) % Div_Tot) + (((j - 1) * Div_X + i) % Block_Size);

						// Forward Lower Right Neighbor
						Cell_Grid[index].Neighbor_List[Neighbor_Index][11] = ((Base_Block + Block_Size) % Div_Tot) + (((j - 1) * Div_X + ((i + 1) % Div_X)) % Block_Size);
					}

					// Forward Center Neighbor
					Cell_Grid[index].Neighbor_List[Neighbor_Index][12] = ((Base_Block + Block_Size) % Div_Tot) + ((j * Div_X + i) % Block_Size);

				}
			}
		}
	}
}

void Build_Neighbor_Offset_List(int index)
{
	int Base_Block, Neighbor_Index;

	int Div_X = Cell_Grid[index].Divisions[0];
	int Div_Y = Cell_Grid[index].Divisions[1];
	int Div_Z = Cell_Grid[index].Divisions[2];
	int Div_Tot = Cell_Grid[index].Total_Divisions;

	int Block_Size = Div_X * Div_Y;

	if (Dimension == 2)
	{
		for (int i = 0; i < Div_X; i++)
		{
			for (int j = 0; j < Div_Y; j++)
			{
				Neighbor_Index = j * Div_X + i;

				// Right Neighbor
				if (i == (Div_X - 1))
				{
					Cell_Grid[index].Neighbor_X_Offset[Neighbor_Index][0] = -Region_Size[0];
					Cell_Grid[index].Neighbor_Y_Offset[Neighbor_Index][0] = 0;
				}
				else
				{
					Cell_Grid[index].Neighbor_X_Offset[Neighbor_Index][0] = 0;
					Cell_Grid[index].Neighbor_Y_Offset[Neighbor_Index][0] = 0;
				}
				
				// Upper Right Neighbor
				if (i == (Div_X - 1))
				{
					if (j == (Div_Y - 1))
					{
						Cell_Grid[index].Neighbor_X_Offset[Neighbor_Index][1] = -Region_Size[0];
						Cell_Grid[index].Neighbor_Y_Offset[Neighbor_Index][1] = -Region_Size[1];
					}
					else
					{
						Cell_Grid[index].Neighbor_X_Offset[Neighbor_Index][1] = -Region_Size[0];
						Cell_Grid[index].Neighbor_Y_Offset[Neighbor_Index][1] = 0;
					}
				}
				else
				{
					if (j == (Div_Y - 1))
					{
						Cell_Grid[index].Neighbor_X_Offset[Neighbor_Index][1] = 0;
						Cell_Grid[index].Neighbor_Y_Offset[Neighbor_Index][1] = -Region_Size[1];
					}
					else
					{
						Cell_Grid[index].Neighbor_X_Offset[Neighbor_Index][1] = 0;
						Cell_Grid[index].Neighbor_Y_Offset[Neighbor_Index][1] = 0;
					}
				}

				// Upper Neighbor
				if (j == (Div_Y - 1))
				{
					Cell_Grid[index].Neighbor_X_Offset[Neighbor_Index][2] = 0;
					Cell_Grid[index].Neighbor_Y_Offset[Neighbor_Index][2] = -Region_Size[1];
				}
				else
				{
					Cell_Grid[index].Neighbor_X_Offset[Neighbor_Index][2] = 0;
					Cell_Grid[index].Neighbor_Y_Offset[Neighbor_Index][2] = 0;
				}

				// Upper Left Neighbor
				if (i == 0)
				{
					if (j == (Div_Y - 1))
					{
						Cell_Grid[index].Neighbor_X_Offset[Neighbor_Index][3] = Region_Size[0];
						Cell_Grid[index].Neighbor_Y_Offset[Neighbor_Index][3] = -Region_Size[1];
					}
					else
					{
						Cell_Grid[index].Neighbor_X_Offset[Neighbor_Index][3] = Region_Size[0];
						Cell_Grid[index].Neighbor_Y_Offset[Neighbor_Index][3] = 0;
					}
				}
				else
				{
					if (j == (Div_Y - 1))
					{
						Cell_Grid[index].Neighbor_X_Offset[Neighbor_Index][3] = 0;
						Cell_Grid[index].Neighbor_Y_Offset[Neighbor_Index][3] = -Region_Size[1];
					}
					else
					{
						Cell_Grid[index].Neighbor_X_Offset[Neighbor_Index][3] = 0;
						Cell_Grid[index].Neighbor_Y_Offset[Neighbor_Index][3] = 0;
					}
				}
			}
		}
	}
	else
	{
		for (int i = 0; i < Div_X; i++)
		{
			for (int j = 0; j < Div_Y; j++)
			{
				for (int k = 0; k < Div_Z; k++)
				{
					Base_Block = k * Block_Size;

					Neighbor_Index = Base_Block + j * Div_X + i;
					
					// effect on z similar for all neighbors
					for (int z = 0; z < 4; z++)
					{
						Cell_Grid[index].Neighbor_Z_Offset[Neighbor_Index][z] = 0;
					}
					
					for (int z = 4; z < 13; z++)
					{
						if (k == (Div_Z - 1))
						{
							Cell_Grid[index].Neighbor_Z_Offset[Neighbor_Index][z] = -Region_Size[2];
						}
						else
						{
							Cell_Grid[index].Neighbor_Z_Offset[Neighbor_Index][z] = 0;
						}
					}

					// Right Neighbor
					if (i == (Div_X - 1))
					{
						Cell_Grid[index].Neighbor_X_Offset[Neighbor_Index][0] = -Region_Size[0];
					}
					else
					{
						Cell_Grid[index].Neighbor_X_Offset[Neighbor_Index][0] = 0;
					}

					Cell_Grid[index].Neighbor_Y_Offset[Neighbor_Index][0] = 0;

					// Upper Right Neighbor
					if (i == (Div_X - 1))
					{
						Cell_Grid[index].Neighbor_X_Offset[Neighbor_Index][1] = -Region_Size[0];
					}
					else
					{
						Cell_Grid[index].Neighbor_X_Offset[Neighbor_Index][1] = 0;
					}

					if (j == (Div_Y - 1))
					{
						Cell_Grid[index].Neighbor_Y_Offset[Neighbor_Index][1] = -Region_Size[1];
					}
					else
					{
						Cell_Grid[index].Neighbor_Y_Offset[Neighbor_Index][1] = 0;
					}

					// Upper Neighbor
					Cell_Grid[index].Neighbor_X_Offset[Neighbor_Index][2] = 0;
					
					if (j == (Div_Y - 1))
					{
						Cell_Grid[index].Neighbor_Y_Offset[Neighbor_Index][2] = -Region_Size[1];
					}
					else
					{
						Cell_Grid[index].Neighbor_Y_Offset[Neighbor_Index][2] = 0;
					}

					// Upper Left Neighbor
					if (i == 0)
					{
						Cell_Grid[index].Neighbor_X_Offset[Neighbor_Index][3] = Region_Size[0];
					}
					else
					{
						Cell_Grid[index].Neighbor_X_Offset[Neighbor_Index][3] = 0;
					}

					if (j == (Div_Y - 1))
					{
						Cell_Grid[index].Neighbor_Y_Offset[Neighbor_Index][3] = -Region_Size[1];
					}
					else
					{
						Cell_Grid[index].Neighbor_Y_Offset[Neighbor_Index][3] = 0;
					}
				
					// Forward Right Neighbor
					if (i == (Div_X - 1))
					{
						Cell_Grid[index].Neighbor_X_Offset[Neighbor_Index][4] = -Region_Size[0];
					}
					else
					{
						Cell_Grid[index].Neighbor_X_Offset[Neighbor_Index][4] = 0;
					}

					Cell_Grid[index].Neighbor_Y_Offset[Neighbor_Index][4] = 0;

					// Forward Upper Right Neighbor
					if (i == (Div_X - 1))
					{
						Cell_Grid[index].Neighbor_X_Offset[Neighbor_Index][5] = -Region_Size[0];
					}
					else
					{
						Cell_Grid[index].Neighbor_X_Offset[Neighbor_Index][5] = 0;
					}

					if (j == (Div_Y - 1))
					{
						Cell_Grid[index].Neighbor_Y_Offset[Neighbor_Index][5] = -Region_Size[1];
					}
					else
					{
						Cell_Grid[index].Neighbor_Y_Offset[Neighbor_Index][5] = 0;
					}

					// Forward Upper Neighbor
					Cell_Grid[index].Neighbor_X_Offset[Neighbor_Index][6] = 0;

					if (j == (Div_Y - 1))
					{
						Cell_Grid[index].Neighbor_Y_Offset[Neighbor_Index][6] = -Region_Size[1];
					}
					else
					{
						Cell_Grid[index].Neighbor_Y_Offset[Neighbor_Index][6] = 0;
					}

					// Forward Upper Left Neighbor
					if (i == 0)
					{
						Cell_Grid[index].Neighbor_X_Offset[Neighbor_Index][7] = Region_Size[0];
					}
					else
					{
						Cell_Grid[index].Neighbor_X_Offset[Neighbor_Index][7] = 0;
					}

					if (j == (Div_Y - 1))
					{
						Cell_Grid[index].Neighbor_Y_Offset[Neighbor_Index][7] = -Region_Size[1];
					}
					else
					{
						Cell_Grid[index].Neighbor_Y_Offset[Neighbor_Index][7] = 0;
					}

					// Forward Left Neighbor
					if (i == 0)
					{
						Cell_Grid[index].Neighbor_X_Offset[Neighbor_Index][8] = Region_Size[0];
					}
					else
					{
						Cell_Grid[index].Neighbor_X_Offset[Neighbor_Index][8] = 0;
					}

					Cell_Grid[index].Neighbor_Y_Offset[Neighbor_Index][8] = 0;

					// Forward Lower Left Neighbor
					if (i == 0)
					{
						Cell_Grid[index].Neighbor_X_Offset[Neighbor_Index][9] = Region_Size[0];
					}
					else
					{
						Cell_Grid[index].Neighbor_X_Offset[Neighbor_Index][9] = 0;
					}

					if (j == 0)
					{
						Cell_Grid[index].Neighbor_Y_Offset[Neighbor_Index][9] = Region_Size[1];
					}
					else
					{
						Cell_Grid[index].Neighbor_Y_Offset[Neighbor_Index][9] = 0;
					}

					// Forward Lower Neighbor
					Cell_Grid[index].Neighbor_X_Offset[Neighbor_Index][10] = 0;

					if (j == 0)
					{
						Cell_Grid[index].Neighbor_Y_Offset[Neighbor_Index][10] = Region_Size[1];
					}
					else
					{
						Cell_Grid[index].Neighbor_Y_Offset[Neighbor_Index][10] = 0;
					}

					// Forward Lower Right Neighbor
					if (i == (Div_X - 1))
					{
						Cell_Grid[index].Neighbor_X_Offset[Neighbor_Index][11] = -Region_Size[0];
					}
					else
					{
						Cell_Grid[index].Neighbor_X_Offset[Neighbor_Index][11] = 0;
					}

					if (j == 0)
					{
						Cell_Grid[index].Neighbor_Y_Offset[Neighbor_Index][11] = Region_Size[1];
					}
					else
					{
						Cell_Grid[index].Neighbor_Y_Offset[Neighbor_Index][11] = 0;
					}

					// Forward Center Neighbor
					Cell_Grid[index].Neighbor_X_Offset[Neighbor_Index][12] = 0;
					Cell_Grid[index].Neighbor_Y_Offset[Neighbor_Index][12] = 0;
				}
			}
		}
	}
}

void Build_Memory()
{
	Particles = new Particle_Structure[Number_Particles + Number_Colloids];

	double Cell_Volume;
	double Particle_Volume;
	
	for (int i = 0; i < 2; i++)
	{
		// calculate how many particles can fit in the smallest grid
		// depending on dimension, then only allocate that many spots (with some to spare)
		// assume particles are smallest size
		if (Dimension == 2)
		{
			Cell_Volume = (Region_Size[0] * Region_Size[1]) / (Cell_Grid[i].Divisions[0] * Cell_Grid[i].Divisions[1]);
			Particle_Volume = Particle_Diameter * Particle_Diameter * 0.75;

			Cell_Grid[i].Max_Per_Cell = (int)floor(Cell_Volume / Particle_Volume * 2);
		}
		else
		{
			Cell_Volume = (Region_Size[0] * Region_Size[1] * Region_Size[2]) / (Cell_Grid[i].Divisions[0] * Cell_Grid[i].Divisions[1] * Cell_Grid[i].Divisions[2]);
			Particle_Volume = Particle_Diameter * Particle_Diameter * Particle_Diameter * 0.5;

			Cell_Grid[i].Max_Per_Cell = (int)floor(Cell_Volume / Particle_Volume * 2);
		}

		Cell_Grid[i].Particle_Cell_List = new int* [Cell_Grid[i].Total_Divisions];
		Cell_Grid[i].Colloid_Cell_List = new int* [Cell_Grid[i].Total_Divisions];

		Cell_Grid[i].Neighbor_List = new int* [Cell_Grid[i].Total_Divisions];

		Cell_Grid[i].Neighbor_X_Offset = new double* [Cell_Grid[i].Total_Divisions];
		Cell_Grid[i].Neighbor_Y_Offset = new double* [Cell_Grid[i].Total_Divisions];
		Cell_Grid[i].Neighbor_Z_Offset = new double* [Cell_Grid[i].Total_Divisions];

		Cell_Grid[i].Particles_Per_Cell = new int[Cell_Grid[i].Total_Divisions];
		Cell_Grid[i].Colloids_Per_Cell = new int[Cell_Grid[i].Total_Divisions];
		Cell_Grid[i].Neighbor_Has_Colloid = new bool[Cell_Grid[i].Total_Divisions];

		for (int j = 0; j < Cell_Grid[i].Total_Divisions; j++)
		{
			Cell_Grid[i].Particle_Cell_List[j] = new int[Cell_Grid[i].Max_Per_Cell];
			Cell_Grid[i].Colloid_Cell_List[j] = new int[Cell_Grid[i].Max_Per_Cell];

			Cell_Grid[i].Neighbor_List[j] = new int[Neighbor_Count];

			Cell_Grid[i].Neighbor_X_Offset[j] = new double[Neighbor_Count];
			Cell_Grid[i].Neighbor_Y_Offset[j] = new double[Neighbor_Count];
			Cell_Grid[i].Neighbor_Z_Offset[j] = new double[Neighbor_Count];

			Cell_Grid[i].Particles_Per_Cell[j] = 0;
		}
	}
}

void Clear_Memory()
{
	delete[] Particles;

	for (int i = 0; i < 2; i++)
	{
		delete[] Cell_Grid[i].Particles_Per_Cell;
		delete[] Cell_Grid[i].Colloids_Per_Cell;
		delete[] Cell_Grid[i].Neighbor_Has_Colloid;

		for (int j = 0; j < Cell_Grid[i].Total_Divisions; j++)
		{
			delete[] Cell_Grid[i].Particle_Cell_List[j];
			delete[] Cell_Grid[i].Colloid_Cell_List[j];
			delete[] Cell_Grid[i].Neighbor_List[j];
			delete[] Cell_Grid[i].Neighbor_X_Offset[j];
			delete[] Cell_Grid[i].Neighbor_Y_Offset[j];
			delete[] Cell_Grid[i].Neighbor_Z_Offset[j];
		}

		delete[] Cell_Grid[i].Particle_Cell_List;
		delete[] Cell_Grid[i].Colloid_Cell_List;
		delete[] Cell_Grid[i].Neighbor_List;
		delete[] Cell_Grid[i].Neighbor_X_Offset;
		delete[] Cell_Grid[i].Neighbor_Y_Offset;
		delete[] Cell_Grid[i].Neighbor_Z_Offset;
	}
}