// Particle_Simulation.h : Include file for standard system include files,
// or project specific include files.

#pragma once

#include <iostream>

struct Cell_Structure
{
	int Total_Divisions;

	int Max_Per_Cell;

	int Divisions[3];

	int** Particle_Cell_List;
	int** Colloid_Cell_List;

	int** Neighbor_List;

	double** Neighbor_X_Offset;
	double** Neighbor_Y_Offset;
	double** Neighbor_Z_Offset;

	int* Particles_Per_Cell;
	int* Colloids_Per_Cell;

	bool* Neighbor_Has_Colloid;

};

struct Particle_Structure
{
	double Position_X;
	double Position_Y;
	double Position_Z;
	double Force_X;
	double Force_Y;
	double Force_Z;
	double Size;
	int Type;

	double Effective_Range;
	double Diffusion;

	double Delta_X;
	double Delta_Y;
	double Delta_Z;

	bool Inside_Reaction;
};

void Build_Memory();
void Clear_Memory();

void Assign_Particle_Data();
void Initialize_Positions();

void Build_Neighbor_List(int);
void Build_Neighbor_Offset_List(int);

void Build_Cell_List_2D(int);
void Build_Cell_List_3D(int);

void Calculate_Forces_2D();
void Calculate_Forces_3D();

void Calculate_Particle_Forces_2D(int, int, double, double);
void Calculate_Particle_Forces_3D(int, int, double, double, double);

void Calculate_Colloid_Forces_2D(int, int, double, double);
void Calculate_Colloid_Forces_3D(int, int, double, double, double);

bool Determine_Proper_Seperation_2D(double, double, double, double, double);
bool Determine_Proper_Seperation_3D(double, double, double, double, double, double, double);