// SelfPropelledProfiles.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <stdio.h>
#include <vector>
#include <math.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include <omp.h>
#include <random>
#include <string>  

#define M_PI 3.14159265359




//-----���������_�������----------------------


int Number_of_particles = 2000;
double ParticlesPerUnit = 20.;
int RadiusSystem = int(sqrt(Number_of_particles / ParticlesPerUnit)) + 1;
int Accuracy = 150, TimesRadius = RadiusSystem;
int RadiusSystem2 = RadiusSystem*RadiusSystem;
double radius_step = (double)RadiusSystem / TimesRadius;
double RadiusInteraction = 1.;
double dt = 1. / 60;
double scalar_noise_amp_init = 0.002,
		vector_noise_amp_init = 0.002,
		eps = 1. / sqrt(Number_of_particles);

void evolution(int, double, double);

int main(void)
{
	for (double scalar_noise = scalar_noise_amp_init + 0.1 ; scalar_noise < 0.6; scalar_noise += 0.1){
		for (double vector_noise = vector_noise_amp_init ; vector_noise < 0.6; vector_noise += 0.1)
			evolution(Number_of_particles, scalar_noise, vector_noise);
	}
	system("PAUSE");
	return 0;
}

void evolution(int ParticleNum, double scalar_noise_amp, double vector_noise_amp){
	
	std::ofstream result("D:\\csd\\result_csd scalar_noise_amp= " + std::to_string(scalar_noise_amp) + " ,scalar_noise_amp= " + std::to_string(vector_noise_amp) + ".txt"),
		profiles("D:\\csd\\profiles_csd scalar_noise_amp= " + std::to_string(scalar_noise_amp) + ", scalar_noise_amp= " + std::to_string(vector_noise_amp) + ".txt"),
		parametrs("D:\\csd\\InitialParametrs_csd scalar_noise_amp= " + std::to_string(scalar_noise_amp) + ", scalar_noise_amp= " + std::to_string(vector_noise_amp) + ".txt");

	parametrs << "Number Of Particles " << Number_of_particles << '\n' <<
		"Radius System " << RadiusSystem << '\n' <<
		"Density " << ParticlesPerUnit << '\n' <<
		"Scalar noise amp " << scalar_noise_amp << '\n' <<
		"Vector noise amp " << vector_noise_amp << '\n';

	double sx = 0., sy = 0., Vx = 0., Vy = 0.;
	int Iteration = 0;

	std::vector<double> X(ParticleNum);
	std::vector<double> Y(ParticleNum);
	std::vector<double> VX(ParticleNum);
	std::vector<double> VY(ParticleNum);
	std::vector<double> OrtsVX(ParticleNum);
	std::vector<double> OrtsVY(ParticleNum);

	double* Amount_Particle = new double[TimesRadius];
	double* MODULE_VELOCITY_HIST = new double[TimesRadius];
	double* MOMENTUM_VELOCITY_HIST = new double[TimesRadius];


	result << "x\ty" << std::endl;
	profiles << "time\tstep\tdensv\tmodule\tmomentum\n";
	std::default_random_engine genereator;
	std::chi_squared_distribution<double> distribution(4.0);
	std::srand((unsigned)time(NULL));
	//-----initial distribution
	for (int i = 0; i < ParticleNum; i++){
		double phi = 2. * M_PI * ((double)(rand()) / RAND_MAX);
		double psi = 2. * M_PI * ((double)(rand()) / RAND_MAX);
		double e = ((double)(rand()) / RAND_MAX);
		X[i] = 0.95*RadiusSystem * e * cos(phi);
		Y[i] = 0.95*RadiusSystem * e * sin(phi);

		double ABS = distribution(genereator) + 2.; 
		VX[i] = ABS*cos(psi);
		VY[i] = ABS*sin(psi);
	}


	int *block = new int[RadiusSystem2 + 1];
	int *lscl = new int[ParticleNum];
	int **cellb = new int*[RadiusSystem2];


	for (int i = 0; i < RadiusSystem2; i++)
	{
		cellb[i] = new int[9];
	}

	int cellv0[9][2] =
	{ { RadiusSystem - 1, RadiusSystem - 1 }, { 0, RadiusSystem - 1 }, { 1, RadiusSystem - 1 }, { RadiusSystem - 1, 0 }, { 0, 0 }, { 1, 0 }, { RadiusSystem - 1, 1 }, { 0, 1 }, { 1, 1 } };

	for (int b = 0; b < RadiusSystem2; b++){
		for (int j = 0; j < 9; j++){
			int bx = cellv0[j][0] + b % RadiusSystem; int by = (cellv0[j][1] + b / RadiusSystem);
			cellb[b][j] = bx%RadiusSystem + (by%RadiusSystem)*RadiusSystem;

		}
	}

	double L0 = 0.,
		L100 = 1e5,
		time_elapsed = 0.;
	int Number_of_steps = 0,n = 0;

	while (fabs(L0 - L100) > eps && Iteration <= Number_of_particles){
		time_elapsed += dt;
		Iteration++;

		if (Number_of_steps%Accuracy == (Accuracy - 1)) {
			L100 = 0.;
			for (int i = 0; i < ParticleNum; i++){
				L100 += fabs((X[i] * VY[i] - Y[i] * VX[i]) / (X[i] * X[i] + Y[i] * Y[i]));
			}
		}

		if (Number_of_steps%Accuracy == 0){
			L0 = 0.;
			for (int i = 0; i < ParticleNum; i++){
				L0 += fabs((X[i] * VY[i] - Y[i] * VX[i]) / (X[i] * X[i] + Y[i] * Y[i]));
			}
		}
		for (int k = 0; k < TimesRadius; k++){
			Amount_Particle[k] = 0.;
			MODULE_VELOCITY_HIST[k] = 0.;
			MOMENTUM_VELOCITY_HIST[k] = 0.;
		}

		Number_of_steps++;

		for (int i = 0; i < RadiusSystem*RadiusSystem + 1; i++){
			block[i] = -1;
		}

		for (int i = 0; i < ParticleNum; i++){
			OrtsVX[i] = VX[i] / hypot(VX[i], VY[i]);
			OrtsVY[i] = VY[i] / hypot(VX[i], VY[i]);
		}

		for (int i = 0; i < ParticleNum; i++){
			int b = int(floor(fabs(X[i]))) + RadiusSystem * int(floor(fabs(Y[i])));
			lscl[i] = block[b];
			block[b] = i;
		}


		for (int b0 = 0; b0 < RadiusSystem2; b0++){
			
			for (int i = block[b0]; i != -1; i = lscl[i]) {
				#pragma omp parallel for
				for (int b1 = 0; b1 < 9; b1++){
					int b2 = cellb[b0][b1];

					for (int j = block[b2]; j != -1; j = lscl[j]) {
						double dx = abs(X[i] - X[j]);
						double dy = abs(Y[i] - Y[j]);
						if ((hypot(dx, dy) < RadiusInteraction) && j != i)
						{
							sx += OrtsVX[j]; sy += OrtsVY[j];
							n++;
						}
					}
				}
				double noise_action_scal = scalar_noise_amp * M_PI *(1 - 2 * ((double)rand() / (RAND_MAX)));
				double noise_action_vect_x = vector_noise_amp * M_PI *(1 - 2 * ((double)rand() / (RAND_MAX)));
				double noise_action_vect_y = vector_noise_amp * M_PI *(1 - 2 * ((double)rand() / (RAND_MAX)));

				double sx2 = sx + noise_action_vect_x * n;
				double sy2 = sy + noise_action_vect_y * n;
				if (abs(sx2) < 1e-5 && abs(sy2) < 1e-5){
					Vx = OrtsVX[i];
					Vy = OrtsVY[i];
					OrtsVX[i] = Vx * cos(noise_action_scal) - Vy * sin(noise_action_scal);
					OrtsVY[i] = Vx * sin(noise_action_scal) + Vy * cos(noise_action_scal);
				}
				else {
					OrtsVX[i] = (sx2 * cos(noise_action_scal) - sy2 * sin(noise_action_scal)) / hypot(sx2, sy2);
					OrtsVY[i] = (sx2 * sin(noise_action_scal) + sy2 * cos(noise_action_scal)) / hypot(sx2, sy2);
				}

				sx = 0., sy = 0.;
				n = 0;
				double mod = hypot(VX[i], VY[i]);
				VX[i] = OrtsVX[i] * mod;
				VY[i] = OrtsVY[i] * mod;

				//--------���� ������ 

				double r = radius_step;
				for (int q = 0; q < TimesRadius; q++){
					if ((hypot(X[i], Y[i]) < r) && (hypot(X[i], Y[i]) > (r - radius_step))) {
						Amount_Particle[q] += (1 / (double)Number_of_particles);
						MODULE_VELOCITY_HIST[q] += ((VX[i] * VX[i] + VY[i] * VY[i]) / (double)Number_of_particles);
						MOMENTUM_VELOCITY_HIST[q] += fabs((X[i] * VY[i] - Y[i] * VX[i]) / ((X[i] * X[i] + Y[i] * Y[i]) * (double)Number_of_particles));
						break;
					}
					r += radius_step;
				}
				result << X[i] << '\t' << Y[i] << '\n';
				double x1 = X[i] + VX[i] * dt; double y1 = Y[i] + VY[i] * dt;
				if (x1*x1 + y1*y1 > (double)(RadiusSystem*RadiusSystem)){
					double SQRT = pow(((X[i] * VX[i] + Y[i] * VY[i]) / (VX[i] * VX[i] + VY[i] * VY[i])), 2) + (((double)(RadiusSystem*RadiusSystem) - Y[i] * Y[i] - X[i] * X[i]) / (VX[i] * VX[i] + VY[i] * VY[i]));
					double dt1 = fabs(sqrt(fabs(SQRT)) - ((X[i] * VX[i] + Y[i] * VY[i]) / (VX[i] * VX[i] + VY[i] * VY[i])));

					double x0 = X[i] + VX[i] * dt1;
					double y0 = Y[i] + VY[i] * dt1;

					Vx = VX[i];
					Vy = VY[i];

					VX[i] = (Vx * ((double)(RadiusSystem*RadiusSystem) - 2. * x0 *x0) - 2. * Vy * x0 * y0) / (double)(RadiusSystem*RadiusSystem);
					VY[i] = (Vy * ((double)(RadiusSystem*RadiusSystem) - 2. * y0 *y0) - 2. * Vx * x0 * y0) / (double)(RadiusSystem*RadiusSystem);

					X[i] = x0 + VX[i] * (dt - dt1);
					Y[i] = y0 + VY[i] * (dt - dt1);

					if (sqrt(X[i] * X[i] + Y[i] * Y[i]) > (double)RadiusSystem){
						X[i] = x0;
						Y[i] = y0;
					}

				}
				else{
					X[i] = x1;
					Y[i] = y1;
				}
			}
		}
		for (int q = 0; q < TimesRadius; q++){
			profiles << time_elapsed << '\t' << q
				<< '\t' << Amount_Particle[q] << '\t' << MODULE_VELOCITY_HIST[q]
				<< '\t' << MOMENTUM_VELOCITY_HIST[q] << '\n';
		}
		}
		parametrs << "Iteration Number= " << Iteration << '\n';
		delete[] block;
		delete[] lscl;
		delete[] Amount_Particle;
		delete[] MODULE_VELOCITY_HIST;
		delete[] MOMENTUM_VELOCITY_HIST;
}
