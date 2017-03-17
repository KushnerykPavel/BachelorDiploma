#include "stdafx.h"
#include "evolution.h"
#include <vector>
#include <math.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>


#define M_PI 3.14159265359

void Evolution(int ParticleNum,int RadiusSystem, int ThreadNumber, double noise_amp){
	std::string Noise = std::to_string(noise_amp), outNoise = "";
	for (int i = 0; i < 3; i++) outNoise += Noise[i];
	std::ofstream profiles("D:\\profiles Thread= " + std::to_string(ThreadNumber) + " Noise_amp= " + outNoise + ".txt");

	int Number_of_steps = 0,
		Accuracy = 100, 
		TimesRadius = 10;

	double L0 = 0.,
		L100 = 1e5,
		eps = 1., 
		phi, psi, e,
		time_elapsed = 0., r=0.,
		sx = 0., sy = 0., SQRT = 0.,
		Vx = 0., Vy = 0., dt1 = 0.,
		mod = 0, x0 = 0., y0 = 0., 
		noise_action =0., dx , dy ,
		dt =1./30, RadiusInteraction = 1.,
		radius_step = (double)RadiusSystem / TimesRadius;


	std::vector<double> X(ParticleNum);
	std::vector<double> Y(ParticleNum);
	std::vector<double> VX(ParticleNum);
	std::vector<double> VY(ParticleNum);
	std::vector<double> OrtsVX(ParticleNum);
	std::vector<double> OrtsVY(ParticleNum);

	double* Amount_Particle = new double[TimesRadius];
	double* MODULE_VELOCITY_HIST = new double[TimesRadius];
	double* MOMENTUM_VELOCITY_HIST = new double[TimesRadius];


	profiles << "time\tstep\tdensv\tmodule\tmomentum\n";

	std::srand((unsigned)time(NULL));
	//-----initial distribution
	for (int i = 0; i < ParticleNum; i++){
		phi = 2. * M_PI * ((double)(rand()) / RAND_MAX);
		psi = 2. * M_PI * ((double)(rand()) / RAND_MAX);
		e = ((double)(rand()) / RAND_MAX);
		X[i] = 0.05 + 0.95 * RadiusSystem * e * cos(phi);
		Y[i] = 0.05 + 0.95 * RadiusSystem * e * sin(phi);
		if(i<=ParticleNum/4){
			VX[i] = 10.*cos(psi);
			VY[i] = 10.*sin(psi);
		}
		else{
			VX[i] = cos(psi);
			VY[i] = sin(psi);
		}

	}


	int *block = new int[RadiusSystem*RadiusSystem + 1];
	int *lscl = new int[ParticleNum];
	int **cellb = new int*[RadiusSystem*RadiusSystem];


	for (int i = 0; i < RadiusSystem*RadiusSystem; i++)
	{
		cellb[i] = new int[9];
	}

	int cellv0[9][2] =
	{ { RadiusSystem - 1, RadiusSystem - 1 }, { 0, RadiusSystem - 1 }, { 1, RadiusSystem - 1 }, { RadiusSystem - 1, 0 }, { 0, 0 }, { 1, 0 }, { RadiusSystem - 1, 1 }, { 0, 1 }, { 1, 1 } };

	for (int b = 0; b < RadiusSystem*RadiusSystem; b++){
		for (int j = 0; j < 9; j++){
			int bx = cellv0[j][0] + b % RadiusSystem; int by = (cellv0[j][1] + b / RadiusSystem);
			cellb[b][j] = bx%RadiusSystem + (by%RadiusSystem)*RadiusSystem;

		}
	}

	

	while (fabs(L0 - L100) > eps){

		time_elapsed += dt;
		
		if (Number_of_steps%Accuracy == (Accuracy - 1)) {
			L100 = 0.;
			for (int i = 0; i < ParticleNum; i++){
				L100 += fabs((X[i] * VY[i] - Y[i] * VX[i]));
			}
		}

		if (Number_of_steps%Accuracy == 0){
			L0 = 0.;
			for (int i = 0; i < ParticleNum; i++){
				L0 += fabs((X[i] * VY[i] - Y[i] * VX[i]));
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


		for (int b0 = 0; b0 < RadiusSystem*RadiusSystem; b0++){
			int i = block[b0];
			while (i != -1){
				for (int b1 = 0; b1 < 9; b1++){
					int b2 = cellb[b0][b1];
					int j = block[b2];
					while (j != -1){
						dx = abs(X[i] - X[j]);
						dy = abs(Y[i] - Y[j]);
						if ((hypot(dx, dy) < RadiusInteraction) && j != i)
						{
							sx += OrtsVX[j]; sy += OrtsVY[j];
						}
						j = lscl[j];
					}
				}
				noise_action = noise_amp * M_PI *(1 - 2 * ((double)rand() / (RAND_MAX)));

				if (abs(sx) < 1e-5 && abs(sy) < 1e-5){
					Vx = OrtsVX[i];
					Vy = OrtsVY[i];
					OrtsVX[i] = Vx * cos(noise_action) - Vy * sin(noise_action);
					OrtsVY[i] = Vx * sin(noise_action) + Vy * cos(noise_action);
				}
				else {
					OrtsVX[i] = (sx * cos(noise_action) - sy * sin(noise_action)) / hypot(sx, sy);
					OrtsVY[i] = (sx * sin(noise_action) + sy * cos(noise_action)) / hypot(sx, sy);
				}

				sx = 0., sy = 0.;
				mod = hypot(VX[i], VY[i]);
				VX[i] = OrtsVX[i] * mod;
				VY[i] = OrtsVY[i] * mod;

				
				r = radius_step;
				for (int q = 0; q < TimesRadius; q++){
					if ((hypot(X[i], Y[i]) < r) && (hypot(X[i], Y[i]) > (r - radius_step))) {
						Amount_Particle[q] += (1 / (double)ParticleNum);
						MODULE_VELOCITY_HIST[q] += ((VX[i] * VX[i] + VY[i] * VY[i]) / (double)ParticleNum);
						MOMENTUM_VELOCITY_HIST[q] += fabs((X[i] * VY[i] - Y[i] * VX[i]) / ((X[i] * X[i] + Y[i] * Y[i]) * (double)ParticleNum));
						break;
					}
					r += radius_step;
				}

				double x1 = X[i] + VX[i] * dt; double y1 = Y[i] + VY[i] * dt;
				if (x1*x1 + y1*y1 > (double)(RadiusSystem*RadiusSystem)){
					SQRT = pow(((X[i] * VX[i] + Y[i] * VY[i]) / (VX[i] * VX[i] + VY[i] * VY[i])), 2) + (((double)(RadiusSystem*RadiusSystem) - Y[i] * Y[i] - X[i] * X[i]) / (VX[i] * VX[i] + VY[i] * VY[i]));
					dt1 = fabs(sqrt(fabs(SQRT)) - ((X[i] * VX[i] + Y[i] * VY[i]) / (VX[i] * VX[i] + VY[i] * VY[i])));

					x0 = X[i] + VX[i] * dt1;
					y0 = Y[i] + VY[i] * dt1;

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

				i = lscl[i];
			}
		}
		for (int q = 0; q < TimesRadius; q++){
			profiles << time_elapsed << '\t' << q
				<< '\t' << Amount_Particle[q] << '\t' << MODULE_VELOCITY_HIST[q]
				<< '\t' << MOMENTUM_VELOCITY_HIST[q] << '\n';
		}
	}
	delete[] block;
	delete[] lscl;
	delete[] Amount_Particle;
	delete[] MODULE_VELOCITY_HIST;
	delete[] MOMENTUM_VELOCITY_HIST;
}