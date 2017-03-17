// SelfParallel.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <math.h>
#include <thread>
#include "evolution.h"

#define CORENUM 4 //Количество потоков

int main(int argc, _TCHAR* argv[])
{
	std::thread Threads[CORENUM];
	
	int Number_of_particles = 5000; //Кол-во частиц 
	double ParticlesPerUnit = 3, //Плотность 
		noise_amp = 0.; //Амплитуда шума
	int RadiusSystem = int(sqrt(Number_of_particles / ParticlesPerUnit)) + 1;

	for (int ThreadNumber = 0; ThreadNumber < CORENUM; ThreadNumber++){
		Threads[ThreadNumber] = std::thread(Evolution, Number_of_particles,RadiusSystem,ThreadNumber,noise_amp);
	}
	for (int ThreadNumber = 0; ThreadNumber < CORENUM; ThreadNumber++){
		Threads[ThreadNumber].join();
	}
	system("PAUSE");
	return 0;
}

