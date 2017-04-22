/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "particle_filter.h"

 using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 50;
	default_random_engine gen;
	normal_distribution<double> N_x(x,std[0]);
	normal_distribution<double> N_y(y,std[1]);
	normal_distribution<double> N_theta(theta,std[2]);

	for(int i = 0 ; i < num_particles ; i++){
		Particle particle;
		particle.id = i;
		particle.x = N_x(gen);
		particle.y = N_y(gen);
		particle.theta = N_theta(gen);
		particle.weight = 1;
		particles.push_back(particle);
		weights.push_back(1);
	}
	is_initialized = true;





}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;
	double new_x,new_y,new_theta;

	for(int i = 0 ; i < num_particles ; i++){
		if(yaw_rate == 0){
			new_x = particles[i].x + velocity * delta_t * cos(particles[i].theta);
			new_y = particles[i].y + velocity * delta_t * sin(particles[i].theta);
			new_theta = particles[i].theta;
		}else{
			new_x = particles[i].x + velocity/yaw_rate*(sin(particles[i].theta+yaw_rate*delta_t)-sin(particles[i].theta));
			new_y = particles[i].y + velocity/yaw_rate*(cos(particles[i].theta)-cos(particles[i].theta+yaw_rate*delta_t));
			new_theta = particles[i].theta + yaw_rate*delta_t;
		}

		normal_distribution<double> N_x(new_x,std_pos[0]);
		normal_distribution<double> N_y(new_y,std_pos[1]);
		normal_distribution<double> N_theta(new_theta,std_pos[2]);

		particles[i].x = N_x(gen);
		particles[i].y = N_y(gen);
		particles[i].theta = N_theta(gen);
	}

}

void ParticleFilter::dataAssociation(vector<LandmarkObs>& predicted, vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.


	//for(int p=0;p<num_particles;p++){
	//	for(int i=0;i<observations.size();i++){
	//		LandmarkObs trans_obs;
	//		trans_obs.x = particles[p].x + observations[i].x * cos(particles[p].theta) - observations[i].y * sin(particles[i].theta);
	//		trans_obs.y = particles[p].y + observations[i].x * sin(particles[p].theta) + observations[i].y * cos(particles[i].theta);
	//		predicted.push_back(trans_obs);
	//	}
	//}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html

	for (int p=0;p<num_particles;p++){
		vector<LandmarkObs> trans_observations;

		for(int i=0; i<observations.size();i++){
			//calculate each of the components
			LandmarkObs trans_obs;
			trans_obs.x = particles[p].x + observations[i].x * cos(particles[p].theta) - observations[i].y * sin(particles[i].theta);
			trans_obs.y = particles[p].y + observations[i].x * sin(particles[p].theta) + observations[i].y * cos(particles[i].theta);
			//add to the vector
			trans_observations.push_back(trans_obs);
		}
		particles[p].weight = 1.0;

		//find the closest landmark
		for(int i=0; i<trans_observations.size();i++){
			double closest = sensor_range;
			int association = 0;
			for (int j=0;j<map_landmarks.landmark_list.size();j++){
				double landmark_x = map_landmarks.landmark_list[j].x_f;
				double landmark_y = map_landmarks.landmark_list[j].y_f;

				//use function from helper.h to calculate distance
				double distance = dist(trans_observations[i].x,trans_observations[i].y,landmark_x,landmark_y);
				//update the closest landmark
				if(distance<closest){
					closest = distance;
					association = j;
				}
			}
			if(association!=0){
				double meas_x = trans_observations[i].x;
				double meas_y = trans_observations[i].y;
				double mu_x = map_landmarks.landmark_list[association].x_f;
				double mu_y = map_landmarks.landmark_list[association].y_f;
				long double mult_variate = 1.0/(2.0*M_PI*std_landmark[0]*std_landmark[1])*exp(-((pow(meas_x-mu_x,2.0)/(2*pow(std_landmark[0],2.0)))+(pow(meas_y-mu_y,2.0))/(2*pow(std_landmark[1],2.0))));

				//cout<<mult_variate<<endl;

				if(mult_variate != 0){
				particles[p].weight *= mult_variate;
				}

			}
		}
		weights[p] = particles[p].weight;

	}


	//for(int i = 0; i<observations.size();i++){
	//std::cout<< observations[i].x<<std::endl;
	//}

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	default_random_engine gen;
	discrete_distribution<int> distribution(weights.begin(),weights.end());

	vector<Particle> resample_particles;

	for(int i=0;i<num_particles;i++){
		resample_particles.push_back(particles[distribution(gen)]);
	}
	particles = resample_particles;



}

void ParticleFilter::write(string filename) {
	// You don't need to modify this file.
	ofstream dataFile;
	dataFile.open(filename, ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
