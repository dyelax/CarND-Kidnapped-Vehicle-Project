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
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

  num_particles = 200;
  particles.resize(num_particles);
  weights.resize(num_particles, 1);

  float std_x = std[0];
  float std_y = std[1];
  float std_theta = std[2];
  
  // Create normal distributions for x, y and theta
  default_random_engine gen;
  normal_distribution<double> dist_x(x, std_x);
  normal_distribution<double> dist_y(y, std_y);
  normal_distribution<double> dist_theta(theta, std_theta);
  
  for (int i = 0; i < num_particles; i++) {
    float sample_x = dist_x(gen);
    float sample_y = dist_y(gen);
    float sample_theta = dist_theta(gen);
    
    // Create new particle
    Particle p;
    p.id = i;
    p.x = sample_x;
    p.y = sample_y;
    p.theta = sample_theta;
    p.weight = weights[i];
    
    particles[i] = p;
  }
  
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
  
  float std_x = std_pos[0];
  float std_y = std_pos[1];
  float std_theta = std_pos[2];
  
  default_random_engine gen;
  for (int i = 0; i < num_particles; i++) {
    Particle p = particles[i];
    
    // Add gaussian noise
    normal_distribution<double> dist_x(p.x, std_x);
    normal_distribution<double> dist_y(p.y, std_y);
    normal_distribution<double> dist_theta(p.theta, std_theta);
    
    float noisy_x = dist_x(gen);
    float noisy_y = dist_y(gen);
    float noisy_theta = dist_theta(gen);
    
    // Compute predicted location
    float new_x, new_y, new_theta;
    if (abs(yaw_rate) < 0.0000001) {
      new_x = noisy_x + velocity * delta_t * cos(new_theta);
      new_y = noisy_y + velocity * delta_t * sin(new_theta);
      new_theta = noisy_theta;
    } else {
      new_x =
        noisy_x +
        (velocity / yaw_rate) * (sin(noisy_theta + (yaw_rate * delta_t)) -
                                 sin(noisy_theta));
      
      new_y =
        noisy_y +
        (velocity / yaw_rate) * (cos(noisy_theta) -
                                 cos(noisy_theta + (yaw_rate * delta_t)));
      
      new_theta = noisy_theta + yaw_rate * delta_t;
    }
    
    // Create new particle with predicted position
    // TODO: check if we need to do this, or can just modify p.
    Particle new_p;
    new_p.x = new_x;
    new_p.y = new_y;
    new_p.theta = new_theta;
    new_p.associations = p.associations;
    new_p.sense_x = p.sense_x;
    new_p.sense_y = p.sense_y;
    
    particles[i] = new_p;
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

  
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
  float std_x = std_landmark[0];
  float std_y = std_landmark[1];
  float denominator = 2 * M_PI * std_x * std_y;
  
  for (int p_i = 0; p_i < num_particles; p_i++) {
    Particle p = particles[p_i];
    
    float new_weight = 1.;
    for (int o_i = 0; o_i < observations.size(); o_i++) {
      LandmarkObs o = observations[o_i];
      
      // Transform observation to map coordinates
      float xt = o.x * cos(p.theta) - o.y * sin(p.theta) + p.x;
      float yt = o.x * sin(p.theta) + o.y * cos(p.theta) + p.y;
      
      // Associate transformed observation with nearest landmark
      float min_dist = sensor_range;
      float x_diff = sensor_range;
      float y_diff = sensor_range;
      for (int l_i = 0; l_i < map_landmarks.landmark_list.size(); l_i++) {
        Map::single_landmark_s landmark = map_landmarks.landmark_list[l_i];
        
        float ol_dist = dist(xt, yt, landmark.x_f, landmark.y_f);
        if (ol_dist < min_dist) {
          min_dist = ol_dist;
          
          x_diff = xt - landmark.x_f;
          y_diff = yt - landmark.y_f;
        }
      }
      
      // Weight particle based on multivariate gaussian of obs vs associated landmark
      float power = -((pow(x_diff, 2) / (2 * pow(std_x, 2))) +
                      (pow(y_diff, 2) / (2 * pow(std_y, 2))));
      new_weight *= exp(power) / denominator;
    }
    weights[p_i] = new_weight;
    particles[p_i].weight = new_weight;
  }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  vector<Particle> new_p(num_particles);

  random_device rd;
  mt19937 gen(rd());
  for (int i = 0; i < num_particles; i++) {
    // Sample a particle
    discrete_distribution<> dist(weights.begin(), weights.end());
    int sample_i = dist(gen);
    new_p[i] = particles[sample_i];
  }
  
  particles = new_p;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
