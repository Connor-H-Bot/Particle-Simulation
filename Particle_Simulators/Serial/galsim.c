// 2D particle simulation by Connor Harris
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double plummer_sphere = 1e-3; // Plummer sphere as a part of the force calculation. 

// Structure for each individual particle
struct particles {
    double *position_x, *position_y, *mass, *velocity_x, *velocity_y, *brightness, *accel_x, *accel_y; 
};

// Function Declarations
int main(int argc, char const *argv[]);
double force(struct particles* all_particles, char axis, int star_count, double gravity, int index_1);
double distance(struct particles* all_particles, int index_1, int index_2);
double vector(struct particles* all_particles, int index_1, int index_2, char axis);

// Calculate force on a particle. Equation: Force(i) = -gravity * mass(i) * sum(mass(j) / (distance(i, j) + plummer_sphere)^3 * vector(i, j))
double force(struct particles* all_particles, char axis, int star_count, double gravity, int index_1) {
    double summation = 0.0;
    double dist;
    for (int i = 0; i < star_count; i++) { 
        if (index_1 != i) {
            dist = distance(all_particles, index_1, i) + plummer_sphere;
            summation += ((all_particles->mass[i]) / (dist * dist * dist)) * vector(all_particles, index_1, i, axis);
        }
    }
    return (-gravity * all_particles->mass[index_1] * summation);
}

// Calculate the distance of two points. Distance = Sqrt((x2 - x1)^2 + (y2-y1)^2)
double distance(struct particles* all_particles, int index_1, int index_2) {
    return sqrt(pow((all_particles->position_x[index_2] - all_particles->position_x[index_1]), 2.0) + pow((all_particles->position_y[index_2] - all_particles->position_y[index_1]), 2.0));
}

// Calculate the vector of point A relative to point B in its respective axis. Vector = (x1 - x2, y1 - y2)
double vector(struct particles* all_particles, int index_1, int index_2, char axis){
    if (axis == 'x') {
        return (all_particles->position_x[index_1] - all_particles->position_x[index_2]);
    }
    else {
        return (all_particles->position_y[index_1] - all_particles->position_y[index_2]);
    }
}

// Main function
int main(int argc, char const *argv[]) { 
    // Check there are 5 parameters, but assume the input is correct (program, star_count, filename, timesteps, delta_t, graphics)
    if (argc != 6) {
        printf("Error: File requires 5 input parameters");
        return 1;
    }

    // Assign user input to variables, and declare other variables used. 
    int star_count = atoi(argv[1]); 
    char *file_name = malloc(strlen(argv[2]) + 1);
    strcpy(file_name, argv[2]);
    FILE *file = fopen(file_name, "rb");
    int time_steps = atoi(argv[3]);
    double delta_time = atof(argv[4]);
    // int graphics = atoi(argv[5]); // Commented out when graphics is not used
    double gravity = (100.0/star_count);
    char X = 'x', Y = 'y';

    // Allocate memory for the structure of arrays.
    struct particles all_particles;
    all_particles.position_x = malloc(star_count * sizeof(double));
    all_particles.position_y = malloc(star_count * sizeof(double));
    all_particles.mass = malloc(star_count * sizeof(double));
    all_particles.velocity_x = malloc(star_count * sizeof(double));
    all_particles.velocity_y = malloc(star_count * sizeof(double));
    all_particles.brightness = malloc(star_count * sizeof(double));
    all_particles.accel_x = malloc(star_count * sizeof(double));
    all_particles.accel_y = malloc(star_count * sizeof(double));

    // Open the input file and use it to populate the particles in the array of structures
    for (int index = 0; index < star_count; index++) {
        if (fread(&all_particles.position_x[index], sizeof(double), 1, file) != 1
            || fread(&all_particles.position_y[index], sizeof(double), 1, file) != 1
            || fread(&all_particles.mass[index], sizeof(double), 1, file) != 1
            || fread(&all_particles.velocity_x[index], sizeof(double), 1, file) != 1
            || fread(&all_particles.velocity_y[index], sizeof(double), 1, file) != 1
            || fread(&all_particles.brightness[index], sizeof(double), 1, file) != 1) {
            fprintf(stderr, "Error reading file\n");
            break;
        }
    }

    // Main loop that runs for the ammount of time steps. Inside the loop are accelleration, velocity, and updated position calculations. 
    for (int step_iteration = 0; step_iteration < time_steps; step_iteration++) { 

        for (int index = 0; index < star_count; index++) { 
            all_particles.accel_x[index] = (force(&all_particles, X, star_count, gravity, index) / all_particles.mass[index]); 
            all_particles.accel_y[index] = (force(&all_particles, Y, star_count, gravity, index) / all_particles.mass[index]);
        }

        for (int index = 0; index < star_count; index++) { // Calculate Velocity and then position

            // Calculate velocity on all particles with equation: Velocity(i) + Delta time * Acceleration(i)
            all_particles.velocity_x[index] += (delta_time * all_particles.accel_x[index]); // X direction
            all_particles.velocity_y[index] += (delta_time * all_particles.accel_y[index]); // Y direction

            // Update position of particles with equation: Position(i) + Delta Time * Velocity(i)
            all_particles.position_x[index] += (delta_time * all_particles.velocity_x[index]); // X directions
            all_particles.position_y[index] += (delta_time * all_particles.velocity_y[index]); // Y directions
        }
    }

    // Declare and populate the output file with the results
    FILE *output = fopen("result.gal", "wb");
    for (int output_line = 0; output_line < star_count; output_line++) {
        double position_x = all_particles.position_x[output_line]; 
        double position_y = all_particles.position_y[output_line]; 
        double mass = all_particles.mass[output_line]; 
        double velocity_x = all_particles.velocity_x[output_line]; 
        double velocity_y = all_particles.velocity_y[output_line]; 
        double brightness = all_particles.brightness[output_line];
        fwrite(&position_x, sizeof(double), 1, output);
        fwrite(&position_y, sizeof(double), 1, output);
        fwrite(&mass, sizeof(double), 1, output);
        fwrite(&velocity_x, sizeof(double), 1, output);
        fwrite(&velocity_y, sizeof(double), 1, output);
        fwrite(&brightness, sizeof(double), 1, output);
    }

    fclose(output);
    free(file_name);
    free(all_particles.accel_x);
    free(all_particles.accel_y);
    free(all_particles.brightness);
    free(all_particles.mass);
    free(all_particles.position_x);
    free(all_particles.position_y);
    free(all_particles.velocity_x);
    free(all_particles.velocity_y);
    return 0;
}