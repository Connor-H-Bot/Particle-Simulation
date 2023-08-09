// 2D particle simulation by Connor Harris. 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>

double plummer_sphere = 1e-3; // Plummer sphere as a part of the force calculation.

// Structure for each individual particle
struct particles {
    double *position_x, *position_y, *mass, *velocity_x, *velocity_y, *brightness, *accel_x, *accel_y; 
};

// Function Declarations
int main(int argc, char const *argv[]);
double force(int idx, char axis, struct particles* all_particles, int star_count, double gravity);
double distance(int idx1, int idx2, struct particles* all_particles);
double vector(int idx1, int idx2, char axis, struct particles* all_particles);

// Calculate force on a particle. Equation: Force(i) = -gravity * mass(i) * sum(mass(j) / (distance(i, j) + plummer_sphere)^3 * vector(i, j))
double force(int idx, char axis, struct particles* all_particles, int star_count, double gravity) {
    double summation = 0.0;
    for (int i = 0; i < star_count; i++) { 
        if (idx != i) {
            summation += ((all_particles->mass[i]) / pow((distance(idx, i, all_particles) + plummer_sphere), 3.0)) * vector(idx, i, axis, all_particles);
        }
    }
    return (-gravity * all_particles->mass[idx] * summation);
}


// Calculate the distance of two points. Distance = Sqrt((x2 - x1)^2 + (y2-y1)^2)
double distance(int idx1, int idx2, struct particles* all_particles) {
    return sqrt(pow((all_particles->position_x[idx2] - all_particles->position_x[idx1]), 2.0) + pow((all_particles->position_y[idx2] - all_particles->position_y[idx1]), 2.0));
}

// Calculate the vector of point A relative to point B in its respective axis. Vector = (x1 - x2, y1 - y2)
double vector(int idx1, int idx2, char axis, struct particles* all_particles){
    if (axis == 'x') {
        return (all_particles->position_x[idx1] - all_particles->position_x[idx2]);
    }
    else {
        return (all_particles->position_y[idx1] - all_particles->position_y[idx2]);
    }
}

// Wall time for timing func
double get_wall_seconds() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
    return seconds;
}

// Main function
int main(int argc, char const *argv[]) { 
    // Check there are 6 parameters, but assume the input is correct (program, star_count, filename, timesteps, delta_t, graphics)
    if (argc != 7) {
        printf("Error: File requires 5 input parameters");
        return 1;
    }

    double start_overall_timer = get_wall_seconds(); // Start timer before I/O

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
    int desired_threads = atoi(argv[6]);
    omp_set_num_threads(desired_threads); // Set the number of threads for OpenMP

    // Allocation of memory for particles
    struct particles all_particles;
    all_particles.position_x = (double *) malloc(star_count * sizeof(double));
    all_particles.position_y = (double *) malloc(star_count * sizeof(double));
    all_particles.mass = (double *) malloc(star_count * sizeof(double));
    all_particles.velocity_x = (double *) malloc(star_count * sizeof(double));
    all_particles.velocity_y = (double *) malloc(star_count * sizeof(double));
    all_particles.brightness = (double *) malloc(star_count * sizeof(double));
    all_particles.accel_x = (double *) malloc(star_count * sizeof(double));
    all_particles.accel_y = (double *) malloc(star_count * sizeof(double));

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
        #pragma omp parallel // Parallelised section
        {
            int thread_id = omp_get_thread_num(); // Assign a thread its unique ID from 1-N
            int total_threads = omp_get_num_threads(); // Variable to hold desired_threads indexed with 0
            int chunk_size = star_count / total_threads; // Divides the star (particle) count by the threads 
            int start_index = thread_id * chunk_size; // Assigns a start index for this unqie threads chunk
            int end_index = (thread_id == total_threads - 1) ? star_count : (thread_id + 1) * chunk_size; // and its end (accounts for a remainder)

            // Calculate acceleration on all particles with equation: Force(i) / Mass(i)
            for (int index = start_index; index < end_index; index++) {
                all_particles.accel_x[index] = (force(index, X, &all_particles, star_count, gravity) / all_particles.mass[index]);
                all_particles.accel_y[index] = (force(index, Y, &all_particles, star_count, gravity) / all_particles.mass[index]);
            }

            // Ensure all acceleration calculations are done
            #pragma omp barrier 

            // Calculate Velocity and then position
            for (int index = start_index; index < end_index; index++) {
                // Calculate velocity on all particles with equation: Velocity(i) + Delta time * Acceleration(i)
                all_particles.velocity_x[index] += (delta_time * all_particles.accel_x[index]); // X direction
                all_particles.velocity_y[index] += (delta_time * all_particles.accel_y[index]); // Y direction

                // Update position of particles with equation: Position(i) + Delta Time * Velocity(i)
                all_particles.position_x[index] += (delta_time * all_particles.velocity_x[index]); // X directions
                all_particles.position_y[index] += (delta_time * all_particles.velocity_y[index]); // Y direction
            }
            // Ensure all position updates are done
            #pragma omp barrier 
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

    // Close files, free memory, gracefully terminate. 
    fclose(output);
    free(file_name);
    free(all_particles.position_x);
    free(all_particles.position_y);
    free(all_particles.mass);
    free(all_particles.velocity_x);
    free(all_particles.velocity_y);
    free(all_particles.brightness);
    free(all_particles.accel_x);
    free(all_particles.accel_y);
    //start_overall_timer = get_wall_seconds() - start_overall_timer; // Timer function. Comment out when not using results
    //printf("\nOverall time: %lf\n\n", start_overall_timer); 
    return 0;
}