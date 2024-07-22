// Galaxy Simulation by Connor Harris.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h> // used for timer

double plummer_sphere = 1e-3; // Plummer sphere as a part of the force calculation.
char X = 'x', Y = 'y';
double gravity;
double delta_time;
int time_steps;
int star_count;

// Structure for each individual particle
struct particles
{
    double *position_x, *position_y, *mass, *velocity_x, *velocity_y, *brightness, *accel_x, *accel_y;
};

// Function Declarations
int main(int argc, char const *argv[]);
int calculate_position(struct particles *all_particles);
int update_accelleration(struct particles *all_particles);
int update_position(struct particles *all_particles);
double force(struct particles *all_particles, char axis, int star_count, double gravity, int index_1);
double distance(struct particles *all_particles, int index_1, int index_2);
double vector(struct particles *all_particles, int index_1, int index_2, char axis);
double force_bidirectional(struct particles *all_particles, char axis, int particle_1, int particle_2);
void reset_array(struct particles *all_particles);

// For each user defined timestep, calculate all particles new position
int calculate_position(struct particles *all_particles)
{
    for (int current_timestep = 0; current_timestep < time_steps; current_timestep++)
    {
        update_accelleration(all_particles);
        update_position(all_particles);
    }
    return 1;
}

/*
    Calculate the new accelleration of all the particles.
        Accelleration = Force(i) / Mass(i)
        Force(i) = -gravity * mass(i) * sum(mass(j) / (distance(i, j) + plummer_sphere)^3 * vector(i, j))
*/
int update_accelleration(struct particles *all_particles)
{
    reset_array(all_particles);
    for (int index = 0; index < star_count; index++)
    {
        for (int i = index; i < star_count; i++)
        {
            double force_x = force_bidirectional(all_particles, X, index, i);
            double force_y = force_bidirectional(all_particles, Y, index, i);

            all_particles->accel_x[index] += (force_x * all_particles->mass[i]);
            all_particles->accel_y[index] += force_y * all_particles->mass[i];
            all_particles->accel_x[i] -= (force_x * all_particles->mass[index]);
            all_particles->accel_y[i] -= (force_y * all_particles->mass[index]);
        }
        all_particles->accel_x[index] = (-gravity * all_particles->mass[index] * all_particles->accel_x[index]) / all_particles->mass[index];
        all_particles->accel_y[index] = (-gravity * all_particles->mass[index] * all_particles->accel_y[index]) / all_particles->mass[index];
    }
    return 1;
}

int update_position(struct particles *all_particles)
{

    for (int index = 0; index < star_count; index++)
    {
        // Calculate velocity on all particles with equation: Velocity(i) + Delta time * Acceleration(i)
        all_particles->velocity_x[index] += (delta_time * all_particles->accel_x[index]); // X direction
        all_particles->velocity_y[index] += (delta_time * all_particles->accel_y[index]); // Y direction

        // Update position of particles with equation: Position(i) + Delta Time * Velocity(i)
        all_particles->position_x[index] += (delta_time * all_particles->velocity_x[index]); // X directions
        all_particles->position_y[index] += (delta_time * all_particles->velocity_y[index]); // Y directions
    }
    return 1;
}

inline void reset_array(struct particles *all_particles)
{ // set the acceleration array back to zero
    for (int i = 0; i < star_count; i++)
    {
        all_particles->accel_x[i] = 0;
        all_particles->accel_y[i] = 0;
    }
}

// Calculate the reciprocal distance * vector(i, j) for two particles
inline double force_bidirectional(struct particles *all_particles, char axis, int particle_1, int particle_2)
{
    // return (dist(i, j) + e3)^3) * vector(i, j)
    double dist = distance(all_particles, particle_1, particle_2) + plummer_sphere;
    double dist_reciprocal = 1.0 / (dist * dist * dist);
    return (dist_reciprocal * vector(all_particles, particle_1, particle_2, axis));
}

// Calculate the distance of two points. Distance = Sqrt((x2 - x1)^2 + (y2-y1)^2)
inline double distance(struct particles *all_particles, int index_1, int index_2)
{
    double x_distance = all_particles->position_x[index_2] - all_particles->position_x[index_1];
    double y_distance = all_particles->position_y[index_2] - all_particles->position_y[index_1];
    return sqrt(x_distance * x_distance + y_distance * y_distance);
}

// Calculate the vector of point A relative to point B in its respective axis. Vector = (x1 - x2, y1 - y2)
inline double vector(struct particles *all_particles, int index_1, int index_2, char axis)
{
    return (axis == X ? (all_particles->position_x[index_1] - all_particles->position_x[index_2])
                      : all_particles->position_y[index_1] - all_particles->position_y[index_2]);
}

// Wall seconds function copied from HPPC labs
double get_wall_seconds()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + (double)tv.tv_usec / 1000000;
}

// Main function
int main(int argc, char const *argv[])
{
    // Check there are 5 parameters, but assume the input is correct (program, star_count, filename, timesteps, delta_t, graphics)
    if (argc != 6)
    {
        printf("Error: File requires 5 input parameters\n");
        return 1;
    }

    // Assign user input to variables, and declare other variables used.
    star_count = atoi(argv[1]);
    char *file_name = malloc(strlen(argv[2]) + 1);
    strcpy(file_name, argv[2]);
    FILE *file = fopen(file_name, "rb");
    time_steps = atoi(argv[3]);
    delta_time = atof(argv[4]);
    // int graphics = atoi(argv[5]); // Commented out when graphics is not used
    gravity = (100.0 / star_count);
    double start_overall_timer = get_wall_seconds(); // Timer for measuring entire process

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
    for (int index = 0; index < star_count; index++)
    {
        if (fread(&all_particles.position_x[index], sizeof(double), 1, file) != 1 ||
            fread(&all_particles.position_y[index], sizeof(double), 1, file) != 1 ||
            fread(&all_particles.mass[index], sizeof(double), 1, file) != 1 ||
            fread(&all_particles.velocity_x[index], sizeof(double), 1, file) != 1 ||
            fread(&all_particles.velocity_y[index], sizeof(double), 1, file) != 1 ||
            fread(&all_particles.brightness[index], sizeof(double), 1, file) != 1)
        {
            fprintf(stderr, "Error reading file\n");
            break;
        }
    }

    // call the function to calculate new positions
    calculate_position(&all_particles);

    // Declare and populate the output file with the results
    FILE *output = fopen("result.gal", "wb");
    for (int output_line = 0; output_line < star_count; output_line++)
    {
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
    start_overall_timer = get_wall_seconds() - start_overall_timer; // Stop timing process for overall runtime
    printf("\nOverall time: %lf\n\n", start_overall_timer);
    return 0;
}
