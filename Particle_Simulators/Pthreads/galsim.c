// 2d Particle Simulation Assignment by Connor Harris. Pthreads version
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <pthread.h>

double plummer_sphere = 1e-3; // Plummer sphere as a part of the force calculation.
char X = 'x', Y = 'y';
double gravity, delta_time;
int time_steps, star_count;

pthread_mutex_t mutex;

struct particles
{
    double *position_x, *position_y, *mass, *velocity_x, *velocity_y, *brightness, *accel_x, *accel_y;
};

struct thread_data
{
    int start_index_load_balanced, end_index_load_balanced, start_index, end_index;
    struct particles *all_particles;
    pthread_barrier_t *barrier;
};

int main(int argc, char const *argv[]);
void load_balance_threads(int *thread_ranges_balanced, int length);
void *calculate_particles(void *threadarg);
void reset_array(struct particles *all_particles, int start_index_load_balanced, int end_index_load_balanced);
double force_bidirectional(struct particles *all_particles, char axis, int particle_1, int particle_2);
double distance(struct particles *all_particles, int current_particle, int next_particle);
double vector(struct particles *all_particles, int current_particle, int next_particle, char axis);

/*
    Higher particle numbers require less work due to the symmetrical force property,
    so this function calculates how many particles each thread should take relative to
    the work required for calcualting the updated accelleration.
*/
inline void load_balance_threads(int *thread_ranges_balanced, int length)
{

    int symmetrical_calculations = (star_count * (star_count + 1) / 2); // sum of natural numbers formula
    int calculations_per_thread = symmetrical_calculations / (length - 1);
    thread_ranges_balanced[0] = 0;
    int increment = 0;
    int decrement = (length - 1);
    int range = 0;

    for (int i = 0; i < star_count; i++)
    {
        range += i;
        if (range >= (calculations_per_thread * increment))
        {
            thread_ranges_balanced[decrement] = (star_count - i);
            increment++;
            decrement--;
            if (increment == (length - 1))
            {
                break;
            }
        }
    }
}

/*
    Symmetrically calculate a particles force and use it to update the accelleration,
    velocity and position. The symmetric force calculations are balanced out so that all threads
    do the same ammount of work, and then they use an even subset of particles
    for the asymmetric calculations.
*/
void *calculate_particles(void *threadarg)
{
    struct thread_data *data;
    data = (struct thread_data *)threadarg;
    int start_index_load_balanced = data->start_index_load_balanced;
    int end_index_load_balanced = data->end_index_load_balanced;
    pthread_barrier_t *barrier = data->barrier;
    int start_index = data->start_index;
    int end_index = data->end_index;

    for (int current_timestep = 0; current_timestep < time_steps; current_timestep++)
    {
        reset_array(data->all_particles, start_index_load_balanced, end_index_load_balanced);
        pthread_barrier_wait(barrier);
        double force_x, force_y, x_i, y_i;
        double local_x_j[star_count];
        double local_y_j[star_count];

        for (int thread_index = start_index_load_balanced; thread_index < end_index_load_balanced; thread_index++)
        {
            x_i = 0.0;
            y_i = 0.0;
            for (int i = thread_index; i < star_count; i++)
            {
                force_x = force_bidirectional(data->all_particles, X, thread_index, i);
                force_y = force_bidirectional(data->all_particles, Y, thread_index, i);
                x_i += (force_x * data->all_particles->mass[i]);
                y_i += (force_y * data->all_particles->mass[i]);
                local_x_j[i] -= (force_x * data->all_particles->mass[thread_index]);
                local_y_j[i] -= (force_y * data->all_particles->mass[thread_index]);
            }

            pthread_mutex_lock(&mutex);
            data->all_particles->accel_x[thread_index] += x_i;
            data->all_particles->accel_y[thread_index] += y_i;

            for (int i = thread_index; i < star_count; i++)
            {
                data->all_particles->accel_x[i] += local_x_j[i];
                data->all_particles->accel_y[i] += local_y_j[i];
                local_x_j[i] = 0.0;
                local_y_j[i] = 0.0;
            }
            pthread_mutex_unlock(&mutex);
        }
        pthread_barrier_wait(barrier);

        for (int thread_index = start_index; thread_index < end_index; thread_index++)
        {
            data->all_particles->accel_x[thread_index] = (-gravity * data->all_particles->mass[thread_index] * data->all_particles->accel_x[thread_index]) / data->all_particles->mass[thread_index];
            data->all_particles->accel_y[thread_index] = (-gravity * data->all_particles->mass[thread_index] * data->all_particles->accel_y[thread_index]) / data->all_particles->mass[thread_index];
        }

        // Calculate a particle segments velocity, then position
        for (int thread_index = start_index; thread_index < end_index; thread_index++)
        {
            data->all_particles->velocity_x[thread_index] += (delta_time * data->all_particles->accel_x[thread_index]);
            data->all_particles->velocity_y[thread_index] += (delta_time * data->all_particles->accel_y[thread_index]);
            data->all_particles->position_x[thread_index] += (delta_time * data->all_particles->velocity_x[thread_index]);
            data->all_particles->position_y[thread_index] += (delta_time * data->all_particles->velocity_y[thread_index]);
        }
        pthread_barrier_wait(barrier);
    }
    pthread_exit(NULL);
    return 0;
}

// Reset the accel values of a given particle subset to zero
void reset_array(struct particles *all_particles, int start_index_load_balanced, int end_index_load_balanced)
{
    for (int i = start_index_load_balanced; i < end_index_load_balanced; i++)
    {
        all_particles->accel_x[i] = 0;
        all_particles->accel_y[i] = 0;
    }
}

// Calculate the reciprocal (dist(i, j) + e3)^3) * vector(i, j) for two particles
inline double force_bidirectional(struct particles *all_particles, char axis, int particle_1, int particle_2)
{
    double dist = distance(all_particles, particle_1, particle_2) + plummer_sphere;
    double dist_reciprocal = 1.0 / (dist * dist * dist);
    return (dist_reciprocal * vector(all_particles, particle_1, particle_2, axis));
}

// Calculate the distance of two points. Distance = Sqrt((x2 - x1)^2 + (y2-y1)^2)
inline double distance(struct particles *all_particles, int current_particle, int next_particle)
{
    double x_distance = all_particles->position_x[next_particle] - all_particles->position_x[current_particle];
    double y_distance = all_particles->position_y[next_particle] - all_particles->position_y[current_particle];
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
    // Check there are 6 parameters, but assume the input is correct (program, star_count, filename, timesteps, delta_t, graphics)
    if (argc != 7)
    {
        printf("Error: File requires 6 input parameters. Try: \n");
        printf("./galsim 1000 input_data/ellipse_N_01000.gal 200 1e-5 0 1\n");
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
    int thread_count = atoi(argv[6]);
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

    // Pthread data declarations
    pthread_t threads[thread_count];
    struct thread_data thread_data_array[thread_count];
    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier, NULL, thread_count);
    pthread_mutex_init(&mutex, NULL);

    // Particle subsets with equal division for asymmetrical calculations, and distributed for symmetrical calculations
    int particles_subset_length = star_count / thread_count;
    int *thread_ranges_balanced = (int *)malloc((thread_count + 1) * sizeof(int));
    load_balance_threads(thread_ranges_balanced, thread_count + 1);

    // TODO: Clean this loop
    for (int current_thread = 0; current_thread < thread_count; current_thread++)
    {
        thread_data_array[current_thread].start_index_load_balanced = thread_ranges_balanced[current_thread];
        thread_data_array[current_thread].end_index_load_balanced = thread_ranges_balanced[current_thread + 1];
        thread_data_array[current_thread].start_index = current_thread * particles_subset_length;
        thread_data_array[current_thread].end_index = (current_thread == thread_count - 1) ? star_count : (current_thread + 1) * particles_subset_length;
        thread_data_array[current_thread].all_particles = &all_particles;
        thread_data_array[current_thread].barrier = &barrier;

        int thread_creation_unsuccessful = pthread_create(&threads[current_thread], NULL, calculate_particles, (void *)&thread_data_array[current_thread]);
        if (thread_creation_unsuccessful)
        {
            printf("Error: return code from pthread_create() is %d\n", thread_creation_unsuccessful);
            exit(1);
        }
    }

    // Join all threads to wait for them to finish
    for (int t = 0; t < thread_count; t++)
    {
        pthread_join(threads[t], NULL);
    }

    pthread_mutex_destroy(&mutex);

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
    printf("\nOverall time: %lf\n", start_overall_timer);
    return 0;
}
