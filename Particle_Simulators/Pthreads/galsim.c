// 2D particle simulation by Connor Harris.  
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <sys/time.h>

// Structure to hold data for each thread
struct thread_data {
    int start_index;
    int end_index;
    struct particle *all_particles;
    int star_count;
    double delta_time;
    double gravity;
    pthread_barrier_t *barrier;
};

// Structure for each individual particle
struct particle {
    double position_x, position_y, mass, velocity_x, velocity_y, brightness, accel_x, accel_y; 
};

// Globals
double plummer_sphere = 1e-3;  
int time_steps = 0;

// Function Prototypes
int main(int argc, char const *argv[]);
double force(struct particle* pointer_1, char axis, struct particle* pointer_2, int star_count, double gravity);
double distance(struct particle* pointer_1, struct particle* pointer_2);
double vector(struct particle* pointer_1, struct particle* pointer_2, char axis);

/* 
 * Calculate force on a particle based on gravitational pull.
 * 
 * Pre: `pointer_1` and `pointer_2` are valid pointers to particle structures.
 *      `axis` is either 'x' or 'y' representing the dimension.
 *      `star_count` is the number of stars in the simulation.
 *      `gravity` is a positive value for gravitational constant.
 * Post: Returns the calculated force on the particle.
 */
double force(struct particle* pointer_1, char axis, struct particle* pointer_2, int star_count, double gravity) {
    double summation = 0.0;
    for (int i = 0; i < star_count; i++) { 
        if (pointer_1 != pointer_2) {
            summation += (pointer_2->mass / pow(distance(pointer_1, pointer_2) + plummer_sphere, 3.0)) * vector(pointer_1, pointer_2, axis);
        }
        pointer_2++;
    }
    return (-gravity * pointer_1->mass * summation);
}

/* 
 * Calculate the Euclidean distance between two particles.
 * 
 * Pre: `pointer_1` and `pointer_2` are valid pointers to particle structures.
 * Post: Returns the distance between the two particles.
 */
double distance(struct particle* pointer_1, struct particle* pointer_2) {
    return sqrt(pow(pointer_2->position_x - pointer_1->position_x, 2.0) + pow(pointer_2->position_y - pointer_1->position_y, 2.0));
}

/* 
 * Calculate the vector between two particles in the specified dimension.
 * 
 * Pre: `pointer_1` and `pointer_2` are valid pointers to particle structures.
 *      `axis` is either 'x' or 'y' representing the dimension.
 * Post: Returns the vector value for the specified dimension.
 */
double vector(struct particle* pointer_1, struct particle* pointer_2, char axis){
    if (axis == 'x') {
        return (pointer_1->position_x - pointer_2->position_x);
    }
    else {
        return (pointer_1->position_y - pointer_2->position_y);
    }
}

/* 
 * Thread function to calculate particle properties.
 * 
 * Pre: `threadarg` is a valid pointer to `thread_data` structure.
 * Post: Updates the particle properties after simulating their motion for given time steps.
 */
void *calculate_particles(void *threadarg) {
    struct thread_data *data;
    data = (struct thread_data *) threadarg;
    
    for (int step_iteration = 0; step_iteration < time_steps; step_iteration++) {
        for (int index = data->start_index; index < data->end_index; index++) {
            struct particle *current_particle = &data->all_particles[index];
            current_particle->accel_x = force(current_particle, 'x', data->all_particles, data->star_count, data->gravity) / current_particle->mass;
            current_particle->accel_y = force(current_particle, 'y', data->all_particles, data->star_count, data->gravity) / current_particle->mass;
        }

        pthread_barrier_wait(data->barrier);

        for (int index = data->start_index; index < data->end_index; index++) {
            struct particle *current_particle = &data->all_particles[index];
            current_particle->velocity_x += (data->delta_time * current_particle->accel_x);
            current_particle->velocity_y += (data->delta_time * current_particle->accel_y);
            current_particle->position_x += (data->delta_time * current_particle->velocity_x);
            current_particle->position_y += (data->delta_time * current_particle->velocity_y);
        }

        pthread_barrier_wait(data->barrier);
    }

    pthread_exit(NULL);
}

/* 
 * Get current wall time in seconds.
 * 
 * Post: Returns the current wall time in seconds.
 */
double get_wall_seconds() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + (double)tv.tv_usec / 1000000;
}

// Main function
int main(int argc, char const *argv[]) { 
    // Check there are 6 parameters, but assume the input is correct (program, star_count, filename, timesteps, delta_t, graphics)
    if (argc != 7) {
        printf("Error: File requires 6 input parameters\n");
        return 1;
    }

    double start_overall_timer = get_wall_seconds(); // Start the timer. Start here to measure read and write computer times with it. 

    // Assign user input to variables, and declare other variables used. 
    int star_count = atoi(argv[1]); 
    char *file_name = malloc(strlen(argv[2]) + 1);
    strcpy(file_name, argv[2]);
    FILE *file = fopen(file_name, "rb");
    time_steps = atoi(argv[3]);
    double delta_time = atof(argv[4]);
    // int graphics = atoi(argv[5]); // Commented out when graphics is not used
    double gravity = (100.0/star_count);
    struct particle *all_particles = malloc(star_count * sizeof(struct particle));

    // Open the input file and use it to populate the particles in the array of structures
    for (int index = 0; index < star_count; index++) {
        if (fread(&all_particles[index].position_x, sizeof(double), 1, file) != 1
            || fread(&all_particles[index].position_y, sizeof(double), 1, file) != 1
            || fread(&all_particles[index].mass, sizeof(double), 1, file) != 1
            || fread(&all_particles[index].velocity_x, sizeof(double), 1, file) != 1
            || fread(&all_particles[index].velocity_y, sizeof(double), 1, file) != 1
            || fread(&all_particles[index].brightness, sizeof(double), 1, file) != 1) {
            fprintf(stderr, "Error reading file\n");
            break;
        }
    }

    // Reading the number of threads from command line arguments
    int thread_count = atoi(argv[6]); 
    pthread_t threads[thread_count];
    struct thread_data thread_data_array[thread_count];
    pthread_barrier_t barrier;

    // Initialize the barrier with the number of threads
    pthread_barrier_init(&barrier, NULL, thread_count);

    // Distribute the work among the threads
    int chunk_size = star_count / thread_count;

    for (int t = 0; t < thread_count; t++) {
        thread_data_array[t].start_index = t * chunk_size;
        thread_data_array[t].end_index = (t == thread_count - 1) ? star_count : (t + 1) * chunk_size;
        thread_data_array[t].all_particles = all_particles;
        thread_data_array[t].star_count = star_count;
        thread_data_array[t].delta_time = delta_time;
        thread_data_array[t].gravity = gravity;
        thread_data_array[t].barrier = &barrier;
        int rc = pthread_create(&threads[t], NULL, calculate_particles, (void*)&thread_data_array[t]);
        if (rc) {
            printf("Error: return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
    }

    // Join all threads to wait for them to finish
    for (int t = 0; t < thread_count; t++) {
        pthread_join(threads[t], NULL);
    }


    // Declare and populate the output file with the results
    FILE *output = fopen("result.gal", "wb");
    for (int output_line = 0; output_line < star_count; output_line++) {
        double position_x = all_particles[output_line].position_x; 
        double position_y = all_particles[output_line].position_y; 
        double mass = all_particles[output_line].mass; 
        double velocity_x = all_particles[output_line].velocity_x; 
        double velocity_y = all_particles[output_line].velocity_y; 
        double brightness = all_particles[output_line].brightness;
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
    free(all_particles);
    //start_overall_timer = get_wall_seconds() - start_overall_timer;  //Commented out when not requiring timing results
    //printf("\nOverall time: %lf\n\n", start_overall_timer);
    return 0;
}
