# Particle-Simulation
This project contains an N-body problem solution which calculates particles in a 2D cartesian plane over a variable ammount of time. 
The purpose of this project is first to create efficient serial code in C, and then parallise it with Pthreads and OpenMP. 

# Usage
To run, compile any of the programs by navigating to their file and running the makefile:
Once inside the directory of the desired program (either using one of the two parallel versions, or the serial) and compile it by running:
```
Make
```
This will compile the code using the -Ofast flag (most aggressivle optimisations) and compiler warnings. This will compile using GCC. Afterwards, the executeable file can be ran using:

```
./galsim number_of_stars filename.gal number_of_timesteps delta_time graphics 
```
And if using a threadded version:
```
./galsim number_of_stars filename.gal number_of_timesteps delta_time graphics number_of_threads
```
So if this program was to be executed in serial with the current binary file (and 200 steps to compare with the provided file):
```
./galsim 100 ../Data/ellipse_N_0100.gal 200 10e3 0 
```

# Checking results 
Place the results.gal binary file inside the "Data" folder, compile compare_gal_files.c and then run:
```
./compare 100 results.gal ellipse_N_00100_after200steps.gal
```

# Running with timers
At the moment there is a timing feature embedded inside the parallel versions of the code. To test how fast it runs, the original C program needs to be edited so that the timing function is not commented out. This will print to terminal after executing the program, with a runtime to show the speedups when altering the thread count. 
The timing function is in the final 3 lines of the program, so change these from comments into code:
```
//start_overall_timer = get_wall_seconds() - start_overall_timer;  //Commented out when not requiring timing results
//printf("\nOverall time: %lf\n\n", start_overall_timer);
```

# Equations used:
In the galaxy simulation, the dynamic behavior of a single particle relies on an intricate process of calculating the gravitational interactions with every other particle. This calculation enables us to update the particle's position using the following sequence of steps:

1. **Calculate the Force on a Particle**: Start with Newton's law of gravitation in two dimensions to compute the force on particle `i` relative to `j`, represented as a function `F(i, j)`:
_Force(i) = -gravity * mass(i) * summation(mass(j) / (distance(i, j) + plummer_sphere)^3 * vector(i, j))_


2. **Compute Acceleration**: Divide the force by the mass of each particle to find the acceleration:
_Force(i) / Mass(i)_


3. **Update Velocity**: With the acceleration calculated, you can then compute the new velocity for each particle using:
_Velocity(i) + Delta time * Acceleration(i)_


4. **Update Position**: Finally, the newly calculated velocity is used to update the position of particle `i` by:
_Position(i) + Delta Time * Velocity(i)_



This sequence of equations brings a particle to life in the simulation, allowing it to move and interact in a manner that adheres to the principles of physics. By following these equations, the simulation provides a portrayal of the forces at work in a galaxy, showing how every particle's motion is influenced by the cumulative effect of the gravitational force exerted by all other particles.


# Issues:
This version of the code has not been implemented with the graphic display features, so the "graphics" setting should always be set to "0". 
