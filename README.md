# Particle-Simulation
This project contains an N-body problem solution which calculates particles in a 2D cartesian plane over a variable ammount of time. 
The purpose of this project is first to create efficient serial code in C, and then parallise it with Pthreads and OpenMP. 

# Usage
To run, compile any of the programs by navigating to their file and running the makefile:
```
Make
```
After which the program becomes an executeable that can be ran as so: (using GCC)

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
In order to update the position for a single particle, its force (relative to every other particle!) needs to be calculated, and afterwards the velocity can be calculated, which is then used to update the position. 

Starting with the Netwons law of gravitation in two dimensions, the force on particle i relative to j (denoted as a function F(i, j)) can be calculated as:
_Force(i) = -gravity * mass(i) * summation(mass(j) / (distance(i, j) + plummer_sphere)^3 * vector(i, j))_

Which is then used to calculate accelleration on each particle:
_Force(i) / Mass(i)_

And after calculating the accelleration on all particles, the following equation is used to calculate the velocity:
_Velocity(i) + Delta time * Acceleration(i)_

Which can finally be used to update the position of i using:
_Position(i) + Delta Time * Velocity(i)_
