/**********************************************************
* Sparse 2015 Project
* Directed by: Euripides Montagne.
*              eurip@cs.ucf.edu
*
* Coded by: Edward Aymerich. 
*           edward.aymerich@knights.ucf.edu
*
* Copyright 2015 - University of Central Florida.
**********************************************************/

#ifndef SP_CLOCK
#define SP_CLOCK

#include <sys/time.h>

/**
 * This structure is used to measure and retrieve
 * execution time.
 */
typedef struct {
	struct timeval begin; // Stores initial time.
	struct timeval end;   // Stores ending time.
	double sec;           // Stores measured time in seconds.
} Clock;

void clock_start(Clock *clock);
void clock_stop(Clock *clock);

/**
 * Begins measuring time.
 */
void clock_start(Clock *clock){
	gettimeofday(&(clock->begin),NULL);
}

/**
 * Stops measuring time.
 * Time measure is store in clock->sec,
 * in seconds.
 */
void clock_stop(Clock *clock){
	gettimeofday(&(clock->end),NULL);
	
	clock->sec = (clock->end).tv_sec;
	clock->sec -= (clock->begin).tv_sec;
	clock->sec *= 1000000.0;
	clock->sec += (clock->end).tv_usec;
	clock->sec -= (clock->begin).tv_usec;
	//clock->sec /= 1000000.0;
}

#endif // SP_CLOCK
