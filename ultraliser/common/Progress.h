/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Marwan Abdellah < marwan.abdellah@epfl.ch >
 *
 * This file is part of Ultraliser < https://github.com/BlueBrain/Ultraliser >
 *
 * This library is free software; you can redistribute it and/or modify it under the terms of the
 * GNU General Public License version 3.0 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the GNU General Public License along with this library;
 * if not, write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
 * MA 02111-1307, USA.
 * You can also find it on the GNU web site < https://www.gnu.org/licenses/gpl-3.0.en.html >
 **************************************************************************************************/

#pragma once

#include <common/Headers.hh>
#include <common/Logging.h>

/**
 * @brief printProgressBar
 * Print the progress bar in a loop.
 *
 * @param current
 * Current count
 * @param total
 * Total count
 * @param barLength
 * The total langth of the progress bar.
 */
void printProgressBar(const size_t& current,
                      const size_t& total,
                      const size_t barLength = 50);

/**
 * @brief printFractionProgressBar
 * Prints a very simple progress bar in a loop that only increments evey 10% of the loop.
 *
 * @param current
 * Current count
 * @param total
 * Total count
 * @param barLength
 * The total langth of the progress bar.
 */
void printFractionProgressBar(const size_t& current,
                              const size_t& total,
                              const size_t barLength = 50);

/**
 * @brief progressUpdate
 * Update the progress of the current thread in a loop.
 *
 * @param progressValue
 * The progress value to be updated.
 */
void progressUpdate(size_t& progressValue);

// Setting a counter
#define LOOP_COUNTER_SET size_t COUNTER = 0
#define LOOP_COUNTER_RESET COUNTER = 0

// Prints a simple message before starting the loop
#define LOOP_STARTS(MESSAGE) (printf("\t%s \n", MESSAGE));

// Print the progress in a loop
#ifdef ENABLE_PROGRESS_BAR
#define LOOP_PROGRESS(PROGRESS, TOTAL) (printProgressBar(PROGRESS, TOTAL))
#else
#define LOOP_PROGRESS(PROGRESS, TOTAL) { }
#endif

// Print the progress in a loop only in fractions
#ifdef ENABLE_PROGRESS_BAR
#define LOOP_PROGRESS_FRACTION(PROGRESS, TOTAL) (printFractionProgressBar(PROGRESS, TOTAL))
#else
#define LOOP_PROGRESS_FRACTION(PROGRESS, TOTAL) {}
#endif

// Print the status after the loop is done
#ifdef ENABLE_PROGRESS_BAR
#define LOOP_DONE { LOOP_PROGRESS(100, 100); printf(" \n"); }
#else
#define LOOP_DONE { }
#endif

// The progress variable itself
#ifdef ENABLE_PROGRESS_BAR
#define PROGRESS ULTRALISER_PROGRESS
#else
#define PROGRESS
#endif

// Set the progress to zero
#ifdef ENABLE_PROGRESS_BAR
#define PROGRESS_SET size_t ULTRALISER_PROGRESS = 0
#else
#define PROGRESS_SET
#endif

// Set the progress at a specific starting value
#ifdef ENABLE_PROGRESS_BAR
#define PROGRESS_SET_AT_VALUE( VALUE ) size_t ULTRALISER_PROGRESS = VALUE;
#else
#define PROGRESS_SET_AT_VALUE( VALUE ) { }
#endif

// Reset the progress
#ifdef ENABLE_PROGRESS_BAR
#define PROGRESS_RESET (ULTRALISER_PROGRESS = 0)
#else
#define PROGRESS_RESET { }
#endif

// Update the progress bar
#ifdef ENABLE_PROGRESS_BAR
#define PROGRESS_UPDATE (progressUpdate(ULTRALISER_PROGRESS))
#else
#define PROGRESS_UPDATE { }
#endif

