/***************************************************************************************************
 * Copyright (c) 2016 - 2024
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

#ifdef ULTRALISER_USE_OPENMP
#include <omp.h>
#endif

#ifdef ULTRALISER_USE_OPENMP
#define OMP_PARALLEL_FOR _Pragma("omp parallel for")
#else
#define OMP_PARALLEL_FOR
#endif

static omp_lock_t createOMPLock() {omp_lock_t LOCK; omp_init_lock(&LOCK); return LOCK; }

#ifdef ULTRALISER_USE_OPENMP
#define CREATE_OMP_LOCK omp_lock_t ULTRALISER_OMP_LOCK = createOMPLock();
#else
#define CREATE_OMP_LOCK
#endif

#ifdef ULTRALISER_USE_OPENMP
#define OMP_SET_LOCK omp_set_lock(&ULTRALISER_OMP_LOCK );
#else
#define OMP_SET_LOCK
#endif

#ifdef ULTRALISER_USE_OPENMP
#define OMP_UNSET_LOCK omp_unset_lock(&ULTRALISER_OMP_LOCK );
#else
#define OMP_UNSET_LOCK
#endif

