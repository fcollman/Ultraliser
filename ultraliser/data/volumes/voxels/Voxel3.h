/***************************************************************************************************
 * Copyright (c) 2016 - 2023
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

namespace Ultraliser
{

/**
 * @brief The Voxel3 struct
 */
struct Voxel3
{
    /**
     * @brief Voxel3
     * @param i
     * @param j
     * @param k
     */
    Voxel3(size_t i, size_t j, size_t k) : x(i), y(j), z(k) { }

    /**
     * @brief x
     */
    size_t x;

    /**
     * @brief y
     */
    size_t y;

    /**
     * @brief z
     */
    size_t z;
};

/**
 * @brief Voxels3
 */
typedef std::vector< Voxel3 > Voxels3;

/**
 * @brief Voxel3Bucket
 */
typedef std::vector< Voxels3 > Voxel3Bucket;

}
