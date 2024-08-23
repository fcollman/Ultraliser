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

namespace Ultraliser
{

/**
 * @brief The Pixel2 struct
 */
struct Pixel2
{
public:

    /**
     * @brief Pixel2
     * @param i
     * @param j
     */
    Pixel2(int64_t i, int64_t j) : x(i), y(j) { }

public:

    /**
     * @brief x
     */
    int64_t x = -1;

    /**
     * @brief y
     */
    int64_t y = -1;

    /**
     * @brief operator +
     * @param rhs
     * @return
     */
    Pixel2 operator +(Pixel2 rhs) { return Pixel2(x + rhs.x, y + rhs.y); }
};

/**
 * @brief Pixels2
 */
typedef std::vector< Pixel2 > Pixels2;

}
