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

#include <data/images/Image.h>
#include <data/volumes/voxels/Voxel3.h>

namespace Ultraliser
{

/**
 * @brief The FloodFiller class
 */
class FloodFiller
{
public:

     /**
     * @brief fill
     * Flood-fill the image.
     *
     * @param image
     * An image to be flood-filled.
     * @param nx
     * Number of pixels on the x-axis.
     * @param ny
     * Number of pixels on the y-axis.
     * @param x
     * Seed pixel on the x-axis.
     * @param y
     * Seed pixel on the y-axis.
     * @param background
     * The default background color.
     * @param fillColor
     * The filling color.
     */
    static void fill(Image* image, const int64_t &nx, const int64_t &ny,
                     const int64_t &x, const int64_t &y,
                     PIXEL_COLOR background, PIXEL_COLOR fillColor);

    /**
     * @brief fill
     * @param image
     * @param filledPixels
     * @param nx
     * @param ny
     * @param x
     * @param y
     * @param background
     * @param fillColor
     */
    static void fill(Image* image, Pixels2& filledPixels,
                     const int64_t &nx, const int64_t &ny,
                     const int64_t &x, const int64_t &y,
                     PIXEL_COLOR background, PIXEL_COLOR fillColor);

    /**
     * @brief fillComponent
     * @param component
     * @param sliceWidth
     * @param sliceHeight
     * @return
     */
    static Pixels2 fillComponent(const std::vector< Pixel2 >& component,
                                 const size_t& sliceWidth, const size_t& sliceHeight);
};

}
