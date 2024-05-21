/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Marwan Abdellah < marwan.abdellah@epfl.ch >
 *      Juan Jose Garcia Cantero < juanjose.garcia@epfl.ch>
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

#include "FloodFiller.h"
#include <stack>

namespace Ultraliser
{

void FloodFiller::fill(Image* image,
                       const int64_t &nx, const int64_t &ny,
                       const int64_t &x, const int64_t &y,
                       PIXEL_COLOR background, PIXEL_COLOR fillColor)
{
    if (background != image->getPixelColorWOBC(x, y))
        return;

    /// Bottom
    int64_t y1 = y;
    while (y1 < ny && background == image->getPixelColorWOBC(x, y1))
        image->setPixelColor(x, y1++, fillColor);

    /// Top
    y1 = y - 1;
    while (y1 >= 0 && background == image->getPixelColorWOBC(x, y1))
        image->setPixelColor(x, y1--, fillColor);

    /// Left
    y1 = y;
    while (y1 < ny && fillColor == image->getPixelColorWOBC(x, y1))
    {
        if (x > 0 && background == image->getPixelColorWOBC(x - 1, y1))
            fill(image, nx, ny, x - 1, y1, background, fillColor);
        y1++;
    }

    y1 = y - 1;
    while (y1 >= 0 && fillColor == image->getPixelColorWOBC(x, y1))
    {
        if (x > 0 && background == image->getPixelColorWOBC(x - 1, y1))
            fill(image, nx, ny, x - 1, y1, background, fillColor);
        y1--;
    }

    /// Right
    y1 = y;
    while (y1 < ny && fillColor == image->getPixelColorWOBC(x, y1))
    {
        if (x < nx - 1 && background == image->getPixelColorWOBC(x + 1, y1))
            fill(image, nx, ny, x + 1, y1, background, fillColor);
        y1++;
    }

    y1 = y - 1;
    while (y1 >= 0 && fillColor == image->getPixelColorWOBC(x, y1))
    {
        if (x < nx - 1 && background == image->getPixelColorWOBC(x + 1, y1))
            fill(image, nx, ny, x + 1, y1, background, fillColor);
        y1--;
    }
}

void FloodFiller::fill(Image* image, Pixels2& filledPixels,
                       const int64_t &nx, const int64_t &ny,
                       const int64_t &x, const int64_t &y,
                       PIXEL_COLOR background, PIXEL_COLOR fillColor)
{
    if (background != image->getPixelColorWOBC(x, y))
        return;

    /// Bottom
    int64_t y1 = y;
    while (y1 < ny && background == image->getPixelColorWOBC(x, y1))
    {
        image->setPixelColor(x, y1, fillColor);
        filledPixels.push_back(Pixel2(x, y1));
        y1++;
    }

    /// Top
    y1 = y - 1;
    while (y1 >= 0 && background == image->getPixelColorWOBC(x, y1))
    {
        image->setPixelColor(x, y1, fillColor);
        filledPixels.push_back(Pixel2(x, y1));
        y1--;
    }

    /// Left
    y1 = y;
    while (y1 < ny && fillColor == image->getPixelColorWOBC(x, y1))
    {
        if (x > 0 && background == image->getPixelColorWOBC(x - 1, y1))
            fill(image, filledPixels, nx, ny, x - 1, y1, background, fillColor);
        y1++;
    }

    y1 = y - 1;
    while (y1 >= 0 && fillColor == image->getPixelColorWOBC(x, y1))
    {
        if (x > 0 && background == image->getPixelColorWOBC(x - 1, y1))
            fill(image, filledPixels, nx, ny, x - 1, y1, background, fillColor);
        y1--;
    }

    /// Right
    y1 = y;
    while (y1 < ny && fillColor == image->getPixelColorWOBC(x, y1))
    {
        if (x < nx - 1 && background == image->getPixelColorWOBC(x + 1, y1))
            fill(image, filledPixels, nx, ny, x + 1, y1, background, fillColor);
        y1++;
    }

    y1 = y - 1;
    while (y1 >= 0 && fillColor == image->getPixelColorWOBC(x, y1))
    {
        if (x < nx - 1 && background == image->getPixelColorWOBC(x + 1, y1))
            fill(image, filledPixels, nx, ny, x + 1, y1, background, fillColor);
        y1--;
    }
}

Pixels2 FloodFiller::fillComponent(const std::vector< Pixel2 >& component,
                                   const size_t& sliceWidth, const size_t& sliceHeight)
{
    Pixels2 filledPixels;

    // Single pixel
    if (component.size() == 1)
    {
        filledPixels.push_back(Pixel2(component[0].x, component[0].y));
        return filledPixels;
    }

    int64_t xMin = std::numeric_limits<int64_t>::max();
    int64_t xMax = -1;
    int64_t yMin = std::numeric_limits<int64_t>::max();
    int64_t yMax = -1;

    // Get the bounds of the component
    for (const auto& p : component)
    {
        if (p.x > xMax) xMax = p.x; if (p.x < xMin) xMin = p.x;
        if (p.y > yMax) yMax = p.y; if (p.y < yMin) yMin = p.y;
    }

    const size_t width = xMax - xMin + 1;
    const size_t height = yMax - yMin + 1;

    // Create an image
    const size_t gap = 8;
    const size_t componentWidth = width + 2 * gap;
    const size_t componentHeight = height + 2 * gap;
    Image componentImage(componentWidth, componentHeight);
    componentImage.fill(WHITE);

    int xStart = 0;
    int yStart = 0;

    if (xMin == 0)
    {
        for (size_t j = 0; j < componentHeight; ++j)
        {
            componentImage.setPixelColor(gap - 1, j, GRAY);
        }

        xStart = gap;
    }

    if (xMax >= sliceWidth - 1)
    {
        for (size_t j = 0; j < componentHeight; ++j)
        {
            componentImage.setPixelColor(componentWidth - gap - 1, j, GRAY);
        }
    }

    if (yMin == 0)
    {
        for (size_t i = 0; i < componentWidth; ++i)
        {
            componentImage.setPixelColor(i, gap - 1, GRAY);
        }

        yStart = gap;
    }

    if (yMax >= sliceHeight - 1)
    {
        for (size_t i = 0; i < componentWidth; ++i)
        {
            componentImage.setPixelColor(i, componentHeight - gap - 1, GRAY);
        }
    }

    // Fill the image with the component object
    for (const auto& p : component)
    {
        int xComponent = p.x - xMin;
        int yComponent = p.y - yMin;
        componentImage.setPixelColor(gap + xComponent , gap + yComponent , PIXEL_COLOR::GRAY);
    }

    // Flood Filler
    PIXEL_COLOR newColor = WHITE;
    PIXEL_COLOR oldColor = BLACK;

    // Apply the flood filling algorithm
    FloodFiller::fill(&componentImage, componentImage.getWidth(),componentImage.getHeight(),
                      xStart, yStart, newColor, oldColor);

    // Get the filled pixels
    for (int64_t i = 0; i < width; ++i)
    {
        for (int64_t j = 0; j < height; ++j)
        {
            // Component color
            const auto componentColor = componentImage.getPixelColor(i + gap, j + gap);
            const int64_t xPixel = xMin +  i ;
            const int64_t yPixel = yMin +  j ;

            bool filled = componentImage.isFilled(i + gap, j + gap);

            if (xPixel > 0 && xPixel < sliceWidth && yPixel > 0 && yPixel < sliceHeight)
            {
                if (filled) { filledPixels.push_back(Pixel2(xPixel, yPixel)); }
            }
        }
    }

    return filledPixels;
}

}
