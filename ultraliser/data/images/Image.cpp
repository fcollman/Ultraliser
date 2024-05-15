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

#include "Image.h"
#include <utilities/Utilities.h>
#include <algorithms/floodfill/FloodFiller.h>
#include <stack>

namespace Ultraliser
{

Image::Image(const size_t &width, const size_t &height)
    : _width(width)
    , _height(height)
    , _numberPixels(width * height)
{
    // Allocate the memory of the image
    _allocateMemory();
}

Image::Image(const std::string& imagePath, bool readMask)
{
    // Read the image
    if (readMask)
    {
        _readMask(imagePath);
    }
    else
    {
        _readPPM(imagePath);
    }
}

void Image::_allocateMemory()
{
    // Allocate the array that should contain the data, and initialize all the elements to Zero
    _data = new uint8_t[_numberPixels]();
}

void Image::_freeMemory()
{
    delete [] _data;
}

size_t Image::dimension(const size_t& i) const
{
    if (i == 0)
    {
        return _width;
    }
    else if (i == 1)
    {
        return _height;
    }
    else
    {
        LOG_WARNING("Image::dimension accepts ONLY 0 or 1!");
        return 0;
    }
}

// Function to skip over comments in the PPM file
void skipComments(std::ifstream& file)
{
    char c = file.peek();
    while (c == '#' || c == '\n' || c == '\r')
    {
        // Skip until the end of the line
        file.ignore(std::numeric_limits< std::streamsize >::max(), '\n');
        c = file.peek();
    }
}

void Image::_readPPM(const std::string &imagePath)
{
    struct Pixel {
        unsigned char red, green, blue;
    };

    std::ifstream file(imagePath, std::ios::binary);
    if (!file)
    {
        LOG_ERROR("Error: Unable to open file [%s]", imagePath.c_str());
        return;
    }

    std::string magic;
    file >> magic;
    LOG_INFO("Magic %s", magic.c_str());
    if (magic != "P6")
    {
        LOG_ERROR("Error: Not a valid PPM file [%s]", imagePath.c_str());
        return;
    }

    // Skip comments after the magic number
    skipComments(file);

    file >> _width >> _height;
    LOG_INFO("Dimensions %d %d", _width, _height);
    int maxColorValue;
    file >> maxColorValue;

    if (maxColorValue != 255)
    {
        LOG_ERROR("Error: Unsupported max color value. "
                  "This program only supports 8-bit PPM files.");
        return;
    }

    // Allocate the data
    _data = new uint8_t[_width * _height];

    uint8_t* pixelsData = new uint8_t[_width * _height * 3];
    file.read(reinterpret_cast< char* >(pixelsData), _width * _height * 3);

    int count = 0;
    size_t pIndex = 0;
    for (size_t i = 0; i < _width * _height; ++i)
    {
        uint8_t r = unsigned(pixelsData[pIndex]); pIndex++;
        uint8_t g = unsigned(pixelsData[pIndex]); pIndex++;
        uint8_t b = unsigned(pixelsData[pIndex]); pIndex++;

        // float pixelValue = (0.2126f * p.red) + (0.7152f * p.green) + (0.0722f * p.blue);
        const uint8_t value = std::max(std::max(r, g), b);
        if (value > 0) { count++; }

        // const int64_t index = _width * _height - ii;
        _data[i] = (value > 0) ? 255 : 0;
    }

    delete [] pixelsData;

    LOG_INFO("Count : %d", count);
}

void Image::_readMask(const std::string& imagePath)
{
    std::ifstream file(imagePath);
    if (!file)
    {
        LOG_ERROR("Error: Unable to open file [%s]", imagePath.c_str());
        return;
    }

    file >> _width >> _height;

    // Allocate the data
    _data = new uint8_t[_width * _height];

    int count = 0;
    for (size_t i = 0; i < _width * _height; ++i)
    {
        int value;
        file >> value;

        if (value > 0) { count++; }

        // const int64_t index = _width * _height - ii;
        _data[i] = value;
    }

    LOG_INFO("Count : %d", count);
}

void Image::writeMask(const std::string &prefix) const
{
    std::stringstream stream;
    stream << prefix << ".mask";

    std::ofstream file(stream.str(), std::ios::app);

    // Write the header
    file << _width << " " << _height << "\n";

    int count = 0;
    for (size_t i = 0; i < _width * _height; ++i)
    {
        int index = i;
        size_t value = _data[index];
        if (value > 0) { count++; }

        if (PIXEL_COLOR(_data[index]) == WHITE)
            value = 255;
        else if (PIXEL_COLOR(_data[index]) == GRAY)
            value = 128;
        else
            value = 0;
        file << value << "\n";
    }

    file.close();
}

void Image::writePPM(const std::string &prefix) const
{
    std::stringstream stream;
    stream << prefix << ".ppm";

    std::ofstream file(stream.str(), std::ios::binary);

    // Write the header
    file << "P6\n" << _width << " " << _height << "\n255\n";

    int count = 0;
    for (size_t i = 0; i < _width * _height; ++i)
    {
        int index = i;
        uint8_t value = _data[index];
        if (value > 0) { count++; }

//        if (PIXEL_COLOR(_data[index]) == WHITE)
//            value = 255;
//        else if (PIXEL_COLOR(_data[index]) == GRAY)
//            value = 128;
//        else
//            value = 0;

        file.put(static_cast<unsigned char>(value)); // Red
        file.put(static_cast<unsigned char>(value)); // Green
        file.put(static_cast<unsigned char>(value)); // Blue
    }

    file.close();
    // LOG_INFO("Count : %d", count);
}

std::vector< std::vector< Image::ImagePixel > > Image::_getComponents()
{
    size_t* labels = new size_t[_width * _height]();

    int labelCount = 1;

    // Define 8-connectivity offsets
    std::vector<ImagePixel> neighbors = {
        ImagePixel(-1, -1), ImagePixel(-1, 0), ImagePixel(-1, 1),
        ImagePixel( 0, -1),                    ImagePixel( 0, 1),
        ImagePixel( 1, -1), ImagePixel( 1, 0), ImagePixel( 1, 1) };

    std::vector< std::vector< ImagePixel > > components;

    // Iterate through each pixel of the image
    for (int x = 0; x < _width; ++x)
    {
        for (int y = 0; y < _height; ++y)
        {
            const size_t index = x + _width * y;

            // Skip background pixels and pixels already labeled
            if (_data[index] == 0 || labels[index] != 0) { continue; }

            // Start a new connected component
            std::stack< ImagePixel > _stack;
            std::vector< ImagePixel > _component;

            _stack.push(ImagePixel(x, y));

            // Explore the connected pixels using a stack
            while (!_stack.empty())
            {
                ImagePixel current = _stack.top();
                _stack.pop();
                _component.push_back(current);

                const size_t currentIndex = current.x + _width * current.y;
                labels[currentIndex] = labelCount;

                // Check 8-connectivity neighbors
                for (const ImagePixel& offset : neighbors)
                {
                    ImagePixel neighbor = current + offset;
                    size_t nIndex = neighbor.x + _width * neighbor.y;

                    const size_t neighbourIndex = neighbor.x + _width * neighbor.y;
                    if (neighbor.x >= 0 && neighbor.x < _width &&
                        neighbor.y >= 0 && neighbor.y < _height &&
                        _data[nIndex] != 0 && labels[neighbourIndex] == 0)
                    {
                        _stack.push(neighbor);
                    }
                }
            }

            components.push_back(_component);
            _component.clear();
            labelCount++;
        }
    }

    delete [] labels;

    neighbors.clear();
    neighbors.shrink_to_fit();

    // Return the components
    return components;
}

void Image::floodFill()
{
    for (size_t i = 0; i < _width * _height; ++i)
    {
        if (_data[i] == 255) { _data[i] = 128; }
    }

    for (size_t i = 0; i < _width * _height; ++i)
    {
        if (_data[i] == 0) { _data[i] = 255; }
    }

    // Flood Filler
    PIXEL_COLOR newColor = WHITE;
    PIXEL_COLOR oldColor = BLACK;
    FloodFiller::fill(this, _width, _height, 0, 0, newColor, oldColor);
}

void Image::_fillComponent(const std::vector< ImagePixel >& component,
                           size_t componentIndex, size_t imageIndex)
{
    // If the component is composed of a single pixel, return, there is nothing to be filled
    if (component.size() == 1)
    {
        setPixelColor(component[0].x, component[0].y, WHITE);
        return;
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

    if (xMax >= _width - 1)
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

    if (yMax >= _height - 1)
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

    FloodFiller::fill(&componentImage, componentImage.getWidth(),componentImage.getHeight(), xStart, yStart, newColor, oldColor);

    for (int64_t i = 0; i < width; ++i)
    {
        for (int64_t j = 0; j < height; ++j)
        {
            // Component color
            const auto componentColor = componentImage.getPixelColor(i + gap, j + gap);
            const int64_t xPixel = xMin +  i ;
            const int64_t yPixel = yMin +  j ;

            bool filled = componentImage.isFilled(i + gap, j + gap);

            if (xPixel > 0 && xPixel < _width && yPixel > 0 && yPixel < _height)
            {
                if (filled)
                {
                    setPixelColor(xPixel, yPixel, PIXEL_COLOR::WHITE);
                }
            }
        }
    }
}

Pixels2 Image::_getFilledPixelsAfterFloodFilling(const std::vector< ImagePixel >& component,
                                                 size_t componentIndex,
                                                 size_t imageIndex)
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

    if (xMax >= _width - 1)
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

    if (yMax >= _height - 1)
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

            if (xPixel > 0 && xPixel < _width && yPixel > 0 && yPixel < _height)
            {
                if (filled) { filledPixels.push_back(Pixel2(xPixel, yPixel)); }
            }
        }
    }

    return filledPixels;
}

Pixels2 Image::getFilledPixelsAfterFloodFilling(size_t imageIndex)
{
    Pixels2 filledPixels;

    auto components = _getComponents();
    for (size_t i = 0; i < components.size(); ++i)
    {
        auto& component = components[i];
        if (component.size() == 0) { continue; }
        else
        {
            auto filledPixelsInComponent =
                    _getFilledPixelsAfterFloodFilling(component, i, imageIndex);
            filledPixels.insert(filledPixels.end(),
                                filledPixelsInComponent.begin(), filledPixelsInComponent.end());
        }
    }

    return filledPixels;
}

void Image::fillComponents(size_t imageIndex)
{
    auto components = _getComponents();
    for (size_t i = 0; i < components.size(); ++i)
    {
        auto& component = components[i];
        if (component.size() > 0)
        {
            _fillComponent(component, i, imageIndex);
        }
    }
}

Image::~Image()
{
    _freeMemory();
}

}
