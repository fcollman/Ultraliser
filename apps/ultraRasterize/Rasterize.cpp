/***************************************************************************************************
 * Copyright (c) 2016 - 2022
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

#include <Ultraliser.h>
#include <AppCommon.h>
#include <AppArguments.h>
#include <data/images/Image.h>

namespace Ultraliser
{

void run(int argc , const char** argv)
{
    // Read the image
    std::string maskPath = "/home/abdellah/testing-components/shapes.mask";
    std::string outputPath = "/home/abdellah/testing-components/shapes_output";

    Image* image = new Image(maskPath, true);

    TIMER_SET;
    image->fillComponents();
    LOG_STATUS_IMPORTANT("Algorithm  fillComponents Stats.");
    LOG_STATS(GET_TIME_SECONDS);

    image->writePPM(outputPath + "-fillComponents");
    delete image;
    image = nullptr;

    image = new Image(maskPath, true);

    TIMER_RESET;
    image->floodFill();
    LOG_STATUS_IMPORTANT("Algorithm  floodFill Stats.");
    LOG_STATS(GET_TIME_SECONDS);

    image->writePPM(outputPath + "-floodFill");
    delete image;
    image = nullptr;
}
}

int main(int argc , const char** argv)
{
    TIMER_SET;

    Ultraliser::run(argc, argv);

    LOG_STATUS_IMPORTANT("Ultralization Stats.");
    LOG_STATS(GET_TIME_SECONDS);

    ULTRALISER_DONE;
}

