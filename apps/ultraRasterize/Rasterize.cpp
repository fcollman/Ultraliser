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
    // std::string inputPath = "/home/abdellah/testing-components/slice_9_before.ppm";
    std::string maskPathI = "/home/abdellah/testing-components/slice_9_before.mask";
    std::string outputPathI = "/home/abdellah/testing-components/slice_9_mask";


    for (int i = 1; i < 180; ++i)
    {
        std::stringstream mP, oP;
        mP << "/home/abdellah/testing-components/slice_" << i << "_before.mask";
        Image* image = new Image(mP.str(), true);

        TIMER_SET;
        image->fillComponents(i);
        // LOG_STATUS_IMPORTANT("Algorithm  fillComponents Stats.");
        // LOG_STATS(GET_TIME_SECONDS);

        oP << "/home/abdellah/testing-components/slice_" << i << "_mask";
        image->writePPM(oP.str() + "-fillComponents");
        delete image;

        Image* image2 = new Image(mP.str(), true);

        TIMER_RESET;
        image2->floodFill();
        // LOG_STATUS_IMPORTANT("Algorithm  floodFill Stats.");
        // LOG_STATS(GET_TIME_SECONDS);

        image2->writePPM(oP.str() + "-floodFill");
        delete image2;
    }
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

