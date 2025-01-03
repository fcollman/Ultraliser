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

namespace Ultraliser
{
AppOptions* parseArguments(const int& argc , const char** argv)
{
    std::unique_ptr< AppArguments > args = std::make_unique< AppArguments >(
        argc, argv, COPYRIGHT
        "Optimize multi region mesh.");

    // Get all the options
    AppOptions* options = args->getOptions();

    LOG_TITLE("Creating Context");

    // Initialize context, once everything is in place and all the options are verified
    options->initializeContext();

    // Return the executable options
    return options;
}

void run(int argc , const char** argv)
{
    // Parse the arguments and get the tool options
    // auto options = parseArguments(argc, argv);

    // Get the list of meshes
    const std::string meshPathD1 = "/Users/abdellah/z43/meshes/test_64-dmc.obj";
    const std::string meshPathD2 = "/Users/abdellah/z43/meshes/test_255-dmc.obj";
    
    // Construct the domains 
    auto domain1 = new Mesh(meshPathD1);
    auto domain2 = new Mesh(meshPathD2);

    // Build the KdTree of both 
    auto d1PointCloud = domain1->constructPointCloud();
    auto d1KdTree = KdTree::from(d1PointCloud);
    
    auto d2PointCloud = domain2->constructPointCloud();
    auto d2KdTree = KdTree::from(d2PointCloud);




    size_t numberSharedPoints = 0;
    for (size_t i = 0; i < d1PointCloud.size(); ++i)
    {
        auto& p1 = d1PointCloud[i];
        auto nPoint = d2KdTree.findNearestPoint(p1);
        if (nPoint.distance < 0.000000000001)
        {
            // Simply increment 
            numberSharedPoints++;
            
            // Add the vertex here to the vertex markers 
            domain1->setVertexMarkers(i);
        } 
    }
    
    std::cout << "Number Points: " << numberSharedPoints << std::endl;

    numberSharedPoints = 0;
    for (size_t i = 0; i < d2PointCloud.size(); ++i)
    {
        auto& p2 = d2PointCloud[i];
        auto nPoint = d1KdTree.findNearestPoint(p2);
        if (nPoint.distance < 0.000000000001)
        {
            // Simply increment 
            numberSharedPoints++;
            
            // Add the vertex here to the vertex markers 
            domain2->setVertexMarkers(i);
        } 
    }

    // Impement the smoothing 
    for (size_t i = 0; i < 5; ++i)
    {
        domain1->smooth();
        domain2->smooth();
        
        // inputMesh->smoothNormals(false);
    }

    // Export the output mesh 
    domain1->exportMesh("/Users/abdellah/z43/meshes/test_64-dmc-snap1", true, true, true, true);
    
    domain2->exportMesh("/Users/abdellah/z43/meshes/test_255-dmc-snap2", true, true, true, true);

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
