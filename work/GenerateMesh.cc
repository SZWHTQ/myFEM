#include <gmsh.h>

#include <iostream>
#include <vector>

#include "GenerateMesh.h"

int generate_mesh(std::vector<double>& nodeCoord,
                  std::vector<size_t>& elementNodeTags,
                  std::vector<size_t>& boundaryNodeTags, double L, double B,
                  double a, double b, double lc, double refinementFactor,
                  bool isSerendipity, int meshAlgorithm) {
    // Initialize the Gmsh library
    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 0);
    // gmsh::option::setNumber("Geometry.OCCKernel", 1);
    try {
        // Start a new model
        gmsh::model::add("RectangleWithHole");

        // Outer rectangle points
        std::vector<int> pointsTag;
        pointsTag.push_back(gmsh::model::occ::addPoint(0, 0, 0, lc));
        pointsTag.push_back(gmsh::model::occ::addPoint(L / 2, 0, 0, lc));
        pointsTag.push_back(gmsh::model::occ::addPoint(L / 2, B / 2, 0, lc));
        pointsTag.push_back(gmsh::model::occ::addPoint(0, B / 2, 0, lc));
        pointsTag.push_back(
            gmsh::model::occ::addPoint(a, 0, 0, lc / refinementFactor));
        pointsTag.push_back(
            gmsh::model::occ::addPoint(0, b, 0, lc / refinementFactor));

        // Outer rectangle lines
        std::vector<int> linesTag;
        linesTag.push_back(
            gmsh::model::occ::addLine(pointsTag[4], pointsTag[1]));
        linesTag.push_back(
            gmsh::model::occ::addLine(pointsTag[1], pointsTag[2]));
        linesTag.push_back(
            gmsh::model::occ::addLine(pointsTag[2], pointsTag[3]));
        linesTag.push_back(
            gmsh::model::occ::addLine(pointsTag[3], pointsTag[5]));
        linesTag.push_back(gmsh::model::occ::addEllipseArc(
            pointsTag[5], pointsTag[0], pointsTag[4], pointsTag[4]));

        // Outer loop (line loop)
        int outerCurveLoop = gmsh::model::occ::addCurveLoop(linesTag);

        // Define the first and second half of the ellipse
        int ellipseArcTag = gmsh::model::occ::addEllipseArc(
            pointsTag[4], pointsTag[0], pointsTag[4], pointsTag[5]);
        int ellipseBoundary1Tag =
            gmsh::model::occ::addLine(pointsTag[5], pointsTag[0]);
        int ellipseBoundary2Tag =
            gmsh::model::occ::addLine(pointsTag[0], pointsTag[4]);

        // Create a curve loop and plane surface
        int innerCurveLoopTag = gmsh::model::occ::addCurveLoop(
            {ellipseArcTag, ellipseBoundary1Tag, ellipseBoundary2Tag});

        // Plane surface with inclusion
        gmsh::model::occ::addPlaneSurface({outerCurveLoop}, 1);

        gmsh::model::occ::addPlaneSurface({innerCurveLoopTag}, 2);

        // Synchronize the model
        gmsh::model::occ::synchronize();

        gmsh::model::addPhysicalGroup(2, {1}, MATRIX_SURFACE_TAG,
                                      "Matrix Surface");
        gmsh::model::addPhysicalGroup(2, {2}, INCLUSION_SURFACE_TAG,
                                      "Inclusion Surface");

        linesTag.push_back(ellipseBoundary1Tag);
        linesTag.push_back(ellipseBoundary2Tag);
        gmsh::model::addPhysicalGroup(1, linesTag, BOUNDARY_TAG, "Boundary");
        gmsh::model::addPhysicalGroup(1, {ellipseArcTag}, INTERFACE_TAG,
                                      "Interface");
        gmsh::model::addPhysicalGroup(0, pointsTag, VERTEX_TAG,
                                      "Matrix Vertex");

        // gmsh::model::occ::removeAllDuplicates();

        // gmsh::option::setNumber("Geometry.AutoCoherence", 1);

        // Set mesh options for second-order elements
        gmsh::option::setNumber("Mesh.ElementOrder", 2);
        // serendipity
        gmsh::option::setNumber("Mesh.SecondOrderIncomplete", isSerendipity);
        gmsh::option::setNumber("Mesh.RecombinationAlgorithm", 2);
        gmsh::option::setNumber("Mesh.RecombineAll", 1);
        // Set mesh algorithm
        gmsh::option::setNumber("Mesh.Algorithm", meshAlgorithm);

        // Generate mesh
        gmsh::model::mesh::generate(2);
        gmsh::model::mesh::removeDuplicateNodes();
        // gmsh::model::mesh::removeDuplicateElements();

        // Save mesh to file
#if WRITE_INP
        gmsh::write("PlateWithInclusion.inp");
#endif
#if WRITE_MSH
        gmsh::write("PlateWithInclusion.msh");
#endif

        // Retrieve the node coordinates
        std::vector<size_t> nodeTags;
        std::vector<double> parametricCoord;
        gmsh::model::mesh::getNodes(nodeTags, nodeCoord, parametricCoord);

        // Retrieve the element nodes
        std::vector<size_t> elementTags;
        int serendipityTag =
            gmsh::model::mesh::getElementType("Quadrangle", 2, isSerendipity);
        gmsh::model::mesh::getElementsByType(serendipityTag, elementTags,
                                             elementNodeTags);
        // int triangleTag =
        //     gmsh::model::mesh::getElementType("Triangle", 2, isSerendipity);

        // Retrieve the boundary nodes
        // // Vertex
        // std::vector<int> vertexTags;
        // gmsh::model::getEntitiesForPhysicalGroup(0, VERTEX_TAG, vertexTags);
        // boundaryNodeTags.insert(boundaryNodeTags.end(), vertexTags.begin(),
        //                         vertexTags.end());
        // Boundary
        std::vector<int> lineTags;
        gmsh::model::getEntitiesForPhysicalGroup(1, INTERFACE_TAG, lineTags);
        std::cout << "lineTags.size() = " << lineTags.size() << std::endl;
        for (const auto& line : lineTags) {
            std::vector<size_t> nodeTags;
            std::vector<double> nodeCoords;
            std::vector<double> parametricCoords;

            // Get nodes on the line
            gmsh::model::mesh::getNodes(nodeTags, nodeCoords, parametricCoords,
                                        1, line);

            // Append boundary node coordinates to the output vector
            boundaryNodeTags.insert(boundaryNodeTags.end(), nodeTags.begin(),
                                    nodeTags.end());
        }
        // gmsh::fltk::run();
    } catch (const std::exception& e) {
        std::cerr << e.what();
        return -2;
    }
    // // Finalize the Gmsh library
    // gmsh::finalize();

    return 0;
}