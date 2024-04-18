#include <gmsh.h>

#include <iostream>
#include <vector>

#include "GenerateMesh.h"

int generate_mesh(std::vector<double>& nodeCoord,
                  std::vector<size_t>& elementNodeTags,
                  std::vector<size_t>& elementMaterialTags,
                  std::vector<size_t>& interfaceNodeTags,
                  toml::table settings) {
    // Define geometry
    double L = settings["Rectangle"]["L"].value_or(2);
    double B = settings["Rectangle"]["B"].value_or(0.9);
    double a = B / settings["Ellipse"]["ksi"].value_or(3.0);
    double b = a / settings["Ellipse"]["a_b"].value_or(1. / 3);
    bool convertToSquare =
        settings["Ellipse"]["convertToSquare"].value_or(false);
    double lc = settings["Mesh"]["size"].value_or(0.02);
    double refinementFactor = settings["Mesh"]["refinementFactor"].value_or(8);
    int meshAlgorithm = settings["Mesh"]["Algorithm"].value_or(8);
    bool isSerendipity = settings["Mesh"]["Serendipity"].value_or(true);
    bool isPlaneStress = settings["Mesh"]["planeStress"].value_or(true);
    bool writeInp = settings["Mesh"]["writeInp"].value_or(false);
    bool writeMsh = settings["Mesh"]["writeMsh"].value_or(false);
    bool runFltk = settings["Mesh"]["runFltk"].value_or(false);

    int err = generate_mesh(nodeCoord, elementNodeTags, elementMaterialTags,
                            interfaceNodeTags, L, B, a, b, lc, refinementFactor,
                            isSerendipity, meshAlgorithm, convertToSquare,
                            writeInp, writeMsh, runFltk);
    return err;
}

int generate_mesh(std::vector<double>& nodeCoord,
                  std::vector<size_t>& elementNodeTags,
                  std::vector<size_t>& elementMaterialTags,
                  std::vector<size_t>& interfaceNodeTags, double L, double B,
                  double a, double b, double lc, double refinementFactor,
                  bool isSerendipity, int meshAlgorithm, bool convertToSquare,
                  bool writeInp, bool writeMsh, bool runFltk) {
    gmsh::option::setNumber("General.Terminal", 0);
    try {
        // Start a new model
        gmsh::model::add("RectangleWithHole");

        // Outer rectangle points
        std::vector<int> pointsTag;
        pointsTag.push_back(gmsh::model::occ::addPoint(0, 0, 0, lc));
        pointsTag.push_back(gmsh::model::occ::addPoint(L / 2, 0, 0, lc));
        pointsTag.push_back(gmsh::model::occ::addPoint(L / 2, B / 2, 0, lc));
        pointsTag.push_back(gmsh::model::occ::addPoint(0, B / 2, 0, lc));

        // Outer loop (line loop)
        std::vector<int> linesTag;
        std::vector<int> interfaceCurveTags = {pointsTag[4]};
        int boundary1Tag = 0;
        int boundary2Tag = 0;
        int outerCurveLoop = 0;
        int innerCurveLoopTag = 0;
        if (convertToSquare) {
            double squareLength = std::sqrt(M_PI * a * b) * 0.5;
            pointsTag.push_back(gmsh::model::occ::addPoint(
                squareLength, 0, 0, lc / refinementFactor));
            pointsTag.push_back(gmsh::model::occ::addPoint(
                0, squareLength, 0, lc / refinementFactor));
            pointsTag.push_back(gmsh::model::occ::addPoint(
                squareLength, squareLength, 0, lc / refinementFactor));

            // Outer loop (line loop)
            linesTag.push_back(
                gmsh::model::occ::addLine(pointsTag[4], pointsTag[1]));
            linesTag.push_back(
                gmsh::model::occ::addLine(pointsTag[1], pointsTag[2]));
            linesTag.push_back(
                gmsh::model::occ::addLine(pointsTag[2], pointsTag[3]));
            linesTag.push_back(
                gmsh::model::occ::addLine(pointsTag[3], pointsTag[5]));
            linesTag.push_back(
                gmsh::model::occ::addLine(pointsTag[5], pointsTag[6]));
            linesTag.push_back(
                gmsh::model::occ::addLine(pointsTag[6], pointsTag[4]));
            outerCurveLoop = gmsh::model::occ::addCurveLoop(linesTag);

            // Inner loop (line loop)
            int interfaceLine1Tag = gmsh::model::occ::addLine(
                pointsTag[4], pointsTag[6]);  // interfaceCurve
            int interfaceLine2Tag = gmsh::model::occ::addLine(
                pointsTag[6], pointsTag[5]);  // interfaceCurve
            boundary1Tag =
                gmsh::model::occ::addLine(pointsTag[5], pointsTag[0]);
            boundary2Tag =
                gmsh::model::occ::addLine(pointsTag[0], pointsTag[4]);

            // Create a curve loop
            innerCurveLoopTag = gmsh::model::occ::addCurveLoop(
                {interfaceLine1Tag, interfaceLine2Tag, boundary1Tag,
                 boundary2Tag});

            // Interface curve
            interfaceCurveTags = {linesTag[4], linesTag[5], interfaceLine1Tag,
                                  interfaceLine2Tag};
        } else {
            pointsTag.push_back(
                gmsh::model::occ::addPoint(a, 0, 0, lc / refinementFactor));
            pointsTag.push_back(
                gmsh::model::occ::addPoint(0, b, 0, lc / refinementFactor));

            // Outer loop (line loop)
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
            outerCurveLoop = gmsh::model::occ::addCurveLoop(linesTag);

            // Inner loop (line loop)
            int ellipseArcTag = gmsh::model::occ::addEllipseArc(
                pointsTag[4], pointsTag[0], pointsTag[4], pointsTag[5]);
            boundary1Tag =
                gmsh::model::occ::addLine(pointsTag[5], pointsTag[0]);
            boundary2Tag =
                gmsh::model::occ::addLine(pointsTag[0], pointsTag[4]);

            // Create a curve loop
            innerCurveLoopTag = gmsh::model::occ::addCurveLoop(
                {ellipseArcTag, boundary1Tag, boundary2Tag});

            // Interface curve
            interfaceCurveTags = {linesTag[4], ellipseArcTag};
        }

        // Plane surface with inclusion
        gmsh::model::occ::addPlaneSurface({outerCurveLoop}, 1);

        gmsh::model::occ::addPlaneSurface({innerCurveLoopTag}, 2);

        // Synchronize the model
        gmsh::model::occ::synchronize();

        gmsh::model::addPhysicalGroup(2, {1}, MATRIX_SURFACE_TAG,
                                      "Matrix Surface");
        gmsh::model::addPhysicalGroup(2, {2}, INCLUSION_SURFACE_TAG,
                                      "Inclusion Surface");

        gmsh::model::addPhysicalGroup(1, interfaceCurveTags, INTERFACE_TAG,
                                      "Interface");

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
        if (writeInp) {
            gmsh::write("PlateWithInclusion.inp");
        }
        if (writeMsh) {
            gmsh::write("PlateWithInclusion.msh");
        }

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

        // Retrieve the element materials by physical group
        {
            elementMaterialTags.resize(elementTags.size());
            std::vector<size_t> eTags;
            std::vector<size_t> nTags;
            // matrix
            std::vector<int> surfaces;
            gmsh::vectorpair dimTags;
            gmsh::model::getEntities(dimTags, 2);
            gmsh::model::getEntitiesForPhysicalGroup(2, MATRIX_SURFACE_TAG,
                                                     surfaces);
            bool isFirst = true;
            size_t startTag = 0;
            for (const auto& s : surfaces) {
                gmsh::model::mesh::getElementsByType(serendipityTag, eTags,
                                                     nTags, s);
                if (isFirst) {
                    isFirst = false;
                    startTag = eTags[0];
                }
                for (size_t i = 0; i < eTags.size(); ++i) {
                    elementMaterialTags[eTags[i] - startTag] = 0;
                }
            }

            // inclusion
            surfaces.clear();
            gmsh::model::getEntitiesForPhysicalGroup(2, INCLUSION_SURFACE_TAG,
                                                     surfaces);
            for (const auto& s : surfaces) {
                gmsh::model::mesh::getElementsByType(serendipityTag, eTags,
                                                     nTags, s);
                for (size_t i = 0; i < eTags.size(); ++i) {
                    elementMaterialTags[eTags[i] - startTag] = 1;
                }
            }
        }

        // Retrieve the interface nodes
        {
            // interface
            std::vector<int> lines;
            gmsh::model::getEntitiesForPhysicalGroup(1, INTERFACE_TAG, lines);
            for (const auto& line : lines) {
                std::vector<size_t> nodeTags_;
                std::vector<double> nodeCoords_;
                std::vector<double> parametricCoords_;

                // Get nodes on the line
                gmsh::model::mesh::getNodes(nodeTags_, nodeCoords_,
                                            parametricCoords_, 1, line);

                // Append boundary node coordinates to the output vector
                interfaceNodeTags.insert(interfaceNodeTags.end(),
                                         nodeTags_.begin(), nodeTags_.end());
            }
            interfaceNodeTags.push_back(pointsTag[4]);
            interfaceNodeTags.push_back(pointsTag[5]);
            if (convertToSquare) {
                interfaceNodeTags.push_back(pointsTag[6]);
            }
        }
        if (runFltk) {
            gmsh::fltk::run();
        }
    } catch (const std::exception& e) {
        std::cerr << e.what();
        return -2;
    }
    // // Finalize the Gmsh library
    // gmsh::finalize();

    return 0;
}