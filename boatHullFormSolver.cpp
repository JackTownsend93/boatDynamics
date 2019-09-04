/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2017 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at 
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include "palabos3D.h"
#include "palabos3D.hh"

#include "pypal/headers3D.h"
#include "pypal/headers3D.hh"

#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <string>

using namespace plb;

typedef double T;

#define WAVE_ABSORPTION

#ifdef WAVE_ABSORPTION
#define DESCRIPTOR descriptors::AbsorbingWaveD3Q19Descriptor
#else
#define DESCRIPTOR descriptors::D3Q19Descriptor
#endif

// Descriptor used for post-processing smoothing
#define SM_DESCRIPTOR descriptors::AdvectionDiffusionD3Q7Descriptor

#define OUT     0
#define IN      1

struct SimulationParameters {

    /*
     * Parameters set by the user.
     * All user input variables and all data in external input files must be in the same system of units.
     */

    T yVel;
    std::vector<T> xDomain;                     // Extent in the x-direction of the physical simulation domain.
    std::vector<T> yDomain;                     // Extent in the y-direction of the physical simulation domain.
    std::vector<T> zDomain;                     // Extent in the z-direction of the physical simulation domain.
    T inletAbsorbingZoneWidth;                  // Absorbing zone widths.
    T outletAbsorbingZoneWidth;
    T lateralAbsorbingZoneWidth;
    T topAbsorbingZoneWidth;
    T fluidHeight;                              // Initial height of the fluid.

    std::string boatStl;                        // Name of the file with the boat hull form geometry.

    T rho;                                      // Fluid density in physical units.
    T nu;                                       // Fluid kinematic viscosity in physical units.
    T surfaceTension;                           // Surface tension coefficient in physical units.

    T inflationParameter;                       // Parameter for the voxelizer.

    T characteristicLength;                     // Length to define dx.
    plint resolution;                           // Total number of lattice nodes in the characteristic length.

    T inletVelocity;                            // Inlet BC.

    T u_Ref;                                    // Reference velocity.
    T u_LB;                                     // Lattice velocity.

    T A;                                        // Wave generation force parameters.
    T P;
    std::vector<T> xWaveDomain;
    std::vector<T> yWaveDomain;
    std::vector<T> zWaveDomain;

    plint maxIter;                              // Maximum number of iterations.

    T cSmago;                                   // Smagorinsky parameter.

    bool strongRepelling;                       // Parameter for the Immersed Boundary method.

    T ambientPressure;                          // Absolute pressure at infinity.

    std::string outDir;                         // Output directory.
    plint statIter;                             // Number of iterations for terminal output.
    plint outIter;                              // Number of iterations for disk output.
    plint cpIter;                               // Number of iterations for checkpointing.
    plint abIter;                               // Number of iterations for checking for user-driven program abortion.

    bool excludeInteriorForOutput;              // Exclude the interior of the boat for pressure and force output or not?
    int numPresLaplaceIter;                     // Number of Laplacian smoothing iterations for the pressure post-processing.

    bool outputInDomain;                        // Save data on disk in a volume domain or not?
    Cuboid<T> outputCuboid;                     // Volume domain for disk output.

    bool outputOnSlices;                        // Save data on disk on a set of slices or not?
    std::vector<T> xPositions;                  // Positions of the x-slices for output.
    std::vector<T> xyRange;                     // y range of the x-slices.
    std::vector<T> xzRange;                     // z range of the x-slices.
    std::vector<T> yPositions;                  // Positions of the y-slices for output.
    std::vector<T> yzRange;                     // z range of the y-slices.
    std::vector<T> yxRange;                     // x range of the y-slices.
    std::vector<T> zPositions;                  // Positions of the z-slices for output.
    std::vector<T> zxRange;                     // x range of the z-slices.
    std::vector<T> zyRange;                     // y range of the z-slices.

    std::string abortFileName;                  // File for signaling program abortion.
    std::string xmlContinueFileName;            // XML file for restarting.
    std::string baseFileName;                   // Basename of the checkpoint files.
    bool useParallelIO;                         // For a desktop PC this should be "false", for a cluster "true".

    /*
     * Parameters NOT set by the user.
     */

    T yVel_LB;
    plint nx, ny, nz;
    plint fluidHeight_LB;
    Box3D fullDomain;
    Box3D bottom, top;
    Box3D initialFluidDomain;
    T dx;
    T dt;
    T rho_LB;
    Array<plint,6> numAbsorbingCells;
    plint totalNumAbsorbingCells;
    T surfaceTension_LB;
    T inletVelocity_LB;
    T initialVelocity_LB;
    T gravity_LB;
    T A_LB;
    T P_LB;
    Box3D waveDomain;
    T omega;
    bool incompressibleModel;
    Array<T,3> physicalLocation;
    RawConnectedTriangleMesh<T>* connectedMesh;
    std::vector<Array<T,3> > vertices;
    std::vector<T> areas;
    std::vector<int> flags;
    Box3D outputDomain;
    std::vector<Box3D> xSlices;
    std::vector<Box3D> ySlices;
    std::vector<Box3D> zSlices;
    std::vector<Box3D> allOutputDomains;
    std::vector<std::string> allOutputDomainNames;
    bool saveDynamicContent;
    plint fileNamePadding;
} param;

T toPhys(T lbVal, plint direction, T dx, Array<T,3> const& location)
{
    PLB_ASSERT(direction >= 0 && direction <= 2);
    return (lbVal * dx + location[direction]);
}

Array<T,3> toPhys(Array<T,3> const& lbVal, T dx, Array<T,3> const& location)
{
    return (lbVal * dx + location);
}

T toLB(T physVal, plint direction, T dx, Array<T,3> const& location)
{
    PLB_ASSERT(direction >= 0 && direction <= 2);
    return (physVal - location[direction]) / dx;
}

Array<T,3> toLB(Array<T,3> const& physVal, T dx, Array<T,3> const& location)
{
    return (physVal - location) / dx;
}

void setParameters(std::string xmlInputFileName)
{
    XMLreader document(xmlInputFileName);

    document["geometry"]["simulationDomain"]["x"].read(param.xDomain);
    plbIOError(param.xDomain.size() != 2 || util::lessEqual(param.xDomain[1], param.xDomain[0]),
            "The x-extent of the simulation domain is wrong.");
    document["geometry"]["simulationDomain"]["y"].read(param.yDomain);
    plbIOError(param.yDomain.size() != 2 || util::lessEqual(param.yDomain[1], param.yDomain[0]),
            "The y-extent of the simulation domain is wrong.");
    document["geometry"]["simulationDomain"]["z"].read(param.zDomain);
    plbIOError(param.zDomain.size() != 2 || util::lessEqual(param.zDomain[1], param.zDomain[0]),
            "The z-extent of the simulation domain is wrong.");

    document["geometry"]["inletAbsorbingZoneWidth"].read(param.inletAbsorbingZoneWidth);
    document["geometry"]["outletAbsorbingZoneWidth"].read(param.outletAbsorbingZoneWidth);
    document["geometry"]["lateralAbsorbingZoneWidth"].read(param.lateralAbsorbingZoneWidth);
    document["geometry"]["topAbsorbingZoneWidth"].read(param.topAbsorbingZoneWidth);

    document["geometry"]["fluidHeight"].read(param.fluidHeight);

    document["geometry"]["boatStl"].read(param.boatStl);

    document["fluid"]["rho"].read(param.rho);
    document["fluid"]["nu"].read(param.nu);
    document["fluid"]["surfaceTension"].read(param.surfaceTension);

    document["solver"]["yVel"].read(param.yVel);
    document["solver"]["inflationParameter"].read(param.inflationParameter);
    plbIOError(param.inflationParameter < (T) 0 || param.inflationParameter > (T) 1,
            "The inflationParameter must take values between 0 and 1.");

    document["solver"]["characteristicLength"].read(param.characteristicLength);
    document["solver"]["resolution"].read(param.resolution);

    document["solver"]["inletVelocity"].read(param.inletVelocity);

    document["solver"]["uRef"].read(param.u_Ref);
    document["solver"]["uLB"].read(param.u_LB);

    document["solver"]["A"].read(param.A);
    document["solver"]["P"].read(param.P);
    document["solver"]["waveDomain"]["x"].read(param.xWaveDomain);
    plbIOError(param.xWaveDomain.size() != 2 || util::lessEqual(param.xWaveDomain[1], param.xWaveDomain[0]),
            "The x-extent of the wave domain is wrong.");
    document["solver"]["waveDomain"]["y"].read(param.yWaveDomain);
    plbIOError(param.yWaveDomain.size() != 2 || util::lessEqual(param.yWaveDomain[1], param.yWaveDomain[0]),
            "The y-extent of the wave domain is wrong.");
    document["solver"]["waveDomain"]["z"].read(param.zWaveDomain);
    plbIOError(param.zWaveDomain.size() != 2 || util::lessEqual(param.zWaveDomain[1], param.zWaveDomain[0]),
            "The z-extent of the wave domain is wrong.");

    document["solver"]["maxIter"].read(param.maxIter);

    document["solver"]["cSmago"].read(param.cSmago);

    document["solver"]["strongRepelling"].read(param.strongRepelling);

    document["solver"]["ambientPressure"].read(param.ambientPressure);

    std::string outDir;
    document["output"]["outDir"].read(outDir);
    if (outDir[outDir.size() - 1] != '/') {
        outDir += '/';
    }
    param.outDir = outDir;
    abortIfCannotCreateFileInDir(param.outDir, "plb-checkfile.txt");

    document["output"]["statIter"].read(param.statIter);
    document["output"]["outIter"].read(param.outIter);
    document["output"]["cpIter"].read(param.cpIter);
    document["output"]["abIter"].read(param.abIter);

    document["output"]["excludeInteriorForOutput"].read(param.excludeInteriorForOutput);
    document["output"]["numPresLaplaceIter"].read(param.numPresLaplaceIter);

    document["output"]["outputInDomain"].read(param.outputInDomain);
    if (param.outputInDomain) {
        std::vector<plint> x, y, z;
        document["output"]["outputDomain"]["x"].read(x);
        plbIOError(x.size() != 2 || util::lessEqual(x[1], x[0]), "The x-extent of the outputDomain is wrong");
        document["output"]["outputDomain"]["y"].read(y);
        plbIOError(y.size() != 2 || util::lessEqual(y[1], y[0]), "The y-extent of the outputDomain is wrong");
        document["output"]["outputDomain"]["z"].read(z);
        plbIOError(z.size() != 2 || util::lessEqual(z[1], z[0]), "The z-extent of the outputDomain is wrong");
        param.outputCuboid.lowerLeftCorner[0] = x[0];
        param.outputCuboid.lowerLeftCorner[1] = y[0];
        param.outputCuboid.lowerLeftCorner[2] = z[0];
        param.outputCuboid.upperRightCorner[0] = x[1];
        param.outputCuboid.upperRightCorner[1] = y[1];
        param.outputCuboid.upperRightCorner[2] = z[1];
    }

    document["output"]["outputOnSlices"].read(param.outputOnSlices);
    if (param.outputOnSlices) {
        document["output"]["outputSlices"]["xSlices"]["xPositions"].read(param.xPositions);
        document["output"]["outputSlices"]["xSlices"]["yRange"].read(param.xyRange);
        plbIOError(param.xyRange.size() != 2 || util::lessEqual(param.xyRange[1], param.xyRange[0]),
                "The y-range of the x-slices is wrong");
        document["output"]["outputSlices"]["xSlices"]["zRange"].read(param.xzRange);
        plbIOError(param.xzRange.size() != 2 || util::lessEqual(param.xzRange[1], param.xzRange[0]),
                "The z-range of the x-slices is wrong");

        document["output"]["outputSlices"]["ySlices"]["yPositions"].read(param.yPositions);
        document["output"]["outputSlices"]["ySlices"]["zRange"].read(param.yzRange);
        plbIOError(param.yzRange.size() != 2 || util::lessEqual(param.yzRange[1], param.yzRange[0]),
                "The z-range of the y-slices is wrong");
        document["output"]["outputSlices"]["ySlices"]["xRange"].read(param.yxRange);
        plbIOError(param.yxRange.size() != 2 || util::lessEqual(param.yxRange[1], param.yxRange[0]),
                "The x-range of the y-slices is wrong");

        document["output"]["outputSlices"]["zSlices"]["zPositions"].read(param.zPositions);
        document["output"]["outputSlices"]["zSlices"]["xRange"].read(param.zxRange);
        plbIOError(param.zxRange.size() != 2 || util::lessEqual(param.zxRange[1], param.zxRange[0]),
                "The x-range of the z-slices is wrong");
        document["output"]["outputSlices"]["zSlices"]["yRange"].read(param.zyRange);
        plbIOError(param.zyRange.size() != 2 || util::lessEqual(param.zyRange[1], param.zyRange[0]),
                "The x-range of the z-slices is wrong");
    }

    document["output"]["abortFileName"].read(param.abortFileName);
    document["output"]["xmlContinueFileName"].read(param.xmlContinueFileName);
    document["output"]["baseFileName"].read(param.baseFileName);
    document["output"]["useParallelIO"].read(param.useParallelIO);
}

void computeOutputDomain(Cuboid<T> const& cuboid, Box3D& box)
{
    if (!param.outputInDomain) {
        return;
    }

    Array<T,3> llc = cuboid.lowerLeftCorner;
    Array<T,3> urc = cuboid.upperRightCorner;

    plint x0 = util::roundToInt(toLB(llc[0], 0, param.dx, param.physicalLocation));
    plint y0 = util::roundToInt(toLB(llc[1], 1, param.dx, param.physicalLocation));
    plint z0 = util::roundToInt(toLB(llc[2], 2, param.dx, param.physicalLocation));

    plint x1 = util::roundToInt(toLB(urc[0], 0, param.dx, param.physicalLocation));
    plint y1 = util::roundToInt(toLB(urc[1], 1, param.dx, param.physicalLocation));
    plint z1 = util::roundToInt(toLB(urc[2], 2, param.dx, param.physicalLocation));

    PLB_ASSERT(x1 >= x0 && y1 >= y0 && z1 >= z0);

    box = Box3D(x0, x1, y0, y1, z0, z1);

#ifdef PLB_DEBUG
    bool intersectionOK =
#endif
        intersect(box, param.fullDomain, box);
    PLB_ASSERT(intersectionOK);
}

void computeOutputSlices()
{
    if (!param.outputOnSlices) {
        return;
    }

    {
        param.xSlices.clear();

        plint y0 = util::roundToInt(toLB(param.xyRange[0], 1, param.dx, param.physicalLocation));
        plint y1 = util::roundToInt(toLB(param.xyRange[1], 1, param.dx, param.physicalLocation));
        plint z0 = util::roundToInt(toLB(param.xzRange[0], 2, param.dx, param.physicalLocation));
        plint z1 = util::roundToInt(toLB(param.xzRange[1], 2, param.dx, param.physicalLocation));
        PLB_ASSERT(y1 >= y0 && z1 >= z0);

        for (size_t i = 0; i < param.xPositions.size(); i++) {
            plint xPos = util::roundToInt(toLB(param.xPositions[i], 0, param.dx, param.physicalLocation));
            plint x0 = xPos - 1;
            plint x1 = xPos + 1;
            param.xSlices.push_back(Box3D(x0, x1, y0, y1, z0, z1));

#ifdef PLB_DEBUG
            bool intersectionOK =
#endif
                intersect(param.xSlices.back(), param.fullDomain, param.xSlices.back());
            PLB_ASSERT(intersectionOK);
        }
    }

    {
        param.ySlices.clear();

        plint z0 = util::roundToInt(toLB(param.yzRange[0], 2, param.dx, param.physicalLocation));
        plint z1 = util::roundToInt(toLB(param.yzRange[1], 2, param.dx, param.physicalLocation));
        plint x0 = util::roundToInt(toLB(param.yxRange[0], 0, param.dx, param.physicalLocation));
        plint x1 = util::roundToInt(toLB(param.yxRange[1], 0, param.dx, param.physicalLocation));

        PLB_ASSERT(z1 >= z0 && x1 >= x0);

        for (size_t i = 0; i < param.yPositions.size(); i++) {
            plint yPos = util::roundToInt(toLB(param.yPositions[i], 1, param.dx, param.physicalLocation));
            plint y0 = yPos - 1;
            plint y1 = yPos + 1;
            param.ySlices.push_back(Box3D(x0, x1, y0, y1, z0, z1));

#ifdef PLB_DEBUG
            bool intersectionOK =
#endif
                intersect(param.ySlices.back(), param.fullDomain, param.ySlices.back());
            PLB_ASSERT(intersectionOK);
        }
    }

    {
        param.zSlices.clear();

        plint x0 = util::roundToInt(toLB(param.zxRange[0], 0, param.dx, param.physicalLocation));
        plint x1 = util::roundToInt(toLB(param.zxRange[1], 0, param.dx, param.physicalLocation));
        plint y0 = util::roundToInt(toLB(param.zyRange[0], 1, param.dx, param.physicalLocation));
        plint y1 = util::roundToInt(toLB(param.zyRange[1], 1, param.dx, param.physicalLocation));

        PLB_ASSERT(x1 >= x0 && y1 >= y0);

        for (size_t i = 0; i < param.zPositions.size(); i++) {
            plint zPos = util::roundToInt(toLB(param.zPositions[i], 2, param.dx, param.physicalLocation));
            plint z0 = zPos - 1;
            plint z1 = zPos + 1;
            param.zSlices.push_back(Box3D(x0, x1, y0, y1, z0, z1));

#ifdef PLB_DEBUG
            bool intersectionOK =
#endif
                intersect(param.zSlices.back(), param.fullDomain, param.zSlices.back());
            PLB_ASSERT(intersectionOK);
        }
    }
}

void orderAllOutputDomains()
{
    param.allOutputDomains.clear();
    param.allOutputDomainNames.clear();

    if (param.outputInDomain) {
        param.allOutputDomains.push_back(param.outputDomain);
        param.allOutputDomainNames.push_back("domain");
    }

    if (param.outputOnSlices) {
        size_t numXdigits = util::val2str(param.xSlices.size()).length();
        for (size_t i = 0; i < param.xSlices.size(); i++) {
            param.allOutputDomains.push_back(param.xSlices[i]);
            param.allOutputDomainNames.push_back(createFileName("slice_x_", i, numXdigits+1));
        }

        size_t numYdigits = util::val2str(param.ySlices.size()).length();
        for (size_t i = 0; i < param.ySlices.size(); i++) {
            param.allOutputDomains.push_back(param.ySlices[i]);
            param.allOutputDomainNames.push_back(createFileName("slice_y_", i, numYdigits+1));
        }

        size_t numZdigits = util::val2str(param.zSlices.size()).length();
        for (size_t i = 0; i < param.zSlices.size(); i++) {
            param.allOutputDomains.push_back(param.zSlices[i]);
            param.allOutputDomainNames.push_back(createFileName("slice_z_", i, numZdigits+1));
        }
    }
}

void setDerivedParameters()
{
    // Derived quantities.
    Cuboid<T> fullCuboid(Array<T,3>(param.xDomain[0], param.yDomain[0], param.zDomain[0]),
                         Array<T,3>(param.xDomain[1], param.yDomain[1], param.zDomain[1]));

    param.dx = param.characteristicLength / (T) (param.resolution - 1);

    T lx = fullCuboid.x1() - fullCuboid.x0();
    T ly = fullCuboid.y1() - fullCuboid.y0();
    T lz = fullCuboid.z1() - fullCuboid.z0();
    param.physicalLocation = Array<T,3>(fullCuboid.x0(), fullCuboid.y0(), fullCuboid.z0());
    param.nx = util::roundToInt(lx / param.dx);
    param.ny = util::roundToInt(ly / param.dx) + 1;
    param.nz = util::roundToInt(lz / param.dx);
    param.fullDomain = Box3D(0, param.nx - 1, 0, param.ny - 1, 0, param.nz - 1);

    param.bottom = Box3D(0, param.nx - 1, 0,            0,            0, param.nz - 1);
    param.top    = Box3D(0, param.nx - 1, param.ny - 1, param.ny - 1, 0, param.nz - 1);

    param.initialFluidDomain = param.fullDomain;
    param.initialFluidDomain.y0++;
    param.fluidHeight_LB = util::roundToInt(param.fluidHeight / param.dx);
    param.initialFluidDomain.y1 = param.fluidHeight_LB - 1;

    param.saveDynamicContent = true;
    param.fileNamePadding = 8;

    param.dt = param.u_LB / param.u_Ref * param.dx;

    param.rho_LB = (T) 1;

    param.yVel_LB = param.yVel * (param.dt / param.dx);

    param.inletVelocity_LB = param.inletVelocity * (param.dt / param.dx);

    param.initialVelocity_LB = param.inletVelocity_LB;

    param.surfaceTension_LB = (param.rho_LB/param.rho) * param.dt*param.dt / (param.dx*param.dx*param.dx) * param.surfaceTension;
    param.gravity_LB = (T) 9.81 * (param.dt*param.dt / param.dx);

    param.numAbsorbingCells[0] = util::roundToInt(param.inletAbsorbingZoneWidth / param.dx);
    param.numAbsorbingCells[1] = util::roundToInt(param.outletAbsorbingZoneWidth / param.dx);
    param.numAbsorbingCells[2] = 0;
    param.numAbsorbingCells[3] = util::roundToInt(param.topAbsorbingZoneWidth / param.dx);
    param.numAbsorbingCells[4] = util::roundToInt(param.lateralAbsorbingZoneWidth / param.dx);
    param.numAbsorbingCells[5] = util::roundToInt(param.lateralAbsorbingZoneWidth / param.dx);

    param.totalNumAbsorbingCells = 0;
    for (int iZone = 0; iZone < 6; iZone++) {
        param.totalNumAbsorbingCells += param.numAbsorbingCells[iZone];
    }

    param.A_LB = param.A;
    param.P_LB = param.P / param.dt;
    param.waveDomain.x0 = util::roundToInt(toLB(param.xWaveDomain[0], 0, param.dx, param.physicalLocation));
    param.waveDomain.x1 = util::roundToInt(toLB(param.xWaveDomain[1], 0, param.dx, param.physicalLocation));
    param.waveDomain.y0 = util::roundToInt(toLB(param.yWaveDomain[0], 1, param.dx, param.physicalLocation));
    param.waveDomain.y1 = util::roundToInt(toLB(param.yWaveDomain[1], 1, param.dx, param.physicalLocation));
    param.waveDomain.z0 = util::roundToInt(toLB(param.zWaveDomain[0], 2, param.dx, param.physicalLocation));
    param.waveDomain.z1 = util::roundToInt(toLB(param.zWaveDomain[1], 2, param.dx, param.physicalLocation));
#ifdef PLB_DEBUG
    bool intersectionOK =
#endif
        intersect(param.waveDomain, param.fullDomain, param.waveDomain);
    PLB_ASSERT(intersectionOK);

    T nu_LB = param.nu * param.dt / (param.dx * param.dx);
    param.omega = (T) 1 / (DESCRIPTOR<T>::invCs2 * nu_LB + (T) 0.5);

    computeOutputDomain(param.outputCuboid, param.outputDomain);
    computeOutputSlices();
    orderAllOutputDomains();
}

T waveForce(T t)
{
    static T pi = std::acos((T) -1);
    return param.A_LB * std::fabs(param.gravity_LB) * std::sin((T) 2 * pi * t / param.P_LB);
}

void setupBoatGeometry(bool operateOnOuterBorder, MultiScalarField3D<int>& tags)
{
    TriangleSet<T>* triangleSet = new TriangleSet<T>(param.boatStl);
    triangleSet->translate(-param.physicalLocation);
    triangleSet->scale((T) 1 / param.dx);

    RawConnectedTriangleMesh<T>* mesh = 0;
    {
        RawTriangleMesh<T> rawMesh = triangleSetToRawTriangleMesh(*triangleSet);
        mesh = new RawConnectedTriangleMesh<T>(MeshConnector<T>(rawMesh).generateConnectedMesh());
    }
    inflate(*mesh, param.inflationParameter);
    plint borderWidth = 1;
    MultiNTensorField3D<int>* voxelMatrix = meshToVoxel(*mesh, tags.getBoundingBox(), borderWidth);
    delete mesh; mesh = 0;

    setToConstant(tags, voxelMatrix->scalarView(), voxelFlag::inside, tags.getBoundingBox(), (int) IN);
    setToConstant(tags, voxelMatrix->scalarView(), voxelFlag::innerBorder, tags.getBoundingBox(), (int) IN);
    if (operateOnOuterBorder) {
        setToConstant(tags, voxelMatrix->scalarView(), voxelFlag::outerBorder, tags.getBoundingBox(), (int) IN);
    }
    delete voxelMatrix; voxelMatrix = 0;

    plint maxRefinements = 100;
    T targetLength = (T) 1;
    bool succeeded = triangleSet->refineRecursively(targetLength, maxRefinements);
    if (!succeeded) {
        pcout << std::endl;
        pcout << "WARNING: The target maximum triangle edge length " << targetLength
              << " for the immersed surface was not reached after " << maxRefinements
              << " refinement iterations." << std::endl;
        pcout << std::endl;
        exit(1);
    }

    {
        RawTriangleMesh<T> rawMesh = triangleSetToRawTriangleMesh(*triangleSet);
        delete triangleSet; triangleSet = 0;
        param.connectedMesh = new RawConnectedTriangleMesh<T>(MeshConnector<T>(rawMesh).generateConnectedMesh());
    }

    plint numVertices = param.connectedMesh->getNumVertices();
    param.vertices.resize(numVertices);
    param.areas.resize(numVertices);
    param.flags.resize(numVertices, (int) 1);
    plint uniqueVertexIdTag = param.connectedMesh->getVertexTag("UniqueID");
    RawConnectedTriangleMesh<T>::PVertexIterator vertexIterator(param.connectedMesh->vertexIterator());
    while (!vertexIterator->end()) {
        RawConnectedTriangleMesh<T>::PVertex vertex = vertexIterator->next();
        plint iVertex = vertex->tag(uniqueVertexIdTag);
        param.vertices[iVertex] = vertex->get();
        param.areas[iVertex] = vertex->area();
    }
    pcout << "The number of vertices on the boat surface (after refinement) is: " << numVertices << std::endl;
}

void initialRhoU(plint iX, plint iY, plint iZ, T& rho, Array<T,3>& u)
{
    rho = param.rho_LB;
    if (iY <= param.fluidHeight_LB) {
        rho = param.rho_LB + DESCRIPTOR<T>::invCs2 * std::fabs(param.gravity_LB) * (param.fluidHeight_LB - iY);
    }
    u = Array<T,3>(param.initialVelocity_LB, (T) 0, (T) 0);
}

template <typename T>
class VelFunction {
public:
    T yVel_LB = param.yVel_LB;
    Array<T,3> operator()(pluint id)
    {
        return Array<T,3>(0.0, yVel_LB, 0.0);
    }
};

void printSimulationParameters()
{
    pcout << "Physical simulation domain: [" << param.xDomain[0] << ", " << param.xDomain[1] << "] x ["
                                             << param.yDomain[0] << ", " << param.yDomain[1] << "] x ["
                                             << param.zDomain[0] << ", " << param.zDomain[1] << "]" << std::endl;
    pcout << "inletAbsorbingZoneWidth = " << param.inletAbsorbingZoneWidth << std::endl;
    pcout << "outletAbsorbingZoneWidth = " << param.outletAbsorbingZoneWidth << std::endl;
    pcout << "lateralAbsorbingZoneWidth = " << param.lateralAbsorbingZoneWidth << std::endl;
    pcout << "topAbsorbingZoneWidth = " << param.topAbsorbingZoneWidth << std::endl;
    pcout << "fluidHeight = " << param.fluidHeight << std::endl;

    pcout << "boatStl = " << param.boatStl << std::endl;

    pcout << "rho = " << param.rho << std::endl;
    pcout << "nu = " << param.nu << std::endl;
    pcout << "surfaceTension = " << param.surfaceTension << std::endl;

    pcout << "inflationParameter = " << param.inflationParameter << std::endl;
    pcout << "characteristicLength = " << param.characteristicLength << std::endl;
    pcout << "resolution = " << param.resolution << std::endl;
    pcout << "inletVelocity = " << param.inletVelocity << std::endl;
    pcout << "u_Ref = " << param.u_Ref << std::endl;
    pcout << "u_LB = " << param.u_LB << std::endl;
    pcout << "maxIter = " << param.maxIter << std::endl;
    pcout << "cSmago = " << param.cSmago << std::endl;
    pcout << "strongRepelling = " << (param.strongRepelling ? "true" : "false") << std::endl;
    pcout << "ambientPressure = " << param.ambientPressure << std::endl;

    pcout << "outDir = " << param.outDir << std::endl;
    pcout << "statIter = " << param.statIter << std::endl;
    pcout << "outIter = " << param.outIter << std::endl;
    pcout << "cpIter = " << param.cpIter << std::endl;
    pcout << "abIter = " << param.abIter << std::endl;
    pcout << "excludeInteriorForOutput = " << (param.excludeInteriorForOutput ? "true" : "false") << std::endl;
    pcout << "numPresLaplaceIter = " << param.numPresLaplaceIter << std::endl;
    pcout << "abortFileName = " << param.abortFileName << std::endl;
    pcout << "xmlContinueFileName = " << param.xmlContinueFileName << std::endl;
    pcout << "baseFileName = " << param.baseFileName << std::endl;
    pcout << "useParallelIO = " << (param.useParallelIO ? "true" : "false") << std::endl;

    pcout << "incompressibleModel = " << (param.incompressibleModel ? "true" : "false") << std::endl;
    pcout << "fluidHeight_LB = " << param.fluidHeight_LB << std::endl;
    pcout << "rho_LB = " << param.rho_LB << std::endl;
    pcout << "surfaceTension_LB = " << param.surfaceTension_LB << std::endl;
    pcout << "gravity_LB = " << param.gravity_LB << std::endl;
    pcout << "inletVelocity_LB = " << param.inletVelocity_LB << std::endl;
    pcout << "initialVelocity_LB = " << param.initialVelocity_LB << std::endl;
    for (int iZone = 0; iZone < 6; iZone++) {
        pcout << "numAbsorbingCells[" << iZone << "] = " << param.numAbsorbingCells[iZone] << std::endl;
    }
    T Re = param.u_Ref * param.characteristicLength / param.nu;
    pcout << "Reynolds number: Re = " << Re << std::endl;
    pcout << "omega = " << param.omega << std::endl;
    pcout << "tau = " << 1.0 / param.omega << std::endl;
    pcout << "dx = " << param.dx << std::endl;
    pcout << "dt = " << param.dt << std::endl;
    pcout << "dt / dx = " << param.dt / param.dx << std::endl;
    pcout << "dt / (dx * dx) = " << param.dt / (param.dx * param.dx) << std::endl;
    pcout << "physicalLocation = (" << param.physicalLocation[0] << ", " << param.physicalLocation[1] << ", "
          << param.physicalLocation[2] << ")" << std::endl;
    pcout << std::endl;
}

void applyAbsorbingZones(FreeSurfaceFields3D<T,DESCRIPTOR>& fields, MultiScalarField3D<T>& targetVolumeFraction,
        MultiTensorField3D<T,3>& targetVelocity)
{
    if (param.totalNumAbsorbingCells == 0) {
        return;
    }

    std::vector<MultiBlock3D*> args;
    args.push_back(&fields.volumeFraction);
    args.push_back(&targetVolumeFraction);
    args.push_back(&fields.rhoBar);
    args.push_back(&fields.mass);
    args.push_back(&fields.j);
    args.push_back(&targetVelocity);
    applyProcessingFunctional(new FreeSurfaceSpongeZone3D<T,DESCRIPTOR>(
                param.nx + 1, param.ny, param.nz + 1, param.numAbsorbingCells,
                param.incompressibleModel),
            fields.volumeFraction.getBoundingBox(), args);
}

void writeVTK(FreeSurfaceFields3D<T,DESCRIPTOR>& fields, MultiTensorField3D<T,3>& force, MultiScalarField3D<T>& smoothVF,
        MultiScalarField3D<T>& pressure, Box3D const& outputDomain, ParallelVtkImageOutput3D<T>& vtkOut)
{
    // CAUTION: If more entries are added to the VTK file in this function, the variable "numEntries" in
    //          the "writeResults" function must be adapted accordingly.

    vtkOut.writeData<float>(*extractSubDomain(pressure, outputDomain), "pressure", 1.0);
    vtkOut.writeData<float>(*copyConvert<int,T>(fields.flag, outputDomain), "flag", 1.0);
    MultiScalarField3D<T>* normalVF = extractSubDomain(fields.volumeFraction, outputDomain).release();
    boundScalarField<T>(*normalVF, (T) 0, (T) 1);
    vtkOut.writeData<float>(*normalVF, "volumeFraction", 1.0);
    delete normalVF; normalVF = 0;
    vtkOut.writeData<float>(*extractSubDomain(smoothVF, outputDomain), "volumeFractionFiltered", 1.0);
    std::auto_ptr<MultiTensorField3D<T,3> > v = freeSurfaceComputeForcedVelocity(fields.lattice, force, fields.flag, outputDomain);
    vtkOut.writeData<3,float>(*v, "velocity", param.dx / param.dt);
}

void writeResults(FreeSurfaceFields3D<T,DESCRIPTOR>& fields, MultiTensorField3D<T,3>& force,
        MultiScalarField3D<int>& tags, plint iIter)
{
    MultiScalarField3D<T>* smoothVF = generateMultiScalarField<T>(fields.volumeFraction, 1).release();
    copy(fields.volumeFraction, *smoothVF, fields.volumeFraction.getBoundingBox());
    lbmSmoothenInPlace<T,SM_DESCRIPTOR>(*smoothVF);
    lbmSmoothenInPlace<T,SM_DESCRIPTOR>(*smoothVF);
    boundScalarField<T>(*smoothVF, (T) 0, (T) 1);

    T pressureScale = param.rho * (param.dx * param.dx) / (param.dt * param.dt) * DESCRIPTOR<T>::cs2;
    T pressureOffset = param.ambientPressure -
        param.rho_LB * param.rho * (param.dx * param.dx) / (param.dt * param.dt) * DESCRIPTOR<T>::cs2;
    std::auto_ptr<MultiScalarField3D<T> > pressure = computeDensity(fields.lattice);
    pressure->periodicity().toggle(0, true);
    pressure->periodicity().toggle(1, false);
    pressure->periodicity().toggle(2, true);

    Box3D smoothDomain(param.fullDomain);
    smoothDomain.y0++;
    smoothDomain.y1--;
    for (int iSmooth = 0; iSmooth < param.numPresLaplaceIter; iSmooth++) {
        lbmSmoothenInPlace<T,SM_DESCRIPTOR>(*pressure, smoothDomain);
    }

    multiplyInPlace(*pressure, pressureScale);
    addInPlace(*pressure, pressureOffset);

    // CAUTION: If more entries are added to the VTK file in the "writeVTK" function, the variable
    //          "numEntries" must be adapted accordingly here.

    plint numEntries = 5;

    for (size_t iDomain = 0; iDomain < param.allOutputDomains.size(); iDomain++) {
        Box3D outputDomain = param.allOutputDomains[iDomain];
        std::string domainName = param.allOutputDomainNames[iDomain];
        std::string fname = createFileName(param.outDir + domainName + "_", iIter, param.fileNamePadding);
        ParallelVtkImageOutput3D<T> vtkOut(fname, numEntries, param.dx, param.physicalLocation);
        writeVTK(fields, force, *smoothVF, *pressure, outputDomain, vtkOut);
    }

    // Use a marching-cube algorithm to reconstruct the free surface and write an STL file.

    std::vector<T> isoLevels;
    isoLevels.push_back((T) 0.5);
    typedef TriangleSet<T>::Triangle Triangle;
    std::vector<Triangle> triangles;
    Box3D marchingCubeDomain(param.fullDomain);
    marchingCubeDomain.x0 += 2;
    marchingCubeDomain.x1 -= 2;
    marchingCubeDomain.y0 += 1;
    marchingCubeDomain.y1 -= 1;
    isoSurfaceMarchingCube(triangles, *smoothVF, isoLevels, marchingCubeDomain);
    delete smoothVF; smoothVF = 0;
    TriangleSet<T> set(triangles);
    set.scale(param.dx);
    set.translate(param.physicalLocation);
    set.writeBinarySTL(createFileName(param.outDir + "interface_", iIter, param.fileNamePadding) + ".stl");

    // Export the pressure on the surface of the boat.

    RawConnectedTriangleMesh<T> mesh(*param.connectedMesh);
    plint envelopeWidth = 1;
    std::auto_ptr<MultiScalarField3D<int> > mask = generateMultiScalarField<int>((MultiBlock3D&) fields.flag, envelopeWidth);
    setToConstant(*mask, mask->getBoundingBox(), (int) 0);
    setToConstant(*mask, fields.flag, (int) freeSurfaceFlag::fluid, fields.flag.getBoundingBox(), (int) 1);
    setToConstant(*mask, fields.flag, (int) freeSurfaceFlag::interface, fields.flag.getBoundingBox(), (int) 1);
    if (param.excludeInteriorForOutput) {
        setToConstant(*mask, tags, (int) IN, tags.getBoundingBox(), (int) 0);
    }
    for (plint iLayer = 0; iLayer < 3; iLayer++) {
        applyProcessingFunctional(new MaskedNTensorNeumannInLayersFunctional3D<T>(0, 1, -1),
                pressure->getBoundingBox(), pressure->nTensorView(), mask->nTensorView());
    }
    nTensorFieldToMesh(pressure->nTensorView(), mesh, "pressure");
    mesh.scale(param.dx);
    mesh.translate(param.physicalLocation);
    std::string fname = createFileName(param.outDir + "boat_pressure_", iIter, param.fileNamePadding) + ".vtk";
    writeVTK(mesh, fname);
}


int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);

    std::cout.precision(10);

    // Command-line arguments

    if (argc != 2 && argc != 3) {
        pcout << "Usage: " << argv[0] << " xml-input-file-name [xml-continue-file-name]" << std::endl;
        exit(1);
    }

    std::string xmlInputFileName;
    xmlInputFileName = std::string(argv[1]);
    abortIfCannotOpenFileForReading(xmlInputFileName);

    std::string xmlRestartFileName;
    bool continueSimulation = false;
    if (argc == 3) {
        xmlRestartFileName = std::string(argv[2]);
        continueSimulation = true;
    }

    int nproc = global::mpi().getSize();

    global::timer("init").start();
    global::timer("totalTime").start();

    // Set the simulation parameters.

    setParameters(xmlInputFileName);
    global::IOpolicy().activateParallelIO(param.useParallelIO);

    setDerivedParameters();

    // Setup the geometry.

    pcout << "Total number of lattice cells: " << param.nx * param.ny * param.nz << std::endl;
    MultiScalarField3D<int> tags(param.nx, param.ny, param.nz, (int) OUT);
    tags.periodicity().toggle(0, true);
    tags.periodicity().toggle(1, false);
    tags.periodicity().toggle(2, true);

    pcout << "Setting up the boat geometry." << std::endl;
    bool operateOnOuterBorder = true;
    setupBoatGeometry(operateOnOuterBorder, tags);

#ifdef PLB_DEBUG
    {
        VtkImageOutput3D<T> vtkOut(param.outDir + "tags", param.dx, param.physicalLocation);
        vtkOut.writeData<float>(*copyConvert<int,T>(tags), "tags", 1.0);
    }
#endif

    // Free-surface blocks.

#ifdef WAVE_ABSORPTION
    Dynamics<T,DESCRIPTOR>* dynamics = new WaveAbsorptionDynamics<T,DESCRIPTOR>(
            new SmagorinskyDynamics<T,DESCRIPTOR>(
                new TruncatedTRTdynamics<T,DESCRIPTOR>(param.omega, 1.1), param.omega, param.cSmago));
    param.incompressibleModel = false;
    pcout << "Dynamics: Wave Absorption Smagorinsky TRT." << std::endl;
#else
    Dynamics<T,DESCRIPTOR>* dynamics = new SmagorinskyDynamics<T,DESCRIPTOR>(
            new TruncatedTRTdynamics<T,DESCRIPTOR>(param.omega, 1.1), param.omega, param.cSmago);
    param.incompressibleModel = false;
    pcout << "Dynamics: Smagorinsky TRT." << std::endl;
#endif

    T contactAngle = (T) -1;
    plint numIBIterations = 4;
    int repelInterface = 0;
    if (param.strongRepelling) {
        repelInterface = 1; // 3;
    }
    FreeSurfaceFields3D<T,DESCRIPTOR> fields(tags.getMultiBlockManagement().getSparseBlockStructure(),
            dynamics->clone(), param.rho_LB, param.surfaceTension_LB, contactAngle, Array<T,0>::zero(),
            numIBIterations, param.vertices, param.areas, param.flags, VelFunction<T>(),
            repelInterface);

    fields.periodicityToggle(0, true);
    fields.periodicityToggle(1, false);
    fields.periodicityToggle(2, true);

    MultiTensorField3D<T,3> force((MultiBlock3D&) tags);
    force.periodicity().toggle(0, true);
    force.periodicity().toggle(1, false);
    force.periodicity().toggle(2, true);
    Array<T,3> bodyForce_LB = Array<T,3>((T) 0, -param.gravity_LB, (T) 0);
    setToConstant<T,3>(force, force.getBoundingBox(), bodyForce_LB);
    if (!util::isZero(param.A_LB)) {
        bodyForce_LB = Array<T,3>((T) 0, -param.gravity_LB + waveForce((T) 0), (T) 0);
        setToConstant<T,3>(force, param.waveDomain, bodyForce_LB);
    }

    // Initialization

    pcout << "Setting up initial condition." << std::endl;

    setToConstant(fields.flag, fields.flag.getBoundingBox(), (int) freeSurfaceFlag::empty);
    setToConstant(fields.flag, param.initialFluidDomain, (int) freeSurfaceFlag::fluid);
    setToConstant(fields.flag, tags, (int) IN, tags.getBoundingBox(), (int) freeSurfaceFlag::empty);
    setToConstant(fields.flag, param.top, (int) freeSurfaceFlag::wall);

    initializeAtEquilibrium(fields.lattice, force, fields.lattice.getBoundingBox(), initialRhoU);
    computeRhoBarJ(fields.lattice, fields.rhoBar, fields.j, fields.lattice.getBoundingBox());

    freeSurfaceAddForceToMomentum<T,DESCRIPTOR>(fields.lattice, fields.rhoBar, fields.j, fields.flag, force,
            fields.lattice.getBoundingBox());

    bool useConstRho = false;
    bool useZeroMomentum = false;
    bool initializeCell = false;
    fields.defaultInitialize(useConstRho, useZeroMomentum, initializeCell);

    Array<bool,3> reflectOnAxis;
    reflectOnAxis[0] = false;
    reflectOnAxis[1] = true;
    reflectOnAxis[2] = false;
    defineDynamics(fields.lattice, param.bottom, new SpecularReflection<T,DESCRIPTOR>(reflectOnAxis, param.rho_LB));
    setToConstant(fields.flag, param.bottom, (int) freeSurfaceFlag::slipWall);

#ifdef WAVE_ABSORPTION
    if (param.totalNumAbsorbingCells != 0) {
        setGenericExternalScalar(fields.lattice, fields.lattice.getBoundingBox(), DESCRIPTOR<T>::ExternalField::sigmaBeginsAt,
                WaveAbsorptionSigmaFunction3D<T>(fields.lattice.getBoundingBox(), param.numAbsorbingCells, param.omega));
#ifdef PLB_DEBUG
        {
            VtkImageOutput3D<T> vtkOut(param.outDir + "sigma", param.dx, param.physicalLocation);
            vtkOut.writeData<float>(*computeExternalScalar(fields.lattice, DESCRIPTOR<T>::ExternalField::sigmaBeginsAt),
                    "sigma", 1.0);
        }
#endif
        setExternalScalar(fields.lattice, fields.lattice.getBoundingBox(), DESCRIPTOR<T>::ExternalField::rhoBarBeginsAt, (T) 0);
        setExternalVector(fields.lattice, fields.lattice.getBoundingBox(), DESCRIPTOR<T>::ExternalField::uBeginsAt,
                Array<T,3>(param.inletVelocity_LB, (T) 0, (T) 0));
    }
#endif

    plint targetEnvelopeWidth = 1;
    std::auto_ptr<MultiScalarField3D<T> > targetVolumeFraction = generateMultiScalarField<T>(
            (MultiBlock3D&) fields.volumeFraction, targetEnvelopeWidth);
    copy(fields.volumeFraction, *targetVolumeFraction, fields.volumeFraction.getBoundingBox());
    std::auto_ptr<MultiTensorField3D<T,3> > targetVelocity = generateMultiTensorField<T,3>(
            (MultiBlock3D&) fields.j, targetEnvelopeWidth);
    setToConstant<T,3>(*targetVelocity, targetVelocity->getBoundingBox(), Array<T,3>(param.initialVelocity_LB, (T) 0, (T) 0));

    applyAbsorbingZones(fields, *targetVolumeFraction, *targetVelocity);

    plint iniIter = 0;

    printSimulationParameters();

    // Checkpointing.

    std::vector<MultiBlock3D*> checkpointBlocks;
    checkpointBlocks.push_back(&fields.lattice);
    checkpointBlocks.push_back(&fields.mass);
    checkpointBlocks.push_back(&fields.flag);
    checkpointBlocks.push_back(&fields.volumeFraction);
    checkpointBlocks.push_back(&fields.outsideDensity);
    checkpointBlocks.push_back(&fields.rhoBar);
    checkpointBlocks.push_back(&fields.j);
    checkpointBlocks.push_back(&tags);
    checkpointBlocks.push_back(&force);
    checkpointBlocks.push_back(targetVolumeFraction.get());
    checkpointBlocks.push_back(targetVelocity.get());

    if (continueSimulation) {
        pcout << "Reading state of the simulation from file: " << xmlRestartFileName << std::endl;
        loadState(checkpointBlocks, iniIter, param.saveDynamicContent, xmlRestartFileName);
        fields.lattice.resetTime(iniIter);
        pcout << std::endl;
    }

    // File preparation.

    FILE* fpEnergy = 0;
    if (global::mpi().isMainProcessor()) {
        std::string fileName = param.outDir + "average_energy.dat";
        fpEnergy = fopen(fileName.c_str(), continueSimulation ? "a" : "w");
        PLB_ASSERT(fpEnergy != 0);
    }

    FILE* fpForce = 0;
    if (global::mpi().isMainProcessor()) {
        std::string fileName = param.outDir + "total_force_on_boat.dat";
        fpForce = fopen(fileName.c_str(), continueSimulation ? "a" : "w");
        PLB_ASSERT(fpForce != 0);
    }

    global::mpi().barrier();
    global::timer("init").stop();

    pcout << "The full initialization phase took " << global::timer("init").getTime() << " seconds on "
          << nproc << " processes." << std::endl;









































    // Starting iterations.

    pcout << std::endl;
    pcout << "Starting simulation." << std::endl;
    pcout << std::endl;
    bool stopExecution = false;
    plint iIter = 0;
    for (iIter = iniIter; iIter < param.maxIter && !stopExecution; iIter++) {
        if (iIter != iniIter && (iIter % param.statIter == 0 || iIter == param.maxIter - 1)) {
            global::timer("io").start();
            pcout << "==== At iteration " << iIter
                  << ", t = " << iIter * param.dt << std::endl;
            T energy = freeSurfaceComputeAverageForcedEnergy(fields.lattice, force, fields.flag) *
                param.rho * (param.dx * param.dx) / (param.dt * param.dt);
            if (!util::isFiniteNumber(energy)) {
                if (global::mpi().isMainProcessor()) {
                    fprintf(stdout, "The simulation is unstable. Aborting ...\n");
                }
                global::mpi().barrier();
                exit(1);
            }
            pcout << "Average kinetic energy: " << energy << std::endl;

            bool isCompressible = !param.incompressibleModel;
            MultiNTensorField3D<T>* stress = pypal_computeStress(fields.lattice, fields.lattice.getBoundingBox(),
                    param.rho_LB, isCompressible);

            plint envelopeWidth = 1;
            std::auto_ptr<MultiScalarField3D<int> > mask = generateMultiScalarField<int>((MultiBlock3D&) fields.flag, envelopeWidth);

            setToConstant(*mask, mask->getBoundingBox(), (int) 0);
            setToConstant(*mask, fields.flag, (int) freeSurfaceFlag::fluid, fields.flag.getBoundingBox(), (int) 1);
            setToConstant(*mask, fields.flag, (int) freeSurfaceFlag::interface, fields.flag.getBoundingBox(), (int) 1);
            if (param.excludeInteriorForOutput) {
                setToConstant(*mask, tags, (int) IN, tags.getBoundingBox(), (int) 0);
            }
            for (plint iLayer = 0; iLayer < 3; iLayer++) {
                applyProcessingFunctional(new MaskedNTensorNeumannInLayersFunctional3D<T>(0, 1, -1),
                        stress->getBoundingBox(), *stress, mask->nTensorView());
            }

            Array<T,3> force_LB = surfaceForceIntegral(*stress, *param.connectedMesh);
            delete stress; stress = 0;

            T forceConversion = param.rho * util::sqr(util::sqr(param.dx)) / util::sqr(param.dt);

            pcout << "Total force on the boat: ";
            pcout << "(" << forceConversion * force_LB[0] << ", "
                         << forceConversion * force_LB[1] << ", "
                         << forceConversion * force_LB[2] << ")" << std::endl;

            if (global::mpi().isMainProcessor()) {
                fprintf(fpEnergy, "% .8e\t% .8e\n", (double) (iIter * param.dt), (double) energy);
                fflush(fpEnergy);
                fprintf(fpForce, "% .8e\t% .8e\t% .8e\t% .8e\n", (double) (iIter * param.dt),
                        (double) (forceConversion * force_LB[0]), (double) (forceConversion * force_LB[1]),
                        (double) (forceConversion * force_LB[2]));
                fflush(fpForce);
            }
            global::timer("io").stop();

            if (iIter != iniIter) {
                pcout << "Time per iteration: " << global::timer("iterations").getTime() / (double) (iIter - iniIter)
                      << std::endl;
                pcout << "Time per iteration considering the last " << param.statIter << " iterations: "
                      << global::timer("statIterations").getTime() / (double) param.statIter << std::endl;
                global::timer("statIterations").reset();
                pcout << "Total run time: " << global::timer("totalTime").getTime() << std::endl;
                pcout << "Total I/O time: " << global::timer("io").getTime() << std::endl;
            }
        }

        if (iIter % param.outIter == 0 || iIter == param.maxIter - 1) {
            global::timer("io").start();
            pcout << "== Output to disk at iteration: " << iIter << std::endl;
            writeResults(fields, force, tags, iIter);
            pcout << std::endl;
            global::timer("io").stop();
        }

        if (param.cpIter > 0 && iIter % param.cpIter == 0 && iIter != iniIter) {
            global::timer("io").start();
            pcout << "Saving the state of the simulation at iteration: " << iIter << std::endl;
            saveState(checkpointBlocks, iIter, param.saveDynamicContent, param.xmlContinueFileName,
                    param.baseFileName, param.fileNamePadding);
            pcout << std::endl;
            global::timer("io").stop();
        }

        if (iIter % param.abIter == 0) {
            stopExecution = abortExecution(param.abortFileName, checkpointBlocks, iIter,
                    param.saveDynamicContent, param.xmlContinueFileName,
                    param.baseFileName, param.fileNamePadding);

            if (stopExecution) {
                pcout << "Aborting execution at iteration: " << iIter << std::endl;
                pcout << std::endl;
            }
        }

        global::timer("iterations").start();
        global::timer("statIterations").start();
        fields.lattice.executeInternalProcessors();
        fields.lattice.evaluateStatistics();
        fields.lattice.incrementTime();

        if (!util::isZero(param.A_LB)) {
            bodyForce_LB = Array<T,3>((T) 0, -param.gravity_LB + waveForce((T) (iIter + 1)), (T) 0);
            setToConstant<T,3>(force, param.waveDomain, bodyForce_LB);
        }
        freeSurfaceAddForceToMomentum<T,DESCRIPTOR>(fields.lattice, fields.rhoBar, fields.j, fields.flag, force,
                fields.lattice.getBoundingBox());

        applyAbsorbingZones(fields, *targetVolumeFraction, *targetVelocity);

        global::timer("iterations").stop();
        global::timer("statIterations").stop();
    }

    // Solver execution statistics.

    global::mpi().barrier();
    global::timer("totalTime").stop();
    global::timer("iterations").stop();
    global::timer("io").stop();

    double totT = global::timer("totalTime").getTime();
    double iniT = global::timer("init").getTime();
    double itT  = global::timer("iterations").getTime();
    double ioT  = global::timer("io").getTime();

    double iniPC = !util::isZero(totT) ? iniT / totT * 100.0 : 0.0;
    double itPC  = !util::isZero(totT) ? itT  / totT * 100.0 : 0.0;
    double ioPC  = !util::isZero(totT) ? ioT  / totT * 100.0 : 0.0;

    pcout << "The solver finished executing successfully with " << nproc << " processes." << std::endl;
    pcout << "The total running time of the solver, was " << totT << " seconds." << std::endl;
    pcout << "The " << iIter - iniIter << " iterations of the solver, took " << itT << " seconds (" << itPC << "%)." << std::endl;
    pcout << "The total I/O time, was " << ioT << " seconds (" << ioPC << "%)." << std::endl;

    FILE *fpExecStat = 0;
    if (global::mpi().isMainProcessor()) {
        std::string fileName = param.outDir + "solver_execution_statistics.txt";
        fpExecStat = fopen(fileName.c_str(), continueSimulation ? "a" : "w");
        PLB_ASSERT(fpExecStat != 0);
        if (continueSimulation) {
            fprintf(fpExecStat, "\nSummary (after restarting) of execution of the solver: %s\n\n", argv[0]);
        } else {
            fprintf(fpExecStat, "Summary of execution of the solver: %s\n\n", argv[0]);
        }
        fprintf(fpExecStat, "Number of processes: %d\n", nproc);
        fprintf(fpExecStat, "Total execution time                             : %g s\n", totT);
        fprintf(fpExecStat, "Total time of the initialization phase           : %g s,\t%g%% of the total time.\n", iniT, iniPC);
        fprintf(fpExecStat, "Total time of the pure solution phase (no output): %g s,\t%g%% of the total time.\n", itT, itPC);
        fprintf(fpExecStat, "Total time of I/O                                : %g s,\t%g%% of the total time.\n", ioT, ioPC);
        fclose(fpExecStat);
    }

    if (global::mpi().isMainProcessor()) {
        fclose(fpForce);
        fclose(fpEnergy);
    }
    
    // Create a finish file on completion (Jack).
    // std::string emptyFile = "finished";
    // fopen(emptyFile.c_str(), "w" );
    
    delete param.connectedMesh;
    delete dynamics;

    return 0;
}
