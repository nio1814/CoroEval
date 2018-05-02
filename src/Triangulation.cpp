/* 
 * CoroEval
 *
 * An evaluation tool for coronary artery reconstructions.
 *
 * Copyright © 2014:
 *
 * Christoph Forman, Universität Erlangen-Nürnberg
 * christoph.forman@cs.fau.de
 *
 * Chris Schwemmer, Universität Erlangen-Nürnberg
 * chris.schwemmer@cs.fau.de
 *
 * Jens Wetzl, Universität Erlangen-Nürnberg
 * jens.wetzl@cs.fau.de
 *
 * CoroEval is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CoroEval is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CoroEval.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "Triangulation.h"
#include "settings.h"
#include "measurementpoint.h"

#include "wm3/Wm3Vector3.h"
#include "wm3/Wm3Matrix3.h"

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <fstream>
#include <limits>


Triangulation::Triangulation(std::vector<MeasurementPoint *> &measurementPoints)
    : m_measurementPoints(measurementPoints)
    , m_mesh(0)
{}

Triangulation::~Triangulation()
{
    if (m_mesh)
    {
        delete m_mesh;
    }
}

void Triangulation::triangulate()
{
    if (m_mesh)
    {
        delete m_mesh;
    }

    m_mesh = new Mesh();

    const size_t numPoints = m_measurementPoints.size();
    const size_t numNormalsPerPoint = Settings::getSettings().getNumNormalsAtPoint();

    // Determine point correspondences

    std::vector<size_t> offsets;
    offsets.resize(numPoints - 1);

    for (size_t i = 0; i < numPoints - 1; ++i)
    {
        // Brute force: try all offsets
        // If this becomes a bottleneck (e.g. numNormalsPerPoint is
        // increased dramatically), use the torsion of the B-spline

        size_t bestOffset;
        float bestOffsetValue = std::numeric_limits<float>::max();

        for (size_t offset = 0; offset < 2 * numNormalsPerPoint; ++offset)
        {
            float combinedDistance = 0.0f;

            for (size_t j = 0; j < 2 * numNormalsPerPoint; ++j)
            {
                Wm3::Vector3d &pt1 = m_surfacePoints[i][j];
                Wm3::Vector3d &pt2 = m_surfacePoints[i+1][(j+offset) % (2 * numNormalsPerPoint)];

                combinedDistance += (pt1 - pt2).SquaredLength();
            }

            if (bestOffsetValue > combinedDistance)
            {
                bestOffsetValue = combinedDistance;
                bestOffset = offset;
            }
        }

        offsets[i] = bestOffset;
    }

    // Assemble the mesh
    std::vector<Mesh::VertexHandle> face_vhandles, previous_vhandles, current_vhandles;
    previous_vhandles.resize(2 * numNormalsPerPoint);
    current_vhandles.resize(2 * numNormalsPerPoint);

    // Initialize for first iteration
    for (size_t j = 0; j < 2 * numNormalsPerPoint; ++j)
    {
        Wm3::Vector3d &v = m_surfacePoints[0][j];
        current_vhandles[j] = m_mesh->add_vertex(Mesh::Point(v[0], v[1], v[2]));
    }

    for (size_t i = 0; i < numPoints - 1; ++i)
    {
        // Copy current to old handles and load new ones
        for (size_t j = 0; j < 2 * numNormalsPerPoint; ++j)
        {
            previous_vhandles[j] = current_vhandles[j];
            Wm3::Vector3d &v = m_surfacePoints[i+1][j];
            current_vhandles[j] = m_mesh->add_vertex(Mesh::Point(v[0], v[1], v[2]));
        }

        // Create tris for current point pair i and i+1
        for (size_t j = 0; j < 2 * numNormalsPerPoint; ++j)
        {
            Mesh::VertexHandle &a = previous_vhandles[j];
            Mesh::VertexHandle &b =  current_vhandles[(j + offsets[i]) % (2 * numNormalsPerPoint)];
            Mesh::VertexHandle &c = previous_vhandles[(j + 1) % (2 * numNormalsPerPoint)];
            Mesh::VertexHandle &d =  current_vhandles[(j + offsets[i] + 1) % (2 * numNormalsPerPoint)];

            face_vhandles.clear();
            face_vhandles.push_back(a);
            face_vhandles.push_back(b);
            face_vhandles.push_back(d);
            m_mesh->add_face(face_vhandles);

            face_vhandles.clear();
            face_vhandles.push_back(a);
            face_vhandles.push_back(d);
            face_vhandles.push_back(c);
            m_mesh->add_face(face_vhandles);
        }
    }
}

bool Triangulation::saveMesh(const std::string &fileName)
{
    this->extractSurfacePoints();
    this->triangulate();

    try
    {
        if (!OpenMesh::IO::write_mesh(*m_mesh, fileName))
        {
            std::cerr << "Cannot write mesh to file '" << fileName << "'" << std::endl;
            return false;
        }
    }
    catch (std::exception &x)
    {
        std::cerr << x.what() << std::endl;
        return false;
    }

    return true;
}

bool Triangulation::saveTubeMesh(const std::string &fileName)
{
    this->extractTubePoints();
    this->triangulate();

    try
    {
        if (!OpenMesh::IO::write_mesh(*m_mesh, fileName))
        {
            std::cerr << "Cannot write mesh to file '" << fileName << "'" << std::endl;
            return false;
        }
    }
    catch (std::exception &x)
    {
        std::cerr << x.what() << std::endl;
        return false;
    }

    return true;
}

bool Triangulation::saveSurfacePoints(const std::string &fileName)
{
    this->extractSurfacePoints();

    std::ofstream outputStream(fileName.c_str());
    if (!outputStream.good())
    {
        outputStream.close();
        return false;
    }

    for (size_t i = 0; i < m_surfacePoints.size(); ++i)
    {
        for (size_t j = 0; j < m_surfacePoints[i].size(); ++j)
        {
            outputStream << m_surfacePoints[i][j] << std::endl;
        }
    }

    outputStream.close();
    return true;
}

void Triangulation::extractSurfacePoints()
{
    const size_t numPoints = m_measurementPoints.size();
    const size_t numNormalsPerPoint = Settings::getSettings().getNumNormalsAtPoint();

    // Save points on vessel circumference
    m_surfacePoints.clear();
    m_surfacePoints.resize(numPoints);

    std::vector<bool> outlierStatus;
    outlierStatus.resize(2 * numNormalsPerPoint);

    for (size_t i = 0; i < m_measurementPoints.size(); ++i)
    {
        MeasurementPoint *mp = m_measurementPoints[i];
        mp->updateDiameterOutlierStatus();

        // Extract points and outlier status
        for (size_t j = 0; j < numNormalsPerPoint; ++j)
        {
            m_surfacePoints[i].push_back(mp->get50PointInVolume(j, Interval::Left));
            outlierStatus[j] = mp->isPointDiameterOutlier(j, Interval::Left);
        }

        for (size_t j = 0; j < numNormalsPerPoint; ++j)
        {
            m_surfacePoints[i].push_back(mp->get50PointInVolume(j, Interval::Right));
            outlierStatus[numNormalsPerPoint + j] = mp->isPointDiameterOutlier(j, Interval::Right);
        }

        size_t numOutliers = 0;
        for (size_t j = 0; j < 2 * numNormalsPerPoint; ++j)
            numOutliers += static_cast<size_t>(outlierStatus[j]);
        std::cout << i << ": " << numOutliers << " outliers" << std::endl;

        // Replace outliers
        for (size_t j = 0; j < 2 * numNormalsPerPoint; ++j)
        {
            if (!outlierStatus[j])
                continue;

            // If two points before and two points after this aren't outliers,
            // least-squares fit a parabola to replace the outlier
            if (!outlierStatus[(j+1) % (2 * numNormalsPerPoint)] &&
                !outlierStatus[(j+2) % (2 * numNormalsPerPoint)] &&
                !outlierStatus[(j-1+2*numNormalsPerPoint) % (2 * numNormalsPerPoint)] &&
                !outlierStatus[(j-2+2*numNormalsPerPoint) % (2 * numNormalsPerPoint)])
            {
                // Rotate points so that the point to be replaced lies on the y axis
                Wm3::Vector3d p  = mp->getPositionInVolume();
                Wm3::Vector3d x0 = m_surfacePoints[i][(j+1) % (2 * numNormalsPerPoint)] - p;
                Wm3::Vector3d x1 = m_surfacePoints[i][(j+2) % (2 * numNormalsPerPoint)] - p;
                Wm3::Vector3d x2 = m_surfacePoints[i][(j-1+2*numNormalsPerPoint) % (2 * numNormalsPerPoint)] - p;
                Wm3::Vector3d x3 = m_surfacePoints[i][(j-2+2*numNormalsPerPoint) % (2 * numNormalsPerPoint)] - p;

                Wm3::Vector3d u, v;
                mp->getPlane(u, v);

                double angle = 0.0;
                if (std::abs(v.Dot(mp->getSamplingDirection(j%numNormalsPerPoint) / mp->getSamplingDirection(j%numNormalsPerPoint).Length() / v.Length())) < 1.0)
                    angle = std::acos(v.Dot(mp->getSamplingDirection(j%numNormalsPerPoint)) / mp->getSamplingDirection(j%numNormalsPerPoint).Length() / v.Length());

                Wm3::Matrix3d mat(mp->getTangent(), angle);

                if (std::abs(std::abs(v.Dot(mat * mp->getSamplingDirection(j%numNormalsPerPoint))) - 1.0) > 0.001)
                {
                    angle = -angle;
                    mat = Wm3::Matrix3d(mp->getTangent(), angle);
                }

                x0 = mat * x0;
                x1 = mat * x1;
                x2 = mat * x2;
                x3 = mat * x3;
                Wm3::Vector3d xorig = mat * (m_surfacePoints[i][j] - p);

                double x[4];
                double y[4];

                x[0] = u.Dot(x0);
                y[0] = v.Dot(x0);
                x[1] = u.Dot(x1);
                y[1] = v.Dot(x1);
                x[2] = u.Dot(x2);
                y[2] = v.Dot(x2);
                x[3] = u.Dot(x3);
                y[3] = v.Dot(x3);

                double x0pow2 = x[0] * x[0];
                double x1pow2 = x[1] * x[1];
                double x2pow2 = x[2] * x[2];
                double x3pow2 = x[3] * x[3];

                double x0pow3 = x[0] * x[0] * x[0];
                double x1pow3 = x[1] * x[1] * x[1];
                double x2pow3 = x[2] * x[2] * x[2];
                double x3pow3 = x[3] * x[3] * x[3];

                double x0pow4 = x[0] * x[0] * x[0] * x[0];
                double x1pow4 = x[1] * x[1] * x[1] * x[1];
                double x2pow4 = x[2] * x[2] * x[2] * x[2];
                double x3pow4 = x[3] * x[3] * x[3] * x[3];

                // Code generated with Maple
//                double coeffa = -(-x[1] * x2pow3 * y[3] + 2.0 * x0pow2 * x3pow2 * y[1] + 2.0 * x0pow3 * x[2] * y[2] - 2.0 * x0pow2 * x2pow2 * y[2] - 2.0 * x1pow2 * x2pow2 * y[2] - 2.0 * x3pow2 * x2pow2 * y[2] + 2.0 * x0pow3 * x[1] * y[1] - 2.0 * x0pow2 * x3pow2 * y[0] - 2.0 * x0pow2 * x2pow2 * y[0] - 2.0 * x0pow2 * x1pow2 * y[0] + 2.0 * x[0] * x3pow3 * y[0] - 2.0 * x0pow2 * x1pow2 * y[1] + 2.0 * x2pow3 * x[3] * y[3] + 2.0 * x[2] * x3pow3 * y[2] + 2.0 * x1pow3 * x[2] * y[2] + 2.0 * x1pow2 * x3pow2 * y[0] + 2.0 * x1pow3 * x[3] * y[3] - 2.0 * x1pow2 * x3pow2 * y[3] - 2.0 * x0pow2 * x3pow2 * y[3] + 2.0 * x[1] * x2pow3 * y[1] - 2.0 * x1pow2 * x3pow2 * y[1] + 2.0 * x0pow2 * x3pow2 * y[2] + 2.0 * x3pow2 * x2pow2 * y[1] + 2.0 * x[0] * x1pow3 * y[0] + 2.0 * x0pow2 * x2pow2 * y[1] + 2.0 * x0pow2 * x1pow2 * y[3] - 2.0 * x3pow2 * x2pow2 * y[3] + 2.0 * x3pow2 * x2pow2 * y[0] + 2.0 * x0pow3 * x[3] * y[3] + 2.0 * x[0] * x2pow3 * y[0] + 2.0 * x0pow2 * x1pow2 * y[2] + 2.0 * x1pow2 * x2pow2 * y[3] + 2.0 * x1pow2 * x3pow2 * y[2] + 2.0 * x1pow2 * x2pow2 * y[0] + 2.0 * x0pow2 * x2pow2 * y[3] - 2.0 * x1pow2 * x2pow2 * y[1] + 2.0 * x[1] * x3pow3 * y[1] - x0pow3 * x[3] * y[1] - x1pow3 * x[3] * y[2] - x1pow3 * x[3] * y[0] - x[1] * x3pow3 * y[0] - x2pow3 * x[3] * y[0] - x1pow3 * x[2] * y[3] - x1pow3 * x[2] * y[0] - x0pow3 * x[1] * y[2] - x0pow3 * x[1] * y[3] - x[1] * x2pow3 * y[0] - x[2] * x3pow3 * y[0] - x0pow3 * x[3] * y[2] - x[0] * x1pow3 * y[3] - x[0] * x1pow3 * y[2] - x[1] * x3pow3 * y[2] - x[2] * x3pow3 * y[1] - x2pow3 * x[3] * y[1] - x[1] * x[2] * x0pow2 * y[2] - x[0] * x[2] * x3pow2 * y[0] - x[0] * x[2] * x3pow2 * y[2] - x[0] * x[2] * x1pow2 * y[0] - x[0] * x[2] * x1pow2 * y[2] - x[0] * x[3] * x2pow2 * y[0] - x[2] * x[3] * x0pow2 * y[2] - x[2] * x[3] * x0pow2 * y[3] - x[2] * x[3] * x1pow2 * y[2] - x[2] * x[3] * x1pow2 * y[3] - x[0] * x[3] * x2pow2 * y[3] - x[0] * x[3] * x1pow2 * y[0] - x[0] * x[3] * x1pow2 * y[3] - x[1] * x[3] * x2pow2 * y[1] - x[1] * x[3] * x2pow2 * y[3] - x[1] * x[3] * x0pow2 * y[1] - x[1] * x[3] * x0pow2 * y[3] - x[0] * x[1] * x2pow2 * y[0] - x[0] * x[1] * x2pow2 * y[1] - x[0] * x[1] * x3pow2 * y[0] - x[0] * x[1] * x3pow2 * y[1] - x[1] * x[2] * x3pow2 * y[1] - x[1] * x[2] * x3pow2 * y[2] - x[1] * x[2] * x0pow2 * y[1] - x[0] * x2pow3 * y[1] - x[0] * x2pow3 * y[3] - x0pow3 * x[2] * y[1] - x0pow3 * x[2] * y[3] - x[0] * x3pow3 * y[1] - x[0] * x3pow3 * y[2] + 2.0 * x[2] * x[3] * x0pow2 * y[0] + 2.0 * x[2] * x[3] * x1pow2 * y[1] + 2.0 * x[0] * x[2] * x1pow2 * y[1] + 2.0 * x[0] * x[2] * x3pow2 * y[3] + 2.0 * x[0] * x[3] * x1pow2 * y[1] + 2.0 * x[0] * x[3] * x2pow2 * y[2] + 2.0 * x[1] * x[3] * x0pow2 * y[0] + 2.0 * x[1] * x[3] * x2pow2 * y[2] + 2.0 * x[0] * x[1] * x2pow2 * y[2] + 2.0 * x[0] * x[1] * x3pow2 * y[3] + 2.0 * x[1] * x[2] * x0pow2 * y[0] + 2.0 * x[1] * x[2] * x3pow2 * y[3]) / (x3pow4 * x2pow2 + x1pow4 * x3pow2 + x1pow4 * x2pow2 + x3pow2 * x2pow4 + x1pow2 * x2pow4 + x1pow2 * x3pow4 + x0pow4 * x1pow2 + x0pow4 * x2pow2 + x0pow4 * x3pow2 + x0pow2 * x1pow4 + x0pow2 * x2pow4 + x0pow2 * x3pow4 - 2.0 * x2pow3 * x3pow3 - 2.0 * x0pow3 * x1pow3 - 2.0 * x1pow3 * x2pow3 - 2.0 * x0pow3 * x3pow3 - 2.0 * x1pow3 * x3pow3 - 2.0 * x0pow3 * x2pow3 - 3.0 * x1pow2 * x2pow2 * x3pow2 - 3.0 * x0pow2 * x1pow2 * x3pow2 - 3.0 * x0pow2 * x1pow2 * x2pow2 - 3.0 * x0pow2 * x2pow2 * x3pow2 + x1pow3 * x[2] * x3pow2 + x[0] * x1pow3 * x2pow2 + x[1] * x3pow3 * x2pow2 - x[0] * x[3] * x1pow4 - x[0] * x[3] * x2pow4 + x[0] * x3pow3 * x1pow2 + x0pow3 * x[1] * x2pow2 + x0pow3 * x[1] * x3pow2 + x0pow3 * x[2] * x1pow2 + x0pow3 * x[2] * x3pow2 + x0pow3 * x[3] * x1pow2 + x0pow3 * x[3] * x2pow2 + x2pow3 * x[3] * x1pow2 + x[2] * x3pow3 * x0pow2 - x[2] * x[3] * x0pow4 - x[2] * x[3] * x1pow4 + x2pow3 * x[3] * x0pow2 + x[2] * x3pow3 * x1pow2 + x[0] * x2pow3 * x1pow2 - x[0] * x[2] * x3pow4 + x[0] * x2pow3 * x3pow2 - x[0] * x[2] * x1pow4 + x1pow3 * x[2] * x0pow2 + x[1] * x2pow3 * x0pow2 - x[1] * x[2] * x0pow4 + x[1] * x2pow3 * x3pow2 + x[0] * x1pow3 * x3pow2 + x[1] * x3pow3 * x0pow2 - x[0] * x[1] * x2pow4 - x[1] * x[2] * x3pow4 - x[1] * x[3] * x2pow4 + x1pow3 * x[3] * x2pow2 - x[1] * x[3] * x0pow4 + x1pow3 * x[3] * x0pow2 - x[0] * x[1] * x3pow4 + x[0] * x3pow3 * x2pow2) / 2.0;
//                double coeffb = (2.0 * x1pow4 * x[3] * y[3] - 2.0 * x0pow3 * x3pow2 * y[3] - 2.0 * x1pow2 * x3pow3 * y[1] + 2.0 * x0pow4 * x[3] * y[3] - 2.0 * x1pow3 * x2pow2 * y[2] + 2.0 * x[1] * x2pow4 * y[1] - 2.0 * x0pow3 * x1pow2 * y[1] + 2.0 * x[0] * x2pow4 * y[0] + 2.0 * x0pow4 * x[1] * y[1] - 2.0 * x0pow2 * x3pow3 * y[0] + 2.0 * x1pow4 * x[2] * y[2] + 2.0 * x[0] * x3pow4 * y[0] + 2.0 * x0pow4 * x[2] * y[2] + 2.0 * x[2] * x3pow4 * y[2] - 2.0 * x0pow2 * x2pow3 * y[0] + 2.0 * x[1] * x3pow4 * y[1] - 2.0 * x1pow3 * x3pow2 * y[3] - 2.0 * x0pow3 * x2pow2 * y[2] - 2.0 * x2pow2 * x3pow3 * y[2] - 2.0 * x0pow2 * x1pow3 * y[0] - 2.0 * x1pow2 * x2pow3 * y[1] + 2.0 * x2pow4 * x[3] * y[3] - 2.0 * x2pow3 * x3pow2 * y[3] + 2.0 * x[0] * x1pow4 * y[0] - 2.0 * x0pow2 * x1pow2 * x[3] * y[3] - 2.0 * x0pow2 * x[1] * x2pow2 * y[1] - 2.0 * x0pow2 * x[1] * x3pow2 * y[1] - 2.0 * x0pow2 * x2pow2 * x[3] * y[3] - 2.0 * x0pow2 * x[2] * x3pow2 * y[2] - 2.0 * x[0] * x1pow2 * x2pow2 * y[0] - 2.0 * x[0] * x1pow2 * x3pow2 * y[0] - 2.0 * x[0] * x2pow2 * x3pow2 * y[0] - 2.0 * x1pow2 * x2pow2 * x[3] * y[3] - x0pow4 * x[1] * y[2] + x0pow3 * x1pow2 * y[2] + x0pow3 * x2pow2 * y[1] - 2.0 * x1pow2 * x[2] * x3pow2 * y[2] - 2.0 * x[1] * x2pow2 * x3pow2 * y[1] - 2.0 * x0pow2 * x1pow2 * x[2] * y[2] + x0pow2 * x1pow3 * y[2] + x0pow3 * x3pow2 * y[2] - x0pow4 * x[2] * y[3] - x[1] * x2pow4 * y[3] - x[1] * x3pow4 * y[0] - x[1] * x3pow4 * y[2] - x2pow4 * x[3] * y[0] - x2pow4 * x[3] * y[1] + x2pow3 * x3pow2 * y[0] + x2pow3 * x3pow2 * y[1] + x2pow2 * x3pow3 * y[0] + x2pow2 * x3pow3 * y[1] - x[2] * x3pow4 * y[0] - x[2] * x3pow4 * y[1] + x1pow2 * x2pow3 * y[0] + x1pow2 * x2pow3 * y[3] + x1pow2 * x3pow3 * y[0] + x1pow2 * x3pow3 * y[2] - x[1] * x2pow4 * y[0] + x0pow2 * x3pow3 * y[1] + x0pow2 * x3pow3 * y[2] - x[0] * x1pow4 * y[2] - x[0] * x1pow4 * y[3] - x[0] * x2pow4 * y[1] - x[0] * x2pow4 * y[3] - x[0] * x3pow4 * y[1] - x[0] * x3pow4 * y[2] - x1pow4 * x[2] * y[0] - x1pow4 * x[2] * y[3] - x1pow4 * x[3] * y[0] - x1pow4 * x[3] * y[2] + x1pow3 * x2pow2 * y[0] + x1pow3 * x2pow2 * y[3] + x1pow3 * x3pow2 * y[0] + x1pow3 * x3pow2 * y[2] + x0pow2 * x1pow3 * y[3] - x0pow4 * x[1] * y[3] + x1pow2 * x[2] * x3pow2 * y[1] + x1pow2 * x[2] * x3pow2 * y[3] + x[1] * x2pow2 * x3pow2 * y[2] + x[1] * x2pow2 * x3pow2 * y[3] + x0pow2 * x1pow2 * x[2] * y[0] + x0pow2 * x1pow2 * x[2] * y[1] + x0pow2 * x1pow2 * x[3] * y[0] + x0pow2 * x1pow2 * x[3] * y[1] + x0pow2 * x[1] * x2pow2 * y[0] + x0pow2 * x[1] * x2pow2 * y[2] + x0pow2 * x[1] * x3pow2 * y[0] + x0pow2 * x[1] * x3pow2 * y[3] + x0pow2 * x2pow2 * x[3] * y[0] + x0pow2 * x2pow2 * x[3] * y[2] + x0pow2 * x[2] * x3pow2 * y[0] + x0pow2 * x[2] * x3pow2 * y[3] + x[0] * x1pow2 * x2pow2 * y[1] + x[0] * x1pow2 * x2pow2 * y[2] + x[0] * x1pow2 * x3pow2 * y[1] + x[0] * x1pow2 * x3pow2 * y[3] + x[0] * x2pow2 * x3pow2 * y[2] + x[0] * x2pow2 * x3pow2 * y[3] + x1pow2 * x2pow2 * x[3] * y[1] + x1pow2 * x2pow2 * x[3] * y[2] + x0pow2 * x2pow3 * y[1] - x0pow4 * x[2] * y[1] + x0pow3 * x2pow2 * y[3] + x0pow3 * x1pow2 * y[3] - x0pow4 * x[3] * y[2] + x0pow3 * x3pow2 * y[1] - x0pow4 * x[3] * y[1] + x0pow2 * x2pow3 * y[3]) / (x3pow4 * x2pow2 + x1pow4 * x3pow2 + x1pow4 * x2pow2 + x3pow2 * x2pow4 + x1pow2 * x2pow4 + x1pow2 * x3pow4 + x0pow4 * x1pow2 + x0pow4 * x2pow2 + x0pow4 * x3pow2 + x0pow2 * x1pow4 + x0pow2 * x2pow4 + x0pow2 * x3pow4 - 2.0 * x2pow3 * x3pow3 - 2.0 * x0pow3 * x1pow3 - 2.0 * x1pow3 * x2pow3 - 2.0 * x0pow3 * x3pow3 - 2.0 * x1pow3 * x3pow3 - 2.0 * x0pow3 * x2pow3 - 0.3e1 * x1pow2 * x2pow2 * x3pow2 - 0.3e1 * x0pow2 * x1pow2 * x3pow2 - 0.3e1 * x0pow2 * x1pow2 * x2pow2 - 0.3e1 * x0pow2 * x2pow2 * x3pow2 + x1pow3 * x[2] * x3pow2 + x[0] * x1pow3 * x2pow2 + x[1] * x3pow3 * x2pow2 - x[0] * x[3] * x1pow4 - x[0] * x[3] * x2pow4 + x[0] * x3pow3 * x1pow2 + x0pow3 * x[1] * x2pow2 + x0pow3 * x[1] * x3pow2 + x0pow3 * x[2] * x1pow2 + x0pow3 * x[2] * x3pow2 + x0pow3 * x[3] * x1pow2 + x0pow3 * x[3] * x2pow2 + x2pow3 * x[3] * x1pow2 + x[2] * x3pow3 * x0pow2 - x[2] * x[3] * x0pow4 - x[2] * x[3] * x1pow4 + x2pow3 * x[3] * x0pow2 + x[2] * x3pow3 * x1pow2 + x[0] * x2pow3 * x1pow2 - x[0] * x[2] * x3pow4 + x[0] * x2pow3 * x3pow2 - x[0] * x[2] * x1pow4 + x1pow3 * x[2] * x0pow2 + x[1] * x2pow3 * x0pow2 - x[1] * x[2] * x0pow4 + x[1] * x2pow3 * x3pow2 + x[0] * x1pow3 * x3pow2 + x[1] * x3pow3 * x0pow2 - x[0] * x[1] * x2pow4 - x[1] * x[2] * x3pow4 - x[1] * x[3] * x2pow4 + x1pow3 * x[3] * x2pow2 - x[1] * x[3] * x0pow4 + x1pow3 * x[3] * x0pow2 - x[0] * x[1] * x3pow4 + x[0] * x3pow3 * x2pow2) / 2.0;
                double coeffc = (x2pow2 * x[0] * x3pow3 * y[0] + x2pow2 * x[0] * x3pow3 * y[2] - x3pow4 * x[0] * x[1] * y[0] - x3pow4 * x[0] * x[1] * y[1] - x3pow4 * x[1] * x[2] * y[1] - x3pow4 * x[1] * x[2] * y[2] - x3pow4 * x[0] * x[2] * y[0] + x0pow2 * x[1] * x3pow3 * y[1] + x0pow2 * x1pow3 * x[3] * y[0] + x0pow2 * x1pow3 * x[3] * y[3] + x0pow2 * x[1] * x2pow3 * y[0] + x3pow4 * x0pow2 * y[1] + x3pow4 * x0pow2 * y[2] + x3pow4 * x1pow2 * y[0] + x3pow4 * x1pow2 * y[2] + x3pow4 * x2pow2 * y[0] + x3pow4 * x2pow2 * y[1] + x0pow4 * x2pow2 * y[1] + x0pow4 * x2pow2 * y[3] + x0pow4 * x3pow2 * y[1] + x0pow2 * x[1] * x2pow3 * y[1] + x0pow2 * x1pow3 * x[2] * y[0] + x0pow2 * x1pow3 * x[2] * y[2] + x0pow2 * x2pow3 * x[3] * y[0] - x1pow4 * x[2] * x[3] * y[2] - x1pow4 * x[2] * x[3] * y[3] - x1pow4 * x[0] * x[3] * y[0] - x1pow4 * x[0] * x[3] * y[3] + x1pow3 * x[3] * x2pow2 * y[2] + x1pow3 * x[3] * x2pow2 * y[3] + x1pow3 * x[0] * x2pow2 * y[2] + x1pow3 * x[0] * x3pow2 * y[3] + x1pow3 * x[0] * x2pow2 * y[0] + x1pow3 * x[0] * x3pow2 * y[0] + x1pow3 * x[2] * x3pow2 * y[3] + x1pow3 * x[2] * x3pow2 * y[2] - x1pow4 * x[0] * x[2] * y[0] - x1pow4 * x[0] * x[2] * y[2] + x1pow2 * x2pow3 * x[3] * y[1] + x1pow2 * x2pow3 * x[3] * y[3] + x1pow2 * x[2] * x3pow3 * y[1] + x1pow2 * x[2] * x3pow3 * y[2] - x0pow4 * x[2] * x[3] * y[2] - x0pow4 * x[2] * x[3] * y[3] + x0pow3 * x[3] * x2pow2 * y[3] + x0pow3 * x[3] * x1pow2 * y[3] - x0pow4 * x[1] * x[3] * y[1] - x0pow4 * x[1] * x[3] * y[3] + x0pow3 * x[1] * x2pow2 * y[2] + x0pow3 * x[1] * x3pow2 * y[3] + x0pow3 * x[1] * x2pow2 * y[1] + x0pow3 * x[1] * x3pow2 * y[1] - x0pow4 * x[1] * x[2] * y[1] - x0pow4 * x[1] * x[2] * y[2] + x0pow3 * x[2] * x1pow2 * y[1] + x0pow3 * x[2] * x3pow2 * y[3] + x0pow3 * x[2] * x3pow2 * y[2] + x0pow3 * x[2] * x1pow2 * y[2] + x0pow3 * x[3] * x1pow2 * y[1] + x0pow3 * x[3] * x2pow2 * y[2] + x0pow2 * x[1] * x3pow3 * y[0] + x1pow4 * x0pow2 * y[2] + x1pow4 * x0pow2 * y[3] + x1pow4 * x2pow2 * y[3] + x2pow4 * x1pow2 * y[0] + x2pow4 * x1pow2 * y[3] + x2pow4 * x3pow2 * y[0] + x2pow4 * x3pow2 * y[1] - 0.2e1 * x1pow3 * x2pow3 * y[0] - 0.2e1 * x0pow3 * x3pow3 * y[1] - 0.2e1 * x1pow3 * x3pow3 * y[2] - 0.2e1 * x0pow3 * x3pow3 * y[2] - 0.2e1 * x0pow3 * x1pow3 * y[3] - 0.2e1 * x2pow3 * x3pow3 * y[0] - 0.2e1 * x0pow3 * x2pow3 * y[1] - 0.2e1 * x0pow3 * x2pow3 * y[3] - 0.2e1 * x2pow3 * x3pow3 * y[1] - 0.2e1 * x1pow3 * x3pow3 * y[0] - 0.2e1 * x1pow3 * x2pow3 * y[3] - 0.2e1 * x0pow3 * x1pow3 * y[2] - 0.2e1 * x0pow2 * x3pow2 * x2pow2 * y[0] - 0.2e1 * x0pow2 * x3pow2 * x2pow2 * y[3] - 0.2e1 * x1pow2 * x3pow2 * x2pow2 * y[2] - 0.2e1 * x1pow2 * x3pow2 * x2pow2 * y[1] - 0.2e1 * x1pow2 * x3pow2 * x2pow2 * y[3] - 0.2e1 * x0pow2 * x1pow2 * x2pow2 * y[2] - 0.2e1 * x0pow2 * x1pow2 * x3pow2 * y[3] - 0.2e1 * x0pow2 * x1pow2 * x2pow2 * y[0] - 0.2e1 * x0pow2 * x1pow2 * x2pow2 * y[1] - 0.2e1 * x0pow2 * x1pow2 * x3pow2 * y[0] - 0.2e1 * x0pow2 * x1pow2 * x3pow2 * y[1] - 0.2e1 * x0pow2 * x3pow2 * x2pow2 * y[2] + x0pow2 * x2pow3 * x[3] * y[3] + x0pow2 * x[2] * x3pow3 * y[0] + x0pow2 * x[2] * x3pow3 * y[2] + x2pow4 * x0pow2 * y[1] + x1pow4 * x2pow2 * y[0] + x1pow4 * x3pow2 * y[0] + x1pow4 * x3pow2 * y[2] + x2pow4 * x0pow2 * y[3] + x0pow4 * x3pow2 * y[2] + x0pow4 * x1pow2 * y[2] + x0pow4 * x1pow2 * y[3] - x3pow4 * x[0] * x[2] * y[2] + x1pow2 * x[0] * x2pow3 * y[0] + x1pow2 * x[0] * x2pow3 * y[1] + x1pow2 * x[0] * x3pow3 * y[0] + x1pow2 * x[0] * x3pow3 * y[1] - x2pow4 * x[0] * x[3] * y[3] - x2pow4 * x[1] * x[3] * y[1] - x2pow4 * x[1] * x[3] * y[3] - x2pow4 * x[0] * x[1] * y[0] - x2pow4 * x[0] * x[1] * y[1] + x2pow3 * x[1] * x3pow2 * y[3] + x2pow3 * x[1] * x3pow2 * y[1] + x2pow3 * x[0] * x3pow2 * y[3] + x2pow3 * x[0] * x3pow2 * y[0] - x2pow4 * x[0] * x[3] * y[0] + x2pow2 * x[1] * x3pow3 * y[1] + x2pow2 * x[1] * x3pow3 * y[2]) / (x3pow4 * x2pow2 + x1pow4 * x3pow2 + x1pow4 * x2pow2 + x3pow2 * x2pow4 + x1pow2 * x2pow4 + x1pow2 * x3pow4 + x0pow4 * x1pow2 + x0pow4 * x2pow2 + x0pow4 * x3pow2 + x0pow2 * x1pow4 + x0pow2 * x2pow4 + x0pow2 * x3pow4 - 0.2e1 * x2pow3 * x3pow3 - 0.2e1 * x0pow3 * x1pow3 - 0.2e1 * x1pow3 * x2pow3 - 0.2e1 * x0pow3 * x3pow3 - 0.2e1 * x1pow3 * x3pow3 - 0.2e1 * x0pow3 * x2pow3 - 0.3e1 * x1pow2 * x2pow2 * x3pow2 - 0.3e1 * x0pow2 * x1pow2 * x3pow2 - 0.3e1 * x0pow2 * x1pow2 * x2pow2 - 0.3e1 * x0pow2 * x2pow2 * x3pow2 + x1pow3 * x[2] * x3pow2 + x[0] * x1pow3 * x2pow2 + x[1] * x3pow3 * x2pow2 - x[0] * x[3] * x1pow4 - x[0] * x[3] * x2pow4 + x[0] * x3pow3 * x1pow2 + x0pow3 * x[1] * x2pow2 + x0pow3 * x[1] * x3pow2 + x0pow3 * x[2] * x1pow2 + x0pow3 * x[2] * x3pow2 + x0pow3 * x[3] * x1pow2 + x0pow3 * x[3] * x2pow2 + x2pow3 * x[3] * x1pow2 + x[2] * x3pow3 * x0pow2 - x[2] * x[3] * x0pow4 - x[2] * x[3] * x1pow4 + x2pow3 * x[3] * x0pow2 + x[2] * x3pow3 * x1pow2 + x[0] * x2pow3 * x1pow2 - x[0] * x[2] * x3pow4 + x[0] * x2pow3 * x3pow2 - x[0] * x[2] * x1pow4 + x1pow3 * x[2] * x0pow2 + x[1] * x2pow3 * x0pow2 - x[1] * x[2] * x0pow4 + x[1] * x2pow3 * x3pow2 + x[0] * x1pow3 * x3pow2 + x[1] * x3pow3 * x0pow2 - x[0] * x[1] * x2pow4 - x[1] * x[2] * x3pow4 - x[1] * x[3] * x2pow4 + x1pow3 * x[3] * x2pow2 - x[1] * x[3] * x0pow4 + x1pow3 * x[3] * x0pow2 - x[0] * x[1] * x3pow4 + x[0] * x3pow3 * x2pow2) / 2.f;

                // We only need the constant parameter of the parabola because
                // we just want to evaluate it at x=0

                m_surfacePoints[i][j] = Wm3::Matrix3d(mp->getTangent(), -angle) * (coeffc * v) + mp->getPositionInVolume();

                // Debug code: generate TikZ visualization of parabola fit
//                for (size_t k = 0; k < 2 * numNormalsPerPoint; ++k)
//                {
//                    if (k == j) continue;
//                    std::cout << "\\draw[fill] (" << u.Dot(mat * (m_surfacePoints[i][k] - p)) << ", " << v.Dot(mat * (m_surfacePoints[i][k] - p)) << ") circle (2pt);" << std::endl;
//                }
//                std::cout << "\\draw[blue] (-2, " << 4 * coeffa - 2 * coeffb + coeffc << ")";
//                for (double k = -2.0; k <= 2.0; k += 0.1)
//                {
//                    std::cout << " -- (" << k << ", " << k*k*coeffa + k*coeffb + coeffc << ")";
//                }
//                std::cout << ";" << std::endl;
//                std::cout << "\\draw[fill=red] (0, " << v.Dot(xorig) << ") circle (2pt);" << std::endl;
//                std::cout << "\\draw[fill=green] (0, " << coeffc << ") circle (2pt);" << std::endl;

//                std::cout << m_surfacePoints[i][j] << std::endl;
//                for (size_t k = 0; k < 2 * numNormalsPerPoint; ++k)
//                {
//                    if (j == k) continue;
//                    std::cout << m_surfacePoints[i][k] << std::endl;
//                }

            }
            else // Replace with circle fit
            {
                // Compute intersections between fitted circle and sampling normal

                Wm3::Vector3d u, v;
                mp->getPlane(u, v);

                double mx = -mp->getImprovedPositionInPlane().X();
                double my = -mp->getImprovedPositionInPlane().Y();
                double sx = u.Dot(mp->getSamplingDirection(j % numNormalsPerPoint));
                double sy = v.Dot(mp->getSamplingDirection(j % numNormalsPerPoint));
                double r  = mp->getDiameterVoxels() / 2.0;

                double tau1 = (-(mx * sx) - (my * sy) + sqrt((2 * mx * sx * my * sy - sy * sy * mx * mx + sy * sy * r * r + sx * sx * r * r - sx * sx * my * my))) / (sy * sy + sx * sx);
                double tau2 = (-(mx * sx) - (my * sy) - sqrt((2 * mx * sx * my * sy - sy * sy * mx * mx + sy * sy * r * r + sx * sx * r * r - sx * sx * my * my))) / (sy * sy + sx * sx);

                Wm3::Vector3d p1 = mp->getPositionInVolume() + tau1 * sx * u + tau1 * sy * v;
                Wm3::Vector3d p2 = mp->getPositionInVolume() + tau2 * sx * u + tau2 * sy * v;

                // Replace with the point closer to the original one

                double dist1 = (p1 - m_surfacePoints[i][j]).SquaredLength();
                double dist2 = (p2 - m_surfacePoints[i][j]).SquaredLength();

//                Wm3::Vector3d oldPoint = m_surfacePoints[i][j];
                m_surfacePoints[i][j] = (dist1 < dist2) ? p1 : p2;

                // Debug code: Matlab visualization of 3D points
//                std::cout << "pt0=[" << oldPoint << "];" << std::endl;
//                std::cout << "pts=[";
//                std::cout << m_surfacePoints[i][j] << std::endl;
//                for (size_t k = 0; k < 2 * numNormalsPerPoint; ++k)
//                {
//                    if (j == k) continue;
//                    std::cout << m_surfacePoints[i][k] << std::endl;
//                }
//                std::cout << "];" << std::endl;
//                std::cout << "pts2=[";
//                for (size_t k = 0; k < 2 * numNormalsPerPoint; ++k)
//                {
//                    std::cout << mp->getImprovedPositionInVolume() + ((k < numNormalsPerPoint) ? -1.0 : 1.0) * mp->getSamplingDirection(k % numNormalsPerPoint) * mp->getDiameterVoxels() / 2.0 << std::endl;
//                }
//                std::cout << "];" << std::endl;
//                std::cout << "pts3=[";
//                for (double k = -2.0; k <= 2.0; k += 0.2)
//                {
//                    std::cout << mp->getPositionInVolume() + k * mp->getSamplingDirection(j % numNormalsPerPoint) << std::endl;
//                }
//                std::cout << "];" << std::endl;
//                std::cout << "pts2 = pts2 - ones(20,1)*mean(pts);" << std::endl;
//                std::cout << "pts3 = pts3 - ones(21,1)*mean(pts);" << std::endl;
//                std::cout << "pt0 = pt0 - mean(pts);" << std::endl;
//                std::cout << "pts = pts - ones(20,1)*mean(pts);" << std::endl;
//                std::cout << "figure, plot3(pts(2:end,1),pts(2:end,2),pts(2:end,3), 'bx'), hold on, plot3(pts(1,1), pts(1,2), pts(1,3), 'rx'), plot3(pts2(:,1),pts2(:,2),pts2(:,3),'gx'), plot3(pt0(1,1),pt0(1,2),pt0(1,3),'o'), plot3(pts3(:,1),pts3(:,2),pts3(:,3),'go');" << std::endl;
//                std::cout << std::endl << "---" << std::endl;
            }
        }
    }
}

void Triangulation::extractTubePoints()
{
    const size_t numPoints = m_measurementPoints.size();
    const size_t numNormalsPerPoint = Settings::getSettings().getNumNormalsAtPoint();

    // Save points on vessel circumference
    m_surfacePoints.clear();
    m_surfacePoints.resize(numPoints);

    for (size_t i = 0; i < m_measurementPoints.size(); ++i)
    {
        MeasurementPoint *mp = m_measurementPoints[i];
        Wm3::Vector3d center = mp->getImprovedPositionInVolume();
        double radius = mp->getDiameter() / 2.0;

        // Extract points
        for (size_t j = 0; j < numNormalsPerPoint; ++j)
        {
            Wm3::Vector3d direction = mp->getSamplingDirection(j);
            m_surfacePoints[i].push_back(center + direction * radius);
        }

        for (size_t j = 0; j < numNormalsPerPoint; ++j)
        {
            Wm3::Vector3d direction = mp->getSamplingDirection(j);
            m_surfacePoints[i].push_back(center - direction * radius);
        }
    }
}

