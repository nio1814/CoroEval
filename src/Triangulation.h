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

#include "wm3/Wm3Vector3.h"

#include <vector>
#include <string>


class MeasurementPoint;
namespace OpenMesh
{
    struct DefaultTraits;
    template <class Traits>
    class TriMesh_ArrayKernelT;
}
typedef OpenMesh::TriMesh_ArrayKernelT<OpenMesh::DefaultTraits> Mesh;


class Triangulation
{
public:
    Triangulation(std::vector<MeasurementPoint*> &measurementPoints);
    ~Triangulation();

    bool saveMesh(const std::string &fileName);
    bool saveTubeMesh(const std::string &fileName);
    bool saveSurfacePoints(const std::string &fileName);

private:
    void extractSurfacePoints();
    void extractTubePoints();
    void triangulate();

    std::vector<MeasurementPoint*> &m_measurementPoints;
    std::vector< std::vector< Wm3::Vector3d > > m_surfacePoints;
    Mesh *m_mesh;

};

