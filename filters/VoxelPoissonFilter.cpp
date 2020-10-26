/******************************************************************************
 * Copyright (c) 2020, Bradley J Chambers (brad.chambers@gmail.com)
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following
 * conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in
 *       the documentation and/or other materials provided
 *       with the distribution.
 *     * Neither the name of Hobu, Inc. or Flaxen Geo Consulting nor the
 *       names of its contributors may be used to endorse or promote
 *       products derived from this software without specific prior
 *       written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGE.
 ****************************************************************************/

#include "VoxelPoissonFilter.hpp"

namespace pdal
{

static StaticPluginInfo const s_info{
    "filters.voxelpoisson", "Voxelized Poisson Sampling",
    "http://pdal.io/stages/filters.voxelpoisson.html"};

CREATE_STATIC_STAGE(VoxelPoissonFilter, s_info)

VoxelPoissonFilter::VoxelPoissonFilter() {}

std::string VoxelPoissonFilter::getName() const
{
    return s_info.name;
}

void VoxelPoissonFilter::addArgs(ProgramArgs& args)
{
    args.add("cell", "Cell size", m_cell, 1.0);
    // consider specifying radius instead, relationship being that the cell
    // size should be radius/sqrt(2), or radius should be sqrt(2)*cell
}

void VoxelPoissonFilter::ready(PointTableRef)
{
    m_populatedVoxels.clear();
    m_radius = std::sqrt(2.0) * m_cell * 0.5;
    log()->get(LogLevel::Debug)
        << "cell " << m_cell << ", radius " << m_radius << std::endl;
}

PointViewSet VoxelPoissonFilter::run(PointViewPtr view)
{
    PointViewPtr output = view->makeNew();
    PointRef point(*view);
    for (PointId id = 0; id < view->size(); ++id)
    {
        point.setPointId(id);
        if (voxelize(point))
            output->appendPoint(*view, id);
    }

    PointViewSet viewSet;
    viewSet.insert(output);
    return viewSet;
}

bool VoxelPoissonFilter::voxelize(PointRef& point)
{
    /*
     * Calculate the voxel coordinates for the incoming point.
     * gx, gy, gz will be the global coordinates from (0, 0, 0).
     */
    double x = point.getFieldAs<double>(Dimension::Id::X);
    double y = point.getFieldAs<double>(Dimension::Id::Y);
    double z = point.getFieldAs<double>(Dimension::Id::Z);
    if (m_populatedVoxels.empty())
    {
        m_originX = x - (m_cell / 2);
        m_originY = y - (m_cell / 2);
        m_originZ = z - (m_cell / 2);
    }

    // Offset by origin.
    x -= m_originX;
    y -= m_originY;
    z -= m_originZ;

    Voxel v = std::make_tuple((int)(std::floor(x / m_cell)),
                              (int)(std::floor(y / m_cell)),
                              (int)(std::floor(z / m_cell)));

    // this is silly, works for now
    double xd = point.getFieldAs<double>(Dimension::Id::X);
    double yd = point.getFieldAs<double>(Dimension::Id::Y);
    double zd = point.getFieldAs<double>(Dimension::Id::Z);

    // get all points in m_populateVoxels +/- 1 voxel away from v in each
    // dimension compute minimum distance to any of these points (there is
    // probably an optimization we can make here), and only insert point if
    // greater than radius away

    double mindist = std::numeric_limits<double>::max();
    for (int xi = std::get<0>(v) - 1; xi < std::get<0>(v) + 2; ++xi)
    {
        for (int yi = std::get<1>(v) - 1; yi < std::get<1>(v) + 2; ++yi)
        {
            for (int zi = std::get<2>(v) - 1; zi < std::get<2>(v) + 2; ++zi)
            {
                auto it = m_populatedVoxels.find(std::make_tuple(xi, yi, zi));
                if (it == m_populatedVoxels.end())
                    continue;
                std::vector<std::tuple<double, double, double>> coords =
                    m_populatedVoxels[std::make_tuple(xi, yi, zi)];
                for (std::tuple<double, double, double> const& coord : coords)
                {
                    double xv = std::get<0>(coord);
                    double yv = std::get<1>(coord);
                    double zv = std::get<2>(coord);
                    double dist = std::sqrt((xv - xd) * (xv - xd) +
                                            (yv - yd) * (yv - yd) +
                                            (zv - zd) * (zv - zd));
                    if (dist < mindist)
                        mindist = dist;
                }
            }
        }
    }

    if (mindist > m_radius)
    {
        std::tuple<double, double, double> vd = std::make_tuple(xd, yd, zd);
        auto it = m_populatedVoxels.find(v);
        if (it != m_populatedVoxels.end())
        {
            m_populatedVoxels[v].push_back(vd);
        }
        else
        {
            std::vector<std::tuple<double, double, double>> coord;
            coord.push_back(vd);
            m_populatedVoxels.emplace(v, coord);
        }
        return true;
    }

    return false;
}

bool VoxelPoissonFilter::processOne(PointRef& point)
{
    return voxelize(point);
}

} // namespace pdal
