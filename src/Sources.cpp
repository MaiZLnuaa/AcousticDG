/*
 *  @File          Source.cpp
 *
 *  @Author        MaiZLnuaa <mai-zl@nuaa.edu.cn>
 *  @Date          Tue May 27 2025 16:26:11
 *
 *  @Description    
 */


#include "Sources.h"

namespace sources
{
std::vector<int> findSourceNodesForOne(Mesh& mesh, Sources& source)
{
    std::vector<int> indices;
    for (int n = 0; n < mesh.getNumNodes(); ++n)
    {
        std::vector<double> coord, paramCoord;
        int dim, tag;
        gmsh::model::mesh::getNode(mesh.getElNodeTags()[n], coord, paramCoord, dim, tag);

        double dx = coord[0] - source.source[1];
        double dy = coord[1] - source.source[2];
        double dz = coord[2] - source.source[3];
        double r2 = dx * dx + dy * dy + dz * dz;
        if (r2 < std::pow(source.source[4], 2))
        {
            indices.push_back(n);
        }
    }
    return indices;
}

void applySourceForOne(std::vector<std::vector<double>>& u,
                       const std::vector<int>& nodeIndices,
                       Sources& source,
                       double t)
{
    if (source.formula == "" && source.data.empty())
    {
        double amp = source.source[5];
        double freq = source.source[6];
        double phase = source.source[7];
        double duration = source.source[8];
        if (t < duration)
        {
            for (int n : nodeIndices)
                u[0][n] = amp * sin(2 * M_PI * freq * t + phase);
        }
    }
    else if (source.data.empty())
    {
        double duration = source.source[5];
        if (t < duration)
        {
            for (int n : nodeIndices)
                u[0][n] = source.value(t);
        }
    }
    else
    {
        for (int n : nodeIndices)
            u[0][n] = source.interpolate_value(t);
    }
}

    
} // namespace sources



