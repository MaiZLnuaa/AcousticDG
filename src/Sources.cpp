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



    void getElSourceVector(Config &config, Mesh &mesh, int eq, int el, double* elSourceVector, double t, const std::vector<std::vector<int>> &srcIndices)
    {
        int elNumNodes = mesh.getElNumNodes();
        std::fill(elSourceVector, elSourceVector + elNumNodes, 0.0); // 清零

        if (eq != 0 && eq != 1 && eq != 2) return;

        const auto &elNodeTags = mesh.getElNodeTags();  // 当前单元的全局节点编号列表

        for (int src = 0; src < (int)config.sources.size(); ++src)
        {
            auto &s = config.sources[src];
            if (s.formula == "" && s.data.empty())
            {
                int src_type = static_cast<int>(s.source[0]);
                if ((src_type == 0 && eq != 0) ||
                    (src_type == 1 && eq != 1) ||
                    (src_type == 2 && (eq != 1 && eq != 2)))
                {
                    continue; // 当前源项不作用于这个 eq
                }

                double x = s.source[1], y = s.source[2], z = s.source[3];
                double amp = s.source[5];
                double freq = s.source[6];
                double phase = s.source[7];
                double duration = s.source[8];
                if (t >= duration) continue;

                for (int n = 0; n < elNumNodes; ++n)
                {
                    int globalIndex = el * elNumNodes + n; // 当前单元的节点在全局节点列表中的索引
                    int globalNode = elNodeTags[globalIndex]; // 当前单元的节点在全局节点列表中的编号

                    // 判断该节点是否在源项索引列表内
                    if (std::find(srcIndices[src].begin(), srcIndices[src].end(), globalIndex) == srcIndices[src].end())
                        continue;

                    std::vector<double> coord, paramCoord;
                    int _dim, _tag;
                    gmsh::model::mesh::getNode(globalNode, coord, paramCoord, _dim, _tag);

                    double dx = coord[0] - x;
                    double dy = coord[1] - y;
                    double dz = coord[2] - z;
                    double r2 = dx * dx + dy * dy + dz * dz;

                    double ALPHA = std::log(2) / 2.0;
                    double spatial = std::exp(-ALPHA * r2);

                    if (src_type == 0 && eq == 0)
                    {
                        elSourceVector[n] += amp * sin(2 * M_PI * freq * t + phase) * spatial;
                    }
                    else if (src_type == 1 && eq == 1)
                    {
                        elSourceVector[n] += amp * cos(M_PI / 10 * dx) * exp((-log(2) / 5)*dy * dy) * sin(2 * M_PI * freq * t);
                    }
                    else if (src_type == 2)
                    {
                        if (eq == 1)
                        {
                            elSourceVector[n] += amp * sin(M_PI / 20 * dx) * exp((-log(2) / 5)*dy * dy) * sin(2 * M_PI * freq * t);
                        }
                        else if (eq == 2)
                        {
                            elSourceVector[n] += -amp * sin(M_PI / 10 * dy) * exp((-log(2) / 5)*dx * dx) * sin(2 * M_PI * freq * t);
                        }
                    }

                    int elNumIntPts = mesh.getElNumIntPts();
                    double Sum = 0.0;
                    for (int g = 0; g < elNumIntPts; g++)
                    {
                        Sum += mesh.elBasisFct(g, n) * mesh.getElWeight(g) * mesh.elJacobianDet(el, g); 
                    }
                    elSourceVector[n] *= Sum;
                }
            }
        }
    }
    
} // namespace sources



