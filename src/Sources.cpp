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
    void getElSourceVector(Config &config, Mesh &mesh, const size_t eq, const size_t el, double* elSourceVector, double t)
    {
        int jId;
        int m_elNumNodes = mesh.getElNumNodes();
        int m_elNumIntPts = mesh.getElNumIntPts();
        double x = config.sources[0].source[1];
        double y = config.sources[0].source[2];
        double z = config.sources[0].source[3];
        double size = config.sources[0].source[4];
        double amp = config.sources[0].source[5];
        double freq = config.sources[0].source[6];
        double phase = config.sources[0].source[7];
        double duration = config.sources[0].source[8];
        if (eq == 0)
        {
            for (int i = 0; i < m_elNumNodes; i++)
            {
                elSourceVector[i] = 0.0;
                std::vector<double> coord, paramCoord;
                int _dim, _tag;
                gmsh::model::mesh::getNode(mesh.getElNodeTags()[el * m_elNumNodes + i], coord, paramCoord, _dim, _tag);
                for (int j = 0; j < m_elNumNodes; j++)
                {
                    jId = el * m_elNumNodes + j;
                    for (int g = 0; g < m_elNumIntPts; g++)
                    {
                        elSourceVector[i] += mesh.elBasisFct(g, i) * mesh.elBasisFct(g, j) * mesh.getElWeight(g) * mesh.elJacobianDet(el, g) * amp * exp(-((coord[0] - x) * (coord[0] - x) + (coord[1] - y) * (coord[1] - y) + (coord[2] - z) * (coord[2] - z)) / size) * cos(2 * M_PI * freq * t + phase) ;
                    }
                }
            }
        }
        else
        {
            for (int i = 0; i < m_elNumNodes; i++)
            {
                elSourceVector[i] = 0.0;
            }
        }
    }
    
} // namespace sources



