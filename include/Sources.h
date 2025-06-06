/*
 *  @File          Sources.h
 *
 *  @Author        MaiZLnuaa <mai-zl@nuaa.edu.cn>
 *  @Date          Tue May 27 2025 16:23:53
 *
 *  @Description    
 */

#pragma once

#include "Mesh.h"
#include <gmsh.h>

namespace sources
{
    std::vector<int> findSourceNodesForOne(Mesh& mesh, Sources& source);

    void applySourceForOne(std::vector<std::vector<double>>& u,
                        const std::vector<int>& nodeIndices,
                        Sources& source,
                        double t);

    void getElSourceVector(Config &config, Mesh &mesh, int eq, int el, double* elSourceVector, double t, const std::vector<std::vector<int>> &srcIndices);
    
} // namespace sources
