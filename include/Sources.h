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
    void getElSourceVector(Config &config, Mesh &mesh, const size_t eq, const size_t el, double* elSourceVector, double t);
    
} // namespace sources
