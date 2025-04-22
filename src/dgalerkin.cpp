#include <cstdio>
#include <errno.h>
#include <gmsh.h>
#include <iostream>
#include <omp.h>

#include <parallel/algorithm>
#include <parallel/settings.h>

#include "Mesh.h"
#include "configParser.h"
#include "solver.h"

int main(int argc, char **argv)
{
    /**
     * The DGarlekin solver requires 2 arguments
     * 1 : the Mesh file (.msh)
     * 2 : the config file (.conf)
     *
     * e.g. ./dgarlerkin mymesh.msh myconfig.conf
     */

    // __gnu_parallel::_Settings s;
    // s.algorithm_strategy = __gnu_parallel::force_parallel;
    // __gnu_parallel::_Settings::set(s);

    auto start = std::chrono::system_clock::now();

    if (argc != 2)
    {
        return E2BIG;
    }
    std::string config_name = argv[1];

    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 1.0);

    Config config;

    if (fileExtension(config_name) == "conf")
        config = config::parseConfig(config_name);
    if (fileExtension(config_name) == "json")
        config = config::parseJSON(config_name);    

    gmsh::logger::write("Config loaded : " + config_name);

    Mesh mesh(config);

    int restartSimulation = config.restartSimulation;

    if (restartSimulation == 1)
    {
        // std::string resartFileName = "highorder_results/result0.vtu";
        std::string resartFileName = config.resartFileName;


        std::vector<double> pressure_value;
        // std::vector<std::vector<double>> velocity_value;
        std::vector<double> velocity_value;
        RestartSimulation::restartsimulation(resartFileName, pressure_value, velocity_value);

        std::vector<std::vector<double>> u(4, std::vector<double>(mesh.getNumNodes(), 0));

        
        for (int n = 0; n < mesh.getNumNodes(); n++)
        {
            u[0][n] = pressure_value[mesh.getElNodeTags()[n]-1];
            // u[1][n] = velocity_value[mesh.getElNodeTags()[n]][0];
            // u[2][n] = velocity_value[mesh.getElNodeTags()[n]][1];
            // u[3][n] = velocity_value[mesh.getElNodeTags()[n]][2];
            u[1][n] = velocity_value[(mesh.getElNodeTags()[n]-1)*3];
            u[2][n] = velocity_value[(mesh.getElNodeTags()[n]-1)*3 + 1];
            u[3][n] = velocity_value[(mesh.getElNodeTags()[n]-1)*3 + 2];
        }
            

        /**
        * Start solver
        */
        if (config.timeIntMethod == "Euler1")
            solver::forwardEuler(u, mesh, config);
        else if (config.timeIntMethod == "Runge-Kutta")
            solver::rungeKutta(u, mesh, config);
        else Fatal_Error("Time integration method error")    
    }else
    {
        /**
         * Initialize the solution:
         */
        std::vector<std::vector<double>> u(4, std::vector<double>(mesh.getNumNodes(), 0));
        for (int i = 0; i < config.initConditions.size(); ++i)
        {
            double x = config.initConditions[i][1];
            double y = config.initConditions[i][2];
            double z = config.initConditions[i][3];
            double size = config.initConditions[i][4];
            double amp = config.initConditions[i][5];

            for (int n = 0; n < mesh.getNumNodes(); n++)
            {
                std::vector<double> coord, paramCoord;
                int _dim, _tag;
                gmsh::model::mesh::getNode(mesh.getElNodeTags()[n], coord, paramCoord, _dim, _tag);
                u[0][n] += amp * exp(-((coord[0] - x) * (coord[0] - x) +
                                    (coord[1] - y) * (coord[1] - y) +
                                    (coord[2] - z) * (coord[2] - z)) /
                                    size);
            }
        }

        /**
        * Start solver
        */
        if (config.timeIntMethod == "Euler1")
            solver::forwardEuler(u, mesh, config);
        else if (config.timeIntMethod == "Runge-Kutta")
            solver::rungeKutta(u, mesh, config);
        else Fatal_Error("Time integration method error")    
    }
    


    

    // /**
    //  * Initialize the solution:
    //  */
    // std::vector<std::vector<double>> u(4, std::vector<double>(mesh.getNumNodes(), 0));
    // for (int i = 0; i < config.initConditions.size(); ++i)
    // {
    //     double x = config.initConditions[i][1];
    //     double y = config.initConditions[i][2];
    //     double z = config.initConditions[i][3];
    //     double size = config.initConditions[i][4];
    //     double amp = config.initConditions[i][5];

    //     for (int n = 0; n < mesh.getNumNodes(); n++)
    //     {
    //         std::vector<double> coord, paramCoord;
    //         int _dim, _tag;
    //         gmsh::model::mesh::getNode(mesh.getElNodeTags()[n], coord, paramCoord, _dim, _tag);
    //         u[0][n] += amp * exp(-((coord[0] - x) * (coord[0] - x) +
    //                                (coord[1] - y) * (coord[1] - y) +
    //                                (coord[2] - z) * (coord[2] - z)) /
    //                              size);
    //     }
    // }

    // /**
    //  * Start solver
    //  */
    // if (config.timeIntMethod == "Euler1")
    //     solver::forwardEuler(u, mesh, config);
    // else if (config.timeIntMethod == "Runge-Kutta")
    //     solver::rungeKutta(u, mesh, config);
    // else Fatal_Error("Time integration method error")    

    // mesh.writePVD("results.pvd");
    mesh.writePVD_highOrder("results_highorder.pvd");

    int i_m_elDim;
    std::vector<int> i_m_elType;

    i_m_elDim = gmsh::model::getDimension();
    gmsh::model::mesh::getElementTypes(i_m_elType, i_m_elDim);

    int _numPrimaryNodes = 0;
    std::string m_elName;
    int m_elDim, m_elOrder, m_elNumNodes;
    std::vector<double> m_elParamCoord;

    gmsh::model::mesh::getElementProperties(i_m_elType[0], m_elName, m_elDim,
                                            m_elOrder, m_elNumNodes, m_elParamCoord, _numPrimaryNodes);

    gmsh::logger::write("==================================================");
   
    gmsh::logger::write("Element dimension : " + std::to_string(m_elDim));
    gmsh::logger::write("Element Type : " + m_elName);
    gmsh::logger::write("Element Order : " + std::to_string(m_elOrder));
                                     
    screen_display::write_string("Calculation finished", GREEN);
    

    gmsh::finalize();

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    screen_display::write_value("Total time:", elapsed.count() * 1.0e-6, "s", BLUE);

    return EXIT_SUCCESS;
}
