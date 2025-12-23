#include <chrono>
#include <gmsh.h>
#include <iostream>
#include <omp.h>
#include <sstream>
#include <utils.h>
#include <vector>

#include "Mesh.h"
#include "configParser.h"

#include "Sources.h"

namespace solver
{

    /**
     * Common variables to all solver
     */
    int elNumNodes;
    int numNodes;
    std::vector<std::string> g_names;
    std::vector<int> elTags;
    std::vector<double> elFlux;
    std::vector<double> elStiffvector;
    std::vector<double> elZKHvector;
    std::vector<double> eldampingPvector;
    std::vector<double> elAuxiliaryTerm1Vector;
    std::vector<double> elAuxiliaryTerm2Vector;
    std::vector<double> eljcaPhiVector;
    std::vector<double> eljcaAuxiliaryTerm1Vector;
    std::vector<double> eljcaAuxiliaryTerm2Vector;
    std::vector<std::vector<std::vector<double>>> Flux;

    std::vector<std::vector<float>> data4wave;

    std::vector<double> elSourceVector;

    /**
     * Perform a numerical step: u[t+1] = dt*M^-1*(S[u[t]]-F[u[t]]) + beta*u[t]
     * for all elements in mesh object.
     *
     * @param mesh Mesh object
     * @param config Configuration file
     * @param u Nodal solution vector
     * @param Flux Nodal physical Flux
     * @param beta double coefficient
     */
    void numStep(Mesh &mesh, Config config, std::vector<std::vector<double>> &u,
                 std::vector<std::vector<std::vector<double>>> &Flux, double beta)
    {

        for (int eq = 0; eq < 4; ++eq)
        {
            // mesh.precomputeFlux(u[eq], Flux[eq], eq);
            mesh.newprecomputeFlux(u[eq], Flux[eq], eq, u[1], u[2], u[3]);

#pragma omp parallel for schedule(static) firstprivate(elFlux, elStiffvector) num_threads(config.numThreads)
            for (int el = 0; el < mesh.getElNum(); ++el)
            {

                mesh.getElFlux(el, elFlux.data());
                mesh.getElStiffVector(el, Flux[eq], u[eq], elStiffvector.data()); // 获得 S_k
                eigen::minus(elStiffvector.data(), elFlux.data(), elNumNodes); // S_k - F_k
                eigen::linEq(&mesh.elMassMatrix(el), &elStiffvector[0], &u[eq][el * elNumNodes],
                             config.timeStep, beta, elNumNodes); // 求 u[t+1] = beta * u[t] + dt * M^-1 * (S_k - F_k)
            }
        }
    }

    void pmlnumStep(Mesh &mesh, Config config, std::vector<std::vector<double>> &u, 
                    std::vector<std::vector<double>> &pml_phi, std::vector<std::vector<std::vector<double>>> &Flux, double beta)
    {

        std::vector<std::vector<double>> u_old(4, std::vector<double>(mesh.getNumNodes(), 0.0));
        for (int eq = 0; eq < 4; ++eq)
        {
            // mesh.precomputeFlux(u[eq], Flux[eq], eq);
            mesh.newprecomputeFlux(u[eq], Flux[eq], eq, u[1], u[2], u[3]);
            u_old[eq] = u[eq];

#pragma omp parallel for schedule(static) firstprivate(elFlux, elStiffvector, eldampingPvector) num_threads(config.numThreads)
            for (int el = 0; el < mesh.getElNum(); ++el)
            {

                mesh.getElFlux(el, elFlux.data());
                mesh.getElStiffVector(el, Flux[eq], u[eq], elStiffvector.data()); // 获得 S_k
                mesh.getDampingPressureVector(eq, el, u[eq], pml_phi, eldampingPvector.data()); // 获得 dampingP_k
                eigen::minus(elStiffvector.data(), elFlux.data(), elNumNodes); // S_k - F_k
                eigen::minus(elStiffvector.data(), eldampingPvector.data(), elNumNodes); // S_k - F_k - dampingP_k
                eigen::linEq(&mesh.elMassMatrix(el), &elStiffvector[0], &u[eq][el * elNumNodes],
                             config.timeStep, beta, elNumNodes); // 求 u[t+1] = beta * u[t] + dt * M^-1 * (S_k - F_k-dampingP_k)

            }
        }

        for (int eq = 0; eq < 3; eq++)
        {
            int velEq = eq + 1;
            #pragma omp parallel for schedule(static) firstprivate(elAuxiliaryTerm1Vector, elAuxiliaryTerm2Vector) num_threads(config.numThreads)
            for (int el = 0; el < mesh.getElNum(); el++)
            {
                if (!mesh.isPML(mesh.elTag(el)))
                {
                    continue;
                }
                mesh.getAuxiliaryEquationTerm1(eq, el, pml_phi, elAuxiliaryTerm1Vector.data());
                mesh.getAuxiliaryEquationTerm2(eq, el, u[velEq], u_old[velEq], pml_phi, elAuxiliaryTerm2Vector.data());
                eigen::minusFromZero(elAuxiliaryTerm1Vector.data(), elAuxiliaryTerm1Vector.data(), elNumNodes); // - Term1_k
                eigen::plus(elAuxiliaryTerm1Vector.data(), elAuxiliaryTerm2Vector.data(), elNumNodes); // + Term2_k
                eigen::linEq(&mesh.elMassMatrix(el), &elAuxiliaryTerm1Vector[0], &pml_phi[eq][el * elNumNodes],
                             config.timeStep, beta, elNumNodes); // 求 pml_phi[t+1] = beta * pml_phi[t] + dt * M^-1 * Term1_k
            }
            
        }
        
    }
        
    void zknumStep(Mesh &mesh, Config config, std::vector<std::vector<double>> &u,
                 std::vector<std::vector<std::vector<double>>> &Flux, double beta)
    {

        // porous media parameters
        double porosity = config.porousParams[0][0];
        double tortuosity = config.porousParams[0][1];
        double resistivity = config.porousParams[0][2];
        double gamma = config.porousParams[0][3];
        for (int eq = 0; eq < 4; ++eq)
        {
            // mesh.precomputeFlux(u[eq], Flux[eq], eq);
            mesh.newprecomputeFlux(u[eq], Flux[eq], eq, u[1], u[2], u[3]);
#pragma omp parallel for schedule(static) firstprivate(elFlux, elStiffvector, elZKHvector) num_threads(config.numThreads)
            for (int el = 0; el < mesh.getElNum(); ++el)
            {
                mesh.getElFlux(el, elFlux.data());
                mesh.getElStiffVector(el, Flux[eq], u[eq], elStiffvector.data()); // 获得 S_k
                eigen::minus(elStiffvector.data(), elFlux.data(), elNumNodes); // S_k - F_k

                mesh.getZKPorousHVector(el, eq, u[eq], elZKHvector.data(), porosity, tortuosity, resistivity, gamma); // 获得 H_k
                eigen::minus(elStiffvector.data(), elZKHvector.data(), elNumNodes); // S_k - F_k - H_k

                eigen::linEq(&mesh.elMassMatrix(el), &elStiffvector[0], &u[eq][el * elNumNodes],
                             config.timeStep, beta, elNumNodes); // 求 u[t+1] = beta * u[t] + dt * M^-1 * (S_k - F_k - H_k)
            }
        }
    }

    void jcanumStep(Mesh &mesh, Config config, std::vector<std::vector<double>> &u, 
                    std::vector<std::vector<double>> &jca_phi, std::vector<std::vector<std::vector<double>>> &Flux, double beta)
    {

        // porous media parameters
        double porosity = config.jcaPorousParams[0][0];
        double tortuosity = config.jcaPorousParams[0][1];
        double resistivity = config.jcaPorousParams[0][2];
        double gamma = config.jcaPorousParams[0][3];
        double viscousLength = config.jcaPorousParams[0][4];
        double thermalLength = config.jcaPorousParams[0][5];
        double prandtlNumber = config.jcaPorousParams[0][6];

        std::vector<std::vector<double>> u_old(4, std::vector<double>(mesh.getNumNodes(), 0.0));
        
        for (int eq = 0; eq < 4; ++eq)
        {
            // mesh.precomputeFlux(u[eq], Flux[eq], eq);
            mesh.newprecomputeFlux(u[eq], Flux[eq], eq, u[1], u[2], u[3]);
            u_old[eq] = u[eq];

#pragma omp parallel for schedule(static) firstprivate(elFlux, elStiffvector, eljcaPhiVector) num_threads(config.numThreads)
            for (int el = 0; el < mesh.getElNum(); ++el)
            {
                mesh.getElFlux(el, elFlux.data());
                mesh.getElStiffVector(el, Flux[eq], u[eq], elStiffvector.data()); // 获得 S_k
                eigen::minus(elStiffvector.data(), elFlux.data(), elNumNodes); // S_k - F_k

                mesh.getZKPorousHVector(el, eq, u[eq], elZKHvector.data(), porosity, tortuosity, resistivity, gamma); // 获得 H_k
                eigen::minus(elStiffvector.data(), elZKHvector.data(), elNumNodes); // S_k - F_k - H_k

                mesh.getjcaPorousPhiVector(el, eq, jca_phi, eljcaPhiVector.data()); // 获得 jca_phi_k
                eigen::minus(elStiffvector.data(), eljcaPhiVector.data(), elNumNodes); // S_k - F_k - H_k - jca_phi_k

                eigen::linEq(&mesh.elMassMatrix(el), &elStiffvector[0], &u[eq][el * elNumNodes],
                             config.timeStep, beta, elNumNodes); // 求 u[t+1] = beta * u[t] + dt * M^-1 * (S_k - F_k - H_k - jca_phi_k)
            }

        }

        for (int eq = 0; eq < 4; eq++)
        {
#pragma omp parallel for schedule(static) firstprivate(eljcaAuxiliaryTerm1Vector, eljcaAuxiliaryTerm2Vector) num_threads(config.numThreads)
            for (int el = 0; el < mesh.getElNum(); el++)
            {
                mesh.getjcaAuxiliaryTerm1(eq, el, jca_phi, eljcaAuxiliaryTerm1Vector.data());
                mesh.getjcaAuxiliaryTerm2(eq, el, u[eq], u_old[eq], jca_phi, eljcaAuxiliaryTerm2Vector.data());
                eigen::minusFromZero(eljcaAuxiliaryTerm2Vector.data(), eljcaAuxiliaryTerm2Vector.data(), elNumNodes); // - Term2_k
                eigen::plus(eljcaAuxiliaryTerm1Vector.data(), eljcaAuxiliaryTerm2Vector.data(), elNumNodes); // + Term1_k
                eigen::linEq(&mesh.elMassMatrix(el), &eljcaAuxiliaryTerm1Vector[0], &jca_phi[eq][el * elNumNodes],
                             config.timeStep, beta, elNumNodes); // 求 jca_phi[t+1] = beta * jca_phi[t] + dt * M^-1 * (Term1_k - Term2_k)
            }
            
        }
        
    }


    void pmlzknumStep(Mesh &mesh, Config config, std::vector<std::vector<double>> &u, std::vector<std::vector<double>> &pml_phi, 
                    std::vector<std::vector<std::vector<double>>> &Flux, double beta)
    {
        std::vector<std::vector<double>> u_old(4, std::vector<double>(mesh.getNumNodes(), 0.0));

        // porous media parameters
        double porosity = config.porousParams[0][0];
        double tortuosity = config.porousParams[0][1];
        double resistivity = config.porousParams[0][2];
        double gamma = config.porousParams[0][3];
        for (int eq = 0; eq < 4; ++eq)
        {
            // mesh.precomputeFlux(u[eq], Flux[eq], eq);
            mesh.newprecomputeFlux(u[eq], Flux[eq], eq, u[1], u[2], u[3]);
            u_old[eq] = u[eq];
#pragma omp parallel for schedule(static) firstprivate(elFlux, elStiffvector, elZKHvector, eldampingPvector) num_threads(config.numThreads)
            for (int el = 0; el < mesh.getElNum(); ++el)
            {
                mesh.getElFlux(el, elFlux.data());
                mesh.getElStiffVector(el, Flux[eq], u[eq], elStiffvector.data()); // 获得 S_k
                eigen::minus(elStiffvector.data(), elFlux.data(), elNumNodes); // S_k - F_k

                mesh.getZKPorousHVector(el, eq, u[eq], elZKHvector.data(), porosity, tortuosity, resistivity, gamma); // 获得 H_k
                eigen::minus(elStiffvector.data(), elZKHvector.data(), elNumNodes); // S_k - F_k - H_k

                mesh.getDampingPressureVector(eq, el, u[eq], pml_phi, eldampingPvector.data());
                eigen::minus(elStiffvector.data(), eldampingPvector.data(), elNumNodes); // S_k - F_k - H_k - dampingP_k
                eigen::linEq(&mesh.elMassMatrix(el), &elStiffvector[0], &u[eq][el * elNumNodes],
                             config.timeStep, beta, elNumNodes); // 求 u[t+1] = beta * u[t] + dt * M^-1 * (S_k - F_k - H_k - dampingP_k)
            }
        }

        for (int eq = 0; eq < 3; eq++)
        {
            int velEq = eq + 1;
#pragma omp parallel for schedule(static) firstprivate(elAuxiliaryTerm1Vector, elAuxiliaryTerm2Vector) num_threads(config.numThreads)
            for (int el = 0; el < mesh.getElNum(); el++)
            {
                if (!mesh.isPML(mesh.elTag(el)))
                {
                    continue;
                }
                mesh.getAuxiliaryEquationTerm1(eq, el, pml_phi, elAuxiliaryTerm1Vector.data());
                mesh.getAuxiliaryEquationTerm2(eq, el, u[velEq], u_old[velEq], pml_phi, elAuxiliaryTerm2Vector.data());
                eigen::minusFromZero(elAuxiliaryTerm1Vector.data(), elAuxiliaryTerm1Vector.data(), elNumNodes); // - Term1_k
                eigen::plus(elAuxiliaryTerm1Vector.data(), elAuxiliaryTerm2Vector.data(), elNumNodes); // + Term2_k
                eigen::linEq(&mesh.elMassMatrix(el), &elAuxiliaryTerm1Vector[0], &pml_phi[eq][el * elNumNodes],
                                config.timeStep, beta, elNumNodes); // 求 pml_phi[t+1] = beta * pml_phi[t] + dt * M^-1 * Term1_k
            }
            
        }
    }

    void pmljcanumStep(Mesh &mesh, Config config, std::vector<std::vector<double>> &u, 
                    std::vector<std::vector<double>> &jca_phi, std::vector<std::vector<double>> &pml_phi, std::vector<std::vector<std::vector<double>>> &Flux, double beta)
    {

        // porous media parameters
        double porosity = config.jcaPorousParams[0][0];
        double tortuosity = config.jcaPorousParams[0][1];
        double resistivity = config.jcaPorousParams[0][2];
        double gamma = config.jcaPorousParams[0][3];
        double viscousLength = config.jcaPorousParams[0][4];
        double thermalLength = config.jcaPorousParams[0][5];
        double prandtlNumber = config.jcaPorousParams[0][6];

        std::vector<std::vector<double>> u_old(4, std::vector<double>(mesh.getNumNodes(), 0.0));
        
        for (int eq = 0; eq < 4; ++eq)
        {
            // mesh.precomputeFlux(u[eq], Flux[eq], eq);
            mesh.newprecomputeFlux(u[eq], Flux[eq], eq, u[1], u[2], u[3]);
            u_old[eq] = u[eq];

#pragma omp parallel for schedule(static) firstprivate(elFlux, elStiffvector, eljcaPhiVector) num_threads(config.numThreads)
            for (int el = 0; el < mesh.getElNum(); ++el)
            {
                mesh.getElFlux(el, elFlux.data());
                mesh.getElStiffVector(el, Flux[eq], u[eq], elStiffvector.data()); // 获得 S_k
                eigen::minus(elStiffvector.data(), elFlux.data(), elNumNodes); // S_k - F_k

                mesh.getZKPorousHVector(el, eq, u[eq], elZKHvector.data(), porosity, tortuosity, resistivity, gamma); // 获得 H_k
                eigen::minus(elStiffvector.data(), elZKHvector.data(), elNumNodes); // S_k - F_k - H_k

                mesh.getjcaPorousPhiVector(el, eq, jca_phi, eljcaPhiVector.data()); // 获得 jca_phi_k
                eigen::minus(elStiffvector.data(), eljcaPhiVector.data(), elNumNodes); // S_k - F_k - H_k - jca_phi_k

                mesh.getDampingPressureVector(eq, el, u[eq], pml_phi, eldampingPvector.data());
                eigen::minus(elStiffvector.data(), eldampingPvector.data(), elNumNodes); // S_k - F_k - H_k - jca_phi_k - dampingP_k

                eigen::linEq(&mesh.elMassMatrix(el), &elStiffvector[0], &u[eq][el * elNumNodes],
                             config.timeStep, beta, elNumNodes); // 求 u[t+1] = beta * u[t] + dt * M^-1 * (S_k - F_k - H_k - jca_phi_k)
            }

        }

        for (int eq = 0; eq < 4; eq++)
        {
#pragma omp parallel for schedule(static) firstprivate(eljcaAuxiliaryTerm1Vector, eljcaAuxiliaryTerm2Vector) num_threads(config.numThreads)
            for (int el = 0; el < mesh.getElNum(); el++)
            {
                mesh.getjcaAuxiliaryTerm1(eq, el, jca_phi, eljcaAuxiliaryTerm1Vector.data());
                mesh.getjcaAuxiliaryTerm2(eq, el, u[eq], u_old[eq], jca_phi, eljcaAuxiliaryTerm2Vector.data());
                eigen::minusFromZero(eljcaAuxiliaryTerm2Vector.data(), eljcaAuxiliaryTerm2Vector.data(), elNumNodes); // - Term2_k
                eigen::plus(eljcaAuxiliaryTerm1Vector.data(), eljcaAuxiliaryTerm2Vector.data(), elNumNodes); // + Term1_k
                eigen::linEq(&mesh.elMassMatrix(el), &eljcaAuxiliaryTerm1Vector[0], &jca_phi[eq][el * elNumNodes],
                             config.timeStep, beta, elNumNodes); // 求 jca_phi[t+1] = beta * jca_phi[t] + dt * M^-1 * (Term1_k - Term2_k)
            }
            
        }

        for (int eq = 0; eq < 3; eq++)
        {
            int velEq = eq + 1;
#pragma omp parallel for schedule(static) firstprivate(elAuxiliaryTerm1Vector, elAuxiliaryTerm2Vector) num_threads(config.numThreads)
            for (int el = 0; el < mesh.getElNum(); el++)
            {
                if (!mesh.isPML(mesh.elTag(el)))
                {
                    continue;
                }
                mesh.getAuxiliaryEquationTerm1(eq, el, pml_phi, elAuxiliaryTerm1Vector.data());
                mesh.getAuxiliaryEquationTerm2(eq, el, u[velEq], u_old[velEq], pml_phi, elAuxiliaryTerm2Vector.data());
                eigen::minusFromZero(elAuxiliaryTerm1Vector.data(), elAuxiliaryTerm1Vector.data(), elNumNodes); // - Term1_k
                eigen::plus(elAuxiliaryTerm1Vector.data(), elAuxiliaryTerm2Vector.data(), elNumNodes); // + Term2_k
                eigen::linEq(&mesh.elMassMatrix(el), &elAuxiliaryTerm1Vector[0], &pml_phi[eq][el * elNumNodes],
                                config.timeStep, beta, elNumNodes); // 求 pml_phi[t+1] = beta * pml_phi[t] + dt * M^-1 * Term1_k
            }
            
        }
        
    }

    void sourcenumStep(Mesh &mesh, Config config, std::vector<std::vector<double>> &u,
                 std::vector<std::vector<std::vector<double>>> &Flux, double beta, double t)
    {

        for (int eq = 0; eq < 4; ++eq)
        {
            // mesh.precomputeFlux(u[eq], Flux[eq], eq);
            mesh.newprecomputeFlux(u[eq], Flux[eq], eq, u[1], u[2], u[3]);

#pragma omp parallel for schedule(static) firstprivate(elFlux, elStiffvector, elSourceVector) num_threads(config.numThreads)
            for (int el = 0; el < mesh.getElNum(); ++el)
            {

                mesh.getElFlux(el, elFlux.data());
                mesh.getElStiffVector(el, Flux[eq], u[eq], elStiffvector.data()); // 获得 S_k
                sources::getElSourceVector(config, mesh, eq, el, elSourceVector.data(), t); // 获得 source_k
                eigen::minus(elStiffvector.data(), elFlux.data(), elNumNodes); // S_k - F_k
                eigen::plus(elStiffvector.data(),elSourceVector.data(),elNumNodes); // S_k - F_k + source_k
                eigen::linEq(&mesh.elMassMatrix(el), &elStiffvector[0], &u[eq][el * elNumNodes],
                             config.timeStep, beta, elNumNodes); // 求 u[t+1] = beta * u[t] + dt * M^-1 * (S_k - F_k + source_k)
            }
        }
    }

    /**
     * Solve using forward explicit scheme. O(h)
     *
     * @param u initial nodal solution vector
     * @param mesh
     * @param config
     */
    void forwardEuler(std::vector<std::vector<double>> &u, Mesh &mesh, Config config, std::vector<std::vector<double>> &pml_phi, std::vector<std::vector<double>> &jca_phi)
    {

        /** Memory allocation */
        elNumNodes = mesh.getElNumNodes();
        numNodes = mesh.getNumNodes();
        elTags = std::vector<int>(&mesh.elTag(0), &mesh.elTag(0) + mesh.getElNum());
        elFlux.resize(elNumNodes);
        elStiffvector.resize(elNumNodes);
        elZKHvector.resize(elNumNodes, 0.0);
        eldampingPvector.resize(elNumNodes, 0.0);
        elAuxiliaryTerm1Vector.resize(elNumNodes, 0.0);
        elAuxiliaryTerm2Vector.resize(elNumNodes, 0.0);
        eljcaAuxiliaryTerm1Vector.resize(elNumNodes, 0.0);
        eljcaAuxiliaryTerm2Vector.resize(elNumNodes, 0.0);
        eljcaPhiVector.resize(elNumNodes, 0.0);
        elSourceVector.resize(elNumNodes, 0.0);
        Flux = std::vector<std::vector<std::vector<double>>>(4, std::vector<std::vector<double>>(mesh.getNumNodes(), std::vector<double>(3)));

        /** Gmsh save init */
        gmsh::model::list(g_names);
        int gp_viewTag = gmsh::view::add("Pressure");
        int gv_viewTag = gmsh::view::add("Velocity");
        int grho_viewTag = gmsh::view::add("Density");
        std::vector<std::vector<double>> g_p(mesh.getElNum(), std::vector<double>(elNumNodes));
        std::vector<std::vector<double>> g_rho(mesh.getElNum(), std::vector<double>(elNumNodes));
        std::vector<std::vector<double>> g_v(mesh.getElNum(), std::vector<double>(3 * elNumNodes));

        /** Precomputation */
        mesh.precomputeMassMatrix();

        /** ComputeSigma */
        mesh.computeSigma();

        /** Source */
        std::vector<std::vector<int>> srcIndices;
        for (int i = 0; i < config.sources.size(); ++i)
        {
            std::vector<int> indice;
            for (int n = 0; n < mesh.getNumNodes(); n++)
            {
                std::vector<double> coord, paramCoord;
                int _dim, _tag;
                gmsh::model::mesh::getNode(mesh.getElNodeTags()[n], coord, paramCoord, _dim, _tag);
                if (pow(coord[0] - config.sources[i].source[1], 2) +
                        pow(coord[1] - config.sources[i].source[2], 2) +
                        pow(coord[2] - config.sources[i].source[3], 2) <
                    pow(config.sources[i].source[4], 2))
                {
                    indice.push_back(n);
                }
            }
            srcIndices.push_back(indice);
        }

        /** Observer */
        std::vector<std::vector<int>> obsIndices;
        std::vector<std::vector<double>> obsPtDistance;
        for (int i = 0; i < config.observers.size(); ++i)
        {
            std::vector<int> indice;
            std::vector<double> dist;
            for (int n = 0; n < mesh.getNumNodes(); n++)
            {
                std::vector<double> coord, paramCoord;
                int _dim, _tag;
                gmsh::model::mesh::getNode(mesh.getElNodeTags()[n], coord, paramCoord, _dim, _tag);
                double distance = sqrt(pow(coord[0] - config.observers[i][0], 2) +
                                       pow(coord[1] - config.observers[i][1], 2) +
                                       pow(coord[2] - config.observers[i][2], 2));
                if (distance < config.observers[i][3])
                {
                    indice.push_back(n);
                    dist.push_back(distance);
                }
            }
            obsIndices.push_back(indice);
            obsPtDistance.push_back(dist);
        }

        /**
         * Main Loop : Time iteration
         */
        std::ofstream outfile("residuals.csv");
        outfile << "time;res_p;res_rho;res_vx;res_vy;res_vz;elapsed_time" << std::endl;

        std::vector<std::ofstream> obs_outfile(config.observers.size());
        data4wave.clear();
        data4wave.resize(config.observers.size());
        for (int obs = 0; obs < config.observers.size(); ++obs)
        {
            std::string filename = "results/observers" + std::to_string(obs + 1) + ".txt";
            obs_outfile[obs].open(filename.c_str());
            obs_outfile[obs] << "time;pressure;density;velocity_x;velocity_y;velocity_z" << std::endl;
        }

        for (int i = 0; i < obsIndices.size(); i++)
        {
            std::cout << "Observer " << i
                    << " @ (" << config.observers[i][0] << ", "
                    << config.observers[i][1] << ", "
                    << config.observers[i][2] << ")"
                    << "  radius = " << config.observers[i][3]
                    << "  → contains " << obsIndices[i].size() << " nodes" << std::endl;

            for (int j = 0; j < obsIndices[i].size(); j++)
            {
                int nodeIndex = obsIndices[i][j];
                std::vector<double> coord, paramCoord;
                int dim, tag;
                gmsh::model::mesh::getNode(mesh.getElNodeTags()[nodeIndex], coord, paramCoord, dim, tag);

                std::cout << "    Node " << nodeIndex
                        << " : (" << coord[0] << ", "
                        << coord[1] << ", "
                        << coord[2] << ")"
                        << "  distance = " << obsPtDistance[i][j] << std::endl;
            }
            std::cout << std::endl;
        }
        std::cout << "Press Enter to start the simulation..." << std::endl;
        getchar(); // 暂停程序，按回车继续

        auto start = std::chrono::system_clock::now();
        for (double t = config.timeStart, step = int(config.timeStart / config.timeStep), tDisplay = 0; t <= config.timeEnd + config.timeStep / 2;
             t += config.timeStep, tDisplay += config.timeStep, ++step)
        {

            auto start_time = std::chrono::system_clock::now();
            std::vector<double> residual(5, 0.0);
            /**
             *  Savings and prints
             */

            if (tDisplay >= config.timeRate - config.timeStep / 2 || step == 0)
            {
                tDisplay = 0;

/** [1] Copy solution to match GMSH format */
#pragma omp parallel for schedule(static) num_threads(config.numThreads)
                for (int el = 0; el < mesh.getElNum(); ++el)
                {
                    for (int n = 0; n < mesh.getElNumNodes(); ++n)
                    {
                        int elN = el * elNumNodes + n;
                        g_p[el][n] = u[0][elN];
                        g_rho[el][n] = u[0][elN] / (config.c0 * config.c0);
                        g_v[el][3 * n + 0] = u[1][elN];
                        g_v[el][3 * n + 1] = u[2][elN];
                        g_v[el][3 * n + 2] = u[3][elN];
                    }
                }

                /** [2] Print and compute iteration time */
                auto end = std::chrono::system_clock::now();
                auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
                gmsh::logger::write("[" + std::to_string(t) + "/" + std::to_string(config.timeEnd) + "s] Step number : " + std::to_string((int)step) + ", Elapsed time: " + std::to_string(elapsed.count()) + "s");
                // screen_display::write_string("time\t\tres_p\t\tres_rho\t\tres_vx\t\tres_vy\t\tres_vz\t\telapsed time", BOLDBLUE);
                // mesh.writeVTK("result.vtk");
                // std::string vtu_filename = "results/result" + std::to_string((int)step) + ".vtu";
                // mesh.writeVTUb(vtu_filename, u);
                std::string vtu_filename_highOrder = "highorder_results/result" + std::to_string((int)step) + ".vtu";
                mesh.writeVTUb_highOrder(vtu_filename_highOrder, u);
                screen_display::write_string("time\t\tres_p\t\tres_rho\t\tres_vx\t\tres_vy\t\tres_vz\t\telapsed time", BOLDBLUE);
            }

            // /**
            //  * Update Source
            //  */

            // for (int src = 0; src < config.sources.size(); ++src)
            // {
            //     if (config.sources[src].formula == "" && config.sources[src].data.empty())
            //     {
            //         double amp = config.sources[src].source[5];
            //         double freq = config.sources[src].source[6];
            //         double phase = config.sources[src].source[7];
            //         double duration = config.sources[src].source[8];
            //         if (t < duration)
            //             for (int n = 0; n < srcIndices[src].size(); ++n)
            //                 u[0][srcIndices[src][n]] = amp * sin(2 * M_PI * freq * t + phase);
            //     }
            //     else
            //     {
            //         if (config.sources[src].data.empty())
            //         {
            //             double duration = config.sources[src].source[5];
            //             if (t < duration)
            //                 for (int n = 0; n < srcIndices[src].size(); ++n)
            //                     u[0][srcIndices[src][n]] = config.sources[src].value(t);
            //         }
            //         else
            //         {
            //             for (int n = 0; n < srcIndices[src].size(); ++n)
            //                 u[0][srcIndices[src][n]] = config.sources[src].interpolate_value(t);
            //         }
            //     }
            // }

            /**
             * First Order Euler
             */
            mesh.updateFlux(u, Flux, config.v0, config.c0, config.rho0);
            // mesh.updatezkFlux(u, Flux, config.v0, config.c0, config.rho0, config.porousParams[0][0], config.porousParams[0][1], config.porousParams[0][2], config.porousParams[0][3]);
            // numStep(mesh, config, u, Flux, 1);
            // pmlnumStep(mesh, config, u, pml_phi, Flux, 1);
            // pmlzknumStep(mesh, config, u, pml_phi, Flux, 1);
            // jcanumStep(mesh, config, u, jca_phi, Flux, 1);
            // pmljcanumStep(mesh, config, u, jca_phi, pml_phi, Flux, 1);
            sourcenumStep(mesh, config, u, Flux, 1, t);

            /**
             * Compute residuals
             */
#pragma omp parallel for schedule(static) num_threads(config.numThreads)
            for (int el = 0; el < mesh.getElNum(); ++el)
            {
                for (int n = 0; n < mesh.getElNumNodes(); ++n)
                {
                    int elN = el * elNumNodes + n;
#pragma omp atomic update
                    residual[0] += pow(g_p[el][n] - u[0][elN], 2);
#pragma omp atomic update
                    residual[1] += pow(g_rho[el][n] - u[0][elN] / (config.c0 * config.c0), 2);
#pragma omp atomic update
                    residual[2] += pow(g_v[el][3 * n + 0] - u[1][elN], 2);
#pragma omp atomic update
                    residual[3] += pow(g_v[el][3 * n + 1] - u[2][elN], 2);
#pragma omp atomic update
                    residual[4] += pow(g_v[el][3 * n + 2] - u[3][elN], 2);
                }
            }
            outfile << t + config.timeStep << ";";
            std::cout << std::scientific << t + config.timeStep << "\t";
            auto end_time = std::chrono::system_clock::now();
            auto elapsed_time = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
            for (int eq = 0; eq < residual.size(); ++eq)
            {
                residual[eq] /= (mesh.getElNum() * mesh.getElNumNodes());
                std::cout << std::scientific << residual[eq] << "\t";
                outfile << residual[eq] << ";";
            }
            std::cout << elapsed_time.count() * 1.0e-6 << " s" << std::endl;
            outfile << elapsed_time.count() * 1.0e-6 << std::endl;

            /**
             * get observers value
             * Franke-Little interpolation method
             */
            for (int obs = 0; obs < config.observers.size(); ++obs)
            {
                double p(0), rho(0), w_sum(0);
                std::vector<double> v = {0, 0, 0};
                for (int n = 0; n < obsIndices[obs].size(); ++n)
                {
                    double R = config.observers[obs][3];                        //! influence sphere
                    double w = 1.0 / (pow(obsPtDistance[obs][n], 2) + 1.0e-12); // fmax(1.0-obsPtDistance[obs][n]/R,0.0);//1.0 / pow(obsPtDistance[obs][n],2);
                    p += u[0][obsIndices[obs][n]] * w;
                    v[0] += u[1][obsIndices[obs][n]] * w;
                    v[1] += u[2][obsIndices[obs][n]] * w;
                    v[2] += u[3][obsIndices[obs][n]] * w;
                    w_sum += w;
                }
                p /= w_sum;
                rho = p / pow(config.c0, 2);
                v[0] /= w_sum;
                v[1] /= w_sum;
                v[2] /= w_sum;
                data4wave[obs].push_back(p);
                obs_outfile[obs] << t << ";" << p << ";" << rho << ";" << v[0] << ";" << v[1] << ";" << v[2] << std::endl;
            }
        }
        for (int obs = 0; obs < config.observers.size(); ++obs)
        {
            io::writeWave(data4wave[obs], "results/observer_" + std::to_string(obs + 1) + ".wav", 1.0 / config.timeStep, 16, 1, 1);
            io::writeFFT(data4wave[obs],config.timeStep,"results/observer_" + std::to_string(obs + 1));
        }

        outfile.close();
        for (int obs = 0; obs < config.observers.size(); ++obs)
            obs_outfile[obs].close();
    }

    /**
     * Solve using explicit Runge-Kutta integration method. O(h^4)
     *
     * @param u initial nodal solution vector
     * @param mesh
     * @param config
     */
    void rungeKutta(std::vector<std::vector<double>> &u, Mesh &mesh, Config config, std::vector<std::vector<double>> &pml_phi, std::vector<std::vector<double>> &jca_phi)
    {

        /** Memory allocation */
        elNumNodes = mesh.getElNumNodes();
        numNodes = mesh.getNumNodes();
        elTags = std::vector<int>(&mesh.elTag(0), &mesh.elTag(0) + mesh.getElNum());
        elFlux.resize(elNumNodes);
        elStiffvector.resize(elNumNodes);
        elZKHvector.resize(elNumNodes, 0.0);
        eldampingPvector.resize(elNumNodes, 0.0);
        elAuxiliaryTerm1Vector.resize(elNumNodes, 0.0);
        elAuxiliaryTerm2Vector.resize(elNumNodes, 0.0);
        eljcaAuxiliaryTerm1Vector.resize(elNumNodes, 0.0);
        eljcaAuxiliaryTerm2Vector.resize(elNumNodes, 0.0);
        eljcaPhiVector.resize(elNumNodes, 0.0);
        elSourceVector.resize(elNumNodes, 0.0);
        std::vector<std::vector<double>> k1, k2, k3, k4;
        std::vector<std::vector<double>> h1, h2, h3, h4;
        std::vector<std::vector<double>> j1, j2, j3, j4;
        Flux = std::vector<std::vector<std::vector<double>>>(4, std::vector<std::vector<double>>(mesh.getNumNodes(), std::vector<double>(3)));

        /** Gmsh save init */
        gmsh::model::list(g_names);
        int gp_viewTag = gmsh::view::add("Pressure");
        int gv_viewTag = gmsh::view::add("Velocity");
        int grho_viewTag = gmsh::view::add("Density");
        std::vector<std::vector<double>> g_p(mesh.getElNum(), std::vector<double>(elNumNodes));
        std::vector<std::vector<double>> g_rho(mesh.getElNum(), std::vector<double>(elNumNodes));
        std::vector<std::vector<double>> g_v(mesh.getElNum(), std::vector<double>(3 * elNumNodes));

        /** Precomputation (constants over time) */
        screen_display::write_string("\t>>> Precomputation", BLUE);
        mesh.precomputeMassMatrix();
        screen_display::write_string("\t>>> precomputeMassMatrix", BLUE);

        /** ComputeSigma */
        mesh.computeSigma();

        /** Source */
        std::vector<std::vector<int>> srcIndices;
        for (int i = 0; i < config.sources.size(); ++i)
        {
            std::vector<int> indice;
            for (int n = 0; n < mesh.getNumNodes(); n++)
            {
                std::vector<double> coord, paramCoord;
                int _dim, _tag;
                gmsh::model::mesh::getNode(mesh.getElNodeTags()[n], coord, paramCoord, _dim, _tag);
                if (pow(coord[0] - config.sources[i].source[1], 2) +
                        pow(coord[1] - config.sources[i].source[2], 2) +
                        pow(coord[2] - config.sources[i].source[3], 2) <
                    pow(config.sources[i].source[4], 2))
                {
                    indice.push_back(n);
                }
            }
            srcIndices.push_back(indice);
        }

        /** Observer */
        std::vector<std::vector<int>> obsIndices;
        std::vector<std::vector<double>> obsPtDistance;
        for (int i = 0; i < config.observers.size(); ++i)
        {
            std::vector<int> indice;
            std::vector<double> dist;
            for (int n = 0; n < mesh.getNumNodes(); n++)
            {
                std::vector<double> coord, paramCoord;
                int _dim, _tag;
                gmsh::model::mesh::getNode(mesh.getElNodeTags()[n], coord, paramCoord, _dim, _tag);
                double distance = sqrt(pow(coord[0] - config.observers[i][0], 2) +
                                       pow(coord[1] - config.observers[i][1], 2) +
                                       pow(coord[2] - config.observers[i][2], 2));
                if (distance < config.observers[i][3])
                {
                    indice.push_back(n);
                    dist.push_back(distance);
                }
            }
            obsIndices.push_back(indice);
            obsPtDistance.push_back(dist);
        }

        for (int i = 0; i < obsIndices.size(); i++)
        {
            std::cout << "Observer " << i
                    << " @ (" << config.observers[i][0] << ", "
                    << config.observers[i][1] << ", "
                    << config.observers[i][2] << ")"
                    << "  radius = " << config.observers[i][3]
                    << "  → contains " << obsIndices[i].size() << " nodes" << std::endl;

            for (int j = 0; j < obsIndices[i].size(); j++)
            {
                int nodeIndex = obsIndices[i][j];
                std::vector<double> coord, paramCoord;
                int dim, tag;
                gmsh::model::mesh::getNode(mesh.getElNodeTags()[nodeIndex], coord, paramCoord, dim, tag);

                std::cout << "    Node " << nodeIndex
                        << " : (" << coord[0] << ", "
                        << coord[1] << ", "
                        << coord[2] << ")"
                        << "  distance = " << obsPtDistance[i][j] << std::endl;
            }
            std::cout << std::endl;
        }
        std::cout << "Press Enter to start the simulation..." << std::endl;
        getchar(); // 暂停程序，按回车继续

        // for (int i = 0; i < obsIndices.size(); i++)
        // {
        //     for (int j = 0; j < obsIndices[i].size(); j++)
        //     {
        //         std::cout << obsIndices[i][j] << "\t";
        //     }
        //     std::cout << std::endl;
        // }
        // getchar();

        /**
         * Main Loop : Time iteration
         */

        std::ofstream outfile("residuals.csv");
        outfile << "time;res_p;res_rho;res_vx;res_vy;res_vz;elapsed_time" << std::endl;

        std::vector<std::ofstream> obs_outfile(config.observers.size());
        data4wave.clear();
        data4wave.resize(config.observers.size());
        for (int obs = 0; obs < config.observers.size(); ++obs)
        {
            std::string filename = "results/observers" + std::to_string(obs + 1) + ".txt";
            obs_outfile[obs].open(filename.c_str());
            obs_outfile[obs] << "time;density;pressure;velocity_x;velocity_y;velocity_z" << std::endl;
        }

        auto start = std::chrono::system_clock::now();
        for (double t = config.timeStart, step = int(config.timeStart / config.timeStep), tDisplay = 0; t <= config.timeEnd + config.timeStep / 2;
             t += config.timeStep, tDisplay += config.timeStep, ++step)
        {
            auto start_time = std::chrono::system_clock::now();
            std::vector<double> residual(5, 0.0);
            /**
             *  Savings and prints
             */
            if (tDisplay >= config.timeRate - config.timeStep / 2 || step == 0)  // 设置 config.timeRate = 0.02 , 但实际上是 0.20000000000000001
            {
                tDisplay = 0;

/** [1] Copy solution to match GMSH format */
// #pragma omp parallel for
#pragma omp parallel for schedule(static) num_threads(config.numThreads)
                for (int el = 0; el < mesh.getElNum(); ++el)
                {
                    for (int n = 0; n < mesh.getElNumNodes(); ++n)
                    {
                        int elN = el * elNumNodes + n;
                        g_p[el][n] = u[0][elN];
                        g_rho[el][n] = u[0][elN] / (config.c0 * config.c0);
                        g_v[el][3 * n + 0] = u[1][elN];
                        g_v[el][3 * n + 1] = u[2][elN];
                        g_v[el][3 * n + 2] = u[3][elN];
                    }
                }
                // gmsh::view::addModelData(gp_viewTag, step, g_names[0], "ElementNodeData", elTags, g_p, t, 1);
                // gmsh::view::addModelData(grho_viewTag, step, g_names[0], "ElementNodeData", elTags, g_rho, t, 1);
                // gmsh::view::addModelData(gv_viewTag, step, g_names[0], "ElementNodeData", elTags, g_v, t, 3);

                /** [2] Print and compute iteration time */
                auto end = std::chrono::system_clock::now();
                auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
                gmsh::logger::write("[" + std::to_string(t) + "/" + std::to_string(config.timeEnd) + "s] Step number : " + std::to_string((int)step) + ", Elapsed time: " + std::to_string(elapsed.count()) + "s");
                // screen_display::write_string("time\t\tres_p\t\tres_rho\t\tres_vx\t\tres_vy\t\tres_vz\t\telapsed time", BOLDBLUE);
                // mesh.writeVTK("result.vtk");
                // std::string vtu_filename = "results/result" + std::to_string((int)step) + ".vtu";
                // mesh.writeVTUb(vtu_filename, u);
                std::string vtu_filename_highOrder = "highorder_results/result" + std::to_string((int)step) + ".vtu";
                mesh.writeVTUb_highOrder(vtu_filename_highOrder, u);
                screen_display::write_string("time\t\tres_p\t\tres_rho\t\tres_vx\t\tres_vy\t\tres_vz\t\telapsed time", BOLDBLUE);
                // mesh.writeVTK("result.vtk",u);
            }

            /**Update Source */
            for (int src = 0; src < config.sources.size(); ++src)
            {
                if (config.sources[src].formula == "" && config.sources[src].data.empty())
                {
                    double amp = config.sources[src].source[5];
                    double freq = config.sources[src].source[6];
                    double phase = config.sources[src].source[7];
                    double duration = config.sources[src].source[8];
                    if (t < duration)
                        for (int n = 0; n < srcIndices[src].size(); ++n)
                            u[0][srcIndices[src][n]] = amp * sin(2 * M_PI * freq * t + phase);
                }
                else
                {
                    if (config.sources[src].data.empty())
                    {
                        double duration = config.sources[src].source[5];
                        if (t < duration)
                            for (int n = 0; n < srcIndices[src].size(); ++n)
                                u[0][srcIndices[src][n]] = config.sources[src].value(t);
                    }
                    else
                    {
                        for (int n = 0; n < srcIndices[src].size(); ++n)
                            u[0][srcIndices[src][n]] = config.sources[src].interpolate_value(t);
                    }
                }
            }


            /**
             * Fourth order Runge-Kutta algorithm
             * u ──> [k1] ──┐
             *             ├─> k2 = u + 0.5*k1 ──┐
             *             │                    ├─> k3 = u + 0.5*k2 ──┐
             *             │                    │                    ├─> k4 = u + 1.0*k3 ──┐
             *             │                    │                    │                     │       
             *             ▼                    ▼                    ▼                     ▼
             *          numStep              numStep              numStep               numStep
             *             │                    │                    │                    │
             *             ▼                    ▼                    ▼                    ▼
             *          Flux[k1]            Flux[k2]            Flux[k3]              Flux[k4]
             *             │                    │                    │                    │
             *             └──────┬─────────────┴────────────────────┴────────────────────┘
             *                    ▼
             *          u += (k1 + 2*k2 + 2*k3 + k4)/6
             */
            k1 = k2 = k3 = k4 = u;
            h1 = h2 = h3 = h4 = pml_phi;
            j1 = j2 = j3 = j4 = jca_phi;
            /** [1] Step R-K */
            mesh.updateFlux(k1, Flux, config.v0, config.c0, config.rho0);
            // numStep(mesh, config, k1, Flux, 0);
            // mesh.updatezkFlux(k1, Flux, config.v0, config.c0, config.rho0, config.porousParams[0][0], config.porousParams[0][1], config.porousParams[0][2], config.porousParams[0][3]);
            // zknumStep(mesh, config, k1, Flux, 0);
            // pmlnumStep(mesh, config, k1, h1, Flux, 0);
            // pmlzknumStep(mesh, config, k1, h1, Flux, 0);
            // jcanumStep(mesh, config, k1, j1, Flux, 0);
            // pmljcanumStep(mesh, config, k1, j1, h1, Flux, 0);
            sourcenumStep(mesh, config, k1, Flux, 0, t);
            for (int eq = 0; eq < u.size(); ++eq)
                eigen::plusTimes(k2[eq].data(), k1[eq].data(), 0.5, numNodes);
            for (int eq = 0; eq < pml_phi.size(); eq++)
                eigen::plusTimes(h2[eq].data(), h1[eq].data(), 0.5, numNodes);
            for (int eq = 0; eq < jca_phi.size(); eq++)
            {
                eigen::plusTimes(j2[eq].data(), j1[eq].data(), 0.5, numNodes);
            }
            
            /** [2] Step R-K */
            mesh.updateFlux(k2, Flux, config.v0, config.c0, config.rho0);
            numStep(mesh, config, k2, Flux, 0);
            // mesh.updatezkFlux(k2, Flux, config.v0, config.c0, config.rho0, config.porousParams[0][0], config.porousParams[0][1], config.porousParams[0][2], config.porousParams[0][3]);
            // zknumStep(mesh, config, k2, Flux, 0);
            // pmlnumStep(mesh, config, k2, h2, Flux, 0);
            // pmlzknumStep(mesh, config, k2, h2, Flux, 0);
            // jcanumStep(mesh, config, k2, j2, Flux, 0);
            // pmljcanumStep(mesh, config, k2, j2, h2, Flux, 0);
            sourcenumStep(mesh, config, k2, Flux, 0, t);
            for (int eq = 0; eq < u.size(); ++eq)
                eigen::plusTimes(k3[eq].data(), k2[eq].data(), 0.5, numNodes);
            for (int eq = 0; eq < pml_phi.size(); eq++)
                eigen::plusTimes(h3[eq].data(), h2[eq].data(), 0.5, numNodes);
            for (int eq = 0; eq < jca_phi.size(); eq++)
            {
                eigen::plusTimes(j3[eq].data(), j2[eq].data(), 0.5, numNodes);
            }
            
            /** [3] Step R-K */
            mesh.updateFlux(k3, Flux, config.v0, config.c0, config.rho0);
            // numStep(mesh, config, k3, Flux, 0);
            // mesh.updatezkFlux(k3, Flux, config.v0, config.c0, config.rho0, config.porousParams[0][0], config.porousParams[0][1], config.porousParams[0][2], config.porousParams[0][3]);
            // zknumStep(mesh, config, k3, Flux, 0);
            // pmlnumStep(mesh, config, k3, h3, Flux, 0);
            // pmlzknumStep(mesh, config, k3, h3, Flux, 0);
            // jcanumStep(mesh, config, k3, j3, Flux, 0);
            // pmljcanumStep(mesh, config, k3, j3, h3, Flux, 0);
            sourcenumStep(mesh, config, k3, Flux, 0, t);
            for (int eq = 0; eq < u.size(); ++eq)
                eigen::plusTimes(k4[eq].data(), k3[eq].data(), 1, numNodes);
            for (int eq = 0; eq < pml_phi.size(); eq++)
                eigen::plusTimes(h4[eq].data(), h3[eq].data(), 1, numNodes);
            for (int eq = 0; eq < jca_phi.size(); eq++)
            {
                eigen::plusTimes(j4[eq].data(), j3[eq].data(), 1, numNodes);
            }
            
            /** [4] Step R-K */
            mesh.updateFlux(k4, Flux, config.v0, config.c0, config.rho0);
            // numStep(mesh, config, k4, Flux, 0);
            // mesh.updatezkFlux(k4, Flux, config.v0, config.c0, config.rho0, config.porousParams[0][0], config.porousParams[0][1], config.porousParams[0][2], config.porousParams[0][3]);
            // zknumStep(mesh, config, k4, Flux, 0);
            // pmlnumStep(mesh, config, k4, h4, Flux, 0);
            // pmlzknumStep(mesh, config, k4, h4, Flux, 0);
            // jcanumStep(mesh, config, k4, j4, Flux, 0);
            // pmljcanumStep(mesh, config, k4, j4, h4, Flux, 0);
            sourcenumStep(mesh, config, k4, Flux, 0, t);
            /** Concat results of R-K iterations */
            // #pragma omp parallel for
            for (int eq = 0; eq < u.size(); ++eq)
            {
                for (int i = 0; i < numNodes; ++i)
                {
                    // #pragma omp atomic
                    u[eq][i] += (k1[eq][i] + 2 * k2[eq][i] + 2 * k3[eq][i] + k4[eq][i]) / 6.0;
                }
            }
            for (int eq = 0; eq < pml_phi.size(); eq++)
            {
                for (int i = 0; i < numNodes; ++i)
                {
                    // #pragma omp atomic
                    pml_phi[eq][i] += (h1[eq][i] + 2 * h2[eq][i] + 2 * h3[eq][i] + h4[eq][i]) / 6.0;
                }
            }
            for (int eq = 0; eq < jca_phi.size(); eq++)
            {
                for (int i = 0; i < numNodes; ++i)
                {
                    // #pragma omp atomic
                    jca_phi[eq][i] += (j1[eq][i] + 2 * j2[eq][i] + 2 * j3[eq][i] + j4[eq][i]) / 6.0;
                }
            }
            

#pragma omp parallel for schedule(static) num_threads(config.numThreads)
            for (int el = 0; el < mesh.getElNum(); ++el)
            {
                for (int n = 0; n < mesh.getElNumNodes(); ++n)
                {
                    int elN = el * elNumNodes + n;
#pragma omp atomic update
                    residual[0] += pow(g_p[el][n] - u[0][elN], 2);
#pragma omp atomic update
                    residual[1] += pow(g_rho[el][n] - u[0][elN] / (config.c0 * config.c0), 2);
#pragma omp atomic update
                    residual[2] += pow(g_v[el][3 * n + 0] - u[1][elN], 2);
#pragma omp atomic update
                    residual[3] += pow(g_v[el][3 * n + 1] - u[2][elN], 2);
#pragma omp atomic update
                    residual[4] += pow(g_v[el][3 * n + 2] - u[3][elN], 2);
                }
            }
            outfile << t + config.timeStep << ";";
            std::cout << std::scientific << t + config.timeStep << "\t";
            auto end_time = std::chrono::system_clock::now();
            auto elapsed_time = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
            for (int eq = 0; eq < residual.size(); ++eq)
            {
                residual[eq] /= (mesh.getElNum() * mesh.getElNumNodes());
                std::cout << std::scientific << residual[eq] << "\t";
                outfile << residual[eq] << ";";
            }
            std::cout << elapsed_time.count() * 1.0e-6 << " s" << std::endl;
            outfile << elapsed_time.count() * 1.0e-6 << std::endl;
            /**
             * get observers value
             * Inverse distance weight interpolation method
             */
            for (int obs = 0; obs < config.observers.size(); ++obs)
            {
                double p(0), rho(0), w_sum(0);
                std::vector<double> v = {0, 0, 0};
                for (int n = 0; n < obsIndices[obs].size(); ++n)
                {
                    double R = config.observers[obs][3];                        //! influence sphere
                    double w = 1.0 / (pow(obsPtDistance[obs][n], 2) + 1.0e-12); // fmax(1.0-obsPtDistance[obs][n]/R,0.0);//1.0 / pow(obsPtDistance[obs][n],2);
                    p += u[0][obsIndices[obs][n]] * w;
                    v[0] += u[1][obsIndices[obs][n]] * w;
                    v[1] += u[2][obsIndices[obs][n]] * w;
                    v[2] += u[3][obsIndices[obs][n]] * w;
                    w_sum += w;
                }
                p /= w_sum;
                rho = p / pow(config.c0, 2);
                v[0] /= w_sum;
                v[1] /= w_sum;
                v[2] /= w_sum;
                data4wave[obs].push_back(p);
                obs_outfile[obs] << t << ";" << rho << ";" << p << ";" << v[0] << ";" << v[1] << ";" << v[2] << std::endl;
            }
        }
        for (int obs = 0; obs < config.observers.size(); ++obs)
        {
            io::writeWave(data4wave[obs], "results/observer_" + std::to_string(obs + 1) + ".wav", 1.0 / config.timeStep, 16, 1, 1);
            io::writeFFT(data4wave[obs],config.timeStep,"results/observer_" + std::to_string(obs + 1));
        }

        outfile.close();
        for (int obs = 0; obs < config.observers.size(); ++obs)
            obs_outfile[obs].close();
    }
}
