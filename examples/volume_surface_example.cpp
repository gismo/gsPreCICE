/** @file bSplineSurface_example.cpp

    @brief Tutorial on gsTensorBSpline class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Speh
*/

#include <cmath>
#include <iostream>

#include <gismo.h>


using namespace gismo;

const double PI = 3.14159265;

int main(int argc, char* argv[])
{
    index_t n = 5;
    index_t m = 5;
    index_t degree = 3;
    std::string output("");

    gsCmdLine cmd("Tutorial on gsTensorBSpline class.");
    cmd.addInt   ("n", "dof1", "Number of basis function in one direction"  , n);
    cmd.addInt   ("m", "dof2", "Number of basis function in other direction", m);
    cmd.addInt   ("d", "degree", "Degree of a surface", degree);
    cmd.addString("o", "output", "Name of the output file.", output);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Adjust values to the minimum required
    degree = math::max( (index_t)(0), degree    );
    n      = math::max(n, degree + 1);
    m      = math::max(m, degree + 1);

    gsInfo << "----------------------\n\n"
              << "n: " << n << "\n\n"
              << "m: " << m << "\n\n"
              << "degree: " << degree << "\n\n"
              << "output: " << output << "\n\n"
              << "----------------------\n\n";

    // 1. construction of a knot vector for each direction
    gsKnotVector<> kv1(0, 1, n - degree - 1, degree + 1);
    gsKnotVector<> kv2(0, 1, m - degree - 1, degree + 1);

    // 2. construction of a basis
    gsTensorBSplineBasis<2, real_t> basis(kv1, kv2);

    // 3. construction of a coefficients
    gsMatrix<> greville = basis.anchors();
    gsMatrix<> coefs (greville.cols(), 3);

    for (index_t col = 0; col != greville.cols(); col++)
    {
        real_t x = greville(0, col);
        real_t y = greville(1, col);

        coefs(col, 0) = x;
        coefs(col, 1) = y;
        coefs(col, 2) = math::sin(x * 2 * PI) * math::sin(y * 2 * PI);
    }

    // 4. putting basis and coefficients toghether
    gsTensorBSpline<2, real_t>  surface(basis, coefs);

    gsExprEvaluator<> ev;
    ev.setIntegrationElements(basis);
    auto G = ev.getMap(surface);

    // Make the volumetric basis
    gsKnotVector<> kv3(0, 1, 0, 3);
    gsTensorBSplineBasis<3, real_t> vbasis(kv1, kv2, kv3);
    gsDebugVar(vbasis);
    gsMatrix<> coefs3D(3*greville.cols(), 3);

    gsFunctionExpr<> thickness("0.1+0.025*sin(2*pi*x)*sin(2*pi*y)",2);
    gsMatrix<> N;
    gsMatrix<> T = thickness.eval(greville);
    gsMatrix<> dcoefs(greville.cols(), 3);
    dcoefs.setZero();//not needed
    for (index_t k=0; k!=greville.cols(); k++)
    {
        N = ev.eval(sn(G).normalized(),greville.col(k));
        dcoefs.row(k).transpose() = N*T(0,k);
    }
    coefs3D.block(0                ,0,greville.cols(),3) = coefs-dcoefs;
    coefs3D.block(greville.cols()  ,0,greville.cols(),3) = coefs;
    coefs3D.block(2*greville.cols(),0,greville.cols(),3) = coefs+dcoefs;



    coefs3D << coefs-dcoefs, coefs, coefs+dcoefs;
    gsDebugVar(coefs3D);
    gsGeometry<>::uPtr volume = vbasis.makeGeometry(coefs3D);
    gsWriteParaview(*volume, "volume",1000,true);


    gsWriteParaviewPoints(coefs3D, "coefs3D");


    // direction to eliminate
    short_t dir = 2;
    // ASSUMES THICKNESS IS ALWAYS IN DIR2
    // for (index_t k=0; k!=volume.coefs().rows()/volume.basis().component(2).size(); k++)

    index_t componentSize = volume->coefs().rows()/3;
    gsMatrix<> coefsSurf(3,componentSize);
    coefsSurf.setZero();
    gsMatrix<> dcoefsSurf(3,componentSize);



    gsDebugVar(componentSize);

    gsDebugVar(coefsSurf.dim());

    gsDebugVar(N);



    for (index_t k=0; k!=componentSize; k++)
    {
        coefsSurf.col(k) = (coefs3D.row(k)+coefs3D.row(k+2*componentSize))/2;
        dcoefsSurf.col(k) = (coefs3D.row(k+2*componentSize)-coefs3D.row(k))/2; // think about the normal
    }
    gsWriteParaviewPoints(coefsSurf, "coefsSurf");
    gsWriteParaviewPoints(dcoefsSurf, "dcoefsSurf");



    // Use the same basis for the mid surface
    gsGeometry<>::uPtr surf_mid = basis.makeGeometry(coefsSurf.transpose());
    gsMultiPatch<> mid_surface_geom;
    mid_surface_geom.addPatch(std::move(surf_mid)); // Use std::move to transfer ownership

    gsGeometry<>::uPtr surfThickness = basis.makeGeometry(dcoefsSurf.transpose());



    gsExprEvaluator<> ev_shell;
    ev_shell.setIntegrationElements(basis);
    auto G_shell = ev_shell.getMap(mid_surface_geom);

    gsMatrix<> N_shell;

    gsMatrix<> T_shell = thickness.eval(greville);

    gsMatrix<> deformed_thickness(greville.cols(), 3);
    dcoefs.setZero();//not needed
    for (index_t k=0; k!=greville.cols(); k++)
    {
        N_shell = ev_shell.eval(sn(G_shell).normalized(),greville.col(k));
        deformed_thickness.row(k).transpose() = N_shell*T_shell(0,k);
    }

    gsDebugVar(deformed_thickness);

    gsMatrix<> deformedcoefs3D(3*greville.cols(), 3);

    gsDebugVar(deformed_thickness.dim());

    gsDebugVar(deformed_thickness.dim());
    gsDebugVar(coefsSurf.dim());

    coefsSurf.transposeInPlace();

    deformedcoefs3D.block(0                ,0,greville.cols(),3) = coefsSurf - deformed_thickness;
    deformedcoefs3D.block(greville.cols()  ,0,greville.cols(),3) = coefsSurf;
    deformedcoefs3D.block(2*greville.cols(),0,greville.cols(),3) = coefsSurf + deformed_thickness;


    deformedcoefs3D << coefsSurf-deformed_thickness, coefsSurf, coefsSurf+deformed_thickness;

    gsDebugVar(deformedcoefs3D);
    gsGeometry<>::uPtr deformed_volume = vbasis.makeGeometry(deformedcoefs3D);
    
    gsWriteParaview(*deformed_volume, "deformed_volume",1000,true);

    coefs3D.transposeInPlace();
    gsWriteParaviewPoints(deformedcoefs3D, "coefs3D");



    gsWriteParaview(*surf_mid, "surf_recon",1000,true);
    gsWriteParaview(*surfThickness, "surfThickness_recon",1000,true);

    std::vector<gsBasis<>::uPtr> bases(2);


    // for (short_t d = 0; d!=vbasis.domainDim(); d++)
    //     if (d!=dir)
    //         bases[d] = vbasis.component(d);
    // gsTensorBSplineBasis<2,real_t> sbasis(*basis[0],*basis[1]);



    // 5. saving surface, basis and control net to a file
    if (output != "")
    {
        std::string out = output + "Geometry";
        gsInfo << "Writing the surface to a paraview file: " << out
                  << "\n\n";

        gsWriteParaview(surface, out, 10000);

        out = output + "Basis";
        gsInfo << "Writing the basis to a paraview file: " << out
                  << "\n\n";

        gsWriteParaview(basis, out);


        out = output + "ContolNet";
        gsInfo << "Writing the control net to a paraview file: " << out
                  << "\n" << "\n";

        gsMesh<> mesh;
        surface.controlNet(mesh);
        gsWriteParaview(mesh, out);

        out = output + "Coefficients";
        gsMatrix <> coefs = surface.coefs();
        coefs.transposeInPlace();
        gsWriteParaviewPoints(coefs, out);

    }
    else
    {
        gsInfo << "Done. No output created, re-run with --output <filename> to get a ParaView "
                  "file containing the solution.\n";
    }

    return 0;
}
