/** @file gsPreCICEUtils.h

    @brief Utilities file for using gsPreCICE extension

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (TU Delft, 2019-2024), J.Li (TU Delft, 2023 - ...)
*/

#pragma once

#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsNurbs/gsTensorNurbsBasis.h>

namespace gismo {

/**
 * @brief      Gets a vector of knot vectors from a basis
 *
 * @param[in]  source      The source basis
 * @param      tensorKnots  The knots in a vector
 *
 * @tparam     DIM         Dimension
 * @tparam     T           Real type
 * @tparam     basis_type  Basis type
 */
template <short_t DIM, class T, template<short_t _DIM, class _T> class basis_type>
inline void getKnots(const gsBasis<T> & source, std::vector<gsKnotVector<T>> & tensorKnots)
{
    if ( const basis_type<DIM,T> * basis = dynamic_cast<const basis_type<DIM,T>*>(&source) )
        for (index_t d=0; d!=DIM; d++)
            tensorKnots[d] = basis->knots(d);
}

/**
 * @brief      Puts all the knots in a vector (for 1D bases only)
 *
 *  Make a matrix with the knot vectors
 *  [[x1_1, x1_2, ..., nan , nan ,..., nan , nan , ...]
 *   [nan , nan , ..., x2_1, x2_2,..., nan , nan , ...]
 *   [nan , nan , ..., nan , nan ,..., x2_1, x2_2, ...]]
 *
 * @param[in]  basis  The basis
 *
 * @tparam     T      Real type
 *
 * @note @hverhelst, @Crazy-Rich-Meghan is this needed?
 */

template<class T>
gsVector<T> knotsToVector(const gsBasis<T> & basis)
{
    const size_t DIM = 1;
    std::vector<gsKnotVector<T>> tensorKnots(DIM);
    switch (DIM)
    {
        case 1:
            getKnots<1,T,gsTensorBSplineBasis>(basis,tensorKnots);
            getKnots<1,T,gsTensorNurbsBasis>(basis,tensorKnots);
            break;
        default:
            GISMO_ERROR("Basis type not understood");
            break;
    }

    std::vector<size_t> sizes(DIM);
    sizes[0] = tensorKnots[0].size();

    gsVector<T> knots(std::accumulate(sizes.begin(),sizes.end(),0));

    for (index_t i = 0; i < knots.size(); ++i)
        knots[i] = tensorKnots[0][i];
    return knots;
}

/** @brief      Puts all the knots in a matrix (for n-D bases)
 *
 *  Make a matrix with the knot vectors
 *  [[x1_1, x1_2, ..., nan , nan ,..., nan , nan , ...]
 *   [nan , nan , ..., x2_1, x2_2,..., nan , nan , ...]
 *   [nan , nan , ..., nan , nan ,..., x2_1, x2_2, ...]]
 *
 * @param[in]  basis  The basis
 *
 * @tparam     T      Real type
 */
template<class T>
gsMatrix<T> knotsToMatrix(const gsBasis<T> & basis)
{
    const size_t DIM = basis.domainDim();
    std::vector<gsKnotVector<T>> tensorKnots(DIM);
    switch (DIM)
    {
        case 1:
            getKnots<1,T,gsTensorBSplineBasis>(basis,tensorKnots);
            getKnots<1,T,gsTensorNurbsBasis>(basis,tensorKnots);
            break;
        case 2:
            getKnots<2,T,gsTensorBSplineBasis>(basis,tensorKnots);
            getKnots<2,T,gsTensorNurbsBasis>(basis,tensorKnots);
            break;
        case 3:
            getKnots<3,T,gsTensorBSplineBasis>(basis,tensorKnots);
            getKnots<3,T,gsTensorNurbsBasis>(basis,tensorKnots);
            break;
        default:
            GISMO_ERROR("Basis type not understood");
            break;
    }

    std::vector<size_t> sizes(DIM);
    std::vector<size_t> strides(DIM);
    for (size_t d=0; d!=DIM; d++)
        sizes[d] = tensorKnots[d].size();

    strides[0]=0;
    for (size_t d=1; d!=DIM; d++)
        strides[d]= strides[d-1]+sizes[d-1];

    gsMatrix<T> knots(DIM,std::accumulate(sizes.begin(),sizes.end(),0));
    knots.setConstant(std::nan("1"));
    for (size_t d=0; d!=DIM; d++)
        knots.block(d,strides[d],1,sizes[d]) = tensorKnots[d].asMatrix();

    return knots;
}

/**
 * @todo @hverhelst, @Crazy-Rich-Meghan
 */
template<class T>
gsMatrix<T> knotVectorUnpack(const gsMatrix<T> & knots, index_t numBoundaries)
{
    gsMatrix<> kv_unpacked;
    kv_unpacked = knots.row(0);
    kv_unpacked.resize(knots.cols()/numBoundaries,numBoundaries);

    return kv_unpacked;
}




// template <class T>
// std::pair<gsMatrix<T>, gsMatrix<T>> packMultiPatch(const gsMultiPatch<T> &mp) {
//     std::vector<gsMatrix<T>> knotMatrices;
//     knotMatrices.reserve(mp.nPatches());
//     std::vector<gsMatrix<T>> coefMatrices;
//     coefMatrices.reserve(mp.nPatches());

//     for (typename gsMultiPatch<T>::const_iterator patch = mp.begin(); patch != mp.end(); ++patch) {
//         // Dereference the patch pointer to access its members
//         knotMatrices.push_back(knotsToMatrix((*patch)->basis()));
//         coefMatrices.push_back((*patch)->coefs().transpose());
//     }

//     index_t knotRows = mp.domainDim();
//     index_t coefRows = mp.targetDim();
//     index_t knotCols = 0;
//     index_t coefCols = 0;

//     for (size_t p = 0; p != mp.nPatches(); ++p) {
//         knotCols += knotMatrices[p].cols();
//         coefCols += coefMatrices[p].cols();
//     }

//     gsMatrix<T> knotMatrix(knotRows, knotCols);
//     gsMatrix<T> coefMatrix(coefRows, coefCols);

//     size_t currentKnotCol = 0;
//     size_t currentCoefCol = 0;

//     for (size_t p = 0; p != mp.nPatches(); ++p) {
//         const auto& km = knotMatrices[p];
//         const auto& cm = coefMatrices[p];

//         // Validate block dimensions before assignment
//         if (currentKnotCol + km.cols() <= knotCols && km.rows() == knotRows &&
//             currentCoefCol + cm.cols() <= coefCols && cm.rows() == coefRows) {
//             knotMatrix.block(0, currentKnotCol, knotRows, km.cols()) = km;
//             coefMatrix.block(0, currentCoefCol, coefRows, cm.cols()) = cm;

//             currentKnotCol += km.cols();
//             currentCoefCol += cm.cols();
//         } else {
//             std::cerr << "Invalid block parameters for patch " << p << std::endl;
//             throw std::logic_error("Block parameters out of range");
//         }
//     }

//     return std::make_pair(knotMatrix, coefMatrix);
// }

/**
 * @brief      Pack the knot and control points matrices of a gsMultiPatch object into a single matrix.
 *
 * @param[in]  mp  The gsMultiPatch object. MultiPatch geometry object.
 *
 * @tparam     T   Real type.
 *
 * @return     A tuple containing the packed knot matrix, the packed ratio matrix, the number of columns in each patch's knot matrix, and the number of columns in each patch's ratio matrix.
 */

template <class T>
std::tuple<gsMatrix<T>, gsMatrix<T>, std::vector<index_t>, std::vector<index_t>> packMultiPatch(const gsMultiPatch<T> &mp)
{
    std::vector<gsMatrix<T>> knotMatrices;
    knotMatrices.reserve(mp.nPatches());
    std::vector<gsMatrix<T>> coefMatrices;
    coefMatrices.reserve(mp.nPatches());
    std::vector<index_t> knotCols;
    std::vector<index_t> coefCols;

    for (typename gsMultiPatch<T>::const_iterator patch = mp.begin(); patch != mp.end(); ++patch) {
        // Dereference the patch pointer to access its members
        knotMatrices.push_back(knotsToMatrix((*patch)->basis()));
        coefMatrices.push_back((*patch)->coefs().transpose());
        knotCols.push_back(knotMatrices.back().cols());
        coefCols.push_back(coefMatrices.back().cols());
    }

    index_t knotRows = mp.domainDim();
    index_t coefRows = mp.targetDim();
    index_t knotColsSum = 0;
    index_t coefColsSum = 0;

    for (size_t p = 0; p != mp.nPatches(); ++p) {
        knotColsSum += knotMatrices[p].cols();
        coefColsSum += coefMatrices[p].cols();
    }

    gsMatrix<T> knotMatrix(knotRows, knotColsSum);
    gsMatrix<T> coefMatrix(coefRows, coefColsSum);

    size_t currentKnotCol = 0;
    size_t currentCoefCol = 0;

    for (size_t p = 0; p != mp.nPatches(); ++p) {
        const auto& km = knotMatrices[p];
        const auto& cm = coefMatrices[p];

        // Validate block dimensions before assignment
        if (currentKnotCol + km.cols() <= knotColsSum && km.rows() == knotRows &&
            currentCoefCol + cm.cols() <= coefColsSum && cm.rows() == coefRows) {
            knotMatrix.block(0, currentKnotCol, knotRows, km.cols()) = km;
            coefMatrix.block(0, currentCoefCol, coefRows, cm.cols()) = cm;

            currentKnotCol += km.cols();
            currentCoefCol += cm.cols();
        } else {
            std::cerr << "Invalid block parameters for patch " << p << std::endl;
            throw std::logic_error("Block parameters out of range");
        }
    }

    return std::make_tuple(knotMatrix, coefMatrix, knotCols, coefCols);
}


/**
 * @brief      Convert the packed matrices into back into a gsMultiPatch<T> object.
 *
 * @param[in]  knotMatrices  The knot matrices
 * @param[in]  coefMatrices  The ratio matrices
 *
 * @tparam     T             Reak type.
 *
 * @return     A gsMultiPatch object.
 */


template <class T>
gsMultiPatch<T> unpackMultiPatch(const gsMatrix<T> &knotMatrix, const gsMatrix<T> &coefMatrix, const std::vector<index_t> &knotCols, const std::vector<index_t> &coefCols) {
    gsMultiPatch<T> mp;
    size_t currentKnotCol = 0;
    size_t currentCoefCol = 0;

    for (size_t p = 0; p < knotCols.size(); ++p) {
        gsMatrix<T> km = knotMatrix.block(0, currentKnotCol, knotMatrix.rows(), knotCols[p]);
        gsMatrix<T> cm = coefMatrix.block(0, currentCoefCol, coefMatrix.rows(), coefCols[p]).transpose();

        // Create the gsBasis object from the knot matrix
        std::shared_ptr<gsBasis<T>> KnotBasis = knotMatrixToBasis(km);

        // Create the gsGeometry object
        auto geom = KnotBasis->makeGeometry(cm);

        // Add the gsGeometry object to gsMultiPatch
        mp.addPatch(std::move(geom));

        currentKnotCol += knotCols[p];
        currentCoefCol += coefCols[p];
    }

    return mp;
}

/**
 * @brief      Unpack the control points matrix of a gsMultiPatch object into separate matrices.
 *
 * @param[in]  controlPoints  The control points matrix
 * @param[in]  kv_unpacked    The unpacked knot matrix
 * @param[in]  knot_index     The index of the knot matrix
 * @param[in]  cp_index       The index of the control points matrix
 *
 *
 * @return     A matrix of control points.
 */

template<class T>
gsMatrix<T> unPackControlPoints(const gsMatrix<T> & controlPoints, const gsMatrix<T> & kv_unpacked, index_t knot_index, index_t cp_index)
{
    // number of cps n = N_knot - p - 1
    int counter = cp_index;
    gsVector<> temp = kv_unpacked.row(knot_index);

    // Calculate the amount of control points based on the knot vector
    for (index_t i = 0; i < temp.size(); ++i)
    {
        if(temp[i] == 0)
            counter ++;
    }
    counter++;

    gsMatrix<T> unpackedCps(controlPoints.rows(), counter - cp_index);


    // Simplified copying of control points (just a placeholder)
    for (index_t i = 0; i < controlPoints.rows(); ++i)
    {
        index_t diff = counter - cp_index;
        index_t startIndex = cp_index - diff + 1;  // Calculate the start index
        for (size_t j = 0; j < diff; ++j)
        {
            unpackedCps(i,j) = controlPoints(i,j + startIndex);
        }
    }
    return unpackedCps;
}

/**
 * @brief      Convert a matrix of knot vectors into a gsBasis object.
 *
 * @param[in]  knots  The matrix of knot vectors
 *
 * @return     A shared pointer to the gsBasis object.
 */

template<class T>
typename gsBasis<T>::Ptr knotMatrixToBasis(const gsMatrix<T> & knots)
{
    gsBasis<> * basis;
    const short_t DIM = knots.rows();

    std::vector<gsKnotVector<T>> KVs(DIM);
    index_t k=0;
    for (size_t d=0; d!=DIM; d++)
    {
        std::vector<T> tmp;
        std::copy_if(knots.row(d).begin(),knots.row(d).end(),
                     std::back_inserter(tmp),
                     [](T a){return !math::isnan(a);});
        KVs[d] = give(gsKnotVector<T>(tmp));


        gsDebug<<"(gsPreCICEUtils.h: There is a memory leak in the line above)\n";
    }

    switch(DIM)
    {
        case 1:
            basis = new  gsBSplineBasis<T>(KVs[0]);
            break;
        case 2:
            basis = new  gsTensorBSplineBasis<2,T>(KVs);
            break;
        case 3:
            basis = new  gsTensorBSplineBasis<3,T>(KVs);
            break;
    }
    return memory::make_shared_not_owned(basis);
}

/**
 * @brief Creates a 3D volume from a 2D surface by extruding along its normals with variable thickness.
 *
 * This function takes a 2D tensor B-spline surface and extrudes it according to a scalar thickness function
 * to create a 3D tensor B-spline volume. The extrusion is performed along the surface normals, with half
 * of the thickness extending in each direction (positive and negative).
 *
 * @tparam T Numeric type (e.g., real_t, double)
 * @param surface Input tensor B-spline surface (2D parametric domain)
 * @param thickness Function that maps from the surface parameter space to thickness values
 *
 * @return A 3D tensor B-spline volume created by extruding the input surface
 *
 * @note The function supports both 2D and 3D target dimensions for the input surface.
 *       For 2D surfaces, the normal vector is assumed to be (0,0,1).
 *       For 3D surfaces, proper surface normals are computed.
 *
 * @throws AssertionError if the thickness function is not defined on a 2D domain
 * @throws AssertionError if the thickness function is not scalar-valued
 * @throws AssertionError if the surface is not defined on a 2D domain
 * @throws Error if the surface's target dimension is not 2D or 3D
 */
template <class T>
gsTensorBSpline<3, T> surfaceToVolume(const gsTensorBSpline<2, T> & surface, const gsFunction<T> & thickness)
{
    // Validate input parameters
    GISMO_ASSERT(thickness.domainDim() == 2,"Thickness function must be defined on a 2D domain.");
    GISMO_ASSERT(thickness.targetDim() == 1,"Thickness function must be scalar valued.");
    GISMO_ASSERT(surface.domainDim() == 2,"Surface must be defined on a 2D domain.");

    // Create a knot vector for the third dimension (thickness direction)
    // Linear basis with 2 elements: [0,0,0.5,1,1] (degree 1, 2 internal knots)
    gsKnotVector<T> kvZ(0, 1, 1, 2);

    // Create a 3D tensor B-spline basis by combining the existing 2D basis with the new Z knot vector
    gsTensorBSplineBasis<3, T> basis(surface.basis().knots(0),
                                      surface.basis().knots(1),
                                      kvZ);

    // Calculate control point count for the new volume
    const index_t nCoefs = surface.coefs().rows();

    // Allocate control points matrix for the volume (3 layers of control points)
    gsMatrix<T> coefs(3*nCoefs, surface.coefs().cols());

    // Get parametric coordinates (anchors) of the surface control points
    gsMatrix<T> anchors = surface.basis().anchors();

    // Matrix to store thickness values at each anchor point
    gsMatrix<T> thicknessValues;

    // Evaluate the thickness function at each anchor point
    thickness.eval_into(anchors, thicknessValues);

    // Matrix to store surface normals
    gsMatrix<T> normals;

    // Handle different dimension cases for the surface
    if (surface.targetDim() == 2)
    {
        // For 2D surfaces, use (0,0,1) as the normal vector
        normals = gsMatrix<T>::Zero(3, anchors.cols());
        normals.row(2).setConstant(1.0); // Set the z-component of the normal to 1
    }
    else if (surface.targetDim() == 3)
    {
        // For 3D surfaces, compute proper surface normals at each anchor point
        gsMapData<T> md(NEED_NORMAL);
        md.points = anchors;
        surface.computeMap(md);
        std::swap(md.normals, normals); // Get the computed normal vectors
    }
    else
    {
        GISMO_ERROR("The surface must be 2D or 3D.");
    }

    // For each control point of the surface
    for (index_t k = 0; k != anchors.cols(); k++)
    {
        // Normalize the normal vector
        const gsVector<T> & normal = normals.col(k).normalized();

        // Get thickness value at this control point
        const T & t = thicknessValues(0,k);

        // Create three layers of control points:
        // Bottom layer: offset by -0.5*thickness along the normal
        coefs.row(k)            = surface.coefs().row(k) - 0.5*t*normal.transpose();
        // Middle layer: original surface control point
        coefs.row(k + nCoefs)   = surface.coefs().row(k);
        // Top layer: offset by +0.5*thickness along the normal
        coefs.row(k + 2*nCoefs) = surface.coefs().row(k) + 0.5*t*normal.transpose();
    }

    // Create and return the 3D tensor B-spline volume
    return gsTensorBSpline<3, T>(basis, coefs);
}

/**
 * @brief Converts a 2D tensor B-spline surface to a 3D tensor B-spline volume.
 *
 * This function creates a volumetric representation of a surface by extruding the
 * surface in the normal direction with a constant thickness.
 *
 * @tparam T Type for coordinates (typically float or double)
 *
 * @param surface The input 2D tensor B-spline surface
 * @param thickness The constant thickness value for the extrusion (default = 1.0)
 *
 * @return A 3D tensor B-spline volume representing the extruded surface
 *
 * @see toVolume(const gsTensorBSpline<2, T> &, const gsFunction<T> &)
 */
template <class T>
gsTensorBSpline<3, T> toVolume(const gsTensorBSpline<2, T> & surface, const T & thickness = 1.0)
{
    return toVolume(surface,gsConstantFunction<T>(thickness, surface.basis().domainDim()));
}

} //namespace gismo
