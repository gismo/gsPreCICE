/** @file gsLookupFunction.h

    @brief Provides declaration of gsLookupFunction class for multipatch support.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J.Li (TU Delft), H. Verhelst
    */

#pragma once

#include <gsCore/gsMultiPatch.h>


namespace gismo
{
template<class T>
class gsLookupFunctionSingle;
/**
 * @brief      Multipatch lookup function evaluator
 * @details    This class enables lookup functions on multipatch geometries where each patch
 *             can have its own lookup function mapping points to data values
 * @param      T     Number format
 */
template <class T>
class gsLookupFunction : public gsFunctionSet<T>
{
public:
    typedef memory::shared_ptr< gsLookupFunction > Ptr;
    typedef memory::unique_ptr< gsLookupFunction > uPtr;
    typedef gsLookupFunctionSingle<T> Piece_t;

    // Using unique_ptr 
    typedef std::vector<memory::unique_ptr<Piece_t>> Container;

    // Constructor
    gsLookupFunction(size_t nPatches = 0)
    : m_container(nPatches)
    {
        // Use n nullptrs to initialize the container
    }
    /// Constructor with points and data for a single patch
    gsLookupFunction(const gsMatrix<T> & points,
                               const gsMatrix<T> & data)
    : gsLookupFunction(1)
    {
        this->add(points, data);
    }

    /// Destructor
    ~gsLookupFunction()
    {
    }

    // clone function for deep copying the lookup function
    private: 
    virtual gsLookupFunction* clone_impl() const override 
    { 
        gsLookupFunction* result = new gsLookupFunction(m_container.size());
        // Deep copy each element
        for (size_t i = 0; i < m_container.size(); ++i) {
            if (m_container[i]) {
                result->m_container[i].reset(new Piece_t(*m_container[i]));
            }
        }
        return result;
    }
    public:
    memory::unique_ptr<gsLookupFunction> clone() const 
    { 
        return memory::unique_ptr<gsLookupFunction>(clone_impl()); 
    }

    /// Number of pieces
    index_t nPieces() const override
    { return m_container.size(); }

    /// Access a piece
    const gsFunction<T> & piece(const index_t k) const override
    {
        GISMO_ASSERT(k < m_container.size(), "Piece index out of bounds");
        GISMO_ASSERT(m_container[k], "Piece at index " + std::to_string(k) + " is not initialized");
        return *m_container[k];
    }

    /// See \a gsFunctionSet
    index_t size() const override
    { return m_container.size(); }

    /// Domain dimension
    short_t domainDim() const override
    {
        if (!m_container.empty() && m_container[0])
            return m_container[0]->domainDim();
        return 0;
    }

    /// Target dimension
    short_t targetDim() const override
    {
        if (!m_container.empty() && m_container[0])
            return m_container[0]->targetDim();
        return 0;
    }

    /// Evaluate function - requires patch index (no default evaluation across all patches)
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override
    {
        GISMO_ERROR("gsLookupFunction requires patch index for evaluation. Use eval_into(patch, u, result)");
    }

    /// Evaluate derivatives - requires patch index (no default evaluation across all patches)
    void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override
    {
        GISMO_ERROR("gsLookupFunction requires patch index for evaluation. Use deriv_into(patch, u, result)");
    }

    /// Evaluate second derivatives - requires patch index (no default evaluation across all patches)
    void deriv2_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override
    {
        GISMO_ERROR("gsLookupFunction requires patch index for evaluation. Use deriv2_into(patch, u, result)");
    }

    /// Evaluate the function at points u for specific patch
    void eval_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        GISMO_ASSERT(patch < m_container.size(), "Patch index out of bounds");
        m_container[patch]->eval_into(u, result);
    }

    /// Evaluate derivatives for specific patch
    void deriv_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        GISMO_ASSERT(patch < m_container.size(), "Patch index out of bounds");
        m_container[patch]->deriv_into(u, result);
    }

    /// Evaluate second derivatives for specific patch
    void deriv2_into(const index_t patch, const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        GISMO_ASSERT(patch < m_container.size(), "Patch index out of bounds");
        m_container[patch]->deriv2_into(u, result);
    }

    /// Get the lookup function container
    const Container & container() const { return m_container; }

    /// Get the lookup function container (non-const)
    Container  & container() { return m_container; }

    /// Add lookup function (appends to the end)
    void add(const gsMatrix<T> & points, const gsMatrix<T> & data)
    {
        // Use unique_ptr - more efficient than shared_ptr
        m_container.push_back(memory::make_unique(new Piece_t(points, data)));
    }

    /// Set lookup function for a specific patch index
    void set(index_t patch, const gsMatrix<T> & points, const gsMatrix<T> & data)
    {
        GISMO_ASSERT(patch < m_container.size(), "Patch index out of bounds");
        // Reset with unique_ptr - old object is automatically destroyed
        m_container[patch].reset(new Piece_t(points, data));
    }

protected:

    Container m_container;
};

/**
 * @brief      Class defining a function that looks up registered data on points.
 * @usage      Combine with gsThinShellAssembler to update the stress on the solid.
 * @details   The gsLookupFunction enables
 *            efficient data lookups based on spatial coordinates. When given a set of points and corresponding
 *            data, it creates a mapping that allows for quick retrieval of the data based on the point
 *            coordinates.
 * @param     T     Number format
 */

template <class T>
class gsLookupFunctionSingle : public gsFunction<T>
{

    struct Compare
    {
        bool operator()(const gsVector<T> & a, const gsVector<T> & b) const
        {
            return std::lexicographical_compare(  a.begin(), a.end(), b.begin(), b.end());
        }
    };

public:
    typedef gsGeometry<T> Base;

    /// Shared pointer for gsLookupFunction
    typedef memory::shared_ptr< gsLookupFunctionSingle > Ptr;

    /// Unique pointer for gsLookupFunction
    typedef memory::unique_ptr< gsLookupFunctionSingle > uPtr;

    /// Default constructor
    gsLookupFunctionSingle() { 
        // Note: m_points and m_data are empty, update() will be called when data is added
    }

    /**
     * @brief      Constructs a new instance of the gsLookupFunctionSingle
     *
     * @param      interface   The precice::SolverInterface (see \a gsPreCICE)
     * @param[in]  meshName      The ID of the mesh on which the data is located
     * @param[in]  dataName      The ID of the data
     * @param[in]  patches     The geometry
     * @param[in]  parametric  Specifies whether the data is defined on the parametric domain or not
     */
    gsLookupFunctionSingle(   const gsMatrix<T> & points,
                              const gsMatrix<T> & data   )
    :
    m_points(points),
    m_data(data)
    {
        this->update();
    }

    /// Constructs a function pointer
    static uPtr make(   const gsMatrix<T> & points,
                        const gsMatrix<T> & data   )
    { return uPtr(new gsLookupFunctionSingle(points, data)); }

    GISMO_CLONE_FUNCTION(gsLookupFunctionSingle)

    /// Access a piece
    const gsLookupFunctionSingle<T> & piece(const index_t) const
    {
        return *this;
    }

    /// See \a gsFunction
    virtual short_t domainDim() const
    { return m_points.rows(); }

    /// Gives the targetDomain, currently only scalar functions (todo)
    virtual short_t targetDim() const
    { return m_data.rows(); }

    /** \brief Evaluate the function at points \a u into \a result.
     * This function evaluates a target function at a given set of input points and
     * stores the results in the provided output matrix. The mapping between input
     * points and corresponding results is defined by the internal mapping table \c m_map.
     * The function ensures that all input points are valid and registered in the table
     * before performing the evaluation.
     *
     * \param[in] u A \c gsMatrix , where each column represents a point in the parameter domain.
     * \param[out] result A \c gsMatrix , where each column will contain the result of
     * evaluating the function at the corresponding input point.
     *
     */
    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override
    {
        const T tolerance = 100 * math::limits::epsilon();

        result.resize(this->targetDim(), u.cols());
        result.setZero();

        for (index_t k = 0; k != u.cols(); ++k)
        {
            auto it = std::find_if(m_map.begin(), m_map.end(),
                [&](const std::pair<gsVector<T>, index_t>& entry) {
                    return (entry.first - u.col(k)).norm() <= tolerance;
                });

            GISMO_ASSERT(std::count_if(m_map.begin(), m_map.end(),[&](const std::pair<gsVector<T>, index_t>& entry) {return (entry.first - u.col(k)).norm() <= tolerance;}) == 1,"gsLookupFunction: Multiple points within tolerance");

            GISMO_ASSERT(it != m_map.end(),
                "Coordinate " + util::to_string(k) + " [" + util::to_string(u.col(k).transpose()) + "] not registered in the table within tolerance.");
            result.col(k) = m_data.col(it->second);
        }
    }

    /// See \a gsFunction
    virtual void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        // This would be nice to have with higher-order (IGA) coupling of precice
        GISMO_NO_IMPLEMENTATION;
    }

    /// See \a gsFunction
    virtual void deriv2_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        // This would be nice to have with higher-order (IGA) coupling of precice
        GISMO_NO_IMPLEMENTATION;
    }

    /// See \a gsFunction
    void evalAllDers_into(const gsMatrix<T> & u, int n,
                          std::vector<gsMatrix<T> > & result) const
    {
        result.resize(1);
        // This would be nice to have with higher-order (IGA) coupling of precice
        gsMatrix<T> tmp;
        this->eval_into(u,tmp);
        result[0]= tmp;
    }

    void update()
    {
        m_map.clear();
        GISMO_ASSERT(m_points.cols()==m_data.cols(),"Points and data must have the same number of columns");
        for (index_t k = 0; k != m_points.cols(); k++)
            m_map.insert({m_points.col(k),k}); // m_map.at(vector) returns the column index of vector
    }

    /// See \a gsFunction
    virtual std::ostream &print(std::ostream &os) const
    {
        os << "gsLookupFunctionSingle\n";
        return os;
    }

protected:

    const gsMatrix<T>& m_points;
    const gsMatrix<T>& m_data;

    std::map<gsVector<T>,index_t,Compare> m_map;

};

} // namespace gismo

