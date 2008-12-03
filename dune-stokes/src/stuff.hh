/**
 *  \file stuff.hh
 *  \brief  contains some stuff
 **/
#ifndef STUFF_HH_INCLUDED
#define STUFF_HH_INCLUDED

#define SEGFAULT int*i=0;*i=9;

#include <iomanip>


namespace Stuff
{

/**
 *  \todo doc me
 **/
template < class ReturnType >
ReturnType fromString(const std::string& s)
{
    std::stringstream ss;
    ss << s;
    ReturnType r;
    ss >> r;
    return r;
}

/**
 *  \todo doc
 **/
template < class ReturnType >
std::string toString(const ReturnType& s)
{
    std::stringstream ss;
    ss << s;
    std::string r;
    ss >> r;
    return r;
}

/**
 *  \brief Only free mem asc. to valid pointer, log warning otherwise
 *
 **/

template < class T >
void safe_delete ( T t )
{
    if (t){
        delete t;
        t=0;
    }
    //else log warning
}

/**
 *  \brief prints a Dune::FieldVector
 *
 *  or anything compatible in terms of Iterators
 *  \tparam should be Dune::FieldVector or compatible
 *  \param  arg
 *          Vector to be printed
 **/
template < class T >
void printFieldVector( T& arg )
{
    std::cout << "\nprinting Dune::FieldVector" << std::endl;
    typedef typename T::ConstIterator
        IteratorType;
    IteratorType itEnd = arg.end();
    for ( IteratorType it = arg.begin(); it != itEnd; ++it ) {
            std::cout << std::setw( 7 ) << std::setprecision( 3 ) << *it;
    }
}

/**
 *  \brief prints a Dune::FieldVector
 *
 *  or anything compatible in terms of Iterators
 *  \tparam should be Dune::FieldVector or compatible
 *  \param  arg
 *          Vector to be printed
 *  \param  name
 *          name to be printed along
 **/
template < class T >
void printFieldVector( T& arg, std::string name )
{
    std::cout << "\nprinting " << name << " (Dune::FieldVector)" << std::endl;
    typedef typename T::ConstIterator
        IteratorType;
    IteratorType itEnd = arg.end();
    for ( IteratorType it = arg.begin(); it != itEnd; ++it ) {
            std::cout << std::setw( 7 ) << std::setprecision( 3 ) << *it;
    }
}

/**
 *  \brief prints a Dune::FieldMatrix
 *
 *  or anything compatible in terms of Iterators
 *  \tparam should be Dune::FieldMatrix or compatible
 *  \param  arg
 *          Matrix to be printed
 **/
template < class T >
void printFieldMatrix( T& arg )
{
    std::cout << "\nprinting Dune::FieldMatrix";
    typedef typename T::ConstRowIterator
        RowIteratorType;
    typedef typename T::row_type::ConstIterator
        VectorInRowIteratorType;
    unsigned int row = 1;
    RowIteratorType rItEnd = arg.end();
    for ( RowIteratorType rIt = arg.begin(); rIt != rItEnd; ++rIt ) {
        std::cout << "\nrow " << row << ":";
        VectorInRowIteratorType vItEnd = rIt->end();
        for (   VectorInRowIteratorType vIt = rIt->begin(); vIt != vItEnd; ++vIt ) {
            std::cout << std::setw( 7 ) << std::setprecision( 3 ) << *vIt;
        }
        row += 1;
    }
}

/**
 *  \brief prints a Dune::FieldMatrix
 *
 *  or anything compatible in terms of Iterators
 *  \tparam should be Dune::FieldMatrix or compatible
 *  \param  arg
 *          Matrix to be printed
 *  \param  name
 *          name to be printed along
 **/
template < class T >
void printFieldMatrix( T& arg, std::string name )
{
    std::cout << "\nprinting " << name << " (Dune::FieldMatrix)";
    typedef typename T::ConstRowIterator
        RowIteratorType;
    typedef typename T::row_type::ConstIterator
        VectorInRowIteratorType;
    unsigned int row = 1;
    RowIteratorType rItEnd = arg.end();
    for ( RowIteratorType rIt = arg.begin(); rIt != rItEnd; ++rIt ) {
        std::cout << "\nrow " << row << ":";
        VectorInRowIteratorType vItEnd = rIt->end();
        for (   VectorInRowIteratorType vIt = rIt->begin(); vIt != vItEnd; ++vIt ) {
            std::cout << std::setw( 7 ) << std::setprecision( 3 ) << *vIt;
        }
        row += 1;
    }
}

} //end namepspace stuff


#endif
