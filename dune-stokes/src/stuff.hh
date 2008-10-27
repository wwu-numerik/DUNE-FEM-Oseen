/**
 *  \file stuff.hh
 *  \brief  contains some stuff
 **/

namespace Stuff
{

/** \todo Please doc me! */
template < class ReturnType >
ReturnType fromString(const std::string& s)
{
    std::stringstream ss;
    ss << s;
    ReturnType r;
    ss >> r;
    return r;
}

/** \todo Please doc me! */
template < class ReturnType >
std::string toString(const ReturnType& s)
{
    std::stringstream ss;
    ss << s;
    std::string r;
    ss >> r;
    return r;
}

}
