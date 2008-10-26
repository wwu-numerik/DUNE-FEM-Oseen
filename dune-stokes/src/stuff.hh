namespace Stuff
{

/** \todo Please doc me! */
template < class ReturnType >
ReturnType fromString(const std::string& s) {
        std::stringstream ss;
        ss << s;
        ReturnType r;
        ss >> r;
        return r;
}

}
