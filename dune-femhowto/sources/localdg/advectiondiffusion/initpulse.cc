/////////////////////////////////////////////////////////////////////
//
//  Rotating Pulse Problem: (page 13)
//
//  P. Bastian. Higher order discontinuous galerkin methods for flow and
//  transport in porous media. In E. Bänsch, editor, Challenges in
//  Scientific Computing - CISC 2002, number 35 in LNCSE, pages 1-22,
//  2003. 
//  
//  link:
//  http://hal.iwr.uni-heidelberg.de/~peter/Papers/cisc2002.pdf
//
/////////////////////////////////////////////////////////////////////
template <class GridType>
class U0 {
 public:
  enum { ConstantVelocity = false };
  enum { dimDomain = GridType::dimensionworld };  
  typedef FieldVector<double,dimDomain> DomainType;
  typedef FieldVector<double,1> RangeType;
  U0(double eps,bool diff_timestep=true) :
    spotmid_(0), epsilon(eps), diff_tstep(diff_timestep) 
  {
    spotmid_[0] = -0.5;
    myName = "RotatingPulse";
  }

  //! end time of simulation 
  double endtime() {
    return 0.25*M_PI;
  }
    
  // initialize function 
  void evaluate(const DomainType& arg, RangeType& res) const {
    evaluate(0,arg,res);
  }
  
  //! exact solution
  void evaluate(const DomainType& arg, const double t, RangeType& res) const 
  {
    evaluate(t,arg,res);
  }

  //! exact solution
  void evaluate(const double t, const DomainType& arg, RangeType& res) const 
  {
    double x = arg[0];
    double y = arg[1];

    double sig2 = 0.004; /* Siehe Paper P.Bastian Gl. 30 */
    double sig2PlusDt4 = sig2+(4.0*epsilon*t);
    double xq = ( x*cos(4.0*t) + y*sin(4.0*t)) - spotmid_[0];
    double yq = (-x*sin(4.0*t) + y*cos(4.0*t)) - spotmid_[1];

    res = (sig2/ (sig2PlusDt4) ) * exp (-( xq*xq + yq*yq ) / sig2PlusDt4 );
  }
  
  //! info for tex file 
  void printmyInfo(string filename)
  {
    std::ostringstream filestream;
    filestream << filename;

    std::ofstream ofs(filestream.str().c_str(), std::ios::app);
  
    ofs << "Problem: " << myName << "\n\n"
        << "Epsilon = " << epsilon << "\n\n";
    ofs.close();
  }

  // not valid, because ConstantVelocity = false
  DomainType constantVelocity() const 
  {
    DomainType v(0);
    return v;
  }

  //! velocity 
  void velocity(const double t, const DomainType& x, DomainType& v) const 
  {
    v[0] = -4.0*x[1];
    v[1] =  4.0*x[0];
    for(int i=2; i<DomainType :: dimension; ++i) v[i] = 0;
  }
  
  DomainType spotmid_;
  
  double epsilon;
  bool diff_tstep;
  string myName;
};
