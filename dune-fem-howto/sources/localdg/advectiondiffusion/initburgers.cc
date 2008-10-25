
template <class GridType>
class U0 {
public:
  enum { ConstantVelocity = true };
  enum { dimDomain = GridType::dimensionworld };  
  typedef FieldVector<double,dimDomain> DomainType;
  typedef FieldVector<double,1> RangeType;
  U0(double eps,bool diff_timestep=true) :
    velocity_(0), epsilon(eps), diff_tstep(diff_timestep) 
  {
    velocity_[0] = 1;
		myName = "Burgers-Diffusion";
	}
  
  double endtime() {
    return 0.4;
  }

  void constantVelocity(DomainType& v) const
  {
    v = velocity_;
  } 

  void velocity(const double t, const DomainType& x, DomainType& v) const
  {
    constantVelocity(v);
  }
  
  void evaluate(const DomainType& arg, RangeType& res) const {
    evaluate(0,arg,res);
  }
  
  void evaluate(const DomainType& arg, const double t, RangeType& res) const 
  {
    evaluate(t,arg,res);
  }
  
  void evaluate(const double t,const DomainType& arg, RangeType& res) const {
    if(epsilon < 1e-9)
      res = ((2.0*arg[0]-1.0)<0.0)? 1.0:-1.0;
    else
      res = -tanh((2.0*arg[0]-1.0)/(4.*epsilon));
  }
	
	void printmyInfo(string filename)
	{
  std::ostringstream filestream;
  filestream << filename;

  std::ofstream ofs(filestream.str().c_str(), std::ios::app);
	
	ofs << "Problem: " << myName << "\n\n"
	    << "Epsilon = " << epsilon << "\n\n"
			<< "Exact solution: $u(x,y,z,t):=\\displaystyle{-\\tanh\\left( \\frac{2x-1}{4\\varepsilon} \\right) }$\\\\\n\n";
	
	ofs.close();
	}
	
	
  DomainType velocity_;
  double epsilon;
  bool diff_tstep;
	string myName;
};
