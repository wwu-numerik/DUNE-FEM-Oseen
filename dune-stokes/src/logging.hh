/**
 *  \file logging.hh
 *  \brief  logging
 **/
 #ifndef LOGGING_HH_INCLUDED
 #define LOGGING_HH_INCLUDED

 #include <fstream>
 #include <ostream>
 #include <sstream>
 #include <ctime>
 #include "stuff.hh"
 #include "parametercontainer.hh"


/** \brief handles three diff logstreams
**/
class Logging
{
    public:

        enum LogLevel{
          MIN = 1,
          ALL = 2,
          ERR = 4
        };

        Logging( )
        {
        }


        ~Logging()
        {

            if ( log_to_file_ ) {
                logfile_ << TimeString() << ": LOG END" << std::endl;
                logfile_.close();
            }
        }

        void Create (LogLevel l, bool log_to_file = 0, std::string logfile = "log" )
        {
            loglevel_ = l;
            log_to_file_ = log_to_file;
            filename_ = logfile;
            if ( log_to_file_ ) {
                logfile_.open ( filename_.c_str() );
            }

        }

        template < class Class >
        void Log( void ( Class::*pf )(std::ostream&) , Class& c)
        {
            //if ( SHOULD I INDEED OUTPUT ANYTHING )
            (c.*pf)( std::cout );
            (c.*pf)( logfile_ );
        }

    private:
        bool log_to_file_;
        std::string filename_;
        std::ofstream logfile_;
        LogLevel loglevel_;

        std::string TimeString()
        {
            const time_t cur_time = time( NULL );
            return ctime ( &cur_time );
        }
};

///global Logging instance
Logging& Logger ()
{
    static Logging log;
    return log;
}

#endif
