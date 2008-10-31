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


/** \brief handles all logging
**/
class Logging
{
    public:

        enum LogFlags{
          LOG_ERR = 1,
          LOG_INFO = 2,
          LOG_DEBUG = 4,
          LOG_CONSOLE = 8,
          LOG_FILE = 16
        };

        Logging( )
        {
        }

        ~Logging()
        {

            if ( ( loglevel_ & LOG_FILE ) != 0 ) {
                logfile_ << TimeString() << ": LOG END" << std::endl;
                logfile_.close();
            }
        }

        void Create (unsigned int l = LOG_CONSOLE | LOG_ERR, std::string logfile = "log" )
        {
            loglevel_ = l;
            filename_ = logfile;
            if ( log_to_file_ ) {
                logfile_.open ( filename_.c_str() );
            }

        }


        template < class Class >
        void LogDebug( void ( Class::*pf )(std::ostream&) , Class& c )
        {
            if ( ( loglevel_ & LOG_DEBUG ) != 0 )
                Log( pf, c );
        }

        template < class Class >
        void LogInfo( void ( Class::*pf )(std::ostream&) , Class& c )
        {
            if ( ( loglevel_ & LOG_INFO ) != 0 )
                Log( pf, c );
        }

        template < class Class >
        void LogErr( void ( Class::*pf )(std::ostream&) , Class& c )
        {
            if ( ( loglevel_ & LOG_ERR ) != 0 )
                Log( pf, c );
        }

        template < class Class >
        void LogDebug( Class& c )
        {
            if ( ( loglevel_ & LOG_DEBUG ) != 0 )
                Log( c );
        }

        template < class Class >
        void LogInfo( Class& c )
        {
            if ( ( loglevel_ & LOG_INFO ) != 0 )
                Log( c );
        }

        template < class Class >
        void LogErr( Class& c )
        {
            if ( ( loglevel_ & LOG_ERR ) != 0 )
                Log( c );
        }

    private:
        bool log_to_file_;
        std::string filename_;
        std::ofstream logfile_;
        unsigned int loglevel_;

        std::string TimeString()
        {
            const time_t cur_time = time( NULL );
            return ctime ( &cur_time );
        }

        template < class Class >
        void Log( void ( Class::*pf )(std::ostream&) , Class& c)
        {
            if ( ( loglevel_ & LOG_CONSOLE ) != 0 )
                (c.*pf)( std::cout );
            if ( ( loglevel_ & LOG_FILE ) != 0 )
                (c.*pf)( logfile_ );
        }

        template < class Class >
        void Log( Class& c )
        {
            if ( ( loglevel_ & LOG_CONSOLE ) != 0 )
                std::cout << c;
            if ( ( loglevel_ & LOG_FILE ) != 0 )
                logfile_ << c;
        }
};

///global Logging instance
Logging& Logger ()
{
    static Logging log;
    return log;
}

#endif
