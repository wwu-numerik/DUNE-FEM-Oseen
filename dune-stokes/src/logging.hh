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
 #include <iomanip>
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

        class LogStream //: virtual public std::ostream
        {
            protected:
                LogFlags loglevel_;
                int logflags_;
                std::stringstream buffer_;
                std::ofstream& logfile_;

            public:
                LogStream( LogFlags loglevel, int logflags, std::ofstream& file )
                    : loglevel_(loglevel), logflags_(logflags),
                      logfile_(file) {}
                ~LogStream(){}
                template < typename T >
                LogStream& operator << ( T in )
                {
                    buffer_ << in;
                    return *this;
                }

                //template < class Class >
                LogStream& operator << ( LogStream& ( *pf )(LogStream&) )
                {
                    buffer_ << pf;
                    return *this;
                }

                LogStream& operator << ( std::ostream& ( *pf )(std::ostream &) )
                {
                    if ( pf == (std::ostream& ( * )(std::ostream&))std::endl )
                    { //flush buffer into stream
                        if ( ( logflags_ & LOG_CONSOLE ) != 0 ) {
                            std::cout << buffer_.str();
                            std::cout.flush();
                        }
                        if ( ( logflags_ & LOG_FILE ) != 0 )
                            logfile_ << "\n" << TimeString() << buffer_.str();
                            logfile_.flush();
                    }
                    else
                        buffer_ << pf;

                    return *this;
                }
        };

        Logging( )
        {
        }

        ~Logging()
        {
            Stuff::safe_delete( stream_info );
            Stuff::safe_delete( stream_dbg );
            Stuff::safe_delete( stream_err );

            if ( ( logflags_ & LOG_FILE ) != 0 ) {
                logfile_ << '\n' << TimeString() << ": LOG END" << std::endl;
                logfile_.close();
            }
        }


        /** \brief setup loglevel, logfilename
            \param logflags any OR'd combination of flags
            \param logfile filename for log, can contain paths, but creation will fail if dir is non-existant
        **/
        void Create (unsigned int logflags = LOG_FILE | LOG_CONSOLE | LOG_ERR, std::string logfile = "dune-stokes.log" )
        {
            logflags_ = logflags;
            filename_ = logfile;
            if ( ( logflags_ & LOG_FILE ) != 0 ) {
                logfile_.open ( filename_.c_str() );
                assert( logfile_.is_open() );
            }
            flagmap[LOG_ERR] = logflags & ( LOG_ERR | LOG_DEBUG | LOG_INFO );
            flagmap[LOG_INFO] = logflags & ( LOG_ERR | LOG_DEBUG | LOG_INFO );
            flagmap[LOG_DEBUG] = logflags & ( LOG_ERR | LOG_DEBUG | LOG_INFO );
            stream_err = new LogStream( LOG_ERR, flagmap[LOG_ERR], logfile_ );
            stream_dbg = new LogStream( LOG_DEBUG, flagmap[LOG_DEBUG], logfile_ );
            stream_info = new LogStream( LOG_INFO, flagmap[LOG_INFO], logfile_ );
        }

         /** \name Log funcs for member-function pointers
         * \{
         */
        template < class Class >
        void LogDebug( void ( Class::*pf )(std::ostream&) , Class& c )
        {
            if ( ( logflags_ & LOG_DEBUG ) != 0 )
                Log( pf, c );
        }

        template < class Class >
        void LogInfo( void ( Class::*pf )(std::ostream&) , Class& c )
        {
            if ( ( logflags_ & LOG_INFO ) != 0 )
                Log( pf, c );
        }

        template < class Class >
        void LogErr( void ( Class::*pf )(std::ostream&) , Class& c )
        {
            if ( ( logflags_ & LOG_ERR ) != 0 )
                Log( pf, c );
        }
        /** \}
        */

         /** \name Log funcs for basic types/classes
         * \{
         */
        template < class Class >
        void LogDebug( Class c )
        {
            if ( ( logflags_ & LOG_DEBUG ) != 0 )
                Log( c );
        }

        template < class Class >
        void LogInfo( Class c )
        {
            if ( ( logflags_ & LOG_INFO ) != 0 )
                Log( c );
        }

        template < class Class >
        void LogErr( Class c )
        {
            if ( ( logflags_ & LOG_ERR ) != 0 )
                Log( c );
        }

        LogStream& Err() { assert( stream_err ); return *stream_err; }
        LogStream& Info() { assert( stream_info ); return *stream_info; }
        LogStream& Dbg() { assert( stream_dbg ); return *stream_dbg; }
        /** \}
        */

        static std::string TimeString()
        {
            const time_t cur_time = time( NULL );
            return ctime ( &cur_time );
        }

    private:
        std::string filename_;
        std::ofstream logfile_;
        typedef std::map<LogFlags,int> FlagMap;
        FlagMap flagmap_;
        int logflags_;
        LogStream* stream_err;
        LogStream* stream_dbg;
        LogStream* stream_info;

        template < class Class >
        void Log( void ( Class::*pf )(std::ostream&) , Class& c)
        {
            if ( ( logflags_ & LOG_CONSOLE ) != 0 )
                (c.*pf)( std::cout );
            if ( ( logflags_ & LOG_FILE ) != 0 )
                (c.*pf)( logfile_ );
        }

        template < class Class >
        void Log( Class c )
        {
            if ( ( logflags_ & LOG_CONSOLE ) != 0 )
                std::cout << c;
            if ( ( logflags_ & LOG_FILE ) != 0 )
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
